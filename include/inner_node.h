//
// Copyright (c) Zhou Zhang.
// Implementation of PLIN inner nodes
//

#pragma once

#include "leaf_node.h"

class InnerNode{

#pragma pack(push)
#pragma pack(1)

    // The size of header should less than NODE_HEADER_SIZE (256B)
    struct InnerNodeMetadata {
        double slope = 0;
        float intercept = 0;
        uint32_t block_number = 0;
        uint64_t record_number = 0;
        _key_t first_key = FREE_FLAG;
        uint32_t level = 0;
        bool is_leaf = true;

        volatile uint32_t locked ;
        inline void get_lock(){
            uint32_t new_value = 0;
            uint32_t old_value = 0;
            do {
                while (true) {
                    old_value = __atomic_load_n(&locked, __ATOMIC_ACQUIRE);
                    if (!(old_value & lockSet)) {
                        old_value &= lockMask;
                        break;
                    }
                }
                new_value = old_value | lockSet;
            } while (!CAS(&locked, &old_value, new_value));
        }
        
        inline bool try_get_lock() {
            uint32_t v = __atomic_load_n(&locked, __ATOMIC_ACQUIRE);
            if (v & lockSet) {
                return false;
            }
            return true;
        }

        inline void release_lock() {
            uint32_t v = locked;
            __atomic_store_n(&locked, v + 1 - lockSet, __ATOMIC_RELEASE);
        }

        inline bool test_lock_set(uint32_t &version) const {
            version = __atomic_load_n(&locked, __ATOMIC_ACQUIRE);
            return (version & lockSet) != 0;
        }
    } inner_node;
    char unused[NODE_HEADER_SIZE - sizeof(InnerNodeMetadata)];

#pragma pack(pop)

    InnerSlot inner_slots[1]; 

public:

    constexpr static uint64_t InnerSlotsPerBlock = BLOCK_SIZE / (sizeof(InnerSlot));

    inline uint32_t predict_block (_key_t key, double slope, float intercept, uint32_t block_number) const {
        int64_t predicted_block = (key * slope + intercept + 0.5) / InnerSlotsPerBlock;
        if (predicted_block < 0)
            return 0;
        else if (predicted_block > block_number - 4)
            return block_number - 2;
        return predicted_block + 1;
    }

    InnerNode(InnerSlot& accelerator, uint64_t block_number, _key_t* keys, InnerSlot* nodes, uint64_t number, uint64_t start_pos, double slope, double intercept, uint32_t level, bool is_leaf, bool rebuild = false) {
        // assert((uint64_t)inner_slots - (uint64_t)&inner_node == NODE_HEADER_SIZE);
        inner_node.block_number = block_number;
        inner_node.slope = slope / INNER_NODE_INIT_RATIO;
        inner_node.intercept = (intercept - start_pos - keys[start_pos] * slope) / INNER_NODE_INIT_RATIO;
        inner_node.record_number = number;
        inner_node.first_key = keys[start_pos];
        inner_node.level = level;
        inner_node.is_leaf = is_leaf;
        do_flush(&inner_node, sizeof(inner_node));
        model_correction(inner_node.slope, inner_node.intercept, (inner_node.block_number - 3) * InnerSlotsPerBlock, keys[start_pos], keys[start_pos + number - 1]);
        // Init
        for (uint64_t i = 0; i < inner_node.block_number; ++i) {
            for (uint8_t j = 0; j < 8; j++) {
                inner_slots[i * InnerSlotsPerBlock + j].min_key = FREE_FLAG;
            }
        }
        // Model-based data placement
        for (uint64_t i = 0; i < number; ++i) {
            data_placement(keys[start_pos + i], nodes[start_pos + i]);
        }
        do_flush(inner_slots, block_number * BLOCK_SIZE);
        // Build accelerator
        accelerator.min_key = inner_node.first_key;
        accelerator.ptr = galc->relative(this);
        accelerator.slope = inner_node.slope;
        accelerator.intercept = inner_node.intercept;
        accelerator.set_block_number(inner_node.block_number);
        accelerator.set_skip(0);
        accelerator.set_type(1);
        accelerator.set_lock();
    }
    

    ~InnerNode() {
        if (inner_node.level > 1) {
            uint64_t slot_number = inner_node.block_number * InnerSlotsPerBlock;
            for (uint64_t slot = 0; slot < slot_number; ++slot) {
                if (inner_slots[slot].min_key != FREE_FLAG) {
                    delete reinterpret_cast<InnerNode*>(galc->absolute(inner_slots[slot].ptr));
                }
            }
        }
        galc->free(this);
    }

    InnerSlot * find_leaf_node (_key_t key, const InnerSlot * accelerator) {
        uint64_t slot = predict_block(key, accelerator->slope, accelerator->intercept, accelerator->block_number()) * InnerSlotsPerBlock;
        // Search left
        if (inner_slots[slot].min_key == FREE_FLAG || key < inner_slots[slot].min_key) {
            while (inner_slots[--slot].min_key == FREE_FLAG) {}
        }
        // Search right
        else {
            while (slot < inner_node.block_number * InnerSlotsPerBlock && inner_slots[slot].min_key != FREE_FLAG && key >= inner_slots[slot].min_key) {
                ++slot;
            }
            --slot;
        }
        if (inner_slots[slot].type()) {
            return reinterpret_cast<InnerNode*>(galc->absolute(inner_slots[slot].ptr))->find_leaf_node(key, &inner_slots[slot]);
        }
        else {
            return &inner_slots[slot];
        }
    }

    uint32_t upsert_node (InnerSlot& node, InnerSlot * accelerator) {//insert new node
        if(!accelerator ->try_get_lock())//
            return 0;
        accelerator ->get_lock();
        uint64_t slot = predict_block(node.min_key, accelerator->slope, accelerator->intercept, accelerator->block_number()) * InnerSlotsPerBlock;
        uint64_t predicted_slot = slot;
        // Search left
        if (inner_slots[slot].min_key == FREE_FLAG || node.min_key < inner_slots[slot].min_key) {
            while (inner_slots[--slot].min_key == FREE_FLAG) {}
        }
        // Search right
        else {
            while (slot < inner_node.block_number * InnerSlotsPerBlock && inner_slots[slot].min_key != FREE_FLAG && node.min_key >= inner_slots[slot].min_key) {
                ++slot;
            }
            --slot;
        }
        if (inner_slots[slot].type()) {
            accelerator ->release_lock();
            return reinterpret_cast<InnerNode*>(galc->absolute(inner_slots[slot].ptr))->upsert_node(node, &inner_slots[slot]);
        }
        // Upsert node
        else {
            if (inner_slots[slot].min_key == node.min_key) {
                LeafNode* old_node = reinterpret_cast<LeafNode*>(galc->absolute(inner_slots[slot].ptr));
                inner_slots[slot] = node;
                do_flush(&inner_slots[slot], sizeof(InnerSlot));
                mfence();
                galc->free(old_node);
                accelerator ->release_lock();
                return 1;
            }
            else if (inner_slots[predicted_slot + InnerSlotsPerBlock - 1].min_key == FREE_FLAG) {
                uint64_t target_slot = predicted_slot;
                if (slot >= predicted_slot) {
                    while (inner_slots[++target_slot].min_key <= node.min_key) {}
                }
                slot = predicted_slot + InnerSlotsPerBlock - 1;
                while (inner_slots[--slot].min_key == FREE_FLAG) {}
                for (++slot; slot > target_slot; --slot) {
                    inner_slots[slot] = inner_slots[slot - 1];
                }
                inner_slots[target_slot] = node;
                do_flush(&inner_slots[predicted_slot], InnerSlotsPerBlock * sizeof(InnerSlot));
                mfence();
                accelerator ->release_lock();
                return 2;
            }
            else {
                accelerator ->release_lock();
                return 3;
            }
        }
    }

    void data_placement (_key_t key, InnerSlot node) {
        uint64_t slot = predict_block(key, inner_node.slope, inner_node.intercept, inner_node.block_number) * InnerSlotsPerBlock;
        for (uint32_t i = 0; i < InnerSlotsPerBlock; ++i) {
            if (inner_slots[slot + i].min_key == FREE_FLAG) {                
                inner_slots[slot + i] = node;
                return;
            }
        }
        for (slot += InnerSlotsPerBlock; slot < inner_node.block_number * InnerSlotsPerBlock; slot++) {
            if (inner_slots[slot].min_key == FREE_FLAG) {                
                inner_slots[slot] = node;
                return;
            }
        }
        assert(slot < inner_node.block_number * InnerSlotsPerBlock);
        return;
    }

    LeafNode* get_leftmost_leaf(){
        uint64_t slot = 0;
        while(inner_slots[slot].min_key == FREE_FLAG) {
            ++slot;
        }
        if (inner_slots[slot].type() == 0) {
            return reinterpret_cast<LeafNode*>(galc->absolute(inner_slots[slot].ptr));
        }
        else {
            return reinterpret_cast<InnerNode*>(galc->absolute(inner_slots[slot].ptr))->get_leftmost_leaf();
        }
    }

};
