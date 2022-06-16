//
// Copyright (c) Zhou Zhang.
// Implementation of PLIN leaf nodes
//

#pragma once

#include "fast&fair.h"

class LeafNode{

#pragma pack(push)
#pragma pack(1)

    // The size of header should less than NODE_HEADER_SIZE (256B)
    struct LeafNodeMetadata {
        double slope = 0;
        float intercept = 0;
        uint32_t block_number = 0;
        LeafNode* prev = NULL;
        LeafNode* next = NULL;
        _key_t first_key = FREE_FLAG;
        uint64_t overflow_number = 0;
        uint64_t number_to_split = 0;
        // locked[0]: split lock; locked[1]: write lock; locked[2]: read lock
        volatile uint32_t locked = 0;
    } leaf_node;
    char unused[NODE_HEADER_SIZE - sizeof(LeafNodeMetadata)];

#pragma pack(pop)

    struct BlockHeader {
        volatile uint32_t locked;   // 1 bit lock & 31 bit block version
        volatile uint32_t local_version;    // 32 bit global version
        void * overflow_tree;

        inline void check_global_version(uint32_t global_version){
            if (global_version != local_version) {
                __atomic_store_n(&locked, 0, __ATOMIC_RELEASE);
                __atomic_store_n(&local_version, global_version, __ATOMIC_RELEASE);
            }
        }

        // lock a block
        inline void get_lock(uint32_t global_version){     
            check_global_version(global_version);
            uint32_t new_value = 0;
            uint32_t old_value = 0;
            do {
                while (true) {
                    old_value = __atomic_load_n(&locked, __ATOMIC_ACQUIRE);
                    if (!(old_value & lockSet)) {
                        break;
                    }
                }
                new_value = old_value | lockSet;
            } while (!CAS(&locked, &old_value, new_value));
        }
        
        // check lock
        inline bool check_lock(uint32_t global_version) {
            check_global_version(global_version);
            uint32_t v = __atomic_load_n(&locked, __ATOMIC_ACQUIRE);
            if (v & lockSet) {
                return false;
            }
            return true;
        }

        inline void release_lock() {
            __atomic_store_n(&locked, (locked + 1) & lockMask, __ATOMIC_RELEASE);
        }

        inline void release_lock_without_change() {
            __atomic_store_n(&locked, locked & lockMask, __ATOMIC_RELEASE);
        }

        // get the block version
        inline void get_version(uint32_t &version) const {
            do {
                version = __atomic_load_n(&locked, __ATOMIC_ACQUIRE);
            } while (version & lockSet);
        }

        // test whether the version has changed, if change, return true
        inline bool check_version_change(uint32_t old_version) const {
            auto value = __atomic_load_n(&locked, __ATOMIC_ACQUIRE);
            return (old_version != value);
        }
    };

    struct DataSlot {
        _key_t key;
        _payload_t payload;
    };

    union LeafSlot {
        BlockHeader header;
        DataSlot data; 
    };
    
    LeafSlot leaf_slots[1];


public:

    constexpr static uint64_t LeafSlotsPerBlock = BLOCK_SIZE / (sizeof(LeafSlot));

    constexpr static uint64_t LeafRealSlotsPerBlock = LeafSlotsPerBlock - 1;

    inline bool try_get_split_lock(){
        uint32_t new_value = 0;
        uint32_t old_value = __atomic_load_n(&leaf_node.locked, __ATOMIC_ACQUIRE);
        if (!(old_value & 1u)) {
            new_value = old_value | 1u;
            return CAS(&leaf_node.locked, &old_value, new_value);
        }
        return false;
    }

    inline void get_write_lock(){
        ADD(&leaf_node.locked,2);
    }

    inline void get_read_lock(){
        ADD(&leaf_node.locked,4);
    }
    
    inline bool check_write_lock(const InnerSlot * accelerator) {
        if (accelerator) {
            return accelerator->check_write_lock();
        }
        uint32_t v = __atomic_load_n(&leaf_node.locked, __ATOMIC_ACQUIRE);
        if (v & 2u) {
            return false;
        }
        return true;
    }

    inline bool check_read_lock(const InnerSlot * accelerator) {
        if (accelerator) {
            return accelerator->check_read_lock();
        }
        uint32_t v = __atomic_load_n(&leaf_node.locked, __ATOMIC_ACQUIRE);
        if (v & 4u) {
            return false;
        }
        return true;
    }

    inline void release_lock() {
        __atomic_store_n(&leaf_node.locked, 0, __ATOMIC_RELEASE);
    }

    static bool leafslot_cmp (LeafSlot a, LeafSlot b) {
        return a.data.key < b.data.key;
    }

    inline uint32_t predict_block (_key_t key, double slope, float intercept, uint32_t block_number) const {
        int64_t predicted_block = (key * slope + intercept + 0.5) / LeafSlotsPerBlock;
        if (predicted_block < 0)
            return 0;
        else if (predicted_block > block_number - 3)
            return block_number - 1;
        return predicted_block + 1;
    }

    LeafNode(InnerSlot& accelerator, uint64_t block_number, _key_t* keys, _payload_t* payloads, uint64_t number, uint64_t start_pos, double slope, double intercept, uint32_t global_version, LeafNode* prev = NULL, LeafNode* next = NULL) {
        // assert((uint64_t)leaf_slots - (uint64_t)&leaf_node == NODE_HEADER_SIZE);
        leaf_node.block_number = block_number;
        leaf_node.slope = slope * LeafSlotsPerBlock / (LEAF_NODE_INIT_RATIO * LeafRealSlotsPerBlock);
        leaf_node.intercept = (intercept - start_pos - keys[start_pos] * slope) * LeafSlotsPerBlock / (LEAF_NODE_INIT_RATIO * LeafRealSlotsPerBlock);
        leaf_node.first_key = keys[start_pos];
        leaf_node.number_to_split = block_number * LeafRealSlotsPerBlock * MAX_OVERFLOW_RATIO;
        leaf_node.prev = galc->relative(prev);
        leaf_node.next = galc->relative(next);
        do_flush(&leaf_node, sizeof(leaf_node));
        model_correction(leaf_node.slope, leaf_node.intercept, (leaf_node.block_number - 2) * LeafSlotsPerBlock, keys[start_pos], keys[start_pos + number - 1]);
        // Init
        for (uint64_t i = 0; i < leaf_node.block_number; ++i) {
            leaf_slots[i * LeafSlotsPerBlock].header.locked = 0;
            leaf_slots[i * LeafSlotsPerBlock].header.local_version = global_version;
            leaf_slots[i * LeafSlotsPerBlock].header.overflow_tree = NULL;
            for (uint8_t j = 1; j < LeafSlotsPerBlock; j++) {
                leaf_slots[i * LeafSlotsPerBlock + j].data.key = FREE_FLAG;
            }
        }
        // Model-based data placement
        for (uint64_t i = 0; i < number; ++i) {
            data_placement(keys[i + start_pos], payloads[i + start_pos]);
        }
        do_flush(leaf_slots, block_number * BLOCK_SIZE);
        // Build accelerator
        accelerator.min_key = leaf_node.first_key;
        accelerator.ptr = galc->relative(this);
        accelerator.slope = leaf_node.slope;
        accelerator.intercept = leaf_node.intercept;
        accelerator.set_block_number(leaf_node.block_number);
        accelerator.set_type(0);
        accelerator.init_lock();
    }

    ~LeafNode () {}
    
    void destroy () {
        galc->free(this);
    }

    void set_prev(LeafNode* prev) {
        if (prev)
            persist_assign(&leaf_node.prev, galc->relative(prev));
    }

    void set_next(LeafNode* next) {
        if (next)
            persist_assign(&leaf_node.next, galc->relative(next));
    }

    LeafNode* get_prev() {
        if (leaf_node.prev)
            return galc->absolute(leaf_node.prev);
        else
            return NULL;
    }

    LeafNode* get_next() {
        if (leaf_node.next)
            return galc->absolute(leaf_node.next);
        else
            return NULL;
    }

    _key_t get_min_key() {
        return leaf_node.first_key;
    }

    void get_info(_key_t& min_key, InnerSlot& accelerator) {
        min_key = leaf_node.first_key;
        accelerator.min_key = min_key;
        accelerator.ptr = galc->relative(this);
        accelerator.slope = leaf_node.slope;
        accelerator.intercept = leaf_node.intercept;
        accelerator.set_block_number(leaf_node.block_number);
        accelerator.set_type(0);
        accelerator.init_lock();
    }

    // 0 : false; 1 : true; 2 :  retry
    uint32_t find (_key_t key, _payload_t & payload, uint32_t global_version, const InnerSlot * accelerator = NULL) {
        if (!check_read_lock(accelerator)) {
            return 2;
        }
        uint32_t block, slot;
        if (accelerator) {
            block = predict_block(key, accelerator->slope, accelerator->intercept, accelerator->block_number());
            if (block >= accelerator->block_number() - 2 && leaf_node.next && get_next()->get_min_key() <= key) {
                return get_next()->find(key, payload, global_version);
            }
            slot = block * LeafSlotsPerBlock;
        }
        else {
            block = predict_block(key, leaf_node.slope, leaf_node.intercept, leaf_node.block_number);
            if (block >= leaf_node.block_number - 2 && leaf_node.next && get_next()->get_min_key() <= key) {
                return get_next()->find(key, payload, global_version);
            }
            slot = block * LeafSlotsPerBlock;
        }

        uint32_t version;
        bool version_changed;
        do {
            // Wait the lock released
            while(!leaf_slots[slot].header.check_lock(global_version)) {}
            leaf_slots[slot].header.get_version(version);
            version_changed = false;
            for (uint32_t i = 1; i < LeafSlotsPerBlock; ++i) {
                if (key == leaf_slots[slot + i].data.key) {
                    payload = leaf_slots[slot + i].data.payload;
                    if (leaf_slots[slot].header.check_version_change(version)) {
                        version_changed = true;
                        break;
                    }
                    else
                        return 1;
                }
            }
        } while (version_changed);

        // Find in overflow block
        if (leaf_slots[slot].header.overflow_tree) {
            btree * overflow_tree = new btree(&leaf_slots[slot].header.overflow_tree, false);
            if(overflow_tree->find(key, payload)){
                return 1;
            }
        }
        return 0;
    }

    // rq_flag = 0: start; = 1: scan; = 2: stop
    void range_query(_key_t lower_bound, _key_t upper_bound, std::vector<std::pair<_key_t, _payload_t>>& answers, uint8_t rq_flag, InnerSlot * accelerator = NULL) {
        uint32_t block = 0;
        uint32_t block_number = 0;
        if (rq_flag == 0) {
            if (accelerator) {
                block = predict_block(lower_bound, accelerator->slope, accelerator->intercept, accelerator->block_number());
                if (block >= accelerator->block_number() - 2 && leaf_node.next && get_next()->get_min_key() <= lower_bound) {
                    return get_next()->range_query(lower_bound, upper_bound, answers, 0);
                }
                block_number = accelerator->block_number();
                rq_flag = 1;
            }
            else {
                block = predict_block(lower_bound, leaf_node.slope, leaf_node.intercept, leaf_node.block_number);
                if (block >= leaf_node.block_number - 2 && leaf_node.next && get_next()->get_min_key() <= lower_bound) {
                    return get_next()->range_query(lower_bound, upper_bound, answers, 0);
                }
                block_number = leaf_node.block_number;
                rq_flag = 1;
            }
        }
        else {
            block_number = leaf_node.block_number;
        }
        while (block < block_number && rq_flag == 1) {
            uint64_t slot = block * LeafSlotsPerBlock;
            for (uint32_t i = 1; i < LeafSlotsPerBlock; ++i) {
                if (leaf_slots[slot + i].data.key >= lower_bound && leaf_slots[slot + i].data.key <= upper_bound) {
                    answers.emplace_back(leaf_slots[slot + i].data.key, leaf_slots[slot + i].data.payload);
                }
                else if (leaf_slots[slot + i].data.key > upper_bound && leaf_slots[slot + i].data.key != FREE_FLAG && leaf_slots[slot + i].data.key != DELETE_FLAG) {
                    rq_flag = 2;
                }
            }
            if (leaf_slots[slot].header.overflow_tree) {
                btree* overflow_tree = new btree(&leaf_slots[slot].header.overflow_tree, false);
                overflow_tree->range_query(lower_bound, upper_bound, answers);
            }
            ++block;
        }
        if (rq_flag == 1 && leaf_node.next) {
            get_next()->range_query(lower_bound, upper_bound, answers, 1);
        }
    }

    // ret = 1 : update in a slot; ret = 2 : insert in a free slot; ret = 3 : update in overflow block; ret = 4 : insert in overflow block; 
    // ret = 5 : insert in overflow block & need to split; ret = 6 : insert in overflow block & need to split orphan node; ret = 7 : the node is locked
    uint32_t upsert (_key_t key, _payload_t payload, uint32_t global_version, LeafNode*& leaf_to_split, const InnerSlot * accelerator = NULL) {

        if (!check_read_lock(accelerator)) {
            return 7;
        }

        uint32_t block, slot;
        if (accelerator) {
            block = predict_block(key, accelerator->slope, accelerator->intercept, accelerator->block_number());
            if (block >= accelerator->block_number() - 2 && leaf_node.next && get_next()->get_min_key() <= key) {
                return get_next()->upsert(key, payload, global_version, leaf_to_split);
            }
            slot = block * LeafSlotsPerBlock;
        }
        else {
            block = predict_block(key, leaf_node.slope, leaf_node.intercept, leaf_node.block_number);
            if (block >= leaf_node.block_number - 2 && leaf_node.next && get_next()->get_min_key() <= key) {
                return get_next()->upsert(key, payload, global_version, leaf_to_split);
            }
            slot = block * LeafSlotsPerBlock;
        }
        
        if (!check_write_lock(accelerator)) {
            leaf_slots[slot].header.release_lock_without_change();
            return 7;
        } 

        // Update
        leaf_slots[slot].header.get_lock(global_version);
        for (uint32_t i = 1; i < LeafSlotsPerBlock; ++i) {
            if (leaf_slots[slot + i].data.key == key) {
                leaf_slots[slot + i].data.payload = payload;
                leaf_slots[slot].header.release_lock();
                do_flush(&leaf_slots[slot + i], sizeof(LeafSlot));
                mfence();
                return 1;
            }
        }

        uint32_t deleted_slot = 0;
        // Insert
        for (uint32_t i = 1; i < LeafSlotsPerBlock; ++i) {
            if (leaf_slots[slot + i].data.key == FREE_FLAG) {      
                leaf_slots[slot + i].data.payload = payload;
                leaf_slots[slot + i].data.key = key;
                leaf_slots[slot].header.release_lock();
                do_flush(&leaf_slots[slot + i], sizeof(LeafSlot));
                mfence();
                return 2;
            }
            else if (leaf_slots[slot + i].data.key == DELETE_FLAG) {
                deleted_slot = i;
            }
        }

        // Upsert into overflow block
        btree* overflow_tree = new btree(&leaf_slots[slot].header.overflow_tree, !leaf_slots[slot].header.overflow_tree);
        uint32_t flag = overflow_tree->upsert(key, payload, deleted_slot);
        if (flag == 2) {
            leaf_slots[slot + deleted_slot].data.payload = payload;
            leaf_slots[slot + deleted_slot].data.key = key;
            leaf_slots[slot].header.release_lock();
            do_flush(&leaf_slots[slot + deleted_slot], sizeof(LeafSlot));
            mfence();
            return 2;
        }
        leaf_slots[slot].header.release_lock();
        if (flag == 4) {
            ++leaf_node.overflow_number;
            if (leaf_node.overflow_number > leaf_node.number_to_split) {
                leaf_to_split = this;
                if (accelerator)
                    flag = 5;
                else        // Orphan node need split
                    flag = 6;
            }
        }
        return flag;
    }

    // 0 : false; 1 : in learned node; 2 : in overflow tree 3 : retry
    uint32_t remove (_key_t key, uint32_t global_version, const InnerSlot * accelerator = NULL) {

        if (!check_read_lock(accelerator)) {
            return 3;
        }

        uint32_t block, slot;
        if (accelerator) {
            block = predict_block(key, accelerator->slope, accelerator->intercept, accelerator->block_number());
            if (block >= accelerator->block_number() - 2 && leaf_node.next && get_next()->get_min_key() <= key) {
                return get_next()->remove(key, global_version);
            }
            slot = block * LeafSlotsPerBlock;
        }
        else {
            block = predict_block(key, leaf_node.slope, leaf_node.intercept, leaf_node.block_number);
            if (block >= leaf_node.block_number - 2 && leaf_node.next && get_next()->get_min_key() <= key) {
                return get_next()->remove(key, global_version);
            }
            slot = block * LeafSlotsPerBlock;
        }
        
        if (!check_write_lock(accelerator)) {
            leaf_slots[slot].header.release_lock_without_change();
            return 3;
        } 

        leaf_slots[slot].header.get_lock(global_version);
        for (uint32_t i = 1; i < LeafSlotsPerBlock; ++i) {
            if (leaf_slots[slot + i].data.key == key) {
                leaf_slots[slot + i].data.key = leaf_slots[slot].header.overflow_tree ? DELETE_FLAG : FREE_FLAG;
                leaf_slots[slot].header.release_lock();
                do_flush(&leaf_slots[slot + i], sizeof(LeafSlot));
                mfence();
                return 1;
            }
        }

        if (leaf_slots[slot].header.overflow_tree) {
            btree* overflow_tree = new btree(&leaf_slots[slot].header.overflow_tree, false);
            if (overflow_tree->remove(key)) {
                leaf_slots[slot].header.release_lock();
                --leaf_node.overflow_number;
                return 2;
            }
        }

        leaf_slots[slot].header.release_lock_without_change();
        return 0;
    }

    void data_placement (_key_t key, _payload_t payload) {
        uint64_t slot = predict_block(key, leaf_node.slope, leaf_node.intercept, leaf_node.block_number) * LeafSlotsPerBlock;
        // Insert
        for (uint32_t i = 1; i < LeafSlotsPerBlock; ++i) {
            if (leaf_slots[slot + i].data.key == FREE_FLAG) {                
                leaf_slots[slot + i].data.payload = payload;
                leaf_slots[slot + i].data.key = key;
                return;
            }
        }
        // Insert into overflow block
        btree* overflow_tree = new btree(&leaf_slots[slot].header.overflow_tree, !leaf_slots[slot].header.overflow_tree);
        overflow_tree->insert(key, payload);
        ++leaf_node.overflow_number;
        return;
    }

    void get_data (std::vector<_key_t>& keys, std::vector<_payload_t>& payloads, uint32_t global_version) {      
        for (uint32_t block = 0; block < leaf_node.block_number; ++block) {
            leaf_slots[block*LeafSlotsPerBlock].header.check_lock(global_version);      // Ensure no write on this block
            std::sort(leaf_slots + (block * LeafSlotsPerBlock + 1), leaf_slots + ((block + 1) * LeafSlotsPerBlock), leafslot_cmp);
            if (leaf_slots[block * LeafSlotsPerBlock].header.overflow_tree) {
                btree* overflow_tree = new btree(&leaf_slots[block * LeafSlotsPerBlock].header.overflow_tree, false);
                std::vector<_key_t> overflowed_keys;
                std::vector<_payload_t> overflowed_payloads;
                overflow_tree->get_data(overflowed_keys, overflowed_payloads);
                uint32_t i = 1, j = 0;
                while (i < LeafSlotsPerBlock && (leaf_slots[block * LeafSlotsPerBlock + i].data.key == DELETE_FLAG || leaf_slots[block * LeafSlotsPerBlock + i].data.key == FREE_FLAG)) {
                    ++i;
                }
                while (i < LeafSlotsPerBlock || j < overflowed_keys.size()) {
                    if (i == LeafSlotsPerBlock || j < overflowed_keys.size() && leaf_slots[block * LeafSlotsPerBlock + i].data.key > overflowed_keys[j]) {
                        keys.emplace_back(overflowed_keys[j]);
                        payloads.emplace_back(overflowed_payloads[j++]);
                    }
                    else {
                        if(leaf_slots[block * LeafSlotsPerBlock + i].data.key < leaf_node.first_key)
                            int k=0;
                        keys.emplace_back(leaf_slots[block * LeafSlotsPerBlock + i].data.key);
                        payloads.emplace_back(leaf_slots[block * LeafSlotsPerBlock + i].data.payload);
                        ++i;
                        if (i < LeafSlotsPerBlock && (leaf_slots[block * LeafSlotsPerBlock + i].data.key == DELETE_FLAG || leaf_slots[block * LeafSlotsPerBlock + i].data.key == FREE_FLAG)) {
                            i = LeafSlotsPerBlock;
                        }
                    }
                }
            }
            else {
                for (uint32_t i = 1; i < LeafSlotsPerBlock; ++i) {
                    if (leaf_slots[block * LeafSlotsPerBlock + i].data.key != FREE_FLAG && leaf_slots[block * LeafSlotsPerBlock + i].data.key != DELETE_FLAG) {
                        keys.emplace_back(leaf_slots[block * LeafSlotsPerBlock + i].data.key);
                        payloads.emplace_back(leaf_slots[block * LeafSlotsPerBlock + i].data.payload);
                    }
                }
            }
        }
    }

};

