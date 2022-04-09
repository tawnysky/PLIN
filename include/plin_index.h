//
// Copyright (c) Zhou Zhang.
// Implementation of PLIN index structure
//

#pragma once

#include "inner_node.h"
#include "piecewise_linear_model.h"

class PlinIndex {
    typedef PlinIndex SelfType;
    // For storing models, to build nodes
    struct Segment {
        _key_t first_key;
        double slope;
        double intercept;
        uint64_t number;
        explicit Segment(const typename OptimalPiecewiseLinearModel<_key_t, size_t>::CanonicalSegment &cs)
                : first_key(cs.get_first_x()),
                number(cs.get_number()) {
            auto [cs_slope, cs_intercept] = cs.get_floating_point_segment(first_key);
            slope = cs_slope;
            intercept = cs_intercept;
        }
    };

    struct split_log {
        LeafNode* leaf_to_split = NULL;
        LeafNode* left_sibling = NULL;
        LeafNode* right_sibling = NULL;
        LeafNode* left_node = NULL;
        LeafNode* right_node = NULL;
        // locked[0]: write lock, locked[1]: valid, locked[2]: orphan
        volatile uint32_t locked = 0;

        inline bool try_get_lock(){
            uint32_t new_value = 0;
            uint32_t old_value = __atomic_load_n(&locked, __ATOMIC_ACQUIRE);
            if (!(old_value & 1u)) {
                new_value = old_value | 1u;
                return CAS(&locked, &old_value, new_value);
            }
            return false;
        }
        
        inline bool check_lock() {
            uint32_t v = __atomic_load_n(&locked, __ATOMIC_ACQUIRE);
            if (v & 1u) {
                return false;
            }
            return true;
        }

        inline bool check_valid() {
            uint32_t v = __atomic_load_n(&locked, __ATOMIC_ACQUIRE);
            if (v & 2u) {
                return true;
            }
            return false;
        }

        inline bool check_orphan() {
            uint32_t v = __atomic_load_n(&locked, __ATOMIC_ACQUIRE);
            if (v & 4u) {
                return true;
            }
            return false;
        }

        inline void set_valid() {
            ADD(&locked, 2);
        }

        inline void set_orphan() {
            ADD(&locked, 4);
        }

        inline void release_lock() {
            __atomic_store_n(&locked, 0, __ATOMIC_RELEASE);
        }
    };

    struct plin_metadata {
        _key_t min_key = 0;
        _key_t max_key = 0;
        btree * left_buffer = NULL;
        btree * right_buffer = NULL;
        uint32_t leaf_number = 0;
        uint32_t orphan_number = 0;
        uint32_t left_buffer_number = 0;
        uint32_t right_buffer_number = 0;
        uint32_t root_number = 0;
        uint32_t level = 0;
        // Used for rebuilding, no SMO in the process of rebuilding, 1 bit lock & 31 bit read counter
        volatile uint32_t smo_lock = 0;
        uint32_t global_version = 0;
        uint32_t rebuilding = 0;
        plin_metadata* old_plin_ = NULL;  
        InnerSlot roots[ROOT_SIZE];
        split_log logs[LOG_NUMBER];
        
        inline bool get_write_lock() {
            uint32_t v;
            do {
                v = __atomic_load_n(&smo_lock, __ATOMIC_ACQUIRE);
                if (v & lockSet)
                    return false;
            } while (!CAS(&smo_lock, &v,  v | lockSet));

            // Wait until the readers all exit the critical section
            do {
                v = __atomic_load_n(&smo_lock, __ATOMIC_ACQUIRE);
            } while (v & lockMask);
            return true;
        }

        inline bool try_get_read_lock() {
            uint32_t v;
            do {
                v = __atomic_load_n(&smo_lock, __ATOMIC_ACQUIRE);
                if (v & lockSet)
                    return false;
            } while (!CAS(&smo_lock, &v,  v + 1));
            return true;
        }

        inline void release_read_lock() { SUB(&smo_lock, 1); }

        inline void release_write_lock() {
            __atomic_store_n(&smo_lock, 0, __ATOMIC_RELEASE);
        }
    };

    plin_metadata * plin_;

public:

    PlinIndex(std::string path, std::string id = "plin", bool recovery = false) {
        if (!recovery) {
            galc = new PMAllocator(path.c_str(), false, id.c_str());
            plin_ = (plin_metadata *) galc->get_root(sizeof(plin_metadata));
            plin_->left_buffer = galc->relative(new btree());   
            plin_->right_buffer = galc->relative(new btree());  
            do_flush(plin_, sizeof(plin_metadata));
            mfence();
        }
        else {
            galc = new PMAllocator(path.c_str(), true, id.c_str());
            plin_ = (plin_metadata *) galc->get_root(sizeof(plin_metadata));
            plin_->global_version++;
            if (plin_->rebuilding > 0) {
                plin_metadata * old_plin_ = galc->absolute(plin_->old_plin_);
                if (plin_->rebuilding > 1) {
                    plin_->root_number = old_plin_->root_number;
                    for (uint32_t i = 0; i < plin_->root_number; ++i)
                        plin_->roots[i] = old_plin_->roots[i];
                    for (uint32_t i = plin_->root_number; i < ROOT_SIZE; ++i)
                        plin_->roots[i].min_key = FREE_FLAG;
                    plin_->level = old_plin_->level;
                    plin_->orphan_number = old_plin_->orphan_number;
                }
                plin_->rebuilding = 0;
                galc->free(old_plin_);
                plin_->old_plin_ = NULL;
                do_flush(plin_, sizeof(plin_metadata));
                mfence();
            }
            for (uint32_t i = 0; i < LOG_NUMBER; ++i) {
                if (!plin_->logs[i].check_lock()) {
                    LeafNode * left_sibling = galc->absolute(plin_->logs[i].left_sibling);
                    LeafNode * right_sibling = galc->absolute(plin_->logs[i].right_sibling);
                    if (plin_->logs[i].check_valid()) {
                        LeafNode * left_node = galc->absolute(plin_->logs[i].left_node);
                        LeafNode * right_node = galc->absolute(plin_->logs[i].right_node);
                        if (left_sibling) {
                            left_sibling->set_next(left_node);
                        }
                        if (right_sibling) {
                            right_sibling->set_prev(right_node);
                        }
                        if(!plin_->logs[i].check_orphan()) {
                            _key_t first_key;
                            InnerSlot accelerator;
                            left_node->get_info(first_key, accelerator);
                            upsert_node(accelerator);
                        }
                    }
                    else {
                        LeafNode * leaf_to_split = galc->absolute(plin_->logs[i].leaf_to_split);
                        leaf_to_split->release_lock();
                        if (!plin_->logs[i].check_orphan()) {
                            uint32_t i = 0;
                            _key_t key = leaf_to_split->get_min_key();
                            while (i < plin_->root_number && key >= plin_->roots[i].min_key) {
                                ++i;
                            }
                            InnerSlot * accelerator = &plin_->roots[--i];
                            if (accelerator->type()) {
                                accelerator = reinterpret_cast<InnerNode*>(galc->absolute(accelerator->ptr))->find_leaf_node(key, accelerator);
                            }
                            accelerator->release_lock();
                        }
                    }
                    if (left_sibling) {
                        left_sibling->release_lock();
                    }
                    if (right_sibling) {
                        right_sibling->release_lock();
                    }
                    plin_->logs[i].release_lock();
                }
            }
            plin_->release_write_lock();
        }
    }

    ~PlinIndex() {}

    void destroy() {
        LeafNode * cur = get_leftmost_leaf();
        LeafNode * next = cur->get_next();
        cur->destroy();
        while (next) {
            cur = next;
            next = cur->get_next();
            cur->destroy();
        }
        for (uint32_t i = 0; i < plin_->root_number; i++) {
            galc->absolute((InnerNode*)(plin_->roots[i].ptr))->destroy();
        }
    }

    void bulk_load(_key_t* keys, _payload_t* payloads, uint64_t number) {
        // Make segmentation for leaves
        std::vector<Segment> segments;
        auto in_fun = [keys](auto i) {
            return std::pair<_key_t, size_t>(keys[i],i);
        };
        auto out_fun = [&segments](auto cs) { segments.emplace_back(cs); };
        uint64_t last_n = make_segmentation(number, EPSILON_LEAF_NODE, in_fun, out_fun);

        // Build leaf nodes
        uint64_t start_pos = 0;
        auto first_keys = new _key_t[last_n];
        auto leaf_nodes = new LeafNode*[last_n];
        auto accelerators = new InnerSlot[last_n];
        LeafNode* prev = NULL;
        for(uint64_t i = 0; i < last_n; ++i) {
            uint64_t block_number = (uint64_t) segments[i].number / (LEAF_NODE_INIT_RATIO * LeafNode::LeafRealSlotsPerBlock) + 3;
            uint64_t node_size_in_byte = block_number * BLOCK_SIZE + NODE_HEADER_SIZE;
            first_keys[i] = segments[i].first_key;
            leaf_nodes[i] = new (galc->malloc(node_size_in_byte)) LeafNode(accelerators[i], block_number, keys, payloads, segments[i].number, start_pos, segments[i].slope, segments[i].intercept, plin_->global_version, prev);
            start_pos += segments[i].number;
            prev = leaf_nodes[i];
        }
        for(uint64_t i = 0; i < last_n - 1; ++i){
            leaf_nodes[i]->set_next(leaf_nodes[i + 1]);
        }
           
        plin_->leaf_number = last_n;
        delete [] leaf_nodes;

        // Build inner nodes recursively
        uint32_t level = 0;
        uint64_t offset = 0;
        auto accelerators_tmp = accelerators;
        auto first_keys_tmp = first_keys;
        auto in_fun_rec = [&first_keys](auto i){
            return std::pair<_key_t, size_t>(first_keys[i],i);
        };
        while (level == 0 || last_n > ROOT_SIZE) {
            level++;
            offset += last_n;
            last_n = make_segmentation(last_n, EPSILON_INNER_NODE, in_fun_rec, out_fun);
            // std::cout<<"Number of inner nodes: "<<last_n<<std::endl;
            start_pos = 0;
            first_keys = new _key_t[last_n];
            accelerators = new InnerSlot[last_n];
            for (uint64_t i = 0; i < last_n; ++i) {
                uint64_t block_number = (uint64_t) segments[offset + i].number / (INNER_NODE_INIT_RATIO * InnerNode::InnerSlotsPerBlock) + EPSILON_INNER_NODE / InnerNode::InnerSlotsPerBlock + 3;
                uint64_t node_size_in_byte = block_number * BLOCK_SIZE + NODE_HEADER_SIZE;
                first_keys[i] = segments[offset + i].first_key;
                new (galc->malloc(node_size_in_byte)) InnerNode(accelerators[i], block_number, first_keys_tmp, accelerators_tmp, segments[offset + i].number, start_pos, segments[offset + i].slope, segments[offset + i].intercept, level);
                start_pos += segments[offset + i].number;
            }
            delete [] accelerators_tmp;
            delete [] first_keys_tmp;
            accelerators_tmp = accelerators;
            first_keys_tmp = first_keys;
        }

        // Build root
        plin_->min_key = keys[0];
        plin_->max_key = keys[number - 1];
        plin_->root_number = last_n;
        plin_->level = level;
        plin_->global_version = 0;
        for (uint32_t i = 0; i < ROOT_SIZE; ++i) {
            plin_->roots[i].min_key = FREE_FLAG;
        }
        for (uint32_t i = 0; i < last_n; ++i) {
            plin_->roots[i] = accelerators_tmp[i];
        }
        do_flush(plin_, sizeof(plin_metadata));
        mfence();
        delete [] accelerators_tmp;
        delete [] first_keys_tmp;
    }

    bool find(_key_t key, _payload_t & payload) {
        if (key >= plin_->min_key && key <= plin_->max_key) {
            do {
                uint32_t i = 0;
                while (i < plin_->root_number && key >= plin_->roots[i].min_key) {
                    ++i;
                }
                InnerSlot * accelerator = &plin_->roots[--i];
                if (accelerator->type()) {
                    accelerator = reinterpret_cast<InnerNode*>(galc->absolute(accelerator->ptr))->find_leaf_node(key, accelerator);
                }
                if (accelerator->check_read_lock()) {
                    uint32_t ret = reinterpret_cast<LeafNode*>(galc->absolute(accelerator->ptr))->find(key, payload, plin_->global_version, accelerator);
                    if (ret == 1) {
                        return true;
                    }
                    else if (ret == 0) {
                        return false;
                    }
                }
            } while (true);
        }
        // Find in buffer
        else if (key < plin_->min_key) {
            return galc->absolute(plin_->left_buffer)->find(key, payload);
        }
        else {
            return galc->absolute(plin_->right_buffer)->find(key, payload);
        }
    }

    void range_query(_key_t lower_bound, _key_t upper_bound, std::vector<std::pair<_key_t, _payload_t>>& answers) {
        if (lower_bound < plin_->min_key) {
            galc->absolute(plin_->left_buffer)->range_query(lower_bound, upper_bound, answers);
            lower_bound = plin_->min_key;
        }
        if (upper_bound > plin_->max_key) {
            galc->absolute(plin_->right_buffer)->range_query(lower_bound, upper_bound, answers);
            upper_bound = plin_->max_key;
        }
        if (upper_bound >= plin_->min_key && lower_bound <= plin_->max_key) {
            uint32_t i = 0;
            while (i < plin_->root_number && lower_bound >= plin_->roots[i].min_key) {
                ++i;
            }
            InnerSlot * accelerator = &plin_->roots[--i];
            if (accelerator->type()) {
                accelerator = reinterpret_cast<InnerNode*>(galc->absolute(accelerator->ptr))->find_leaf_node(lower_bound, accelerator);
            }
            reinterpret_cast<LeafNode*>(galc->absolute(accelerator->ptr))->range_query(lower_bound, upper_bound, answers, 0, accelerator);
        }
    }

    void upsert(_key_t key, _payload_t payload) {
        if (key > plin_->min_key && key < plin_->max_key) {
            uint32_t ret;
            do {
                uint32_t i = 0;
                while (i < plin_->root_number && key >= plin_->roots[i].min_key) {
                    ++i;
                }
                InnerSlot * accelerator = &plin_->roots[--i];
                
                if (accelerator->type()) {
                    accelerator = reinterpret_cast<InnerNode*>(galc->absolute(accelerator->ptr))->find_leaf_node(key, accelerator);
                }
                LeafNode * leaf_to_split;
                // ret = 1 : update in a slot; ret = 2 : insert in a free slot; ret = 3 : update in overflow block; ret = 4 : insert in overflow block; 
                // ret = 5 : insert in overflow block & need to split; ret = 6 : insert in overflow block & need to split orphan node; ret = 7 : the node is locked
                ret = reinterpret_cast<LeafNode*>(galc->absolute(accelerator->ptr))->upsert(key, payload, plin_->global_version, leaf_to_split, accelerator);
                // Split leaf node
                if (ret == 5) {
                    split(leaf_to_split, accelerator);
                }
                else if (ret == 6) {
                    split(leaf_to_split, NULL);
                }
            } while (ret == 7);
        }
        // Upsert in buffer
        else if (key < plin_->min_key) {
            uint32_t ret = galc->absolute(plin_->left_buffer)->upsert(key, payload, 0);
            if (ret == 4) {
                // TODO: merge buffer
                if (++plin_->left_buffer_number > MAX_BUFFER) {}
            }
        }
        else {
            uint32_t ret = galc->absolute(plin_->right_buffer)->upsert(key, payload, 0);
            if (ret == 4) {
                // TODO: merge buffer
                if (++plin_->right_buffer_number > MAX_BUFFER) {}
            }
        }
    }

    void remove (_key_t key) {
        if (key > plin_->min_key && key < plin_->max_key) {
            uint32_t ret;
            do {
                uint32_t i = 0;
                while (i < plin_->root_number && key >= plin_->roots[i].min_key) {
                    ++i;
                }
                InnerSlot * accelerator = &plin_->roots[--i];
                if (accelerator->type()) {
                    accelerator = reinterpret_cast<InnerNode*>(galc->absolute(accelerator->ptr))->find_leaf_node(key, accelerator);
                }
                ret = reinterpret_cast<LeafNode*>(galc->absolute(accelerator->ptr))->remove(key, plin_->global_version, accelerator);
            } while (ret == 3);
        }
        else if (key < plin_->min_key) {
            if (galc->absolute(plin_->left_buffer)->remove(key)) {
                --plin_->left_buffer_number;
            }
        }
        else {
            if (galc->absolute(plin_->right_buffer)->remove(key)) {
                --plin_->right_buffer_number;
            }
        }
    }

    InnerSlot* get_parent(const InnerSlot * node){
        uint32_t i = 0;
        while (i < plin_->root_number && node->min_key >= plin_->roots[i].min_key) {
            ++i;
        }
        --i;
        return &plin_->roots[i];
    }

    // Split & insert nodes
    void split (LeafNode * leaf_to_split, InnerSlot * accelerator) { 

        LeafNode * left_sibling = leaf_to_split->get_prev();;
        LeafNode * right_sibling = leaf_to_split->get_next();

        // Check wether the index is rebuilding, no smo in rebuilding process
        if(!plin_->try_get_read_lock()) {
            return;
            if (!leaf_to_split->try_get_split_lock()) {
                plin_->release_read_lock();
                return;
                // Get split lock of the node, the prev node, and the next node
                if ((!left_sibling) || (!left_sibling->try_get_split_lock())) {
                    leaf_to_split->release_lock();
                    plin_->release_read_lock();
                    return;
                    if ((!right_sibling) || (!right_sibling->try_get_split_lock())) {
                        if (left_sibling)
                            left_sibling->release_lock();
                        leaf_to_split->release_lock();
                        plin_->release_read_lock();
                        return;
                    }
                }
            }
        }

        if (accelerator) {
            accelerator->get_write_lock();
        }
        else {
            leaf_to_split->get_write_lock();
        }

        // Write log
        uint32_t log_number = LOG_NUMBER;
        do {
            for (uint32_t i = 0; i < LOG_NUMBER; ++i) {
                if (plin_->logs[i].try_get_lock()) {
                    log_number = i;
                    break;
                }
            }
        } while (log_number == LOG_NUMBER);
        plin_->logs[log_number].leaf_to_split = galc->relative(leaf_to_split);
        if (left_sibling)
            plin_->logs[log_number].left_sibling = galc->relative(left_sibling);
        else
            plin_->logs[log_number].left_sibling = NULL;
        if (right_sibling)
            plin_->logs[log_number].right_sibling = galc->relative(right_sibling);
        else
            plin_->logs[log_number].right_sibling = NULL;
        
        std::vector<_key_t> keys;
        std::vector<_payload_t> payloads;
        // Merge & sort data in the node
        leaf_to_split->get_data(keys, payloads, plin_->global_version);

        std::vector<Segment> segments;
        auto in_fun = [keys](auto i) {
            return std::pair<_key_t, size_t>(keys[i],i);
        };
        auto out_fun = [&segments](auto cs) { segments.emplace_back(cs); };
        // Train models
        uint64_t last_n = make_segmentation(keys.size(), EPSILON_LEAF_NODE, in_fun, out_fun);
        uint64_t start_pos = 0;
        auto leaf_nodes = new LeafNode*[last_n];
        auto accelerators = new InnerSlot[last_n];
        LeafNode * prev = left_sibling;
        // Build leaf nodes
        for (uint64_t i = 0; i < last_n; ++i) {
            uint64_t block_number = (uint64_t) segments[i].number / (LEAF_NODE_INIT_RATIO * LeafNode::LeafRealSlotsPerBlock) + 3;
            uint64_t node_size_in_byte = block_number * BLOCK_SIZE + NODE_HEADER_SIZE;
            leaf_nodes[i] = new (galc->malloc(node_size_in_byte)) LeafNode(accelerators[i], block_number, &keys[0], &payloads[0], segments[i].number, start_pos, segments[i].slope, segments[i].intercept, plin_->global_version, prev);
            start_pos += segments[i].number;
            prev = leaf_nodes[i];
        }
        for (uint64_t i = 0; i < last_n - 1; ++i){
            leaf_nodes[i]->set_next(leaf_nodes[i + 1]);
        }
        leaf_nodes[last_n - 1]->set_next(right_sibling);

        
        if (accelerator) {
            accelerator->get_read_lock();
        }
        else {
            leaf_to_split->get_read_lock();
            plin_->logs[log_number].set_orphan();
        }
        plin_->logs[log_number].left_node = galc->relative(leaf_nodes[0]);
        plin_->logs[log_number].right_node = galc->relative(leaf_nodes[last_n - 1]);
        plin_->logs[log_number].set_valid();
        do_flush(&plin_->logs[log_number], sizeof(split_log));
        mfence();

        if (left_sibling)
            left_sibling->set_next(leaf_nodes[0]);
        if (right_sibling)
            right_sibling->set_prev(leaf_nodes[last_n - 1]);
        delete [] leaf_nodes;
        leaf_to_split->destroy();
        
        // Insert leaf nodes
        if(accelerator) {
            for (uint64_t i = 0; i < last_n; ++i) {
                upsert_node(accelerators[i]);
            }
        }
        else {
            plin_->orphan_number += last_n - 1;
        }
        delete [] accelerators;

        plin_->release_read_lock();
        if (left_sibling)
            left_sibling->release_lock();
        if (right_sibling)
            right_sibling->release_lock();
        plin_->logs[log_number].release_lock();

        if (double(plin_->orphan_number) / plin_->leaf_number > MAX_ORPHAN_RATIO) {
            #ifdef BACKGROUND_REBUILD
                std::thread rebuild_thread(&SelfType::rebuild_inner_nodes, this);
                rebuild_thread.detach();
            #else
                rebuild_inner_nodes();
            #endif
        }
    }

    void upsert_node (InnerSlot& node) {
        uint32_t i = 0;
        while (i < plin_->root_number && node.min_key >= plin_->roots[i].min_key) {
            ++i;
        }
        --i;
        uint32_t ret = reinterpret_cast<InnerNode*>(galc->absolute(plin_->roots[i].ptr))->upsert_node(node, &plin_->roots[i]);
        if (ret > 1) {
            ++plin_->leaf_number;
        }
        if (ret > 2) {
            ++plin_->orphan_number;
        }
    }

    LeafNode* get_leftmost_leaf() {
        return reinterpret_cast<InnerNode*>(galc->absolute(plin_->roots[0].ptr))->get_leftmost_leaf();
    }

    // Rebuild inner nodes, allow insert and search, no SMO
    void rebuild_inner_nodes() {
        // Get leaf nodes
        if(!plin_->get_write_lock())
            return;

        uint64_t last_n = plin_->leaf_number;
        auto first_keys = new _key_t[last_n];
        auto accelerators = new InnerSlot[last_n];
        LeafNode* node = get_leftmost_leaf();
        
        for (uint64_t i = 0; i < last_n; ++i) {
            node->get_info(first_keys[i], accelerators[i]);
            node = node->get_next();
        }

        // Build inner nodes recursively
        std::vector<Segment> segments;
        uint32_t level = 0;
        uint64_t offset = 0;
        auto accelerators_tmp = accelerators;
        auto first_keys_tmp = first_keys;
        auto in_fun_rec = [&first_keys](auto i){
            return std::pair<_key_t, size_t>(first_keys[i],i);
        };
        auto out_fun = [&segments](auto cs) { segments.emplace_back(cs); };
        while (level == 0 || last_n > ROOT_SIZE) {
            level++;    
            last_n = make_segmentation(last_n, EPSILON_INNER_NODE, in_fun_rec, out_fun);
            // std::cout<<"Number of inner nodes: "<<last_n<<std::endl;
            uint64_t start_pos = 0;
            first_keys = new _key_t[last_n];
            accelerators = new InnerSlot[last_n];
            for (uint64_t i = 0; i < last_n; ++i) {
                uint64_t block_number = (uint64_t) segments[offset + i].number / (INNER_NODE_INIT_RATIO * InnerNode::InnerSlotsPerBlock) + EPSILON_INNER_NODE / InnerNode::InnerSlotsPerBlock + 2;
                uint64_t node_size_in_byte = block_number * BLOCK_SIZE + NODE_HEADER_SIZE;
                first_keys[i] = segments[offset + i].first_key;
                new (galc->malloc(node_size_in_byte)) InnerNode(accelerators[i], block_number, first_keys_tmp, accelerators_tmp, segments[offset + i].number, start_pos, segments[offset + i].slope, segments[offset + i].intercept, level, true);
                start_pos += segments[offset + i].number;
            }
            delete [] accelerators_tmp;
            delete [] first_keys_tmp;
            accelerators_tmp = accelerators;
            first_keys_tmp = first_keys;
            offset += last_n;
        }
        plin_->rebuilding = 1;
        plin_metadata * old_plin_ = (plin_metadata *) galc->malloc(sizeof(plin_metadata));
        plin_->old_plin_ = galc->relative(old_plin_);
        old_plin_->root_number = plin_->root_number;
        for (uint32_t i = 0; i < plin_->root_number; ++i)
            old_plin_->roots[i] = plin_->roots[i];
        old_plin_->level = plin_->level;
        old_plin_->orphan_number = plin_->orphan_number;
        do_flush(old_plin_, sizeof(plin_metadata));
        mfence();

        // Build new root
        plin_->rebuilding = 2;
        plin_->root_number = last_n;
        plin_->level = level;
        plin_->orphan_number = 0;
        for (uint32_t i = 0; i < last_n; ++i) {
            plin_->roots[i] = accelerators_tmp[i];
        }
        for (uint32_t i = last_n; i < ROOT_SIZE; ++i) {
            plin_->roots[i].min_key = FREE_FLAG;
        }
        plin_->release_write_lock();
        for (uint32_t i = 0; i < old_plin_->root_number; i++) {
           galc->absolute((InnerNode*)(old_plin_->roots[i].ptr))->destroy();
        }
        plin_->rebuilding = 0;
        galc->free(old_plin_);
        plin_->old_plin_ = NULL;
        do_flush(plin_, sizeof(plin_metadata));
        mfence();
        delete [] accelerators_tmp;
        delete [] first_keys_tmp;
    }
};

