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
        LeafNode* left_sibling = NULL;
        LeafNode* right_sibling = NULL;
        LeafNode* left_node = NULL;
        LeafNode* right_node = NULL;
        bool valid = false;
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
        split_log log;
        plin_metadata* old_plin_ = NULL;
        InnerSlot roots[ROOT_SIZE];

        uint32_t global_version = 0;
        volatile uint32_t smo_lock; //used for rebuilding, no SMO in the process of rebuilding

        inline bool get_write_lock() {
            uint32_t v = __atomic_load_n(&smo_lock, __ATOMIC_ACQUIRE);
            if ((v & lockSet)) {
                return false;
            }
            uint32_t old_value = v & lockMask;
            uint32_t new_value = old_value | lockSet;

            while (!CAS(&smo_lock, &old_value, new_value)) {
                if ((old_value & lockSet)) {
                    return false;
                }
                old_value = old_value & lockMask;
                new_value = old_value | lockSet;
            }

            // wait until the readers all exit the critical section
            v = __atomic_load_n(&smo_lock, __ATOMIC_ACQUIRE);
            while (v & lockMask) {
                v = __atomic_load_n(&smo_lock, __ATOMIC_ACQUIRE);
            }
            return true;
        }
        inline bool try_get_write_lock() {
            uint32_t v = __atomic_load_n(&smo_lock, __ATOMIC_ACQUIRE);
            uint32_t old_value = v & lockMask;
            uint32_t new_value = old_value | lockSet;

            if (!CAS(&smo_lock, &old_value, new_value)) {
                return false;
            }
        }
        inline bool test_write_lock_set() {
            uint32_t v = __atomic_load_n(&smo_lock, __ATOMIC_ACQUIRE);
            return v & lockSet;
        }

        inline bool get_read_lock() {
            uint32_t v = __atomic_load_n(&smo_lock, __ATOMIC_ACQUIRE);
            if ((v & lockSet))
                return false;

            uint32_t old_value = v & lockMask;
            auto new_value = ((v & lockMask) + 1) & lockMask;
            while (!CAS(&smo_lock, &old_value, new_value)) {
                if ((old_value & lockSet)) {
                    return false;
                }
                old_value = old_value & lockMask;
                new_value = ((old_value & lockMask) + 1) & lockMask;
            }

            return true;
        }

        inline bool try_get_read_lock() {
            uint32_t v = __atomic_load_n(&smo_lock, __ATOMIC_ACQUIRE);
            if ((v & lockSet))
                return false;
            return true;
        }

        inline void release_read_lock() { SUB(&smo_lock, 1); }

        inline void reset_rw_lock() {
            smo_lock = 0;
        }

        inline void release_write_lock() {
            __atomic_store_n(&smo_lock, 0, __ATOMIC_RELEASE);
        }
    };

    plin_metadata * plin_;



public:

    PlinIndex(std::string path, std::string id = "plin") {
        galc = new PMAllocator(path.c_str(), false, id.c_str());
        plin_ = (plin_metadata *) galc->get_root(sizeof(plin_metadata));
        plin_->left_buffer = galc->relative(new btree());//memory leak;
        plin_->right_buffer = galc->relative(new btree());//memory leak;
        do_flush(plin_, sizeof(plin_metadata));
        mfence();
    }

    ~PlinIndex() {}

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
            assert(block_number < 1 << 30);
            uint64_t node_size_in_byte = block_number * BLOCK_SIZE + NODE_HEADER_SIZE;
            first_keys[i] = segments[i].first_key;
            leaf_nodes[i] = new (galc->malloc(node_size_in_byte)) LeafNode(accelerators[i], block_number, keys, payloads, segments[i].number, start_pos, segments[i].slope, segments[i].intercept, prev);
            start_pos += segments[i].number;
            prev = leaf_nodes[i];
        }
        for(uint64_t i = 0; i < last_n - 1; ++i){
            leaf_nodes[i]->set_next(leaf_nodes[i + 1]);
        }
           
        plin_->leaf_number = last_n;
        delete [] leaf_nodes;

        // Build inner nodes recursively
        bool is_leaf = true;
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
                assert(block_number < 1 << 30);
                uint64_t node_size_in_byte = block_number * BLOCK_SIZE + NODE_HEADER_SIZE;
                first_keys[i] = segments[offset + i].first_key;
                new (galc->malloc(node_size_in_byte)) InnerNode(accelerators[i], block_number, first_keys_tmp, accelerators_tmp, segments[offset + i].number, start_pos, segments[offset + i].slope, segments[offset + i].intercept, level, is_leaf);
                start_pos += segments[offset + i].number;
            }
            is_leaf = false;
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
        plin_metadata * plin;
        if (key >= plin_->min_key && key <= plin_->max_key) {            
            do
            {
                uint32_t i = 0;
                plin = plin_;
                while (i < plin->root_number && key >= plin->roots[i].min_key) {
                    ++i;
                }
                const InnerSlot * accelerator = &plin->roots[--i];
                if (accelerator->type()) {
                    accelerator = reinterpret_cast<InnerNode*>(galc->absolute(accelerator->ptr))->find_leaf_node(key, accelerator);
                }
                LeafNode * leaf = reinterpret_cast<LeafNode*>(galc->absolute(accelerator->ptr));
                auto ret = leaf->find(key, payload, accelerator);
                if(ret == 1)
                    return true;
                if(ret == 0)
                    return false;
            } while (1);              
        }
        // Find in buffer
        else if (key < plin_->min_key) {
            return galc->absolute(plin_->left_buffer)->update(key, payload);
        }
        else {
            return galc->absolute(plin_->right_buffer)->update(key, payload);
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

    void upsert(_key_t key, _payload_t payload, bool update = false) {
        double t1, t2, t;
        plin_metadata * plin;
        if (key > plin_->min_key && key < plin_->max_key) {
            
            uint32_t ret ;
            
            Retry_Upsert:
                uint32_t i = 0;
                plin = plin_;
                while (i < plin->root_number && key >= plin->roots[i].min_key) {
                    ++i;
                }
                const InnerSlot * accelerator = &plin->roots[--i];
                
                if (accelerator->type()) {
                    accelerator = reinterpret_cast<InnerNode*>(galc->absolute(accelerator->ptr))->find_leaf_node(key, accelerator);
                }
                // ret = 1 : update in a slot; ret = 2 : insert in a free slot; ret = 3 : update in overflow block; ret = 4 : insert in overflow block; 
                // ret = 5 : insert in overflow block & need to split; ret = 6 : insert in overflow block & need to split orphan node
                // ret = 7 : the leaf is in SMO
            
                ret = reinterpret_cast<LeafNode*>(galc->absolute(accelerator->ptr))->upsert(key, payload, accelerator, update);
                // Split leaf node
                switch (ret)
                {
                case 5:
                    split(accelerator, 0);
                    break;
                case 6:
                    split(accelerator, 1);
                    break;
                case 7:// the insert node is in SMO, retry
                    goto Retry_Upsert;
                default:
                    break;
                }
                           
        }
        // Upsert in buffer
        else if (key < plin_->min_key) {
            uint32_t ret = galc->absolute(plin_->left_buffer)->upsert(key, payload);
            if (ret > 3) {
                // TODO: merge buffer
                if (++plin_->left_buffer_number > MAX_BUFFER) {}
            }
        }
        else {
            uint32_t ret = galc->absolute(plin_->right_buffer)->upsert(key, payload);
            if (ret > 3) {
                // TODO: merge buffer
                if (++plin_->right_buffer_number > MAX_BUFFER) {}
            }
        }
    }

    InnerSlot* get_parent(const InnerSlot * node){//
        uint32_t i = 0;
        uint32_t ret;
        while (i < plin_->root_number && node->min_key >= plin_->roots[i].min_key) {
            ++i;
        }
        --i;
        return &plin_->roots[i];
    }

    // Split & insert nodes
    uint32_t split (const InnerSlot * accelerator, uint32_t type) {//type 0: leaf node split, type 1: orphan node split

        if(!plin_->try_get_read_lock())//check wether the index is rebuilding, no smo in rebuilding process
            return 1;
        plin_->get_read_lock();
        
        if(!accelerator->try_get_lock())
            return 1;
        accelerator->get_lock();

        LeafNode* leaf, *prev, *left_sibling, *right_sibling;
        if(!type){
            leaf = reinterpret_cast<LeafNode*>(galc->absolute(accelerator->ptr));
        }
        else{
            leaf = reinterpret_cast<LeafNode*>(galc->absolute(accelerator->ptr))->get_next();
            accelerator->release_lock();
            if(!leaf->try_get_lock())
                return 1;
            leaf->get_lock();
        }
        prev = leaf->get_prev();
        left_sibling = prev;
        right_sibling = leaf->get_next();

        std::vector<_key_t> keys;
        std::vector<_payload_t> payloads;
        // Merge & sort data in the node
        
        _key_t error_key = 0;

        leaf->get_data(keys, payloads, error_key);

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
        // Build leaf nodes
        for (uint64_t i = 0; i < last_n; ++i) {
            uint64_t block_number = (uint64_t) segments[i].number / (LEAF_NODE_INIT_RATIO * LeafNode::LeafRealSlotsPerBlock) + 3;
            assert(block_number < 1 << 30);
            uint64_t node_size_in_byte = block_number * BLOCK_SIZE + NODE_HEADER_SIZE;
            leaf_nodes[i] = new (galc->malloc(node_size_in_byte)) LeafNode(accelerators[i], block_number, &keys[0], &payloads[0], segments[i].number, start_pos, segments[i].slope, segments[i].intercept, prev);
            start_pos += segments[i].number;
            prev = leaf_nodes[i];
        }
        for (uint64_t i = 0; i < last_n - 1; ++i){
            leaf_nodes[i]->set_next(leaf_nodes[i + 1]);
        }
           
        leaf_nodes[last_n - 1]->set_next(right_sibling);

        if (left_sibling)
            plin_->log.left_sibling = galc->relative(left_sibling);
        else
            plin_->log.left_sibling = NULL;
        if (right_sibling)
            plin_->log.right_sibling = galc->relative(right_sibling);
        else
            plin_->log.right_sibling = NULL;
        plin_->log.left_node = galc->relative(leaf_nodes[0]);
        if (last_n > 1)
            plin_->log.right_node = galc->relative(leaf_nodes[last_n - 1]);
        else
            plin_->log.right_node = NULL;
        plin_->log.valid = true;
        do_flush(&plin_->log, sizeof(split_log));
        mfence();

        if (left_sibling)
            left_sibling->set_next(leaf_nodes[0]);
        if (right_sibling)
            right_sibling->set_prev(leaf_nodes[last_n - 1]);
        
        delete [] leaf_nodes;
        // Insert leaf nodes
        if(!type)
            for (uint64_t i = 0; i < last_n; ++i) {
                upsert_node(accelerators[i]);
            }
        else{
            leaf->release_lock();
            plin_->orphan_number  = plin_->orphan_number + last_n -1;
        }
            

        plin_->log.valid = false;
        do_flush(&plin_->log, sizeof(split_log));
        mfence();

        if (double(plin_->orphan_number) / plin_->leaf_number > MAX_ORPHAN_RATIO) {
            #ifdef BACKGROUND_REBUILD
                std::thread rebuild_thread(&SelfType::rebuild_inner_nodes, this);
                rebuild_thread.detach();
            #else
                rebuild_inner_nodes();

            #endif
        }

        plin_->release_read_lock();
        return 1;
    }

    void upsert_node (InnerSlot& node) {
        uint32_t i = 0;
        uint32_t ret;
        do
        {
            while (i < plin_->root_number && node.min_key >= plin_->roots[i].min_key) {
                ++i;
            }
            --i;
            ret = reinterpret_cast<InnerNode*>(galc->absolute(plin_->roots[i].ptr))->upsert_node(node, &plin_->roots[i]);//将叶子节点插入到中间节点
            if (ret > 1) {
                ++plin_->leaf_number;
            }
            if (ret > 2) {
                ++plin_->orphan_number;
            }
        } while (ret==0);
        
        
    }

    LeafNode* get_leftmost_leaf() {
        return reinterpret_cast<InnerNode*>(galc->absolute(plin_->roots[0].ptr))->get_leftmost_leaf();
    }

    // Rebuild inner nodes
    void rebuild_inner_nodes() {//rebuild allow insert and search, no SMO
        // Get leaf nodes
        if(!plin_->try_get_read_lock())
            return;
        plin_->get_write_lock();
        uint64_t last_n = plin_->leaf_number;
        auto first_keys = new _key_t[last_n];
        auto accelerators = new InnerSlot[last_n];
        LeafNode* node = get_leftmost_leaf();
        
        
        for (uint64_t i = 0; i < last_n; ++i) {
            node->get_info(first_keys[i], accelerators[i]);
            node = node->get_next();//get the info of next node
        }

        // Build inner nodes recursively
        std::vector<Segment> segments;
        bool is_leaf = true;
        uint32_t level = 0;
        uint64_t offset = 0;
        auto accelerators_tmp = accelerators;
        auto first_keys_tmp = first_keys;
        auto in_fun_rec = [&first_keys](auto i){
            return std::pair<_key_t, size_t>(first_keys[i],i);
        };
        auto out_fun = [&segments](auto cs) { segments.emplace_back(cs); };
        while (level == 0 || last_n > ROOT_SIZE) {    
            last_n = make_segmentation(last_n, EPSILON_INNER_NODE, in_fun_rec, out_fun);
            // std::cout<<"Number of inner nodes: "<<last_n<<std::endl;
            uint64_t start_pos = 0;
            first_keys = new _key_t[last_n];
            accelerators = new InnerSlot[last_n];
            for (uint64_t i = 0; i < last_n; ++i) {
                uint64_t block_number = (uint64_t) segments[offset + i].number / (INNER_NODE_INIT_RATIO * InnerNode::InnerSlotsPerBlock) + EPSILON_INNER_NODE / InnerNode::InnerSlotsPerBlock + 2;
                assert(block_number < 1 << 30);
                uint64_t node_size_in_byte = block_number * BLOCK_SIZE + NODE_HEADER_SIZE;
                first_keys[i] = segments[offset + i].first_key;
                new (galc->malloc(node_size_in_byte)) InnerNode(accelerators[i], block_number, first_keys_tmp, accelerators_tmp, segments[offset + i].number, start_pos, segments[offset + i].slope, segments[offset + i].intercept, level, is_leaf,true);
                start_pos += segments[offset + i].number;
            }
            is_leaf = false;
            delete [] accelerators_tmp;
            delete [] first_keys_tmp;
            accelerators_tmp = accelerators;
            first_keys_tmp = first_keys;
            offset += last_n;
            level++;
        }

        
        plin_->old_plin_ = (plin_metadata *) galc->malloc(sizeof(plin_metadata));
        plin_->old_plin_->root_number = plin_->root_number;
        for (uint32_t i = 0; i < plin_->root_number; ++i)
            plin_->old_plin_->roots[i] = plin_->roots[i];
        do_flush(plin_->old_plin_, sizeof(plin_metadata));
        mfence();
        // Build new root
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
        
        galc->free(plin_->old_plin_);
        plin_->old_plin_ = NULL;
        do_flush(plin_, sizeof(void*));
        mfence();

        
        delete [] accelerators_tmp;
        delete [] first_keys_tmp;
    }

};

