/*  
    A wrapper for using PMDK allocator easily (refer to pmwcas)
    Copyright (c) Luo Yongping. THIS SOFTWARE COMES WITH NO WARRANTIES, 
    USE AT YOUR OWN RISK!
*/

#pragma once

#include <sys/stat.h>
#include <cstdio>
#include <omp.h>

#include <libpmemobj.h>
#include <libpmem.h>

#include "flush.h"
#include "spinlock.h"

POBJ_LAYOUT_BEGIN(pmallocator);
POBJ_LAYOUT_TOID(pmallocator, char)
POBJ_LAYOUT_END(pmallocator)

/*
    Persistent Memory Allocator: a wrapper of PMDK allocation lib: https://pmem.io/pmdk/

    PMAllocator allocates persistent memory from a pool file that resides on NVM file system. 
    It uses malloc() and free() as the allocation and reclaiment interfaces. 
    Other public interfaces like get_root(), absolute() and relative() are essential to memory
    management in persistent environment. 
*/


class PMAllocator {
private:
    static const int PIECE_CNT = 64;
    static const size_t ALIGN_SIZE = 256;
    static const int64_t DEFAULT_POOL_SIZE = POOLSIZE * (1024UL * 1024 * 1024); // set the default pool size to be 10GB 
    
    struct MetaType {
        char * buffer[PIECE_CNT];
        size_t blk_per_piece;
        size_t cur_blk;
        // entrance of DS in buffer
        void * entrance;
    };
    MetaType * meta_;

    // volatile domain
    PMEMobjpool *pop_;

    char * buff_[PIECE_CNT];
    char * buff_aligned_[PIECE_CNT];
    size_t piece_size_;
    size_t cur_blk_;
    size_t max_blk_;
    Spinlock alloc_mtx;

    bool file_exist (const char *pool_path) {
        struct stat buffer;
        return (stat(pool_path, &buffer) == 0);
    }

public: 
    /*
     *  Construct a PM allocator, map a pool file into virtual memory
     *  @param filename     pool file name
     *  @param recover      if doing recover, false for the first time 
     *  @param layout_name  ID of a group of allocations (in characters), each ID corresponding to a root entry
     *  @param pool_size    pool size of the pool file, vaild if the file doesn't exist
     */
    PMAllocator(const char *file_name, bool recover, const char *layout_name, int64_t pool_size = DEFAULT_POOL_SIZE) {
        PMEMobjpool *tmp_pool = nullptr;
        pool_size = pool_size + ((pool_size & ((1 << 23) - 1)) > 0 ? (1 << 23) : 0); // align to 8MB
	    if(recover == false) {
            if(file_exist(file_name)) {
                printf("[CAUTIOUS]: The pool file already exists\n");
                printf("Try (1) remove the pool file %s\nOr  (2) set the recover parameter to be true\n", file_name);
                exit(-1);
            }
            pop_ = pmemobj_create(file_name, layout_name, pool_size, S_IWUSR | S_IRUSR);
            meta_ = (MetaType *)pmemobj_direct(pmemobj_root(pop_, sizeof(MetaType)));
            
            // maintain volatile domain
            uint64_t alloc_size = pool_size >> 3; // 1/8 of the pool is used as block alloction
            for(int i = 0; i < PIECE_CNT; i++) {
                buff_[i] = (char *)mem_alloc(alloc_size / PIECE_CNT);
                buff_aligned_[i] = (char *) ((uint64_t)buff_[i] + ((uint64_t) buff_[i] % ALIGN_SIZE == 0 ? 0 : (ALIGN_SIZE - (uint64_t) buff_[i] % ALIGN_SIZE)));
            }
            piece_size_ = (alloc_size / PIECE_CNT) / ALIGN_SIZE - 1;
            cur_blk_ = 0;
            max_blk_ = piece_size_ * PIECE_CNT;
            
            // initialize meta_
            for(int i = 0; i < PIECE_CNT; i++) 
                meta_->buffer[i] = relative(buff_[i]);
            meta_->blk_per_piece = piece_size_;
            meta_->cur_blk = 0;
            meta_->entrance = NULL;
            do_flush(meta_, sizeof(MetaType));
        } else {
            if(!file_exist(file_name)) {
                printf("Pool File Not Exist\n");
		        exit(-1);
	        }
            pop_ = pmemobj_open(file_name, layout_name);
            meta_ = (MetaType *)pmemobj_direct(pmemobj_root(pop_, sizeof(MetaType)));
            // maintain volatile domain
            for(int i = 0; i < PIECE_CNT; i++) {
                buff_[i] = absolute(meta_->buffer[i]);
                buff_aligned_[i] = (char *) ((uint64_t)buff_[i] + ((uint64_t) buff_[i] % ALIGN_SIZE == 0 ? 0 : (ALIGN_SIZE - (uint64_t) buff_[i] % ALIGN_SIZE)));
            }
            cur_blk_ = meta_->cur_blk;
            piece_size_ = meta_->blk_per_piece;
            max_blk_ = piece_size_ * PIECE_CNT;
        }
    }

    ~PMAllocator() {
        pmemobj_close(pop_);
    }

public:
    /*  
     *  Get/allocate the root entry of the allocator.
     *  
     *  The root entry is the entrance of one group of allocation, each group is
     *  identified by the layout_name when constructing it.
     * 
     *  Each group of allocations is a independent, self-contained in-memory structure in the pool
     *  such as b-tree or link-list
     */
    void * get_root(size_t nsize) { // the root of DS stored in buff_ is recorded at meta_->entrance
        if(meta_->entrance == NULL) {
            meta_->entrance = relative(malloc(nsize));
            do_flush(meta_, sizeof(MetaType));
        }
        return absolute(meta_->entrance);
    }

    /*
     *  Allocate a non-root piece of persistent memory from the mapped pool
     *  return the virtual memory address
     */
    void * malloc(size_t nsize) { 
        if(nsize >= (1 << 12)) { // large than 4KB
            void * mem = mem_alloc(nsize + ALIGN_SIZE); // not aligned
            //  |  UNUSED    |HEADER|       memory you can use     |
            // mem             (mem + off)
            uint64_t offset = ALIGN_SIZE - (uint64_t)mem % ALIGN_SIZE;
            // store a header in the front
            uint64_t * header = (uint64_t *)((uint64_t)mem + offset - 8);
            *header = offset;

            return (void *)((uint64_t)mem + offset);
        }
        

        int blk_demand = (nsize + ALIGN_SIZE - 1) / ALIGN_SIZE;
        // case 1: not enough in the buffer
        if(blk_demand + cur_blk_ > max_blk_) {
            printf("run out of memory\n");
            exit(-1);
        }
        // case 2: current piece can not accommdate this allocation
        alloc_mtx.lock();
        int piece_id = cur_blk_ / piece_size_;
        if((cur_blk_ % piece_size_ + blk_demand) > piece_size_) {
            void * mem = buff_aligned_[piece_id + 1]; // allocate from a new peice
            cur_blk_ = piece_size_ * (piece_id + 1) + blk_demand;
            meta_->cur_blk = cur_blk_;
            do_flush(&(meta_->cur_blk), 8);
        
            alloc_mtx.unlock();
            return mem;
        } 
        // case 3: current piece has enough space
        else {
            void * mem = buff_aligned_[piece_id] + ALIGN_SIZE * (cur_blk_ % piece_size_);

            cur_blk_ = cur_blk_ + blk_demand;
            meta_->cur_blk = cur_blk_;
            do_flush(&(meta_->cur_blk), 8);       
            alloc_mtx.unlock();
            return mem;
        }
    }

    void free(void* addr) {
        for(int i = 0; i < PIECE_CNT; i++) {
            uint64_t offset = (uint64_t)addr - (uint64_t)buff_aligned_[i];
            if(offset > 0 && offset < piece_size_ * ALIGN_SIZE) {
                // the addr is in this piece, do not reclaim it
                return ;
            }
        }

        // larger than 4KB, reclaim it 
        uint64_t * header = (uint64_t *)((uint64_t)addr - 8);
        uint64_t offset = *header; 

        alloc_mtx.lock();
        auto oid_ptr = pmemobj_oid((void *)((uint64_t)addr - offset));
        alloc_mtx.unlock();

        TOID(char) ptr_cpy;
        TOID_ASSIGN(ptr_cpy, oid_ptr);
        POBJ_FREE(&ptr_cpy);
        
    }  

    /*
     *  Distinguish from virtual memory address and offset in the pool
     *  Each memory piece allocated from the pool has an in-pool offset, which remains unchanged
     *  until reclaiment. We cannot ensure that the pool file is mapped at the same position at 
     *  any time, so it may locate at different virtual memory addresses next time. 
     *  
     *  So the rule is that, using virtual memory when doing normal operations like to DRAM
     *  space, using offset to store link relationship, for exmaple, next pointer in linklist
     */

    /*
     *  convert an offset to a virtual memory address
     */
    template<typename T>
    inline T *absolute(T *pmem_offset) {
        if(pmem_offset == NULL)
            return NULL;
        return reinterpret_cast<T *>(reinterpret_cast<uint64_t>(pmem_offset) + reinterpret_cast<char *>(pop_));
    }
    
    /*
     *  convert a virtual memory address to an offset
     */
    template<typename T>
    inline T *relative(T *pmem_direct) {
        if(pmem_direct == NULL)
            return NULL;
        return reinterpret_cast<T *>(reinterpret_cast<char *>(pmem_direct) - reinterpret_cast<char *>(pop_));
    }

private:
    void * mem_alloc(size_t nsize) {
        PMEMoid tmp;
        
        alloc_mtx.lock();
        pmemobj_alloc(pop_, &tmp, nsize, TOID_TYPE_NUM(char), NULL, NULL);
        
        
        void * mem = pmemobj_direct(tmp);
        alloc_mtx.unlock();
        assert(mem != nullptr);
        return mem;
    }
};

extern PMAllocator * galc;

