//
// Copyright (c) Zhou Zhang.
//


#pragma once

#include "parameters.h"

#include <omp.h>


// #define DEBUG
#define BACKGROUND_REBUILD

#define CAS(_p, _u, _v)                                                        \
  (__atomic_compare_exchange_n(_p, _u, _v, false, __ATOMIC_ACQUIRE,            \
                               __ATOMIC_ACQUIRE))
// ADD and SUB return the value after add or sub
#define ADD(_p, _v) (__atomic_add_fetch(_p, _v, __ATOMIC_ACQUIRE))
#define SUB(_p, _v) (__atomic_sub_fetch(_p, _v, __ATOMIC_ACQUIRE))
#define LOAD(_p) (__atomic_load_n(_p, __ATOMIC_SEQ_CST))
#define STORE(_p, _v) (__atomic_store_n(_p, _v, __ATOMIC_RELEASE))

#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)


struct InnerSlot {
    _key_t min_key;
    void * ptr;
    double slope;
    float intercept;

    uint32_t info;//info[31]: lock, info[29]: skip, info[28]: type, info[27:0]: block_number
    inline uint32_t block_number() const{
        return info & 0x0fffffffu;
    }
    inline uint32_t type() const{
        return (info & 0x1fffffffu) >>28;
    }
    inline uint32_t skip() const{
        return (info & 0x2fffffffu) >>29;
    }
    inline void set_block_number(uint32_t block_number){
        info = (info & 0xf0000000u) | (block_number & 0x0fffffffu) ;
    }
    inline void set_type(uint32_t type){
        info = (info & (~(1u << 28))) | (type << 28) ;
    }
    inline void set_skip(uint32_t skip){
        info = (info & (~(1u << 29))) | (skip << 29) ;
    }
    inline void set_lock(){
        info = info & (~(1u << 31));
    }

    inline bool get_write_lock() const {
        uint64_t v = __atomic_load_n((uint64_t* ) &ptr, __ATOMIC_ACQUIRE);
        if ((v & lockSet_64)) {
            return false;
        }
        uint64_t old_value = v & lockMask_64;
        uint64_t new_value = old_value | lockSet_64;

        while (!CAS((uint64_t* ) &ptr, &old_value, new_value)) {
            if ((old_value & lockSet_64)) {
                return false;
            }
            old_value = old_value & lockMask_64;
            new_value = old_value | lockSet_64;
        }

        // wait until the readers all exit the critical section
        v = __atomic_load_n((uint64_t* ) &ptr, __ATOMIC_ACQUIRE);
        while (v & lockMask_64) {
            v = __atomic_load_n((uint64_t* ) &ptr, __ATOMIC_ACQUIRE);
        }
        return true;
    }

    inline void get_lock() const {
        uint32_t new_value = 0;
        uint32_t old_value = 0;
        do {
            while (true) {
                old_value = __atomic_load_n(&info, __ATOMIC_ACQUIRE);
                if (!(old_value & lockSet)) {
                    old_value &= lockMask;
                    break;
                }
            }
            new_value = old_value | lockSet;
        } while (!CAS((uint32_t *)&info, &old_value, new_value));
    }
    inline bool try_get_lock() const{
        uint32_t v = __atomic_load_n(&info, __ATOMIC_ACQUIRE);
        if (v & lockSet) {
            return false;
        }
        return true;
    }

    inline void release_lock() const{
        uint32_t v = info;
        __atomic_store_n((uint32_t *)&info, v - lockSet, __ATOMIC_RELEASE);
    }

    inline bool test_lock_set(uint32_t &version) const{
        version = __atomic_load_n(&info, __ATOMIC_ACQUIRE);
        return (version & lockSet) != 0;
    }
};

static void model_correction(double& slope, float& intercept, uint64_t size, _key_t min_key, _key_t max_key) {
    double min_pos = min_key * slope + intercept;
    double max_pos = max_key * slope + intercept + 1;
    if (min_pos > 0 && max_pos < size - 1) {
        slope = (size - 1) / (max_key - min_key);
        intercept = 0 - slope * min_key;
    }
    else if (min_pos > 0) {
        slope = max_pos / (max_key - min_key);
        intercept = 0 - slope * min_key;
    }
    else if (max_pos < size) {
        slope = (size - min_pos - 1) / (max_key - min_key);
        intercept = size - 1 - slope * max_key;
    }
}


inline void fence() { asm volatile("" : : : "memory"); }

inline uint32_t cmpxchgb(uint32_t *object, uint32_t expected,
                               uint32_t desired) {
  asm volatile("lock; cmpxchgb %2,%1"
               : "+a"(expected), "+m"(*object)
               : "r"(desired)
               : "cc");
  fence();
  return expected;
}

