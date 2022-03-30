/*
    Copyright (c) Luo Yongping. THIS SOFTWARE COMES WITH NO WARRANTIES, 
    USE AT YOUR OWN RISK!
*/


#pragma once

#include <x86intrin.h>

#include "parameters.h"

extern uint32_t flush_cnt;

static inline void mfence() {
    asm volatile("sfence" ::: "memory");
}

static inline void flush(void * ptr) {
#ifdef CLWB
    _mm_clwb(ptr);
#elif defined(CLFLUSHOPT)
    _mm_clflushopt(ptr);
#else
    _mm_clflush(ptr);
#endif
    // flush_cnt++;
}

inline void do_flush (void *data, int len) {
#ifdef DOFLUSH
    volatile char *ptr = (char *)((unsigned long long)data &~(CACHELINE_SIZE-1));
    for(; ptr < (char *)data + len; ptr += CACHELINE_SIZE) {
        flush((void *)ptr);
    }
#endif //DOFLUSH
}

inline void do_flush_with_double_fence (void *data, int len, bool fence=true)
{
#ifdef DOFLUSH
    volatile char *ptr = (char *)((unsigned long long)data &~(CACHELINE_SIZE-1));
    // if(fence) mfence();
    for(; ptr < (char *)data + len; ptr += CACHELINE_SIZE){
        flush((void *)ptr);
    }
    if(fence) mfence();
#endif //DOFLUSH
}

template<typename T>
inline void persist_assign(T* addr, const T &v) { // To ensure atomicity, the size of T should be less equal than 8
    *addr = v;
    do_flush(addr, sizeof(T));
}

