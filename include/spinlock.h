/* Copyright(c): Guo Zhongming 
*/

#ifndef __SPINLOCK_H__
#define __SPINLOCK_H__

#include <atomic>
#include <emmintrin.h>
#include <thread>

class Spinlock {
public:
    Spinlock() {
        atomic_val.store(0, std::memory_order_relaxed);
    }

    Spinlock(const Spinlock &) = delete;
    Spinlock & operator = (const Spinlock &) = delete;

public:
    inline void lock() {
        while(atomic_val.exchange(1, std::memory_order_acquire) == 1){
            while(1) {
                _mm_pause(); // delay for 140 cycle

                if(atomic_val.load(std::memory_order_relaxed) == 0) // check the atomic_val
                    break;
                
                std::this_thread::yield(); // delay for 113ns

                if(atomic_val.load(std::memory_order_relaxed) == 0) // check the atomic_val
                    break;
            }

            // if at here, the atomic_val must be just be 0
        }

        return ;
    }

    inline void unlock() {
        atomic_val.exchange(0, std::memory_order_acquire);
    }

    inline bool trylock() {
        return atomic_val.exchange(1, std::memory_order_acquire) == 0;
    }

private:
    std::atomic_short atomic_val;
}; 

#endif // __SPINLOCK_H__