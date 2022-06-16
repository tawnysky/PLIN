//
// Copyright (c) Zhou Zhang.
//

#pragma once

#include <iostream>
#include <cassert>
#include <cstdint>
#include <climits>
#include <cmath>
#include <vector>
#include <algorithm>
#include <atomic>
// #include <chrono>

// #define DEBUG

// #define BACKGROUND_REBUILD
#define BACKGROUND_SPLIT
// #define BACKGROUND_CHECKLOGS

// How to flush
#define CLWB
// #define CLFLUSHOPT
#define DOFLUSH

using _key_t = double;
using _payload_t = uint64_t;

// Size of data structures
constexpr static uint64_t CACHELINE_SIZE = 64;
constexpr static uint64_t POOLSIZE = 32;
constexpr static uint64_t BLOCK_SIZE = 256;

// Header size
constexpr static uint64_t NODE_HEADER_SIZE = 256;

// Root size
constexpr static uint64_t ROOT_SIZE = 8;
constexpr static uint64_t LOG_NUMBER = 16;

// Epsilon for training models
constexpr static uint64_t EPSILON_LEAF_NODE = 256;
constexpr static uint64_t EPSILON_INNER_NODE = 8;

// Initial filling ratio
constexpr static double INNER_NODE_INIT_RATIO = 0.2;
constexpr static double LEAF_NODE_INIT_RATIO = 0.5;

// Max overflow ratio to split node
constexpr static double MAX_OVERFLOW_RATIO = 0.3;

// Max orphan node ratio to rebuild inner nodes
constexpr static double MAX_ORPHAN_RATIO = 0.1;

// Max buffer size
constexpr static uint64_t MAX_BUFFER = (1 << 20);

constexpr static _key_t FREE_FLAG = (std::numeric_limits<_key_t>::max());
constexpr static _key_t DELETE_FLAG = (std::numeric_limits<_key_t>::max() - 1);

// Concurrency control
constexpr static uint32_t lockSet = ((uint32_t)1 << 31);
constexpr static uint32_t lockMask = ((uint32_t)1 << 31) - 1;
constexpr static uint32_t double_bit_lockSet = ((uint32_t)3 << 30);
constexpr static uint32_t double_bit_lockMask = ((uint32_t)1 << 30) - 1;
constexpr static uint32_t second_bit_lockSet = ((uint32_t)1 << 30);
constexpr static uint32_t second_bit_lockMask = ((uint32_t)3 << 30) - 1;
