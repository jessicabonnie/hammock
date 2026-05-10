#ifndef HAMMOCK_XXHASH_HPP
#define HAMMOCK_XXHASH_HPP

#include <cstdint>
#include <string>

namespace xxhash {
    constexpr uint64_t PRIME64_1 = 0x9E3779B185EBCA87ULL;
    constexpr uint64_t PRIME64_2 = 0xC2B2AE3D27D4EB4FULL;
    constexpr uint64_t PRIME64_3 = 0x165667B19E3779F9ULL;
    constexpr uint64_t PRIME64_4 = 0x85EBCA77C2B2AE63ULL;
    constexpr uint64_t PRIME64_5 = 0x27D4EB2F165667C5ULL;

    inline uint64_t rotl64(uint64_t x, int r) {
        return (x << r) | (x >> (64 - r));
    }

    inline uint64_t avalanche(uint64_t h64) {
        h64 ^= h64 >> 33;
        h64 *= PRIME64_2;
        h64 ^= h64 >> 29;
        h64 *= PRIME64_3;
        h64 ^= h64 >> 32;
        return h64;
    }

    uint64_t hash64(const void* input, size_t len, uint64_t seed = 0);
    inline uint64_t hash64(const std::string& str, uint64_t seed = 0) {
        return hash64(str.data(), str.size(), seed);
    }

    // XXH32 — separate algorithm from xxh64. Required to match Python's
    // xxhash.xxh32(...).intdigest() used for Mode B/C subsampling decisions.
    uint32_t hash32(const void* input, size_t len, uint32_t seed = 0);
    inline uint32_t hash32(const std::string& str, uint32_t seed = 0) {
        return hash32(str.data(), str.size(), seed);
    }
}

#endif
