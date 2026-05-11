#ifndef HAMMOCK_XXHASH_HPP
#define HAMMOCK_XXHASH_HPP

#include <cstdint>
#include <cstring>
#include <string>

namespace xxhash {
    constexpr uint64_t PRIME64_1 = 0x9E3779B185EBCA87ULL;
    constexpr uint64_t PRIME64_2 = 0xC2B2AE3D27D4EB4FULL;
    constexpr uint64_t PRIME64_3 = 0x165667B19E3779F9ULL;
    constexpr uint64_t PRIME64_4 = 0x85EBCA77C2B2AE63ULL;
    constexpr uint64_t PRIME64_5 = 0x27D4EB2F165667C5ULL;

    constexpr uint32_t PRIME32_1 = 0x9E3779B1U;
    constexpr uint32_t PRIME32_2 = 0x85EBCA77U;
    constexpr uint32_t PRIME32_3 = 0xC2B2AE3DU;
    constexpr uint32_t PRIME32_4 = 0x27D4EB2FU;
    constexpr uint32_t PRIME32_5 = 0x165667B1U;

    inline uint64_t rotl64(uint64_t x, int r) {
        return (x << r) | (x >> (64 - r));
    }

    inline uint32_t rotl32(uint32_t x, int r) {
        return (x << r) | (x >> (32 - r));
    }

    inline uint64_t avalanche(uint64_t h64) {
        h64 ^= h64 >> 33;
        h64 *= PRIME64_2;
        h64 ^= h64 >> 29;
        h64 *= PRIME64_3;
        h64 ^= h64 >> 32;
        return h64;
    }

    // Bodies live in the header so callers (especially the Mode B/C stride
    // inner loop) can inline them. For short inputs (<32 bytes) the compiler
    // can also DCE the long-input branch when the call site's length is
    // statically bounded.
    inline uint64_t hash64(const void* input, size_t len, uint64_t seed = 0) {
        const uint8_t* p = (const uint8_t*)input;
        const uint8_t* const end = p + len;
        uint64_t h64;

        if (len >= 32) {
            const uint8_t* const limit = end - 32;
            uint64_t v1 = seed + PRIME64_1 + PRIME64_2;
            uint64_t v2 = seed + PRIME64_2;
            uint64_t v3 = seed + 0;
            uint64_t v4 = seed - PRIME64_1;

            do {
                uint64_t k1, k2, k3, k4;
                std::memcpy(&k1, p, 8); p += 8;
                std::memcpy(&k2, p, 8); p += 8;
                std::memcpy(&k3, p, 8); p += 8;
                std::memcpy(&k4, p, 8); p += 8;

                v1 = rotl64(v1 + k1 * PRIME64_2, 31) * PRIME64_1;
                v2 = rotl64(v2 + k2 * PRIME64_2, 31) * PRIME64_1;
                v3 = rotl64(v3 + k3 * PRIME64_2, 31) * PRIME64_1;
                v4 = rotl64(v4 + k4 * PRIME64_2, 31) * PRIME64_1;
            } while (p <= limit);

            h64 = rotl64(v1, 1) + rotl64(v2, 7) + rotl64(v3, 12) + rotl64(v4, 18);

            v1 *= PRIME64_2; v1 = rotl64(v1, 31); v1 *= PRIME64_1;
            h64 ^= v1; h64 = h64 * PRIME64_1 + PRIME64_4;

            v2 *= PRIME64_2; v2 = rotl64(v2, 31); v2 *= PRIME64_1;
            h64 ^= v2; h64 = h64 * PRIME64_1 + PRIME64_4;

            v3 *= PRIME64_2; v3 = rotl64(v3, 31); v3 *= PRIME64_1;
            h64 ^= v3; h64 = h64 * PRIME64_1 + PRIME64_4;

            v4 *= PRIME64_2; v4 = rotl64(v4, 31); v4 *= PRIME64_1;
            h64 ^= v4; h64 = h64 * PRIME64_1 + PRIME64_4;
        } else {
            h64 = seed + PRIME64_5;
        }

        h64 += (uint64_t)len;

        while (p + 8 <= end) {
            uint64_t k1;
            std::memcpy(&k1, p, 8);
            k1 *= PRIME64_2;
            k1 = rotl64(k1, 31);
            k1 *= PRIME64_1;
            h64 ^= k1;
            h64 = rotl64(h64, 27) * PRIME64_1 + PRIME64_4;
            p += 8;
        }

        if (p + 4 <= end) {
            uint32_t k1;
            std::memcpy(&k1, p, 4);
            h64 ^= (uint64_t)k1 * PRIME64_1;
            h64 = rotl64(h64, 23) * PRIME64_2 + PRIME64_3;
            p += 4;
        }

        while (p < end) {
            h64 ^= (*p++) * PRIME64_5;
            h64 = rotl64(h64, 11) * PRIME64_1;
        }

        return avalanche(h64);
    }

    inline uint64_t hash64(const std::string& str, uint64_t seed = 0) {
        return hash64(str.data(), str.size(), seed);
    }

    // Short-input fast path. Caller knows the input is small (<32 bytes) so
    // we drop the wide-lane bulk-mixing loop entirely. The Mode B/C stride
    // loop calls this for every genomic point — input is chr+sep+pos and
    // never gets near 32 bytes in practice. If a caller hits the rare long
    // case, we fall back to the general hash64 (branch predicts well).
    inline uint64_t hash64_short(const void* input, size_t len, uint64_t seed = 0) {
        if (__builtin_expect(len >= 32, 0)) {
            return hash64(input, len, seed);
        }
        const uint8_t* p = (const uint8_t*)input;
        const uint8_t* const end = p + len;
        uint64_t h64 = seed + PRIME64_5 + (uint64_t)len;

        while (p + 8 <= end) {
            uint64_t k1;
            std::memcpy(&k1, p, 8);
            k1 *= PRIME64_2;
            k1 = rotl64(k1, 31);
            k1 *= PRIME64_1;
            h64 ^= k1;
            h64 = rotl64(h64, 27) * PRIME64_1 + PRIME64_4;
            p += 8;
        }

        if (p + 4 <= end) {
            uint32_t k1;
            std::memcpy(&k1, p, 4);
            h64 ^= (uint64_t)k1 * PRIME64_1;
            h64 = rotl64(h64, 23) * PRIME64_2 + PRIME64_3;
            p += 4;
        }

        while (p < end) {
            h64 ^= (*p++) * PRIME64_5;
            h64 = rotl64(h64, 11) * PRIME64_1;
        }

        return avalanche(h64);
    }

    // XXH32 — separate algorithm from xxh64. Required to match Python's
    // xxhash.xxh32(...).intdigest() used for Mode B/C subsampling decisions.
    inline uint32_t hash32(const void* input, size_t len, uint32_t seed = 0) {
        const uint8_t* p = (const uint8_t*)input;
        const uint8_t* const end = p + len;
        uint32_t h32;

        if (len >= 16) {
            const uint8_t* const limit = end - 16;
            uint32_t v1 = seed + PRIME32_1 + PRIME32_2;
            uint32_t v2 = seed + PRIME32_2;
            uint32_t v3 = seed + 0;
            uint32_t v4 = seed - PRIME32_1;

            do {
                uint32_t k1, k2, k3, k4;
                std::memcpy(&k1, p, 4); p += 4;
                std::memcpy(&k2, p, 4); p += 4;
                std::memcpy(&k3, p, 4); p += 4;
                std::memcpy(&k4, p, 4); p += 4;
                v1 = rotl32(v1 + k1 * PRIME32_2, 13) * PRIME32_1;
                v2 = rotl32(v2 + k2 * PRIME32_2, 13) * PRIME32_1;
                v3 = rotl32(v3 + k3 * PRIME32_2, 13) * PRIME32_1;
                v4 = rotl32(v4 + k4 * PRIME32_2, 13) * PRIME32_1;
            } while (p <= limit);

            h32 = rotl32(v1, 1) + rotl32(v2, 7) + rotl32(v3, 12) + rotl32(v4, 18);
        } else {
            h32 = seed + PRIME32_5;
        }

        h32 += static_cast<uint32_t>(len);

        while (p + 4 <= end) {
            uint32_t k1;
            std::memcpy(&k1, p, 4);
            h32 += k1 * PRIME32_3;
            h32 = rotl32(h32, 17) * PRIME32_4;
            p += 4;
        }

        while (p < end) {
            h32 += (*p++) * PRIME32_5;
            h32 = rotl32(h32, 11) * PRIME32_1;
        }

        h32 ^= h32 >> 15;
        h32 *= PRIME32_2;
        h32 ^= h32 >> 13;
        h32 *= PRIME32_3;
        h32 ^= h32 >> 16;
        return h32;
    }

    inline uint32_t hash32(const std::string& str, uint32_t seed = 0) {
        return hash32(str.data(), str.size(), seed);
    }
}

#endif
