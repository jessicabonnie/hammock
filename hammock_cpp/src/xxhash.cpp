#include "../include/xxhash.h"
#include <cstring>

namespace xxhash {

uint64_t hash64(const void* input, size_t len, uint64_t seed) {
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
            memcpy(&k1, p, 8); p += 8;
            memcpy(&k2, p, 8); p += 8;
            memcpy(&k3, p, 8); p += 8;
            memcpy(&k4, p, 8); p += 8;

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
        memcpy(&k1, p, 8);
        k1 *= PRIME64_2;
        k1 = rotl64(k1, 31);
        k1 *= PRIME64_1;
        h64 ^= k1;
        h64 = rotl64(h64, 27) * PRIME64_1 + PRIME64_4;
        p += 8;
    }

    if (p + 4 <= end) {
        uint32_t k1;
        memcpy(&k1, p, 4);
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

}  // namespace xxhash
