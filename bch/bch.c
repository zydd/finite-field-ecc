#include <stdint.h>
#include <string.h>
#include <stdio.h>

// static const uint32_t gen = 0x2;                // 𝑥
// static const uint32_t irr_poly = 0x43;          // 𝑥⁶ + 𝑥 + 1
// static const uint32_t field_charac = (1 << 6);  // 2 ** 6

// (𝑥¹⁸ + 𝑥¹⁷ + 𝑥¹⁶) + (𝑥¹⁵) + (𝑥⁹) + (𝑥⁷ + 𝑥⁶) + (𝑥³ + 𝑥² + 𝑥 + 1)
static const uint32_t generator18 = 0x782cf;
// (𝑥³³ + 𝑥³²) + (𝑥³⁰ + 𝑥²⁹ + 𝑥²⁸) + (𝑥²⁷ + 𝑥²⁶) + (𝑥²³ + 𝑥²² + 𝑥²⁰)
// + (𝑥¹⁵ + 𝑥¹⁴ + 𝑥¹³) + (𝑥¹¹ + 𝑥⁹ + 𝑥⁸) + (𝑥⁶ + 𝑥⁵) + (𝑥² + 𝑥 + 1)
static const uint64_t generator33 = 0x37cd0eb67;

static const uint8_t exp[] = {
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x03, 0x06,
    0x0c, 0x18, 0x30, 0x23, 0x05, 0x0a, 0x14, 0x28,
    0x13, 0x26, 0x0f, 0x1e, 0x3c, 0x3b, 0x35, 0x29,
    0x11, 0x22, 0x07, 0x0e, 0x1c, 0x38, 0x33, 0x25,
    0x09, 0x12, 0x24, 0x0b, 0x16, 0x2c, 0x1b, 0x36,
    0x2f, 0x1d, 0x3a, 0x37, 0x2d, 0x19, 0x32, 0x27,
    0x0d, 0x1a, 0x34, 0x2b, 0x15, 0x2a, 0x17, 0x2e,
    0x1f, 0x3e, 0x3f, 0x3d, 0x39, 0x31, 0x21, 0x01,
};

static const uint8_t log[] = {
    0x00, 0x3f, 0x01, 0x06, 0x02, 0x0c, 0x07, 0x1a,
    0x03, 0x20, 0x0d, 0x23, 0x08, 0x30, 0x1b, 0x12,
    0x04, 0x18, 0x21, 0x10, 0x0e, 0x34, 0x24, 0x36,
    0x09, 0x2d, 0x31, 0x26, 0x1c, 0x29, 0x13, 0x38,
    0x05, 0x3e, 0x19, 0x0b, 0x22, 0x1f, 0x11, 0x2f,
    0x0f, 0x17, 0x35, 0x33, 0x25, 0x2c, 0x37, 0x28,
    0x0a, 0x3d, 0x2e, 0x1e, 0x32, 0x16, 0x27, 0x2b,
    0x1d, 0x3c, 0x2a, 0x15, 0x14, 0x3b, 0x39, 0x3a,
};

static inline uint8_t gf_mul(uint8_t a, uint8_t b) {
    if (a == 0 || b == 0)
        return 0;

    unsigned r = log[a] + log[b];
    if (r >= 63)
        r -= 63;

    return exp[r];
}

static inline uint8_t gf_div(uint8_t a, uint8_t b) {
    unsigned r = log[a] + 63 - log[b];
    if (r >= 63)
        r -= 63;

    return exp[r];
}

static inline uint8_t gf_inv(uint8_t a) {
    return exp[63 - log[1]];
}

uint8_t gf_poly_eval(const uint8_t poly[], const unsigned size, const uint8_t x0) {
    uint8_t r = 0;
    uint8_t x = 1;
    for (unsigned i = 0; i < size; ++i) {
        r ^= gf_mul(poly[i], x);
        x = gf_mul(x, x0);
    }
    return r;
}

uint8_t gf2_poly63_eval(const uint8_t poly[8], const uint8_t x0) {
    uint8_t r = 0;
    uint8_t x = 1;
    for (unsigned i = 1; i < 64; ++i) {
        unsigned byte = i / 8;
        uint8_t coef = (poly[7 - byte] >> (i - byte * 8)) & 1;

        if (coef)
            r ^= x;
        x = gf_mul(x, x0);
    }
    return r;
}

void gf_ex_synth_div(uint8_t a[], unsigned size_a, const uint8_t b[], unsigned size_b) {
    uint8_t normalizer = b[size_b-1];

    for (unsigned i = size_a -1; i > size_b-1 -1; --i) {
        a[i] = gf_div(a[i], normalizer);

        uint8_t coef = a[i];
        if (coef != 0)
            for (int j = 1; j < size_b; j++)
                a[i - j] ^= gf_mul(b[size_b-1 - j], coef);
    }
}

int syndromes(const uint8_t data[8], uint8_t synds[], const unsigned size) {
    int err = 0;

    for (unsigned i = 0; i < size; ++i) {
        uint8_t a = exp[i + 1];
        err |= synds[i] = gf2_poly63_eval(data, a);
    }

    return err;
}

static inline void gf_poly_shift(uint8_t poly[], const unsigned size, const unsigned n) {
    unsigned i = size;
    for (; i > n; --i)
        poly[i-1] = poly[i-1 - n];
    for (; i > 0; --i)
        poly[i-1] = 0;
}

static inline void gf_poly_scale(uint8_t poly[], const unsigned size, const uint8_t a) {
    for (unsigned i = 0; i < size; ++i)
        poly[i] = gf_mul(a, poly[i]);
}

static inline void gf_poly_sub(uint8_t poly_a[], const uint8_t poly_b[], const unsigned size) {
    for (unsigned i = 0; i < size; ++i)
        poly_a[i] ^= poly_b[i];
}

void gf_poly_mul(uint8_t r[],
                 const uint8_t poly_a[], const unsigned size_a,
                 const uint8_t poly_b[], const unsigned size_b) {
    for (int i = 0; i < size_a + size_b; ++i)
        r[i] = 0;

    for (int i = 0; i < size_a; ++i)
        for (int j = 0; j < size_b; ++j)
            r[i + j] ^= gf_mul(poly_a[i], poly_b[j]);
}

static inline void gf_poly_deriv(uint8_t poly[], const unsigned size) {
    for (int i = 0; i < size - 1; ++i)
        poly[i] = (i & 1) ? 0 : poly[i + 1];

    poly[size-1] = 0;
}

unsigned gf_poly_mod(uint8_t a[], const unsigned size_a,
                     const uint8_t b[], unsigned size_b) {
    while (size_b > 0 && b[size_b - 1] == 0)
        --size_b;
    gf_ex_synth_div(a, 6, b, size_b);
    return size_a - size_b;
}

unsigned berlekamp_massey(const uint8_t synds[], uint8_t C[], const unsigned size) {
    // uint8_t C[size]; // correction polynomial
    uint8_t B[12]; // previous C
    uint8_t T[12];

    for (int i = 0; i < size; ++i) {
        C[i] = i == 0;
        B[i] = i == 0;
    }

    unsigned L = 0;
    unsigned m = 1;
    uint8_t b = 1;

    for (int n = 0; n < size; ++n) {
        // discrepancy
        unsigned d = synds[n];
        for (int i = 1; i < L+1; ++i)
            d ^= gf_mul(C[i], synds[n-i]);

        if (d == 0) {
            // discrepancy is zero; annihilation continues
            m = m + 1;
        } else if (2 * L <= n) {
            memcpy(T, C, size); // temporary copy of C

            gf_poly_shift(B, size, m);
            gf_poly_scale(B, size, gf_div(d, b));

            gf_poly_sub(C, B, size);

            L = n + 1 - L;
            memcpy(B, T, size);

            b = d;
            m = 1;
        } else {
            memcpy(T, B, size);

            gf_poly_shift(T, size, m);

            gf_poly_scale(T, size, gf_div(d, b));
            gf_poly_sub(C, T, size);

            m = m + 1;
        }
    }

    return L;
}


unsigned find_err_pos(const uint8_t poly[], const unsigned poly_size, uint8_t err_pos[], const unsigned err_size) {
    unsigned count = 0;

    for (int i = 0; i < 63; ++i) {
        uint8_t xinv = gf_div(1, exp[i]);
        if (gf_poly_eval(poly, poly_size, xinv) == 0)
            err_pos[count++] = i;
    }

    return count;
}

void encode63_45(uint8_t data[8]) {
    const uint32_t poly = generator18 << (32 - 18);
    uint32_t rem = 0;

    for (unsigned i = 0; i < 5; ++i) {
        rem ^= (uint32_t) data[i] << 24;

        for (int j = 0; j < 8; ++j) {
            if (rem & 0x80000000)
                rem = (rem << 1) ^ poly;
            else
                rem = (rem << 1);
        }
    }

    {
        data[5] &= 0xf8;
        rem ^= data[5] << 24;

        for (int i = 0; i < 5; i++) {
            if (rem & 0x80000000)
                rem = (rem << 1) ^ poly;
            else
                rem = (rem << 1);
        }
    }
    rem >>= 32 - 18;

    data[5] |= rem >> 15;
    data[6] = rem >> 7;
    data[7] = (rem << 1) | 1;
}

void encode63_30(uint8_t data[8]) {
    const uint64_t poly = generator33 << (64 - 33);
    uint64_t rem = 0;

    for (unsigned i = 0; i < 3; ++i) {
        rem ^= (uint64_t) data[i] << 56;

        for (int j = 0; j < 8; ++j) {
            if (rem & 0x8000000000000000)
                rem = (rem << 1) ^ poly;
            else
                rem = (rem << 1);
        }
    }

    {
        data[3] &= 0xfc;
        rem ^= (uint64_t) data[3] << 56;

        for (int i = 0; i < 6; i++) {
            if (rem & 0x8000000000000000)
                rem = (rem << 1ull) ^ poly;
            else
                rem = (rem << 1ull);
        }
    }
    rem >>= 64 - 33;

    data[3] |= rem >> 31;
    data[4] = rem >> 23;
    data[5] = rem >> 15;
    data[6] = rem >> 7;
    data[7] = (rem << 1) | 1;
}

int check63_45(const uint8_t data[8]) {
    const uint32_t poly = generator18 << (32 - 18);
    uint32_t rem = 0;

    for (unsigned i = 0; i < 7; ++i) {
        rem ^= data[i] << 24;

        for (int j = 0; j < 8; ++j) {
            if (rem & 0x80000000)
                rem = (rem << 1) ^ poly;
            else
                rem = (rem << 1);
        }
    }

    {
        rem ^= (data[7] & 0xfe) << 24;

        for (int i = 0; i < 8; i++) {
            if (rem & 0x80000000)
                rem = (rem << 1) ^ poly;
            else
                rem = (rem << 1);
        }
    }
    return (rem >> (32 - 18)) == 0;
}


int decode63_45(uint8_t data[8]) {
    uint8_t synds[6];
    uint8_t err_poly[6];

    if (!syndromes(data, synds, 6))
        return 1;

    unsigned error_count = berlekamp_massey(synds, err_poly, 6);

    uint8_t err_pos[3];
    if (error_count == 0 || error_count != find_err_pos(err_poly, 6, err_pos, 3))
        return 0;

    for (unsigned i = 0; i < error_count; ++i) {
        unsigned bit = err_pos[i] + 1;
        unsigned byte = bit / 8;
        data[7 - byte] ^= 1 << (bit - byte * 8);
    }

    return 1;
}

int decode63_30(uint8_t data[8]) {
    uint8_t synds[12];
    uint8_t err_poly[12];

    if (!syndromes(data, synds, 12))
        return 1;

    unsigned error_count = berlekamp_massey(synds, err_poly, 12);

    uint8_t err_pos[6];
    if (error_count == 0 || error_count != find_err_pos(err_poly, 12, err_pos, 6))
        return 0;

    for (unsigned i = 0; i < error_count; ++i) {
        unsigned bit = err_pos[i] + 1;
        unsigned byte = bit / 8;
        data[7 - byte] ^= 1 << (bit - byte * 8);
    }

    return 1;
}
