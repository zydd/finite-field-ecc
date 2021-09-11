#pragma once

#include <cassert>
#include <cstdint>

template<uint8_t primitive, uint16_t poly1>
struct GF {
    static_assert((poly1 & 0x80) == 0); // needed for mul4 and mul8
    static const uint32_t poly4 = (poly1 & 0xff) * 0x01010101;
    static const uint64_t poly8 = (poly1 & 0xff) * 0x0101010101010101;

    std::array<uint8_t, 256> exp_table;
    std::array<uint8_t, 256> log_table;
#ifdef GF_MUL_TABLE
    std::array<std::array<uint8_t, 256>, 256> mul_table;
#endif

    inline GF() {
        uint8_t x = 1;
        for (unsigned i = 0; i < 256; ++i) {
            exp_table[i] = x;
            log_table[x] = i;
            x = _mul(x, primitive);
        }

#ifdef GF_MUL_TABLE
        for (unsigned i = 0; i < 256; ++i) {
            for (unsigned j = 0; j < 256; ++j)
                mul_table[i][j] = _mul(i, j);
        }
#endif
    }

    inline uint8_t _mul(uint8_t a, uint8_t b) {
        uint8_t r = 0;
        for (int i = 7; i >= 0; --i) {
            if (r & 0x80)
                r = (r << 1) ^ poly1;
            else
                r = (r << 1);

            if (a & (1 << i))
                r ^= b;
        }

        return r;
    }

#ifdef GF_MUL_TABLE
    inline uint8_t mul(uint8_t a, uint8_t b) {
        return mul_table[a][b];
    }
#else
    inline uint8_t mul(uint8_t a, uint8_t b) {
        if (a == 0 || b == 0)
            return 0;

        unsigned r = log_table[a] + log_table[b];
        if (r >= 255)
            r -= 255;

        return exp_table[r];
    }
#endif

    inline uint32_t mul(uint32_t a, uint32_t b) {
        uint32_t r = 0;
        for (int i = 7; i >= 0; --i) {
            uint32_t m = r & 0x80808080;
            m = m - (m >> 7);

            r = ((r & 0x7f7f7f7f) << 1) ^ (poly4 & m);

            uint32_t n = (a & (0x01010101 << i)) >> i;
            n = (n << 8) - n;

            r ^= b & n;
        }
        return r;
    }

    inline uint64_t mul(uint64_t a, uint64_t b) {
        uint64_t r = 0;
        for (int i = 7; i >= 0; --i) {
            uint64_t m = r & 0x8080808080808080;
            m = m - (m >> 7);

            r = ((r & 0x7f7f7f7f7f7f7f7f) << 1) ^ (poly8 & m);

            uint64_t n = (a & (0x0101010101010101 << i)) >> i;
            n = (n << 8) - n;

            r ^= b & n;
        }
        return r;
    }

    inline uint8_t inv(uint8_t a) {
        return exp_table[255 - log_table[a]];
    }

    inline uint8_t div(uint8_t a, uint8_t b) {
        if (a == 0)
            return 0;

        unsigned r = log_table[a] + 255 - log_table[b];
        if (r >= 255)
            r -= 255;

        return exp_table[r];
    }

    inline uint8_t exp(uint8_t a) {
        return exp_table[a];
    }

    inline uint8_t log(uint8_t a) {
        return log_table[a];
    }

    inline unsigned ex_synth_div(uint8_t a[], unsigned size_a, const uint8_t b[], unsigned size_b) {
        // uint8_t normalizer = b[0];
        // assert(b[0] == 1);

        if (size_b > size_a)
            return 0;

        for (unsigned i = 0; i < size_a - size_b + 1; ++i) {
            uint8_t coef = a[i];// = div(a[i], normalizer);

            if (coef == 0)
                continue;

            for (unsigned j = 1; j < size_b; ++j)
                a[i + j] ^= mul(b[j], coef);
        }

        return size_a - size_b + 1;
    }

    inline unsigned poly_mod(uint8_t r[], const uint8_t a[], unsigned size_a, const uint8_t b[], unsigned size_b) {
        // assert(b[0] == 1);
        b += 1;

        if (size_b < 2)
            return 0;

        const unsigned size_r = size_b - 1;
        std::copy_n(a, std::min(size_a, size_r), r);

        if (size_a < size_b)
            return size_a;

        for (unsigned i = 0; i < size_a - size_r; ++i) {
            uint8_t c = r[0];
            std::rotate(r, r + 1, r + size_r);
            r[size_r - 1] = a[size_r + i];

            if (c == 0)
                continue;

            for (unsigned j = 0; j < size_r; ++j)
                r[j] ^= mul(b[j], c);
        }

        return size_r;
    }

    inline void poly_mod_x_n(uint8_t rem[], const uint8_t a[], const unsigned size_a, const uint8_t b[], const unsigned size_b) {
        if (size_a >= size_b) {
            std::copy_n(a, size_b, rem);
            for (unsigned i = 0; i < size_a - size_b; ++i) {
                uint8_t c = rem[0];
                rem[0] = a[size_b + i];
                std::rotate(rem, rem + 1, rem + size_b);

                if (c == 0)
                    continue;

                for (unsigned j = 0; j < size_b; ++j)
                    rem[j] ^= mul(b[j], c);
            }
        } else {
            std::memset(rem, 0, size_b - size_a);
            std::copy_n(a, size_a, &rem[size_b - size_a]);
        }

        for (unsigned i = 0; i < size_b; ++i) {
            uint8_t c = rem[0];
            rem[0] = 0;
            std::rotate(rem, rem + 1, rem + size_b);

            if (c == 0)
                continue;

            for (unsigned j = 0; j < size_b; ++j)
                rem[j] ^= mul(b[j], c);
        }
    }

    inline uint8_t poly_eval(const uint8_t poly[], const unsigned size, const uint8_t x) {
        uint8_t r = 0;
        for (unsigned i = 0; i < size; ++i)
            r = mul(r, x) ^ poly[i];
        return r;
    }

    inline uint32_t poly_eval(const uint8_t poly[], const unsigned size, const uint32_t x) {
        uint32_t r = 0;
        for (unsigned i = 0; i < size; ++i)
            r = mul(r, x) ^ (poly[i] * 0x01010101);
        return r;
    }

    inline uint64_t poly_eval(const uint8_t poly[], const unsigned size, const uint64_t x) {
        uint64_t r = 0;
        for (unsigned i = 0; i < size; ++i)
            r = mul(r, x) ^ (poly[i] * 0x0101010101010101);
        return r;
    }

    inline void poly_shift(uint8_t poly[], const unsigned size, const unsigned n) {
        unsigned i = 0;
        for (; i < size - n; ++i)
            poly[i] = poly[i + n];
        for (; i < size; ++i)
            poly[i] = 0;
    }

    inline void poly_scale(uint8_t poly[], const unsigned size, const uint8_t a) {
        for (unsigned i = 0; i < size; ++i)
            poly[i] = mul(a, poly[i]);
    }

    inline void poly_add(uint8_t poly_a[], const uint8_t poly_b[], const unsigned size) {
        for (unsigned i = 0; i < size; ++i)
            poly_a[i] ^= poly_b[i];
    }

    inline unsigned poly_mul(uint8_t r[],
            const uint8_t poly_a[], const unsigned size_a,
            const uint8_t poly_b[], const unsigned size_b) {
        unsigned res_len = size_a + size_b - 1;

        for (unsigned i = 0; i < res_len; ++i)
            r[i] = 0;

        for (unsigned i = 0; i < size_a; ++i)
            for (unsigned j = 0; j < size_b; ++j)
                r[i + j] ^= mul(poly_a[i], poly_b[j]);

        return res_len;
    }

    inline unsigned poly_deriv(uint8_t poly[], unsigned size) {
        if (size & 1) {
            for (unsigned i = 0; i < size - 1; ++i)
                poly[i] = (i & 1) ? 0 : poly[i + 1];

            return size - 2;
        } else {
            for (unsigned i = 1; i < size - 1; i += 2)
                poly[i] = 0;

            return size - 1;
        }
    }
};
