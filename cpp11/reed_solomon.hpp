#pragma once

#include <algorithm>
#include <cstring>
#include <functional>
#include <memory>

#include "galois.hpp"

template<unsigned ecc, uint8_t gf_base = 2, uint16_t gf_poly = 0x11d, typename Word = unsigned>
struct RS {
    static_assert(ecc < 255);

    static const unsigned ecc_w = (ecc / sizeof(Word)) + !!(ecc % sizeof(Word));

#ifndef RS_ENCODE_ONLY
    GF<gf_base, gf_poly> gf;
#endif
    std::array<uint8_t, ecc + 1> generator;
    alignas(Word) std::array<uint8_t, ecc_w * sizeof(Word)> generator_roots;
#if defined(RS_GENERATOR_LUT) || defined(RS_ENCODE_ONLY)
    alignas(Word) std::array<std::array<uint8_t, ecc_w * sizeof(Word)>, 256> generator_lut;
#endif
#ifdef RS_POLY_ROOT_LUT
    alignas(Word) std::array<uint8_t, 256> err_poly_roots;
#endif

    inline RS() {
#ifdef RS_ENCODE_ONLY
        GF<gf_base, gf_poly> gf;
#endif

        std::array<uint8_t, ecc + 1> temp = {};

        auto p1 = (ecc & 1) ? &generator[0] : &temp[0];
        auto p2 = (ecc & 1) ? &temp[0] : &generator[0];

        unsigned len = 1;
        p2[0] = 1;

        for (unsigned i = 0; i < ecc; ++i) {
            uint8_t factor[] = {1, gf.exp(i)};
            len = gf.poly_mul(p1, p2, len, factor, 2);

            auto t = p1;
            p1 = p2;
            p2 = t;
        }

        for (unsigned i = 0; i < ecc; ++i)
            generator_roots[i] = gf.exp(i);

#if defined(RS_GENERATOR_LUT) || defined(RS_ENCODE_ONLY)
        for (unsigned i = 0; i < 256; ++i) {
            uint8_t data[ecc + 1] = {uint8_t(i)};
            gf.ex_synth_div(&data[0], ecc + 1, &generator[0], ecc + 1);

            if (ecc == sizeof(Word)) {
                for (unsigned j = 0; j < ecc; ++j)
                    generator_lut[i][ecc - 1 - j] = data[j + 1];
            } else {
                for (unsigned j = 0; j < ecc; ++j)
                    generator_lut[i][j] = data[j + 1];
            }
        }
#endif

#ifdef RS_POLY_ROOT_LUT
        for (unsigned i = 0; i < 256; ++i)
            err_poly_roots[i] = gf.inv(gf.exp(i));
#endif
    }

#if defined(RS_GENERATOR_LUT) || defined(RS_ENCODE_ONLY)
    inline void encode(uint8_t *data, unsigned size) {
        auto data_len = size - ecc;
        auto rem = &data[data_len];

        if (ecc == sizeof(Word)) {
            auto lut = reinterpret_cast<const Word *>(&generator_lut[0][0]);
            const unsigned shift = (sizeof(Word) - 1) * 8;

            Word w = 0;
            for (unsigned i = 0; i < data_len; ++i) {
                uint8_t pos = (w ^ (Word(data[i]) << shift)) >> shift;

                w = (w << 8) ^ lut[pos];
            }
        
            // std::copy_n(reinterpret_cast<const uint8_t *>(&w), ecc, rem);
            for (unsigned i = 0; i < ecc; ++i)
                rem[i] = w >> (shift - 8 * i);
        } else {
            std::memset(rem, 0, ecc);
            for (unsigned i = 0; i < data_len; ++i) {
                uint8_t pos = rem[0] ^ data[i];
                rem[0] = 0;
                std::rotate(rem, rem + 1, rem + ecc);
                std::transform(rem, rem + ecc,
                        &generator_lut[pos][0],
                        rem, std::bit_xor<uint8_t>());
            }
        }
    }
#else
    inline void encode(uint8_t *data, unsigned size) {
        auto data_len = size - ecc;
        gf.poly_mod_x_n(&data[data_len], data, data_len, &generator[1], ecc);
    }
#endif

#ifndef RS_ENCODE_ONLY
    inline void decode(uint8_t *data, unsigned size) {
        Word synds_w[ecc_w];
        auto synds = reinterpret_cast<uint8_t *>(synds_w);
        auto gen_roots = reinterpret_cast<const Word *>(&generator_roots[0]);

        for (unsigned i = 0; i < ecc_w; ++i)
            synds_w[i] = gf.poly_eval(data, size, gen_roots[i]);

        if (std::all_of(synds, synds + ecc, std::logical_not<uint8_t>()))
            return;

        uint8_t err_poly[ecc];
        auto errors = berlekamp_massey(synds, err_poly);

        uint8_t err_pos[ecc / 2];
        find_roots(&err_poly[ecc-errors-1], errors+1, err_pos, size);

        uint8_t err_mag[ecc / 2];
        forney(synds, &err_poly[ecc-errors-1], err_pos, errors, err_mag);

        for (unsigned i = 0; i < errors; ++i)
            data[size - 1 - err_pos[i]] ^= err_mag[i];
    }

    inline unsigned berlekamp_massey(const uint8_t synds[ecc], uint8_t err_poly[ecc]) {
        uint8_t prev[ecc];
        uint8_t temp[ecc];

        for (unsigned i = 0; i < ecc; ++i) {
            err_poly[i] = 0;
            prev[i] = 0;
        }
        prev[ecc-1] = 1;
        err_poly[ecc-1] = 1;

        unsigned errors = 0;
        unsigned m = 1;
        uint8_t b = 1;

        for (unsigned n = 0; n < ecc; ++n) {
            unsigned d = synds[n]; // discrepancy
            for (unsigned i = 1; i < errors + 1; ++i)
                d ^= gf.mul(err_poly[ecc - 1 - i], synds[n-i]);

            if (d == 0) {  // discrepancy is already zero
                m = m + 1;
            } else if (2 * errors <= n) {
                std::copy_n(err_poly, ecc, temp);

                gf.poly_shift(prev, ecc, m);
                gf.poly_scale(prev, ecc, gf.div(d, b));

                gf.poly_add(err_poly, prev, ecc);

                errors = n + 1 - errors;
                std::copy_n(temp, ecc, prev);

                b = d;
                m = 1;
            } else {
                std::copy_n(prev, ecc, temp);

                gf.poly_shift(temp, ecc, m);

                gf.poly_scale(temp, ecc, gf.div(d, b));
                gf.poly_add(err_poly, temp, ecc);

                m = m + 1;
            }
        }

        return errors;
    }

    inline unsigned find_roots(
            const uint8_t poly[], unsigned poly_size,
            uint8_t roots[], unsigned size)
    {
        unsigned count = 0;

#ifdef RS_POLY_ROOT_LUT
        auto err_poly_roots_64 = reinterpret_cast<const Word *>(&err_poly_roots[0]);
        for (unsigned i = 0; i <= size/sizeof(Word); ++i) {
            Word eval = gf.poly_eval(poly, poly_size, err_poly_roots_64[i]);
            for (unsigned j = 0; j < 8; ++j) {
                if (reinterpret_cast<uint8_t *>(&eval)[j] == 0) {
                    auto pos = i * 8 + j;

                    if (pos < size)
                        roots[count++] = pos;
                    else
                        break;
                }
            }
        }
#else
        for (unsigned i = 0; i < size; ++i) {
            if (gf.poly_eval(poly, poly_size, gf.inv(gf.exp(i))) == 0)
                roots[count++] = i;
        }
#endif

        return count;
    }

    inline void forney(
            const uint8_t synds[ecc], uint8_t err_poly[], const uint8_t err_pos[],
            const unsigned err_count, uint8_t err_mag[])
    {
        uint8_t temp[ecc] = {1};

        uint8_t err_eval[ecc * 2];
        uint8_t synds_rev[ecc];
        std::reverse_copy(synds, &synds[ecc], synds_rev);

        auto err_eval_size = gf.poly_mul(
                err_eval, synds_rev,
                ecc, err_poly, err_count + 1);

        auto err_eval_begin = gf.ex_synth_div(
                err_eval, err_eval_size,
                temp, ecc);

        while (err_eval[err_eval_begin] == 0) {
            err_eval_begin++;
            assert(err_eval_begin < err_eval_size);
        }
        err_eval_size = err_eval_size - err_eval_begin;

        auto err_poly_deriv_size = gf.poly_deriv(err_poly, err_count + 1);

        for (unsigned i = 0; i < err_count; ++i) {
            auto xi = gf.exp(err_pos[i]);
            auto xi_inv = gf.inv(xi);

            auto y = gf.poly_eval(&err_eval[err_eval_begin], err_eval_size, xi_inv);
            auto d = gf.poly_eval(err_poly, err_poly_deriv_size, xi_inv);

            err_mag[i] = gf.mul(xi, gf.div(y, d));
        }
    }

#endif
};
