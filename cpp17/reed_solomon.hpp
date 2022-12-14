#pragma once

#include <algorithm>
#include <functional>

#include "galois.hpp"

namespace detail {
    template<typename T, unsigned Bits>
    struct bit_array {
        static_assert(Bits <= sizeof(T) * 8);

        uint8_t *ptr;

        struct bit_field {
            uint8_t *ptr;
            uint8_t offset;

            void operator=(T const& in) {
                auto ptr = this->ptr;
                auto offset = this->offset;

                for (unsigned i = 0; i < Bits;) {
                    uint8_t val;
                    uint8_t mask = 0xff >> offset;

                    uint8_t rem_bits = Bits - i;
                    uint8_t avail_bits = 8 - offset;

                    if (rem_bits < avail_bits) {
                        val = in << (avail_bits - rem_bits);
                        mask &= 0xff << (avail_bits - rem_bits);
                    } else {
                        val = in >> (rem_bits - avail_bits);
                    }

                    *ptr++ = (*ptr & ~mask) | (val & mask);
                    i += avail_bits;
                    offset = 0;
                }
            }

            operator T() const {
                auto ptr = this->ptr;
                auto offset = this->offset;
                T ret = 0;

                for (unsigned i = 0; i < Bits;) {
                    uint8_t mask = 0xff >> offset;

                    uint8_t rem_bits = Bits - i;
                    uint8_t avail_bits = 8 - offset;

                    if (rem_bits < avail_bits) {
                        mask &= 0xff << (avail_bits - rem_bits);

                        ret <<= rem_bits;
                        ret |= (*ptr++ & mask) >> (avail_bits - rem_bits);
                    } else {
                        ret <<= avail_bits;
                        ret |= (*ptr++ & mask);
                    }

                    i += avail_bits;
                    offset = 0;
                }
                return ret;
            }
        };

        bit_field operator[](int i) {
            auto bit_i = i * Bits;
            return bit_field{ptr + (bit_i >> 3), uint8_t(bit_i & 0x07)};
        }
    };
}


template<typename _GF, unsigned Ecc>
struct rs_base {
    using GF = _GF;
    static constexpr auto ecc = Ecc;
};

template<typename RS>
struct rs_generator {
    static inline constexpr struct sdata_t {
        using GFT = typename RS::GF::Repr;

        GFT generator[RS::ecc + 1] = {};
        GFT roots[RS::ecc] = {};

        inline constexpr sdata_t() {
            GFT temp[RS::ecc + 1] = {};

            auto p1 = (RS::ecc & 1) ? &generator[0] : &temp[0];
            auto p2 = (RS::ecc & 1) ? &temp[0] : &generator[0];

            unsigned len = 1;
            p2[0] = 1;

            for (unsigned i = 0; i < RS::ecc; ++i) {
                roots[i] = RS::GF::exp(i);
                GFT factor[] = {1, RS::GF::sub(0, roots[i])};
                len = RS::GF::poly_mul(p1, p2, len, factor, 2);

                auto t = p1;
                p1 = p2;
                p2 = t;
            }
        }
    } sdata{};
};

template<typename RS>
struct rs_encode_basic {
    static constexpr auto& generator = rs_generator<RS>::sdata.generator;

    template<typename S, typename T>
    static inline void encode(S output, T const& input, unsigned input_size) {
        RS::GF::poly_mod_x_n(output, input, input_size, &generator[1], RS::ecc);

        if constexpr (RS::GF::prime != 2) {
            for (unsigned i = 0; i < RS::ecc; ++i)
                output[i] = RS::GF::sub(0, output[i]);
        }
    }
};

template<typename RS>
struct rs_encode_lut {
    static_assert(RS::GF::prime == 2);

    static constexpr auto& generator = rs_generator<RS>::sdata.generator;

    static inline constexpr struct sdata_t {
        uint8_t generator_lut[RS::GF::charact][RS::ecc] = {};

        constexpr inline sdata_t() {
            for (unsigned i = 0; i < RS::GF::charact; ++i) {
                uint8_t data[RS::ecc + 1] = {0};
                data[0] = uint8_t(i);
                RS::GF::ex_synth_div(&data[0], RS::ecc + 1, &generator[0], RS::ecc + 1);

                for (unsigned j = 0; j < RS::ecc; ++j)
                    generator_lut[i][j] = data[j + 1];
            }
        }
    } sdata{};

    static inline void encode(uint8_t *output, const uint8_t *data, unsigned size) {
        std::fill_n(output, RS::ecc, 0x00);
        for (unsigned i = 0; i < size; ++i) {
            uint8_t pos = output[0] ^ data[i];
            output[0] = 0;
            std::rotate(output, output + 1, output + RS::ecc);
            std::transform(output, output + RS::ecc,
                    &sdata.generator_lut[pos][0],
                    output, std::bit_xor());
        }
    }
};

template<typename Word, unsigned N>
struct rs_encode_slice {
    template<typename RS>
    struct type {
        using GFT = typename RS::GF::Repr;
        static_assert(RS::ecc == sizeof(Word));
        static_assert(RS::GF::prime == 2);

        static constexpr auto& generator = rs_generator<RS>::sdata.generator;

        static inline constexpr struct sdata_t {
            Word generator_lut[N][RS::GF::charact] = {};

            inline constexpr sdata_t() {
                for (unsigned i = 0; i < RS::GF::charact; ++i) {
                    uint8_t data[RS::ecc + 1] = {0};
                    data[0] = uint8_t(i);
                    RS::GF::ex_synth_div(&data[0], RS::ecc + 1, &generator[0], RS::ecc + 1);

                    // endianess dependent
                    for (unsigned j = 0; j < RS::ecc; ++j)
                        generator_lut[0][i] |= Word(data[j + 1]) << j * 8;
                }

                for (unsigned i = 0; i < RS::GF::charact; ++i) {
                    for (unsigned j = 1; j < N; ++j)
                        generator_lut[j][i] =
                                (generator_lut[j - 1][i] >> 8) ^
                                generator_lut[0][generator_lut[j - 1][i] % RS::GF::charact];
                }
            }
        } sdata{};

        static inline void encode(uint8_t *output, const uint8_t *data, unsigned size) {
            unsigned i = 0;
            Word rem = 0;

            if constexpr (N > 1) {
                static_assert(N % RS::ecc == 0);

                for (; size - i >= N; i += N) {
                    auto in = reinterpret_cast<const Word *>(&data[i]);
                    rem ^= *in++;
                    Word t = 0;
                    for (unsigned j = 0; j < RS::ecc; ++j)
                        t ^= sdata.generator_lut[N - j - 1][(rem >> (8 * j)) & 0xff];

                    for (unsigned k = 1; k < N / RS::ecc; ++k) {
                        auto two = *in++;
                        for (unsigned j = 0; j < RS::ecc; ++j)
                            t ^= sdata.generator_lut[N - k * RS::ecc - j - 1][(two >> (8 * j)) & 0xff];
                    }
                    rem = t;
                }
            }

            for (; i < size; ++i)
                rem = (rem >> 8) ^ sdata.generator_lut[0][(rem & 0xff) ^ data[i]];

            // endianess dependent
            for (unsigned i = 0; i < RS::ecc; ++i)
                output[i] = rem >> (8 * i);
        }
    };
};

template<typename RS>
struct rs_synds_basic {
    using GFT = typename RS::GF::Repr;
    using synds_array_t = GFT[RS::ecc];

    static constexpr auto& gen_roots = rs_generator<RS>::sdata.roots;

    template<typename S, typename T>
    static inline void synds(synds_array_t synds, S const& data, unsigned size, T const& rem) {
        for (unsigned i = 0; i < RS::ecc; ++i) {
            auto t = RS::GF::poly_eval(data, size, gen_roots[i]);
            synds[RS::ecc - i - 1] = RS::GF::poly_eval(rem, RS::ecc, gen_roots[i], t);
        }
    }
};


template<typename Word>
struct rs_synds_lut_t {
    template<typename RS>
    struct type {
        static constexpr auto ecc_w = (RS::ecc / sizeof(Word)) + !!(RS::ecc % sizeof(Word));
        static constexpr auto synds_size = ecc_w * sizeof(Word);
        using synds_array_t /* alignas(sizeof(Word)) */ = uint8_t[synds_size];

        static inline constexpr struct sdata_t {
            Word gen_roots[ecc_w] = {};

            inline constexpr sdata_t() {
                for (unsigned i = 0; i < RS::ecc; ++i)
                    // reverse-order syndromes, endianess dependent
                    gen_roots[i / sizeof(Word)] |= Word(RS::GF::exp(RS::ecc - i - 1)) << (i % sizeof(Word)) * 8;
            }
        } sdata{};

        static inline void synds(uint8_t synds[], const uint8_t *data, unsigned size, const uint8_t *rem) {
            auto synds_w = reinterpret_cast<Word *>(synds);
            for (unsigned i = 0; i < ecc_w; ++i) {
                auto t = gf_wide_mul<typename RS::GF, Word>::poly_eval(data, size, sdata.gen_roots[i]);
                synds_w[i] = gf_wide_mul<typename RS::GF, Word>::poly_eval(rem, RS::ecc, sdata.gen_roots[i], t);
            }
        }
    };
};
template<typename RS>
using rs_synds_lut4 = rs_synds_lut_t<uint32_t>::type<RS>;
template<typename RS>
using rs_synds_lut8 = rs_synds_lut_t<uint64_t>::type<RS>;


template<typename RS>
struct rs_roots_eval_basic {
    using GFT = typename RS::GF::Repr;

    static inline unsigned roots(
            const GFT poly[], unsigned poly_size,
            GFT roots[], unsigned size)
    {
        unsigned count = 0;

        for (unsigned i = 0; i < size; ++i) {
            if (RS::GF::poly_eval(poly, poly_size, RS::GF::inv(RS::GF::exp(i))) == 0)
                roots[count++] = i;
        }

        return count;
    }
};

template<typename RS>
struct rs_roots_eval_chien {
    static inline unsigned roots(
            const uint8_t poly[], unsigned poly_size,
            uint8_t roots[], unsigned size)
    {
        unsigned count = 0;

        uint8_t coefs[RS::ecc];
        std::reverse_copy(&poly[0], &poly[poly_size], &coefs[0]);

        for (int i = 254; i >= 0; --i) {
            uint8_t sum = 1;

            for (unsigned j = 1; j < poly_size; ++j) {
                coefs[j] = RS::GF::mul(coefs[j], RS::GF::exp(j));
                sum = RS::GF::add(sum, coefs[j]);
            }

            if (sum == 0) {
                roots[count] = i;
                if (++count >= poly_size-1)
                    break;
            }
        }

        return count;
    }
};

template<typename Word>
struct rs_roots_eval_lut_t {
    template<typename RS>
    struct type {
        static constexpr auto ecc_w = (RS::ecc / sizeof(Word)) + !!(RS::ecc % sizeof(Word));

        static inline constexpr struct sdata_t {
            union {
                std::array<uint8_t, 256> u8;
                std::array<Word, 256 / sizeof(Word)> word;
            } err_poly_roots{};

            inline constexpr sdata_t() {
                for (unsigned i = 0; i < 256; ++i)
                    err_poly_roots.u8[i] = RS::GF::inv(RS::GF::exp(i));
            }
        } sdata{};

        static inline unsigned roots(
                const uint8_t poly[], unsigned poly_size,
                uint8_t roots[], unsigned size)
        {
            unsigned count = 0;

            for (unsigned i = 0; i <= size/sizeof(Word); ++i) {
                Word eval = gf_wide_mul<typename RS::GF, Word>::poly_eval(poly, poly_size, sdata.err_poly_roots.word[i]);
                for (unsigned j = 0; j < sizeof(Word); ++j) {
                    if (reinterpret_cast<uint8_t *>(&eval)[j] == 0) {
                        auto pos = i * sizeof(Word) + j;

                        if (pos < size)
                            roots[count++] = pos;
                        else
                            break;
                    }
                }
            }

            return count;
        }
    };
};

template<typename RS>
using rs_roots_eval_lut4 = rs_roots_eval_lut_t<uint32_t>::type<RS>;
template<typename RS>
using rs_roots_eval_lut8 = rs_roots_eval_lut_t<uint64_t>::type<RS>;

template<typename RS>
struct rs_decode {
    using GFT = typename RS::GF::Repr;

    template<typename T, typename U>
    static inline bool decode(T data, unsigned size, U rem) {
        typename RS::synds_array_t synds;
        RS::synds(synds, &data[0], size, rem);

        if (std::all_of(&synds[0], &synds[RS::ecc], std::logical_not()))
            return true;

        GFT err_poly[RS::ecc];
        auto errors = berlekamp_massey(synds, err_poly);

        GFT err_pos[RS::ecc / 2];
        auto roots = RS::roots(&err_poly[RS::ecc-errors-1], errors+1, err_pos, size + RS::ecc);

        if (errors != roots)
            return false;

        GFT err_mag[RS::ecc / 2];
        forney(synds, &err_poly[RS::ecc-errors-1], err_pos, errors, err_mag);

        for (unsigned i = 0; i < errors; ++i) {
            unsigned pos = size + RS::ecc - 1 - err_pos[i];
            if (pos >= size + RS::ecc)
                return false;

            if (pos < size)
                data[pos] = RS::GF::add(data[pos], err_mag[i]);
            else
                rem[pos - size] = RS::GF::add(rem[pos - size], err_mag[i]);
        }

        return true;
    }

    template<typename T, typename U, typename V>
    static inline bool decode(T data, unsigned size, U rem, const V err_idx, unsigned errors) {
        if (errors > RS::ecc)
            return false;

        typename RS::synds_array_t synds;
        RS::synds(synds, &data[0], size, rem);

        if (std::all_of(&synds[0], &synds[RS::ecc], std::logical_not()))
            return true;

        GFT err_pos[RS::ecc];
        for (unsigned i = 0; i < errors; ++i) {
            if (err_idx[i] > size + RS::ecc - 1)
                return false;

            err_pos[i] = size + RS::ecc - 1 - err_idx[i];
        }

        GFT err_poly[RS::ecc + 1] = {1};
        unsigned err_poly_len = 1;

        {
            GFT temp[RS::ecc + 1] = {1};

            auto p1 = (errors & 1) ? &err_poly[0] : &temp[0];
            auto p2 = (errors & 1) ? &temp[0] : &err_poly[0];

            for (unsigned i = 0; i < errors; ++i) {
                GFT factor[2] = {RS::GF::sub(0, RS::GF::exp(err_pos[i])), 1};
                err_poly_len = RS::GF::poly_mul(p1, p2, err_poly_len, factor, 2);
                std::swap(p1, p2);
            }
        }

        assert(err_poly_len == errors + 1);

        GFT err_mag[RS::ecc];
        forney(synds, err_poly, err_pos, errors, err_mag);

        for (unsigned i = 0; i < errors; ++i) {
            unsigned pos = err_idx[i];
            if (pos >= size + RS::ecc)
                return false;

            if (pos < size)
                data[pos] = RS::GF::add(data[pos], err_mag[i]);
            else
                rem[pos - size] = RS::GF::add(rem[pos - size], err_mag[i]);
        }

        return true;
    }

    static inline unsigned berlekamp_massey(const GFT synds_rev[RS::ecc], GFT err_poly[RS::ecc]) {
        GFT prev[RS::ecc] = {};
        GFT temp[RS::ecc];
        GFT synds[RS::ecc];
        std::reverse_copy(synds_rev, &synds_rev[RS::ecc], synds);
        std::fill_n(err_poly, RS::ecc, 0);

        prev[RS::ecc-1] = 1;
        err_poly[RS::ecc-1] = 1;

        unsigned errors = 0;
        unsigned m = 1;
        GFT b = 1;

        for (unsigned n = 0; n < RS::ecc; ++n) {
            unsigned d = synds[n]; // discrepancy
            for (unsigned i = 1; i < errors + 1; ++i)
                d = RS::GF::add(d, RS::GF::mul(err_poly[RS::ecc - 1 - i], synds[n-i]));

            if (d == 0) {  // discrepancy is already zero
                m = m + 1;
            } else if (2 * errors <= n) {
                std::copy_n(err_poly, RS::ecc, temp);

                RS::GF::poly_shift(prev, RS::ecc, m);
                RS::GF::poly_scale(prev, RS::ecc, RS::GF::div(d, b));

                RS::GF::poly_sub(err_poly, prev, RS::ecc);

                errors = n + 1 - errors;
                std::copy_n(temp, RS::ecc, prev);

                b = d;
                m = 1;
            } else {
                std::copy_n(prev, RS::ecc, temp);

                RS::GF::poly_shift(temp, RS::ecc, m);

                RS::GF::poly_scale(temp, RS::ecc, RS::GF::div(d, b));
                RS::GF::poly_sub(err_poly, temp, RS::ecc);

                m = m + 1;
            }
        }

        return errors;
    }

    static inline void forney(
            const GFT synds_rev[RS::ecc], GFT err_poly[], const GFT err_pos[],
            const unsigned err_count, GFT err_mag[])
    {
        GFT err_eval[RS::ecc * 2];
        auto err_eval_size = RS::GF::poly_mul(err_eval,
                synds_rev, RS::ecc,
                err_poly, err_count + 1);

        GFT x_poly[RS::ecc + 1] = {1};
        auto err_eval_begin = RS::GF::ex_synth_div(
                err_eval, err_eval_size,
                x_poly, RS::ecc + 1);
        while (err_eval[err_eval_begin] == 0) {
            err_eval_begin++;
            assert(err_eval_begin < err_eval_size);
        }
        err_eval_size = err_eval_size - err_eval_begin;

        RS::GF::poly_deriv(err_poly, err_count + 1);
        auto err_poly_deriv = err_poly + 1;
        auto err_poly_deriv_size = err_count;

        for (unsigned i = 0; i < err_count; ++i) {
            auto xi = RS::GF::exp(err_pos[i]);
            auto xi_inv = RS::GF::inv(xi);

            auto n = RS::GF::poly_eval(&err_eval[err_eval_begin], err_eval_size, xi_inv);
            auto d = RS::GF::poly_eval(err_poly_deriv, err_poly_deriv_size, xi_inv);

            err_mag[i] = RS::GF::mul(xi, RS::GF::div(n, d));
        }
    }
};

template<typename GF, unsigned Ecc, template<class>typename...Fs>
struct rs_impl_base : rs_base<GF, Ecc>, Fs<rs_base<GF, Ecc>>... { };

template<typename Impl, template<class>typename...Fs>
struct rs_impl : rs_base<typename Impl::GF, Impl::ecc>, Fs<Impl>... {
    static inline auto static_data_size = detail::get_sdata_size<rs_generator<Impl>, Fs<Impl>...>();
};

template<typename GF, unsigned Ecc, template<class>typename...Fs>
struct RS : rs_impl<rs_impl_base<GF, Ecc, Fs...>, Fs...> {
    static_assert(Ecc < 255);
};

