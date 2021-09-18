#pragma once

#include <algorithm>
#include <cstring>
#include <functional>
#include <memory>

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


template<typename GF_t, unsigned Ecc>
struct rs_base {
    using GF = GF_t;
    static constexpr auto ecc = Ecc;
};

template<typename RS>
struct rs_generator {
    static inline constexpr struct sdata_t {
        std::array<typename RS::GF::Repr, RS::ecc + 1> generator{};

        inline constexpr sdata_t() {
            std::array<typename RS::GF::Repr, RS::ecc + 1> temp{};

            auto p1 = (RS::ecc & 1) ? &generator[0] : &temp[0];
            auto p2 = (RS::ecc & 1) ? &temp[0] : &generator[0];

            unsigned len = 1;
            p2[0] = 1;

            for (unsigned i = 0; i < RS::ecc; ++i) {
                typename RS::GF::Repr factor[] = {1, RS::GF::exp(i)};
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
    static inline void encode(S output, T const& input, unsigned size) {
        auto data_len = size - RS::ecc;
        RS::GF::poly_mod_x_n(output, input, data_len, &generator[1], RS::ecc);
    }
};

template<typename Word>
struct rs_encode_lut_t {
    template<typename RS>
    struct type {
        static constexpr auto ecc_w = (RS::ecc / sizeof(Word)) + !!(RS::ecc % sizeof(Word));
        static constexpr auto& generator = rs_generator<RS>::sdata.generator;

        static inline constexpr struct sdata_t {
            union {
                std::array<std::array<uint8_t, sizeof(Word)>, RS::GF::charact> u8;
                std::array<Word, RS::GF::charact> word;
            } generator_lut{};

            inline constexpr sdata_t() {
                for (unsigned i = 0; i < RS::GF::charact; ++i) {
                    uint8_t data[RS::ecc + 1] = {uint8_t(i)};
                    RS::GF::ex_synth_div(&data[0], RS::ecc + 1, &generator[0], RS::ecc + 1);

                    static_assert(RS::ecc == sizeof(Word));

                    for (unsigned j = 0; j < RS::ecc; ++j)
                        generator_lut.u8[i][RS::ecc - 1 - j] = data[j + 1]; // little endian
                }
            }
        } sdata{};

        static inline void encode(uint8_t *data, unsigned size) {
            auto data_len = size - RS::ecc;
            auto rem = &data[data_len];

            static_assert(RS::ecc == sizeof(Word));

            const auto shift = (sizeof(Word) - 1) * 8;

            Word w = 0;
            for (unsigned i = 0; i < data_len; ++i) {
                auto pos = (w >> shift) ^ data[i];

                w = (w << 8) ^ sdata.generator_lut.word[pos];
            }

            // std::copy_n(reinterpret_cast<const uint8_t *>(&w), RS::ecc, rem);
            for (unsigned i = 0; i < RS::ecc; ++i)
                rem[i] = w >> (shift - 8 * i);
        }
    };
};

template<>
struct rs_encode_lut_t<uint8_t> {
    template<typename RS>
    struct type {
        static constexpr auto& generator = rs_generator<RS>::sdata.generator;

        static inline constexpr struct sdata_t {
            std::array<std::array<uint8_t, RS::ecc>, 256> generator_lut{};

            constexpr inline sdata_t() {
                for (unsigned i = 0; i < 256; ++i) {
                    uint8_t data[RS::ecc + 1] = {uint8_t(i)};
                    RS::GF::ex_synth_div(&data[0], RS::ecc + 1, &generator[0], RS::ecc + 1);

                    for (unsigned j = 0; j < RS::ecc; ++j)
                        generator_lut[i][j] = data[j + 1];
                }
            }
        } sdata{};

        static inline void encode(uint8_t *data, unsigned size) {
            auto data_len = size - RS::ecc;
            auto rem = &data[data_len];
            std::memset(rem, 0, RS::ecc);
            for (unsigned i = 0; i < data_len; ++i) {
                uint8_t pos = rem[0] ^ data[i];
                rem[0] = 0;
                std::rotate(rem, rem + 1, rem + RS::ecc);
                std::transform(rem, rem + RS::ecc,
                        &sdata.generator_lut[pos][0],
                        rem, std::bit_xor<uint8_t>());
            }
        }
    };
};

template<typename RS>
using rs_encode_lut = rs_encode_lut_t<uint8_t>::type<RS>;
template<typename RS>
using rs_encode_lut4 = rs_encode_lut_t<uint32_t>::type<RS>;
template<typename RS>
using rs_encode_lut8 = rs_encode_lut_t<uint64_t>::type<RS>;



template<typename RS>
struct rs_synds_basic {
    using GFT = typename RS::GF::Repr;
    using synds_array_t = GFT[RS::ecc];

    static inline constexpr struct sdata_t {
        std::array<GFT, RS::ecc> gen_roots{};

        inline constexpr sdata_t() {
            for (unsigned i = 0; i < RS::ecc; ++i)
                gen_roots[i] = RS::GF::exp(i);
        }
    } sdata{};

    template<typename S>
    static inline void synds(S const& data, unsigned size, synds_array_t synds) {
        for (unsigned i = 0; i < RS::ecc; ++i)
            synds[i] = RS::GF::poly_eval(data, size, sdata.gen_roots[i]);
    }
};


template<typename Word>
struct rs_synds_lut_t {
    template<typename RS>
    struct type {
        static constexpr auto ecc_w = (RS::ecc / sizeof(Word)) + !!(RS::ecc % sizeof(Word));
        static constexpr auto synds_size = ecc_w * sizeof(Word);
        static constexpr auto synds_align = sizeof(Word);
        using synds_array_t alignas(sizeof(Word)) = uint8_t[synds_size];

        static inline constexpr struct sdata_t {
            union {
                std::array<uint8_t, synds_size> u8;
                std::array<Word, ecc_w> word;
            } gen_roots{};

            inline constexpr sdata_t() {
                for (unsigned i = 0; i < RS::ecc; ++i)
                    gen_roots.u8[i] = Word(RS::GF::exp(i));
            }
        } sdata{};

        static inline void synds(const uint8_t *data, unsigned size, uint8_t synds[]) {
            auto synds_w = reinterpret_cast<Word *>(synds);
            for (unsigned i = 0; i < ecc_w; ++i)
                synds_w[i] = gf_wide_mul<typename RS::GF, Word>::poly_eval(data, size, sdata.gen_roots.word[i]);
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

            for (unsigned j = 1; j < poly_size; ++j)
                sum ^= coefs[j] = RS::GF::mul(coefs[j], RS::GF::exp(j));

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

    template<typename T>
    static inline void decode(T& data, unsigned size) {
        typename RS::synds_array_t synds;
        RS::synds(data, size, synds);

        if (std::all_of(&synds[0], &synds[RS::ecc], std::logical_not()))
            return;

        GFT err_poly[RS::ecc];
        auto errors = berlekamp_massey(synds, err_poly);

        GFT err_pos[RS::ecc / 2];
        RS::roots(&err_poly[RS::ecc-errors-1], errors+1, err_pos, size);

        GFT err_mag[RS::ecc / 2];
        forney(synds, &err_poly[RS::ecc-errors-1], err_pos, errors, err_mag);

        for (unsigned i = 0; i < errors; ++i)
            data[size - 1 - err_pos[i]] ^= err_mag[i];
    }

    static inline unsigned berlekamp_massey(const GFT synds[RS::ecc], GFT err_poly[RS::ecc]) {
        GFT prev[RS::ecc];
        GFT temp[RS::ecc];

        for (unsigned i = 0; i < RS::ecc; ++i) {
            err_poly[i] = 0;
            prev[i] = 0;
        }
        prev[RS::ecc-1] = 1;
        err_poly[RS::ecc-1] = 1;

        unsigned errors = 0;
        unsigned m = 1;
        GFT b = 1;

        for (unsigned n = 0; n < RS::ecc; ++n) {
            unsigned d = synds[n]; // discrepancy
            for (unsigned i = 1; i < errors + 1; ++i)
                d ^= RS::GF::mul(err_poly[RS::ecc - 1 - i], synds[n-i]);

            if (d == 0) {  // discrepancy is already zero
                m = m + 1;
            } else if (2 * errors <= n) {
                std::copy_n(err_poly, RS::ecc, temp);

                RS::GF::poly_shift(prev, RS::ecc, m);
                RS::GF::poly_scale(prev, RS::ecc, RS::GF::div(d, b));

                RS::GF::poly_add(err_poly, prev, RS::ecc);

                errors = n + 1 - errors;
                std::copy_n(temp, RS::ecc, prev);

                b = d;
                m = 1;
            } else {
                std::copy_n(prev, RS::ecc, temp);

                RS::GF::poly_shift(temp, RS::ecc, m);

                RS::GF::poly_scale(temp, RS::ecc, RS::GF::div(d, b));
                RS::GF::poly_add(err_poly, temp, RS::ecc);

                m = m + 1;
            }
        }

        return errors;
    }

    static inline void forney(
            const GFT synds[RS::ecc], GFT err_poly[], const GFT err_pos[],
            const unsigned err_count, GFT err_mag[])
    {
        GFT temp[RS::ecc] = {1};

        GFT err_eval[RS::ecc * 2];
        typename RS::synds_array_t synds_rev;
        std::reverse_copy(synds, &synds[RS::ecc], synds_rev);

        auto err_eval_size = RS::GF::poly_mul(
                err_eval, synds_rev,
                RS::ecc, err_poly, err_count + 1);

        auto err_eval_begin = RS::GF::ex_synth_div(
                err_eval, err_eval_size,
                temp, RS::ecc);

        while (err_eval[err_eval_begin] == 0) {
            err_eval_begin++;
            assert(err_eval_begin < err_eval_size);
        }
        err_eval_size = err_eval_size - err_eval_begin;

        auto err_poly_deriv_size = RS::GF::poly_deriv(err_poly, err_count + 1);

        for (unsigned i = 0; i < err_count; ++i) {
            auto xi = RS::GF::exp(err_pos[i]);
            auto xi_inv = RS::GF::inv(xi);

            auto n = RS::GF::poly_eval(&err_eval[err_eval_begin], err_eval_size, xi_inv);
            auto d = RS::GF::poly_eval(err_poly, err_poly_deriv_size, xi_inv);

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

