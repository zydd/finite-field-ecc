#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>

template<uint8_t Primitive, uint16_t Poly1>
struct gf_base {
    static constexpr auto primitive = Primitive;
    static constexpr auto poly1 = Poly1;
};

template<typename GF>
struct gf_mul_cpu {
    static inline constexpr uint8_t mul(uint8_t a, uint8_t b) {
        uint8_t r = 0;
        for (int i = 7; i >= 0; --i) {
            if (r & 0x80)
                r = (r << 1) ^ GF::poly1;
            else
                r = (r << 1);

            if (a & (1 << i))
                r ^= b;
        }

        return r;
    }
};

template<typename GF>
struct gf_exp_log_lut {
    static inline constexpr struct sdata_t {
        std::array<uint8_t, 256> exp{};
        std::array<uint8_t, 256> log{};

        constexpr inline sdata_t() {
            uint8_t x = 1;
            for (unsigned i = 0; i < 256; ++i) {
                exp[i] = x;
                log[x] = i;
                x = gf_mul_cpu<GF>::mul(x, GF::primitive);
            }
        }
    } sdata{};

    static inline constexpr uint8_t inv(uint8_t a) {
        return sdata.exp[255 - sdata.log[a]];
    }

    static inline constexpr uint8_t div(uint8_t a, uint8_t b) {
        if (a == 0)
            return 0;

        unsigned r = sdata.log[a] + 255 - sdata.log[b];
        if (r >= 255)
            r -= 255;

        return sdata.exp[r];
    }

    static inline constexpr uint8_t exp(uint8_t a) {
        return sdata.exp[a];
    }

    static inline constexpr uint8_t log(uint8_t a) {
        return sdata.log[a];
    }
};

template<typename GF>
struct gf_mul_lut {
    static inline constexpr struct sdata_t {
        std::array<std::array<uint8_t, 256>, 256> mul{};

        constexpr inline sdata_t() {
            for (unsigned i = 0; i < 256; ++i) {
                for (unsigned j = 0; j < 256; ++j)
                    mul[i][j] = gf_mul_cpu<GF>::mul(i, j);
            }
        }
    } sdata{};

    static inline constexpr uint8_t mul(uint8_t a, uint8_t b) {
        return sdata.mul[a][b];
    }
};

template<typename GF>
struct gf_mul_exp_log_lut {
    static constexpr auto& sdata = gf_exp_log_lut<GF>::sdata;

    static inline constexpr uint8_t mul(uint8_t a, uint8_t b) {
        if (a == 0 || b == 0)
            return 0;

        unsigned r = sdata.log[a] + sdata.log[b];
        if (r >= 255)
            r -= 255;

        return sdata.exp[r];
    }
};

template<typename GF, typename T>
class gf_wide_mul {
private:
    static inline constexpr T _w(uint8_t n) {
        T p = 0;
        for (unsigned i = 0; i < sizeof(T); ++i)
            p |= T(n) << (i * 8);

        return p;
    }

public:
    static constexpr T polyw = _w(GF::poly1);

    static inline constexpr T mul(T a, T b) {
        T r = 0;
        for (int i = 7; i >= 0; --i) {
            T m = r & _w(0x80);

            static_assert((GF::poly1 & 0x80) == 0); // needed for the next statement
            m = m - (m >> 7);

            r = ((r & _w(0x7f)) << 1) ^ (polyw & m);

            T n = (a & (_w(0x01) << i)) >> i;
            n = (n << 8) - n;

            r ^= b & n;
        }
        return r;
    }

    static inline constexpr T poly_eval(const uint8_t poly[], const unsigned size, const T x) {
        T r = 0;
        for (unsigned i = 0; i < size; ++i)
            r = mul(r, x) ^ (poly[i] * _w(0x01));
        return r;
    }
};

namespace detail {
    template<typename T, typename E = void>
    struct get_sdata_size_helper { static const auto value = 0; };
    template<typename T>
    struct get_sdata_size_helper<T, std::enable_if_t<sizeof(typename T::sdata_t)>> {
        static const auto value = sizeof(typename T::sdata_t); };

    template<typename T, typename...Ts>
    auto get_sdata_size() {
        if constexpr (sizeof...(Ts) > 0)
            return get_sdata_size_helper<T>::value + get_sdata_size<Ts...>();
        else
            return get_sdata_size_helper<T>::value;
    }
}

template<typename GF>
struct gf_poly {
    static inline constexpr unsigned ex_synth_div(uint8_t a[], unsigned size_a, const uint8_t b[], unsigned size_b) {
        // uint8_t normalizer = b[0];
        // assert(b[0] == 1);

        if (size_b > size_a)
            return 0;

        for (unsigned i = 0; i < size_a - size_b + 1; ++i) {
            uint8_t coef = a[i];// = div(a[i], normalizer);

            if (coef == 0)
                continue;

            for (unsigned j = 1; j < size_b; ++j)
                a[i + j] ^= GF::mul(b[j], coef);
        }

        return size_a - size_b + 1;
    }

    static inline constexpr unsigned poly_mod(uint8_t r[], const uint8_t a[], unsigned size_a, const uint8_t b[], unsigned size_b) {
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
                r[j] ^= GF::mul(b[j], c);
        }

        return size_r;
    }

    static inline constexpr void poly_mod_x_n(uint8_t rem[], const uint8_t a[], const unsigned size_a, const uint8_t b[], const unsigned size_b) {
        if (size_a >= size_b) {
            std::copy_n(a, size_b, rem);
            for (unsigned i = 0; i < size_a - size_b; ++i) {
                uint8_t c = rem[0];
                rem[0] = a[size_b + i];
                std::rotate(rem, rem + 1, rem + size_b);

                if (c == 0)
                    continue;

                for (unsigned j = 0; j < size_b; ++j)
                    rem[j] ^= GF::mul(b[j], c);
                // std::transform(&rem[0], &rem[size_b], &b[0], &rem[0],
                //         [c](uint8_t a, uint8_t b) { return a ^ GF::mul(b, c); });
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
                rem[j] ^= GF::mul(b[j], c);
        }
    }

    static inline constexpr uint8_t poly_eval(const uint8_t poly[], const unsigned size, const uint8_t x) {
        uint8_t r = 0;
        for (unsigned i = 0; i < size; ++i)
            r = GF::mul(r, x) ^ poly[i];
        return r;
    }

    static inline constexpr void poly_shift(uint8_t poly[], const unsigned size, const unsigned n) {
        unsigned i = 0;
        for (; i < size - n; ++i)
            poly[i] = poly[i + n];
        for (; i < size; ++i)
            poly[i] = 0;
    }

    static inline constexpr void poly_scale(uint8_t poly[], const unsigned size, const uint8_t a) {
        for (unsigned i = 0; i < size; ++i)
            poly[i] = GF::mul(a, poly[i]);
    }

    static inline constexpr void poly_add(uint8_t poly_a[], const uint8_t poly_b[], const unsigned size) {
        for (unsigned i = 0; i < size; ++i)
            poly_a[i] ^= poly_b[i];
    }

    static inline constexpr unsigned poly_mul(uint8_t r[],
            const uint8_t poly_a[], const unsigned size_a,
            const uint8_t poly_b[], const unsigned size_b) {
        unsigned res_len = size_a + size_b - 1;

        for (unsigned i = 0; i < res_len; ++i)
            r[i] = 0;

        for (unsigned i = 0; i < size_a; ++i)
            for (unsigned j = 0; j < size_b; ++j)
                r[i + j] ^= GF::mul(poly_a[i], poly_b[j]);

        return res_len;
    }

    static inline constexpr unsigned poly_deriv(uint8_t poly[], unsigned size) {
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

template<uint8_t Primitive, uint16_t Poly1, template<class>typename...Fs>
struct gf_impl_base : gf_base<Primitive, Poly1>, Fs<gf_base<Primitive, Poly1>>... { };

template<typename Impl, template<class>typename...Fs>
struct gf_impl : gf_base<Impl::primitive, Impl::poly1>, Fs<Impl>... {
    static inline auto static_data_size = detail::get_sdata_size<Fs<Impl>...>();
};

template<uint8_t Primitive, uint16_t Poly1, template<class>typename...Fs>
struct GF :
        gf_impl<gf_impl_base<Primitive, Poly1, Fs...>, Fs...>,
        gf_poly<gf_impl<gf_impl_base<Primitive, Poly1, Fs...>, Fs...>>
{ };
