#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>

namespace detail {
    static inline constexpr unsigned ilog2_floor(unsigned a) {
        unsigned r = 0;
        while (a >>= 1)
            r += 1;
        return r;
    }

    static inline constexpr unsigned ipow(unsigned a, unsigned b) {
        unsigned r = 1;
        for (unsigned i = 0; i < b; ++i)
            r *= a;
        return r;
    }
}

template<typename T, T Prime, T Power, T Primitive, T Poly1>
struct gf_base {
    using Repr = T;

    static constexpr auto prime = Prime;
    static constexpr auto power = Power;
    static constexpr auto primitive = Primitive;
    static constexpr auto poly1 = Poly1;
    static constexpr auto charact = detail::ipow(prime, power);
};


template<typename GF>
struct gf_add_ring {
    using T = typename GF::Repr;
    static inline constexpr T add(T const& lhs, T const& rhs) { return (lhs + rhs) % GF::prime; }
    static inline constexpr T sub(T const& lhs, T const& rhs) {
        if (lhs >= rhs)
            return (lhs - rhs) % GF::prime;
        else
            return ((GF::charact - rhs) + lhs) % GF::prime;
    }
};


template<typename GF>
struct gf_add_xor {
    using GFT = typename GF::Repr;
    static inline constexpr GFT add(GFT const& lhs, GFT const& rhs) { return lhs ^ rhs; }
    static inline constexpr GFT sub(GFT const& lhs, GFT const& rhs) { return lhs ^ rhs; }
};

template<typename GF>
struct gf_mul_cpu {
    using GFT = typename GF::Repr;

    template<typename _GF = GF, std::enable_if_t<_GF::prime == 2, bool> = true>
    static inline constexpr GFT mul(GFT const& a, GFT const& b) {
        GFT r = 0;

        constexpr auto iter = detail::ilog2_floor(GF::charact >> 1);
        for (int i = iter; i >= 0; --i) {
            if (r & (GF::charact >> 1))
                r = (r << 1) ^ GF::poly1;
            else
                r = (r << 1);

            if (a & (1 << i))
                r ^= b;
        }

        return r;
    }

    template<typename _GF = GF, std::enable_if_t<_GF::power == 1, bool> = true>
    static inline constexpr GFT mul(GFT const& lhs, GFT const& rhs) { return (lhs * rhs) % GF::prime; }
};

template<typename GF>
struct gf_exp_log_lut {
    using GFT = typename GF::Repr;

    static inline constexpr struct sdata_t {
        std::array<GFT, GF::charact> exp{};
        std::array<GFT, GF::charact> log{};

        constexpr inline sdata_t() {
            GFT x = 1;
            for (unsigned i = 0; i < GF::charact; ++i) {
                exp[i] = x;
                log[x] = i;
                x = gf_mul_cpu<GF>::mul(x, GF::primitive);
            }
        }
    } sdata{};

    static inline constexpr GFT inv(GFT const& a) {
        return sdata.exp[GF::charact-1 - sdata.log[a]];
    }

    static inline constexpr GFT div(GFT const& a, GFT const& b) {
        if (a == 0)
            return 0;

        unsigned r = sdata.log[a] + GF::charact-1 - sdata.log[b];
        if (r >= GF::charact-1)
            r -= GF::charact-1;

        return sdata.exp[r];
    }

    static inline constexpr GFT exp(GFT const& a) {
        return sdata.exp[a];
    }

    static inline constexpr GFT log(GFT const& a) {
        return sdata.log[a];
    }

    static inline constexpr GFT pow(GFT const& a, GFT const& b) {
        return sdata.exp[(sdata.log[a] * b) % (GF::charact - 1)];
    }
};

template<typename GF>
struct gf_mul_lut {
    using GFT = typename GF::Repr;
    static_assert(GF::charact <= 4096); // limit 16MB

    static inline constexpr struct sdata_t {
        GFT mul[GF::charact][GF::charact] = {};

        constexpr inline sdata_t() {
            for (unsigned i = 0; i < GF::charact; ++i) {
                for (unsigned j = 0; j < GF::charact; ++j)
                    mul[i][j] = gf_mul_cpu<GF>::mul(i, j);
            }
        }
    } sdata{};

    static inline constexpr GFT mul(GFT const& a, GFT const& b) {
        return sdata.mul[a][b];
    }
};

template<typename GF>
struct gf_mul_exp_log_lut {
    using GFT = typename GF::Repr;

    static constexpr auto& sdata = gf_exp_log_lut<GF>::sdata;

    static inline constexpr GFT mul(GFT const& a, GFT const& b) {
        if (a == 0 || b == 0)
            return 0;

        unsigned r = sdata.log[a] + sdata.log[b];
        if (r >= GF::charact-1)
            r -= GF::charact-1;

        return sdata.exp[r];
    }
};

template<typename GF, typename Word>
class gf_wide_mul {
    static_assert(std::is_same_v<typename GF::Repr, uint8_t>);

private:
    static inline constexpr Word _w(uint8_t n) {
        Word p = 0;
        for (unsigned i = 0; i < sizeof(Word); ++i)
            p |= Word(n) << (i * 8);

        return p;
    }

public:
    static constexpr Word polyw = _w(GF::poly1);

    static inline constexpr Word mul(Word a, Word b) {
        Word r = 0;
        for (int i = 7; i >= 0; --i) {
            Word m = r & _w(0x80);

            static_assert((GF::poly1 & 0x80) == 0); // needed for the next statement
            m = m - (m >> 7);

            r = ((r & _w(0x7f)) << 1) ^ (polyw & m);

            Word n = (a & (_w(0x01) << i)) >> i;
            n = (n << 8) - n;

            r ^= b & n;
        }
        return r;
    }

    static inline constexpr Word poly_eval(const uint8_t poly[], const unsigned size, const Word x, Word r = 0) {
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
    using GFT = typename GF::Repr;

    template<typename T, typename U>
    static inline constexpr unsigned ex_synth_div(T a[], unsigned size_a, const U b[], unsigned size_b) {
        // T normalizer = b[0];
        // assert(b[0] == 1);

        if (size_b > size_a)
            return 0;

        for (unsigned i = 0; i < size_a - size_b + 1; ++i) {
            T c = a[i];// = div(a[i], normalizer);

            if (c == 0)
                continue;

            for (unsigned j = 1; j < size_b; ++j)
                a[i + j] = GF::sub(a[i + j], GF::mul(b[j], c));
        }

        return size_a - size_b + 1;
    }

    template<typename T>
    static inline constexpr unsigned poly_mod(T r[], const T a[], unsigned size_a, const T b[], unsigned size_b) {
        // assert(b[0] == 1);
        b += 1;

        if (size_b < 2)
            return 0;

        const unsigned size_r = size_b - 1;
        std::copy_n(a, std::min(size_a, size_r), r);

        if (size_a < size_b)
            return size_a;

        for (unsigned i = 0; i < size_a - size_r; ++i) {
            T c = r[0];
            std::rotate(r, r + 1, r + size_r);
            r[size_r - 1] = a[size_r + i];

            if (c == 0)
                continue;

            for (unsigned j = 0; j < size_r; ++j)
                r[j] = GF::sub(r[j], GF::mul(b[j], c));
        }

        return size_r;
    }

    template<typename T, typename U, typename V>
    static inline constexpr void poly_mod_x_n(T rem[],
            const U a[], const unsigned size_a,
            const V b[], const unsigned size_b)
    {
        if (size_a >= size_b) {
            std::copy_n(a, size_b, rem);
            for (unsigned i = 0; i < size_a - size_b; ++i) {
                T c = rem[0];
                rem[0] = a[size_b + i];
                std::rotate(rem, rem + 1, rem + size_b);

                if (c == 0)
                    continue;

                for (unsigned j = 0; j < size_b; ++j)
                    rem[j] = GF::sub(rem[j], GF::mul(b[j], c));
                // std::transform(&rem[0], &rem[size_b], &b[0], &rem[0],
                //         [c](T a, T b) { return GF::sub(a, GF::mul(b, c)); });
            }
        } else {
            std::fill_n(rem, size_b - size_a, 0);
            std::copy_n(a, size_a, &rem[size_b - size_a]);
        }

        for (unsigned i = 0; i < size_b; ++i) {
            T c = rem[0];
            rem[0] = 0;
            std::rotate(rem, rem + 1, rem + size_b);

            if (c == 0)
                continue;

            for (unsigned j = 0; j < size_b; ++j)
                rem[j] = GF::sub(rem[j], GF::mul(b[j], c));
        }
    }

    template<typename T>
    static inline constexpr GFT poly_eval(const T poly[], const unsigned size, GFT const& x, GFT r = 0) {
        for (unsigned i = 0; i < size; ++i)
            r = GF::add(GF::mul(r, x), poly[i]);
        return r;
    }

    template<typename T>
    static inline constexpr void poly_shift(T poly[], const unsigned size, const unsigned n) {
        unsigned i = 0;
        for (; i < size - n; ++i)
            poly[i] = poly[i + n];
        for (; i < size; ++i)
            poly[i] = 0;
    }

    template<typename T>
    static inline constexpr void poly_scale(T poly[], const unsigned size, T const& a) {
        for (unsigned i = 0; i < size; ++i)
            poly[i] = GF::mul(a, poly[i]);
    }

    template<typename T>
    static inline constexpr void poly_add(T poly_a[], const T poly_b[], const unsigned size) {
        for (unsigned i = 0; i < size; ++i)
            poly_a[i] = GF::add(poly_a[i], poly_b[i]);
    }

    template<typename T>
    static inline constexpr void poly_sub(T poly_a[], const T poly_b[], const unsigned size) {
        for (unsigned i = 0; i < size; ++i)
            poly_a[i] = GF::sub(poly_a[i], poly_b[i]);
    }

    template<typename T, typename U>
    static inline constexpr unsigned poly_mul(GFT r[],
            T poly_a[], const unsigned size_a,
            U poly_b[], const unsigned size_b) {
        unsigned res_len = size_a + size_b - 1;

        for (unsigned i = 0; i < res_len; ++i)
            r[i] = 0;

        for (unsigned i = 0; i < size_a; ++i)
            for (unsigned j = 0; j < size_b; ++j)
                r[i + j] = GF::add(r[i + j], GF::mul(poly_a[i], poly_b[j]));

        return res_len;
    }

    template<typename T, typename _GF = GF, std::enable_if_t<_GF::prime == 2, bool> = true>
    static inline constexpr unsigned poly_deriv(T poly[], unsigned size) {
        for (unsigned i = 1; i < size; ++i)
            poly[size - i] = (i & 1) ? poly[size - i - 1] : 0;

        return size - 1;
    }

    template<typename T, typename _GF = GF, std::enable_if_t<_GF::power == 1, bool> = true>
    static inline constexpr unsigned poly_deriv(T poly[], unsigned size) {
        for (unsigned i = 1; i < size; ++i)
            poly[size - i] = GF::mul(poly[size - i - 1], i);
        poly[0] = 0;
        return size - 1;
    }
};

template<typename T, T Prime, T Power, T Primitive, T Poly1, template<class>typename...Fs>
struct gf_impl_base :
        gf_base<T, Prime, Power, Primitive, Poly1>,
        Fs<gf_base<T, Prime, Power, Primitive, Poly1>>... { };

template<typename Impl, template<class>typename...Fs>
struct gf_impl : gf_base<typename Impl::Repr, Impl::prime, Impl::power, Impl::primitive, Impl::poly1>, Fs<Impl>... {
    static inline auto static_data_size = detail::get_sdata_size<Fs<Impl>...>();
    typename Impl::Repr value;

    explicit inline constexpr gf_impl(): value(typename Impl::Repr()) { }
    explicit inline constexpr gf_impl(typename Impl::Repr const& v): value(v) { }

    friend inline constexpr gf_impl operator*(gf_impl const& lhs, gf_impl const& rhs) {
        return gf_impl{gf_impl::mul(lhs.value, rhs.value)};
    }

    friend inline constexpr gf_impl operator+(gf_impl const& lhs, gf_impl const& rhs) {
        return gf_impl{gf_impl::add(lhs.value, rhs.value)};
    }

    friend inline constexpr gf_impl operator-(gf_impl const& lhs, gf_impl const& rhs) {
        return gf_impl{gf_impl::sub(lhs.value, rhs.value)};
    }

    friend constexpr inline bool operator==(gf_impl const& lhs, gf_impl const& rhs) {
        return lhs.value == rhs.value;
    }
};

template<typename T, T Prime, T Power, T Primitive, T Poly1, template<class>typename...Fs>
struct GF :
        gf_impl<gf_impl_base<T, Prime, Power, Primitive, Poly1, Fs...>, Fs...>,
        gf_poly<gf_impl<gf_impl_base<T, Prime, Power, Primitive, Poly1, Fs...>, Fs...>>
{
    using gf_impl<gf_impl_base<T, Prime, Power, Primitive, Poly1, Fs...>, Fs...>::gf_impl;
};
