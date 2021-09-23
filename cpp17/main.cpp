#include <chrono>
#include <cstdio>
#include <iostream>
#include <random>
#include <type_traits>

#include "reed_solomon.hpp"

std::mt19937_64 mersenne;
using hrc = std::chrono::high_resolution_clock;

template<typename T, auto M = 256>
void fill_random(T buf[], unsigned size) {
    for (unsigned i = 0; i < size; ++i)
        buf[i] = mersenne() % M;
}

void benchmark_enc_dec() {
    std::cout << "benchmark_enc_dec" << std::endl;

    const unsigned ecclen = 8;
    const unsigned msglen = 40-ecclen;

    uint8_t buffer[msglen + ecclen];

    using GF = ::GF<uint16_t, 2, 8, 2, 0x11d & 0x1ff, gf_add_xor, gf_exp_log_lut, gf_mul_exp_log_lut>;
    using RS = ::RS<GF, ecclen, rs_encode_basic, rs_synds_basic, rs_roots_eval_basic, rs_decode>;

    std::cout << "sizeof(RS<" << ecclen << ">) = " << sizeof(RS) << std::endl;
    std::cout << "GF::static_data_size: " << GF::static_data_size << std::endl;
    std::cout << "RS::static_data_size: " << RS::static_data_size << std::endl;

    size_t processed = 0;
    auto t_enc = hrc::duration(0);
    auto t_dec = hrc::duration(0);

    auto start = hrc::now();
    while (hrc::now() - start < std::chrono::seconds(2)) {
        fill_random(buffer, msglen + ecclen);

        {
            auto start = hrc::now();
            RS::encode(buffer + msglen, buffer, msglen);
            t_enc += hrc::now() - start;
        }

        uint8_t buffer2[msglen + ecclen];
        std::copy_n(buffer, msglen + ecclen, buffer2);

        // errors
        for (unsigned i = 0; i < ecclen/2; ++i)
            buffer[mersenne() % (msglen + ecclen)] ^= mersenne();

        {
            auto start = hrc::now();
            RS::decode(buffer, msglen + ecclen);
            t_dec += hrc::now() - start;
        }

        assert(std::equal(buffer, buffer + msglen + ecclen, buffer2));
        processed += msglen;
    }

    auto enc_tp = processed / (t_enc.count() / 1e9);
    auto dec_tp = processed / (t_dec.count() / 1e9);

    std::cout << "encode: " << (enc_tp / 1e6) << " MB/s" << std::endl;
    std::cout << "decode: " << (dec_tp / 1e6) << " MB/s" << std::endl;
}

void benchmark_enc_257() {
    std::cout << "benchmark_enc_257" << std::endl;

    const unsigned ecclen = 4;
    const unsigned msglen = 16;

    uint16_t buffer[msglen + ecclen];

    using GF = ::GF<uint16_t, 257, 1, 3, 0, gf_add_ring, gf_mul_cpu, gf_exp_log_lut>;
    using RS = ::RS<GF, ecclen, rs_encode_basic, rs_synds_basic, rs_roots_eval_basic, rs_decode>;

    std::cout << "sizeof(RS<" << ecclen << ">) = " << sizeof(RS) << std::endl;
    std::cout << "GF::static_data_size: " << GF::static_data_size << std::endl;
    std::cout << "RS::static_data_size: " << RS::static_data_size << std::endl;

    size_t processed = 0;
    auto t_enc = hrc::duration(0);
    auto t_dec = hrc::duration(0);

    auto start = hrc::now();
    while (hrc::now() - start < std::chrono::seconds(2)) {
        std::fill(buffer, buffer + msglen + ecclen, 0);
        fill_random(buffer, msglen);
        std::fill(buffer, buffer + msglen, 5);

        {
            auto start = hrc::now();
            for (unsigned i = 0; i < msglen + ecclen; ++i)
                std::printf("%02x ", buffer[i]);
            std::printf("\n");
            RS::encode(buffer + msglen, buffer, msglen);
            auto &gen = rs_generator<RS>::sdata.generator;
            for (unsigned i = 0; i < gen.size(); ++i)
                std::printf("%d ", gen[i]);
            std::printf("\n");
            t_enc += hrc::now() - start;
        }

        uint16_t buffer2[msglen + ecclen];
        std::copy_n(buffer, msglen + ecclen, buffer2);

        // errors
        for (unsigned i = 0; i < ecclen/2; ++i)
            buffer[mersenne() % (msglen + ecclen)] ^= mersenne();

        {
            auto start = hrc::now();
            RS::decode(buffer, msglen + ecclen);
            t_dec += hrc::now() - start;
        }

        assert(std::equal(buffer, buffer + msglen + ecclen, buffer2));
        processed += msglen;
    }

    auto enc_tp = processed / (t_enc.count() / 1e9);
    auto dec_tp = processed / (t_dec.count() / 1e9);

    std::cout << "encode: " << (enc_tp / 1e6) << " MB/s" << std::endl;
    std::cout << "decode: " << (dec_tp / 1e6) << " MB/s" << std::endl;
}


void test_bit_array() {
    constexpr auto arr_size = 10000;
    constexpr auto bits = 30;
    using array_t = uint32_t;

    uint8_t arr[arr_size] = {};
    detail::bit_array<array_t, bits> wrap{arr};

    mersenne.seed(42);
    for (unsigned i = 0; i < arr_size * 8 / bits; ++i) {
        wrap[i] = mersenne() & ((1 << bits) - 1);
    }

    mersenne.seed(42);
    for (unsigned i = 0; i < arr_size * 8 / bits; ++i) {
        array_t expected = mersenne() & ((1 << bits) - 1);
        // std::printf("i %d: %04x %04x\n", i, array_t(wrap[i]), expected);

        if (array_t(wrap[i]) != expected) {
            std::printf("i %d: %04x %04x\n", i, array_t(wrap[i]), expected);
            return;
        }
    }
}

int main(int argc, char *argv[]) {
    using GF = ::GF<uint8_t, 2, 8, 2, 0x11d & 0xff, gf_mul_cpu>;
    static_assert(GF(GF::mul(23, 47)) == GF(23) * GF(47));

    test_bit_array();

    mersenne.seed(6942);
    benchmark_enc_dec();
    benchmark_enc_257();

    return 0;
}
