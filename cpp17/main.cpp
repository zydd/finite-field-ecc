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
    const unsigned msglen = 255-ecclen;
    const unsigned errors = 4;

    uint8_t buffer[msglen + ecclen];
    uint8_t err_pos[ecclen];

    using GF = ::GF<uint8_t, 2, 8, 2, uint8_t(0x11d), gf_add_xor, gf_exp_log_lut, gf_mul_exp_log_lut>;
    using RS = ::RS<GF, ecclen, rs_encode_slice<uint64_t, 16>::type, rs_synds_lut8, rs_roots_eval_chien, rs_decode>;

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
        for (unsigned i = 0; i < errors; ++i) {
            unsigned pos;
            do
                pos = mersenne() % (msglen + ecclen);
            while (std::find(err_pos, err_pos + i, pos) != (err_pos + i));

            err_pos[i] = pos;
            buffer[pos] ^= mersenne();
        }

        {
            auto start = hrc::now();
            RS::decode(&buffer[0], msglen, &buffer[msglen]);
            // RS::decode(&buffer[0], msglen, &buffer[msglen], err_pos, errors);
            t_dec += hrc::now() - start;
        }

        if (! std::equal(buffer, buffer + msglen + ecclen, buffer2)) {
            std::printf("exp: ");
            for (unsigned i = 0; i < msglen + ecclen; ++i)
                std::printf("%02x ", buffer2[i]);
            std::printf("\n");
            std::printf("got: ");
            for (unsigned i = 0; i < msglen + ecclen; ++i)
                std::printf("%02x ", buffer[i]);
            std::printf("\n");
            assert(false);
        }
        processed += msglen;
    }

    auto enc_tp = processed / (t_enc.count() / 1e9);
    auto dec_tp = processed / (t_dec.count() / 1e9);

    std::cout << "encode: " << (enc_tp / 1e6) << " MB/s" << std::endl;
    std::cout << "decode: " << (dec_tp / 1e6) << " MB/s" << std::endl;
}

void benchmark_enc_257() {
    std::cout << "benchmark_enc_257" << std::endl;

    const unsigned ecclen = 8;
    const unsigned msglen = 256-ecclen;
    const unsigned errors = 8;

    uint16_t buffer[msglen + ecclen];
    unsigned err_pos[ecclen];

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
        fill_random(buffer, msglen + ecclen);

        {
            auto start = hrc::now();
            RS::encode(buffer + msglen, buffer, msglen);
            t_enc += hrc::now() - start;
        }

        uint16_t buffer2[msglen + ecclen];
        std::copy_n(buffer, msglen + ecclen, buffer2);


        // errors
        for (unsigned i = 0; i < errors; ++i) {
            unsigned pos;
            do
                pos = mersenne() % (msglen + ecclen);
            while (std::find(err_pos, err_pos + i, pos) != (err_pos + i));

            err_pos[i] = pos;
            buffer[pos] ^= mersenne();
        }

        {
            auto start = hrc::now();
            // RS::decode(&buffer[0], msglen, &buffer[msglen]);
            RS::decode(&buffer[0], msglen, &buffer[msglen], err_pos, errors);
            t_dec += hrc::now() - start;
        }

        if (! std::equal(buffer, buffer + msglen + ecclen, buffer2)) {
            std::printf("exp: ");
            for (unsigned i = 0; i < msglen + ecclen; ++i)
                std::printf("%02x ", buffer2[i]);
            std::printf("\n");
            std::printf("got: ");
            for (unsigned i = 0; i < msglen + ecclen; ++i)
                std::printf("%02x ", buffer[i]);
            std::printf("\n");
            assert(false);
        }
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

    // test_bit_array();

    mersenne.seed(42);
    benchmark_enc_dec();
    benchmark_enc_257();

    return 0;
}
