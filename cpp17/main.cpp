#include <chrono>
#include <cstdio>
#include <iostream>
#include <random>
#include <type_traits>

#include "reed_solomon.hpp"

std::mt19937_64 mersenne;
using hrc = std::chrono::high_resolution_clock;

template<typename T>
void fill_random(T buf[], unsigned size) {
    for (unsigned i = 0; i < size; ++i)
        buf[i] = mersenne();
}

void benchmark_enc_dec() {
    const unsigned ecclen = 8;
    const unsigned msglen = 40-ecclen;

    uint8_t buffer[msglen + ecclen];

    using GF = ::GF<uint8_t, 2, 8, 2, 0x11d & 0xff, gf_add_xor, gf_exp_log_lut, gf_mul_exp_log_lut>;
    using RS = ::RS<GF, ecclen, rs_encode_basic, rs_synds_lut8, rs_roots_eval_basic, rs_decode>;

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
            RS::encode(buffer, msglen + ecclen);
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

int main(int argc, char *argv[]) {
    using GF = GF<uint8_t, 2, 8, 2, 0x11d & 0xff, gf_mul_cpu>;
    static_assert(GF(GF::mul(23, 47)) == GF(23) * GF(47));

    mersenne.seed(42);
    benchmark_enc_dec();


    return 0;
}
