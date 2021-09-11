#define RS_GENERATOR_LUT
#include "reed_solomon.hpp"

static const auto ecclen = 4;
using GF256 = GF<uint8_t, 2, 8, 2, 0x11d & 0xff, gf_add_xor, gf_exp_log_lut, gf_mul_exp_log_lut>;
using RS0 = RS<GF256, ecclen, rs_encode_basic, rs_synds_lut8, rs_roots_eval_basic, rs_decode>;

extern "C" {

void *gf_init() {
    return new RS0;
}

void gf_uninit(void *rs) {
    delete reinterpret_cast<RS0 *>(rs);
}

uint8_t gf_mul(void *rs, uint8_t a, uint8_t b) {
    return RS0::GF::mul(a, b);
}

uint8_t _mul(void *rs, uint8_t a, uint8_t b) {
    return GF<uint8_t, 2, 8, 2, 0x11d & 0xff, gf_mul_cpu>::mul(a, b);
}

uint16_t gf_mul16(void *rs, uint16_t a, uint16_t b) {
    return GF<uint16_t, 2, 16, 2, 0x1002d & 0xffff, gf_mul_cpu>::mul(a, b);
}

uint32_t gf_mul4(void *rs, uint32_t a, uint32_t b) {
    return gf_wide_mul<RS0::GF, uint32_t>::mul(a, b);
}

uint8_t gf_inv(void *rs, uint8_t a) {
    return RS0::GF::inv(a);
}

uint8_t gf_div(void *rs, uint8_t a, uint8_t b) {
    return RS0::GF::div(a, b);
}

unsigned ex_synth_div(void *rs, uint8_t a[], unsigned size_a, const uint8_t b[], unsigned size_b) {
    return RS0::GF::ex_synth_div(a, size_a, b, size_b);
}

unsigned gf_poly_mod(void *rs, uint8_t r[], const uint8_t a[], unsigned size_a, const uint8_t b[], unsigned size_b) {
    return RS0::GF::poly_mod(r, a, size_a, b, size_b);
}

unsigned gf_poly_mul(void *rs, uint8_t r[], const uint8_t a[], unsigned size_a, const uint8_t b[], unsigned size_b) {
    return RS0::GF::poly_mul(r, a, size_a, b, size_b);
}

uint8_t gf_poly_eval(void *rs, uint8_t a[], unsigned size_a, const uint8_t x) {
    return RS0::GF::poly_eval(a, size_a, x);
}

uint32_t gf_poly_eval4(void *rs, uint8_t a[], unsigned size_a, const uint32_t x) {
    return gf_wide_mul<RS0::GF, uint32_t>::poly_eval(a, size_a, x);
}

void encode(void *rs, uint8_t a[], unsigned size) {
    return reinterpret_cast<RS0 *>(rs)->encode(a, size);
}

void decode(void *rs, uint8_t a[], unsigned size) {
    return reinterpret_cast<RS0 *>(rs)->decode(a, size);
}

}
