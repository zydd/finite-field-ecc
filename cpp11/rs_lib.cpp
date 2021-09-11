#define RS_GENERATOR_LUT
#include "reed_solomon.hpp"

using RS0 = RS<4, 2, 0x11d, uint32_t>;

extern "C" {

void *gf_init() {
    return new RS0;
}

void gf_uninit(void *rs) {
    delete reinterpret_cast<RS0 *>(rs);
}

uint8_t gf_mul(void *rs, uint8_t a, uint8_t b) {
    return reinterpret_cast<RS0 *>(rs)->gf.mul(a, b);
}

uint8_t _mul(void *rs, uint8_t a, uint8_t b) {
    return reinterpret_cast<RS0 *>(rs)->gf._mul(a, b);
}

uint8_t gf_inv(void *rs, uint8_t a) {
    return reinterpret_cast<RS0 *>(rs)->gf.inv(a);
}

uint8_t gf_div(void *rs, uint8_t a, uint8_t b) {
    return reinterpret_cast<RS0 *>(rs)->gf.div(a, b);
}

unsigned ex_synth_div(void *rs, uint8_t a[], unsigned size_a, const uint8_t b[], unsigned size_b) {
    return reinterpret_cast<RS0 *>(rs)->gf.ex_synth_div(a, size_a, b, size_b);
}

unsigned gf_poly_mod(void *rs, uint8_t r[], const uint8_t a[], unsigned size_a, const uint8_t b[], unsigned size_b) {
    return reinterpret_cast<RS0 *>(rs)->gf.poly_mod(r, a, size_a, b, size_b);
}

unsigned gf_poly_mul(void *rs, uint8_t r[], const uint8_t a[], unsigned size_a, const uint8_t b[], unsigned size_b) {
    return reinterpret_cast<RS0 *>(rs)->gf.poly_mul(r, a, size_a, b, size_b);
}

uint8_t gf_poly_eval(void *rs, uint8_t a[], unsigned size_a, const uint8_t x) {
    return reinterpret_cast<RS0 *>(rs)->gf.poly_eval(a, size_a, x);
}

uint32_t gf_poly_eval4(void *rs, uint8_t a[], unsigned size_a, const uint32_t x) {
    return reinterpret_cast<RS0 *>(rs)->gf.poly_eval(a, size_a, x);
}

void encode(void *rs, uint8_t a[], unsigned size) {
    return reinterpret_cast<RS0 *>(rs)->encode(a, size);
}

void decode(void *rs, uint8_t a[], unsigned size) {
    return reinterpret_cast<RS0 *>(rs)->decode(a, size);
}

}
