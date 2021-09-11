#!/usr/bin/env python3

import ctypes
import os
import random
import time
import sys

sys.path.append(os.path.abspath('../python'))

import gf
import bch
import rs

def test(fn):
    def _wrapper(*args, **kwargs):
        print(f'{fn.__name__}')
        fn(*args, **kwargs)
    return _wrapper

class RSC:
    def __init__(self, power, prim, poly, ecc_len):
        self.c_lib = ctypes.CDLL('./lib.so')

        self.c_lib.gf_mul.restype       = ctypes.c_uint8
        self.c_lib._mul.restype         = ctypes.c_uint8
        self.c_lib.gf_mul4.restype      = ctypes.c_uint32
        self.c_lib.gf_inv.restype       = ctypes.c_uint8
        self.c_lib.gf_div.restype       = ctypes.c_uint8
        self.c_lib.gf_init.restype      = ctypes.c_void_p
        self.c_lib.gf_mul16.restype     = ctypes.c_uint16
        self.c_lib.ex_synth_div.restype = ctypes.c_uint
        self.c_lib.gf_poly_mul.restype  = ctypes.c_uint
        self.c_lib.gf_poly_eval.restype = ctypes.c_uint8
        self.c_lib.gf_poly_eval4.restype= ctypes.c_uint32

        self.gf_ctx = self.c_lib.gf_init()

    def _mul(self, a, b):
        return self.c_lib._mul(self.gf_ctx, ctypes.c_uint8(a), ctypes.c_uint8(b))

    def gf_mul(self, a, b):
        return self.c_lib.gf_mul(self.gf_ctx, ctypes.c_uint8(a), ctypes.c_uint8(b))

    def gf_mul16(self, a, b):
        return self.c_lib.gf_mul16(self.gf_ctx, ctypes.c_uint16(a), ctypes.c_uint16(b))

    def gf_mul4(self, a, b):
        return self.c_lib.gf_mul4(self.gf_ctx, ctypes.c_uint32(a), ctypes.c_uint32(b))

    def gf_inv(self, a):
        return self.c_lib.gf_inv(self.gf_ctx, ctypes.c_uint8(a))

    def gf_div(self, a, b):
        return self.c_lib.gf_div(self.gf_ctx, ctypes.c_uint8(a), ctypes.c_uint8(b))

    def ex_synth_div(self, a, b):
        res = (ctypes.c_uint8 * len(a))(*a)
        dividend = (ctypes.c_uint8 * len(b))(*b)
        sep = self.c_lib.ex_synth_div(self.gf_ctx, res, len(a), dividend, len(b))
        res = list(res)
        return res[:sep], res[sep:]

    def poly_mod(self, a, b):
        a = (ctypes.c_uint8 * len(a))(*a)
        b = (ctypes.c_uint8 * len(b))(*b)
        r = (ctypes.c_uint8 * len(b))()
        s = self.c_lib.gf_poly_mod(self.gf_ctx, r, a, len(a), b, len(b))
        return list(r[:s])

    def poly_eval(self, a, x):
        poly = (ctypes.c_uint8 * len(a))(*a)
        return self.c_lib.gf_poly_eval(self.gf_ctx, poly, len(a), ctypes.c_uint8(x))

    def poly_eval4(self, a, x):
        poly = (ctypes.c_uint8 * len(a))(*a)
        x = int.from_bytes(bytes(x), 'big')
        r = self.c_lib.gf_poly_eval4(self.gf_ctx, poly, len(a), ctypes.c_uint32(x))
        return list(int(r).to_bytes(4, 'big'))

    def poly_mul(self, a, b):
        pa = (ctypes.c_uint8 * len(a))(*a)
        pb = (ctypes.c_uint8 * len(b))(*b)
        res = (ctypes.c_uint8 * (len(a) + len(b) - 1))()
        self.c_lib.gf_poly_mul(self.gf_ctx, res, pa, len(a), pb, len(b))
        return list(res)

    def encode(self, a):
        res = (ctypes.c_uint8 * len(a))(*a)
        self.c_lib.encode(self.gf_ctx, res, len(a))
        return list(res)

    def decode(self, a):
        res = (ctypes.c_uint8 * len(a))(*a)
        self.c_lib.decode(self.gf_ctx, res, len(a))
        return list(res)

if os.system('g++ -O3 -std=c++17 -Wall -shared -fPIC ./lib.cpp -o lib.so') != 0:
    quit()

ecc_len = 4
GF = rs.GF
poly = GF.poly_to_int(GF.p, GF.poly)
print(f'init field: power: {GF.k} prim: {GF.a} poly: 0x{poly:x}')
RS = RSC(GF.k, GF.a, poly, ecc_len)

GF64k = gf.GF(2, 16, *next(GF.irr_polynomials(2, 16, 2)))
poly = GF64k.poly_to_int(GF64k.p, GF64k.poly)
print(f'init field: power: {GF64k.k} prim: {GF64k.a} poly: 0x{poly:x}')

def assert_eq(a, b, l, r):
    assert l == r, (int(a), int(b), int(l), int(r))

@test
def test_mul():
    for a in range(GF.p ** GF.k):
        for b in range(GF.p ** GF.k):
            assert_eq(a, b, int(GF(a) * GF(b)), RS._mul(a, b))
        print(f'{a}/255    ', end='\r')

@test
def test_gf_mul():
    for a in range(GF.p ** GF.k):
        for b in range(GF.p ** GF.k):
            assert_eq(a, b, int(GF(a) * GF(b)), RS.gf_mul(a, b))
        print(f'{a}/255    ', end='\r')

@test
def test_gf_mul4():
    for a in range(GF.p ** GF.k):
        for b in range(GF.p ** GF.k):
            # d = [int(GF(b + i)) for i in range(4)]
            # ref = int.from_bytes(bytes([int(GF(a) * GF(i)) for i in d]), 'big')
            # print(a, d, ref)
            # d = int.from_bytes(bytes(d), 'big')
            # print(d, RS.gf_mul4(a, d))

            # assert_eq(a, b, ref, RS.gf_mul4(a, d))
            assert_eq(a, b, int(GF(a) * GF(b)), RS.gf_mul4(a, b))
        print(f'{a}/255    ', end='\r')

@test
def test_gf_inv():
    for a in range(1, GF.p ** GF.k):
        assert_eq(a, -1, int(GF(a).inv()), RS.gf_inv(a))

@test
def test_gf_div():
    for a in range(GF.p ** GF.k):
        for b in range(1, GF.p ** GF.k):
            assert_eq(a, b, int(GF(a) // GF(b)), RS.gf_div(a, b))
        print(f'{a}/255    ', end='\r')

@test
def test_gf_mul16():
    for i in range(50000):
        a = random.randrange(GF64k.p ** GF64k.k)
        b = random.randrange(GF64k.p ** GF64k.k)
        assert_eq(a, b, int(GF64k(a) * GF64k(b)), RS.gf_mul16(a, b))
        if i & 0xff == 0:
            print(f'{i}/100000    ', end='\r')

@test
def test_ex_synth_div():
    for len_a in range(32):
        for len_b in range(1, len_a + 2):
            a = [random.randrange(GF.p ** GF.k) for _ in range(len_a)]
            b = [random.randrange(GF.p ** GF.k) for _ in range(len_b)]
            b[0] = 1

            q, r = divmod(gf.P(GF, reversed(a)), gf.P(GF, reversed(b)))
            q, r = list(map(int, reversed(q.x))), list(map(int, reversed(r.x)))
            gfqr = [0] * (len(a) - len(b) - len(q) + 1) + q, [0] * (len(b) - 1 - len(r)) + r

            rsqr = RS.ex_synth_div(a, b)

            if gfqr != rsqr:
                print(f'a: {a} b: {b}')
                print(f'gf: {gfqr}')
                print(f'rs: {rsqr}')
                assert False
        print(f'{len_a}/31    ', end='\r')

@test
def test_poly_mod():
    for len_a in range(100):
        for len_b in range(1, len_a + 2):
            a = [random.randrange(GF.p ** GF.k) for _ in range(len_a)]
            b = [random.randrange(GF.p ** GF.k) for _ in range(len_b)]
            b[0] = 1

            ex = RS.ex_synth_div(a, b)
            # print(a,b,ex)
            mod = RS.poly_mod(a, b)

            assert ex[1] == mod, (ex[1], mod)

        print(f'{len_a}/99    ', end='\r')

@test
def test_poly_eval():
    for i in range(32):
        for x in range(GF.p ** GF.k):
            a = [random.randrange(GF.p ** GF.k) for _ in range(i)]

            gfpoly = gf.P(GF, reversed(a))

            assert int(gfpoly.eval(GF(x))) == RS.poly_eval(a, x)
        print(f'{i}/31    ', end='\r')

@test
def test_poly_eval4():
    for i in range(32):
        a = [random.randrange(GF.p ** GF.k) for _ in range(i)]
        x = [random.randrange(GF.p ** GF.k) for _ in range(4)]

        r1 = [RS.poly_eval(a, i) for i in x]
        r2 = RS.poly_eval4(a, x)

        assert r1 == r2

@test
def test_poly_mul():
    for i in range(20):
        for j in range(20):
            a = [random.randrange(1, GF.p ** GF.k)] + [random.randrange(GF.p ** GF.k) for _ in range(i)]
            b = [random.randrange(1, GF.p ** GF.k)] + [random.randrange(GF.p ** GF.k) for _ in range(j)]

            ref = gf.P(GF, reversed(a)) * gf.P(GF, reversed(b))

            assert list(map(int, ref.x[::-1])) == RS.poly_mul(a,b)

@test
def test_encode():
    gen = rs.rs_generator(ecc_len)
    for _ in range(10):
        a = [random.randrange(GF.p ** GF.k) for _ in range(2)]

        enc = RS.encode(a + [0] * ecc_len)

        ref = rs.rs_encode_systematic(a[::-1], gen)
        ref = [0] * (len(enc) - len(ref.x)) + list(map(int, ref[::-1]))

        if enc != ref:
            print(f'enc: {enc}')
            print(f'ref: {ref}')
            assert False

@test
def test_decode():
    for _ in range(10000):
        a = [random.randrange(GF.p ** GF.k) for _ in range(12)]

        enc = RS.encode(a + [0] * ecc_len)

        # ref = rs.rs_encode_systematic(a[::-1], rs.rs_generator(ecc_len))
        # ref = [0] * (len(enc) - len(ref.x)) + list(map(int, ref[::-1]))

        # assert enc == ref, ref

        for i in range(ecc_len//2):
            e1 = random.randrange(len(enc))
            enc[e1] ^= random.randrange(1, 256)

        dec = RS.decode(enc)

        # assert ref[:len(a)] == a, (ref[:len(a)], a)
        if dec[:len(a)] != a:
            ref = rs.decode(enc, ecc_len)
            print(f'data: {a}')
            print(f'dec:  {dec}')
            print(f'ref:  {ref}')
            assert False

if __name__ == '__main__':
    random.seed(42)
    test_mul()
    test_gf_mul()
    test_gf_mul4()
    test_gf_mul16()
    test_gf_inv()
    test_gf_div()
    test_ex_synth_div()
    test_poly_mod()
    test_poly_eval()
    test_poly_eval4()
    test_poly_mul()
    test_encode()
    test_decode()
