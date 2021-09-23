#!/usr/bin/env python3

import random
import sys

import gf
from gf import P, F, PrimeField
import bch

GF = PrimeField(257, next(PrimeField.primitives(257)))
print(GF)


def rs_generator(n):
    g = P(GF, [1])
    for i in range(n):
        g *= P(GF, [-GF.gen(i), 1])
    return g

def rs_encode(msg, genertor):
    poly = P(GF, msg) * genertor

    return poly #// P(GF, [poly.x[-1]])

def rs_encode_systematic(msg, generator):
    poly = P(GF, [0] * generator.deg() + msg)

    return poly - (poly % generator)

def rs_syndromes(poly, n):
    return [poly.eval(GF.gen(i)) for i in range(n)]

def berlekamp_massey(s):
    C = P(GF, [GF(1)]) # correction polynomial
    B = P(GF, [GF(1)]) # previous C

    L = 0
    m = 1
    b = GF(1)

    for n in range(len(s)):
        # discrepancy

        d = s[n]
        for i in range(1, L+1):
            d += C[i] * s[n-i]

        if int(d) == 0:
            # discrepancy is zero; annihilation continues
            m = m + 1
        elif 2 * L <= n:
            T = P(GF, C) # temporary copy of C

            B.x = [B.N(0)] * m + B.x
            C = C - B.scale(d // b)

            L = n + 1 - L
            B = T
            b = d
            m = 1
        else:
            T = P(GF, B)

            T.x = [T.N(0)] * m + T.x
            C = C - T.scale(d // b)

            m = m + 1

    return C

def forney(synds, err_poly, err_pos):
    # number or possible corrections
    # t = ecc_len // GF.k

    # assert len(synds) == 2 * t
    # print(f'synds: {list(map(int, synds))}')
    print(f'err_poly: {list(map(int, err_poly.x))}')

    # err_eval = (synds * err_poly) % P(GF, [0] * len(synds.x) + [1])
    err_eval = (P(GF, synds) * err_poly) % P(GF, [0] * len(synds) + [1])
    print(f'err_eval: {list(map(int, err_eval.x))}')

    # err_poly_deriv = list(err_poly.x)
    # for i in range(0, len(err_poly_deriv), 2):
    #     err_poly_deriv[i] = 0
    # err_poly_deriv = P(GF, err_poly_deriv) // P(GF, [0, 1])
    # print(f'err_poly_deriv: {list(map(int, err_poly_deriv.x))}')

    err_poly_deriv = err_poly.deriv()
    print(f'err_poly_deriv: {list(map(int, err_poly_deriv.x))}')

    X = [GF.gen(i) for i in err_pos]
    # print(f'X: {list(map(int, X))}')

    errors = []
    for pos, Xi in zip(err_pos, X):
        Xi_i = Xi.inv()

        y = err_eval.eval(Xi_i)
        d = err_poly_deriv.eval(Xi_i)

        # print(f'err_eval[{int(Xi_i)}]: {int(y)}')
        # print(f'd[{int(Xi_i)}]: {int(d)}')

        y = - Xi * y // d

        errors.append((pos, y))

    return errors

# gf.run_tests()

def decode(msg, ecc_len):
    msg_enc = gf.P(GF, msg[::-1])
    synds = rs_syndromes(msg_enc, ecc_len)
    # print(f'synds: {list(map(int, synds))}')
    err_poly = berlekamp_massey(synds)
    # print(f'err_poly: {err_poly}')
    err_pos = []
    for i in range(len(msg)):
        if int(err_poly.eval(GF.gen(i).inv())) == 0:
            err_pos.append(i)
    print(f'err_pos: {err_pos}')
    err_mag = bch.forney(synds, err_poly, err_pos)

    # print(f'err_mag: {list(map(int, list(zip(*err_mag))[1]))}')
    msg_enc.x += [GF(0)] * (len(msg) - len(msg_enc.x))
    for pos, mag in err_mag:
        msg_enc.x[pos] += mag

    return list(map(int, msg_enc[::-1]))

if __name__ == '__main__':
    # random.seed(42)

    ecc_len = 4
    msg_len = 16

    # gf.run_tests()

    message = [(random.randrange(0, 256)*0 + 5) % GF.p ** GF.k for _ in range(msg_len)]
    print(f'message: {message}')

    generator = rs_generator(ecc_len)
    print(f'generator:', list(map(int, generator.x)))

    msg_enc = rs_encode_systematic(message[::-1], generator)
    assert any(map(int, rs_syndromes(msg_enc, ecc_len))) == False

    for i in range(ecc_len // 2):
        e_i = random.randrange(0, len(msg_enc.x))
        # e_i = i
        e = msg_enc.N(random.randrange(0, GF.p ** GF.k))
        msg_enc.x[e_i] += e
        print(f'error: {(e_i, int(e))}')

    synds = rs_syndromes(msg_enc, ecc_len)
    print(f'syndromes: {list(map(int, synds))}')

    err_poly = berlekamp_massey(synds)
    # print(f'err_poly: {err_poly}')

    err_pos = []
    for i in range(msg_len + ecc_len):
        if int(err_poly.eval(GF.gen(i).inv())) == 0:
            err_pos.append(i)

    assert len(err_pos) == err_poly.deg()

    # print(f'err_pos: {err_pos}')

    # errors1 = forney(synds, err_poly, err_pos)
    errors = forney(synds, err_poly, err_pos)
    print('errors', list(map(lambda p: (p[0], int(p[1])), errors)))

    msg_enc.x += [0] * (msg_len + ecc_len - len(msg_enc.x))
    print(f'msg_enc: {list(map(int, msg_enc[::-1]))}')

    for pos, mag in errors:
       msg_enc.x[pos] -= mag

    msg_dec = list(map(int, msg_enc[::-1]))

    # print(f'message: {message}')

    print(f'msg_dec: {msg_dec[:msg_len + ecc_len]}')

    print(message == msg_dec[:msg_len])
