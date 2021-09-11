#!/usr/bin/env python3

import random
import sys
import itertools

# random.seed(42)

import gf
from gf import P

# GF = gf.GF(2, 6, *next(gf.GF.irr_polynomials(2, 6, 2)))
GF = gf.GF256
# print('irr', GF.poly)

# for i in range(GF.p ** GF.k):
#     print(f'{int(GF.gen(i)):02x}', end=' ')
# print()

def bch_generator(corr : int):
    '''
    @param: corr bits to correct
    '''

    mp = set(map(lambda x: x[1],
        itertools.takewhile(lambda x: x[0] < 2 * corr + 1,
            GF.minimal_polynomials())))

    g = P(GF, [1])

    for m in mp:
        g *= m
    return g

def rs_encode(msg, genertor):
    poly = P(GF, msg) * genertor

    return poly #// P(GF, [poly.x[-1]])

def rs_encode_systematic(msg, generator):
    poly = P(GF, [0] * generator.deg() + msg)

    return poly - (poly % generator)

def bch_syndromes(poly, d):
    return [poly.eval(GF.gen(i + 1)) for i in range(d)]

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

generator = bch_generator(4)

ecc_len = generator.deg()
msg_len = GF.p ** GF.k - ecc_len - 1


def test():
    print(f'generator:', generator)
    print(f'msg_len: {msg_len}')
    print(f'ecc_len: {ecc_len}')
    print(f'msg+ecc: {msg_len+ecc_len}')

    message = [random.randrange(0, 2) for _ in range(msg_len)]
    # print(f'message: {message}')

    msg_enc = rs_encode_systematic(message[::-1], generator)
    msg_enc.x += [GF(0)] * (msg_len + ecc_len - len(msg_enc.x))

    assert any(map(int, bch_syndromes(msg_enc, 2*3))) == False

    for i in range(3):
        e_i = random.randrange(0, len(msg_enc.x))
        e = msg_enc.N(1)
        msg_enc.x[e_i] += e
    #     print(f'error: {(e_i, int(e))}')
    # print(f'msg_enc: {list(map(int, msg_enc[::-1]))}')

    synds = bch_syndromes(msg_enc, 2*3)
    # print(f'syndromes: {list(map(int, synds))}')

    err_poly = berlekamp_massey(synds)
    # print(f'err_poly: {err_poly}')

    # find roots of the error locator polynomial
    # could use Chien search
    err_pos = []
    for i in range(msg_len + ecc_len):
        if int(err_poly.eval(GF.gen(i).inv())) == 0:
            err_pos.append(i)

    assert len(err_pos) == err_poly.deg()

    print(f'err_pos: {err_pos}')

    errors = forney(synds, err_poly, err_pos)

    print(f'msg_enc: {"".join(map(lambda x: str(int(x)), msg_enc[::-1][:msg_len]))}')

    for pos, mag in errors:
        msg_enc.x[pos] = int(msg_enc.x[pos] + mag)

    msg_dec = list(map(int, msg_enc[::-1]))

    print(f'message: {"".join(map(str, message))}')

    print(f'msg_dec: {"".join(map(str, msg_dec[:msg_len]))}')

    assert (message == msg_dec[:msg_len])

if __name__ == '__main__':
    test()
