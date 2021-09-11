import ctypes
import os
import random
import time

import gf
import bch

def load_bchlib(recompile=False):
    if (recompile or not os.path.exists('./bch') or
        os.path.getctime('./bch') < os.path.getctime('./bch.c')):
        os.system('gcc -O3 -Wall -shared ./bch.c -o bch')

    c_lib = ctypes.CDLL('./bch')

    c_lib.decode63_45.restype = ctypes.c_bool
    c_lib.decode63_30.restype = ctypes.c_bool
    c_lib.check63_45.restype = ctypes.c_bool
    # c_lib.gf_mul.restype = ctypes.c_uint8
    # c_lib.gf_div.restype = ctypes.c_uint8

    return c_lib

def test33(seed):
    # for i in range(1, 64):
    #     for j in range(1, 64):
    #         assert bchlib.gf_mul(i, j) == int(bch.GF(i) * bch.GF(j)), seed
    #         assert bchlib.gf_div(i, j) == int(bch.GF(i) // bch.GF(j)), seed

    # generator33 = gf.P(bch.GF, map(int, bch.GF.poly_from_int(2, 0x37CD0EB67)))

    random.seed(seed)
    data_orig = random.randrange(0, 2 ** 30)

    # print(f'data_orig: {data_orig:030b}')

    data_tx = ctypes.create_string_buffer((data_orig << (33 + 1)).to_bytes(8, 'big'))
    bchlib.encode63_30(data_tx)
    # assert bchlib.check63_30(data_tx), seed

    data_tx = int.from_bytes(data_tx[:8], 'big')

    # print(f'data_tx:   {data_tx:064b}')

    # data_tx_r = bch.rs_encode_systematic(list(map(int, f'{data_orig:030b}'))[::-1], generator33)
    # print(f'data_tx_r: {gf.GF.poly_to_int(2, data_tx_r) << 1 | 1:064b}')

    data_rx = data_tx

    errors = set(63 - random.randrange(0, 30) for _ in range(6))
    for i in errors:
        data_rx ^= 1 << i
    # print(f'errors:    {sorted(errors)}')

    # print(f'data_rx:   {data_rx:064b}')

    # data_rx_r = gf.P(bch.GF, map(int, gf.GF.poly_from_int(2, data_rx >> 1)))
    # print(f'data_rx_r: {gf.GF.poly_to_int(2, data_rx_r) << 1 | 1:064b}')

    # data = ctypes.create_string_buffer(data_rx.to_bytes(8, 'big'))
    # synds = ctypes.create_string_buffer(12)
    # synds_r = bch.bch_syndromes(data_rx_r, 12)
    # bchlib.syndromes(data, synds, 12)
    # print('synds_r:', list(map(int, synds_r)))
    # print('synds:  ', list(map(int, synds[:12])))
    # assert list(map(int, synds_r)) == list(map(int, synds[:12]))

    # err_poly_r = bch.berlekamp_massey(synds_r)
    # print('err_poly_r:', list(map(int, err_poly_r.x)))


    # err_pos_r = []
    # for i in range(63):
    #     if int(err_poly_r.eval(bch.GF.gen(i).inv())) == 0:
    #         err_pos_r.append(i+1)

    # assert len(err_pos_r) == err_poly_r.deg()

    # err_poly = ctypes.create_string_buffer(12)
    # bchlib.berlekamp_massey_6(synds, err_poly)
    # print('err_poly:  ', list(map(int, err_poly[:12])))

    # err_pos = ctypes.create_string_buffer(6)
    # error_count = bchlib.find_err_pos(err_poly, err_pos, 6)
    # print('error_count:', error_count)
    # print('err_pos: ', list(map(lambda x: x + 1, err_pos[:6])))
    # print(f'err_pos_r: {err_pos_r}')

    data_corr = ctypes.create_string_buffer(data_rx.to_bytes(8, 'big'))

    # assert not bchlib.check63_30(data_corr), seed
    assert bchlib.decode63_30(data_corr), seed
    # assert bchlib.check63_30(data_corr), seed

    data_corr = int.from_bytes(data_corr[:8], 'big')
    data_corr >>= 33 + 1
    # print(f'data_corr: {data_corr:030b}')

    assert data_corr == data_orig, seed

def test18(seed):
    random.seed(seed)
    data_orig = random.randrange(0, 2 ** 45)

    # print(f'data_orig: {data_orig:045b}')

    data_tx = ctypes.create_string_buffer((data_orig << (18 + 1)).to_bytes(8, 'big'))
    bchlib.encode63_45(data_tx)
    # assert bchlib.check63_45(data_tx), seed

    data_tx = int.from_bytes(data_tx[:8], 'big')

    # print(f'data_tx:   {data_tx:064b}')

    data_rx = data_tx

    errors = set(random.randrange(1, 64) for _ in range(3))
    for i in errors:
        data_rx ^= 1 << i
    # print(f'errors:    {errors}')

    # print(f'data_rx:   {data_rx:064b}')

    data_corr = ctypes.create_string_buffer(data_rx.to_bytes(8, 'big'))

    # assert not bchlib.check63_45(data_corr), seed
    assert bchlib.decode63_45(data_corr), seed
    # assert bchlib.check63_45(data_corr), seed

    data_corr = int.from_bytes(data_corr[:8], 'big')
    data_corr >>= 18 + 1
    # print(f'data_corr: {data_corr:045b}')

    # assert data_corr == data_orig, seed

if __name__ == '__main__':
    seed = random.randrange(0, 2 ** 64)

    bchlib = load_bchlib()

    for i in range(100000):
        test33(seed)
        test18(seed)
        seed = random.randrange(0, 2 ** 64)
