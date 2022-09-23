#!/usr/bin/env sage

from Crypto.Util.number import bytes_to_long as bl, long_to_bytes as lb, isPrime, GCD
from os import urandom
from hashlib import sha384

### Given ###

magic_len = 17
magic_num = bl(urandom(magic_len))
s = 248


def hash_(a):
    return int(sha384(lb(a)).hexdigest(), 16)


def magic_dlog(P, E, data):
    assert P >> s == magic_num
    assert isPrime(P)
    assert pow(data, E, P) == hash_(data) % P
    print("success")


### Solution ###


def crt_list_decoding(B, p, d=4):
    """Calculate a smooth number.

    :param B: Pair of the lower and upper bound.
    :param p: List of the first n primes.
    :param d: LLL parameter.
    """

    def fill_zero(l, len_=d * 2):
        """Fills an array of leading zeroes.

        :param l: An array l, indicating coeffs of
        a polynomial from increasing power left to right.
        :param len_: The length to fill the array to.
        """
        return l + [0] * (len_ - len(l))

    R = -B[0]
    B = B[1] - B[0]
    P = product(p)

    Z = PolynomialRing(ZZ, "x")
    gx = [P ** (d - i) * (x - R) ** i for i in range(d)]
    hx = [(x - R) ** d * x**i for i in range(d)]

    L = Matrix(ZZ, [fill_zero(f.subs(x=x * B).list()) for f in gx + hx]).LLL()
    poly = Z(list(map(int, L[0]))).subs(x=x / B)
    return -Integer(poly.roots()[0][0])


lower_bound = magic_num << s
upper_bound = lower_bound + (1 << s) - 1

P = crt_list_decoding((lower_bound, upper_bound), primes_first_n(384))
print(factor(P))
