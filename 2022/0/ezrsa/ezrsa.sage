#!/usr/bin/env sage

from Crypto.Util.number import *
from os import urandom

### Given ###


def xgcd(a, b):
    """Return xa + by = q where q is gcd(a, b).

    :param a: A number.
    :param b: A number.
    """

    if 0 == b:
        return 1, 0, a
    x, y, q = xgcd(b, a % b)
    x, y = y, (x - a // b * y)
    return x, y, q


def gen(p_size, d_size, l_size):
    """generate a public key, private key, and hint.

    :param p_size: Public key size.
    :param d_size: Private key size.
    :param l_size: Hint size.
    """

    while True:
        p = getPrime(p_size)
        q = getPrime(p_size)
        if GCD(p - 1, q - 1) == 2:
            break
    d_p = getPrime(d_size)
    d_q = getPrime(d_size)
    s, t, g = xgcd(p - 1, q - 1)
    if s < 0:
        s += q - 1
    else:
        t += p - 1
    n = p * q
    phi = (p - 1) * (q - 1)
    e = (
        (inverse(d_p, p - 1) * t * (q - 1) + inverse(d_q, q - 1) * s * (p - 1))
        // g
        % phi
    )
    assert e * d_q % (q - 1) == 1
    assert e * d_p % (p - 1) == 1
    k = (e * d_p - 1) // (p - 1)
    l = (e * d_q - 1) // (q - 1)
    return (n, e), (d_p, d_q, p, q), (d_p % (2**l_size), d_q % (2**l_size))


def encrypt(m, pk):
    """Get the ciphertext.

    :param m: Plaintext you want to encrypt.
    :param pk: Public key.
    """

    n, e = pk
    return pow(m, e, n)


def ezrsa():
    """ezrsa challenge. Given the return values, get the flag."""

    from secret import flag

    pk, sk, hint = gen(p_size, d_size, l_size)
    flag = urandom(2 * p_size // 8 - len(flag) - 1) + flag
    enc = encrypt(int(flag.hex(), 16), pk)
    return enc, pk, hint


p_size = 1000
d_size = 105
l_size = 55
# We are given the output of a secret call to ezrsa
enc = 35558284230663313298312684064040643811204702946900174110911295087662938676356112802781671547473910691476600838877279843972105403072929243674403244286458898562457747942651643439624568905004454158744508429126554955023110569348839934098381885098523538078300248638407684468503519326866276798222721018258242443186786917829878515320321445508466038372324063139762003962072922393974710763356236627711414307859950011736526834586028087922704589199885845050751932885698053938070734392371814246294798366452078193195538346218718588887085179856336533576097324041969786125752304133487678308830354729347735644474025828
pk = (
    44774502335951608354043148360684114092901940301155357314508676399067538307546121753785009844275454594381602690061553832466871574728524408152400619047820736137949166290404514747591817206669966103047443912935755873432503095952914080827536130899968275165557303493867755627520568588808534411526058896791373252974606364861105086430729757064078675811147972536205086402752245214343186536177015741922559575575911278873118556603923689408629477875537332177644886701517140711134017511229202430437068095342526435886609381176269251580339549071944830141516001532825295594908434587285225415103472279090325281062442217,
    29624366183227462965645558392954094074485353876807451497147549927093025197118051280445930543762170853769573962200247669305286333212410439624262142109295839433584663989554419810341266820063074908743295553517790354149623873028162282751352613333181218478850463012413786673509078012976454604598813805735677104174112776060905225493357010861225261560490401501912259585922988353328944884443953564154752191932500561561256069872534626325000901099904014035414792860997025614313564862063784602254606240743545483125618939111639728114664995759380293512809125885893543730614962375399353971677980309835647540883700977,
)
hint = (5013415024346389, 4333469053087705)

### Solution ###

N, e = pk
P = PolynomialRing(ZZ, "xp, xq, yp, yq, zp, zq")
(
    xp,
    xq,
    yp,
    yq,
    zp,
    zq,
) = P.gens()
f = xp * yp - xq
g = yp * zp - N * zq
h = N * xp * zq - xq * zp
f1 = xp * yp - xq - e * hint[0]
g1 = yp * zp - N * zq + e * hint[1] * yp
h1 = (
    N * xp * zp
    - xq * zp
    - e**2 * hint[0] * hint[1]
    - e * hint[0] * zp
    - e * hint[1] * xq
)


def get_quo(f, g):
    """f = x * g + r.
    We return x, r.

    :param f: The bigger polynomial.
    :param g: The smaller one.
    """

    r = f.reduce([g])
    return (f - r) // g, r


def replace_with(monomial, replaced, replacer):
    """substitute all the replaced terms in the monomial
    with replacer.

    :param monomial: a polynomial.
    :param replaced: a polynomial.
    :param replacer: a polynomial.
    """

    while 1:
        quotient, remainder = get_quo(monomial, replaced)
        if remainder != 0:
            break
        monomial = quotient * replacer
    return monomial


def trans(f):
    """trans is the operation on a polynomial
    that applies a series of operations on its monomials.
    It creates two copies of the Newton polytope of f
    where one lies in the (xp, yp, zp) plane and the other lies
    in the (xq, yq, zq) plane.

    :param f: the polynomial.
    """

    f_temp = 0
    for coeff, monomial in list(f):
        monomial = replace_with(monomial, yp * yq, N)
        quotient, remainder = get_quo(monomial, yp)
        if quotient != 0:
            monomial = replace_with(monomial, xp, xq + 1)
            monomial = replace_with(monomial, zp, zq - 1)
        else:
            monomial = replace_with(monomial, xp, xq - 1)
            monomial = replace_with(monomial, zp, zq + 1)
        f_temp += coeff * monomial
    return f_temp


# Shift polynomials will have some root (x0, y0, z0) mod e^(2m)
# This is just a random natural number
m = 6
# Calculate M partitions
M1 = []
M2 = []
M3 = []
M4 = []
for a in range(m + 1):
    for c in range(m + 1):
        for b in range(a + c + 1):
            poly = xp**a * yp**b * zp**c
            if a <= c and b <= c - a:
                M1.append(poly)
            elif a > c and b < a - c:
                M2.append(poly)
            elif (a + b + c) % 2 == 0:
                M3.append(poly)
            else:
                M4.append(poly)


def power_of(monomial, term):
    """Get the power of the term in a monomial.

    :param monomial: A monomial.
    :param term: A term.
    """

    i = 0
    while 1:
        quotient, remainder = get_quo(monomial, replaced)
        if remainder != 0:
            break
        monomial = quotient * replacer
        i += 1
    return i


# Define exponent functions
def Ef(a, b, c):
    """Ef.

    :param a: a number.
    :param b: a number.
    :param c: a number.
    """
    f = xp**a + yp**b + zp**c
    if f in M1:
        return 0
    elif f in M2:
        return b
    elif f in M3:
        return (a + b - c) // 2
    elif f in M4:
        return (a + b - c + 1) // 2
    else:
        raise Exception("f should be in M1, M2, M3, or M4")


def Eg(a, b, c):
    f = xp**a + yp**b + zp**c
    if f in M1:
        return b
    elif f in M2:
        return 0
    elif f in M3:
        return (-a + b + c) // 2
    elif f in M4:
        return (-a + b + c - 1) // 2
    else:
        raise Exception("f should be in M1, M2, M3, or M4")


def Eh(a, b, c):
    f = xp**a + yp**b + zp**c
    if f in M1:
        return a
    elif f in M2:
        return c
    elif f in M3:
        return (a - b + c) // 2
    elif f in M4:
        return (a - b + c - 1) // 2
    else:
        raise Exception("f should be in M1, M2, M3, or M4")


def Ex(a, b, c):
    f = xp**a + yp**b + zp**c
    if f in M2:
        return a - b - c
    elif f in M1 or f in M3 or f in M4:
        return 0
    else:
        raise Exception("f should be in M1, M2, M3, or M4")


def Ez(a, b, c):
    f = xp**a + yp**b + zp**c
    if f in M1:
        return -a - b + c
    elif f in M2 or f in M3:
        return 0
    elif f in M4:
        return 1
    else:
        raise Exception("f should be in M1, M2, M3, or M4")


### Testing ###

assert trans(yp * yq) == N
assert trans((yp * yq) ** 2) == N**2
