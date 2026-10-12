# Explicit checks remain active with PYTHONOPTIMIZE=1.
"""Independent exact checks; does not import any submitted implementation.

This checks the four-dissection by the Euler product, all cusp-order signs and
translated denominators, class descriptions, and multiplier arithmetic on a
bounded matrix family. It does not replace the exact principal-part check.
"""
from fractions import Fraction as Q
from math import gcd, isqrt
import json
from pathlib import Path
N = 64 * 19 ** 2
D = 8 * 19
ETA = {1: 5, 2: 19, 4: -6, 19: 18, 361: -1}
F4_ETA = {1: -14, 2: 19, 4: -6}

def order(exponents, denominator):
    return sum((Q(r * gcd(denominator, delta) ** 2, delta) for delta, r in exponents.items()), Q()) / 24

def cusp_representatives():
    result = []
    for v2 in range(7):
        for v19 in range(3):
            c = 2 ** v2 * 19 ** v19
            m = gcd(c, N // c)
            for residue in range(m):
                if gcd(residue, m) != 1:
                    continue
                a = residue
                while gcd(a, c) != 1:
                    a += m
                result.append((a, c))
    return result
cusps = cusp_representatives()
if not len(cusps) == 240:
    raise AssertionError('mathematical verification failed')
critical = []
ordinary_translates = 0
critical_groups = {}
for a, c in cusps:
    counts = {'positive': 0, 'zero': 0, 'negative': 0}
    leading_orders = set()
    for t in range(D):
        target = Q(a, c) + Q(t, D)
        cp = target.denominator
        transformed_order = order(ETA, cp) * Q(c, cp) ** 2
        kind = 'positive' if transformed_order > 0 else 'negative' if transformed_order < 0 else 'zero'
        counts[kind] += 1
        leading_orders.add(transformed_order)
        if c % 361:
            if not cp % 361 != 0:
                raise AssertionError((a, c, t, cp))
            if not transformed_order > 0:
                raise AssertionError((a, c, t, cp, transformed_order))
            ordinary_translates += 1
    if c % 361 == 0:
        critical.append({'cusp': f'{a}/{c}', 'F4_order': str(order(F4_ETA, c)), 'Fpp_order': str(order(ETA, c)), 'translation_counts': counts, 'leading_orders': [str(x) for x in sorted(leading_orders)]})
        critical_groups[str(c)] = str(order(ETA, c))
if not ordinary_translates == 228 * 152 == 34656:
    raise AssertionError('mathematical verification failed')
if not len(critical) == 12:
    raise AssertionError('mathematical verification failed')

def overpartitions_from_product(M):
    out = [1] + [0] * M
    for k in range(1, M + 1):
        for n in range(k, M + 1):
            out[n] += out[n - k]
        for n in range(M, k - 1, -1):
            out[n] += out[n - k]
    return out

def multiply(a, b, M):
    out = [0] * (M + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(min(len(b), M + 1 - i)):
                out[i + j] += x * b[j]
    return out

def power(a, exponent, M):
    out = [1] + [0] * M
    for _ in range(exponent):
        out = multiply(out, a, M)
    return out

def inverse(a, M):
    if not a[0] == 1:
        raise AssertionError('mathematical verification failed')
    out = [1] + [0] * M
    support = [(i, x) for i, x in enumerate(a) if i and x]
    for n in range(1, M + 1):
        out[n] = -sum((x * out[n - i] for i, x in support if i <= n))
    return out
M = 256
pb = overpartitions_from_product(4 * M)
phi = [0] * (M + 1)
phi[0] = 1
for k in range(1, isqrt(M) + 1):
    phi[k * k] = 2
phim = [(-1) ** n * x for n, x in enumerate(phi)]
F4 = multiply(power(phi, 3, M), power(inverse(phim, M), 4, M), M)
if not F4 == [pb[4 * n] for n in range(M + 1)]:
    raise AssertionError('mathematical verification failed')
if not pb[:10] == [1, 2, 4, 8, 14, 24, 40, 64, 100, 154]:
    raise AssertionError('mathematical verification failed')

def legendre(a, p):
    return 0 if a % p == 0 else 1 if pow(a, (p - 1) // 2, p) == 1 else -1
residue_checks = []
for p, ell, r, sign, expected in [(17, 11, 3, -1, [19, 35, 43, 51, 83]), (23, 13, 5, 1, [29, 53, 61, 69, 77, 101])]:
    got = [m for m in range(8 * ell) if m % 8 == r and legendre(m, ell) == sign]
    if not got == expected:
        raise AssertionError('mathematical verification failed')
    residue_checks.append({'p': p, 'classes': got})

def jacobi(a, n):
    n = abs(n)
    if not (n > 0 and n % 2 == 1):
        raise AssertionError('mathematical verification failed')
    a %= n
    value = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                value = -value
        a, n = (n, a)
        if a % 4 == n % 4 == 3:
            value = -value
        a %= n
    return value if n == 1 else 0
matrix_cases = 0
for u in range(1, 13):
    c = N * u
    for d in range(-61, 62, 2):
        if gcd(c, d) != 1:
            continue
        a = pow(d, -1, c)
        b = (a * d - 1) // c
        for t in range(D):
            tp = t * d * d % D
            ap = a + Q(t * c, D)
            dp = d - Q(c * tp, D)
            bp = b + Q(t * d, D) - ap * Q(tp, D)
            if not all((x.denominator == 1 for x in (ap, dp, bp))):
                raise AssertionError('mathematical verification failed')
            ap, dp, bp = (int(ap), int(dp), int(bp))
            if not ap * dp - bp * c == 1:
                raise AssertionError('mathematical verification failed')
            if not (d - dp) % (8 * u) == 0:
                raise AssertionError('mathematical verification failed')
            if not (dp % 4 == d % 4 and jacobi(c, dp) == jacobi(c, d)):
                raise AssertionError('mathematical verification failed')
            matrix_cases += 1
index = Q(N) * Q(3, 2) * Q(20, 19)
if not index == 36480:
    raise AssertionError('mathematical verification failed')
bound_square = Q(35, 12) * index
if not bound_square == 106400:
    raise AssertionError('mathematical verification failed')
if not bound_square / 2 == 53200:
    raise AssertionError('mathematical verification failed')
if not (pow(17, 16, 19) == 5 and pow(17, 33, 19) == pow(17, 15, 19) == 7):
    raise AssertionError('mathematical verification failed')
if not (pow(11, 6, 17) == 8 and pow(13, 9, 23) == 3):
    raise AssertionError('mathematical verification failed')
report = {'status': 'passed', 'cusp_count': len(cusps), 'ordinary_cusps': 228, 'strictly_positive_translated_cusp_orders': ordinary_translates, 'critical_cusps': critical, 'four_dissection_checked': {'n_min': 0, 'n_max': M, 'pbar_max_argument': 4 * M}, 'residue_classes': residue_checks, 'theta_multiplier_matrix_cases': matrix_cases, 'matrix_scope': 'c=23104*u, 1<=u<=12; odd -61<=d<=61 coprime to c; 0<=t<152', 'index': int(index), 'square_bound': int(bound_square), 'half_weight_bound': 53200, 'Hecke_residues': {'17^16 mod19': 5, '17^33 mod19': 7, '11^6 mod17': 8, '13^9 mod23': 3}, 'limitations': ['Finite series and matrix checks do not prove unrestricted identities.', 'Cusp-order checks establish positivity away from 361|c; critical cancellation needs the separately replayed exact cyclotomic certificate.']}
Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
