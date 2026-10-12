# Explicit checks remain active with PYTHONOPTIMIZE=1.
from collections import Counter
from math import gcd
from sympy import divisors, factorint

def phi(d, e, C):
    C[2 * d] += 5 * e
    C[d] -= 2 * e
    C[4 * d] -= 2 * e

def phim(d, e, C):
    C[d] += 2 * e
    C[2 * d] -= e

def Zc():
    C = Counter()
    phim(1, 5, C)
    phi(7, 4, C)
    phi(14, 6, C)
    phi(28, 12, C)
    phi(56, 1, C)
    return C
N = 224
res = []
for name, (a, b) in {'X^2Z': (8, 0), 'XYZ': (4, 4), 'Y^2Z': (0, 8)}.items():
    C = Zc()
    phi(1, a, C)
    phim(1, b, C)
    C = {d: e for d, e in C.items() if e}
    k = sum(C.values()) / 2
    if not all((N % d == 0 for d in C)):
        raise AssertionError('mathematical verification failed')
    c1 = sum((d * e for d, e in C.items())) % 24
    c2 = sum((N // d * e for d, e in C.items())) % 24
    s = 1
    for d, e in C.items():
        s *= d ** e
    from fractions import Fraction
    s = Fraction(1)
    for d, e in C.items():
        s *= Fraction(d) ** e
    sq = {}
    for p, m in factorint(s.numerator).items():
        sq[p] = sq.get(p, 0) + m
    for p, m in factorint(s.denominator).items():
        sq[p] = sq.get(p, 0) - m
    core = 1
    for p, m in sq.items():
        if m % 2:
            core *= p
    from fractions import Fraction
    ords = []
    for c in divisors(N):
        o = sum((Fraction(gcd(c, d) ** 2 * e, gcd(c, N // c) * c * d) for d, e in C.items())) * Fraction(N, 24)
        ords.append(o)
    print(name, 'k=', k, 'sum d r mod24=', c1, 'sum N/d r mod24=', c2, 'char=((-1)^k*', core, '/.)', 'min cusp order', min(ords))
    if k != 18 or c1 or c2 or (core != 14) or (min(ords) < 0):
        raise ValueError('eta modularity/cusp certificate failed')
    print('eta vector', sorted(C.items()), 'cusp orders', list(zip(divisors(N), map(str, ords))))
