# Explicit checks remain active with PYTHONOPTIMIZE=1.
from flint import fmpz_poly, nmod_poly
M = 4000
L = 16 * M + 1

def theta(d, sign, n, mod=None):
    c = [0] * n
    k = 0
    while d * k * k < n:
        c[d * k * k] += (1 if k == 0 else 2) * sign ** k
        k += 1
    return nmod_poly(c, mod) if mod else fmpz_poly(c)

def psi(d, n):
    c = [0] * n
    k = 0
    while d * k * (k + 1) // 2 < n:
        c[d * k * (k + 1) // 2] += 1
        k += 1
    return fmpz_poly(c)

def inv(f, n):
    return fmpz_poly(nmod_poly(f.coeffs(), 2 ** 61 - 1).inverse_series_trunc(n).coeffs()) if False else None

def series_inv(f, n):
    g = fmpz_poly([1])
    m = 1
    while m < n:
        m = min(2 * m, n)
        g = g * (2 - (f * g).__class__(((f * g).coeffs() + [0] * m)[:m]))
        g = fmpz_poly(g.coeffs()[:m])
    return g

def tr(f, n):
    return fmpz_poly(f.coeffs()[:n])

def mul(a, b, n):
    return tr(a * b, n)

def pw(a, e, n):
    r = fmpz_poly([1])
    for _ in range(e):
        r = mul(r, a, n)
    return r

def E(f, n):
    c = f.coeffs()
    return fmpz_poly(c[0::2][:n])

def co(f, n):
    c = f.coeffs()
    return c + [0] * (n - len(c))
pb = co(series_inv(theta(1, -1, L), L), L)
if not pb[:7] == [1, 2, 4, 8, 14, 24, 40]:
    raise AssertionError('mathematical verification failed')
n = M + 1
phi = theta(1, 1, 2 * n)
phim = theta(1, -1, 2 * n)
if not co(mul(theta(1, 1, 4 * n), theta(1, -1, 4 * n), 4 * n), 4 * n) == co(pw(theta(2, -1, 4 * n), 2, 4 * n), 4 * n):
    raise AssertionError('mathematical verification failed')
if not co(pw(theta(1, 1, 4 * n), 2, 4 * n) + pw(theta(1, -1, 4 * n), 2, 4 * n), 4 * n) == co(2 * pw(theta(2, 1, 4 * n), 2, 4 * n), 4 * n):
    raise AssertionError('mathematical verification failed')
if not co(theta(1, 1, 4 * n) + theta(1, -1, 4 * n), 4 * n) == co(2 * theta(4, 1, 4 * n), 4 * n):
    raise AssertionError('mathematical verification failed')
iv = lambda f, m: series_inv(f, m)
S1 = mul(theta(2, 1, n), iv(pw(theta(1, -1, n), 2, n), n), n)
S2 = mul(pw(theta(1, 1, n), 3, n), iv(pw(theta(1, -1, n), 4, n), n), n)
E7 = E(pw(phi, 7, 2 * n), n)
S3 = mul(E7, iv(pw(theta(1, -1, n), 8, n), n), n)
E7b = E(pw(theta(1, 1, 4 * n), 7, 4 * n), 2 * n)
S4 = mul(E(mul(E7b, pw(phi, 8, 2 * n), 2 * n), n), iv(pw(theta(1, -1, n), 16, n), n), n)
for k, S in ((2, S1), (4, S2), (8, S3), (16, S4)):
    d = [i for i in range(n) if co(S, n)[i] != pb[k * i]]
    if not not d:
        raise AssertionError((k, d[:5], len(pb)))
print('exact S1..S4 verified for', n, 'coefficients')
X = pw(theta(1, 1, n), 4, n)
Y = pw(theta(1, -1, n), 4, n)
if not co(E(pw(phi, 8, 2 * n), n), n) == co(8 * mul(X, X, n) - 8 * mul(X, Y, n) + mul(Y, Y, n), n):
    raise AssertionError('mathematical verification failed')
print('E[phi^8] identity verified')
H = mul(mul(mul(X, X, n) - mul(X, Y, n) + mul(Y, Y, n), pw(theta(1, -1, n), 5, n), n), mul(mul(pw(theta(7, 1, n), 4, n), pw(theta(14, 1, n), 6, n), n), pw(theta(28, 1, n), 12, n), n), n)
R = mul(H, iv(pw(theta(56, -1, n), 24, n), n), n)
if not all(((a - b) % 7 == 0 for a, b in zip(co(R, n), [pb[16 * i] for i in range(n)]))):
    raise AssertionError('mathematical verification failed')
print('mod 7 congruence A == H/phi(-q^56)^24 verified for', n, 'coefficients')
from sympy import legendre_symbol as Ls
bad = [i for i in range(n) if i % 8 == 3 and i % 7 and (Ls(i % 7, 7) == 1) and pb[16 * i] % 7]
print('direct conjecture check n<=', M, 'failures', bad)
if bad:
    raise ValueError('direct congruence mismatch')
Ht = mul(H, theta(56, 1, n), n)
R2 = mul(Ht, iv(mul(pw(theta(56, -1, n), 24, n), theta(56, 1, n), n), n), n)
if not all(((a - b) % 7 == 0 for a, b in zip(co(R2, n), [pb[16 * i] for i in range(n)]))):
    raise AssertionError('mathematical verification failed')
print('mod 7 congruence A == Htilde/(phi(-q^56)^24 phi(q^56)) verified')
