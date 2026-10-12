# Explicit checks remain active with PYTHONOPTIMIZE=1.
from flint import nmod_poly
N = 3136
idx = N * 3 // 2 * 8 // 7
B = 18 * idx // 12
if not (N, idx, B) == (3136, 5376, 8064):
    raise AssertionError('mathematical verification failed')
n = B + 1

def theta(d, sign):
    c = [0] * n
    k = 0
    while d * k * k < n:
        c[d * k * k] = (1 if k == 0 else 2) * sign ** k % 7
        k += 1
    return nmod_poly(c, 7)
pw = lambda a, e: a.pow_trunc(e, n)
X = pw(theta(1, 1), 4)
Y = pw(theta(1, -1), 4)
Q = X.mul_low(X, n) - X.mul_low(Y, n) + Y.mul_low(Y, n)
T = theta(56, 1).mul_low(pw(theta(7, 1), 4), n).mul_low(pw(theta(14, 1), 6), n).mul_low(pw(theta(28, 1), 12), n)
H = Q.mul_low(pw(theta(1, -1), 5), n).mul_low(T, n)
c = [int(x) for x in H.coeffs()]
c += [0] * (n - len(c))
QR = {1, 2, 4}
cls = [m for m in range(3, n, 8) if m % 7 in QR]
bad = [m for m in cls if c[m]]
print('N', N, 'index', idx, 'B', B, '| m<=B with m=3 mod 8,(m/7)=1:', len(cls), 'nonzero mod 7:', bad[:10])
print('control, (m/7)=-1 class nonzero count:', sum((1 for m in range(3, n, 8) if m % 7 in {3, 5, 6} and c[m])))
if bad or len(cls) != 432 or c[0] != 1:
    raise ValueError('mod7 Sturm certificate failed')
if sum((1 for m in range(3, n, 8) if m % 7 in {3, 5, 6} and c[m])) != 370:
    raise ValueError('mod7 complementary control failed')
