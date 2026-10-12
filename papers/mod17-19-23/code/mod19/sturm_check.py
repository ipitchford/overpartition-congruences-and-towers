# Explicit checks remain active with PYTHONOPTIMIZE=1.
"""Sturm check of g|T(17^2) == 14 g (mod 19), g = R(F''), weight 35/2, Gamma0(64*361), trivial character.
Coefficients of g are c(n)*[n=3 mod 8]*[(n/19)=1] mod 19, with c(n)=pbar(4n): sum c(n)q^n = phi(q)^3/phi(-q)^4.
T(p^2) (Shimura, psi trivial, lambda=17): b(n)=a(p^2 n)+((-1)^17 n/p) p^16 a(n)+p^33 a(n/p^2).
Sturm bound: (35/2)/12 * [SL2:Gamma0(23104)] = 53200."""
import sys, time
import flint
P = 19
idx = 64 * 3 // 2 * 361 * 20 // 19
B = (35 * idx + 23) // 24
M = 289 * B + 10
t0 = time.time()

def theta(sign):
    co = [0] * (M + 1)
    co[0] = 1
    k = 1
    while k * k <= M:
        co[k * k] = 2 * sign ** k % P
        k += 1
    return flint.nmod_poly(co, P)
phi = theta(1)
phim = theta(-1)
inv = phim.inverse_series_trunc(M + 1)
F = phi.pow_trunc(3, M + 1).mul_low(inv.pow_trunc(4, M + 1), M + 1)
c = [int(x) for x in F.coeffs()]
c += [0] * (M + 1 - len(c))
print('series done', round(time.time() - t0, 1), 's; M =', M, 'B =', B, flush=True)
if not (c[0] == 1 and c[1] == 14 % P and (c[2] == 100 % P)):
    raise AssertionError(c[:3])

def leg(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def G(n):
    if n < 0 or n != int(n):
        return 0
    return c[n] if n % 8 == 3 and leg(n, 19) == 1 else 0
p17 = pow(17, 16, P)
p33 = pow(17, 33, P)
if not (p17 == 5 and p33 == pow(17, 15, P)):
    raise AssertionError('mathematical verification failed')
bad = []
nonzero = 0
for n in range(1, B + 1):
    a2 = G(289 * n)
    a1 = G(n)
    a0 = G(n // 289) if n % 289 == 0 else 0
    b = (a2 + leg(-n, 17) * p17 * a1 + p33 * a0 - 14 * a1) % P
    if a1 or a2:
        nonzero += 1
    if b:
        bad.append(n)
print('checked n <= B =', B, '; n with nonzero data:', nonzero, '; failures:', len(bad), bad[:20])
cnt = 0
fail = 0
for n in range(3, M // 289 + 1, 8):
    if leg(n, 19) == 1 and leg(n, 17) == -1:
        cnt += 1
        if c[289 * n] % P:
            fail += 1
print('conjecture instances checked:', cnt, 'failures:', fail, 'time', round(time.time() - t0, 1))
if bad or fail or B != 53200 or (cnt != 1487) or (nonzero == 0):
    raise ValueError('mod19 Hecke-Sturm certificate failed')
