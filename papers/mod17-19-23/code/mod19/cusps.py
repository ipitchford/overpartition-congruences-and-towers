# Explicit checks remain active with PYTHONOPTIMIZE=1.
"""Exact principal parts of G = sum_{n=3 (8), (n/19)=1} c(n) q^n at every cusp of Gamma0(L).

F4 = eta(2z)^19 / (eta(z)^14 eta(4z)^6) = sum c(n) q^n, c(n) = pbar(4n).
G(z) = sum_t W_t F4(z + t/152).  At a cusp s = a/c (A_s in SL2(Z), c>0) we compute
(-i(c tau+d))^{1/2} * G(A_s tau) as a Puiseux series in tau-exponents (e(E tau)),
using the exact eta transformation law, and report every exponent E <= 0 with its
coefficient (computed to 60 digits; phases are exact rationals).
Also reports the order of E2 = eta(19z)^19/eta(361z) at the same cusp (same units).
"""
import sys
from fractions import Fraction as Fr
from math import gcd
import mpmath as mp
mp.mp.dps = 60
import os, json
ETA_F4 = {int(k): v for k, v in json.loads(os.environ.get('ETAQ', '{"2": 19, "1": -14, "4": -6}')).items()}
WT = Fr(sum(ETA_F4.values()), 2)
E2 = {19: 19, 361: -1}

def e(x):
    return mp.expjpi(2 * mp.mpf(x.numerator) / x.denominator)

def dedekind(d, c):
    s = Fr(0)
    for n in range(1, c):
        x = Fr(n, c)
        y = Fr(d * n % c, c)
        if y != 0:
            s += (x - Fr(1, 2)) * (y - Fr(1, 2))
    return s
_ds = {}

def epsilon(a, b, c, d):
    if not (c > 0 and a * d - b * c == 1):
        raise AssertionError('mathematical verification failed')
    key = (d % c, c)
    if key not in _ds:
        _ds[key] = dedekind(d % c, c)
    return Fr(a + d, 24 * c) - _ds[key] / 2

def egcd(a, b):
    if b == 0:
        return (a, 1, 0) if a >= 0 else (-a, -1, 0)
    g, x, y = egcd(b, a % b)
    return (g, y, x - a // b * y)

def decompose(p, q, r, s):
    """M=[[p,q],[r,s]] (det>0, r>0) = gamma * U, gamma in SL2(Z) with c>0, U=[[u1,u2],[0,u3]], u1>0."""
    g, x0, y0 = egcd(p, r)
    ga, gb, gc, gd = (p // g, -y0, r // g, x0)
    u1 = gd * p - gb * r
    u2 = gd * q - gb * s
    u3 = -gc * q + ga * s
    if not (-gc * p + ga * r == 0 and u1 == g):
        raise AssertionError('mathematical verification failed')
    return ((ga, gb, gc, gd), (u1, u2, u3))

def eta_quot_series(eq, gamma, U, cap):
    a, b, c, d = gamma
    phase = Fr(0)
    mag2 = Fr(1)
    facs = []
    lead = Fr(0)
    u1, u2, u3 = U
    for delta, r in eq.items():
        (ga, gb, gc, gd), (v1, v2, v3) = decompose(delta * a, delta * b, c, d)
        phase += r * epsilon(ga, gb, gc, gd)
        mag2 *= Fr(v1, delta) ** r
        A = v1 * u1
        B = v1 * u2 + v2 * u3
        Z = v3 * u3
        phase += r * Fr(B, 24 * Z)
        lead += Fr(A, 24 * Z) * r
        facs.append((r, A, B, Z))
    budget = cap - lead
    series = {Fr(0): mp.mpc(1)}
    for r, A, B, Z in facs:
        step = Fr(A, Z)
        n = 1
        while n * step <= budget:
            X = e(Fr(n * B, Z))
            fac = {}
            m = 0
            while m * n * step <= budget:
                coef = mp.binomial(r, m) * (-X) ** m
                if coef == 0:
                    break
                fac[m * n * step] = coef
                m += 1
            new = {}
            for E1, w1 in series.items():
                for E2_, w2 in fac.items():
                    EE = E1 + E2_
                    if EE <= budget:
                        new[EE] = new.get(EE, 0) + w1 * w2
            series = new
            n += 1
    return (phase, mag2, {E + lead: v for E, v in series.items()})

def F4_direct(z):

    def eta(z):
        q = mp.exp(2j * mp.pi * z)
        return mp.exp(2j * mp.pi * z / 24) * mp.qp(q)
    return eta(2 * z) ** 19 / (eta(z) ** 14 * eta(4 * z) ** 6)

def G_at_cusp(a, c, cap=Fr(1, 10 ** 6), nt=12):
    """principal part (exponents <= cap) of (-i(c tau+d))^{1/2} G(A_s tau)."""
    g_, d, mb = egcd(a, c)
    if not g_ == 1:
        raise AssertionError('mathematical verification failed')
    d0, b0 = (d, -mb)
    if not a * d0 - b0 * c == 1:
        raise AssertionError('mathematical verification failed')
    chi = {k: 1 if pow(k, 9, 19) == 1 else -1 for k in range(1, 19)}
    gs = mp.mpc(0, mp.sqrt(19))
    total = {}
    global ABS
    ABS = {}
    for j in range(8):
        for k in range(19):
            bk = Fr(9, 19) if k == 0 else None
            W = mp.mpf(1) / 8 * e(Fr(-3 * j, 8))
            W *= mp.mpf(9) / 19 if k == 0 else -(1 + chi[k] * gs) / 38
            t = 19 * j + 8 * k
            p, q, r, s = (152 * a + t * c, 152 * b0 + t * d0, 152 * c, 152 * d0)
            gam, U = decompose(p, q, r, s)
            lam = Fr(152, U[2])
            phase, mag2, ser = eta_quot_series(ETA_F4, gam, U, cap)
            k_w = Fr(-1, 2)
            K = e(phase) * mp.sqrt(mp.mpf(mag2.numerator) / mag2.denominator) * (mp.mpf(lam.numerator) / lam.denominator) ** (mp.mpf(WT.numerator) / WT.denominator)
            for E, v in ser.items():
                if E <= cap:
                    total[E] = total.get(E, 0) + W * K * v
                    ABS[E] = ABS.get(E, 0) + abs(W * K * v)
    return total

def E2_order(a, c):
    g_, d, mb = egcd(a, c)
    gam = (a, -mb, c, d)
    o = Fr(0)
    for delta, r in E2.items():
        (ga, gb, gc, gd), (v1, v2, v3) = decompose(delta * a, delta * -mb, c, d)
        o += r * Fr(v1, 24 * v3)
    return o

def cusps(L):
    out = []
    for c in range(1, L + 1):
        if L % c:
            continue
        m = gcd(c, L // c)
        seen = set()
        for a in range(1, 10 * L):
            if gcd(a, c) != 1:
                continue
            cls = a % m
            if cls in seen or (gcd(cls, m) != 1 and m > 1):
                continue
            seen.add(cls)
            out.append((a, c))
            if len(seen) == sum((1 for x in range(m) if gcd(x, m) == 1)):
                break
    return out
if __name__ == '__main__':
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 64 * 361
    for a, c in cusps(L):
        tot = G_at_cusp(a, c)
        polar = {E: v for E, v in sorted(tot.items()) if abs(v) > mp.mpf(10) ** (-40)}
        neg = {E: v for E, v in polar.items() if E < 0}
        print(f"cusp {a}/{c}: E2ord={E2_order(a, c)} polar={[(str(E), mp.nstr(abs(v), 8)) for E, v in neg.items()]} const={('nonzero' if Fr(0) in polar else '0')}", flush=True)
