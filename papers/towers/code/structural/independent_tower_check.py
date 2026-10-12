# Explicit checks remain active with PYTHONOPTIMIZE=1.
"""Independent audit: numerical basis construction, operator algebra and classification.
No submitted checker is imported. Arithmetic is integral or reduced modulo p.
"""
from math import isqrt, gcd
from fractions import Fraction
from pathlib import Path
import json, time, platform, sys
import numpy as np
import sympy as sp
T0 = time.perf_counter()

def conv(a, b, N, p=None):
    if p is None:
        out = np.convolve(np.array(a, dtype=object), np.array(b, dtype=object))[:N]
        return [int(x) for x in out] + [0] * max(0, N - len(out))
    out = np.convolve(np.array(a, dtype=np.int64), np.array(b, dtype=np.int64))[:N] % p
    return out.tolist() + [0] * max(0, N - len(out))

def ppow(a, e, N, p=None):
    z = [1] + [0] * (N - 1)
    for _ in range(e):
        z = conv(z, a, N, p)
    return z

def theta(N):
    a = [0] * N
    a[0] = 1
    for k in range(1, isqrt(N - 1) + 1):
        a[k * k] = 2
    return a

def eta_F(N, p):
    a = [1] + [0] * (N - 1)
    for d in range(2, N, 2):
        for _ in range(4):
            for n in range(d, N):
                a[n] = (a[n] + a[n - d]) % p
    for d in range(4, N, 4):
        for _ in range(8):
            for n in range(N - 1, d - 1, -1):
                a[n] = (a[n] - a[n - d]) % p
    return [0] + a[:-1]

def Fdiv(N):
    return [sum((d for d in range(1, isqrt(n) + 1) if n % d == 0 for d in ([d, n // d] if d * d != n else [d]))) if n % 2 else 0 for n in range(N)]

def basis(lam, N, p=None):
    ph = theta(N)
    f = Fdiv(N)
    if p:
        f = [x % p for x in f]
    return [conv(ppow(f, a, N, p), ppow(ph, 2 * lam + 1 - 4 * a, N, p), N, p) for a in range((2 * lam + 1) // 4 + 1)]

def solve(im, B, p=None):
    a = []
    r = im[:]
    for j in range(len(B)):
        x = r[j]
        a.append(x if p is None else x % p)
        r = [r[n] - a[-1] * B[j][n] for n in range(len(r))]
        if p:
            r = [x % p for x in r]
    return a

def combine(v, B, n, p=None):
    a = [sum((v[j] * B[j][i] for j in range(len(B)))) for i in range(n)]
    return [x % p for x in a] if p else a

def leg(n, l):
    if n % l == 0:
        return 0
    return 1 if pow(n % l, (l - 1) // 2, l) == 1 else -1

def act(v, U, p):
    return [sum((x * row[j] for x, row in zip(v, U))) % p for j in range(len(v))]

def orbit(v, U, p):
    seen = {}
    vv = []
    while tuple(v) not in seen:
        seen[tuple(v)] = len(vv)
        vv.append(v)
        v = act(v, U, p)
    return (vv, seen[tuple(v)])

def min_annihilator(v, U, p):
    rows = []
    seen = []
    w = v[:]
    for degree in range(len(v) + 1):
        augmented = [[seen[j][i] for j in range(len(seen))] + [w[i]] for i in range(len(v))]
        if degree:
            nr = len(augmented)
            nc = degree
            r = 0
            piv = []
            for c in range(nc):
                rr = next((k for k in range(r, nr) if augmented[k][c] % p), None)
                if rr is None:
                    continue
                augmented[r], augmented[rr] = (augmented[rr], augmented[r])
                inv = pow(augmented[r][c], -1, p)
                augmented[r] = [z * inv % p for z in augmented[r]]
                for k in range(nr):
                    if k != r:
                        mult = augmented[k][c]
                        augmented[k] = [(a - mult * b) % p for a, b in zip(augmented[k], augmented[r])]
                piv.append(c)
                r += 1
            consistent = all((any((row[c] for c in range(nc))) or row[-1] == 0 for row in augmented))
            if consistent and len(piv) == nc:
                coeff = [0] * nc
                for i, c in enumerate(piv):
                    coeff[c] = -augmented[i][-1] % p
                return coeff + [1]
        seen.append(w)
        w = act(w, U, p)
    raise AssertionError('no annihilator')
reports = {}
for p, lam, xc, cases in [(13, 5, [1, 4, 1], [(7, 1, 7)]), (17, 7, [1, 13, 13, 0], [(11, 1, 3), (29, -1, 3)]), (23, 10, [1, 9, 5, 14, 17, 20], [(13, 1, 5)])]:
    K = len(xc)
    sb = 32 * (lam + 1)
    N = max(sb + 1, 841 * (K - 1) + 1)
    B = basis(lam, N, p)
    F = Fdiv(240)
    if not eta_F(240, p) == [x % p for x in F]:
        raise AssertionError('mathematical verification failed')
    k = (p - 1) // 2
    b = k // 2
    phpow = ppow(theta(p * b + 1), p - 1, p * b + 1, p)
    if not conv(combine(xc, B, b + 1, p), theta(b + 1), b + 1, p) == [phpow[p * n] for n in range(b + 1)]:
        raise AssertionError('mathematical verification failed')
    EB = basis(lam, 4 * (lam + 1) + 1)
    U = []
    for bb in EB:
        im = bb[::4][:lam + 2]
        v = solve(im, EB)
        if not im == combine(v, EB, lam + 2):
            raise AssertionError('mathematical verification failed')
        U.append(v)
    vv, start = orbit(xc, U, p)
    whole = []
    for i, v in enumerate(vv):
        f = combine(v, B, sb + 1, p)
        if all((f[n] == 0 for n in range(7, sb + 1, 8))):
            whole.append(i)
    rr = {'p': p, 'lambda': lam, 'U_exact': U, 'U_mod': [[z % p for z in row] for row in U], 'charpoly_mod': str(sp.Poly(sp.Matrix(U).charpoly().as_expr(), modulus=p).as_expr()), 'c_orbit_length': len(vv), 'c_orbit_return': start, 'c_annihilator_ascending': min_annihilator(xc, U, p), 'c_class7_zero_indices': whole, 'cases': []}
    for l, chi, cls in cases:
        images = []
        for bb in B:
            im = [(bb[l * l * n] + leg((-1) ** lam * n, l) * pow(l, lam - 1, p) * bb[n] + (pow(l, 2 * lam - 1, p) * bb[n // (l * l)] if n % (l * l) == 0 else 0)) % p for n in range(K)]
            images.append(solve(im, B, p))
        h = [(a - chi * pow(l, lam - 1, p) * z) % p for a, z in zip(act(xc, images, p), xc)]
        ho, hi = orbit(h, U, p)
        for v in ho:
            f = combine(v, B, sb + 1, p)
            if not all((f[n] == 0 for n in range(cls, sb + 1, 8))):
                raise AssertionError('mathematical verification failed')
        rr['cases'].append({'ell': l, 'T_mod': images, 'h': h, 'orbit': ho, 'return': hi, 'annihilator_ascending': min_annihilator(h, U, p), 'projection_sturm_bound': sb})
    if p == 13:
        extra = []
        for l in [5, 29]:
            images = []
            for bb in B:
                im = [(bb[l * l * n] + leg(-n, l) * pow(l, 4, p) * bb[n] + (pow(l, 9, p) * bb[n // (l * l)] if n % (l * l) == 0 else 0)) % p for n in range(K)]
                images.append(solve(im, B, p))
            h = [(a + pow(l, 4, p) * z) % p for a, z in zip(act(xc, images, p), xc)]
            kk = act(h, U, p)
            f = combine(kk, B, sb + 1, p)
            if not all((f[n] == 0 for n in range(1, sb + 1, 4))):
                raise AssertionError('mathematical verification failed')
            if l == 5:
                if not all((f[2 * n] == 0 for n in [1, 3, 5])):
                    raise AssertionError('mathematical verification failed')
            extra.append({'ell': l, 'h': h, 'U_h': kk, 'IV_a_bound': sb, 'IV_b_bound': 6 if l == 5 else None})
        rr['IV'] = extra
        rr['c_orbit'] = vv
        exponents = {2: -4, 4: 8}
        cusp = {}
        for d in [1, 2, 4]:
            cusp[str(d)] = str(Fraction(4, 24 * gcd(d, 4 // d) * d) * sum((Fraction(gcd(d, delta) ** 2 * e, delta) for delta, e in exponents.items())))
        rr['F_eta_cusp_orders'] = cusp
        if not all((Fraction(v) >= 0 for v in cusp.values())):
            raise AssertionError('mathematical verification failed')
        if not sum((d * e for d, e in exponents.items())) % 24 == 0:
            raise AssertionError('mathematical verification failed')
        if not sum((4 // d * e for d, e in exponents.items())) % 24 == 0:
            raise AssertionError('mathematical verification failed')
    reports[p] = rr
    print(json.dumps(rr, sort_keys=True))
summary = {'python': sys.version, 'platform': platform.platform(), 'numpy': np.__version__, 'sympy': sp.__version__, 'elapsed_seconds': time.perf_counter() - T0, 'reports': reports, 'status': 'all assertions passed'}
Path(__file__).with_name('independent_tower_results.json').write_text(json.dumps(summary, indent=2))
print('ALL INDEPENDENT ASSERTIONS PASSED', summary['elapsed_seconds'])
