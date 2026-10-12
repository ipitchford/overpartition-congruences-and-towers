# Explicit checks remain active with PYTHONOPTIMIZE=1.
"""Reconstruct exact U(4) eigenforms underlying every tower correction."""
from pathlib import Path
from math import gcd, lcm, isqrt
from fractions import Fraction
import json, time, sympy as s, numpy as np
T0 = time.perf_counter()
root = Path(__file__).parent
r = json.loads((root / 'independent_tower_results.json').read_text())['reports']
out = {}
for p, lam, alpha, cls in [(13, 5, 16, 7), (17, 7, -64, 3), (23, 10, -512, 5)]:
    d = r[str(p)]
    U = s.Matrix(d['U_exact'])
    v = (U.T - alpha * s.eye(len(d['U_exact']))).nullspace()[0]
    den = lcm(*[z.q for z in v])
    v = [int(z * den) for z in v]
    gg = gcd(*v)
    v = [z // gg for z in v]
    if next((z for z in v if z)) < 0:
        v = [-z for z in v]
    if not list(U.T * s.Matrix(v)) == [alpha * z for z in v]:
        raise AssertionError('mathematical verification failed')
    scale = next((h * pow(z % p, -1, p) % p for h, z in zip(d['cases'][0]['h'], v) if z % p))
    if not [scale * z % p for z in v] == d['cases'][0]['h']:
        raise AssertionError('mathematical verification failed')
    N = 32 * (lam + 1) + 1
    ph = [0] * N
    ph[0] = 1
    for k in range(1, isqrt(N - 1) + 1):
        ph[k * k] = 2
    F = [sum((k for k in range(1, n + 1) if n % k == 0)) if n % 2 else 0 for n in range(N)]

    def mul(a, b):
        return list(np.convolve(np.array(a, dtype=object), np.array(b, dtype=object))[:N])

    def pw(a, k):
        v = [1] + [0] * (N - 1)
        for _ in range(k):
            v = mul(v, a)
        return v
    B = [mul(pw(F, a), pw(ph, 2 * lam + 1 - 4 * a)) for a in range(len(v))]
    H = [sum((v[a] * B[a][n] for a in range(len(v)))) for n in range(N)]
    witnesses = [(n, int(H[n])) for n in range(cls, N, 8) if H[n]]
    if not not witnesses:
        raise AssertionError('mathematical verification failed')
    result = {'lambda': lam, 'alpha_integer': alpha, 'alpha_mod_p': alpha % p, 'primitive_integer_eigenvector': v, 'h_mod_scale': scale, 'projection_class': cls, 'exact_projection_bound': N - 1, 'nonzero_coefficients': witnesses, 'characteristic_polynomial': str(s.factor(U.charpoly(s.Symbol('x')).as_expr()))}
    if p in (13, 17):
        exponents = {1: 2, 2: 7, 4: 2} if p == 13 else {1: -6, 2: 27, 4: -6}
        result['eta_exponents'] = exponents
        result['cusp_orders'] = {str(d): str(Fraction(4, 24 * gcd(d, 4 // d) * d) * sum((Fraction(gcd(d, e) ** 2 * a, e) for e, a in exponents.items()))) for d in (1, 2, 4)}
        if not all((Fraction(x) >= 0 for x in result['cusp_orders'].values())):
            raise AssertionError('mathematical verification failed')
        # Integral eta-product coefficients from q E'/E, including its shift.
        shift = sum(e * a for e, a in exponents.items()) // 24
        if sum(e * a for e, a in exponents.items()) != 24 * shift:
            raise ValueError('nonintegral eta shift')
        logder = [0] * N
        for e, a in exponents.items():
            for k in range(e, N, e):
                for j in range(k, N, k):
                    logder[j] -= a * k
        coeff = [1] + [0] * (N - 1)
        for n in range(1, N):
            total = sum(logder[k] * coeff[n-k] for k in range(1, n+1))
            if total % n:
                raise ValueError('nonintegral eta coefficient')
            coeff[n] = total // n
        eta = [0] * shift + coeff[:N-shift]
        if eta != H:
            raise ValueError('eta identity mismatch')
        result['eta_identity_coefficients_checked'] = N
    out[p] = result
    print(json.dumps(result), flush=True)
(root / 'exact_correction_eigenvectors.json').write_text(json.dumps(out, indent=2))
print('ALL EXACT EIGENFORM CHECKS PASSED', time.perf_counter() - T0)
