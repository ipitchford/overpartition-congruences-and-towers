"""Exact (cyclotomic) verification that the principal part AND constant term of
(-i(c tau+d))^{-k} G''(A_s tau) vanish at every cusp s=a/c of Gamma0(23104) with 361 | c,
where G'' = R(F''), F'' = F4 * E1 * E2 (eta exponents given by ETAQ).
Every contribution is rho * q * e(r), rho>0 with rho^2 rational, q rational, r rational.
We group by exponent E and by rho^2 squarefree class, and test the root-of-unity sum exactly
in Q(zeta_n) via reduction modulo the cyclotomic polynomial (flint fmpq_poly)."""
import sys; sys.path.insert(0, '.')
from cusps import *
import flint
from math import lcm

def sqfree_split(fr):
    # fr = s * m^2 with s squarefree-ish integer ratio; return (s_num, s_den) canonical and m
    def split(n):
        s, m, p = 1, 1, 2
        while p * p <= n:
            while n % (p * p) == 0:
                n //= p * p; m *= p
            if n % p == 0:
                n //= p; s *= p
            p += 1
        return s * n, m
    sn, mn = split(fr.numerator); sd, md = split(fr.denominator)
    # sqrt(sn/sd) = sqrt(sn*sd)/sd
    return sn * sd, Fr(mn, md * sd)

chi = {k: (1 if pow(k, 9, 19) == 1 else -1) for k in range(1, 19)}

def exact_terms(a, c):
    g_, d0, mb = egcd(a, c); b0 = -mb
    terms = {}  # (E, sqfree) -> list of (q Fraction, r Fraction)
    for j in range(8):
        for k in range(19):
            # W = (1/8) e(-3j/8) * b_k ; b_0 = 9/19 ; b_k = -(1 + chi(k) g)/38, g = sum_m chi(m) e(m/19)
            Wlist = [(Fr(9, 19 * 8), Fr(-3 * j, 8))] if k == 0 else \
                    [(Fr(-1, 38 * 8), Fr(-3 * j, 8))] + [(Fr(-chi[k] * chi[m], 38 * 8), Fr(-3 * j, 8) + Fr(m, 19)) for m in range(1, 19)]
            t = 19 * j + 8 * k
            p, q, r, s = 152 * a + t * c, 152 * b0 + t * d0, 152 * c, 152 * d0
            gam, U = decompose(p, q, r, s)
            lam = Fr(152, U[2])
            # exact series: redo eta_quot_series symbolically
            A_, B_, C_, D_ = gam
            phase = Fr(0); mag2 = Fr(1); lead = Fr(0); facs = []
            u1, u2, u3 = U
            for delta, rr in ETA_F4.items():
                (ga, gb, gc, gd), (v1, v2, v3) = decompose(delta * A_, delta * B_, C_, D_)
                phase += rr * epsilon(ga, gb, gc, gd)
                mag2 *= Fr(v1, delta) ** rr
                AA = v1 * u1; BB = v1 * u2 + v2 * u3; ZZ = v3 * u3
                phase += rr * Fr(BB, 24 * ZZ)
                lead += Fr(AA, 24 * ZZ) * rr
                facs.append((rr, AA, BB, ZZ))
            budget = -lead
            if budget < 0:
                continue
            # series terms: dict exponent -> list of (int coeff, phase)
            ser = {Fr(0): [(1, Fr(0))]}
            from math import comb
            def binom(rr, m):
                if rr >= 0: return comb(rr, m)
                return (-1) ** m * comb(-rr + m - 1, m)
            for rr, AA, BB, ZZ in facs:
                step = Fr(AA, ZZ); n = 1
                while n * step <= budget:
                    fac = {}
                    m = 0
                    while m * n * step <= budget:
                        cf = binom(rr, m) * (-1) ** m
                        if cf == 0: break
                        fac[m * n * step] = (cf, Fr(m * n * BB, ZZ))
                        m += 1
                    new = {}
                    for E1, lst in ser.items():
                        for E2_, (cf, ph) in fac.items():
                            EE = E1 + E2_
                            if EE <= budget:
                                new.setdefault(EE, []).extend((x * cf, y + ph) for x, y in lst)
                    ser = new; n += 1
            rho2 = mag2 * lam ** (2 * WT)  # rho = sqrt(mag2) * lam^WT
            sf, mfac = sqfree_split(rho2)
            for E, lst in ser.items():
                for (x, ph) in lst:
                    for (wq, wph) in Wlist:
                        terms.setdefault((E + lead, sf), []).append((wq * x * mfac, wph + ph + phase))
    return terms

def cyclo_zero(lst):
    n = 1
    for q, r in lst:
        n = lcm(n, (r % 1).denominator)
    coeffs = {}
    for q, r in lst:
        e_ = int((r % 1) * n)
        coeffs[e_] = coeffs.get(e_, 0) + q
    poly = flint.fmpq_poly([flint.fmpq(Fr(coeffs.get(i, 0)).numerator, Fr(coeffs.get(i, 0)).denominator) for i in range(n)])
    phi_n = flint.fmpq_poly(flint.fmpz_poly.cyclotomic(n).coeffs()) if hasattr(flint.fmpz_poly, 'cyclotomic') else None
    rem = poly % phi_n
    return rem == 0, n

L = 23104
allok = True
counts = [0, 0, 0]
for a, c in cusps(L):
    if c % 361: continue
    counts[0] += 1
    T = exact_terms(a, c)
    res = []
    for (E, sf), lst in sorted(T.items()):
        if E > 0: continue
        z, n = cyclo_zero(lst)
        counts[1] += 1
        counts[2] += len(lst)
        res.append((str(E), sf, len(lst), n, z))
        allok &= z
    print(f"{a}/{c}:", res, flush=True)
print("ALL principal parts and constant terms vanish exactly at 361|c cusps:", allok)
if not allok:
    raise ValueError('nonzero critical cusp coefficient')
if counts != [12, 28, 32928]:
    raise ValueError(f'incomplete cusp coverage: {counts}')
print('COVERAGE cusps/groups/summands:', counts)
