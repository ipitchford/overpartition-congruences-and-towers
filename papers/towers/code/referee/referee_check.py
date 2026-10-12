# Explicit checks remain active with PYTHONOPTIMIZE=1.
"""Independent standard-library verification at mathematically sufficient bounds."""
from math import isqrt

def mul(a, b, N, p=None):
    out = [0] * N
    nz = [(j, x) for j, x in enumerate(b) if x]
    for i, x in enumerate(a):
        if not x:
            continue
        for j, y in nz:
            if i + j >= N:
                break
            out[i + j] += x * y
    return [x % p for x in out] if p else out

def basis(lam, N, p=None):
    phi = [0] * N
    phi[0] = 1
    for j in range(1, isqrt(N - 1) + 1):
        phi[j * j] = 2
    F = [0] * N
    for d in range(1, N, 2):
        for j in range(d, N, 2 * d):
            F[j] += d
    powers = [[1] + [0] * (N - 1)]
    for j in range(2 * lam + 1):
        powers.append(mul(powers[-1], phi, N, p))
    fp = powers[0]
    B = []
    for a in range((2 * lam + 1) // 4 + 1):
        B.append(mul(fp, powers[2 * lam + 1 - 4 * a], N, p))
        if a + 1 < (2 * lam + 1) // 4 + 1:
            fp = mul(fp, F, N, p)
    return (B, phi, powers)

def coords(f, B, p=None):
    r = list(f)
    x = []
    for a in range(len(B)):
        v = r[a]
        x.append(v)
        r = [u - v * b for u, b in zip(r, B[a])]
        if p:
            r = [u % p for u in r]
    return x

def comb(v, B, N, p=None):
    ans = [sum((x * b[n] for x, b in zip(v, B))) for n in range(N)]
    return [x % p for x in ans] if p else ans

def apply(v, A, p):
    return [sum((v[a] * A[a][b] for a in range(len(v)))) % p for b in range(len(v))]

def leg(n, l):
    v = pow(n % l, (l - 1) // 2, l)
    return -1 if v == l - 1 else v

def main(p, lam, xc, cases):
    K = len(xc)
    bound = 32 * (lam + 1)
    N = max(bound + 1, max((l * l * (K - 1) + 1 for l, _, _ in cases)))
    B, phi, _ = basis(lam, N, p)
    weight = (p - 1) // 2
    sb = weight // 2
    M = p * sb + 1
    ph = [0] * M
    ph[0] = 1
    for j in range(1, isqrt(M - 1) + 1):
        ph[j * j] = 2
    power = [1] + [0] * (M - 1)
    for _ in range(p - 1):
        power = mul(power, ph, M, p)
    product = mul(comb(xc, B, sb + 1, p), phi[:sb + 1], sb + 1, p)
    if not product == [power[p * n] for n in range(sb + 1)]:
        raise AssertionError('mathematical verification failed')
    print('p', p, 'initial identity PASS through integral-weight Sturm bound', sb)
    ub = lam + 1
    EB, ep, _ = basis(lam, 4 * ub + 1)
    U = []
    for b in EB:
        image = [b[4 * n] for n in range(ub + 1)]
        v = coords(image, EB)
        if not comb(v, EB, ub + 1) == image:
            raise AssertionError('mathematical verification failed')
        U.append(v)
    print('exact U4 rows', U, 'PASS through bound', ub)
    for l, chi, classes in cases:
        T = []
        for b in B:
            im = [(b[l * l * n] + leg((-1) ** lam * n, l) * pow(l, lam - 1, p) * b[n] + (pow(l, 2 * lam - 1, p) * b[n // (l * l)] if n % (l * l) == 0 else 0)) % p for n in range(K)]
            T.append(coords(im, B, p))
        h = [(a - chi * pow(l, lam - 1, p) * c) % p for a, c in zip(apply(xc, T, p), xc)]
        print('l', l, 'T', T, 'h', h)
        if classes is not None:
            orbit = []
            v = h
            while v not in orbit:
                orbit.append(v)
                v = apply(v, U, p)
            return_index = orbit.index(v)
            for v in orbit:
                f = comb(v, B, bound + 1, p)
                if not all((f[n] == 0 for n in range(bound + 1) if n % 8 in classes)):
                    raise AssertionError('mathematical verification failed')
            print('orbit length', len(orbit), 'return', return_index, 'projection PASS through', bound)
        else:
            k = apply(h, U, p)
            f = comb(k, B, bound + 1, p)
            if not all((f[n] == 0 for n in range(bound + 1) if n % 4 == 1)):
                raise AssertionError('mathematical verification failed')
            print('IV-a k', k, 'PASS through', bound)
            if l == 5:
                if not all((f[2 * n] == 0 for n in range(1, 7, 2))):
                    raise AssertionError('mathematical verification failed')
                print('IV-b w PASS through bound 6')
main(17, 7, [1, 13, 13, 0], [(11, 1, {3}), (29, -1, {3})])
main(23, 10, [1, 9, 5, 14, 17, 20], [(13, 1, {5})])
main(13, 5, [1, 4, 1], [(7, 1, {7}), (5, -1, None), (29, -1, None)])
print('ALL CERTIFYING CHECKS PASSED')
