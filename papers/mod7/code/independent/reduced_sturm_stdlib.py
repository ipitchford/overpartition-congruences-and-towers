# Explicit checks remain active with PYTHONOPTIMIZE=1.
"""Independent stdlib finite check for the sharper twisting level 3136.

This verifies a reviewer-proposed shorter proof, not a replacement of the
submitted full-bound replay. Sparse theta factors are multiplied directly;
no FLINT, SymPy, fast polynomial routines, or submitted scripts are imported.
"""
from pathlib import Path
from math import isqrt
import hashlib
import json
import sys
import time
OUT = Path(__file__).resolve().parent
N = 3136
B = 8064
LENGTH = B + 1

def multiply_theta(coefficients, d, sign):
    result = coefficients.copy()
    for k in range(1, isqrt(B // d) + 1):
        degree = d * k * k
        factor = -2 if sign < 0 and k % 2 else 2
        for n in range(LENGTH - degree):
            result[n + degree] += factor * coefficients[n]
    return [v % 7 for v in result]

def build_z():
    result = [1] + [0] * B
    for d, sign, exponent in [(1, -1, 5), (7, 1, 4), (14, 1, 6), (28, 1, 12), (56, 1, 1)]:
        for _ in range(exponent):
            result = multiply_theta(result, d, sign)
    return result
start = time.monotonic()
z = build_z()
monomials = []
for positive_exponent, negative_exponent in [(8, 0), (4, 4), (0, 8)]:
    term = z.copy()
    for _ in range(positive_exponent):
        term = multiply_theta(term, 1, 1)
    for _ in range(negative_exponent):
        term = multiply_theta(term, 1, -1)
    monomials.append(term)
H = [(a - b + c) % 7 for a, b, c in zip(*monomials)]
relevant = [n for n in range(LENGTH) if n % 56 in (11, 43, 51)]
nonzero = [n for n in relevant if H[n]]
control = [n for n in range(LENGTH) if n % 8 == 3 and n % 7 in (3, 5, 6) and H[n]]
if not N == 2 ** 6 * 7 ** 2:
    raise AssertionError('mathematical verification failed')
if not B == 18 * (N * 3 // 2 * 8 // 7) // 12:
    raise AssertionError('mathematical verification failed')
if not len(relevant) == 432:
    raise AssertionError('mathematical verification failed')
if not not nonzero:
    raise AssertionError('mathematical verification failed')
if not (H[0] == 1 and control):
    raise AssertionError('mathematical verification failed')
coefficient_file = OUT / 'H-mod7-through-8064.bin'
coefficient_file.write_bytes(bytes(H))
receipt = {'status': 'passed', 'role': 'independent arithmetic for reviewer-proposed sharper level', 'level': N, 'weight': 18, 'index': 5376, 'bound_inclusive': B, 'target_count': len(relevant), 'nonzero_target_indices': nonzero, 'nonzero_control_count': len(control), 'H0': H[0], 'coefficient_encoding': 'byte n is coefficient of q^n modulo 7', 'coefficient_sha256': hashlib.sha256(bytes(H)).hexdigest(), 'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'elapsed_seconds': time.monotonic() - start, 'python': sys.version}
OUT.joinpath('reduced_sturm_stdlib.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt, indent=2))
