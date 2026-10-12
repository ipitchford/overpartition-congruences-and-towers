# exact check: sum_t W_t e(n t/152) = [n=3 mod 8]*[(n/19)=1] for all n mod 152 (exact cyclotomic arithmetic)
import sys; sys.path.insert(0,'.')
from fractions import Fraction as Fr
import flint
chi={k:(1 if pow(k,9,19)==1 else -1) for k in range(1,19)}
phi=flint.fmpq_poly(flint.fmpz_poly.cyclotomic(152).coeffs())
ok=True
for n in range(152):
    co=[Fr(0)]*152
    for j in range(8):
        for k in range(19):
            t=19*j+8*k
            W=[(Fr(9,152),Fr(-3*j,8))] if k==0 else [(Fr(-1,304),Fr(-3*j,8))]+[(Fr(-chi[k]*chi[m],304),Fr(-3*j,8)+Fr(m,19)) for m in range(1,19)]
            for q,r in W:
                e=int(((r+Fr(n*t,152))%1)*152); co[e]+=q
    p=flint.fmpq_poly([flint.fmpq(x.numerator,x.denominator) for x in co])%phi
    want=1 if (n%8==3 and n%19 and chi[n%19]==1) else 0
    ok&=(p==want)
print("indicator weights exact:",ok)
if not ok:
    raise ValueError('projection weights do not select the intended coefficients')
