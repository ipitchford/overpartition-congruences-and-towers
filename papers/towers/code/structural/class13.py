"""Complete twelve-state certificate for whole-class vanishing; no assert gates."""
import json
from pathlib import Path
from math import isqrt
import flint

root=Path(__file__).parent
r=json.loads((root/'independent_tower_results.json').read_text())['reports']['13']
U=r['U_exact']; v=[1,4,1]; orbit=[]
for i in range(12):
    if v in orbit: raise ValueError('early orbit return')
    orbit.append(v)
    v=[sum(v[a]*U[a][b] for a in range(3))%13 for b in range(3)]
if v!=orbit[0]: raise ValueError('no twelve-state return')
N=193
ph=[0]*N;ph[0]=1
for k in range(1,isqrt(N-1)+1): ph[k*k]=2
F=[0]*N
for d in range(1,N,2):
    for n in range(d,N,2*d):F[n]+=d
ph=flint.nmod_poly(ph,13);F=flint.nmod_poly(F,13)
B=[F.pow_trunc(a,N).mul_low(ph.pow_trunc(11-4*a,N),N) for a in range(3)]
zeros=[];witness=[]
for i,v in enumerate(orbit):
    f=sum((b*x for b,x in zip(B,v)),flint.nmod_poly([],13))
    vals=[int(f[n]) for n in range(7,N,8)]
    if not any(vals):zeros.append(i)
    witness.append(vals[0])
    if vals[0]!=(7+4*pow(5,i,13))%13:raise ValueError('witness formula failed')
if zeros!=[3,7,11]:raise ValueError('wrong complete vanishing states')
out={'U_rows':U,'orbit':orbit,'bound_inclusive':192,'zero_states':zeros,
     'witness_at_7':witness,'status':'passed'}
(root/'class13.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
