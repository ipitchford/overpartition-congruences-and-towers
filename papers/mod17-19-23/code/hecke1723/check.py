# Independent checks: (a) class descriptions; (b) Sturm-range vanishing of h on the 8-class with bounds;
# (c) direct overpartition check mod p for small n via flint.
import flint
N=4000
def series(p):
    phi=[0]*N; k=0
    while k*k<N: phi[k*k]+=1 if k==0 else 2; k+=1
    F=[0]*N
    for n in range(1,N,2): F[n]=sum(d for d in range(1,n+1) if n%d==0)
    R=lambda a:flint.nmod_poly([x%p for x in a],p)
    return R(phi),R(F)
def trunc(f): return flint.nmod_poly(f.coeffs()[:N],f.modulus())
def pw(f,e):
    r=flint.nmod_poly([1],f.modulus())
    for _ in range(e): r=trunc(r*f)
    return r
def coef(f,n):
    c=f.coeffs(); return int(c[n]) if n<len(c) else 0
QR=lambda a,l: 0 if a%l==0 else (1 if pow(a,(l-1)//2,l)==1 else -1)
for (p,ell,lam,r,cls,hc,sgn) in [(17,11,7,3,[19,35,43,51,83],[0,12,12,0],-1),(23,13,10,5,[29,53,61,69,77,101],[0,16,16,14,7,15],1)]:
    A=8*ell
    exp=sorted(m for m in range(A) if m%8==r and QR(m,ell)==sgn)
    print(p,"classes",exp,"==",cls,exp==cls)
    if exp != cls:
        raise ValueError('incorrect residue classes')
    phi,F=series(p)
    h=flint.nmod_poly([0],p)
    for a,x in enumerate(hc): h=h+trunc(pw(F,a)*pw(phi,2*lam+1-4*a))*x
    g=flint.nmod_poly([coef(h,n) if n%8==r else 0 for n in range(N)],p)
    gphi=trunc(g*phi)
    bound=(lam+1)*384//12
    print(" Sturm bound (weight %d, Gamma0(256)) = %d; g*phi zero up to %d:"%(lam+1,bound,N-1), all(coef(gphi,n)==0 for n in range(N)))
    if N <= bound or any(coef(gphi,n) for n in range(N)):
        raise ValueError('class projection failed')
