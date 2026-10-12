"""Failure-sensitive verification entry point; run from any directory.

Each subprocess runs in its own source directory; checks stay active under
PYTHONOPTIMIZE=1. Logs and exact versions are retained under verification/.
"""
from pathlib import Path
import os, sys, subprocess, time, json, hashlib, platform
import importlib.metadata

ROOT = Path(__file__).resolve().parent
PLANS = {
    'mod7': ['code/ghn_check.py', 'code/verify_identities.py', 'code/sturm_check.py',
             'code/independent/reduced_sturm_stdlib.py'],
    'mod17-19-23': ['code/hecke1723/check.py', 'code/independent/independent_checks.py',
                    'code/mod19/check_weights.py', 'code/mod19/exact361.py',
                    'code/mod19/sturm_check.py'],
    'towers': ['code/referee/referee_check.py',
               'code/structural/independent_tower_check.py',
               'code/structural/eigenform_reconstruction.py', 'code/structural/class13.py'],
}

def main():
    choice = sys.argv[1] if len(sys.argv)>1 else 'all'
    if choice not in (*PLANS, 'all'):
        raise SystemExit('use: python verify.py [all|mod7|mod17-19-23|towers]')
    dest = ROOT/'verification'
    dest.mkdir(exist_ok=True)
    result = {'python':sys.version, 'platform':platform.platform(),
              'optimised':bool(sys.flags.optimize), 'commands':[],
              'versions':{x:importlib.metadata.version(x) for x in ['python-flint','sympy','mpmath','numpy']}}
    for name in (PLANS if choice=='all' else [choice]):
        for rel in PLANS[name]:
            script = ROOT/'papers'/name/rel
            env = dict(os.environ, ETAQ='{"1":5,"2":19,"4":-6,"19":18,"361":-1}',
                       OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1')
            log = dest/(name+'-'+script.stem+('-opt' if sys.flags.optimize else '')+'.log')
            started=time.monotonic()
            print('RUN',name,rel,flush=True)
            with log.open('w') as f:
                p=subprocess.run([sys.executable,str(script)],cwd=script.parent,env=env,
                                 stdout=f,stderr=subprocess.STDOUT)
            result['commands'].append({'paper':name,'script':rel,
                'sha256':hashlib.sha256(script.read_bytes()).hexdigest(),
                'exitCode':p.returncode,'seconds':time.monotonic()-started,
                'log':str(log.relative_to(ROOT)), 'logSha256':hashlib.sha256(log.read_bytes()).hexdigest()})
            (dest/(choice+('-opt' if sys.flags.optimize else '')+'.json')).write_text(json.dumps(result,indent=2)+'\n')
            if p.returncode:
                print(log.read_text()[-4000:]); raise SystemExit(p.returncode)
    print('ALL REQUESTED OBLIGATIONS PASSED',flush=True)

if __name__=='__main__': main()
