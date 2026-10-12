"""Mutate real checker parameters in isolated copies; require semantic rejection."""
from pathlib import Path
import tempfile, shutil, subprocess, sys, os, json
ROOT=Path(__file__).resolve().parent
CASES=[
 ('mod7','code/sturm_check.py','432','433'),
 ('mod17-19-23','code/mod19/exact361.py','[12, 28, 32928]','[12, 29, 32928]'),
 ('towers','code/structural/eigenform_reconstruction.py','(13, 5, 16, 7)','(13, 5, 16, 3)'),
 ('towers','code/structural/class13.py','[3,7,11]','[2,7,11]'),
]
rows=[]
for mode in ('0','1'):
 for paper,rel,old,new in CASES:
  with tempfile.TemporaryDirectory(prefix='overpartition-control-') as td:
   dst=Path(td)/paper
   shutil.copytree(ROOT/'papers'/paper,dst)
   f=dst/rel;s=f.read_text()
   if old not in s:raise ValueError('mutation target absent: '+rel)
   f.write_text(s.replace(old,new,1))
   env=dict(os.environ,PYTHONOPTIMIZE=mode,ETAQ='{"1":5,"2":19,"4":-6,"19":18,"361":-1}')
   p=subprocess.run([sys.executable,str(f)],cwd=f.parent,env=env,capture_output=True,text=True)
   if p.returncode==0 or not any(x in p.stderr for x in ('ValueError','AssertionError')):
    raise ValueError('mutation not semantically rejected: '+rel+' '+p.stderr[-500:])
   rows.append({'script':rel,'mutation':[old,new],'optimize':mode,'exitCode':p.returncode,'diagnostic':p.stderr.splitlines()[-1]})
print(json.dumps(rows,indent=2))
print('ALL EIGHT SEMANTIC NEGATIVE CONTROLS REJECTED')
