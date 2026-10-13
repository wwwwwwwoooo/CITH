"""Rerun the four archived numerical commands in a disposable working copy."""
from pathlib import Path
import shutil,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parents[1]
def main():
 with tempfile.TemporaryDirectory(prefix='cith_v73_') as td:
  work=Path(td)/'v7_3';shutil.copytree(ROOT/'archive/v7_3',work)
  for name in ['validate_additions.py','validate_v73.py','audit_archived_inputs.py','make_v73_phase.py']:
   p=subprocess.run([sys.executable,str(work/'code'/name)],cwd=work,text=True,capture_output=True,timeout=240)
   if p.returncode:raise RuntimeError(name+'\n'+p.stderr)
   (work/'results'/(Path(name).stem+'_rerun_stdout.txt')).write_text(p.stdout+p.stderr)
   print(name,'PASS',flush=True)
  target=ROOT/'results/archive_rerun';shutil.copytree(work/'results',target,dirs_exist_ok=True)
if __name__=='__main__':main()
