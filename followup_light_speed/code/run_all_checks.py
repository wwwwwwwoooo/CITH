"""Run the twelve CITH-LS validators and preserve literal stdout."""
from pathlib import Path
import json,re,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]

def main():
 results=ROOT/'results';results.mkdir(exist_ok=True)
 logs=[];summary={};pcount=fcount=0
 for script in sorted((ROOT/'code').glob('validate_ls_*.py')):
  proc=subprocess.run([sys.executable,str(script)],cwd=ROOT.parent,text=True,capture_output=True,timeout=900)
  out=proc.stdout+proc.stderr
  (results/(script.stem+'_stdout.txt')).write_text(out)
  logs.append(f'=== {script.name}; exit={proc.returncode} ===\n'+out)
  counts=[(int(p),int(f)) for p,f in re.findall(r'PASS=(\d+) FAIL=(\d+)',out)]
  p=sum(x[0] for x in counts);f=sum(x[1] for x in counts);pcount+=p;fcount+=f
  summary[script.name]={'exit_code':proc.returncode,'pass':p,'fail':f}
  print(script.name,proc.returncode,p,f,flush=True)
  if proc.returncode or f:raise RuntimeError(f'Validation failure in {script.name}; see saved stdout')
 (results/'ALL_LS_CHECKS_OUTPUT.txt').write_text('\n'.join(logs))
 payload={'scripts':summary,'script_count':len(summary),'total_assertions_pass':pcount,'total_assertions_fail':fcount,'interpretation':'Counts are assertions, not independent random draws.'}
 (results/'ALL_LS_CHECKS_SUMMARY.json').write_text(json.dumps(payload,indent=2));print('TOTAL',pcount,fcount)
if __name__=='__main__':main()
