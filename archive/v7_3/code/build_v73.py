"""Build article, line-numbered article, supplement, and expanded-source check."""
from pathlib import Path
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]

def compile_tex(directory: Path, stem: str, passes: int=3, outdir: Path|None=None) -> None:
    logdir=ROOT/'results/build_logs';logdir.mkdir(exist_ok=True)
    label=stem if directory==ROOT else directory.name+'_'+stem
    combined=[]
    for _ in range(passes):
        cmd=['pdflatex','-interaction=nonstopmode','-halt-on-error','-file-line-error']
        if outdir:
            outdir.mkdir(exist_ok=True,parents=True);cmd+=['-output-directory',str(outdir)]
        cmd+=[stem+'.tex']
        result=subprocess.run(cmd,cwd=directory,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        combined.append(result.stdout)
        (logdir/(label+'.build.log')).write_text('\n'.join(combined))
        if result.returncode:
            raise RuntimeError('TeX build failed: '+label+'\n'+result.stdout[-3500:])
    built=outdir if outdir else directory
    shutil.copy2(built/(stem+'.log'),logdir/(label+'.log'))


def main() -> None:
    if not shutil.which('pdflatex'):
        raise RuntimeError('pdflatex is required to build the requested PDFs.')
    # Bibliography is unchanged; preserve and use the supplied .bbl files.
    for stem in ['manuscript','review']:
        compile_tex(ROOT/'manuscript',stem)
    compile_tex(ROOT/'supplement','supplement')
    subprocess.run([sys.executable,str(ROOT/'code/export_release.py')],check=True)
    compile_tex(ROOT,'CITH_V7_3_full',outdir=ROOT/'results/full_build')
    for source,dest in [('manuscript/manuscript.pdf','CITH_V7_3_revised.pdf'),
                        ('manuscript/review.pdf','CITH_V7_3_line_numbered.pdf'),
                        ('supplement/supplement.pdf','CITH_V7_3_supplement.pdf')]:
        shutil.copy2(ROOT/source,ROOT/dest)
    print('Built all PDFs and compiled the expanded LaTeX source.')

if __name__=='__main__':main()
