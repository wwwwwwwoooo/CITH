"""Synchronize v7.3 modular sources, inclusion fragments, full source."""
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
MAIN=ROOT/'manuscript'
SUP=ROOT/'supplement'

def expand_inputs(text: str, directory: Path) -> str:
    def replace(m: re.Match) -> str:
        p=directory/m.group(1)
        if not p.suffix:p=p.with_suffix('.tex')
        return '\n% BEGIN expanded '+p.name+'\n'+expand_inputs(p.read_text(),p.parent)+'\n% END expanded '+p.name+'\n'
    return re.sub(r'\\input\{([^}]+)\}',replace,text)

def subsection(body: str, label: str) -> str:
    pos=body.index(r'\label{'+label+'}')
    start=body.rfind(r'\subsection{',0,pos)
    after=re.search(r'\n\\(?:subsection|section)\{',body[pos:])
    end=pos+after.start() if after else len(body)
    return body[start:end].strip()+'\n'

def main() -> None:
    pre=(MAIN/'manuscript.tex').read_text()
    body=(MAIN/'body.tex').read_text()
    full=expand_inputs(pre,MAIN).replace(r'\graphicspath{{figures/}{../figures/}}',
                                          r'\graphicspath{{manuscript/figures/}{figures/}}')
    full=full.replace(r'\bibliographystyle{cithplainnat}','')
    full=full.replace(r'\bibliography{references}',(MAIN/'manuscript.bbl').read_text())
    assert r'\input{' not in full and r'\bibliography{' not in full
    (ROOT/'CITH_V7_3_full.tex').write_text(full)
    abstract=pre.split(r'\begin{abstract}',1)[1].split(r'\par\smallskip',1)[0].strip()+'\n'
    for p in [ROOT/'new_abstract.tex',SUP/'abstract.tex']:p.write_text(abstract)
    sections=ROOT/'sections';sections.mkdir(exist_ok=True)
    for p in sections.glob('*.tex'):p.unlink()
    chunks={}
    matches=list(re.finditer(r'^\\section\*?\{([^}]+)\}',body,re.M))
    for i,m in enumerate(matches):
        end=matches[i+1].start() if i+1<len(matches) else len(body)
        slug=re.sub(r'[^a-z0-9]+','_',m.group(1).lower()).strip('_')
        chunk=body[m.start():end]
        chunks[i+1]=chunk
        (sections/f'{i+1:02d}_{slug}.tex').write_text(chunk)
    (SUP/'section4_replacement.tex').write_text(chunks[4])
    (SUP/'section7_replacement.tex').write_text(chunks[7])
    (SUP/'conclusion_replacement.tex').write_text(chunks[8])
    production=subsection(body,'sec:production')
    (SUP/'production_addition.tex').write_text(production)
    persistence=subsection(body,'prop:persistence')
    curvature=subsection(body,'prop:curvature')
    (SUP/'A1_A2_proofs.tex').write_text(persistence+'\n'+curvature)
    start=body.index('The canonical unit points')
    end=body.index('The parameter order is',start)
    (SUP/'generator_specification.tex').write_text(body[start:end])
    (ROOT/'new_proofs_full.tex').write_text(
        '% Inclusion fragment. Compile CITH_V7_3_full.tex for complete context.\n'
        +curvature+'\n'+persistence+'\n'+production)
    cor=(SUP/'v73_production_corollary.tex').read_text()
    (ROOT/'CITH_V7_3_new_corollary.tex').write_text(cor)
    print('Synchronized full source, sections, abstract and proofs.')

if __name__=='__main__':main()
