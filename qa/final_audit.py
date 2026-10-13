"""Audit the assembled source, count main text, verify frozen files, split PDFs.

Requires PyMuPDF for PDF work and texcount on PATH for word counting.
It validates mechanical properties, not mathematical truth.
"""
from pathlib import Path
import hashlib,json,re,subprocess,sys
import fitz
R=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 source=(R/'CITH_V8_full.tex').read_text()
 body=source.split('% BEGIN MAIN TEXT\n',1)[1].split('% END MAIN TEXT',1)[0]
 abstract=re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}',source,re.S).group(1)
 (R/'qa/abstract_for_count.tex').write_text(abstract)
 def count(path):
  text=subprocess.check_output(['texcount','-brief',str(path)],text=True)
  m=re.search(r'(\d+)\+(\d+)\+(\d+)',text)
  if not m:raise RuntimeError(text)
  return {'text':int(m[1]),'headers':int(m[2]),'captions':int(m[3]),'literal_output':text.strip()}
 wc=count(R/'latex/main_text.tex');ac=count(R/'qa/abstract_for_count.tex')
 labels=re.findall(r'\\label\{([^}]+)\}',source)
 refs=[x for m in re.finditer(r'\\(?:[cC]ref|ref|eqref)\*?\{([^}]+)\}',source) for x in m[1].split(',')]
 old=(R/'archive/CITH_V7_3_full.tex').read_text()
 oldlabels=set(re.findall(r'\\label\{([^}]+)\}',old))
 log=(R/'CITH_V8_full.log').read_text()
 warnings=[s for s in log.splitlines() if re.search(r'Warning:|pdfTeX warning|Overfull|undefined references|multiply defined',s,re.I)]
 aux=(R/'CITH_V8_full.aux').read_text()
 m=re.search(r'\\newlabel\{page:supplement-start\}\{\{[^}]*\}\{(\d+)\}',aux)
 if not m:raise RuntimeError('Missing supplement start label')
 start=int(m[1])-1
 pdf=fitz.open(R/'CITH_V8_full.pdf')
 assert 'Supplementary Material' in pdf[start].get_text()
 for name,lo,hi in [('CITH_V8_Manuscript.pdf',0,start-1),('CITH_V8_Supplement.pdf',start,len(pdf)-1)]:
  out=fitz.open();out.insert_pdf(pdf,from_page=lo,to_page=hi);out.set_metadata({'title':name.replace('_',' ').replace('.pdf',''),'author':'Hyeonjun Joo'});out.save(R/name);out.close()
 pages=[]
 for i,p in enumerate(pdf):
  text=p.get_text(); outside=[]
  for b in p.get_text('blocks'):
   if b[0]<-1 or b[1]<-1 or b[2]>p.rect.width+1 or b[3]>p.rect.height+1:outside.append(b[:4])
  pages.append({'page':i+1,'characters':len(text),'text_outside_page':outside})
 # Check figures against the original archive. No archive files are modified here.
 figs={f'fig{i}.pdf':sha(R/f'figures/fig{i}.pdf')==sha(R/f'archive/v7_3/manuscript/figures/fig{i}.pdf') for i in range(1,5)}
 # Statement index makes each declared result discoverable for human review.
 statements=[]
 for m in re.finditer(r'\\begin\{(proposition|lemma|corollary|productioncorollary)\}(?:\[([^\]]*)\])?(.*?)\\end\{\1\}',source,re.S):
  ids=re.findall(r'\\label\{([^}]+)\}',m[3]);statements.append({'environment':m[1],'title':m[2],'labels':ids})
 report={'main_word_count':wc,'abstract_word_count':ac,'labels':len(labels),'reference_uses':len(refs),'duplicate_labels':sorted({x for x in labels if labels.count(x)>1}),'unresolved_references':sorted(set(refs)-set(labels)),'old_labels_not_retained':sorted(oldlabels-set(labels)),'latex_warnings':warnings,'pdf_pages':len(pdf),'manuscript_pages':start,'supplement_pages':len(pdf)-start,'source_figure_hash_matches':figs,'statement_count':len(statements),'statements':statements,'page_checks':pages,'scope':'Mechanical checks and word counts; mathematical review recorded separately.'}
 (R/'results/final_audit.json').write_text(json.dumps(report,indent=2))
 (R/'results/word_count.txt').write_text('Main text: '+wc['literal_output']+'\nAbstract: '+ac['literal_output']+'\nMain includes only sections 1-9; excludes abstract, declarations, appendices, references, supplement.\n')
 assert 9000<=wc['text']<=11000
 assert ac['text']<=250
 assert not report['duplicate_labels'] and not report['unresolved_references'] and not warnings
 assert all(figs.values())
 assert not any(p['text_outside_page'] for p in pages)
 print(json.dumps({k:v for k,v in report.items() if k not in ['statements','page_checks']},indent=2))
if __name__=='__main__':main()
