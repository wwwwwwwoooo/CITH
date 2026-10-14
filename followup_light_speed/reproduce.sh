#!/bin/sh
set -eu
cd "$(dirname "$0")"
python code/run_all_checks.py
python code/make_figures.py
for i in 1 2 3; do pdflatex -interaction=nonstopmode -halt-on-error CITH_LS_draft.tex > /dev/null; done
python code/audit_draft.py
