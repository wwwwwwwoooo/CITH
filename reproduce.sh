#!/bin/sh
set -eu
cd "$(dirname "$0")"
mkdir -p results
python code/validate_additions.py > results/additions_stdout.json
python code/validate_v73.py > results/v73_checks_stdout.json
python code/audit_archived_inputs.py > results/archive_audit_stdout.json
python code/make_v73_phase.py > results/v73_phase_stdout.json
printf '%s\n' 'Numerical reproduction completed. Run python code/build_v73.py to rebuild PDFs.'
