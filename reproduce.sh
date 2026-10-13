#!/bin/sh
set -eu
cd "$(dirname "$0")"
python code/run_all_checks.py
python code/rerun_archive.py
