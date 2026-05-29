#!/bin/bash
cd build
source ../../../../tools/benchmark/.venv/bin/activate
python ../../../../tools/check_speedup.py -lab_path ../ -num_runs 3 > ../speedup.txt
