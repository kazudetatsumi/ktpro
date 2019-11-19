#!/bin/bash
export OMP_NUM_THREADS="1"
time /home/kazu/ktpro/optbinwidth01d_wholefort.py 2>&1 | tee 01d.log &
time /home/kazu/ktpro/optbinwidth2d1_wholefort.py 2>&1 | tee 2d1.log &
time /home/kazu/ktpro/optbinwidth2d2_wholefort.py 2>&1 | tee 2d2.log &
time /home/kazu/ktpro/optbinwidth2d3_wholefort.py 2>&1 | tee 2d3.log &
time /home/kazu/ktpro/optbinwidth2d4_wholefort.py 2>&1 | tee 2d4.log &
time /home/kazu/ktpro/optbinwidth2d5_wholefort.py 2>&1 | tee 2d5.log &
time /home/kazu/ktpro/optbinwidth2d6_wholefort.py 2>&1 | tee 2d6.log &
export OMP_NUM_THREADS="2"
time /home/kazu/ktpro/optbinwidth3d1_wholefort.py 2>&1 | tee 3d1.log &
time /home/kazu/ktpro/optbinwidth3d2_wholefort.py 2>&1 | tee 3d2.log &
time /home/kazu/ktpro/optbinwidth3d3_wholefort.py 2>&1 | tee 3d3.log &
time /home/kazu/ktpro/optbinwidth3d4_wholefort.py 2>&1 | tee 3d4.log &
export OMP_NUM_THREADS="8"
time /home/kazu/ktpro/optbinwidth4d-first.py 2>&1 | tee 4d-first.log &
time /home/kazu/ktpro/optbinwidth4d-second.py 2>&1 | tee 4d-second.log &
