#!/bin/bash

mkdir -p work && cd work
python ../scan_semicirc.py --skipshow
python ../scan_semicirc.py --mu --skipshow

# 84.79user 288.73system 7:22.89elapsed 84%CPU (0avgtext+0avgdata 132344maxresident)k
# 776inputs+115816outputs (0major+35129707minor)pagefaults 0swaps
