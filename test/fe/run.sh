#!/bin/bash

mkdir -p dft && cd dft
${COMRISB_BIN}/../ComBin/cif2matdelab.py ../bcc.cif -k 3
mpirun -np 8 ${COMRISB_BIN}/rspflapw.exe
cd ..
mkdir -p dftg/u0j0/lowh
cp comrisb.ini dftg/u0j0/
cp bcc.cif ginit.json dftg/u0j0/lowh/.
cd dftg/u0j0/lowh
${COMRISB_BIN}/init_grisb.py -u eV -s 1
cd ..
${COMRISB_BIN}/comrisb.py -c
