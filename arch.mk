compfl = -O2 -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback

PF90 = h5pfc
F90 = h5pfc

FPPFLAGS += -DUSE_HDF5
LAPACK_LIB = -mkl

FIX_FORM = -fixed
FREE_FORM = -free

# C/C++ compiler
CC = icc
C++ = icpc

# C compiler options.
CFLAGS = -O2
