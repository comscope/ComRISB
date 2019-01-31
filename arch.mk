# USE_HDF5 =true


#compfl = -debug -g -CB -check bounds -traceback -check uninit -fp-model precise
compfl = -O3

PF90 = mpiifort
F90 = ifort

ifdef USE_HDF5
    FPPFLAGS += -DUSE_HDF5
    PF90 = h5pfc
endif

LAPACK_LIB = -mkl

#### CTQMC ######################

BASE_CPPFLAGS = -DNDEBUG
BASE_LIBS = -lm

CXX = CC
CXX_MPI = CC -DHAVE_MPI
CXXFLAGS = -std=c++11 -fexceptions -Wall -O3

#### ComRISB ######################
# will always use hdf5 internally, the interface, however, as of now, 
# only works if the others are compiled without USE_HDF5.
PF90_risb = h5pfc

FIX_FORM = -fixed
FREE_FORM = -free

# C++ compiler
C++ = CC

# C compiler options.
CFLAGS = -O2
