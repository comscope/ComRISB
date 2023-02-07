include ../../../../Makefile.in

DCPP = -Dconstraint_s2

#FFLAGS = -O0 -check bound -check pointer -check arg_temp_created \
		-traceback 
#DCPP = -Ddebug_mode

obj = GPREC.o GTIME.o GCONSTANT.o GHDF5.o GPRIMME.o \
	  UTIL.o SPARSE.o GHDF5_SM.o GSPCI.o GSPCI_LIB.o

LIB = ../../../lib/libprimme.a
INCLUDE=-I../../../include

all: ${obj} MAIN.o
	${PF90} ${FFLAGS} -cpp -Dreal_version -o exe_spcir MAIN.o \
		${obj} ${HDF5_LIB} ${MKL_LIB} ${LIB}

GPRIMME.o: ../../GPRIMME.f90
	${F90} ${FIX_FORM} ${FFLAGS} -cpp -Dreal_version \
			-c $< ${HDF5_INCLUDE} $(INCLUDE)

%.o: ../%.f90
	${PF90} ${FREE_FORM} ${FFLAGS} -cpp ${DCPP} -Dreal_version \
			-c $< ${HDF5_INCLUDE} $(INCLUDE)

%.o: ../../%.f90
	${PF90} ${FREE_FORM} ${FFLAGS} -cpp  ${DCPP} -Dreal_version \
			-c $< ${HDF5_INCLUDE} $(INCLUDE)

# special fix for ifort
GSPCI_LIB.o: ../../GSPCI_LIB.f90
	${PF90} ${FREE_FORM} ${FFLAGS} -O0 -cpp  ${DCPP} -Dreal_version \
			-c $< ${HDF5_INCLUDE} $(INCLUDE)

clean:
	rm -f *o *mod *out

clean_all:
	rm -f *o *mod *out exe_spcir

install:
	cp exe_spcir ${DESTDIR}
