#i&&&&&&&&&&&&&&&&&&&&&&&&&&&&        UNIX MPICH MPI       &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
include ../arch.mk

wan90_dir = ../wannier90_2.1

comp = $(PF90) $(FPPFLAGS) -c $(compfl)
link = $(PF90)  $(compfl)
LIB = $(LAPACK_LIB) -L$(wan90_dir) -lwannier

.SUFFIXES: .mod .F .o 

o_files = \
  comwann_mod.o \
 beging.o \
   bessr.o \
   dqdall.o \
 cross.o \
 cft.o \
 diagonalization_lapack.o \
 det3.o \
 dif.o \
 fac2.o \
 fft3.o \
 fifi_j_prepare.o \
 integral_band_pw_band_mt.o \
 integral_band_pw_band.o \
 output_u_f0f2f4f6.o \
 size_shift_par.o \
 tcgstb.o \
 unport.o \
getlimj.o \
 gbs_number.o \
 lget.o \
 prd1.o \
 prd2.o \
 rotmat_wan_proj.o \
 rotate.o \
 sphharm.o \
 solid.o \
 sym_z_0.o \
 sym_a_0.o \
 wan_pb_rs.o \
 wannier_mmn.o \
 wannier_amn_from_win.o \
 wannier_interface.o \
 qft.o \
 init.o \
 start_par_solid.o \
 read_restart_file.o \
 set_memory.o \
 read_wan_bni.o \
 read_wan_ini.o \
 

mod_files = \
 comwann_mod.mod

ifdef USE_HDF5
  o_files +=  read_wan_bni_hdf5.o
  o_files +=  read_restart_file_hdf5.o
  o_files +=  write_wannier_dat_hdf5.o
  o_files +=  hdf5io_mod.o
  mod_files += hdf5io_mod.mod
endif	

TARGET = ComWann

$(TARGET): $(mod_files) $(o_files)
	$(link) $(o_files) -o $(TARGET) $(LIB)
	mkdir -p ../bin
	cp $(TARGET) ../bin

.f90.o:
	$(comp) $<
.F.o:
	$(comp) $<

%.mod: %.F
	$(comp) $<

%.mod: %.f90
	$(comp) $<
clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf ComWann

include ./make.depends
