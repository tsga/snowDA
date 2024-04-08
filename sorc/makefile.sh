#! /bin/sh -l
set -ex

#-----------------------------------------------------
#-use standard module.
#-----------------------------------------------------
#
#export jedi_mods=/scratch1/NCEPDEV/da/Tseganeh.Gichamo/GDASApp/modulefiles
#module use $jedi_mods
##module use /scratch1/NCEPDEV/da/Tseganeh.Gichamo/GDASApp/modulefiles
#module load GDAS/hera
#module list

export ufs_mods=/scratch1/NCEPDEV/da/Tseganeh.Gichamo/global-workflow/sorc/ufs_model.fd/modulefiles

module purge
module use $ufs_mods
module load ufs_hera.intel
module load netcdf-hdf5parallel/4.7.4
module list

export FCMP=mpiifort
#export FCMP=${FCMP:-ifort}

###export DEBUG='-ftrapuv -check all -check nooutput_conversion -fp-stack-check -fstack-protector -traceback -g'
export NETCDF_INCLUDE=-I{$NETCDF}/include
export INCS=${NETCDF_INCLUDE} #$IP_INCS
export FFLAGS="-O3 -fp-model precise -r8 -init=huge  -convert big_endian -heap-arrays -g -traceback" 
#use -O0 for debug
export OMPFLAG=-qopenmp
export LDFLG=-qopenmp
export LDFLAGS="-L${NETCDF}/lib -L/scratch1/NCEPDEV/da/Tseganeh.Gichamo/APPS/lib"
export LDLIBS="-lnetcdf -lnetcdff -llapack -lblas"  # -lmpifort"
export LIBSM="${W3EMC_LIBd} ${BACIO_LIB4} ${IP_LIBd} ${SP_LIBd} ${LDLIBS}"   #NETCDF_LDFLAGS_F}" ${IP_LIBd}
#export LIBSM="${LDLIBS}" 

make -f Makefile clean
make -f Makefile
make -f Makefile install
make -f Makefile clean
