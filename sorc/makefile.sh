#!/bin/ksh
set -x

#-----------------------------------------------------
#-use standard module.
#-----------------------------------------------------
#
#export FCMP=mpiifort
export FCMP=${FCMP:-ifort}

###export DEBUG='-ftrapuv -check all -check nooutput_conversion -fp-stack-check -fstack-protector -traceback -g'
export NETCDF_INCLUDE=' -I/scratch2/NCEPDEV/nwprod/NCEPLIBS/src/netcdf_parallel2/include -I/apps/hdf5/1.10.6/intel_seq/18.0.5/include'
export INCS="-I$IP_INCd ${NETCDF_INCLUDE}"
export FFLAGS="$INCS -O3 -fp-model precise -r8 -convert big_endian -heap-arrays -g -traceback" 
export OMPFLAG=-qopenmp
export LDFLG=-qopenmp

export NETCDF_LDFLAGS_F='-L/apps/hdf5/1.10.6/intel_seq/18.0.5/lib -L/scratch2/NCEPDEV/nwprod/NCEPLIBS/src/netcdf_parallel2/lib -lnetcdff'
export LIBSM="${W3NCO_LIBd} ${BACIO_LIB4} ${IP_LIBd} ${SP_LIBd} ${NETCDF_LDFLAGS_F}"

make -f Makefile clean
make -f Makefile
make -f Makefile install
make -f Makefile clean
