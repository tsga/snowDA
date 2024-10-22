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
#module use $ufs_mods
#module load ufs_hera.intel
module load intel/2022.1.2  impi/2022.1.2
module load netcdf-hdf5parallel/4.7.4
#module load bacio/2.0.3 w3emc/2.4.0  #bacio/2.4.1 w3emc/2.10.0
#module load ip/3.0.2 sp/2.0.3  # ip/4.3.0 sp/2.3.3 
module load nco/5.1.6
module list

export FCMP=mpiifort
#export FCMP=${FCMP:-ifort}

###export DEBUG='-ftrapuv -check all -check nooutput_conversion -fp-stack-check -fstack-protector -traceback -g'
export NETCDF_INCLUDE=-I$NETCDF/include
export INCS=$NETCDF_INCLUDE  #$IP_INCS"
export FFLAGS="-O0 -fp-model precise -r8 -init=huge  -convert big_endian -heap-arrays -g -traceback" 
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
