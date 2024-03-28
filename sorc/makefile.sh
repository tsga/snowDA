#! /bin/sh -l
set -x

#-----------------------------------------------------
#-use standard module.
#-----------------------------------------------------
#
export jedi_mods=/scratch1/NCEPDEV/da/Tseganeh.Gichamo/GDASApp/modulefiles
module use $jedi_mods
#module use /scratch1/NCEPDEV/da/Tseganeh.Gichamo/GDASApp/modulefiles
module load GDAS/hera
module list

export FCMP=mpiifort
#export FCMP=${FCMP:-ifort}

###export DEBUG='-ftrapuv -check all -check nooutput_conversion -fp-stack-check -fstack-protector -traceback -g'
export NETCDF_INCLUDE=' -I/scratch2/NCEPDEV/nwprod/NCEPLIBS/src/netcdf_parallel2/include -I//scratch1/NCEPDEV/nems/role.epic/spack-stack/spack-stack-1.6.0/envs/unified-env-rocky8/install/intel/2021.5.0/hdf5-1.14.0-lixiejp/include'  # -I/apps/oneapi/mpi/2021.5.1/include'
export INCS=${NETCDF_INCLUDE}   #$IP_INCS
export FFLAGS="$INCS -O3 -fp-model precise -r8 -convert big_endian -heap-arrays -g -traceback" 
export OMPFLAG=-qopenmp
export LDFLG=-qopenmp
export NETCDF_LDFLAGS_F='-L/scratch1/NCEPDEV/nems/role.epic/spack-stack/spack-stack-1.6.0/envs/unified-env-rocky8/install/intel/2021.5.0/hdf5-1.14.0-lixiejp/include/lib -L/scratch2/NCEPDEV/nwprod/NCEPLIBS/src/netcdf_parallel2/lib' #-L/apps/oneapi/mpi/2021.5.1/lib'
export LDLIBS="-lnetcdff -llapack -lblas"  # -lmpifort"
export LIBSM="${W3EMC_LIBd} ${BACIO_LIB4} ${SP_LIBd} ${LDLIBS}"   #NETCDF_LDFLAGS_F}" ${IP_LIBd}
#export LIBSM="${LDLIBS}" 

make -f Makefile clean
make -f Makefile
#make -f Makefile install
#make -f Makefile clean
