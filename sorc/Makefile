SHELL=	/bin/sh
FC      =$(FCMP)
#LDFLAGS_F  =$(NETCDF_LDFLAGS_F)
#-L/scratch1/NCEPDEV/da/Tseganeh.Gichamo/apps/lib
#LIBS    =$(LIBSM)
#LIBS = $(LDLIBS) -llapack  # -lblas"
OBJS	= ufs_land_update_snow.o scf_density.o snow_da.o snow_analysis.o driver_snowda.o

CMD     =driver_snowda
$(CMD):	$(OBJS)
	${FC} ${DEBUG} $(FFLAGS) $(LDFLG) $(OBJS) $(LDFLAGS) $(LIBSM)  -o $(CMD)

snow_da.o: snow_da.f90
	$(FC)  $(FFLAGS) $(INCS)  -c snow_da.f90
	
ufs_land_update_snow.o: ufs_land_update_snow.f90
	 $(FC)  $(FFLAGS) $(INCS) -c ufs_land_update_snow.f90

scf_density.o: scf_density.f90
	 $(FC)  $(FFLAGS) $(INCS)  -c scf_density.f90

snow_analysis.o: snow_da.o ufs_land_update_snow.o scf_density.o snow_analysis.f90
	$(FC)  $(FFLAGS) $(INCS) -c snow_analysis.f90

driver_snowda.o: snow_analysis.o driver_snowda.f90
	$(FC)  $(FFLAGS) $(INCS)  -c driver_snowda.f90

install:
	-cp -p $(CMD) ../../bin/.
clean:
	-rm -f $(OBJS) *.mod $(CMD)
