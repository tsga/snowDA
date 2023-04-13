Program ResampleForcing

    USE NETCDF

    Use, Intrinsic :: IEEE_ARITHMETIC	
    
    IMPLICIT NONE
    ! include 'mpif.h'
	
    integer, parameter :: dp = kind(1.d0)

    REAL, ALLOCATABLE     :: Src_Lat(:), Src_Lon(:), Dest_Lat(:), Dest_Lon(:)
    Integer, ALLOCATABLE  :: Index_Source_atDest(:)

    INTEGER             :: ERROR, NCID, NCID_Dest, ID_DIM, ID_VAR, indx, io, NDIM, &
                           NDIM_i, NDIM_o, iv
    Integer             :: n_forc
    LOGICAL             :: file_exists
    CHARACTER(LEN=500)  :: namelist_name, Src_Static_File, Dest_Static_File, DIM_NAME, TIM_NAME, &
                           Src_Lat_Name, Src_Lon_Name, Dest_Lat_Name, Dest_Lon_Name, &
                           Src_Forc_File, Dest_Forc_File    !,forc_inp_path

    Real, allocatable   :: Src_forcArray(:, :), Dest_forcArray(:, :)

    character(len=*), parameter, dimension (8)   :: forc_var_list = &
          (/"precipitation_bilinear", "solar_radiation",   &
            "longwave_radiation", "temperature",           &
            "wind_speed", "specific_humidity",             &
            "surface_pressure", "precipitation_conserve"/) 

    character(len=*), parameter, dimension (8)   :: dest_forc_var_list = &
            (/"precipitation", "solar_radiation",   &
              "longwave_radiation", "temperature",           &
              "wind_speed", "specific_humidity",             &
              "surface_pressure", "precipitation_conserve"/) 

    ! integer            :: header_buffer_val = 16384

    REAL      :: latlon_tol != 0.5

    NAMELIST /NAMSNO/ Src_Static_File, Dest_Static_File, DIM_NAME, TIM_NAME, &
            Src_Lat_Name, Src_Lon_Name, Dest_Lat_Name, Dest_Lon_Name, &
            n_forc, Src_Forc_File, Dest_Forc_File, latlon_tol

    DATA namelist_name/""/
    DATA Src_Static_File/"/scratch2/BMC/gsienkf/Tseganeh.Gichamo/ufsoi/static/ufs-land_C96_static_fields.nc"/ 
    DATA Dest_Static_File/"/scratch2/BMC/gsienkf/Tseganeh.Gichamo/ufsoi/test2/restart/SNOTEL/ufs-land_C96_static_fields.nc"/ 
    DATA DIM_NAME/"location"/
    DATA TIM_NAME/"time"/
    DATA Src_Lat_Name/"latitude"/
    DATA Src_Lon_Name/"longitude"/
    DATA Dest_Lat_Name/"latitude"/
    DATA Dest_Lon_Name/"longitude"/
    Data n_forc/7/
    Data latlon_tol/0.5/
    ! Data forc_inp_path/"./"/
    Data Src_Forc_File/"/scratch2/BMC/gsienkf/Tseganeh.Gichamo/ufsoi/forcing/C96_GDAS_forcing_2019-12-15.nc"/ 
    Data Dest_Forc_File/"./snotel_forcing_2015-12-15.nc"/ 
    
    print *, "starting resampling forcing"

    call get_command_argument(1, namelist_name)
    if(namelist_name == "") then 
            print *,  "no namelist given. Using default values "
            ! stop 10  
    endif

    open(30, file=namelist_name, form="formatted")
    read(30, NAMSNO)
    close(30)

    WRITE(6, NAMSNO)
   
    INQUIRE(FILE=trim(Src_Static_File), EXIST=file_exists)
    if (.not. file_exists) then
            print *, 'ResampleForcing error, file does not exist', &
                    trim(Src_Static_File), ' exiting'
            stop 10
    endif
    ERROR=NF90_OPEN(TRIM(Src_Static_File),NF90_NOWRITE,NCID)
    CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(Src_Static_File) )    
    ERROR=NF90_INQ_DIMID(NCID, TRIM(DIM_NAME), ID_DIM)
    CALL NETCDF_ERR(ERROR, 'ERROR READING src Dimension' )    
    ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM_o)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Size of src Dimension' )
      
    ALLOCATE(Src_Lat(NDIM_o))
    ALLOCATE(Src_Lon(NDIM_o)) 
    
    ERROR=NF90_INQ_VARID(NCID, TRIM(Src_Lat_Name), ID_VAR)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Src_Lat_Name ID' )
    ERROR=NF90_GET_VAR(NCID, ID_VAR, Src_Lat)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Src_Lat RECORD' )
    ERROR=NF90_INQ_VARID(NCID, TRIM(Src_Lon_Name), ID_VAR)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Src_Lon_Name ID' )
    ERROR=NF90_GET_VAR(NCID, ID_VAR, Src_Lon)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Src_Lon RECORD' )
    ! ERROR=NF90_INQ_VARID(NCID, TRIM(Eval_VAR_NAME), ID_VAR)
    ! CALL NETCDF_ERR(ERROR, 'ERROR READING Index_Obs_EvalPts ID' )
    ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Index_Obs_EvalPts)
    ! CALL NETCDF_ERR(ERROR, 'ERROR READING Index_Obs_EvalPts RECORD' )
    ERROR = NF90_CLOSE(NCID)
    CALL NETCDF_ERR(ERROR, 'ERROR closing src file' )

    INQUIRE(FILE=trim(Dest_Static_File), EXIST=file_exists)
    if (.not. file_exists) then
            print *, 'ResampleForcing error, file does not exist', &
                    trim(Dest_Static_File), ' exiting'
            stop 10
    endif
    ERROR=NF90_OPEN(TRIM(Dest_Static_File),NF90_NOWRITE,NCID)
    CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(Dest_Static_File) )    
    ERROR=NF90_INQ_DIMID(NCID, TRIM(DIM_NAME), ID_DIM)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Dest Dimension' )    
    ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM_i)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dest Dimension' )
      
    ALLOCATE(Dest_Lat(NDIM_i))
    ALLOCATE(Dest_Lon(NDIM_i)) 
    ALLOCATE(Index_Source_atDest(NDIM_i)) 

    ERROR=NF90_INQ_VARID(NCID, TRIM(Dest_Lat_Name), ID_VAR)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Dest_Lat_Name ID' )
    ERROR=NF90_GET_VAR(NCID, ID_VAR, Dest_Lat)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Dest_Lat RECORD' )
    ERROR=NF90_INQ_VARID(NCID, TRIM(Dest_Lon_Name), ID_VAR)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Dest_Lon_Name ID' )
    ERROR=NF90_GET_VAR(NCID, ID_VAR, Dest_Lon)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Dest_Lon RECORD' )
    ! ERROR=NF90_INQ_VARID(NCID, TRIM(Eval_VAR_NAME), ID_VAR)
    ! CALL NETCDF_ERR(ERROR, 'ERROR READING Index_Obs_EvalPts ID' )
    ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Index_Obs_EvalPts)
    ! CALL NETCDF_ERR(ERROR, 'ERROR READING Index_Obs_EvalPts RECORD' )
    ERROR = NF90_CLOSE(NCID)
    CALL NETCDF_ERR(ERROR, 'ERROR closing Dest file' )
    
    Where(Src_Lon > 180) Src_Lon = Src_Lon - 360
    Where(Dest_Lon > 180) Dest_Lon = Dest_Lon - 360

    ! print *, "Src_Lon ", Src_Lon
    ! print *
    ! print *, "Dest_Lon", Dest_Lon
    ! print *
    ! print *, "Src_Lat ", Src_Lat
    ! print *
    ! print *, "Dest_Lat", Dest_Lat

    Index_Source_atDest = -999

    Do indx=1, NDIM_i
        Do io = 1, NDIM_o    
            ! if((.NOT. IEEE_IS_NAN(SND_GHCND(io))) .AND. &
            if((abs(Src_Lat(io) - Dest_Lat(indx)) < latlon_tol) .AND. &
               (abs(Src_Lon(io) - Dest_Lon(indx)) < latlon_tol)  &
              ) then            
                Index_Source_atDest(indx) = io
                Exit
            endif
            if (io .eq. NDIM_o) print*, "no matching values to dest indx = ", indx 
        End do
    Enddo

    print *, "done finding resampling indices"

    INQUIRE(FILE=trim(Src_Forc_File), EXIST=file_exists)   
    if (.not. file_exists) then
            print *, 'ResampleForcing error, file does not exist', &
                    trim(Src_Forc_File), ' exiting'
            stop 10
    endif
    ERROR=NF90_OPEN(TRIM(Src_Forc_File), NF90_NOWRITE, NCID)
    CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(Src_Forc_File) ) 
    ERROR=NF90_INQ_DIMID(NCID, TRIM(TIM_NAME), ID_DIM)
    CALL NETCDF_ERR(ERROR, 'ERROR READING src time Dimension' )    
    ERROR=NF90_INQUIRE_DIMENSION(NCID, ID_DIM, LEN=NDIM)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Size of src time Dimension' )   
   
    ALLOCATE(Src_forcArray(NDIM_o, NDIM))
    ALLOCATE(Dest_forcArray(NDIM_i, NDIM))

    INQUIRE(FILE=trim(Dest_Forc_File), EXIST=file_exists)
    if (.not. file_exists) then
            print *, 'ResampleForcing error, file does not exist', &
                    trim(Dest_Forc_File) , ' exiting'
            stop 10
    endif
    ERROR=NF90_OPEN(TRIM(Dest_Forc_File),NF90_WRITE, NCID_Dest)
    CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(Dest_Forc_File) )  
    
    Do iv=1, n_forc
        print *, "writing "//TRIM(forc_var_list(iv)) 

        Dest_forcArray = -999

        ERROR=NF90_INQ_VARID(ncid, TRIM(forc_var_list(iv)), ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING src '//TRIM(forc_var_list(iv))//' ID' )
        ERROR=NF90_GET_VAR(ncid, ID_VAR, Src_forcArray)
        CALL NETCDF_ERR(ERROR, 'ERROR READING src RECORD for '//TRIM(forc_var_list(iv)) )

        Do indx=1, NDIM_i
            if(Index_Source_atDest(indx) > 0) then 
                Dest_forcArray(indx, :) = Src_forcArray(Index_Source_atDest(indx), :)
            Endif
        Enddo
        
        ERROR=NF90_INQ_VARID(NCID_Dest, TRIM(dest_forc_var_list(iv)), ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING dest '//TRIM(dest_forc_var_list(iv))//' ID' )
        ERROR = nf90_put_var(NCID_Dest, ID_VAR, Dest_forcArray)  !, start = (/1/), count = (/NDIM/))
        CALL NETCDF_ERR(ERROR, 'ERROR Writing Dest RECORD '//TRIM(dest_forc_var_list(iv)) ) 
    Enddo

    ERROR = NF90_CLOSE(NCID)
    CALL NETCDF_ERR(ERROR, 'ERROR closing '//TRIM(Src_Forc_File) )
    ERROR = NF90_CLOSE(NCID_Dest)
    CALL NETCDF_ERR(ERROR, 'ERROR closing '//TRIM(Dest_Forc_File) )

    DEALLOCATE(Src_Lat, Src_Lon, Dest_Lat, Dest_Lon, Index_Source_atDest)
    DEALLOCATE(Src_forcArray, Dest_forcArray) 

    print*, "Done!"

        
End Program ResampleForcing

SUBROUTINE NETCDF_ERR( ERR, STRING )
    
    !--------------------------------------------------------------
    ! IF AT NETCDF CALL RETURNS AN ERROR, PRINT OUT A MESSAGE
    ! AND STOP PROCESSING.
    !--------------------------------------------------------------
    USE NETCDF
        IMPLICIT NONE
    
        ! include 'mpif.h'
    
        INTEGER, INTENT(IN) :: ERR
        CHARACTER(LEN=*), INTENT(IN) :: STRING
        CHARACTER(LEN=80) :: ERRMSG
    
        IF( ERR == NF90_NOERR )RETURN
        ERRMSG = NF90_STRERROR(ERR)
        PRINT*,''
        PRINT*,'FATAL ERROR: ', TRIM(STRING), ': ', TRIM(ERRMSG)
        PRINT*,'STOP.'
    
        RETURN
END SUBROUTINE NETCDF_ERR