Program BlackListEval

    USE NETCDF

    Use, Intrinsic :: IEEE_ARITHMETIC	
    
    IMPLICIT NONE
    ! include 'mpif.h'
	
    integer, parameter :: dp = kind(1.d0)

    REAL, ALLOCATABLE     :: SND_GHCND(:), Lat_EvalPts(:), Lon_EvalPts(:), &
                             STN_Lat(:), STN_Lon(:)
    Integer, ALLOCATABLE  :: Index_Obs_EvalPts(:)

    INTEGER             :: ERROR, NCID, grp_ncid, ID_DIM, ID_VAR, indx, io, NDIM, NDIM_i
    LOGICAL             :: file_exists
    CHARACTER(LEN=500)  :: namelist_name, STN_OBS_File, STN_DIM_NAME, STN_VAR_NAME, &
                           Eval_OBS_File, Eval_DIM_NAME, &
                           Eval_Lat_NAME, Eval_Lon_NAME    !Eval_VAR_NAME, 

    Integer             :: dim_eval, id_indxobseval
    integer            :: header_buffer_val = 16384

    REAL,Parameter      :: latlon_tol = 0.00001

    NAMELIST /NAMSNO/ STN_OBS_File, STN_DIM_NAME, STN_VAR_NAME, &
        Eval_OBS_File, Eval_DIM_NAME, Eval_Lat_NAME, Eval_Lon_NAME

    DATA namelist_name/""/
    DATA STN_OBS_File/"ghcn_snwd_ioda_20191215.nc"/ 
    DATA STN_DIM_NAME/"nlocs"/
    DATA STN_VAR_NAME/"totalSnowDepth"/
    DATA Eval_OBS_File/"obs_exclude_list_2019.nc"/
    DATA Eval_DIM_NAME/"eval_points"/
    ! DATA Eval_VAR_NAME/"Index_Obs_EvalPts"/
    DATA Eval_Lat_NAME/"latitude"/
    DATA Eval_Lon_NAME/"longitude"/

    call get_command_argument(1, namelist_name)
    if(namelist_name == "") then 
            print *,  "add namelist to the command line: "
            stop 10  
    endif

    open(30, file=namelist_name, form="formatted")
    read(30, NAMSNO)
    close(30)
   
    INQUIRE(FILE=trim(Eval_OBS_File), EXIST=file_exists)
    if (.not. file_exists) then
            print *, 'Blacklist obs error, file does not exist', &
                    trim(Eval_OBS_File), ' exiting'
            stop 10
    endif
    ERROR=NF90_OPEN(TRIM(Eval_OBS_File),NF90_WRITE,NCID)
    CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(Eval_OBS_File) )    
    ERROR=NF90_INQ_DIMID(NCID, TRIM(Eval_DIM_NAME), ID_DIM)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Eval Dimension' )    
    ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM_i)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Eval Dimension' )
    
    ALLOCATE(Index_Obs_EvalPts(NDIM_i))   
    ALLOCATE(Lat_EvalPts(NDIM_i))
    ALLOCATE(Lon_EvalPts(NDIM_i)) 
    
    ERROR=NF90_INQ_VARID(NCID, TRIM(Eval_Lat_NAME), ID_VAR)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Eval_Lat_NAME ID' )
    ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_EvalPts)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Eval_Lat_NAME RECORD' )
    ERROR=NF90_INQ_VARID(NCID, TRIM(Eval_Lon_NAME), ID_VAR)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Eval_Lon_NAME ID' )
    ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_EvalPts)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Eval_Lon_NAME RECORD' )
    ! ERROR=NF90_INQ_VARID(NCID, TRIM(Eval_VAR_NAME), ID_VAR)
    ! CALL NETCDF_ERR(ERROR, 'ERROR READING Index_Obs_EvalPts ID' )
    ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Index_Obs_EvalPts)
    ! CALL NETCDF_ERR(ERROR, 'ERROR READING Index_Obs_EvalPts RECORD' )
    ERROR = NF90_CLOSE(NCID)
    CALL NETCDF_ERR(ERROR, 'ERROR closing file with eval index' )

    INQUIRE(FILE=trim(STN_OBS_File), EXIST=file_exists)
    if (.not. file_exists) then
            print *, 'Observation_Read_GHCND_Tile_excNaN erro,,file does not exist', &
                    trim(STN_OBS_File) , ' exiting'
            stop 10
    endif
    ERROR=NF90_OPEN(TRIM(STN_OBS_File),NF90_WRITE,NCID)
    CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(STN_OBS_File) )  
    
    ! define mode
    ERROR=NF90_REDEF(NCID)
    CALL NETCDF_ERR(ERROR, 'ERROR define mode' )  
     !eval obs indx 
    error = nf90_def_dim(ncid, 'obs_excluded', NDIM_i, dim_eval)
    call netcdf_err(error, 'DEFINING obs_excluded' )
    ! index of Obs at eval points 
    error = nf90_def_var(ncid, 'Index_Obs_Excluded', NF90_INT, dim_eval, id_indxobseval)
    call netcdf_err(error, 'DEFINING Index_Obs_Excluded' )
    error = nf90_put_att(ncid, id_indxobseval, "long_name", "Index of Observations Excluded")
    call netcdf_err(error, 'DEFINING Index_Obs_Excluded LONG NAME')
    error = nf90_put_att(ncid, id_indxobseval, "units", "-")
    call netcdf_err(error, 'DEFINING Index_Obs_Excluded UNITS')

    error = nf90_enddef(ncid, header_buffer_val,4,0,4)
    call netcdf_err(error, 'End DEFINE HEADER' )

    ERROR=NF90_INQ_DIMID(NCID, TRIM(STN_DIM_NAME), ID_DIM)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension' )    
    ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension' )

    ALLOCATE(SND_GHCND(NDIM))    
    ALLOCATE(STN_Lat(NDIM))  
    ALLOCATE(STN_Lon(NDIM))  
    
    ERROR=nf90_inq_ncid(ncid, "MetaData", grp_ncid)
    CALL NETCDF_ERR(ERROR, 'ERROR MetaData ID' )
    ERROR=NF90_INQ_VARID(grp_ncid, TRIM(Eval_Lat_NAME), ID_VAR)
    CALL NETCDF_ERR(ERROR, 'ERROR READING STN_Lat ID' )
    ERROR=NF90_GET_VAR(grp_ncid, ID_VAR, STN_Lat)
    CALL NETCDF_ERR(ERROR, 'ERROR READING STN_Lat RECORD' )
    ERROR=NF90_INQ_VARID(grp_ncid, TRIM(Eval_Lon_NAME), ID_VAR)
    CALL NETCDF_ERR(ERROR, 'ERROR READING STN_Lon ID' )
    ERROR=NF90_GET_VAR(grp_ncid, ID_VAR, STN_Lon)
    CALL NETCDF_ERR(ERROR, 'ERROR READING STN_Lon RECORD' )

    ERROR=nf90_inq_ncid(ncid, "ObsValue", grp_ncid)
    CALL NETCDF_ERR(ERROR, 'ERROR ObsValue ID' )
    ERROR=NF90_INQ_VARID(grp_ncid, TRIM(STN_VAR_NAME), ID_VAR)
    CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD ID' )
    ERROR=NF90_GET_VAR(grp_ncid, ID_VAR, SND_GHCND)
    CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD RECORD' )
    
    Index_Obs_EvalPts = -999

    Do indx=1, NDIM_i
        Do io = 1, NDIM    
            if((.NOT. IEEE_IS_NAN(SND_GHCND(io))) .AND. &
               (abs(Lat_EvalPts(indx) - STN_Lat(io)) < latlon_tol) .AND. &
               (abs(Lon_EvalPts(indx) - STN_Lon(io)) < latlon_tol)  &
              ) then            
                Index_Obs_EvalPts(indx) = io
                !print*, "no matching values!"
                cycle
            endif
        End do
    Enddo

    Do indx=1, NDIM_i
        if(Index_Obs_EvalPts(indx) > 0) then 
            SND_GHCND(Index_Obs_EvalPts(indx)) = -999.9
        Endif
    Enddo

    ERROR = nf90_put_var(grp_ncid, ID_VAR, SND_GHCND)  !, start = (/1/), count = (/NDIM/))
    CALL NETCDF_ERR(ERROR, 'ERROR Writing SNWD RECORD' ) 

    error = nf90_put_var(ncid, id_indxobseval, Index_Obs_EvalPts)
    call netcdf_err(error, 'WRITING Index_Obs_Excluded RECORD' )

    ERROR = NF90_CLOSE(NCID)
    CALL NETCDF_ERR(ERROR, 'ERROR closing SND File' )

    DEALLOCATE(Index_Obs_EvalPts)
    DEALLOCATE(SND_GHCND)
    DEALLOCATE(Lat_EvalPts, Lon_EvalPts, STN_Lat, STN_Lon) 
        
End Program BlackListEval

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