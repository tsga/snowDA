Program BlackListEval

    IMPLICIT NONE

    USE NETCDF

    Use, Intrinsic :: IEEE_ARITHMETIC	

    include 'mpif.h'
	
	integer, parameter :: dp = kind(1.d0)

    INTEGER               :: NDIM
    REAL, ALLOCATABLE)    :: SND_GHCND(:)
    Integer, ALLOCATABLE  :: Index_Obs_EvalPts(:)

    INTEGER                :: ERROR, NCID, grp_ncid, ID_DIM, ID_VAR, indx, NDIM, NDIM_i
    LOGICAL                :: file_exists
    CHARACTER(LEN=500)  :: namelist_name, STN_OBS_File, STN_DIM_NAME, STN_VAR_NAME, &
                           Eval_OBS_File, Eval_DIM_NAME, Eval_VAR_NAME

    NAMELIST /NAMSNO/ STN_OBS_File, STN_DIM_NAME, STN_VAR_NAME, &
        Eval_OBS_File, Eval_DIM_NAME, Eval_VAR_NAME

    DATA namelist_name/""/
    DATA STN_OBS_File/'/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND_IODA/ghcn_snwd_ioda_20191215.nc'/ 
    DATA STN_DIM_NAME/"nlocs"/
    DATA STN_VAR_NAME/"totalSnowDepth"/
    DATA Eval_OBS_File/"/scratch2/BMC/gsienkf/Tseganeh.Gichamo/ufsoi/test2/restart/ENSRF_SND_5pEV/SNOANLEnSRF.20191027_18.nc"/
    DATA Eval_DIM_NAME/"eval_points"/
    DATA Eval_VAR_NAME/"Index_Obs_EvalPts"/

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
            print *, 'Observation_Read_GHCND_Tile_excNaN erro,,file does not exist', &
                    trim(Eval_OBS_File) , ' exiting'
            call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
    endif
    ERROR=NF90_OPEN(TRIM(Eval_OBS_File),NF90_NOWRITE,NCID)
    CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(Eval_OBS_File) )    
    ERROR=NF90_INQ_DIMID(NCID, TRIM(Eval_DIM_NAME), ID_DIM)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Eval Dimension' )    
    ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM_i)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Eval Dimension' )
    
    ALLOCATE(Index_Obs_EvalPts(NDIM_i))    

    ERROR=NF90_INQ_VARID(NCID, TRIM(Eval_VAR_NAME), ID_VAR)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Index_Obs_EvalPts ID' )
    ERROR=NF90_GET_VAR(NCID, ID_VAR, Index_Obs_EvalPts)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Index_Obs_EvalPts RECORD' )
    ERROR = NF90_CLOSE(NCID)
    CALL NETCDF_ERR(ERROR, 'ERROR closing file with eval index' )

    INQUIRE(FILE=trim(STN_OBS_File), EXIST=file_exists)
    if (.not. file_exists) then
            print *, 'Observation_Read_GHCND_Tile_excNaN erro,,file does not exist', &
                    trim(STN_OBS_File) , ' exiting'
            call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
    endif
    ERROR=NF90_OPEN(TRIM(STN_OBS_File),NF90_WRITE,NCID)
    CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(STN_OBS_File) )    
    ERROR=NF90_INQ_DIMID(NCID, TRIM(STN_DIM_NAME), ID_DIM)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension' )    
    ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM)
    CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension' )

    ALLOCATE(SND_GHCND(NDIM))    

    ERROR=nf90_inq_ncid(ncid, "ObsValue", grp_ncid)
    CALL NETCDF_ERR(ERROR, 'ERROR ObsValue ID' )
    ERROR=NF90_INQ_VARID(grp_ncid, TRIM(STN_VAR_NAME), ID_VAR)
    CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD ID' )
    ERROR=NF90_GET_VAR(grp_ncid, ID_VAR, SND_GHCND)
    CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD RECORD' )

    Do indx=1, NDIM_i
        SND_GHCND(Index_Obs_EvalPts(indx)) = -999.9
    Enddo

    ERROR = nf90_put_var(grp_ncid, ID_VAR, SND_GHCND)  !, start = (/1/), count = (/NDIM/))
    CALL NETCDF_ERR(ERROR, 'ERROR Writing SNWD RECORD' ) 

    ERROR = NF90_CLOSE(NCID)
    CALL NETCDF_ERR(ERROR, 'ERROR closing SND File' )

    DEALLOCATE(Index_Obs_EvalPts)
    DEALLOCATE(SND_GHCND)
                
    RETURN
        
End Program BlackListEval