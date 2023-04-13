Program ims_scfi

    IMPLICIT NONE

    USE NETCDF

    Use, Intrinsic :: IEEE_ARITHMETIC	

    include 'mpif.h'
	
	integer, parameter :: dp = kind(1.d0)

    Integer   :: IVEGSRC, land_type, grid_type
	INTEGER   :: IY, IM, ID, IH, IDIM, JDIM, NUM_TILES, MAX_TASKS
    Integer   :: NPROCS, MYRANK, LENSFC

	INTEGER			     :: IERR
	
	CHARACTER(len=250)   :: ims_inp_file, ims_inp_file_indices
	CHARACTER(len=250) 	 :: out_file
	CHARACTER(len=4)     :: y_str, m_str, d_Str, h_str, fvs_tile
    Character(LEN=3)     :: rank_str
    logical              :: fv3_index, print_debg_info

    
    
    ! for mpi par
    INTEGER   :: mp_start, mp_end, N_sA, N_sA_Ext
    Integer   :: LENSFC_proc
 
    Real, Allocatable	:: SNCOV_IMS(:), SNO_IMS_at_Grid(:)  
    Real, Allocatable	:: RLA(:), RLO(:), OROG(:)
    Integer             :: num_subgrd_ims_cels
    CHARACTER(LEN=500)  :: IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, &
                           vector_restart_prefix, static_prefix, output_path, &          

    NAMELIST/NAMSNO/ IDIM, JDIM, NUM_TILES, &
    IY, IM, ID, IH, MAX_TASKS    &
    IVEGSRC, land_type, grid_type, fv3_index, &
    num_subgrd_ims_cels, IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH,  &
    vector_restart_prefix, static_prefix, output_path, &
    print_debg_info, PRINTRANK
!   
    DATA land_type/2/
    Data grid_type/2/
    Data IVEGSRC/2/
    DATA IDIM,JDIM,NUM_TILES/96,96,6/
    Data LENSFC/18360/
    DATA IY,IM,ID,IH,FH/2019,12,15,18/
    DATA MAX_TASKS/99999/
    DATA num_subgrd_ims_cels/627/   !30  ! (max) number of ims subcells within a tile grid cell
    DATA IMS_SNOWCOVER_PATH/'./'/
    DATA IMS_INDEXES_PATH/'./'/
    DATA print_debg_info/.false./
    Data fv3_index/.false./
    DATA PRINTRANK/4/
    Data vector_restart_prefix/"./"/
    DATA static_prefix/"./"/
    DATA output_path/"./"/

    CALL MPI_INIT(IERR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, IERR)

    PRINT*,"STARTING Snow analysis on RANK "

    CALL BAOPENR(360, "imsfsca.nml", IERR)             !"snowDA.nml"   !
    READ(360, NML=NAMSNO)
    IF (MYRANK==PRINTRANK) WRITE(6, NAMSNO)
    
    if (grid_type == 1) then
        if (LENSFC <> IDIM * JDIM) then
            print*, "inconsistent grid sizes"
            stop 
        endif
    endif

    N_sA = LENSFC/NPROCS    !/ Np_til  ! sub array length per proc
    N_sA_Ext = MOD(LENSFC, NPROCS)  !LENSFC - N_sA * Np_til ! extra grid cells
    if(MYRANK < N_sA_Ext) then   !p_tRank== 0) then 
        !  mp_start = 1
        mp_start = MYRANK * (N_sA + 1) + 1 
        mp_end = mp_start + N_sA   
        LENSFC_proc = N_sA + 1
    else
        mp_start = MYRANK * N_sA + N_sA_Ext + 1 
        mp_end = mp_start + N_sA - 1    
        LENSFC_proc = N_sA 
    endif

    Allocate(SNO_IMS_at_Grid(LENSFC_proc))
    Allocate(SNCOV_IMS(LENSFC_proc))
    Allocate(RLA(LENSFC_proc))
    Allocate(RLO(LENSFC_proc))
    Allocate(OROG(LENSFC_proc))
   
    Call get_ims_depth_scf(MYRANK, IDIM, JDIM, IY, IM, ID, IH, IVEGSRC, &
            num_subgrd_ims_cels, IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, & 
            vector_restart_prefix, static_prefix, fv3_index, print_debg_info, &
            mp_start, mp_end, LENSFC_proc, &
            RLA, RLO, OROG, SNCOV_IMS, SNO_IMS_at_Grid)
	
	Write(rank_str, '(I3.3)') (MYRANK+1)
    
    write(y_str, "(I4)") IY
    write(m_str, "(I0.2)") IM
    write(d_str, "(I0.2)") ID
    write(h_str, "(I0.2)") IH
    write(fvs_tile, "(I0.2)") IDIM

    ! IMSscf.20191215.C48.nc
    out_file = output_path// &
               "IMSscf."//TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//".C"//fvs_tile//".nc"	
    Call Write_IMS_SND_fSCA_vector(myrank, NPROCS, out_file, lensfc, LENSFC_proc, &
            N_sA, N_sA_Ext, RLA, RLO, OROG, SNCOV_IMS, SNO_IMS_at_Grid) 

    Deallocate(SNCOV_IMS)
    Deallocate(SNO_IMS_at_Grid)
    DEAllocate(RLA)
    DEAllocate(RLO)
    DEAllocate(OROG)

	PRINT*,'Finished'

    CALL MPI_FINALIZE(IERR)

    STOP

END Program ims_scfi


subroutine get_ims_depth_scf(MYRANK, IDIM, JDIM, IY, IM, ID, IH, &
                IVEGSRC, num_subgrd_ims_cels, &
                IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, & 
                vector_restart_prefix, static_prefix, &
                fv3_index, print_debg_info, &
                mp_start, mp_end, LENSFC_proc, &
                RLA, RLO, OROG, SNCOV_IMS, SNO_IMS_at_Grid)

    !----------------------------------------------------------------------
    IMPLICIT NONE
    !
    USE M_SCF_Den

    include 'mpif.h'

    integer, parameter :: dp = kind(1.d0)

    INTEGER, intent(in)    :: MYRANK, IDIM, JDIM, &  
                              IY, IM, ID, IH, LENSFC, IVEGSRC

    CHARACTER(LEN=*), Intent(In)   :: IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH
    CHARACTER(LEN=*), Intent(In)   :: vector_restart_prefix, static_prefix
         
    INTEGER, intent(in) :: num_subgrd_ims_cels
    LOGICAL             :: fv3_index, print_debg_info
    ! for mpi par
    INTEGER   :: mp_start, mp_end
    Integer   :: LENSFC_proc

    Real, intent(out)   :: SNO_IMS_at_Grid(LENSFC_proc), SNCOV_IMS(LENSFC_proc)
    REAL, intent(out)   :: RLA(LENSFC_proc), RLO(LENSFC_proc), OROG(LENSFC_proc)  !, OROG_UF(LENSFC)

    CHARACTER(LEN=5)    :: TILE_NUM
        ! Character(LEN=3)    :: rank_str
    INTEGER	         :: IERR	
    REAL             :: SNDFCS(LENSFC_proc), SWEFCS(LENSFC_proc), STC(LENSFC_proc)
    
    
    REAL            :: VETFCS(LENSFC_proc)  !, SNUP_Array(LENSFC_proc)   
    INTEGER         :: LANDMASK(LENSFC_proc)  !, DAMASK(LENSFC)        
    Integer         :: tile_xy(LENSFC_proc), Idim_xy(LENSFC_proc), Jdim_xy(LENSFC_proc)

    CHARACTER(len=250)   :: ims_inp_file, ims_inp_file_indices
    CHARACTER(len=4)     :: y_str, m_str, d_Str, h_str, fvs_tile
    REAL                 :: lat_min, lat_max, lon_min, lon_max      
    ! Real                 :: SNCOV_IMS(LENSFC_proc)  ! ims resampled at each grid 

    CHARACTER(len=250)       :: forc_inp_file, da_out_file, noda_inp_path
    CHARACTER(LEN=4)         :: RANKCH 
    CHARACTER(LEN=500)       :: static_filename

    Real               :: snodens, SNODENS_Grid(LENSFC_proc)
    Integer            :: veg_type_landice  ! 10.21.20: no assmn over land ice

    INTEGER            :: send_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, &
                          dest_Aoffset_end, arLen, pindex
    INTEGER            :: mpiReal_size, rsize, isize, mpiInt_size, ixy, ie
    REAL               :: tmp

    Real, parameter    :: nodata_val = -9999.9

    integer, parameter :: lsm_type = 2
    
    ! noah models specific? Needed to ID glaciers.
    if (IVEGSRC == 2) then   ! sib
            veg_type_landice=13
    else
            veg_type_landice=15
    endif
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

    write(fvs_tile, "(I0.2)") IDIM
    ! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE 
    static_filename = trim(static_prefix)//"/ufs-land_C"//trim(fvs_tile)//"_static_fields.nc"
    Call ReadTileInfo(static_filename, LENSFC_proc, veg_type_landice, mp_start, mp_end, &
            tile_xy, Idim_xy, Jdim_xy, RLA, RLO, OROG, VETFCS, LANDMASK) 
         
    write(y_str, "(I4)") IY
    write(m_str, "(I0.2)") IM
    write(d_str, "(I0.2)") ID
    write(h_str, "(I0.2)") IH

    WRITE(RANKCH, '(I0.1)') myrank 

    forc_inp_file=TRIM(vector_restart_prefix)//"/ufs_land_restart."// &
        TRIM(y_str)//"-"//TRIM(m_str)//"-"//TRIM(d_str)//"_"//TRIM(h_str)//"-00-00.nc"
    Call read_NoahMP_restart(myrank, forc_inp_file, mp_start, LENSFC_proc, &
                             SNDFCS, SWEFCS, STC)

    ! grid snow density
    call calc_density(LENSFC_proc, lsm_type, LANDMASK, SWEFCS, SNDFCS, &
                        STC, SNODENS_Grid)

    ims_inp_file = TRIM(IMS_SNOWCOVER_PATH)//"IMS.SNCOV."// &
                    TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//TRIM(h_str)//".nc"                !                   
    ims_inp_file_indices = TRIM(IMS_INDEXES_PATH)//"C"//TRIM(ADJUSTL(fvs_tile))//&
                                ".IMS.Indices."//trim(ADJUSTL(RANKCH))//".nc"                                                       
    if (MYRANK==PRINTRANK) PRINT *, 'reading IMS file', trim(ims_inp_file) 
    ! Call Observation_Read_IMS_Full(ims_inp_file, ims_inp_file_indices, &
    !                         MYRANK, JDIM, IDIM, num_subgrd_ims_cels, SNCOV_IMS)
    Call Read_IMS_and_Resample_toTarget(fv3_index, ims_inp_file, IMS_INDEXES_PATH, &
            LENSFC_proc, num_subgrd_ims_cels, IDIM, JDIM, &
            tile_xy, Idim_xy, Jdim_xy, SNCOV_IMS)
    if(myrank==PRINTRANK .and. print_debg_info) then
        PRINT*, "SNCOV from rank: ", MYRANK
        PRINT*, SNCOV_IMS
    endif 
    if (myrank==PRINTRANK) PRINT*,'Finished reading SNCOV, converting to snow depth' 

    call calcSD_noahmp(LENSFC_proc, SNCOV_IMS, VETFCS, &
                            SNODENS_Grid, SNO_IMS_at_Grid)

    if (myrank==PRINTRANK .and. print_debg_info ) then
        PRINT*, "SNCOV derived snodepth at each grid cell from rank: ", MYRANK
        PRINT*, SNO_IMS_at_Grid
    endif
    if (myrank==PRINTRANK) PRINT*,'Finished converting SNCOV observations'

999 CONTINUE
    !PRINT*,'Finished OI DA ON RANK: ', MYRANK
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

    !STOP

    RETURN

END subroutine get_ims_depth_scf

    subroutine ReadTileInfo(filename, LENSFC, veg_type_landice, mp_start, mp_end, &
                            tile_xy, Idim_xy, Jdim_xy, RLA, RLO, OROG, VETFCS, LANDMASK) 

        use netcdf

        ! type(noahmp_type) :: noahmp
        CHARACTER(LEN=*)  :: filename
        integer           :: LENSFC, vector_length, veg_type_landice, mp_start, mp_end
        integer           :: ncid, dimid, varid, status, IDIM
        real              :: RLA(LENSFC), RLO(LENSFC), OROG(LENSFC), VETFCS(LENSFC)
        Integer           :: tile_xy(LENSFC), Idim_xy(LENSFC), Jdim_xy(LENSFC), &
                             LANDMASK(LENSFC)
        LOGICAL             :: file_exists
        ! CHARACTER(LEN=500)  :: filename
        ! Character(LEN=4)    :: grid_str
        Integer             :: Dummy(LENSFC), i
        
        vector_length = mp_end - mp_start + 1 
    
        INQUIRE(FILE=trim(filename), EXIST=file_exists)
        if (.not. file_exists) then 
            print *, 'error,file does not exist', &   
                    trim(filename) , ' exiting'
            call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif           

        status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
        CALL NETCDF_ERR(status, 'opening file '//trim(filename))

        status = nf90_inq_varid(ncid, "cube_tile", varid)
        ! if(status /= nf90_noerr) call handle_err(status)
        CALL NETCDF_ERR(status, 'reading cube_tile var id' )
        status = nf90_get_var(ncid, varid, tile_xy, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading cube tile value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "cube_i", varid)
        CALL NETCDF_ERR(status, 'reading cube i vari id' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, Idim_xy, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading cube i value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "cube_j", varid)
        CALL NETCDF_ERR(status, 'reading cube j vari id' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, Jdim_xy, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading cube j value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "vegetation_category", varid)
        CALL NETCDF_ERR(status, 'reading vegetation_category vari id' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, Dummy, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading vegetation_category value' )
        ! if(status /= nf90_noerr) call handle_err(status)
        VETFCS = float(Dummy)

        status = nf90_inq_varid(ncid, "latitude", varid)
        CALL NETCDF_ERR(status, 'reading latitude var ID' )
        ! if(status /= nf90_noerr) call NETCDF_ERR(status)
        status = nf90_get_var(ncid, varid, RLA, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading latitude value' )
        ! if(status /= nf90_noerr) call NETCDF_ERR(status)

        status = nf90_inq_varid(ncid, "longitude", varid)
        CALL NETCDF_ERR(status, 'reading longitude var ID' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, RLO, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading longitude value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "elevation", varid)
        CALL NETCDF_ERR(status, 'reading elevation var ID' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, OROG, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading elevation value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "land_mask", varid)
        CALL NETCDF_ERR(status, 'reading land_mask var ID' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, LANDMASK, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading land_mask value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        do i = 1, LENSFC 
            ! if land, but not land ice, set mask to 1.
            if (NINT(VETFCS(i)) /=  veg_type_landice) then 
                 LANDMASK(i) = 1 
            else 
                 LANDMASK(i) = 0
            endif 
         enddo

        status = nf90_close(ncid)
        CALL NETCDF_ERR(status, 'closing file' )

    end subroutine ReadTileInfo

    subroutine read_NoahMP_restart(myrank, forc_inp_file, startIndx, vector_length, &
                                 SNDFCS, SWEFCS, STC)

        use netcdf

        CHARACTER(LEN=*)  :: forc_inp_file
        integer           :: myrank, startIndx, vector_length  !, LENSFC
        Real              :: SNDFCS(vector_length), SWEFCS(vector_length), &
                             STC(vector_length)

        integer           :: ncid, dimid, varid, status
        LOGICAL           :: file_exists
        !vector_length =  endIndx - startIndx + 1

        INQUIRE(FILE=trim(forc_inp_file), EXIST=file_exists)
        if (.not. file_exists) then 
            print *, 'error,file does not exist', &   
                    trim(forc_inp_file) , ' exiting'
            call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif        
    
        status = nf90_open(trim(forc_inp_file), NF90_NOWRITE, ncid)

        status = nf90_inq_varid(ncid, "weasd", varid)
        status = nf90_get_var(ncid, varid, SWEFCS,      &
            start = (/startIndx,1/), count = (/vector_length, 1/))

        status = nf90_inq_varid(ncid, "snwdph", varid)
        status = nf90_get_var(ncid, varid, SNDFCS,     &
            start = (/startIndx,1/), count = (/vector_length, 1/))

        status = nf90_inq_varid(ncid, "stc", varid)
        status = nf90_get_var(ncid, varid, STC, &
            start = (/startIndx, 1, 1/), &
            count = (/vector_length, 1, 1/))

        status = nf90_close(ncid)
     

    end subroutine read_NoahMP_restart

    !This reads the whole IMS file and uses a-priori prepared indices to sample those wihin the grid cel
     SUBROUTINE Read_IMS_and_Resample_toTarget(fv3_index, inp_file, IMS_INDEXES_PATH, &
                    LENSFC, num_sub, IDIM, JDIM, tile_xy, Idim_xy, Jdim_xy, SNCOV_IMS)
                    ! Ele               &
                    
        IMPLICIT NONE
    
        include 'mpif.h'                  
    
        !ToDO: Can you use variable length char array ?
        logical                        :: fv3_index
        CHARACTER(LEN=*), Intent(In)   :: inp_file, IMS_INDEXES_PATH !, dim_name
        INTEGER, Intent(In)            :: LENSFC, num_sub, IDIM, JDIM
        Integer, Intent(In)      ::  tile_xy(LENSFC), Idim_xy(LENSFC), Jdim_xy(LENSFC)
        ! ToDO: ims snow cover is of 'byte' type (Chcek the right one)  
        Real, Intent(Out)       :: SNCOV_IMS(LENSFC)     
    
        CHARACTER(LEN=250)      :: inp_file_indices !, dim_name
        CHARACTER(LEN=3)        :: fvs_tile
        INTEGER, ALLOCATABLE    :: SNCOV_IMS_2D_full(:,:)    !SNCOV_IMS_1D(:), 
        Integer                 :: data_grid_ims_ind(num_sub, LENSFC) 
        ! Integer                 :: data_grid_ims_ind_in(num_sub, JDIM, IDIM) 
        ! Real                    :: grid_dat(LENSFC)        
        INTEGER                :: ERROR, NCID, ID_DIM, ID_VAR, DIM_LEN_lat, DIM_LEN_lon
        LOGICAL                :: file_exists
        integer                :: myindx, ixy  
        integer                :: DUMMY1(6, num_sub, JDIM, IDIM)
        CHARACTER(LEN=1)       :: RANKCH
        
        write(fvs_tile, "(I3)") IDIM
        If(fv3_index) then
!index files from fv3 type files
            Do myindx = 1, 6 
                WRITE(RANKCH, '(I1.1)') myindx
                inp_file_indices = TRIM(IMS_INDEXES_PATH)//"/C"//TRIM(ADJUSTL(fvs_tile)) &
                                            //".IMS.Indices.tile"//RANKCH//".nc"
                INQUIRE(FILE=trim(inp_file_indices), EXIST=file_exists)
                if (.not. file_exists) then 
                    print *, 'error,file does not exist', &   
                            trim(inp_file_indices) , ' exiting'
                    call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
                endif           

                ERROR=NF90_OPEN(TRIM(inp_file_indices), NF90_NOWRITE,NCID)
                CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_file_indices) )

                ERROR=NF90_INQ_VARID(NCID, "IMS_Indices", ID_VAR)
                CALL NETCDF_ERR(ERROR, 'READING IMS_Indices ID' )
                ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy1(myindx,:, :,:),  &
                        start = (/ 1, 1, 1 /), count = (/ num_sub, JDIM, IDIM/))
                CALL NETCDF_ERR(ERROR, 'READING snwdph' )
                
                ERROR = NF90_CLOSE(NCID)    
                CALL NETCDF_ERR(ERROR, 'closing file' )
            End do 
            ! 7.20.21 map fv3 to land
            Do ixy=1, LENSFC 
                data_grid_ims_ind(:,ixy)=dummy1(tile_xy(ixy),:,Jdim_xy(ixy),Idim_xy(ixy))
            end do
        else
! read index file for grid/vector array           
            inp_file_indices = TRIM(IMS_INDEXES_PATH)//"C"//TRIM(ADJUSTL(fvs_tile))//&
                                                    ".IMS.Indices"//".nc" 
            INQUIRE(FILE=trim(inp_file_indices), EXIST=file_exists)

            if (.not. file_exists) then
                    print *, 'Observation_Read_IMS_Full error, index file does not exist', &
                            trim(inp_file_indices) , ' exiting'
                    call MPI_ABORT(MPI_COMM_WORLD, 10) 
            endif
        
            ERROR=NF90_OPEN(TRIM(inp_file_indices),NF90_NOWRITE, NCID)
            CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_file_indices) )
        
            ERROR=NF90_INQ_VARID(NCID, 'IMS_Indices', ID_VAR)
            CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV Indices ID' )
            ERROR=NF90_GET_VAR(NCID, ID_VAR, data_grid_ims_ind, start = (/ 1, 1 /), &
                                    count = (/ num_sub, LENSFC/))
            CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV Indices' )
        
            ERROR = NF90_CLOSE(NCID)
        Endif
        ! print*, "IMS Indices 3D"
        ! print*, data_grid_ims_ind
        ! print*

! IMS raw data
        INQUIRE(FILE=trim(inp_file), EXIST=file_exists)

        if (.not. file_exists) then
                print *, 'iObservation_Read_IMS_Full error,file does not exist', &
                        trim(inp_file) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif
    
        ERROR=NF90_OPEN(TRIM(inp_file),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_file) )
    
        ERROR=NF90_INQ_DIMID(NCID, 'lat', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension lat' )
        
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=DIM_LEN_lat)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension Lat' )
        
        ERROR=NF90_INQ_DIMID(NCID, 'lon', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension lon' )
    
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=DIM_LEN_lon)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension Lon' )
    
        ALLOCATE(SNCOV_IMS_2D_full(DIM_LEN_lon, DIM_LEN_lat))   
        ! print*, "initial IMS array size (lon, lat)= ", DIM_LEN_lon, " ",DIM_LEN_lat
        !ALLOCATE(Ele(DIM_LEN_lat, DIM_LEN_lon))
        ERROR=NF90_INQ_VARID(NCID, 'Band1', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, SNCOV_IMS_2D_full, start = (/ 1, 1 /), &
                                count = (/ DIM_LEN_lon, DIM_LEN_lat/))
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV RECORD' )
        ! ERROR=NF90_INQ_VARID(NCID, 'lat', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Lat var ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_IMS_1D)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )
        ! ERROR=NF90_INQ_VARID(NCID, 'lon', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Lon var ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_IMS_1D)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )
    
        ! need to read corresponding elevation values 
        ! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
        
        ERROR = NF90_CLOSE(NCID)
        ! print*, "IMS 2D"
        ! print*, SNCOV_IMS_2D_full
        ! print*
        
        ! DIM_LEN = DIM_LEN_lat * DIM_LEN_lon
        ! ALLOCATE(SNCOV_IMS_1D(DIM_LEN))
        ! SNCOV_IMS_1D = Reshape(SNCOV_IMS_2D_full, (/DIM_LEN/))
    
        ! print*, "IMS 1D"
        ! print*, SNCOV_IMS_1D
        ! print*
        Where(SNCOV_IMS_2D_full /= 4) SNCOV_IMS_2D_full = 0
        Where(SNCOV_IMS_2D_full == 4) SNCOV_IMS_2D_full = 1

        call resample_interpolate_toTarget(SNCOV_IMS_2D_full, data_grid_ims_ind, &
                                           DIM_LEN_lat, DIM_LEN_lon, LENSFC, num_sub, &  !myrank, &
                                           SNCOV_IMS)
    
        ! SNCOV_IMS = Reshape(grid_dat, (/n_lat * n_lon/))
    
        DEALLOCATE(SNCOV_IMS_2D_full)
                  
        RETURN
        
     End SUBROUTINE Read_IMS_and_Resample_toTarget

    subroutine resample_interpolate_toTarget(data_grid_ims, data_grid_ims_ind, &
                                            nlat_ims, nlon_ims, LENSFC, num_sub, &  !myrank, &
                                                grid_dat)
                                                
        Use, Intrinsic :: IEEE_ARITHMETIC
    
        Implicit None
    
        Integer, Intent(In)     :: nlat_ims, nlon_ims, LENSFC, num_sub 
        Integer, Intent(In)     :: data_grid_ims(nlon_ims, nlat_ims), data_grid_ims_ind(num_sub, LENSFC) 
        Real, Intent(Out)       :: grid_dat(LENSFC)
    
        Integer   :: jc, num_loc_counter, ixy !jy, ix, 
        Integer   :: lonlatcoord_ims, loncoord_ims, latcoord_ims
        
        grid_dat = IEEE_VALUE(grid_dat, IEEE_QUIET_NAN)
    
        ! Do jy=1, n_lat
        ! !print*, "process: ", myrank, "loop ", indx
    
        !     Do ix=1, n_lon
        Do ixy = 1, LENSFC

            num_loc_counter = data_grid_ims_ind(1, ixy)
            if (num_loc_counter < 1) then 
                !print*, "no matching values!"
                cycle
            end if
            !print*, "jy ", jy, " ix ", ix
            grid_dat(ixy) = 0.
            Do jc = 2, num_loc_counter+1
                lonlatcoord_ims = data_grid_ims_ind(jc, ixy) - 1 
                latcoord_ims = lonlatcoord_ims / nlon_ims + 1
                loncoord_ims = mod(lonlatcoord_ims, nlon_ims) + 1
                if(latcoord_ims > nlat_ims) then
                    latcoord_ims = nlat_ims
                    print*, "Warning! lat coordinate outside domain boundary"
                endif
                if(loncoord_ims > nlon_ims) then
                    loncoord_ims = nlon_ims
                    print*, "Warning! lon coordinate outside domain boundary"
                endif
                grid_dat(ixy) =  grid_dat(ixy) + data_grid_ims(loncoord_ims, latcoord_ims)              
            End do

            grid_dat(ixy) =  grid_dat(ixy) / num_loc_counter ! first location, num obs

            ! End do
    
        End do
    
        return !grid_dat
    
    End subroutine resample_interpolate_toTarget

 !This reads the whole IMS file and uses a-priori prepared indices to sample those wihin the grid cel
 SUBROUTINE Observation_Read_IMS_Full(inp_file, inp_file_indices, &
				n_lat, n_lon, num_sub, &
				SNCOV_IMS)
				! Ele		&
				
	IMPLICIT NONE

	include 'mpif.h'		  

	!ToDO: Can you use variable length char array ?
	CHARACTER(LEN=*), Intent(In)   :: inp_file, inp_file_indices !, dim_name
	INTEGER, Intent(In)            :: n_lat, n_lon, num_sub
	! ToDO: ims snow cover is of 'byte' type (Chcek the right one)	
	Real, Intent(Out)       :: SNCOV_IMS(n_lat * n_lon) 	

	INTEGER, ALLOCATABLE    :: SNCOV_IMS_2D_full(:,:)    !SNCOV_IMS_1D(:), 
	Integer                 :: data_grid_ims_ind(num_sub, n_lon, n_lat) 
	Real                    :: grid_dat(n_lon, n_lat)
	
	INTEGER                :: ERROR, NCID, ID_DIM, ID_VAR, DIM_LEN, DIM_LEN_lat, DIM_LEN_lon

	ERROR=NF90_OPEN(TRIM(inp_file),NF90_NOWRITE,NCID)
	CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_file) )

	ERROR=NF90_INQ_DIMID(NCID, 'lat', ID_DIM)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension lat' )
	
	ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=DIM_LEN_lat)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension Lat' )
	
	ERROR=NF90_INQ_DIMID(NCID, 'lon', ID_DIM)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension lon' )

	ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=DIM_LEN_lon)
	CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension Lon' )

    ALLOCATE(SNCOV_IMS_2D_full(DIM_LEN_lon, DIM_LEN_lat))	
	! print*, "initial IMS array size (lon, lat)= ", DIM_LEN_lon, " ",DIM_LEN_lat
	!ALLOCATE(Ele(DIM_LEN_lat, DIM_LEN_lon))
	ERROR=NF90_INQ_VARID(NCID, 'Band1', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, SNCOV_IMS_2D_full, start = (/ 1, 1 /), &
							count = (/ DIM_LEN_lon, DIM_LEN_lat/))
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV RECORD' )
	! ERROR=NF90_INQ_VARID(NCID, 'lat', ID_VAR)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Lat var ID' )
	! ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_IMS_1D)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )
	! ERROR=NF90_INQ_VARID(NCID, 'lon', ID_VAR)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Lon var ID' )
	! ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_IMS_1D)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )

	! need to read corresponding elevation values 
	! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
	! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele)
	! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
	
	ERROR = NF90_CLOSE(NCID)

	ERROR=NF90_OPEN(TRIM(inp_file_indices),NF90_NOWRITE, NCID)
	CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_file_indices) )

	ERROR=NF90_INQ_VARID(NCID, 'IMS_Indices', ID_VAR)
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV Indices ID' )
	ERROR=NF90_GET_VAR(NCID, ID_VAR, data_grid_ims_ind, start = (/ 1, 1, 1 /), &
							count = (/ num_sub, n_lon, n_lat/))
	CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV Indices' )

	ERROR = NF90_CLOSE(NCID)

	! print*, "IMS 2D"
	! print*, SNCOV_IMS_2D_full
	! print*
	
	! DIM_LEN = DIM_LEN_lat * DIM_LEN_lon
	! ALLOCATE(SNCOV_IMS_1D(DIM_LEN))
	! SNCOV_IMS_1D = Reshape(SNCOV_IMS_2D_full, (/DIM_LEN/))

	! print*, "IMS 1D"
	! print*, SNCOV_IMS_1D
	! print*
    Where(SNCOV_IMS_2D_full /= 4) SNCOV_IMS_2D_full = 0
	Where(SNCOV_IMS_2D_full == 4) SNCOV_IMS_2D_full = 1
	

	! print*, "IMS binary"
	! print*, SNCOV_IMS_1D
	! print*

	! print*, "IMS Indices 3D"
	! print*, data_grid_ims_ind
	! print*

	call resample_to_model_tiles_intrp(SNCOV_IMS_2D_full, data_grid_ims_ind, &
	                                   DIM_LEN_lat, DIM_LEN_lon, n_lat, n_lon, num_sub, &  !myrank, &
		                               grid_dat)

	SNCOV_IMS = Reshape(grid_dat, (/n_lat * n_lon/))

	DEALLOCATE(SNCOV_IMS_2D_full)
			  
	RETURN
	
 End SUBROUTINE Observation_Read_IMS_Full

 subroutine resample_to_model_tiles_intrp(data_grid_ims, data_grid_ims_ind, &
                                            nlat_ims, nlon_ims, n_lat, n_lon, num_sub, &  !myrank, &
											grid_dat)
											
    Use, Intrinsic :: IEEE_ARITHMETIC

    Implicit None

    Integer, Intent(In)     :: nlat_ims, nlon_ims, n_lat, n_lon, num_sub 
    Integer, Intent(In)     :: data_grid_ims(nlon_ims, nlat_ims), data_grid_ims_ind(num_sub, n_lon, n_lat) 
    Real, Intent(Out)       :: grid_dat(n_lon, n_lat)

	Integer   :: jc, jy, ix, num_loc_counter
	Integer   :: lonlatcoord_ims, loncoord_ims, latcoord_ims
    
	grid_dat = IEEE_VALUE(grid_dat, IEEE_QUIET_NAN)

    Do jy=1, n_lat
    !print*, "process: ", myrank, "loop ", indx

        Do ix=1, n_lon
            
            num_loc_counter = data_grid_ims_ind(1, ix, jy)
            if (num_loc_counter < 1) then 
                !print*, "no matching values!"
                cycle
            end if
            !print*, "jy ", jy, " ix ", ix
            grid_dat(ix, jy) = 0.
			Do jc = 2, num_loc_counter+1
				lonlatcoord_ims = data_grid_ims_ind(jc, ix, jy) - 1 
				latcoord_ims = lonlatcoord_ims / nlon_ims + 1
				loncoord_ims = mod(lonlatcoord_ims, nlon_ims) + 1
				if(latcoord_ims > nlat_ims) then
					latcoord_ims = nlat_ims
					print*, "Warning! lat coordinate outside domain boundary"
				endif
				if(loncoord_ims > nlon_ims) then
					loncoord_ims = nlon_ims
					print*, "Warning! lon coordinate outside domain boundary"
				endif
                grid_dat(ix, jy) =  grid_dat(ix, jy) + data_grid_ims(loncoord_ims, latcoord_ims)              
            End do

            grid_dat(ix, jy) =  grid_dat(ix, jy) / num_loc_counter ! first location, num obs

        End do

    End do

    return !grid_dat

End subroutine resample_to_model_tiles_intrp

  Subroutine Write_IMS_SND_fSCA_vector(myrank, NPROCS, output_file, &
            lensfc, LENSFC_land, N_sA, N_sA_Ext,    &    
            Lat_atObs, Lon_atObs, OROG_at_stn, SNCOV_IMS, SND_IMS) 

        implicit none
        include "mpif.h"
        
        CHARACTER(LEN=*), Intent(In)    :: output_file
        integer, intent(in)             :: myrank, NPROCS, lensfc, LENSFC_land
        integer, intent(in)             :: N_sA, N_sA_Ext
    
        Real, intent(in)    :: Lat_atObs(LENSFC_land), Lon_atObs(LENSFC_land), &
                               OROG_at_stn(LENSFC_land)

        Real, intent(in)    :: SNCOV_IMS(LENSFC_land), SND_IMS(LENSFC_land)  

        Real    :: RLA(LENSFC), RLO(LENSFC), OROG(LENSFC), SNCOV(LENSFC), SND(LENSFC)  

        integer            :: fsize=65536, inital=0
        integer            :: header_buffer_val = 16384
        integer            :: error, i, ncid
        integer            :: dim_x
        integer            :: id_sn
        integer            :: id_lon, id_lat, id_oro, id_imscov, id_imssnd
        integer            :: mpiReal_size, rsize

        Integer     :: dest_Aoffset, dest_Aoffset_end, arLen, IERR, ixy, jndx

        rsize = SIZEOF(Lat_atObs)/LENSFC_land

        Call MPI_TYPE_SIZE(MPI_REAL, mpiReal_size, IERR) 
        If (rsize == 4 ) then 
            mpiReal_size = MPI_REAL4
        elseif (rsize == 8 ) then 
            mpiReal_size = MPI_REAL8
        elseif (rsize == 16 ) then 
            mpiReal_size = MPI_REAL16
        else
            PRINT*," Possible mismatch between Fortran Real ", rsize," and Mpi Real ", mpiReal_size
            Stop
        endif

        if (MYRANK > 0) then
            call MPI_SEND(Lat_atObs, LENSFC_land, mpiReal_size, 0,   &
                                    MYRANK+1, MPI_COMM_WORLD, IERR) 
            call MPI_SEND(Lon_atObs, LENSFC_land, mpiReal_size, 0,   &
                                    100*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(OROG_at_stn, LENSFC_land, mpiReal_size, 0,   &
                                    200*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(SNCOV_IMS, LENSFC_land, mpiReal_size, 0,   &
                                    300*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(SND_IMS, LENSFC_land, mpiReal_size, 0,   &
                                   400*(MYRANK+1), MPI_COMM_WORLD, IERR)
        Endif
        if (MYRANK == 0) then
            RLA(1:1+N_sA) = Lat_atObs
            RLO(1:1+N_sA) = Lon_atObs
            OROG(1:1+N_sA) = OROG_at_stn
            SNCOV(1:1+N_sA) = SNCOV_IMS
            SND(1:1+N_sA) = SND_IMS
            Do ixy =  1, NPROCS - 1   ! sender proc index within tile group
                if(ixy < N_sA_Ext) then   !p_tRank== 0) then 
                    dest_Aoffset = ixy * (N_sA + 1) + 1    ! 1
                    dest_Aoffset_end = dest_Aoffset + N_sA    !+ 1 
                    arLen = N_sA + 1  
                else
                    dest_Aoffset = ixy * N_sA + N_sA_Ext + 1  !N_sA_Ext*(N_sA+1) + (MYRANK-N_sA_Ext)*N_sA     
                    dest_Aoffset_end = dest_Aoffset + N_sA - 1 
                    arLen = N_sA
                endif
                call MPI_RECV(RLA(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                ixy+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(RLO(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                100*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(OROG(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                200*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(SNCOV(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                300*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(SND(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                400*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            End do
        ! print*,"Process ", myrank, "writing output data to: ",trim(output_file)
        !--- create the file
        ! if (myrank ==) then 
            error = NF90_CREATE(output_file, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), ncid, initialsize=inital, chunksize=fsize)
            call netcdf_err(error, 'CREATING FILE='//trim(output_file) )     
            ! print*,"Process ", myrank, "created file: ",trim(output_file)
            !--- define dimensions
            error = nf90_def_dim(ncid, 'numobs', lensfc, dim_x)
            call netcdf_err(error, 'DEFINING location DIMENSION' )

            !--- define fields
            ! error = nf90_def_var(ncid, 'location', NF90_FLOAT, dim_x, id_x)
            ! call netcdf_err(error, 'DEFINING location FIELD' )
            ! error = nf90_put_att(ncid, id_x, "long_name", "location")
            ! call netcdf_err(error, 'DEFINING location LONG NAME' )
            ! error = nf90_put_att(ncid, id_x, "units", "none")
            ! call netcdf_err(error, 'DEFINING location UNITS' )
            ! error = nf90_put_att(ncid, id_x, "cartesian_axis", "X")
            ! call netcdf_err(error, 'WRITING location FIELD' )
            
            error = nf90_def_var(ncid, 'lon', NF90_DOUBLE, dim_x, id_lon)
            call netcdf_err(error, 'DEFINING lon' )
            error = nf90_put_att(ncid, id_lon, "long_name", "longitude")
            call netcdf_err(error, 'DEFINING lon LONG NAME' )

            error = nf90_def_var(ncid, 'lat', NF90_DOUBLE, dim_x, id_lat)
            call netcdf_err(error, 'DEFINING lat' )
            error = nf90_put_att(ncid, id_lat, "long_name", "latitude")
            call netcdf_err(error, 'DEFINING lat LONG NAME' )

            error = nf90_def_var(ncid, 'oro', NF90_DOUBLE, dim_x, id_oro)
            call netcdf_err(error, 'DEFINING oro' )
            error = nf90_put_att(ncid, id_oro, "long_name", "orography")
            call netcdf_err(error, 'DEFINING oro LONG NAME' )

            error = nf90_def_var(ncid, 'IMSscf', NF90_DOUBLE, dim_x, id_imscov)
            call netcdf_err(error, 'DEFINING IMSscf' )
            error = nf90_put_att(ncid, id_imscov, "long_name", "IMS snow covered fraction")
            call netcdf_err(error, 'DEFINING IMSscf LONG NAME' )
            error = nf90_put_att(ncid, id_imscov, "units", "-")
            call netcdf_err(error, 'DEFINING IMSscf UNITS' )

            error = nf90_def_var(ncid, 'IMSsnd', NF90_DOUBLE, dim_x, id_imssnd)
            call netcdf_err(error, 'DEFINING IMSsnd' )
            error = nf90_put_att(ncid, id_imssnd, "long_name", "IMS snow depth")
            call netcdf_err(error, 'DEFINING IMSsnd LONG NAME' )
            error = nf90_put_att(ncid, id_imssnd, "units", "mm")
            call netcdf_err(error, 'DEFINING IMSsnd UNITS' )

            error = nf90_enddef(ncid, header_buffer_val,4,0,4)
            call netcdf_err(error, 'DEFINING HEADER' )

        
            error = nf90_put_var( ncid, id_lat, RLA)
            call netcdf_err(error, 'WRITING id_lat RECORD' )

            error = nf90_put_var( ncid, id_lon, RLO)
            call netcdf_err(error, 'WRITING id_lon RECORD' ) 

            error = nf90_put_var( ncid, id_oro, OROG)
            call netcdf_err(error, 'WRITING id_oro RECORD' )

            error = nf90_put_var( ncid, id_imscov, SNCOV)
            call netcdf_err(error, 'WRITING id_imscov RECORD' )

            error = nf90_put_var( ncid, id_imssnd, SND)
            call netcdf_err(error, 'WRITING id_imssnd RECORD' )
     
            error = nf90_close(ncid)
            call netcdf_err(error, 'closing DA file' )

        endif    
             
 End subroutine Write_IMS_SND_fSCA_vector

!  subroutine read_forc_noah()
    
!     forc_inp_file = TRIM(SFC_FORECAST_PREFIX)//TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"."//TRIM(h_str)// &
!                             "0000.sfc_data."//TILE_NUM//".nc"

!     if (MYRANK==PRINTRANK) PRINT *, 'snowDA: reading model backgroundfile', trim(forc_inp_file)                               
!     Call READ_Forecast_Data_atPath(forc_inp_file, veg_type_landice, LENSFC, SWEFCS, SNDFCS, &
!                                 VETFCS, LANDMASK)

!     SNODENS_Grid = SWEFCS/SNDFCS
!     WHERE (SNODENS_Grid < 0.0001) SNODENS_Grid = 0.

!     if (COUNT (LANDMASK==1 .and. SNDFCS> 0.01) > 0) then 
!         ! mean density over snow-covered land
!         snodens = SUM(SNODENS_Grid, Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
!                     / COUNT (LANDMASK==1 .and. SNDFCS> 0.01) 
!         if (MYRANK==PRINTRANK) & 
!                 PRINT *, 'snowDA: mean density ', snodens 
!     else 
!         snodens = 0.1  ! default value if have no snow in tile
!         if (MYRANK==PRINTRANK) & 
!                 PRINT *, 'snowDA: no snow in current tiles, using default density ', snodens 
!     endif
!     ! PRINT *, 'snowDA:density ', MYRANK, snodens 
!     ! for grid cells with no valid density, fill in the average snodens
!     Where(.not.(LANDMASK==1 .and. SNDFCS>0.01 )) SNODENS_Grid = snodens
!     If (p_tRank==0)  print*, "Tile ", p_tN, ' mean snow density', snodens
!     tmp = SUM(SWEFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
!                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.01)
!     If (p_tRank==0)  print*, "Tile ", p_tN,  ' mean SWE', tmp
!     tmp = SUM(SNDFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
!                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.01)
!     If (p_tRank==0)  print*, "Tile ", p_tN,  ' mean SND', tmp  
    
    ! if (land_type == 1) then 
    !     ! SNUP array will be used later, to test whether SCF > 100%
    !     call CalcSWEFromSnowCover(SNCOV_IMS, VETFCS, LENSFC_proc, SNO_IMS_at_Grid, SNUP_Array)
    !         SNO_IMS_at_Grid = SNO_IMS_at_Grid/SNODENS_Grid ! convert SWE to SND
    ! else

!     end subroutine read_forc_noah

 