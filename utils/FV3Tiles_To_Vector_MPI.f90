 Subroutine FV3Tiles_To_Vector(LENS_OUT, vector_rand_ens)
    ! PROGRAM FV3Tiles_To_Vector(LENS_OUT, vector_rand_ens)

    IMPLICIT NONE
    !
    include 'mpif.h'
    
    Integer, intent(in)   :: LENS_OUT
    Real, intent(out)     :: vector_rand_ens(LENS_OUT)
    INTEGER :: IDIM, JDIM, NUM_TILES, IY, IM, ID, IH
    REAL    :: FH, DELTSFC
    INTEGER :: IERR
    INTEGER :: NPROCS, MYRANK, NUM_THREADS, NUM_PARTHDS, MAX_TASKS
    REAL                :: horz_len_scale, ver_len_scale, temp_len_scale 
    Integer             :: ens_size, t_indx
    ! LOGICAL             :: rcov_localize, ens_inflate
    CHARACTER(LEN=500)  :: static_filename, fv3_prefix, vector_prefix
    CHARACTER(LEN=50)   :: rand_var
    integer             :: PRINTRANK
    LOGICAL             :: print_debg_info


! noah mp updte   
    Integer                :: LENSFC, LENSFC_land 
    INTEGER, allocatable   :: tile_xy(:), Idim_xy(:), Jdim_xy(:), VETFCS(:)
    Real, allocatable      :: RLA_land(:), RLO_land(:), OROG_land(:)

    ! Real, allocatable      :: rand_Ta(:), rand_Prec(:), rand_sRad(:), rand_lRad(:)
    ! Real, allocatable      :: rand_sHum(:), rand_wSpeed 
    Real, allocatable      :: rand_Ta3D(:,:,:)
    ! Real, allocatable      :: vector_rand_Ta(:)

    CHARACTER(LEN=500)     :: fv3_inp_file, vector_inp_file
    CHARACTER(LEN=2)       :: RANKCH 
    
    Integer                :: ixy, vector_size, arr_indx
    ! real,allocatable       :: DUMMY1(:,:), DUMMY1f(:,:), DUMMY2(:,:), DUMMY3(:,:)   !, DUMMY4(:)
    ! for mpi par
        ! INTEGER            :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end
    INTEGER            :: N_sA, N_sA_Ext, mp_start, mp_end
    INTEGER            :: mpiReal_size, mpiInt_size, rsize, isize, dest_Aoffset, dest_Aoffset_end, arLen
    integer            :: group_world, comm_tile, tile_group !comm_world, 
    integer, allocatable :: tile_members(:)  

    NAMELIST/NAMSNO/ IDIM, JDIM, NUM_TILES, IY, IM, ID, IH, FH, DELTSFC, &
                    horz_len_scale, ver_len_scale, temp_len_scale, ens_size, t_indx, &
                    static_filename, fv3_prefix, vector_prefix,  rand_var, &
                    PRINTRANK, print_debg_info             
    !
    DATA IDIM,JDIM,NUM_TILES/96,96,6/ 
    DATA IY,IM,ID,IH,FH/1997,8,2,0,0./
    DATA DELTSFC/0.0/, MAX_TASKS/99999/
    DATA horz_len_scale/55.0/
    DATA ver_len_scale/800./
    DATA temp_len_scale/24./
    DATA ens_size/20/
    DATA t_indx/2/
    DATA print_debg_info/.false./
    DATA PRINTRANK/4/
    DATA static_filename/""/
    DATA fv3_prefix/"./"/
    DATA vector_prefix/""/
    DATA rand_var/"smc"/
    ! DATA obs_srch_rad/250.0/   
    ! DATA stdev_obsv_depth/40.0/
    ! DATA stdev_obsv_sncov/80.0/
    ! DATA stdev_back/30.0/
    ! DATA num_assim_steps/1/  ! For multiple time steps of assimilation
    ! DATA dT_Asssim/24.0/     ! hrs. For multiple time steps of assimilation
    
    CALL MPI_INIT(IERR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, IERR)

    NUM_THREADS = NUM_PARTHDS()

    PRINT*
    PRINT*,"STARTING CYCLE PROGRAM ON RANK ", MYRANK
    PRINT*,"RUNNING WITH ", NPROCS, "TASKS"
    PRINT*,"AND WITH ", NUM_THREADS, " THREADS."

    IF (MYRANK==0) PRINT*,"READING NAMELIST."

    CALL BAOPENR(360, "tiles_to_vector.nml", IERR)             !"snowDA.nml"   !
    READ(360, NML=NAMSNO)
    IF (MYRANK==0) WRITE(6, NAMSNO)
    LENSFC = IDIM*JDIM ! TOTAL NUMBER OF POINTS FOR THE CUBED-SPHERE TILE

    ALLOCATE(SNOFCS(LENSFC))
    ALLOCATE(LANDMASK(LENSFC))
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    ! IF (MAX_TASKS < 99999 .AND. MYRANK > (MAX_TASKS - 1)) THEN
    !     PRINT*,"USER SPECIFIED MAX NUMBER OF TASKS: ", MAX_TASKS
    !     PRINT*,"WILL NOT RUN CYCLE PROGRAM ON RANK: ", MYRANK
    !     GOTO 333
    ! ENDIF

    !   comm_world = comm_world
    Allocate(tile_members(NUM_TILES))
    do ixy = 0, NUM_TILES - 1   
        tile_members(ixy+1) = ixy
    enddo
    call MPI_Comm_group(MPI_COMM_WORLD, group_world, IERR)        
    call MPI_Group_incl(group_world, NUM_TILES, tile_members, tile_group, IERR)
    call MPI_Comm_create(MPI_COMM_WORLD, tile_group, comm_tile, IERR)
    print*, "Proc ",myrank, " done creating new tile group and comm"

    IF ( MYRANK > (NUM_TILES - 1)) THEN
        PRINT*," WILL NOT RUN ON RANK: ", MYRANK
        GOTO 444
    ENDIF
      
    !   static_filename = "ufs-land_C"// trim(str(NDIM)) //"_static_fields.nc"
    call ReadTileInfo(trim(static_filename), vector_size, tile_xy, Idim_xy, Jdim_xy, &
                            RLA_land, RLO_land, OROG_land, VETFCS)              
    ! IF (MYRANK==0) then
        PRINT*," proc ", myrank, "finished reading tile info"
        !  print*, "tile_xy", tile_xy
        !  print*, "Idim_xy", Idim_xy
        !  print*, "Jdim_xy", Jdim_xy
        ! print*, "vegetation_category", VETFCS_land   
        ! print*, "RLA", RLA_land
        ! print*, "RLO", RLO_land
        ! print*, "OROG", OROG_land
    ! Endif

    ! allocate(vector_rand_Ta(vector_size))

    ! allocate(rand_Ta(LENSFC), rand_Prec(LENSFC), rand_sRad(LENSFC))
    ! allocate(rand_lRad(LENSFC),rand_sHum(LENSFC), rand_wSpeed(LENSFC))
    allocate(rand_Ta3D(6, IDIM, JDIM))
    
    ! Read fv3    
    !workg_T162_984x488.tile03.nc 
    Do ixy = 1, NUM_TILES       !vector_length 
        WRITE(RANKCH, '(I2.2)') ixy     !myrank
        fv3_inp_file = TRIM(fv3_prefix)//"tile"//RANKCH//".nc"
        IF (MYRANK==0) PRINT*," reading file ", fv3_inp_file
        Call Read_FV3_File(fv3_inp_file, rand_var, IDIM, JDIM, t_indx, rand_Ta3D(ixy,:,:))
        IF (MYRANK==0) PRINT*,"finished reading fv3"
    end do 
    
    Do ixy=1, vector_length 
        vector_rand_ens = rand_Ta3D(tile_xy(ixy), Jdim_xy(ixy), Idim_xy(ixy))
        ! vector_rand_Ta(ixy) = rand_Ta3D(tile_xy(ixy), Jdim_xy(ixy), Idim_xy(ixy))
    end do
    vector_inp_file

! Real data type size corresponding to mpi
    ! rsize = SIZEOF(DELTSFC)
    ! Call MPI_TYPE_SIZE(MPI_REAL, mpiReal_size, IERR) 
    ! If (rsize == 4 ) then 
    !         mpiReal_size = MPI_REAL4
    ! elseif (rsize == 8 ) then 
    !         mpiReal_size = MPI_REAL8
    ! elseif (rsize == 16 ) then 
    !         mpiReal_size = MPI_REAL16
    ! else
    !         PRINT*," Possible mismatch between Fortran Real ", rsize," and Mpi Real ", mpiReal_size
    !         Stop
    ! endif
    ! isize = SIZEOF(ens_size) 
    ! Call MPI_TYPE_SIZE(MPI_INTEGER, mpiInt_size, IERR) 
    ! If (isize == 2 ) then 
    !     mpiInt_size = MPI_INTEGER2
    ! elseif (isize == 4 ) then 
    !     mpiInt_size = MPI_INTEGER4
    ! elseif (isize == 8 ) then 
    !     mpiInt_size = MPI_INTEGER8
    ! else
    !     PRINT*," Possible mismatch between Fortran Int ", isize," and Mpi Int ", mpiInt_size
    !     Stop
    ! endif
    ! if (MYRANK > 0) then
    !     call MPI_SEND(SNDANL, LENSFC, mpiReal_size, 0,   &
    !                     MYRANK+1, MPI_COMM_WORLD, IERR) 
    ! endif
    ! if (MYRANK == 0) then 
    !     DUMMY1(1,:) = SNDANL     
    !     Do ixy =  1, NUM_TILES-1   ! sender proc index within tile group
    !         call MPI_RECV(DUMMY1(ixy+1,:), LENSFC, mpiReal_size, ixy,      &
    !                         ixy+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
    !     End do
    ! endif
    ! IF (MYRANK==0) then
    !     PRINT*," root received analyis outputs; broadcasting..."
    ! Endif
    !  ! 7.20.21 map fv3 to land
    ! Do ixy=1, vector_length 
    !     vector_rand_Ta(ixy) = rand_Ta2D(Jdim_xy(ixy), Idim_xy(ixy))
    ! end do

    ! 7.20.21 map fv3 to land
    ! Do ixy=0, LENSFC_land - 1 
    !     arr_indx = JDIM * (Idim_xy(mp_start + ixy) - 1) + Jdim_xy(mp_start + ixy)
    !     if(arr_indx > LENSFC) then
    !         print*, "erroneous index ", arr_indx
    !         stop
    !     endif
    !     anl_out(ixy+1) = dummy1(tile_xy(mp_start + ixy), arr_indx)
    ! end do
 
    Deallocate(tile_members)
    deallocate(tile_xy)
    deallocate(Idim_xy)
    deallocate(Jdim_xy)
    deallocate(RLA_land)
    deallocate(RLO_land)
    deallocate(OROG_land)
    deallocate(VETFCS)
    ! deallocate(rand_Ta, rand_Prec, rand_sRad)
    ! deallocate(rand_lRad, rand_sHum, rand_wSpeed)    
    deallocate(rand_Ta3D)
    ! deallocate(vector_rand_Ta)

444 Continue
! 333 CONTINUE

    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    PRINT*
    PRINT*,'COMPLETED NORMALLY ON RANK: ', MYRANK

    CALL MPI_FINALIZE(IERR)

    STOP

 END subroutine FV3Tiles_To_Vector

 subroutine ReadTileInfo(filename, vector_length, tile_xy, Idim_xy, Jdim_xy, &
                        RLA, RLO, OROG, VETFCS) ! vegetation_category)

        use netcdf

        ! type(noahmp_type) :: noahmp
        CHARACTER(LEN=*)  :: filename
        integer           :: vector_length
        integer           :: ncid, dimid, varid, status
        real, allocatable :: RLA(:), RLO(:), OROG(:)
        Integer, allocatable      :: tile_xy(:), Idim_xy(:), Jdim_xy(:), VETFCS(:) !, vegetation_category(:)
        LOGICAL                   :: file_exists

        INQUIRE(FILE=trim(filename), EXIST=file_exists)
        if (.not. file_exists) then 
            print *, 'error,file does not exist', &   
                    trim(filename) , ' exiting'
            call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif           

        status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
        CALL NETCDF_ERR(status, 'opening file '//trim(filename))

        status = nf90_inq_dimid(ncid, "location", dimid)
        status = nf90_inquire_dimension(ncid, dimid, len = vector_length)

        allocate(tile_xy(vector_length))
        allocate(Idim_xy(vector_length))
        allocate(Jdim_xy(vector_length))
        allocate(VETFCS(vector_length))
        allocate(RLA(vector_length))
        allocate(RLO(vector_length))
        allocate(OROG(vector_length))

        status = nf90_inq_varid(ncid, "cube_tile", varid)
        ! if(status /= nf90_noerr) call handle_err(status)
        CALL NETCDF_ERR(status, 'reading cube_tile var id' )
        status = nf90_get_var(ncid, varid, tile_xy, &
            start = (/1/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading cube tile value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "cube_i", varid)
        CALL NETCDF_ERR(status, 'reading cube i vari id' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, Idim_xy, &
            start = (/1/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading cube i value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "cube_j", varid)
        CALL NETCDF_ERR(status, 'reading cube j vari id' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, Jdim_xy, &
            start = (/1/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading cube j value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "vegetation_category", varid)
        CALL NETCDF_ERR(status, 'reading vegetation_category vari id' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, VETFCS, &
            start = (/1/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading vegetation_category value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "latitude", varid)
        ! if(status /= nf90_noerr) call NETCDF_ERR(status)
        status = nf90_get_var(ncid, varid, RLA, &
            start = (/1/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading latitude value' )
        ! if(status /= nf90_noerr) call NETCDF_ERR(status)

        status = nf90_inq_varid(ncid, "longitude", varid)
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, RLO, &
            start = (/1/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading longitude value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "elevation", varid)
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, OROG, &
            start = (/1/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading elevation value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        ! status = nf90_inq_varid(ncid, "land_mask", varid)
        ! if(status /= nf90_noerr) call handle_err(status)
        ! status = nf90_get_var(ncid, varid, this%land_mask, &
        !     start = (/namelist%begsub/), count = (/namelist%lensub/))
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_close(ncid)
        CALL NETCDF_ERR(status, 'closing file' )

    end subroutine ReadTileInfo

    subroutine Read_FV3_File(filename, var_name, IDIM, JDIM, t_indx, Var_Out)

        use netcdf

        ! type(noahmp_type) :: noahmp
        CHARACTER(LEN=*)  :: filename, var_name
        integer           :: IDIM, JDIM, t_indx
        Real              :: Var_Out(IDIM, JDIM)   !, DUMMY(IDIM, JDIM)
        ! real, allocatable :: RLA(:), RLO(:), OROG(:)
        LOGICAL                   :: file_exists
        integer           :: ncid, dimid, varid, status

        INQUIRE(FILE=trim(filename), EXIST=file_exists)
        if (.not. file_exists) then 
            print *, 'error,file does not exist', &   
                    trim(filename) , ' exiting'
            call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif           

        status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
        CALL NETCDF_ERR(status, 'opening file '//trim(filename))

        status = nf90_inq_varid(ncid, trim(var_name), varid)
        ! if(status /= nf90_noerr) call handle_err(status)
        CALL NETCDF_ERR(status, 'reading '//trim(var_name)//' var id' )
        status = nf90_get_var(ncid, varid, Var_Out, &
            start = (/1, 1, t_indx/), count = (/IDIM, JDIM, 1/))
        CALL NETCDF_ERR(status, 'reading cube tile value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        ! Var_Out = RESHAPE(DUMMY, (/IDIM * JDIM/))   

        status = nf90_close(ncid)
        CALL NETCDF_ERR(status, 'closing file '//trim(filename) )

    end subroutine Read_FV3_File

    SUBROUTINE NETCDF_ERR( ERR, STRING )

    !--------------------------------------------------------------
    ! IF AT NETCDF CALL RETURNS AN ERROR, PRINT OUT A MESSAGE
    ! AND STOP PROCESSING.
    !--------------------------------------------------------------

        IMPLICIT NONE

        include 'mpif.h'

        INTEGER, INTENT(IN) :: ERR
        CHARACTER(LEN=*), INTENT(IN) :: STRING
        CHARACTER(LEN=80) :: ERRMSG

        IF( ERR == NF90_NOERR )RETURN
        ERRMSG = NF90_STRERROR(ERR)
        PRINT*,''
        PRINT*,'FATAL ERROR: ', TRIM(STRING), ': ', TRIM(ERRMSG)
        PRINT*,'STOP.'
        CALL MPI_ABORT(MPI_COMM_WORLD, 999)

        RETURN
    END SUBROUTINE NETCDF_ERR