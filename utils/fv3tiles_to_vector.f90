 Subroutine fv3tiles_to_vector(LENS_OUT, num_ens, vector_rand_ens)
    ! PROGRAM FV3Tiles_To_Vector(LENS_OUT, vector_rand_ens)

    IMPLICIT NONE
    !
    include 'mpif.h'
    
    Integer, intent(in)   :: LENS_OUT, num_ens
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
  
    Integer                :: LENSFC, LENSFC_land 
    INTEGER, allocatable   :: tile_xy(:), Idim_xy(:), Jdim_xy(:), VETFCS(:)
    Real, allocatable      :: RLA_land(:), RLO_land(:), OROG_land(:)
    ! Real, allocatable      :: rand_Ta(:), rand_Prec(:), rand_sRad(:), rand_lRad(:)
    ! Real, allocatable      :: rand_sHum(:), rand_wSpeed 
    Real, allocatable      :: rand_Ta3D(:,:,:), sfc_wts_out(:,:,:)
    ! Real, allocatable      :: vector_rand_Ta(:)
    CHARACTER(LEN=500)     :: fv3_inp_file, vector_inp_file
    CHARACTER(LEN=2)       :: RANKCH 

    CHARACTER(len=250)     :: forc_inp_file, forc_inp_path
    CHARACTER(len=500)     :: forc_inp_file_ens
    CHARACTER(LEN=2)       :: ensCH 
    
    Integer                :: ixy, vector_size, arr_indx, iproc, ie
    ! real,allocatable       :: DUMMY1(:,:)
    integer, allocatable :: tile_members(:), n_surf_vars  
    INTEGER            :: mpiReal_size, rsize, isize, mpiInt_size

    LOGICAL     :: file_exists
    Integer     :: ncid, error, varid
    character(len=128)  :: forc_var_list(6)
    Real                :: forcArray(LENS_OUT)

    NAMELIST/NAMSNO/ IDIM, JDIM, NUM_TILES, IY, IM, ID, IH, FH, DELTSFC, &
                    horz_len_scale, ver_len_scale, temp_len_scale, ens_size, t_indx, &
                    static_filename, fv3_prefix, vector_prefix,  rand_var, &
                    PRINTRANK, print_debg_info, n_surf_vars             
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
    Data n_surf_vars/6/
    ! DATA obs_srch_rad/250.0/   
    ! DATA stdev_obsv_depth/40.0/
    ! DATA stdev_obsv_sncov/80.0/
    ! DATA stdev_back/30.0/
    ! DATA num_assim_steps/1/  ! For multiple time steps of assimilation
    ! DATA dT_Asssim/24.0/     ! hrs. For multiple time steps of assimilation
    
    CALL MPI_INIT(IERR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, IERR)

    PRINT*," fv3tiles_to_vector RANK ", MYRANK, " WITH ", NPROCS, "TASKS"
 
    rsize = SIZEOF(horz_len_scale)
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
    isize = SIZEOF(vector_size) 
    Call MPI_TYPE_SIZE(MPI_INTEGER, mpiInt_size, IERR) 
    If (isize == 2 ) then 
        mpiInt_size = MPI_INTEGER2
    elseif (isize == 4 ) then 
        mpiInt_size = MPI_INTEGER4
    elseif (isize == 8 ) then 
        mpiInt_size = MPI_INTEGER8
    else
        PRINT*," Possible mismatch between Fortran Int ", isize," and Mpi Int ", mpiInt_size
        Stop
    endif

    PRINT*,"READING NAMELIST."
    ! CALL BAOPENR(360, "tiles_to_vector.nml", IERR)             !"snowDA.nml"   !
    open(360, file="tiles_to_vector.nml", form="formatted")
    read(360, NAMSNO)
    ! READ(360, NML=NAMSNO)
    close(360)

    WRITE(6, NAMSNO)
    LENSFC = IDIM*JDIM ! TOTAL NUMBER OF POINTS FOR THE CUBED-SPHERE TILE

    vector_size = LENS_OUT
    allocate(tile_xy(vector_size))
    allocate(Idim_xy(vector_size))
    allocate(Jdim_xy(vector_size))
    allocate(VETFCS(vector_size))
    allocate(RLA_land(vector_size))
    allocate(RLO_land(vector_size))
    allocate(OROG_land(vector_size))
    !   static_filename = "ufs-land_C"// trim(str(NDIM)) //"_static_fields.nc"
    call ReadTileInfo(trim(static_filename), vector_size, tile_xy, Idim_xy, Jdim_xy, &
                            RLA_land, RLO_land, OROG_land, VETFCS)              
    ! IF (MYRANK==0) then
        PRINT*, "finished reading tile info" !" proc ", myrank, 
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
    allocate(sfc_wts_out(n_surf_vars, IDIM, JDIM))

    Do ie = 1, num_ens
        IF (MYRANK==0) PRINT*," calling stand alone stochy, ens ", ie
        call standalone_stochy_sfc(IDIM, JDIM, n_surf_vars, sfc_wts_out)
        IF (MYRANK==0) PRINT*," Done stand alone stochy, ens ", ie     

        Do ixy = 1, n_surf_vars        !vector_length 
            if (MYRANK /= 0) then  
                call MPI_SEND(sfc_wts_out(ixy,:,:),IDIM*JDIM, mpiReal_size, 0, &
                                10*(MYRANK+1), MPI_COMM_WORLD, IERR) 
            else
                rand_Ta3D(1,:,:) = sfc_wts_out(ixy,:,:)
                Do iproc = 1, NUM_TILES-1 
                    call MPI_RECV(rand_Ta3D(iproc+1,:,:), IDIM*JDIM, mpiReal_size, &
                    iproc, 10*(iproc+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                Enddo

                Do iv=1, vector_size 
                    vector_rand_ens(iv) = rand_Ta3D(tile_xy(iv), Jdim_xy(iv), Idim_xy(iv))
                end do

                WRITE(ensCH, '(I2.2)') ie
                forc_inp_file_ens=TRIM(forc_inp_path)//"/ens"//ensCH//"/"//TRIM(forc_inp_file)

                ! REstart file
                INQUIRE(FILE=trim(forc_inp_file_ens), EXIST=file_exists)
                if (.not. file_exists) then 
                    print *, 'error,file does not exist', &   
                            trim(forc_inp_file_ens) , ' exiting'
                    call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
                endif
                error = nf90_open(trim(forc_inp_file_ens), NF90_WRITE, ncid)
                call netcdf_err(error, 'opening restart file' )

                ! Start writing restart file
                error = nf90_inq_varid(ncid, trim(forc_var_list(ixy)), varid)
                call netcdf_err(error, 'getting varid '//trim(forc_var_list(ixy)) )
                !read
                ERROR=NF90_GET_VAR(NCID, varid, forcArray)
                CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )
                
                forcArray = forcArray * vector_rand_ens

                error = nf90_put_var(ncid, varid , forcArray   , &
                    start = (/1,1/), count = (/vector_size, 1/))
                call netcdf_err(error, 'writing '//trim(forc_var_list(ixy)) )

                error = nf90_close(ncid)
                call netcdf_err(error, 'closing '//trim(forc_var_list(ixy)) )
                
            Endif
            
        end do 
        ! each of the 6 proc receives one variable
        ! Do ixy = 0, NUM_TILES-1       !vector_length 
        !     if (MYRANK /= ixy) then  
        !         call MPI_SEND(sfc_wts_out(ixy+1,:,:),IDIM*JDIM, mpiReal_size, ixy, &
        !                         10*(MYRANK+1), MPI_COMM_WORLD, IERR) 
        !     else
        !         rand_Ta3D(ixy+1,:,:) = sfc_wts_out(ixy+1,:,:)
        !         Do iproc = 0, NUM_TILES-1 
        !             if (iproc /= ixy) then 
        !                 call MPI_RECV(rand_Ta3D(iproc+1,:,:), IDIM*JDIM, mpiReal_size, &
        !                 iproc, 10*(iproc+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
        !             endif
        !         Enddo
        !     Endif
        ! end do 
        
        Do ixy=1, vector_size 
            vector_rand_ens(ie, ixy) = rand_Ta3D(tile_xy(ixy), Jdim_xy(ixy), Idim_xy(ixy))
            ! vector_rand_Ta(ixy) = rand_Ta3D(tile_xy(ixy), Jdim_xy(ixy), Idim_xy(ixy))
        end do
        PRINT*,"proc ", myrank, " finished copying rand for ens ", ie
    Enddo

    ! Do ixy = 1, NUM_TILES       !vector_length 
        !     WRITE(RANKCH, '(I2.2)') ixy     !myrank
        !     fv3_inp_file = TRIM(fv3_prefix)//RANKCH//".nc"
        !     ! IF (MYRANK==0) 
        !     PRINT*," reading file ", fv3_inp_file
        !     Call Read_FV3_File(trim(fv3_inp_file), trim(rand_var), IDIM, JDIM, t_indx, &
        !                         rand_Ta3D(ixy,:,:))
        !     ! IF (MYRANK==0) 
        !     PRINT*,"finished reading ", fv3_inp_file
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
    
    
    Call ReadRestartNoahMP(myrank, TRIM(forc_inp_file), noda_inp_path, mp_start, mp_end, &
                               LENSFC, noahmp, SNDnoDA, SWEnoDA)  !, SNDFCS, SWEFCS)
    
    

    ! Deallocate(tile_members)
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
    deallocate(sfc_wts_out)

    PRINT*,'COMPLETED FV3Tiles_To_Vector NORMALLY'

 END subroutine FV3Tiles_To_Vector

!  subroutine write_rand_pattern(file_name)
   
!     ! call get_outfile(fname)
!     ! write(strid,'(I2.2)') my_id+1
!     ! if (ntile_out.EQ.0) write_this_tile=.true.
!     ! if ((my_id+1).EQ.ntile_out) write_this_tile=.true.
!     ! print*,trim(fname)//'.tile'//strid//'.nc',write_this_tile
!     CHARACTER(LEN=500)  :: rand_pattern_filename
!     if (write_this_tile) then
        
!         rand_pattern_filename = trim(fname)//'.tile'//strid//'.nc'
 
!         fid=30+my_id
!         ierr=nf90_create(trim(fname)//'.tile'//strid//'.nc',cmode=NF90_CLOBBER,ncid=ncid)
!         ierr=NF90_DEF_DIM(ncid,"grid_xt",nx,xt_dim_id)
!         ierr=NF90_DEF_DIM(ncid,"grid_yt",ny,yt_dim_id)
!         if (do_skeb)ierr=NF90_DEF_DIM(ncid,"p_ref",nlevs,zt_dim_id)
!         ierr=NF90_DEF_DIM(ncid,"time",NF90_UNLIMITED,time_dim_id)
!        !> - Define the dimension variables.
!         ierr=NF90_DEF_VAR(ncid,"grid_xt",NF90_FLOAT,(/ xt_dim_id /), xt_var_id)
!         ierr=NF90_PUT_ATT(ncid,xt_var_id,"long_name","T-cell longitude")
!         ierr=NF90_PUT_ATT(ncid,xt_var_id,"cartesian_axis","X")
!         ierr=NF90_PUT_ATT(ncid,xt_var_id,"units","degrees_E")
!         ierr=NF90_DEF_VAR(ncid,"grid_yt",NF90_FLOAT,(/ yt_dim_id /), yt_var_id)
!         ierr=NF90_PUT_ATT(ncid,yt_var_id,"long_name","T-cell latitude")
!         ierr=NF90_PUT_ATT(ncid,yt_var_id,"cartesian_axis","Y")
!         ierr=NF90_PUT_ATT(ncid,yt_var_id,"units","degrees_N")
!         ierr=NF90_DEF_VAR(ncid,"grid_lat",NF90_FLOAT,(/ xt_dim_id, yt_dim_id, time_dim_id /), var_id_lat)
!         ierr=NF90_PUT_ATT(ncid,var_id_lat,"long_name","T-cell latitudes")
!         ierr=NF90_PUT_ATT(ncid,var_id_lat,"units","degrees_N")
!         ierr=NF90_PUT_ATT(ncid,var_id_lat,"missing_value",undef)
!         ierr=NF90_PUT_ATT(ncid,var_id_lat,"_FillValue",undef)
!         ierr=NF90_DEF_VAR(ncid,"grid_lon",NF90_FLOAT,(/ xt_dim_id, yt_dim_id, time_dim_id /), var_id_lon)
!         ierr=NF90_PUT_ATT(ncid,var_id_lon,"long_name","T-cell longitudes")
!         ierr=NF90_PUT_ATT(ncid,var_id_lon,"units","degrees_N")
!         ierr=NF90_PUT_ATT(ncid,var_id_lon,"missing_value",undef)
!         ierr=NF90_PUT_ATT(ncid,var_id_lon,"_FillValue",undef)
!         ierr=NF90_DEF_VAR(ncid,"tile_num",NF90_FLOAT,(/ xt_dim_id, yt_dim_id, time_dim_id /), var_id_tile)
!         ierr=NF90_PUT_ATT(ncid,var_id_tile,"long_name","tile number")
!         ierr=NF90_PUT_ATT(ncid,var_id_tile,"missing_value",undef)
!         ierr=NF90_PUT_ATT(ncid,var_id_tile,"_FillValue",undef)
!         if (do_skeb)then
!            ierr=NF90_DEF_VAR(ncid,"p_ref",NF90_FLOAT,(/ zt_dim_id /), zt_var_id)
!            ierr=NF90_PUT_ATT(ncid,zt_var_id,"long_name","reference pressure")
!            ierr=NF90_PUT_ATT(ncid,zt_var_id,"cartesian_axis","Z")
!            ierr=NF90_PUT_ATT(ncid,zt_var_id,"units","Pa")
!         endif
!         ierr=NF90_DEF_VAR(ncid,"time",NF90_FLOAT,(/ time_dim_id /), time_var_id)
!         ierr=NF90_DEF_VAR(ncid,"time",NF90_FLOAT,(/ time_dim_id /), time_var_id)
!         ierr=NF90_PUT_ATT(ncid,time_var_id,"long_name","time")
!         ierr=NF90_PUT_ATT(ncid,time_var_id,"units","hours since 2014-08-01 00:00:00")
!         ierr=NF90_PUT_ATT(ncid,time_var_id,"cartesian_axis","T")
!         ierr=NF90_PUT_ATT(ncid,time_var_id,"calendar_type","JULIAN")
!         ierr=NF90_PUT_ATT(ncid,time_var_id,"calendar","JULIAN")
!         if (do_sppt)then
!            ierr=NF90_DEF_VAR(ncid,"sppt_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), varid1)
!            ierr=NF90_PUT_ATT(ncid,varid1,"long_name","sppt pattern")
!            ierr=NF90_PUT_ATT(ncid,varid1,"units","None")
!            ierr=NF90_PUT_ATT(ncid,varid1,"missing_value",undef)
!            ierr=NF90_PUT_ATT(ncid,varid1,"_FillValue",undef)
!            ierr=NF90_PUT_ATT(ncid,varid1,"cell_methods","time: point")
!         endif
!         if (do_shum)then
!            ierr=NF90_DEF_VAR(ncid,"shum_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), varid2)
!            ierr=NF90_PUT_ATT(ncid,varid2,"long_name","shum pattern")
!            ierr=NF90_PUT_ATT(ncid,varid2,"units","None")
!            ierr=NF90_PUT_ATT(ncid,varid2,"missing_value",undef)
!            ierr=NF90_PUT_ATT(ncid,varid2,"_FillValue",undef)
!            ierr=NF90_PUT_ATT(ncid,varid2,"cell_methods","time: point")
!         endif
!         if (do_skeb)then
!            ierr=NF90_DEF_VAR(ncid,"skebu_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,zt_dim_id,time_dim_id/), varid3)
!            ierr=NF90_DEF_VAR(ncid,"skebv_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,zt_dim_id,time_dim_id/), varid4)
!            ierr=NF90_PUT_ATT(ncid,varid3,"long_name","skeb u pattern")
!            ierr=NF90_PUT_ATT(ncid,varid3,"units","None")
!            ierr=NF90_PUT_ATT(ncid,varid3,"missing_value",undef)
!            ierr=NF90_PUT_ATT(ncid,varid3,"_FillValue",undef)
!            ierr=NF90_PUT_ATT(ncid,varid3,"cell_methods","time: point")
!            ierr=NF90_PUT_ATT(ncid,varid4,"long_name","skeb v pattern")
!            ierr=NF90_PUT_ATT(ncid,varid4,"units","None")
!            ierr=NF90_PUT_ATT(ncid,varid4,"missing_value",undef)
!            ierr=NF90_PUT_ATT(ncid,varid4,"_FillValue",undef)
!            ierr=NF90_PUT_ATT(ncid,varid4,"cell_methods","time: point")
!         endif
!         if (lndp_type > 0)then
!            do l=1,n_var_lndp
!               ierr=NF90_DEF_VAR(ncid,lndp_var_list(l),NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), varidl(l))
!               ierr=NF90_PUT_ATT(ncid,varidl(l),"long_name",lndp_var_list(l)//" pattern")
!               ierr=NF90_PUT_ATT(ncid,varidl(l),"units","None")
!               ierr=NF90_PUT_ATT(ncid,varidl(l),"missing_value",undef)
!               ierr=NF90_PUT_ATT(ncid,varidl(l),"_FillValue",undef)
!               ierr=NF90_PUT_ATT(ncid,varidl(l),"cell_methods","time: point")
!            enddo
!         endif
!         ierr=NF90_ENDDEF(ncid)
!         ierr=NF90_PUT_VAR(ncid,xt_var_id,grid_xt)
!         ierr=NF90_PUT_VAR(ncid,yt_var_id,grid_yt)

!          ! put lat lon and tile number
!     !ierr=NF90_PUT_VAR(ncid,var_id_lon,transpose(xlon(isc:iec,jsc:iec)),(/1,1,1/))
!     !ierr=NF90_PUT_VAR(ncid,var_id_lat,transpose(xlat(isc:iec,jsc:iec)),(/1,1,1/))
!     ierr=NF90_PUT_VAR(ncid,var_id_lon,transpose(xlon(:,:)),(/1,1,1/))
!     ierr=NF90_PUT_VAR(ncid,var_id_lat,transpose(xlat(:,:)),(/1,1,1/))
!     tile_number=my_id+1
       
!     ierr=NF90_PUT_VAR(ncid,var_id_tile,tile_number,(/1,1,1/))

!            do l=1,n_var_lndp
!               do j=1,ny
!                  workg(:,j)=sfc_wts(j,:,l)
!               enddo
!               ierr=NF90_PUT_VAR(ncid,varidl(l),workg,(/1,1,tpt/))
!            enddo
!         endif
!         ierr=NF90_PUT_VAR(ncid,time_var_id,ts,(/tpt/))
!         endif
!         tpt=tpt+1
!      endif
!   enddo
!   if (write_this_tile) ierr=NF90_CLOSE(ncid)


!  end subroutine write_rand_pattern

 subroutine ReadTileInfo(filename, vector_length, tile_xy, Idim_xy, Jdim_xy, &
                        RLA, RLO, OROG, VETFCS) ! vegetation_category)

        use netcdf

        ! type(noahmp_type) :: noahmp
        CHARACTER(LEN=*)  :: filename
        integer           :: vector_length, vector_length_in
        integer           :: ncid, dimid, varid, status
        real              :: RLA(vector_length), RLO(vector_length), OROG(vector_length)
        Integer  :: tile_xy(vector_length), Idim_xy(vector_length), &
                    Jdim_xy(vector_length), VETFCS(vector_length) !, vegetation_category(:)
        LOGICAL                   :: file_exists

        INQUIRE(FILE=trim(filename), EXIST=file_exists)
        if (.not. file_exists) then 
            print *, 'error,file does not exist', &   
                    trim(filename) , ' exiting'
            ! call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
            STOP
        endif           

        status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
        CALL NETCDF_ERR(status, 'opening file '//trim(filename))

        status = nf90_inq_dimid(ncid, "location", dimid)
        status = nf90_inquire_dimension(ncid, dimid, len = vector_length_in)
        
        if (vector_length /= vector_length_in) then
            print*, "wrong vector size, stop"
            stop
        endif
        ! allocate(tile_xy(vector_length))
        ! allocate(Idim_xy(vector_length))
        ! allocate(Jdim_xy(vector_length))
        ! allocate(VETFCS(vector_length))
        ! allocate(RLA(vector_length))
        ! allocate(RLO(vector_length))
        ! allocate(OROG(vector_length))

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
            ! call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
            STOP
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
        
        use netcdf

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
        ! CALL MPI_ABORT(MPI_COMM_WORLD, 999)
        STOP

        RETURN
    END SUBROUTINE NETCDF_ERR