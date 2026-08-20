      module hdf5_readers
      use hdf5
      implicit none
      contains 

      subroutine read1d_32(group_id, name, x, n, ierr)

        use hdf5
        implicit none

        integer(HID_T), intent(in) :: group_id
        character(len=*), intent(in) :: name
        integer, intent(out) :: ierr

        real(4), allocatable, intent(out) :: x(:)
        integer, intent(out) :: n

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T) :: dims(1)

        call h5dopen_f(group_id, name, dset_id, ierr)
        if (ierr /= 0) return

        call h5dget_space_f(dset_id, space_id, ierr)
        call h5sget_simple_extent_dims_f(space_id, dims, dims, ierr)

        n = dims(1)
        allocate(x(n))

        call h5dread_f(dset_id, H5T_NATIVE_REAL, x, dims, ierr)

        call h5sclose_f(space_id, ierr)
        call h5dclose_f(dset_id, ierr)

      end subroutine read1d_32

      subroutine read1d_64(group_id, name, x, n, ierr)

        use hdf5
        implicit none

        integer(HID_T), intent(in) :: group_id
        character(len=*), intent(in) :: name
        integer, intent(out) :: ierr
        integer, intent(out) :: n

        real(8), allocatable, intent(out) :: x(:)

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T) :: dims(1)

        call h5dopen_f(group_id, name, dset_id, ierr)
        if (ierr /= 0) return

        call h5dget_space_f(dset_id, space_id, ierr)
        call h5sget_simple_extent_dims_f(space_id, dims, dims, ierr)

        n = dims(1)
        allocate(x(n))

        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, x, dims, ierr)

        call h5sclose_f(space_id, ierr)
        call h5dclose_f(dset_id, ierr)

      end subroutine read1d_64

      subroutine readint1d_32(group_id, name, x, n, ierr)

        use hdf5
        implicit none

        integer(HID_T), intent(in) :: group_id
        character(len=*), intent(in) :: name
        integer, intent(out) :: ierr

        integer(4), allocatable, intent(out) :: x(:)
        integer, intent(out) :: n

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T) :: dims(1)

        call h5dopen_f(group_id, name, dset_id, ierr)
        if (ierr /= 0) return

        call h5dget_space_f(dset_id, space_id, ierr)
        call h5sget_simple_extent_dims_f(space_id, dims, dims, ierr)

        n = dims(1)
        allocate(x(n))

        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, x, dims, ierr)

        call h5sclose_f(space_id, ierr)
        call h5dclose_f(dset_id, ierr)

      end subroutine readint1d_32

      subroutine readint1d_64(group_id, name, x, n, ierr)

        use hdf5
        implicit none

        integer(HID_T), intent(in) :: group_id
        character(len=*), intent(in) :: name
        integer, intent(out) :: ierr
        integer, intent(out) :: n

        integer(8), allocatable, intent(out) :: x(:)

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T) :: dims(1)

        call h5dopen_f(group_id, name, dset_id, ierr)
        if (ierr /= 0) return

        call h5dget_space_f(dset_id, space_id, ierr)
        call h5sget_simple_extent_dims_f(space_id, dims, dims, ierr)

        n = dims(1)
        allocate(x(n))

        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, x, dims, ierr)

        call h5sclose_f(space_id, ierr)
        call h5dclose_f(dset_id, ierr)

      end subroutine readint1d_64

      end module hdf5_readers

      module readgridlevels
      implicit none 
      contains

*********************************************************************
       subroutine read_grid_levels(filename_params, level, ngrids)
*********************************************************************
*      reads data from ENZO hdf5 inputs
*********************************************************************

       implicit none 

       character(len=*), intent(in) :: filename_params
       integer, allocatable, intent(out) :: level(:)
       integer, intent(out) :: ngrids

       integer unit, ios 
       integer igrid, nx, ny, nz 
       integer capacity 
       integer i, j, k
       real*8 :: left(3), right(3)
       real*8 dx 
       character(len=512) :: filename, line
       real*8, allocatable :: tmp_dx(:)
       integer, allocatable :: tmp_level(:)
       real*8, allocatable :: dxgrid(:)

       filename = trim(filename_params) // '.hierarchy'

       open(unit, file=filename, status='old', action='read',
     &      iostat=ios)

       if (ios.ne.0) then 
          write(*,*) 'Error opening file: ', trim(filename)
          stop
       end if

       ngrids = 0 
       capacity = 1000 

       allocate(dxgrid(capacity))
       
       do ! main loop, undetermined length (number of grids)
        ! Look for next grid
        read(unit, '(A)', iostat=ios) line
        if (ios.ne.0) exit ! end of file

        if (index(line, 'Grid').eq.1) then
          read(line(index(line, '=')+1:),*) igrid

          ! If necessary, increase the size of the array 
          if (igrid.gt.capacity) then 
            allocate(tmp_dx(capacity))
            tmp_dx = 0.d0 
            tmp_dx = dxgrid
            deallocate(dxgrid)
            capacity = capacity * 2
            allocate(dxgrid(capacity))
            dxgrid = 0.d0
            dxgrid(1:int(capacity/2)) = tmp_dx
            deallocate(tmp_dx)
          end if

          ! Read until GridDimension
          do ! loop until we reach GridDimension
            read(unit, '(A)', iostat=ios) line
            if (ios.ne.0) exit ! end of file

            if (index(line, 'GridDimension').eq.1) then
              read(line(index(line, '=')+1:),*) nx, ny, nz
              exit
            end if
          end do ! loop until we reach GridDimension

          ! Read until GridLeftEdge
          do ! loop until we reach GridLeftEdge
            read(unit, '(A)', iostat=ios) line
            if (ios.ne.0) exit ! end of file

            if (index(line, 'GridLeftEdge').eq.1) then
              read(line(index(line, '=')+1:),*) 
     &                         left(1),left(2),left(3)
              exit
            end if
          end do ! loop until we reach GridLeftEdge

          ! Read until GridRightEdge
          do ! loop until we reach GridRightEdge
            read(unit, '(A)', iostat=ios) line
            if (ios.ne.0) exit ! end of file

            if (index(line, 'GridRightEdge').eq.1) then
              read(line(index(line, '=')+1:),*) 
     &                         right(1),right(2),right(3)
              exit
            end if
          end do ! loop until we reach GridRightEdge

          ! Remove 3 ghost zones 
          nx = nx - 6 
          dx = (right(1) - left(1)) / dble(nx)
          dxgrid(igrid) = dx
          ngrids = max(ngrids, igrid)

        end if
       end do ! main loop, undetermined length (number of grids)

       close(unit)

*      Reside to the number of grids actually read
       allocate(tmp_dx(ngrids))
       tmp_dx = dxgrid(1:ngrids)
       deallocate(dxgrid)
       allocate(dxgrid(ngrids))
       dxgrid = tmp_dx
       deallocate(tmp_dx)

*      Now, assign the level of each grid based on its dx
       allocate(level(ngrids))
       do i = 1, ngrids
         level(i) = int(log(dxgrid(1)/dxgrid(i))/log(2.0) + 0.5)
       end do

       end subroutine read_grid_levels

      end module readgridlevels

*********************************************************************
       subroutine read_enzo_hdf5(iter,nx,ny,nz,t,zeta,
     &             u2dm,u3dm,u4dm,masap,rxpa,rypa,rzpa,oripa,
     &             n_dm,var,n_st,uv,um,hubble_littleh,
     &             filespersnap)
*********************************************************************
*      reads data from ENZO hdf5 inputs
*********************************************************************

       use hdf5
       use hdf5_readers
       use iso_c_binding
       use readgridlevels
       implicit none

       include 'input_files/asohf_parameters.dat'

       integer nx,ny,nz,iter,ndxyz,filespersnap
       real*4 t,aaa,bbb,ccc,map,zeta,uv,um,hubble_littleh
       integer var !(=1: only dm; =2: dm+stars)

       integer n_dm,n_st,nbas,nst0

       real*4 u2dm(parti_read)
       real*4 u3dm(parti_read)
       real*4 u4dm(parti_read)
       real*4 masap(parti_read)
       real*4 rxpa(parti_read)
       real*4 rypa(parti_read)
       real*4 rzpa(parti_read)
       integer oripa(parti_read)

       real cio_mass,cio_speed,cio_length,cio_alpha,cio_xc,cio_yc,cio_zc
       common /conv_io/ cio_mass,cio_speed,cio_length,cio_alpha,cio_xc,
     &                 cio_yc,cio_zc
       real fact_mass,fact_speed,fact_length

       integer i,j,n,icpu,nread,i1,i2,nthis
       integer nmembers,objtype,imember,igrid
       character(len=64) :: obj_name
       character(len=64) :: gridname
       logical :: exists, verbose
       integer, allocatable :: level(:)
       integer ngrids

*     io variables
       real*4,allocatable::scr41(:)
       real*4,allocatable::scr42(:,:)
       integer*4,allocatable::scrint41(:)
       integer*8,allocatable::scrint81(:)
       integer*4,allocatable::scrint1(:) !delete!
       real*8,allocatable::scr81(:)
       real*8,allocatable::scr82(:,:)
       
       type(c_ptr) :: redshift_ptr
       real(c_double), target :: redshift

       character*4 iter_string
       character*4 cpu_string

       character*200 postfix,enzofolderprefix,enzoprefix
       character*200 enzoparamspostfix
       common /POSTFIX/ postfix,enzofolderprefix,enzoprefix,
     &                  enzoparamspostfix
       character*200 filename,filename_params,filenamebase

       character(len=512) :: line
       character(len=128)  :: key
       
       ! HDF5 variables
       integer(HID_T) :: file_id, dset_id, dtype_id, space_id, attr_id
       integer(HID_T) :: group_id
       integer :: ierr
       integer(hsize_t), dimension(1) :: dims1d, maxdims1d
       integer(hsize_t), dimension(2) :: dims2d, maxdims2d
       integer(hsize_t) :: precision

***** what i need to output/load is *******
*      rxpa,rypa,rzpa,u2dm,u3dm,u4dm,masap,oripa
*      n_dm, n_st, zeta
*****************************************
      write(*,*) '!!!!'
      write(*,*) 'ENZO reader, warning!'
      write(*,*) 'Not implemented for stars yet'
      write(*,*) '!!!!'

      verbose = .false.

*     READING DATA
      WRITE(ITER_STRING, '(I4.4)') ITER !For saving files to disk
      WRITE(*,*) ' > Reading iter ',ITER_STRING

      filenamebase = './simulation/' // trim(enzofolderprefix) //
     &       trim(ITER_STRING) // '/' // trim(enzoprefix) // 
     &       trim(ITER_STRING) // trim(postfix)

      if (var.eq.2) then
        write(*,*) 'WARNING: reading stars not implemented yet'
        write(*,*) 'set VAR=1 to read only DM particles'
        stop
      end if
      n_st = 0 !!! Since we are ignoring stars for now, we set n_st=0

      ! Obtain the redshift of the snapshot
      filename_params = './simulation/' // trim(enzofolderprefix) //
     &       trim(ITER_STRING) // '/' // trim(enzoprefix) // 
     &       trim(ITER_STRING) // trim(enzoparamspostfix)
      open(10,file=trim(filename_params),status='old',action='read')
      do 
        read(10,'(A)',iostat=ierr) line
        if (ierr /= 0) exit
        if (index(line, 'CosmologyCurrentRedshift') > 0) then
          j = index(line, '=') + 1
          read(line(j:),*) zeta
          exit
        end if
      end do
      close(10)
      write(*,*) 'Redshift: zeta=', zeta

      write(*,*) 'filespersnap=',filespersnap
      ! check all the cpu files exist
      do icpu=0,filespersnap-1
       WRITE(CPU_STRING, '(I4.4)') icpu
       filename = trim(filenamebase) // '.cpu' // trim(CPU_STRING)
       inquire(file=trim(filename), exist=exists, iostat=ierr)
       if (.not.exists) then
        write(*,*) 'File ',trim(filename),' does not exist'
        stop
       end if
      end do

      call read_grid_levels(filename_params, level, ngrids)
      write(*,*) 'Number of grids: ', ngrids
      do i=0,maxval(level)
        write(*,*) 'Number of grids at level ', i, ': ', count(level==i)
      end do

      ! Now, we have to loop over CPUs and over Grids, and if they have
      !  particles (and only in that case), we read them and store them.
      call h5open_f(ierr)
      nread = 0

      do icpu=0,filespersnap-1 ! ======================================
        write(cpu_string, '(I4.4)') icpu
        filename = trim(filenamebase) // '.cpu' // trim(cpu_string) 

        ! Open file to retrieve all the grids in it
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr)

        ! number of members in the group (number of grids plus extra possible data)
        call h5gn_members_f(file_id, "/", nmembers, ierr)

        if (verbose)
     & write(*,*) 'file ',trim(filename),' has ',nmembers,' members'

        do imember=0,nmembers-1 ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call h5gget_obj_info_idx_f(file_id, "/", imember, obj_name, 
     &                                  objtype, ierr)
          if (objtype.ne.H5G_GROUP_F) cycle ! skip if not a group 
          if (index(trim(obj_name), "Grid") == 0) cycle ! skip if not a grid
          gridname = trim(obj_name)

          if (verbose) then
            write(*,*) 'reading grid ',trim(gridname)
          end if

          ! Check if the grid contains particles, else, skip it
          call h5lexists_f(file_id, 
     &     trim(gridname) // "/particle_position_x", exists, ierr)

          if (.not.exists) then
            if (verbose) then
              write(*,*) 'grid ',trim(gridname),
     &                   ' has no particles, skipping'
            end if
            cycle
          end if

          if (verbose) then
            write(*,*) 'grid ',trim(gridname),
     &               ' has particles, reading them'
          end if

          ! get the grid number from the name 
          read(gridname(5:),*) igrid

*         1st, get number of particles 
          call h5gopen_f(file_id, trim(gridname), group_id, ierr)
          call h5dopen_f(group_id, "particle_position_x", dset_id, ierr)
          call h5dget_space_f(dset_id, space_id, ierr)
          call h5sget_simple_extent_dims_f(space_id, dims1d, maxdims1d, 
     &                                     ierr)

          nthis = dims1d(1)
          if (verbose)
     &     write(*,*) 'nread=',nread,' nthis=',nthis

          ! get the datatype and its precision (assume it's the same for all)
          call h5dget_type_f(dset_id, dtype_id, ierr)
          call h5tget_precision_f(dtype_id, precision, ierr)  ! bits
          call h5sclose_f(space_id, ierr)
          call h5dclose_f(dset_id, ierr)
          
          ! Determine the indices 
          i1 = nread + 1
          i2 = nread + nthis
          nread = nread + nthis

          if (verbose) then
            write(*,*) 'Precision of data in file: ', precision, ' bits'
          end if
          if (precision.ne.32 .and. precision.ne.64) then
            write(*,*) 'Unsupported precision: ', precision
            stop
          end if

          if (nread.gt.parti_read) then
            write(*,*) 'warning: bad dimensioning of parti_read',
     &                  nread,'>',parti_read
            stop
          end if
          
          !! TO DO
          if (precision == 32) then 
            call read1d_32(group_id, "particle_position_x", scr41, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            rxpa(i1:i2)=scr41(1:nthis)
          
            call read1d_32(group_id, "particle_position_y", scr41, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            rypa(i1:i2)=scr41(1:nthis)
          
            call read1d_32(group_id, "particle_position_z", scr41, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            rzpa(i1:i2)=scr41(1:nthis)
          
            call read1d_32(group_id, "particle_velocity_x", scr41, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            u2dm(i1:i2)=scr41(1:nthis)
          
            call read1d_32(group_id, "particle_velocity_y", scr41, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            u3dm(i1:i2)=scr41(1:nthis)
          
            call read1d_32(group_id, "particle_velocity_z", scr41, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            u4dm(i1:i2)=scr41(1:nthis)
          
            call read1d_32(group_id, "particle_mass", scr41, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            masap(i1:i2)=scr41(1:nthis) / 8.0**level(igrid) ! correct for the level of the grid
          
            call readint1d_32(group_id, "particle_id", scrint41, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            oripa(i1:i2)=scrint41(1:nthis)
            deallocate(scr41, scrint41)
          
          else if (precision == 64) then 
          
            call read1d_64(group_id, "particle_position_x", scr81, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            rxpa(i1:i2)=scr81(1:nthis)
          
            call read1d_64(group_id, "particle_position_y", scr81, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            rypa(i1:i2)=scr81(1:nthis)
          
            call read1d_64(group_id, "particle_position_z", scr81, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            rzpa(i1:i2)=scr81(1:nthis)
          
            call read1d_64(group_id, "particle_velocity_x", scr81, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            u2dm(i1:i2)=scr81(1:nthis)
          
            call read1d_64(group_id, "particle_velocity_y", scr81, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            u3dm(i1:i2)=scr81(1:nthis)
          
            call read1d_64(group_id, "particle_velocity_z", scr81, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            u4dm(i1:i2)=scr81(1:nthis)
          
            call read1d_64(group_id, "particle_mass", scr81, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            masap(i1:i2)=scr81(1:nthis) / 8.0**level(igrid) ! correct for the level of the grid
          
            call readint1d_32(group_id, "particle_index", scrint41, 
     &                     n, ierr)
            if (n.ne.nthis) then
             write(*,*) 'Error: number of particles read != nthis'
             write(*,*) 'n=',n,' nthis=',nthis
             stop
            end if
            oripa(i1:i2)=scrint41(1:nthis)
            deallocate(scr81, scrint41)
          end if

          call h5gclose_f(group_id, ierr)
        end do ! imember=0,nmembers-1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call h5fclose_f(file_id, ierr)
      end do ! icpu=0,filespersnap-1 ==================================

      n_dm = nread

       write(*,*)
       write(*,*) 'Total number of particles read: ', nread
       write(*,*) 'input. dm x positions (min,max):',
     &            minval(rxpa(1:n_dm)),maxval(rxpa(1:n_dm))
       write(*,*) 'input. dm y positions (min,max):',
     &            minval(rypa(1:n_dm)),maxval(rypa(1:n_dm))
       write(*,*) 'input. dm z positions (min,max):',
     &            minval(rzpa(1:n_dm)),maxval(rzpa(1:n_dm))
       write(*,*) 'input. dm x velocities (min,max):',
     &            minval(u2dm(1:n_dm)),maxval(u2dm(1:n_dm))
       write(*,*) 'input. dm y velocities (min,max):',
     &            minval(u3dm(1:n_dm)),maxval(u3dm(1:n_dm))
       write(*,*) 'input. dm z velocities (min,max):',
     &            minval(u4dm(1:n_dm)),maxval(u4dm(1:n_dm))
       write(*,*) 'input. dm masses (min,max):',
     &            minval(masap(1:n_dm)),maxval(masap(1:n_dm))
       write(*,*) 'input. dm unique ids (min,max):',
     &            minval(oripa(1:n_dm)),maxval(oripa(1:n_dm))
       write(*,*)
       write(*,*) 'total dm particles in iter=',n_dm


       if (var.eq.2) then
        write(*,*)
        write(*,*) 'input. st x positions (min,max):',
     &     minval(rxpa(n_dm+1:n_dm+n_st)),maxval(rxpa(n_dm+1:n_dm+n_st))
        write(*,*) 'input. st y positions (min,max):',
     &     minval(rypa(n_dm+1:n_dm+n_st)),maxval(rypa(n_dm+1:n_dm+n_st))
        write(*,*) 'input. st z positions (min,max):',
     &     minval(rzpa(n_dm+1:n_dm+n_st)),maxval(rzpa(n_dm+1:n_dm+n_st))
        write(*,*) 'input. st x velocities (min,max):',
     &     minval(u2dm(n_dm+1:n_dm+n_st)),maxval(u2dm(n_dm+1:n_dm+n_st))
        write(*,*) 'input. st y velocities (min,max):',
     &     minval(u3dm(n_dm+1:n_dm+n_st)),maxval(u3dm(n_dm+1:n_dm+n_st))
        write(*,*) 'input. st z velocities (min,max):',
     &     minval(u4dm(n_dm+1:n_dm+n_st)),maxval(u4dm(n_dm+1:n_dm+n_st))
        write(*,*) 'input. st masses (min,max):',
     &   minval(masap(n_dm+1:n_dm+n_st)),maxval(masap(n_dm+1:n_dm+n_st))
        write(*,*) 'input. st unique ids (min,max):',
     &   minval(oripa(n_dm+1:n_dm+n_st)),maxval(oripa(n_dm+1:n_dm+n_st))
        write(*,*)
        write(*,*) 'total stellar particles in iter=',n_st
       end if

       fact_mass=cio_mass/um
       fact_speed=(cio_speed/uv)*(1+zeta)**(cio_alpha-1.0)
       fact_length=cio_length

!$omp parallel do shared(n_dm,n_st,rxpa,rypa,rzpa,u2dm,u3dm,u4dm,masap,
!$omp+                   fact_length,fact_speed,fact_mass,cio_xc,cio_yc,
!$omp+                   cio_zc),
!$omp+            private(i),
!$omp+            default(none)
       do i=1,n_dm+n_st
        rxpa(i)=(rxpa(i)-cio_xc)*fact_length
        rypa(i)=(rypa(i)-cio_yc)*fact_length
        rzpa(i)=(rzpa(i)-cio_zc)*fact_length

        u2dm(i)=u2dm(i)*fact_speed
        u3dm(i)=u3dm(i)*fact_speed
        u4dm(i)=u4dm(i)*fact_speed

        masap(i)=masap(i)*fact_mass
       end do

       write(*,*)
       write(*,*) 'after unit conversion...'
       write(*,*) 'x positions (min,max), in Mpc:',
     &     minval(rxpa(1:n_dm+n_st)),maxval(rxpa(1:n_dm+n_st))
       write(*,*) 'y positions (min,max), in Mpc:',
     &     minval(rypa(1:n_dm+n_st)),maxval(rypa(1:n_dm+n_st))
       write(*,*) 'z positions (min,max), in Mpc:',
     &     minval(rzpa(1:n_dm+n_st)),maxval(rzpa(1:n_dm+n_st))
       write(*,*) 'x velocities (min,max), in c:',
     &     minval(u2dm(1:n_dm+n_st)),maxval(u2dm(1:n_dm+n_st))
       write(*,*) 'y velocities (min,max), in c:',
     &     minval(u3dm(1:n_dm+n_st)),maxval(u3dm(1:n_dm+n_st))
       write(*,*) 'z velocities (min,max), in c:',
     &     minval(u4dm(1:n_dm+n_st)),maxval(u4dm(1:n_dm+n_st))
       write(*,*) 'masses (min,max), in internal units:',
     &   minval(masap(1:n_dm+n_st)),maxval(masap(1:n_dm+n_st))
       write(*,*) 'unique ids (min,max):',
     &   minval(oripa(1:n_dm+n_st)),maxval(oripa(1:n_dm+n_st))
       write(*,*)
       return
       end