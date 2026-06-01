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

*********************************************************************
       subroutine read_enzo_hdf5(iter,nx,ny,nz,t,zeta,
     &             u2dm,u3dm,u4dm,masap,rxpa,rypa,rzpa,oripa,
     &             n_dm,var,n_st,uv,um,hubble_littleh)
*********************************************************************
*      reads data from ENZO hdf5 inputs
*********************************************************************

       use hdf5
       use hdf5_readers
       use iso_c_binding
       implicit none

       include 'input_files/asohf_parameters.dat'

       integer nx,ny,nz,iter,ndxyz
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

       integer i,j,n

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

       character*200 postfix,enzofolderprefix,enzoprefix
       character*200 enzoparamspostfix
       common /POSTFIX/ postfix,enzofolderprefix,enzoprefix,
     &                  enzoparamspostfix
       character*200 filename,filename_params

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
      write(*,*) 'Only implemented for uniform grid case, and no stars'
      write(*,*) 'And everything written in the cpu0000 file'
      write(*,*) '!!!!'

*     READING DATA
      WRITE(ITER_STRING, '(I4.4)') ITER !For saving files to disk
      WRITE(*,*) ' > Reading iter ',ITER_STRING

      call h5open_f(ierr)
      filename = './simulation/' // trim(enzofolderprefix) //
     &       trim(ITER_STRING) // '/' // trim(enzoprefix) // 
     &       trim(ITER_STRING) // trim(postfix)
      write(*,*) 'filename=',filename
      call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr)

*     1st, get number of particles
      call h5gopen_f(file_id, "/Grid00000001", group_id, ierr)
      call h5dopen_f(group_id, "particle_position_x", dset_id, ierr)

      call h5dget_space_f(dset_id, space_id, ierr)
      call h5sget_simple_extent_dims_f(space_id, dims1d, maxdims1d, 
     &                                 ierr)
      write(*,*) 'Number of DM particles in iter=',dims1d
      n_dm = dims1d(1)

      ! get the datatype and its precision (assume it's the same for all)
      call h5dget_type_f(dset_id, dtype_id, ierr)
      call h5tget_precision_f(dtype_id, precision, ierr)  ! bits
      call h5sclose_f(space_id, ierr)
      call h5dclose_f(dset_id, ierr)

      if (var.eq.2) then 
C       call h5dopen_f(file_id, "/PartType4/Coordinates", dset_id, ierr)
C       call h5dget_space_f(dset_id, space_id, ierr)
C       call h5sget_simple_extent_dims_f(space_id, dims2d, maxdims2d, 
C     &                                 ierr)
C       write(*,*) 'Number of stellar particles in iter=',dims2d
C       n_st = dims2d(2)
C       call h5sclose_f(space_id, ierr)
C       call h5dclose_f(dset_id, ierr)
        write(*,*) 'WARNING: reading stars not implemented yet'
        stop
      else 
       n_st = 0
      end if

      write(*,*) 'Precision of data in file: ', precision, ' bits'
      if (precision.ne.32 .and. precision.ne.64) then
       write(*,*) 'Unsupported precision: ', precision
       stop
      end if

      if (n_dm+n_st.gt.parti_read) then
       write(*,*) 'warning: bad dimensioning of parti_read',
     &            n_dm+n_st,'>',parti_read
       stop
      end if

      call h5gclose_f(group_id, ierr)
      call h5fclose_f(file_id, ierr)

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

*     2nd read DM (for now only first grid)
      call h5open_f(ierr)
      call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr)
      call h5gopen_f(file_id, "/Grid00000001", group_id, ierr)

      if (precision == 32) then 
        call read1d_32(group_id, "particle_position_x", scr41, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        rxpa(1:n_dm)=scr41(1:n_dm)

        call read1d_32(group_id, "particle_position_y", scr41, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        rypa(1:n_dm)=scr41(1:n_dm)

        call read1d_32(group_id, "particle_position_z", scr41, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        rzpa(1:n_dm)=scr41(1:n_dm)

        call read1d_32(group_id, "particle_velocity_x", scr41, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        u2dm(1:n_dm)=scr41(1:n_dm)

        call read1d_32(group_id, "particle_velocity_y", scr41, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        u3dm(1:n_dm)=scr41(1:n_dm)

        call read1d_32(group_id, "particle_velocity_z", scr41, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        u4dm(1:n_dm)=scr41(1:n_dm)

        call read1d_32(group_id, "particle_mass", scr41, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        masap(1:n_dm)=scr41(1:n_dm)

        call readint1d_32(group_id, "particle_id", scrint41, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        oripa(1:n_dm)=scrint41(1:n_dm)
        deallocate(scr41, scrint41)

      else if (precision == 64) then 

        call read1d_64(group_id, "particle_position_x", scr81, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        rxpa(1:n_dm)=scr81(1:n_dm)

        call read1d_64(group_id, "particle_position_y", scr81, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        rypa(1:n_dm)=scr81(1:n_dm)

        call read1d_64(group_id, "particle_position_z", scr81, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        rzpa(1:n_dm)=scr81(1:n_dm)

        call read1d_64(group_id, "particle_velocity_x", scr81, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        u2dm(1:n_dm)=scr81(1:n_dm)

        call read1d_64(group_id, "particle_velocity_y", scr81, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        u3dm(1:n_dm)=scr81(1:n_dm)

        call read1d_64(group_id, "particle_velocity_z", scr81, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        u4dm(1:n_dm)=scr81(1:n_dm)

        call read1d_64(group_id, "particle_mass", scr81, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        masap(1:n_dm)=scr81(1:n_dm)

        call readint1d_64(group_id, "particle_index", scrint81, n, ierr)
        if (n.ne.n_dm) then
         write(*,*) 'Error: number of particles read != n_dm'
         stop
        end if
        oripa(1:n_dm)=scrint81(1:n_dm)
        deallocate(scr81, scrint81)
      end if

*      3rd read stars (if necessary)
       if (var.eq.2.and. n_st.gt.0) then
        write(*,*) 'WARNING: reading stars not implemented yet'
        stop
       end if
         
       write(*,*)
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