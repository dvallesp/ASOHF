*********************************************************************
       subroutine read_gizmo_hdf5(iter,nx,ny,nz,t,zeta,
     &             u2dm,u3dm,u4dm,masap,rxpa,rypa,rzpa,oripa,
     &             n_dm,var,n_st,uv,um,hubble_littleh)
*********************************************************************
*      reads data from gadget unformatted inputs
*********************************************************************

       use hdf5
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

       integer i,j

*     io variables
       real*4,allocatable::scr41(:)
       real*4,allocatable::scr42(:,:)
       integer*4,allocatable::scrint1(:)
       real*8,allocatable::scr81(:)
       real*8,allocatable::scr82(:,:)
       
       type(c_ptr) :: redshift_ptr
       real(c_double), target :: redshift

       character*3 iter_string

       character*200 postfix
       common /POSTFIX/ postfix
       
       ! HDF5 variables
       integer(HID_T) :: file_id, dset_id, dtype_id, space_id, attr_id
       integer :: ierr
       integer(hsize_t), dimension(1) :: dims1d
       integer(hsize_t), dimension(2) :: dims2d, maxdims2d
       integer(hsize_t) :: precision

***** what i need to output/load is *******
*      rxpa,rypa,rzpa,u2dm,u3dm,u4dm,masap,oripa
*      n_dm, n_st, zeta
*****************************************

*     READING DATA
      WRITE(ITER_STRING, '(I3.3)') ITER !For saving files to disk
      WRITE(*,*) ' > Reading iter',ITER_STRING

      call h5open_f(ierr)
      call h5fopen_f(trim('./simulation/snap_'//ITER_STRING//postfix),
     &               H5F_ACC_RDONLY_F, file_id, ierr)

*     1st, get number of particles
      call h5dopen_f(file_id, "/PartType1/Coordinates", dset_id, ierr)
      call h5dget_space_f(dset_id, space_id, ierr)
      call h5sget_simple_extent_dims_f(space_id, dims2d, maxdims2d, 
     &                                 ierr)
      write(*,*) 'Number of DM particles in iter=',dims2d
      n_dm = dims2d(2)

       ! get the datatype and its precision (assume it's the same for all)
      call h5dget_type_f(dset_id, dtype_id, ierr)
      call h5tget_precision_f(dtype_id, precision, ierr)  ! bits
      call h5sclose_f(space_id, ierr)
      call h5dclose_f(dset_id, ierr)

      if (var.eq.2) then 
       call h5dopen_f(file_id, "/PartType4/Coordinates", dset_id, ierr)
       call h5dget_space_f(dset_id, space_id, ierr)
       call h5sget_simple_extent_dims_f(space_id, dims2d, maxdims2d, 
     &                                 ierr)
       write(*,*) 'Number of stellar particles in iter=',dims2d
       n_st = dims2d(2)
       call h5sclose_f(space_id, ierr)
       call h5dclose_f(dset_id, ierr)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   TO DO: should we add parttype 2 and 3 ??
!!!!   This may increase a lot the number of particles read, most of them being outside the domain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Obtain the redshift of the snapshot
      ! Redshift is an attribute of the header group, with name 'Redshift'
      call h5gopen_f(file_id, '/Header', dset_id, ierr)
      call h5aopen_f(dset_id, 'Redshift', attr_id, ierr)
      redshift_ptr = c_loc(redshift)
      call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, redshift_ptr, ierr)
      zeta = redshift
      ! Print the result
      print *, 'Redshift = ', zeta
      ! Close everything
      call h5aclose_f(attr_id, ierr)
      call h5dclose_f(dset_id, ierr)


*     2nd read DM (for now only PartType1)
       if (precision == 32) then
         call h5dopen_f(file_id, "/PartType1/Coordinates", dset_id, 
     &                  ierr)
         allocate(scr42(3,n_dm))
         dims2d(1) = 3
         dims2d(2) = n_dm
         call h5dread_f(dset_id, H5T_NATIVE_REAL, scr42, dims2d, ierr)
         rxpa(1:n_dm)=scr42(1,1:n_dm)
         rypa(1:n_dm)=scr42(2,1:n_dm)
         rzpa(1:n_dm)=scr42(3,1:n_dm)
         call h5dclose_f(dset_id, ierr)
         
         call h5dopen_f(file_id, "/PartType1/Velocities", dset_id, 
     &                  ierr)
         call h5dread_f(dset_id, H5T_NATIVE_REAL, scr42, dims2d, ierr)
         u2dm(1:n_dm)=scr42(1,1:n_dm)
         u3dm(1:n_dm)=scr42(2,1:n_dm)     
         u4dm(1:n_dm)=scr42(3,1:n_dm)
         call h5dclose_f(dset_id, ierr)
         deallocate(scr42)

         allocate(scr41(n_dm))
         dims1d(1) = n_dm
         call h5dopen_f(file_id, "/PartType1/Masses", dset_id, ierr)
         call h5dread_f(dset_id, H5T_NATIVE_REAL, scr41, dims1d, ierr)
         masap(1:n_dm)=scr41(1:n_dm)
         call h5dclose_f(dset_id, ierr)
         deallocate(scr41)

         allocate(scrint1(n_dm))
         call h5dopen_f(file_id, "/PartType1/ParticleIDs",
     &               dset_id, ierr)
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER, scrint1, dims1d,
     &               ierr)  
         oripa(1:n_dm)=scrint1(1:n_dm)
         call h5dclose_f(dset_id, ierr)
         deallocate(scrint1)

       else if (precision == 64) then
         call h5dopen_f(file_id, "/PartType1/Coordinates", dset_id, 
     &                  ierr)
         allocate(scr82(3,n_dm))
         dims2d(1) = 3
         dims2d(2) = n_dm
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, scr82, dims2d, 
     &                  ierr)
         rxpa(1:n_dm)=scr82(1,1:n_dm)
         rypa(1:n_dm)=scr82(2,1:n_dm)
         rzpa(1:n_dm)=scr82(3,1:n_dm)
         call h5dclose_f(dset_id, ierr)

         call h5dopen_f(file_id, "/PartType1/Velocities", dset_id, 
     &                  ierr)
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, scr82, dims2d, 
     &                  ierr)
         u2dm(1:n_dm)=scr82(1,1:n_dm)
         u3dm(1:n_dm)=scr82(2,1:n_dm)     
         u4dm(1:n_dm)=scr82(3,1:n_dm)
         call h5dclose_f(dset_id, ierr)
         deallocate(scr82)

         allocate(scr81(n_dm))
         dims1d(1) = n_dm
         call h5dopen_f(file_id, "/PartType1/Masses", dset_id, ierr)
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, scr81, dims1d, 
     &                  ierr)
         masap(1:n_dm)=scr81(1:n_dm)
         call h5dclose_f(dset_id, ierr)
         deallocate(scr81)

         allocate(scrint1(n_dm))
         call h5dopen_f(file_id, "/PartType1/ParticleIDs",
     &               dset_id, ierr)
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER, scrint1, dims1d,
     &               ierr)
         oripa(1:n_dm)=scrint1(1:n_dm)
         call h5dclose_f(dset_id, ierr)
         deallocate(scrint1)
       end if

*      3rd read stars (if necessary)
       if (var.eq.2.and. n_st.gt.0) then
        if (precision == 32) then
          call h5dopen_f(file_id, "/PartType4/Coordinates", dset_id, 
     &                  ierr)
          allocate(scr42(3,n_st))
          dims2d(1) = 3
          dims2d(2) = n_st
          call h5dread_f(dset_id, H5T_NATIVE_REAL, scr42, dims2d, ierr)
          rxpa(n_dm+1:n_dm+n_st)=scr42(1,1:n_st)
          rypa(n_dm+1:n_dm+n_st)=scr42(2,1:n_st)
          rzpa(n_dm+1:n_dm+n_st)=scr42(3,1:n_st)
          call h5dclose_f(dset_id, ierr)
          
          call h5dopen_f(file_id, "/PartType4/Velocities", dset_id, 
     &                  ierr)
          call h5dread_f(dset_id, H5T_NATIVE_REAL, scr42, dims2d, ierr)
          u2dm(n_dm+1:n_dm+n_st)=scr42(1,1:n_st)
          u3dm(n_dm+1:n_dm+n_st)=scr42(1,1:n_st)     
          u4dm(n_dm+1:n_dm+n_st)=scr42(1,1:n_st)
          call h5dclose_f(dset_id, ierr)
          deallocate(scr42)
 
          allocate(scr41(n_st))
          dims1d(1) = n_st
          call h5dopen_f(file_id, "/PartType4/Masses", dset_id, ierr)
          call h5dread_f(dset_id, H5T_NATIVE_REAL, scr41, dims1d, ierr)
          masap(n_dm+1:n_dm+n_st)=scr41(1:n_st)
          call h5dclose_f(dset_id, ierr)
          deallocate(scr41)
 
          allocate(scrint1(n_st))
          call h5dopen_f(file_id, "/PartType4/ParticleIDs",
     &                dset_id, ierr)
          call h5dread_f(dset_id, H5T_NATIVE_INTEGER, scrint1, dims1d,
     &                ierr)  
          oripa(n_dm+1:n_dm+n_st)=scrint1(1:n_st)
          call h5dclose_f(dset_id, ierr)
          deallocate(scrint1)
        else if (precision == 64) then
          call h5dopen_f(file_id, "/PartType4/Coordinates", dset_id, 
     &                  ierr)
          allocate(scr82(3,n_st))
          dims2d(1) = 3
          dims2d(2) = n_st
          call h5dread_f(dset_id, H5T_NATIVE_REAL, scr82, dims2d, ierr)
          rxpa(n_dm+1:n_dm+n_st)=scr82(1,1:n_st)
          rypa(n_dm+1:n_dm+n_st)=scr82(2,1:n_st)
          rzpa(n_dm+1:n_dm+n_st)=scr82(3,1:n_st)
          call h5dclose_f(dset_id, ierr)
          
          call h5dopen_f(file_id, "/PartType4/Velocities", dset_id, 
     &                  ierr)
          call h5dread_f(dset_id, H5T_NATIVE_REAL, scr82, dims2d, ierr)
          u2dm(n_dm+1:n_dm+n_st)=scr82(1,1:n_st)
          u3dm(n_dm+1:n_dm+n_st)=scr82(1,1:n_st)     
          u4dm(n_dm+1:n_dm+n_st)=scr82(1,1:n_st)
          call h5dclose_f(dset_id, ierr)
          deallocate(scr82)
 
          allocate(scr81(n_st))
          dims1d(1) = n_st
          call h5dopen_f(file_id, "/PartType4/Masses", dset_id, ierr)
          call h5dread_f(dset_id, H5T_NATIVE_REAL, scr81, dims1d, ierr)
          masap(n_dm+1:n_dm+n_st)=scr81(1:n_st)
          call h5dclose_f(dset_id, ierr)
          deallocate(scr81)
 
          allocate(scrint1(n_st))
          call h5dopen_f(file_id, "/PartType4/ParticleIDs",
     &                dset_id, ierr)
          call h5dread_f(dset_id, H5T_NATIVE_INTEGER, scrint1, dims1d,
     &                ierr)  
          oripa(n_dm+1:n_dm+n_st)=scrint1(1:n_st)
          call h5dclose_f(dset_id, ierr)
          deallocate(scrint1)
        end if
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