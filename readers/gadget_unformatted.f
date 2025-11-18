*********************************************************************
      subroutine read_gadget_unformatted(iter,nx,ny,nz,t,zeta,
     &             u2dm,u3dm,u4dm,masap,rxpa,rypa,rzpa,oripa,
     &             n_dm,var,n_st,uv,um,hubble_littleh)
*********************************************************************
*      reads data from gadget unformatted inputs
*********************************************************************

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

*     local variables
      integer ntot,igas0,igas1,idm0,idm1,ist0,ist1,ibh0,ibh1,i,j
      integer igas0m,igas1m,idm0m,idm1m,ist0m,ist1m,ibh0m,ibh1m ! indices for the mass arrays (may be different from the position, velocity arrays!!)

*     io variables
      character*4 blocklabel
      integer*4 blocksize,npp(6)
      real*8 mass_arr(6)
      real*8 time8,zeta8
      real*8 caca(5)
      real*8 boxsize8,omega_m8,omega_lambda8,hubble_param8

      real*4,allocatable::scr41(:)
      real*4,allocatable::scr42(:,:)
      integer*4,allocatable::scrint1(:)

      character*3 iter_string
***** fix unit conversions for gadget *******
      !cio_mass=1.e10/hubble_littleh
      !cio_speed=299792.458
      !cio_length=1.0/hubble_littleh
      !cio_alpha=0.5
*********************************************

*     reading data
      write(iter_string, '(i3.3)') iter !for saving files to disk
      write(*,*) ' > reading iter',iter_string

      open(11, file='./simulation/snap_'//iter_string,
     &     status='unknown',action='read', form='unformatted')

*      read the header ************************************************
       read(11) blocklabel,blocksize
       !write(*,*) 'found block ', blocklabel, ' with length', blocksize
       read(11) npp,mass_arr,time8,zeta8,caca,boxsize8,omega_m8,
     &          omega_lambda8,hubble_param8

       zeta=zeta8

       !write(*,*) npp
       !write(*,*) mass_arr
       !write(*,*) 'redshift:', zeta8
       !write(*,*) 'box size (ckpc/h):', boxsize8
       !write(*,*) 'om, olambda, h:', omega_m8, omega_lambda8,
       !&            hubble_param8

       ntot=sum(npp)

       igas0=1
       igas1=npp(1)
       idm0=igas1+1
       idm1=sum(npp(1:4))
       ist0=idm1+1
       ist1=sum(npp(1:5))
       ibh0=ist1+1
       ibh1=sum(npp(1:6))
       !write(*,*) 'gas indices',igas0,igas1,igas1-igas0+1
       !write(*,*) 'dm indices ',idm0,idm1,idm1-idm0+1
       !write(*,*) 'st indices ',ist0,ist1,ist1-ist0+1
       !write(*,*) 'bh indices ',ibh0,ibh1,ibh1-ibh0+1

       if (mass_arr(2).gt.0.0) then
        ! no DM particle masses written explicitlym they are in the header
        write(*,*) 'dm masses from the header!'
        igas0m=igas0
        igas1m=igas1 
        idm0m=-1
        idm1m=-1
        ist0m=igas1m+1 
        ist1m=igas1m+npp(5)
        ibh0m=ist1m+1
        ibh1m=ist1m+npp(6)
       else
        write(*,*) 'dm masses from the file content!'
        igas0m=igas0
        igas1m=igas1
        idm0m=idm0
        idm1m=idm1
        ist0m=ist0
        ist1m=ist1
        ibh0m=ibh0
        ibh1m=ibh1
       end if

       n_dm=idm1-idm0+1
       if (var.eq.2) then
        n_st=ist1-ist0+1
       else
        n_st=0
       end if

*      read the particle positions ************************************
       read(11) blocklabel,blocksize
       write(*,*) 'found block ', blocklabel, ' with length', blocksize
       allocate(scr42(3,ntot))
       read(11) ((scr42(j,i),j=1,3),i=1,ntot) ! all particles
       !write(*,*) '-x-', minval(scr42(1,:)), maxval(scr42(1,:))
       !write(*,*) '-y-', minval(scr42(2,:)), maxval(scr42(2,:))
       !write(*,*) '-z-', minval(scr42(3,:)), maxval(scr42(3,:))
       rxpa(1:n_dm)=scr42(1,idm0:idm1)
       rypa(1:n_dm)=scr42(2,idm0:idm1)
       rzpa(1:n_dm)=scr42(3,idm0:idm1)
       if (var.eq.2) then
        rxpa(n_dm+1:n_dm+n_st)=scr42(1,ist0:ist1)
        rypa(n_dm+1:n_dm+n_st)=scr42(2,ist0:ist1)
        rzpa(n_dm+1:n_dm+n_st)=scr42(3,ist0:ist1)
       end if

*      read the particle velocities************************************
       read(11) blocklabel,blocksize
       write(*,*) 'found block ', blocklabel, ' with length', blocksize
       read(11) ((scr42(j,i),j=1,3),i=1,ntot) ! all particles
       !write(*,*) '-vx-', minval(scr42(1,:)), maxval(scr42(1,:))
       !write(*,*) '-vy-', minval(scr42(2,:)), maxval(scr42(2,:))
       !write(*,*) '-vz-', minval(scr42(3,:)), maxval(scr42(3,:))
       u2dm(1:n_dm)=scr42(1,idm0:idm1)
       u3dm(1:n_dm)=scr42(2,idm0:idm1)
       u4dm(1:n_dm)=scr42(3,idm0:idm1)
       if (var.eq.2) then
        u2dm(n_dm+1:n_dm+n_st)=scr42(1,ist0:ist1)
        u3dm(n_dm+1:n_dm+n_st)=scr42(2,ist0:ist1)
        u4dm(n_dm+1:n_dm+n_st)=scr42(3,ist0:ist1)
       end if
       deallocate(scr42)

*      read the particle ids****************************************
       read(11) blocklabel,blocksize
       write(*,*) 'found block ', blocklabel, ' with length', blocksize
       allocate(scrint1(ntot))
       read(11) (scrint1(i),i=1,ntot) ! all particles
       !write(*,*) '-id-', minval(scrint1(:)), maxval(scrint1(:))
       oripa(1:n_dm)=scrint1(idm0:idm1)
       if (var.eq.2) then
        oripa(n_dm+1:n_dm+n_st)=scrint1(ist0:ist1)
       end if
       deallocate(scrint1)

       ! Compatibility with 64-bit IDs: I will just number the particles in the input file
       write(*,*) 'WARNING / NOTE: the unique ids of the particles' 
       write(*,*) ' are not read from the input file, but just numbered'
       write(*,*) ' from 1 to n_dm+n_st'
       oripa(1:n_dm)=(/ (idm0+i-1,i=1,n_dm) /)
       if (var.eq.2) then
        oripa(n_dm+1:n_dm+n_st)=(/ (ist0+i-1,i=1,n_st) /)
       end if

*      read the particle masses****************************************
       read(11) blocklabel,blocksize
       write(*,*) 'found block ', blocklabel, ' with length', blocksize
       allocate(scr41(ntot))
       if (idm0m.gt.0) then
        read(11) (scr41(i),i=1,ntot) ! all particles
       else
        read(11) (scr41(i),i=1,ntot-n_dm) ! all particles
       end if
       !write(*,*) '-m-', minval(scr41(:)), maxval(scr41(:))
       if (idm0m.gt.0) then
        masap(1:n_dm)=scr41(idm0m:idm1m)
       else
        masap(1:n_dm)=mass_arr(2) ! all DM particles have the same mass
       end if
       if (var.eq.2) then
        masap(n_dm+1:n_dm+n_st)=scr41(ist0m:ist1m)
       end if
       deallocate(scr41)

      close(11)

      if (n_dm+n_st.gt.parti_read) then
       write(*,*) 'warning: bad dimensioning of parti_read',
     &            n_dm+n_st,'>',parti_read
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
      write(*,*) 'x positions (min,max), in mpc:',
     &     minval(rxpa(1:n_dm+n_st)),maxval(rxpa(1:n_dm+n_st))
      write(*,*) 'y positions (min,max), in mpc:',
     &     minval(rypa(1:n_dm+n_st)),maxval(rypa(1:n_dm+n_st))
      write(*,*) 'z positions (min,max), in mpc:',
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
