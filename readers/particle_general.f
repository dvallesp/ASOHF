*********************************************************************
      subroutine read_particles_general(iter,nx,ny,nz,t,zeta,
     &             u2dm,u3dm,u4dm,masap,rxpa,rypa,rzpa,oripa,
     &             n_dm,var,n_st,uv,um,hubble_littleh)
*********************************************************************
*      reads data from generic input format
*********************************************************************

      implicit none

      include 'input_files/asohf_parameters.dat'

      integer nx,ny,nz,iter,ndxyz
      real*4 t,aaa,bbb,ccc,map,zeta,uv,um,hubble_littleh
      integer var !(=1: only dm; =2: dm+stars)

      integer i,j,k,ix,nl,ir,irr,n_dm,n_st,nbas,nst0

      real*4 u2dm(parti_read)
      real*4 u3dm(parti_read)
      real*4 u4dm(parti_read)
      real*4 masap(parti_read)
      real*4 rxpa(parti_read)
      real*4 rypa(parti_read)
      real*4 rzpa(parti_read)

      integer oripa(parti_read)

      real*4 ubas(0:parti_read)
      integer ubasint(0:parti_read)

      real cio_mass,cio_speed,cio_length,cio_alpha,cio_xc,cio_yc,cio_zc
      common /conv_io/ cio_mass,cio_speed,cio_length,cio_alpha,cio_xc,
     &                 cio_yc,cio_zc
      real fact_mass,fact_speed,fact_length

      character*5 iter_string

*     reading data
      write(iter_string, '(i5.5)') iter !for saving files to disk
      write(*,*) 'reading iter',iter

**    dark matter
      open (32,file='./simulation/particles'//iter_string,
     &         status='unknown',action='read',form='unformatted')
      read(32) zeta
      read(32) n_dm

      read(32) (ubas(i),i=1,n_dm)
      rxpa(1:n_dm)=ubas(1:n_dm)
      read(32) (ubas(i),i=1,n_dm)
      rypa(1:n_dm)=ubas(1:n_dm)
      read(32) (ubas(i),i=1,n_dm)
      rzpa(1:n_dm)=ubas(1:n_dm)
      read(32) (ubas(i),i=1,n_dm)
      u2dm(1:n_dm)=ubas(1:n_dm)
      read(32) (ubas(i),i=1,n_dm)
      u3dm(1:n_dm)=ubas(1:n_dm)
      read(32) (ubas(i),i=1,n_dm)
      u4dm(1:n_dm)=ubas(1:n_dm)
      read(32) (ubas(i),i=1,n_dm)
      masap(1:n_dm)=ubas(1:n_dm)
      read(32) (ubasint(i),i=1,n_dm)
      oripa(1:n_dm)=ubasint(1:n_dm)

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
       read(32) n_st
       if (n_dm+n_st.gt.parti_read) then
        write(*,*) 'warning: bad dimensioning of parti_read',
     &              n_dm+n_st,'>',parti_read
        stop
       end if

       read(32) (ubas(i),i=1,n_st)
       rxpa(n_dm+1:n_dm+n_st)=ubas(1:n_st)
       read(32) (ubas(i),i=1,n_st)
       rypa(n_dm+1:n_dm+n_st)=ubas(1:n_st)
       read(32) (ubas(i),i=1,n_st)
       rzpa(n_dm+1:n_dm+n_st)=ubas(1:n_st)
       read(32) (ubas(i),i=1,n_st)
       u2dm(n_dm+1:n_dm+n_st)=ubas(1:n_st)
       read(32) (ubas(i),i=1,n_st)
       u3dm(n_dm+1:n_dm+n_st)=ubas(1:n_st)
       read(32) (ubas(i),i=1,n_st)
       u4dm(n_dm+1:n_dm+n_st)=ubas(1:n_st)
       read(32) (ubas(i),i=1,n_st)
       masap(n_dm+1:n_dm+n_st)=ubas(1:n_st)
       read(32) (ubasint(i),i=1,n_st)
       oripa(n_dm+1:n_dm+n_st)=ubasint(1:n_st)

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
      else
       n_st=0
      end if

      close(32)

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