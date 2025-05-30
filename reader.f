#if reader == 1
#include "readers/masclet.f"
#endif

#if reader == 0
#include "readers/particle_general.f"
#endif

#if reader == 2
#include "readers/gadget_unformatted.f"
#endif

#if reader == 3
#include "readers/gizmo_hdf5.f"
#endif

*********************************************************************
       subroutine sort_dm_particles(n_dm,npart_esp,n_st,ir_kern_stars)
*********************************************************************
*      reorders dm particles by species (assumes there are n_esp
*       especies, each 8 times lighter than the previous one)
*********************************************************************
       use particles
       implicit none
       include 'input_files/asohf_parameters.dat'

       integer n_dm
       integer npart_esp(0:n_esp-1)
       integer n_st,ir_kern_stars

       integer i,j,k,n,iesp,conta
       real mlow,mhigh,bas,maxmass,minmass
       integer,allocatable::indices(:)
       real,allocatable::scr(:,:)
       integer,allocatable::scrint(:,:)

       write(*,*) 'sorting particles by mass'

       maxmass=maxval(masap(1:n_dm))
       minmass=minval(masap(1:n_dm))
       npart_esp=0

       conta=0
       allocate(indices(1:n_dm))
       do iesp=0,n_esp-1
        mhigh=2*maxmass/8.0**iesp
        mlow=0.5*maxmass/8.0**iesp
        if (mhigh.lt.minmass) exit

        do i=1,n_dm
         bas=masap(i)
         if (bas.lt.mhigh) then
         if (bas.gt.mlow) then
          conta=conta+1
          indices(conta)=i
         end if
         end if
        end do
        if (iesp.eq.0) then
         npart_esp(iesp)=conta
        else
         npart_esp(iesp)=conta-sum(npart_esp(0:iesp-1))
        end if
        if (npart_esp(iesp).gt.0) then
         write(*,*) 'of species',iesp,', no. particles:',npart_esp(iesp)
        end if
       end do

       if (conta.ne.n_dm.or.sum(npart_esp(0:n_esp-1)).ne.n_dm) then
        write(*,*) 'wrong sorting, cannot continue',conta,n_dm
        stop
       end if

       allocate(scr(7,n_dm),scrint(1,n_dm))

!$omp parallel do shared(scr,scrint,rxpa,rypa,rzpa,u2dm,u3dm,u4dm,masap,
!$omp+                   oripa,indices,n_dm),
!$omp+            private(i),
!$omp+            default(none)
       do i=1,n_dm
        scr(1,i)=rxpa(indices(i))
        scr(2,i)=rypa(indices(i))
        scr(3,i)=rzpa(indices(i))
        scr(4,i)=u2dm(indices(i))
        scr(5,i)=u3dm(indices(i))
        scr(6,i)=u4dm(indices(i))
        scr(7,i)=masap(indices(i))
        scrint(1,i)=oripa(indices(i))
       end do

       deallocate(indices)

!$omp parallel do shared(scr,scrint,rxpa,rypa,rzpa,u2dm,u3dm,u4dm,masap,
!$omp+                   oripa,indices,n_dm),
!$omp+            private(i),
!$omp+            default(none)
       do i=1,n_dm
        rxpa(i)=scr(1,i)
        rypa(i)=scr(2,i)
        rzpa(i)=scr(3,i)
        u2dm(i)=scr(4,i)
        u3dm(i)=scr(5,i)
        u4dm(i)=scr(6,i)
        masap(i)=scr(7,i)
        oripa(i)=scrint(1,i)
       end do

       deallocate(scr,scrint)

       if (n_st.gt.0) then
        
        if (ir_kern_stars.gt.n_esp-1) then
         write(*,*) 'warning: ir_kern_stars > n_esp-1',
     &              ir_kern_stars,n_esp-1
         write(*,*) 'fix n_esp to ir_kern_stars+1, at least'
         stop
        end if

        npart_esp(ir_kern_stars)=npart_esp(ir_kern_stars)+n_st
        write(*,*) 'stars: of species',ir_kern_stars,
     &             ', no. particles:',npart_esp(ir_kern_stars)
       end if

c       write(*,*) 'checking...'
*      check
c       bas=masap(1)
c       do i=2,n_dm
c        if (masap(i).gt.1.0001*bas) then
c         write(*,*) 'wrong, i=',i
c         stop
c        end if
c        bas=masap(i)
c       end do
c       stop

       return
       end

*********************************************************************
       subroutine sort_dm_particles_localdensity(n_dm,npart_esp,n_st,
     &                                           ir_kern_stars,rodo,re0)
*********************************************************************
*      assigns dm particles in species, according to local density,
*       and reorders them accordingly.
*********************************************************************
       use particles
       implicit none
       include 'input_files/asohf_parameters.dat'

       integer n_dm
       integer npart_esp(0:n_esp-1)
       integer n_st,ir_kern_stars
       real*4 rodo,re0

       real*4 dx,dy,dz
       common /espaciado/ dx,dy,dz

       integer i,j,k,n,iesp,conta,nx,ny,nz,ix,jy,kz,plev
       real mlow,mhigh,bas,maxmass,minmass
       real xmin,ymin,zmin
       integer,allocatable::indices(:)
       real,allocatable::scr(:,:)
       integer,allocatable::scrint(:,:)
       real,allocatable::dens(:,:,:)
       integer,allocatable::mocklevel(:)

       nx=nmax
       ny=nmay
       nz=nmaz

       allocate(dens(nx,ny,nz))
!$omp parallel do shared(nx,ny,nz,dens),
!$omp+            private(ix,jy,kz),
!$omp+            default(none)
       do kz=1,nz
       do jy=1,ny
       do ix=1,nx
        dens(ix,jy,kz)=0.0
       end do
       end do
       end do

       xmin=-nx*dx/2.0
       ymin=-ny*dy/2.0
       zmin=-nz*dz/2.0
!$omp parallel do shared(n_dm,rxpa,rypa,rzpa,xmin,ymin,zmin,dx,dy,dz,
!$omp+                   masap,nx,ny,nz),
!$omp+            private(ix,jy,kz,i),
!$omp+            reduction(+:dens), default(none)
       do i=1,n_dm
        ix=int((rxpa(i)-xmin)/dx)+1
        jy=int((rypa(i)-ymin)/dy)+1
        kz=int((rzpa(i)-zmin)/dz)+1
        if (ix.lt.1) ix=1
        if (ix.gt.nx) ix=nx
        if (jy.lt.1) jy=1
        if (jy.gt.ny) jy=ny
        if (kz.lt.1) kz=1
        if (kz.gt.nz) kz=nz
        dens(ix,jy,kz)=dens(ix,jy,kz)+masap(i)
       end do

       bas=dx*dy*dz*rodo*re0**3
!$omp parallel do shared(nx,ny,nz,dens,bas),
!$omp+            private(ix,jy,kz),
!$omp+            default(none)
       do kz=1,nz
       do jy=1,ny
       do ix=1,nx
        dens(ix,jy,kz)=dens(ix,jy,kz)/bas
       end do
       end do
       end do

       allocate(mocklevel(n_dm))

!$omp parallel do shared(n_dm,rxpa,rypa,rzpa,xmin,ymin,zmin,dx,dy,dz,
!$omp+                   nx,ny,nz,dens,mocklevel),
!$omp+            private(ix,jy,kz,i,bas),
!$omp+            default(none)
       do i=1,n_dm
        ix=int((rxpa(i)-xmin)/dx)+1
        jy=int((rypa(i)-ymin)/dy)+1
        kz=int((rzpa(i)-zmin)/dz)+1
        if (ix.lt.1) ix=1
        if (ix.gt.nx) ix=nx
        if (jy.lt.1) jy=1
        if (jy.gt.ny) jy=ny
        if (kz.lt.1) kz=1
        if (kz.gt.nz) kz=nz
        bas=dens(ix,jy,kz)
        if (bas.gt.0.0) then
         mocklevel(i)=max(min(int(log(bas)/log(8.0)),n_esp-1),0)
        else
         mocklevel(i)=0
        end if
       end do

       deallocate(dens)

       write(*,*) 'sorting particles by local density'

       npart_esp=0
       conta=0

       allocate(indices(1:n_dm))
       do iesp=0,n_esp-1
        do i=1,n_dm
         plev=mocklevel(i)
         if (plev.eq.iesp) then
          conta=conta+1
          indices(conta)=i
         end if
        end do
        if (iesp.eq.0) then
         npart_esp(iesp)=conta
        else
         npart_esp(iesp)=conta-sum(npart_esp(0:iesp-1))
        end if
        if (npart_esp(iesp).gt.0) then
         write(*,*) 'of species',iesp,', no. particles:',npart_esp(iesp)
        end if
       end do

       deallocate(mocklevel)

       if (conta.ne.n_dm.or.sum(npart_esp(0:n_esp-1)).ne.n_dm) then
        write(*,*) 'wrong sorting, cannot continue',conta,n_dm
        stop
       end if

       allocate(scr(7,n_dm),scrint(1,n_dm))

!$omp parallel do shared(scr,scrint,rxpa,rypa,rzpa,u2dm,u3dm,u4dm,masap,
!$omp+                   oripa,indices,n_dm),
!$omp+            private(i),
!$omp+            default(none)
       do i=1,n_dm
        scr(1,i)=rxpa(indices(i))
        scr(2,i)=rypa(indices(i))
        scr(3,i)=rzpa(indices(i))
        scr(4,i)=u2dm(indices(i))
        scr(5,i)=u3dm(indices(i))
        scr(6,i)=u4dm(indices(i))
        scr(7,i)=masap(indices(i))
        scrint(1,i)=oripa(indices(i))
       end do

       deallocate(indices)

!$omp parallel do shared(scr,scrint,rxpa,rypa,rzpa,u2dm,u3dm,u4dm,masap,
!$omp+                   oripa,indices,n_dm),
!$omp+            private(i),
!$omp+            default(none)
       do i=1,n_dm
        rxpa(i)=scr(1,i)
        rypa(i)=scr(2,i)
        rzpa(i)=scr(3,i)
        u2dm(i)=scr(4,i)
        u3dm(i)=scr(5,i)
        u4dm(i)=scr(6,i)
        masap(i)=scr(7,i)
        oripa(i)=scrint(1,i)
       end do

       deallocate(scr,scrint)

       if (n_st.gt.0) then

        if (ir_kern_stars.gt.n_esp-1) then
         write(*,*) 'warning: ir_kern_stars > n_esp-1',
     &              ir_kern_stars,n_esp-1
         write(*,*) 'fix n_esp to ir_kern_stars+1, at least'
         stop
        end if


        npart_esp(ir_kern_stars)=npart_esp(ir_kern_stars)+n_st
        write(*,*) 'stars: of species',ir_kern_stars,
     &             ', no. particles:',npart_esp(ir_kern_stars)
       end if

c       write(*,*) 'checking...'
*      check
c       bas=masap(1)
c       do i=2,n_dm
c        if (masap(i).gt.1.0001*bas) then
c         write(*,*) 'wrong, i=',i
c         stop
c        end if
c        bas=masap(i)
c       end do
c       stop

       return
       end

*********************************************************************
       subroutine sort_dm_particles_localdensity_and_mass(n_dm,
     &                   npart_esp,n_st,ir_kern_stars,rodo,re0)
*********************************************************************
*      assigns dm particles in species, according to local density,
*       and reorders them accordingly.
*********************************************************************
       use particles
       implicit none
       include 'input_files/asohf_parameters.dat'

       integer n_dm
       integer npart_esp(0:n_esp-1)
       integer n_st,ir_kern_stars
       real*4 rodo,re0

       real*4 dx,dy,dz
       common /espaciado/ dx,dy,dz

       integer i,j,k,n,iesp,conta,nx,ny,nz,ix,jy,kz,plev,nn,ii,jj,kk
       integer fac_grid,maxlev
       real mlow,mhigh,bas,maxmass,minmass,fac,levmass,levdens
       real xmin,ymin,zmin,dxpa,dypa,dzpa
       integer,allocatable::indices(:)
       real,allocatable::scr(:,:)
       integer,allocatable::scrint(:,:)
       real,allocatable::dens(:,:,:)
       integer,allocatable::mocklevel(:)

       maxmass=maxval(masap(1:n_dm))
       minmass=minval(masap(1:n_dm))
       maxlev=1 ! particles at this level and above are pointlike, lower bigger cloud
       fac_grid=2 ! we consider this factor the base grid


       nx=nmax*fac_grid
       ny=nmay*fac_grid
       nz=nmaz*fac_grid

       dxpa=dx/fac_grid
       dypa=dy/fac_grid
       dzpa=dz/fac_grid

       allocate(dens(nx,ny,nz))
!$omp parallel do shared(nx,ny,nz,dens),
!$omp+            private(ix,jy,kz),
!$omp+            default(none)
       do kz=1,nz
       do jy=1,ny
       do ix=1,nx
        dens(ix,jy,kz)=0.0
       end do
       end do
       end do

       xmin=-nx*dxpa/2.0
       ymin=-ny*dypa/2.0
       zmin=-nz*dzpa/2.0
!$omp parallel do shared(n_dm,rxpa,rypa,rzpa,xmin,ymin,zmin,dxpa,dypa,
!$omp+                   dzpa,masap,nx,ny,nz,maxmass,maxlev),
!$omp+            private(ix,jy,kz,i,plev,ii,jj,kk,nn,fac),
!$omp+            reduction(+:dens), default(none)
       do i=1,n_dm
        ix=int((rxpa(i)-xmin)/dxpa)+1
        jy=int((rypa(i)-ymin)/dypa)+1
        kz=int((rzpa(i)-zmin)/dzpa)+1
        if (ix.lt.1) ix=1
        if (ix.gt.nx) ix=nx
        if (jy.lt.1) jy=1
        if (jy.gt.ny) jy=ny
        if (kz.lt.1) kz=1
        if (kz.gt.nz) kz=nz
        plev=int(log(maxmass/masap(i))+0.01)
        nn=2**max(maxlev-plev,0)-1
        fac=1./float((2*nn+1)**3)
        do kk=-nn,nn
        do jj=-nn,nn
        do ii=-nn,nn
         dens(ix,jy,kz)=dens(ix,jy,kz)+fac*masap(i)
        end do
        end do
        end do
       end do

       bas=dx*dy*dz*rodo*re0**3
!$omp parallel do shared(nx,ny,nz,dens,bas),
!$omp+            private(ix,jy,kz),
!$omp+            default(none)
       do kz=1,nz
       do jy=1,ny
       do ix=1,nx
        dens(ix,jy,kz)=dens(ix,jy,kz)/bas
       end do
       end do
       end do

       allocate(mocklevel(n_dm))

!$omp parallel do shared(n_dm,rxpa,rypa,rzpa,xmin,ymin,zmin,dxpa,dypa,
!$omp+                   dzpa,nx,ny,nz,dens,mocklevel,maxmass,masap),
!$omp+            private(ix,jy,kz,i,bas,levmass,levdens),
!$omp+            default(none)
       do i=1,n_dm
        ix=int((rxpa(i)-xmin)/dxpa)+1
        jy=int((rypa(i)-ymin)/dypa)+1
        kz=int((rzpa(i)-zmin)/dzpa)+1
        if (ix.lt.1) ix=1
        if (ix.gt.nx) ix=nx
        if (jy.lt.1) jy=1
        if (jy.gt.ny) jy=ny
        if (kz.lt.1) kz=1
        if (kz.gt.nz) kz=nz
        bas=dens(ix,jy,kz)
        if (bas.gt.0.0) then
         levmass=log(maxmass/masap(i))/log(8.0)
         levdens=log(bas)/log(8.0)
         mocklevel(i)=max(min(
     &                   int(max(levmass,levdens)+0.01),
     &                n_esp-1),0)
        else
         mocklevel(i)=0
        end if
       end do

       deallocate(dens)

       write(*,*) 'sorting particles by local density + particle mass'

       npart_esp=0
       conta=0

       allocate(indices(1:n_dm))
       do iesp=0,n_esp-1
        do i=1,n_dm
         plev=mocklevel(i)
         if (plev.eq.iesp) then
          conta=conta+1
          indices(conta)=i
         end if
        end do
        if (iesp.eq.0) then
         npart_esp(iesp)=conta
        else
         npart_esp(iesp)=conta-sum(npart_esp(0:iesp-1))
        end if
        if (npart_esp(iesp).gt.0) then
         write(*,*) 'of species',iesp,', no. particles:',npart_esp(iesp)
        end if
       end do

       deallocate(mocklevel)

       if (conta.ne.n_dm.or.sum(npart_esp(0:n_esp-1)).ne.n_dm) then
        write(*,*) 'wrong sorting, cannot continue',conta,n_dm
        stop
       end if

       allocate(scr(7,n_dm),scrint(1,n_dm))

!$omp parallel do shared(scr,scrint,rxpa,rypa,rzpa,u2dm,u3dm,u4dm,masap,
!$omp+                   oripa,indices,n_dm),
!$omp+            private(i),
!$omp+            default(none)
       do i=1,n_dm
        scr(1,i)=rxpa(indices(i))
        scr(2,i)=rypa(indices(i))
        scr(3,i)=rzpa(indices(i))
        scr(4,i)=u2dm(indices(i))
        scr(5,i)=u3dm(indices(i))
        scr(6,i)=u4dm(indices(i))
        scr(7,i)=masap(indices(i))
        scrint(1,i)=oripa(indices(i))
       end do

       deallocate(indices)

!$omp parallel do shared(scr,scrint,rxpa,rypa,rzpa,u2dm,u3dm,u4dm,masap,
!$omp+                   oripa,indices,n_dm),
!$omp+            private(i),
!$omp+            default(none)
       do i=1,n_dm
        rxpa(i)=scr(1,i)
        rypa(i)=scr(2,i)
        rzpa(i)=scr(3,i)
        u2dm(i)=scr(4,i)
        u3dm(i)=scr(5,i)
        u4dm(i)=scr(6,i)
        masap(i)=scr(7,i)
        oripa(i)=scrint(1,i)
       end do

       deallocate(scr,scrint)

       if (n_st.gt.0) then
        
        if (ir_kern_stars.gt.n_esp-1) then
         write(*,*) 'warning: ir_kern_stars > n_esp-1',
     &              ir_kern_stars,n_esp-1
         write(*,*) 'fix n_esp to ir_kern_stars+1, at least'
         stop
        end if


        npart_esp(ir_kern_stars)=npart_esp(ir_kern_stars)+n_st
        write(*,*) 'stars: of species',ir_kern_stars,
     &             ', no. particles:',npart_esp(ir_kern_stars)
       end if

c       write(*,*) 'checking...'
*      check
c       bas=masap(1)
c       do i=2,n_dm
c        if (masap(i).gt.1.0001*bas) then
c         write(*,*) 'wrong, i=',i
c         stop
c        end if
c        bas=masap(i)
c       end do
c       stop

       return
       end
