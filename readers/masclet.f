*********************************************************************
       subroutine read_masclet(var,iter,nx,ny,nz,ndxyz,t,zeta,nl,npatch,
     &           pare,patchnx,patchny,patchnz,patchx,patchy,patchz,
     &           patchrx,patchry,patchrz,
     &           u2dm,u3dm,u4dm,masap,npart,rxpa,rypa,rzpa,oripa,n_dm)
*********************************************************************
*      reads masclet data: grids, gas density (clus files) and
*      dm particles information.
*      must be checked depending on the version/flavour of masclet
*      the simulation has been run with
*********************************************************************

       implicit none

       include 'input_files/asohf_parameters.dat'

       integer nx,ny,nz,iter,ndxyz,n_dm
       real*4 t,aaa,bbb,ccc,map,zeta

       integer i,j,k,ix,nl,ir,irr,var,n1,n2,n3

*      variables
       real*4 u1(nmax,nmay,nmaz)
       real*4 u11(namrx,namry,namrz,npalev)
       common /varia/ u1,u11  !!,u12,u13,u14

       integer npatch(0:nlevels)
       integer pare(npalev)
       integer patchnx(npalev)
       integer patchny(npalev)
       integer patchnz(npalev)
       integer patchx(npalev)
       integer patchy(npalev)
       integer patchz(npalev)
       real*4  patchrx(npalev)
       real*4  patchry(npalev)
       real*4  patchrz(npalev)

       integer npart(0:nlevels)
       real*4 u2dm(parti_read)
       real*4 u3dm(parti_read)
       real*4 u4dm(parti_read)
       real*4 masap(parti_read)
       real*4 rxpa(parti_read)
       real*4 rypa(parti_read)
       real*4 rzpa(parti_read)

       integer oripa(parti_read)

       real*4, allocatable::scr(:,:,:)

       real*4 ubas(0:parti_read)
       integer ubas2(0:parti_read),conta,low1,low2

       character*5 iter_string


*      reading data
       write(iter_string, '(i5.5)') iter !for saving files to disk
       write(*,*) 'reading iter',iter

       open (33,file='./simu_masclet/grids'//iter_string,
     &       status='unknown',action='read')
       open (31,file='./simu_masclet/clus'//iter_string,
     &       status='unknown',action='read',form='unformatted')
       open (32,file='./simu_masclet/cldm'//iter_string,
     &       status='unknown',action='read',form='unformatted')


*      grid data
       read(33,*) irr,t,nl,map
       read(33,*) zeta
       read(33,*) ir,ndxyz
       write(*,*) 'ir,nl,ndxyz,map', ir,nl,ndxyz,map

       do ir=1,nl
       read(33,*) irr,npatch(ir), npart(ir)
       write(*,*) 'npatch(ir), npart(ir)',npatch(ir), npart(ir)
       read(33,*)

       if (ir.ne.irr) write(*,*)'warning: fail in restart'
       low1=sum(npatch(0:ir-1))+1
       low2=sum(npatch(0:ir))
*       do i=1,npatch(ir)
       do i=low1,low2
        read(33,*) patchnx(i),patchny(i),patchnz(i)
        read(33,*) patchx(i),patchy(i),patchz(i)
        read(33,*) aaa,bbb,ccc
        patchrx(i)=aaa
        patchry(i)=bbb
        patchrz(i)=ccc
        read(33,*) pare(i)
       end do
       end do
       close(33)
       npart(0)=ndxyz

c       if (var.eq.1) then
*      baryonic
       read(31)
       ir=0
        n1=nx
        n2=ny
        n3=nz
        read(31) !(((u1g(i,j,k),i=1,n1),j=1,n2),k=1,n3)
        read(31) !u2
        read(31) !u3
        read(31) !u4
        read(31) !pres
        read(31) !pot
        read(31) !!opot
        read(31) !!caution with this line!! depends on masclet version: t
        read(31) !!new: metalicity!! depends on masclet version!: tracer
        read(31) !!abs(cr0amr(1:n1,1:n2,1:n3)-1)
        read(31) !bx
        read(31) !by
        read(31) !bz

       allocate(scr(namrx,namry,namrz))
       scr=0.0
       do ir=1,nl
       low1=sum(npatch(0:ir-1))+1
       low2=sum(npatch(0:ir))
       do i=low1,low2
        n1=patchnx(i)
        n2=patchny(i)
        n3=patchnz(i)
        read(31) !(((scr(ix,j,k),ix=1,n1),j=1,n2),k=1,n3)
           !u11g(1:n1,1:n2,1:n3,i)=scr(1:n1,1:n2,1:n3)
        read(31) !!(((scr(ix,j,k),ix=1,n1),j=1,n2),k=1,n3)
           !!u12(1:n1,1:n2,1:n3,i)=scr(1:n1,1:n2,1:n3)
        read(31) !!(((scr(ix,j,k),ix=1,n1),j=1,n2),k=1,n3)
           !!u13(1:n1,1:n2,1:n3,i)=scr(1:n1,1:n2,1:n3)
        read(31) !!(((scr(ix,j,k),ix=1,n1),j=1,n2),k=1,n3)
           !!u14(1:n1,1:n2,1:n3,i)=scr(1:n1,1:n2,1:n3)
        read(31) !pres21
        read(31) !pot1
        read(31) !opot
        read(31) !!caution with this line!! depends on masclet version: t
        read(31) !!new: metalicity!! depends on masclet version! tracer
        read(31) !!abs(cr0amr11(:,:,:,i)-1) (celdas=0 se eliminan del nivel por estar refinadas)
        read(31) !!abs(solapst(:,:,:,i)-1): se guarda solapst para saber que celdas estan solapadas
        read(31) !bx
        read(31) !by
        read(31) !bz


       end do
       end do
       deallocate(scr)

c      end if
      close(31)

***       if (var.eq.2) then
**     dark matter
       read(32)
       ir=0
        n1=nx
        n2=ny
        n3=nz

        read(32) (((u1(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        read(32) (rxpa(i),i=1,ndxyz)
        read(32) (rypa(i),i=1,ndxyz)
        read(32) (rzpa(i),i=1,ndxyz)
c        write(*,*)'hola2', maxval(rxpa(1:ndxyz)), maxval(rypa(1:ndxyz))

        read(32) (u2dm(i),i=1,ndxyz)
        read(32) (u3dm(i),i=1,ndxyz)
        read(32) (u4dm(i),i=1,ndxyz)
        read(32) (oripa(i),i=1,ndxyz)    !particle id
        conta=ndxyz
        masap(1:ndxyz)=map
        write(*,*) 'oripa=',maxval(oripa(1:ndxyz)),
     &                       minval(oripa(1:ndxyz))
        write(*,*) 'npart(0)=',ir, npart(ir),conta


       allocate(scr(namrx,namry,namrz))
       scr=0.0
       do ir=1,nl
       low1=sum(npatch(0:ir-1))+1
       low2=sum(npatch(0:ir))
       do i=low1,low2
        n1=patchnx(i)
        n2=patchny(i)
        n3=patchnz(i)

        read(32) (((scr(ix,j,k),ix=1,n1),j=1,n2),k=1,n3)
           u11(1:n1,1:n2,1:n3,i)=scr(1:n1,1:n2,1:n3)
       end do

        ubas=0.0
        ubas2=0
        read(32) (ubas(ix),ix=1,npart(ir))
        rxpa(conta+1:conta+npart(ir))=ubas(1:npart(ir))
        read(32) (ubas(ix),ix=1,npart(ir))
        rypa(conta+1:conta+npart(ir))=ubas(1:npart(ir))
        read(32) (ubas(ix),ix=1,npart(ir))
        rzpa(conta+1:conta+npart(ir))=ubas(1:npart(ir))
        read(32) (ubas(ix),ix=1,npart(ir))
        u2dm(conta+1:conta+npart(ir))=ubas(1:npart(ir))
        read(32) (ubas(ix),ix=1,npart(ir))
        u3dm(conta+1:conta+npart(ir))=ubas(1:npart(ir))
        read(32) (ubas(ix),ix=1,npart(ir))
        u4dm(conta+1:conta+npart(ir))=ubas(1:npart(ir))
        read(32) (ubas(ix),ix=1,npart(ir))
        masap(conta+1:conta+npart(ir))=ubas(1:npart(ir))

        read(32) (ubas2(ix),ix=1,npart(ir))
        if (npart(ir).gt.0)
     &   oripa(conta+1:conta+npart(ir))=ubas2(1:npart(ir))

        if (npart(ir).gt.0) then
        write(*,*) 'oripa=',maxval(oripa(conta+1:conta+npart(ir))),
     &                       minval(oripa(conta+1:conta+npart(ir)))
        end if

        conta=conta+npart(ir)
        write(*,*) 'npart(ir)=',ir,npart(ir),conta


       end do

       deallocate(scr)

       close(32)
***       end if

       write(*,*) 'total particles in iter=',conta
       if (conta.ne.n_dm) then
        write(*,*) 'warning: conta != n_dm',conta,n_dm
        stop
       end if

       return
       end

*********************************************************************
       subroutine read_particles_masclet(iter,nx,ny,nz,t,zeta,
     &             u2dm,u3dm,u4dm,masap,rxpa,rypa,rzpa,oripa,
     &             n_dm,var,n_st)
*********************************************************************
*      reads masclet data: grids and dm particles information.
*      must be checked depending on the version/flavour of masclet
*      the simulation has been run with
*********************************************************************

       implicit none

       include 'input_files/asohf_parameters.dat'

       integer nx,ny,nz,iter,ndxyz
       real*4 t,aaa,bbb,ccc,map,zeta
       integer var !(=1: only dm; =2: dm+stars)

       real bas
       integer i,j,k,ix,nl,ir,irr,n1,n2,n3,n_dm,n_st,nbas,are_bh,nst0
       integer consider_bh

       integer,allocatable::npatch(:)

       integer,allocatable::npart(:),npartst(:),npartbh(:)
       real*4 u2dm(parti_read)
       real*4 u3dm(parti_read)
       real*4 u4dm(parti_read)
       real*4 masap(parti_read)
       real*4 rxpa(parti_read)
       real*4 rypa(parti_read)
       real*4 rzpa(parti_read)

       integer oripa(parti_read)

       real*4 ubas(0:parti_read)
       integer ubas2(0:parti_read),conta,low1,low2

       integer max_oripa_dm
       common /oripastcorrect/ max_oripa_dm

       character*5 iter_string

       are_bh=1 ! depends on masclet version (are there bhs??)
       consider_bh=0

*      reading data
       write(iter_string, '(i5.5)') iter !for saving files to disk
       write(*,*) 'reading iter',iter

*      grid data
       open (33,file='./simu_masclet/grids'//iter_string,
     &       status='unknown',action='read')
       read(33,*) irr,t,nl,map
       read(33,*) zeta
       read(33,*) ir,ndxyz,nst0
       !write(*,*) 'ir,nl,ndxyz,map', ir,nl,ndxyz,map

       allocate(npatch(0:nl),npart(0:nl),npartst(0:nl),npartbh(0:nl))
       npatch=0
       npart=0
       npartst=0
       npartbh=0

       do ir=1,nl
       if (are_bh.eq.0) then
        read(33,*) irr,npatch(ir),npart(ir),npartst(ir)!,npartbh(ir)
       else
        read(33,*) irr,npatch(ir),npart(ir),npartst(ir),npartbh(ir)
       end if
       !write(*,*) 'npatch(ir), npart(ir)',npatch(ir), npart(ir)
       read(33,*)

       if (ir.ne.irr) write(*,*)'warning: fail in restart'
       low1=sum(npatch(0:ir-1))+1
       low2=sum(npatch(0:ir))
*       do i=1,npatch(ir)
       do i=low1,low2
        read(33,*) !patchnx(i),patchny(i),patchnz(i)
        read(33,*) !patchx(i),patchy(i),patchz(i)
        read(33,*) !patchrx(i),patchry(i),patchrz(i)
        read(33,*) !pare(i)
       end do
       end do
       npart(0)=ndxyz
       npartst(0)=nst0

       close(33)

**     dark matter
       open (32,file='./simu_masclet/cldm'//iter_string,
     &       status='unknown',action='read',form='unformatted')
       read(32)
       !ir=0
       read(32) !(((u1(i,j,k),i=1,nx),j=1,ny),k=1,nz)
       read(32) (rxpa(i),i=1,ndxyz)
       read(32) (rypa(i),i=1,ndxyz)
       read(32) (rzpa(i),i=1,ndxyz)
       read(32) (u2dm(i),i=1,ndxyz)
       read(32) (u3dm(i),i=1,ndxyz)
       read(32) (u4dm(i),i=1,ndxyz)
       read(32) (oripa(i),i=1,ndxyz)
       conta=ndxyz
       masap(1:ndxyz)=map

c       write(*,*) 'oripa=',maxval(oripa(1:ndxyz)),
c     &                      minval(oripa(1:ndxyz))
       write(*,*) 'npart(ir)=',0, npart(0),conta

       do ir=1,nl
        low1=sum(npatch(0:ir-1))+1
        low2=sum(npatch(0:ir))

        do i=low1,low2
         read(32) !(((scr(ix,j,k),ix=1,n1),j=1,n2),k=1,n3)
        end do

        ubas=0.0
        ubas2=0
        read(32) (ubas(ix),ix=1,npart(ir))
        rxpa(conta+1:conta+npart(ir))=ubas(1:npart(ir))
        read(32) (ubas(ix),ix=1,npart(ir))
        rypa(conta+1:conta+npart(ir))=ubas(1:npart(ir))
        read(32) (ubas(ix),ix=1,npart(ir))
        rzpa(conta+1:conta+npart(ir))=ubas(1:npart(ir))
        read(32) (ubas(ix),ix=1,npart(ir))
        u2dm(conta+1:conta+npart(ir))=ubas(1:npart(ir))
        read(32) (ubas(ix),ix=1,npart(ir))
        u3dm(conta+1:conta+npart(ir))=ubas(1:npart(ir))
        read(32) (ubas(ix),ix=1,npart(ir))
        u4dm(conta+1:conta+npart(ir))=ubas(1:npart(ir))
        read(32) (ubas(ix),ix=1,npart(ir))
        masap(conta+1:conta+npart(ir))=ubas(1:npart(ir))

        read(32) (ubas2(ix),ix=1,npart(ir))
        oripa(conta+1:conta+npart(ir))=ubas2(1:npart(ir))

        conta=conta+npart(ir)
        write(*,*) 'npart(ir)=',ir,npart(ir),conta
       end do

       close(32)

       write(*,*) 'total dm particles in iter=',conta
       n_dm=sum(npart(0:nl))
       if (conta.ne.n_dm) then
        write(*,*) 'warning: n_dm rewritten:',n_dm,'-->',conta
        n_dm=conta
       end if

*      fix oripas: particles of the heavier species get negative
       bas=maxval(masap(1:n_dm))
       if (bas.gt.4.0*minval(masap(1:n_dm))) then
        bas=0.9*bas
!$omp parallel do shared(n_dm,bas,oripa,masap),
!$omp+            private(i),
!$omp+            default(none)
        do i=1,n_dm
         if (masap(i).gt.bas) oripa(i)=-abs(oripa(i))
        end do
       end if
*      end fix oripas

       if (var.eq.2) then
        max_oripa_dm = maxval(oripa(1:n_dm))
        n_st=sum(npartst(0:nl))+sum(npartbh(0:nl))

        if (n_dm+n_st.gt.parti_read) then
         write(*,*) 'warning: bad dimensioning of parti_read',
     &               n_dm+n_st,'>',parti_read
         stop
        end if

        open (34,file='./simu_masclet/clst'//iter_string,
     &        status='unknown',action='read',form='unformatted')

        read(34) !iter,t4,zeta
        !ir=0
        ! note: in principle we assume no stellar particles will be
        !  located at l=0. this is an assumption in masclet in its
        !  present version.
        nbas=npartst(0)+npartbh(0)
        read(34) !(((u1st(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        read(34) (ubas(i),i=1,nbas)
        rxpa(conta+1:conta+nbas)=ubas(1:nbas)
        read(34) (ubas(i),i=1,nbas)
        rypa(conta+1:conta+nbas)=ubas(1:nbas)
        read(34) (ubas(i),i=1,nbas)
        rzpa(conta+1:conta+nbas)=ubas(1:nbas)
        read(34) (ubas(i),i=1,nbas)
        u2dm(conta+1:conta+nbas)=ubas(1:nbas)
        read(34) (ubas(i),i=1,nbas)
        u3dm(conta+1:conta+nbas)=ubas(1:nbas)
        read(34) (ubas(i),i=1,nbas)
        u4dm(conta+1:conta+nbas)=ubas(1:nbas)
        read(34) (ubas(i),i=1,nbas)
        masap(conta+1:conta+nbas)=ubas(1:nbas)
        read(34) !(ubas(i),i=1,ndxyz)
        read(34) !(ubas(i),i=1,ndxyz)
        conta=conta+nbas

        write(*,*) 'nst(ir)=',0,nbas,conta

        do ir=1,nl
         low1=sum(npatch(0:ir-1))+1
         low2=sum(npatch(0:ir))
         do i=low1,low2
          read(34) !(((scr(ix,j,k),ix=1,n1),j=1,n2),k=1,n3)
         end do

         nbas=npartst(ir)
         read(34) (ubas(ix),ix=1,nbas)
         rxpa(conta+1:conta+nbas)=ubas(1:nbas)
         read(34) (ubas(ix),ix=1,nbas)
         rypa(conta+1:conta+nbas)=ubas(1:nbas)
         read(34) (ubas(ix),ix=1,nbas)
         rzpa(conta+1:conta+nbas)=ubas(1:nbas)
         read(34) (ubas(ix),ix=1,nbas)
         u2dm(conta+1:conta+nbas)=ubas(1:nbas)
         read(34) (ubas(ix),ix=1,nbas)
         u3dm(conta+1:conta+nbas)=ubas(1:nbas)
         read(34) (ubas(ix),ix=1,nbas)
         u4dm(conta+1:conta+nbas)=ubas(1:nbas)
         read(34) (ubas(ix),ix=1,nbas)
         masap(conta+1:conta+nbas)=ubas(1:nbas)
         read(34)
         read(34)
         read(34) (ubas2(ix),ix=1,nbas)
         oripa(conta+1:conta+nbas)=max_oripa_dm+ubas2(1:nbas)

         conta=conta+nbas

         nbas=npartbh(ir)
         if (are_bh.eq.1) then
          if (consider_bh.eq.1) then
           read(34) (ubas(ix),ix=1,nbas)
           rxpa(conta+1:conta+nbas)=ubas(1:nbas)
           read(34) (ubas(ix),ix=1,nbas)
           rypa(conta+1:conta+nbas)=ubas(1:nbas)
           read(34) (ubas(ix),ix=1,nbas)
           rzpa(conta+1:conta+nbas)=ubas(1:nbas)
           read(34) (ubas(ix),ix=1,nbas)
           u2dm(conta+1:conta+nbas)=ubas(1:nbas)
           read(34) (ubas(ix),ix=1,nbas)
           u3dm(conta+1:conta+nbas)=ubas(1:nbas)
           read(34) (ubas(ix),ix=1,nbas)
           u4dm(conta+1:conta+nbas)=ubas(1:nbas)
           read(34) (ubas(ix),ix=1,nbas)
           masap(conta+1:conta+nbas)=ubas(1:nbas)
           read(34)
           read(34) (ubas2(ix),ix=1,nbas)
           oripa(conta+1:conta+nbas)=ubas2(1:nbas)
           conta=conta+nbas
          else
           read(34) !(ubas(ix),ix=1,nbas)
           !rxpa(conta+1:conta+nbas)=ubas(1:nbas)
           read(34) !(ubas(ix),ix=1,nbas)
           !rypa(conta+1:conta+nbas)=ubas(1:nbas)
           read(34) !(ubas(ix),ix=1,nbas)
           !rzpa(conta+1:conta+nbas)=ubas(1:nbas)
           read(34) !(ubas(ix),ix=1,nbas)
           !u2dm(conta+1:conta+nbas)=ubas(1:nbas)
           read(34) !(ubas(ix),ix=1,nbas)
           !u3dm(conta+1:conta+nbas)=ubas(1:nbas)
           read(34) !(ubas(ix),ix=1,nbas)
           !u4dm(conta+1:conta+nbas)=ubas(1:nbas)
           read(34) !(ubas(ix),ix=1,nbas)
           !masap(conta+1:conta+nbas)=ubas(1:nbas)
           read(34) !
           read(34) !(ubas2(ix),ix=1,nbas)
           !oripa(conta+1:conta+nbas)=ubas2(1:nbas)
           !conta=conta+nbas
          end if
         end if

         write(*,*) 'nst(ir)=',ir,npartst(ir)+npartbh(ir),conta
        end do

        close(34)
       else
        n_st=0
       end if !(var.eq.2)

       deallocate(npatch,npart,npartst,npartbh)

       return
       end