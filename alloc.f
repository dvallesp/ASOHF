************************************************************************
      SUBROUTINE INIT_OUTVARS(MASA,RADIO,CLUSRX,CLUSRY,CLUSRZ,CLUSRXCM,
     &           CLUSRYCM,CLUSRZCM,MSUB,RSUB,PATCHCLUS,REALCLUS,
     &           HALBORDERS,LEVHAL,CONCENTRA,ANGULARM,VMAXCLUS,IPLIP,VX,
     &           VY,VZ,MEAN_VR,VCMAX,MCMAX,RCMAX,M200C,M500C,M2500C,
     &           R200C,R500C,R2500C,M200M,M500M,M2500M,R200M,R500M,
     &           R2500M,PROFILES,VELOCITY_DISPERSION,RMAXSIGMA,
     &           MMAXSIGMA,KINETIC_E,POTENTIAL_E,FSUB,NSUBS,
     &           INDCS_PARTICLES_PER_HALO,DMPCLUS,NHALLEV,SUBS_LEV,
     &           SUBHALOS,EIGENVAL,INERTIA_TENSOR,NCLUS)
************************************************************************

       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

*      ---HALOS AND SUBHALOS---
       REAL*4 MASA(MAXNCLUS), RADIO(MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
       REAL*4 CLUSRXCM(MAXNCLUS),CLUSRYCM(MAXNCLUS),CLUSRZCM(MAXNCLUS)
       REAL*4 MSUB(MAXNCLUS),RSUB(MAXNCLUS)
       INTEGER PATCHCLUS(MAXNCLUS),REALCLUS(MAXNCLUS)
       INTEGER HALBORDERS(MAXNCLUS),LEVHAL(MAXNCLUS)
       REAL*4 CONCENTRA(NMAXNCLUS)
       REAL*4 ANGULARM(3,NMAXNCLUS)
       REAL*4 VMAXCLUS(NMAXNCLUS)
       INTEGER IPLIP(NMAXNCLUS)
       REAL*4 VX(NMAXNCLUS),VY(NMAXNCLUS),VZ(NMAXNCLUS)
       REAL*4 MEAN_VR(NMAXNCLUS)
       REAL*4 VCMAX(NMAXNCLUS),MCMAX(NMAXNCLUS),RCMAX(NMAXNCLUS)
       REAL*4 M200C(NMAXNCLUS),R200C(NMAXNCLUS)
       REAL*4 M500C(NMAXNCLUS),R500C(NMAXNCLUS)
       REAL*4 M2500C(NMAXNCLUS),R2500C(NMAXNCLUS)
       REAL*4 M200M(NMAXNCLUS),R200M(NMAXNCLUS)
       REAL*4 M500M(NMAXNCLUS),R500M(NMAXNCLUS)
       REAL*4 M2500M(NMAXNCLUS),R2500M(NMAXNCLUS)
       REAL*4 PROFILES(NBINS,2,NMAXNCLUS)
       REAL*4 VELOCITY_DISPERSION(NMAXNCLUS)
       REAL*4 RMAXSIGMA(NMAXNCLUS),MMAXSIGMA(NMAXNCLUS)
       REAL*4 KINETIC_E(NMAXNCLUS),POTENTIAL_E(NMAXNCLUS)
       REAL*4 FSUB(NMAXNCLUS) ! fraction of mass in substructures
       INTEGER NSUBS(NMAXNCLUS) ! number of substructures
       INTEGER INDCS_PARTICLES_PER_HALO(2,NMAXNCLUS)
       INTEGER DMPCLUS(MAXNCLUS)
       INTEGER SUBHALOS(NMAXNCLUS)
       REAL*4 EIGENVAL(3,NMAXNCLUS)
       REAL*4 INERTIA_TENSOR(6,NMAXNCLUS)
       INTEGER NHALLEV(0:NLEVELS),SUBS_LEV(0:NLEVELS)
       INTEGER NCLUS

       INTEGER NMAXNCLUSBAS,I

       NMAXNCLUSBAS=MAXNCLUS
!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,CLUSRX,CLUSRY,CLUSRZ,CLUSRXCM,
!$OMP+                   CLUSRYCM,CLUSRZCM,MASA,RADIO,MSUB,RSUB,
!$OMP+                   PATCHCLUS,REALCLUS,HALBORDERS,LEVHAL),
!$OMP+            PRIVATE(I)
       DO I=1,NMAXNCLUSBAS
        CLUSRX(I)=0.0
        CLUSRY(I)=0.0
        CLUSRZ(I)=0.0
        CLUSRXCM(I)=0.0
        CLUSRYCM(I)=0.0
        CLUSRZCM(I)=0.0
        MASA(I)=0.0
        RADIO(I)=0.0
        MSUB(I)=0.0
        RSUB(I)=0.0
        PATCHCLUS(I)=0
        REALCLUS(I)=0
        HALBORDERS(I)=0
        LEVHAL(I)=0
       END DO

       NMAXNCLUSBAS=NMAXNCLUS
!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,VCMAX,MCMAX,RCMAX,DMPCLUS,M200C,
!$OMP+                   M500C,M2500C,M200M,M500M,M2500M,R200C,R500C,
!$OMP+                   R2500C,R200M,R500M,R2500M,IPLIP,REALCLUS,
!$OMP+                   LEVHAL,EIGENVAL,INERTIA_TENSOR,MEAN_VR,
!$OMP+                   VELOCITY_DISPERSION,RMAXSIGMA,MMAXSIGMA,
!$OMP+                   KINETIC_E,POTENTIAL_E,FSUB,NSUBS,
!$OMP+                   INDCS_PARTICLES_PER_HALO,CONCENTRA,ANGULARM,
!$OMP+                   VX,VY,VZ,PROFILES,VMAXCLUS,SUBHALOS),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
       DO I=1,NMAXNCLUSBAS
        VCMAX(I)=0.0
        MCMAX(I)=0.0
        RCMAX(I)=0.0
        M200C(I)=0.0
        M500C(I)=0.0
        M2500C(I)=0.0
        M200M(I)=0.0
        M500M(I)=0.0
        M2500M(I)=0.0
        R200C(I)=0.0
        R500C(I)=0.0
        R2500C(I)=0.0
        R200M(I)=0.0
        R500M(I)=0.0
        R2500M(I)=0.0
        IPLIP(I)=0
        LEVHAL(I)=0
        DMPCLUS(I)=0
        REALCLUS(I)=0    !de momento no hay halos
        EIGENVAL(:,I)=0.0
        INERTIA_TENSOR(:,I)=0.0
        MEAN_VR(I)=0.0
        VELOCITY_DISPERSION(I)=0.0
        RMAXSIGMA(I)=0.0
        MMAXSIGMA(I)=0.0
        KINETIC_E(I)=0.0
        POTENTIAL_E(I)=0.0
        FSUB(I)=0.0
        NSUBS(I)=0
        INDCS_PARTICLES_PER_HALO(:,I)=0
        CONCENTRA(I)=0.0
        ANGULARM(:,I)=0.0
        VMAXCLUS(I)=0.0
        VX(I)=0.0
        VY(I)=0.0
        VZ(I)=0.0
        PROFILES(:,:,I)=0
        SUBHALOS(I)=0
       END DO

       NHALLEV=0
       SUBS_LEV=0
       NCLUS=0

       RETURN
      END
************************************************************************
      SUBROUTINE INIT_GRIDVARS(NX,NY,NZ,PATCHNX,PATCHNY,PATCHNZ,PATCHX,
     &           PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,PARE,NPATCH,U1,
     &           U11,ROTE,RETE)
************************************************************************
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NX,NY,NZ
       INTEGER NPATCH(0:NLEVELS),PARE(NPALEV)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
       REAL*4  PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
       REAL*4 U1(NMAX,NMAY,NMAZ),U11(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL ROTE,RETE

       INTEGER I,J,K,IX,JY,KZ,N1,NPBAS

       PATCHNX=0
       PATCHNY=0
       PATCHNZ=0
       PATCHX=0
       PATCHY=0
       PATCHZ=0
       PATCHRX=0.0
       PATCHRY=0.0
       PATCHRZ=0.0
       NPATCH=0
       PARE=0

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U1),PRIVATE(I,J,K)
       DO K=1,NZ
       DO J=1,NY
       DO I=1,NX
        U1(I,J,K)=-1.0
       END DO
       END DO
       END DO

**********
       N1=NAMRX
       NPBAS=NPALEV
**********

!$OMP PARALLEL DO SHARED(NPBAS,N1,U11),
!$OMP+        PRIVATE(IX,JY,KZ,I)
       DO I=1,NPBAS
        DO KZ=1,N1
        DO JY=1,N1
        DO IX=1,N1
         U11(IX,JY,KZ,I)=-1.0
        END DO
        END DO
        END DO
       END DO


       ROTE=0.0
       RETE=0.0

       RETURN
      END

************************************************************************
      SUBROUTINE READ_AND_ALLOC_PARTICLES(ITER,NX,
     &           NY,NZ,T,ZETA,N_DM,VAR,N_ST,N_PARTICLES,UV,UM,
     &           HUBBLE_LITTLEH,LADO0,LADO)
************************************************************************
       USE PARTICLES
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER ITER,NX,NY,NZ
       REAL T,ZETA
       INTEGER N_DM,VAR,N_ST,N_PARTICLES
       REAL UV,UM,HUBBLE_LITTLEH,LADO0,LADO

       REAL*4 U2DM_R(PARTI_READ),U3DM_R(PARTI_READ),U4DM_R(PARTI_READ)
       REAL*4 MASAP_R(PARTI_READ)
       REAL*4 RXPA_R(PARTI_READ),RYPA_R(PARTI_READ),RZPA_R(PARTI_READ)
       INTEGER ORIPA_R(PARTI_READ),KEEP(PARTI_READ)
       INTEGER I

*      Domain decomposition
       REAL CIO_MASS,CIO_SPEED,CIO_LENGTH,CIO_ALPHA,CIO_XC,CIO_YC,CIO_ZC
       COMMON /CONV_IO/ CIO_MASS,CIO_SPEED,CIO_LENGTH,CIO_ALPHA,CIO_XC,
     &                  CIO_YC,CIO_ZC
       INTEGER DO_DOMDECOMP
       REAL DDXL,DDXR,DDYL,DDYR,DDZL,DDZR
       COMMON /DOM_DECOMP/ DO_DOMDECOMP,DDXL,DDXR,DDYL,DDYR,DDZL,DDZR
       REAL X1,X2,Y1,Y2,Z1,Z2,FACT_LENGTH,BASX,BASY,BASZ
       real basxc,basyc,baszc,x3,x4,y3,y4,z3,z4
       INTEGER PARTIST,II

       REAL*4 DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ
       REAL*4  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
       COMMON /GRID/   RADX,RADY,RADZ

       REAL DDXC,DDYC,DDZC,SHIFT_X,SHIFT_Y,SHIFT_Z
       real shift_x0,shift_y0,shift_z0

       integer period_x,period_y,period_z 
       common /period_correct/ period_x,period_y,period_z

!$OMP PARALLEL DO SHARED(RXPA_R,RYPA_R,RZPA_R,U2DM_R,U3DM_R,U4DM_R,
!$OMP+                   MASAP_R,ORIPA_R),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
       DO I=1,PARTI_READ
        RXPA_R(I)=0.0
        RYPA_R(I)=0.0
        RZPA_R(I)=0.0
        U2DM_R(I)=0.0
        U3DM_R(I)=0.0
        U4DM_R(I)=0.0
        MASAP_R(I)=0.0
        ORIPA_R(I)=0
       END DO

#if reader == 1
        CALL READ_PARTICLES_MASCLET(ITER,NX,NY,NZ,T,ZETA,
     &                              U2DM_R,U3DM_R,U4DM_R,MASAP_R,
     &                              RXPA_R,RYPA_R,RZPA_R,ORIPA_R,
     &                              N_DM,VAR,N_ST)
#elif reader == 0
         CALL READ_PARTICLES_GENERAL(ITER,NX,NY,NZ,T,ZETA,
     &                               U2DM_R,U3DM_R,U4DM_R,MASAP_R,
     &                               RXPA_R,RYPA_R,RZPA_R,ORIPA_R,
     &                               N_DM,VAR,N_ST,
     &                               UV,UM,HUBBLE_LITTLEH)
#elif reader == 2
         CALL READ_GADGET_UNFORMATTED(ITER,NX,NY,NZ,T,ZETA,
     &                                U2DM_R,U3DM_R,U4DM_R,MASAP_R,
     &                                RXPA_R,RYPA_R,RZPA_R,ORIPA_R,
     &                                N_DM,VAR,N_ST,
     &                                UV,UM,HUBBLE_LITTLEH)
#elif reader == 3
         CALL READ_GIZMO_HDF5(ITER,NX,NY,NZ,T,ZETA,
     &                        U2DM_R,U3DM_R,U4DM_R,MASAP_R,
     &                        RXPA_R,RYPA_R,RZPA_R,ORIPA_R,
     &                        N_DM,VAR,N_ST,
     &                        UV,UM,HUBBLE_LITTLEH)
#else
         WRITE(*,*) 'ERROR: reader not defined'
         STOP 1
#endif

       N_PARTICLES=N_DM+N_ST
       WRITE(*,*) 'DM, stars, total particles:',N_DM,N_ST,N_PARTICLES

****************************
*      DOMAIN DECOMPOSITION
****************************

       X1=(DDXL-CIO_XC)*CIO_LENGTH
       X2=(DDXR-CIO_XC)*CIO_LENGTH
       Y1=(DDYL-CIO_YC)*CIO_LENGTH
       Y2=(DDYR-CIO_YC)*CIO_LENGTH
       Z1=(DDZL-CIO_ZC)*CIO_LENGTH
       Z2=(DDZR-CIO_ZC)*CIO_LENGTH

       if (do_domdecomp.eq.2) then ! we automatically select the region,
              ! based on the particles read, with some extra buffer 
        ! For periodicity correction
        shift_x0 = 0. 
        shift_y0 = 0. 
        shift_z0 = 0.
        period_x = 0
        period_y = 0
        period_z = 0

        x1 = minval(rxpa_r(1:n_dm+n_st))
        x2 = maxval(rxpa_r(1:n_dm+n_st))
        y1 = minval(rypa_r(1:n_dm+n_st))
        y2 = maxval(rypa_r(1:n_dm+n_st))
        z1 = minval(rzpa_r(1:n_dm+n_st))
        z2 = maxval(rzpa_r(1:n_dm+n_st))

        ! Detect if the halo goes through the periodic boundary
        if (x2-x1 > 0.9*lado0) then ! particles span almost the full box
         x3 = x1 + 0.25*(x2-x1)
         x4 = x2 - 0.25*(x2-x1)
         if (count(rxpa_r(1:n_dm+n_st).ge.x3.and.
     &             rxpa_r(1:n_dm+n_st).le.x4).eq.0) then ! there's a gap between the two ends
          write(*,*) 'Warning: correcting for periodic boundary in x'
!$omp parallel do shared(rxpa_r,n_dm,n_st,lado0), 
!$omp+            private(i), default(none)
          do i=1,n_dm+n_st
           if (rxpa_r(i).lt.0.) then 
            rxpa_r(i)=rxpa_r(i)+lado0/2.
           else
            rxpa_r(i)=rxpa_r(i)-lado0/2.
           end if
          end do
          shift_x0 = lado0/2.

          x1 = minval(rxpa_r(1:n_dm+n_st))
          x2 = maxval(rxpa_r(1:n_dm+n_st))
          period_x = 1
         end if ! count(...)
        end if !(x2-x1 > 0.9*lado0)

        if (y2-y1 > 0.9*lado0) then ! particles span almost the full box
         y3 = y1 + 0.25*(y2-y1)
         y4 = y2 - 0.25*(y2-y1)
         if (count(rypa_r(1:n_dm+n_st).ge.y3.and.
     &             rypa_r(1:n_dm+n_st).le.y4).eq.0) then ! there's a gap between the two ends
          write(*,*) 'Warning: correcting for periodic boundary in y'
!$omp parallel do shared(rypa_r,n_dm,n_st,lado0), 
!$omp+            private(i), default(none)
          do i=1,n_dm+n_st
           if (rypa_r(i).lt.0.) then 
            rypa_r(i)=rypa_r(i)+lado0/2.
           else
            rypa_r(i)=rypa_r(i)-lado0/2.
           end if
          end do
          shift_y0 = lado0/2.

          y1 = minval(rypa_r(1:n_dm+n_st))
          y2 = maxval(rypa_r(1:n_dm+n_st))
          period_y = 1
         end if ! count(...)
        end if !(y2-y1 > 0.9*lado0)

        if (z2-z1 > 0.9*lado0) then ! particles span almost the full box
         z3 = z1 + 0.25*(z2-z1)
         z4 = z2 - 0.25*(z2-z1)
         if (count(rzpa_r(1:n_dm+n_st).ge.z3.and.
     &             rzpa_r(1:n_dm+n_st).le.z4).eq.0) then ! there's a gap between the two ends
          write(*,*) 'Warning: correcting for periodic boundary in z'
!$omp parallel do shared(rzpa_r,n_dm,n_st,lado0), 
!$omp+            private(i), default(none)
          do i=1,n_dm+n_st
           if (rzpa_r(i).lt.0.) then 
            rzpa_r(i)=rzpa_r(i)+lado0/2.
           else
            rzpa_r(i)=rzpa_r(i)-lado0/2.
           end if
          end do
          shift_z0 = lado0/2.

          z1 = minval(rzpa_r(1:n_dm+n_st))
          z2 = maxval(rzpa_r(1:n_dm+n_st))
          period_z = 1
         end if ! count(...)
        end if !(z2-z1 > 0.9*lado0)

        basx = (x2 - x1) * 1.1 / 2.
        basy = (y2 - y1) * 1.1 / 2.
        basz = (z2 - z1) * 1.1 / 2.
        basxc = (x2 + x1) / 2.
        basyc = (y2 + y1) / 2.
        baszc = (z2 + z1) / 2.

        x1 = basxc - basx
        x2 = basxc + basx
        y1 = basyc - basy
        y2 = basyc + basy
        z1 = baszc - basz
        z2 = baszc + basz

        write(*,*) 'Domain decomposition automatically selected!'
       end if ! (do_domdecomp.eq.2)

!$OMP PARALLEL DO SHARED(KEEP,N_PARTICLES), PRIVATE(I), DEFAULT(NONE)
       DO I=1,N_PARTICLES
        KEEP(I)=1
       END DO

       IF (DO_DOMDECOMP.EQ.0) THEN
        PARTI=N_PARTICLES
        PARTIST=N_ST
       ELSE ! options 1 and 2
        PARTI=0
        PARTIST=0
!$OMP PARALLEL DO SHARED(N_PARTICLES,RXPA_R,RYPA_R,RZPA_R,X1,Y1,Z1,X2,
!$OMP+                   Y2,Z2,KEEP,N_DM),
!$OMP+            PRIVATE(I,BASX,BASY,BASZ),
!$OMP+            REDUCTION(+:PARTI,PARTIST)
!$OMP+            DEFAULT(NONE)
        DO I=1,N_PARTICLES
         BASX=(RXPA_R(I)-X1)*(X2-RXPA_R(I))
         IF (BASX.LT.0.0) THEN
          KEEP(I)=0
          CYCLE
         END IF

         BASY=(RYPA_R(I)-Y1)*(Y2-RYPA_R(I))
         IF (BASY.LT.0.0) THEN
          KEEP(I)=0
          CYCLE
         END IF

         BASZ=(RZPA_R(I)-Z1)*(Z2-RZPA_R(I))
         IF (BASZ.LT.0.0) THEN
          KEEP(I)=0
          CYCLE
         END IF

         PARTI=PARTI+1
         IF (I.GT.N_DM) PARTIST=PARTIST+1
        END DO

        WRITE(*,*) 'Performing domain decomposition in'
        WRITE(*,*) '... (input):',DDXL,DDXR,DDYL,DDYR,DDZL,DDZR
        WRITE(*,*) '... (internal)',X1,X2,Y1,Y2,Z1,Z2
        WRITE(*,*) 'Particles (read --> keep)',N_PARTICLES,'-->',PARTI
        WRITE(*,*) '... Of which DM',N_DM,'-->',PARTI-PARTIST
        WRITE(*,*) '... Of which stellar',N_ST,'-->',PARTIST
       END IF ! (DO_DOMDECOMP.EQ.0) THEN; else

****   ALLOCATE PARTICLES FROM PARTICLES MODULE
       ALLOCATE(RXPA(PARTI), RYPA(PARTI), RZPA(PARTI))
       ALLOCATE(U2DM(PARTI), U3DM(PARTI), U4DM(PARTI))
       ALLOCATE(MASAP(PARTI))
       ALLOCATE(ORIPA(PARTI))
       PARTI_PPH = PARTI
       ALLOCATE(PARTICLES_PER_HALO(PARTI_PPH))
***********************************************

       II=0
       DO I=1,N_PARTICLES
        IF (KEEP(I).EQ.0) CYCLE
        II=II+1
        RXPA(II)=RXPA_R(I)
        RYPA(II)=RYPA_R(I)
        RZPA(II)=RZPA_R(I)
        U2DM(II)=U2DM_R(I)
        U3DM(II)=U3DM_R(I)
        U4DM(II)=U4DM_R(I)
        MASAP(II)=MASAP_R(I)
        ORIPA(II)=ORIPA_R(I)
        PARTICLES_PER_HALO(II)=0
       END DO

       IF (II.NE.PARTI) THEN
        WRITE(*,*) 'ERROR in domain decompose, II.NE.PARTI',II,PARTI
        STOP 1
       END IF

       N_DM=PARTI-PARTIST
       N_ST=PARTIST
       N_PARTICLES=PARTI

       IF (do_domdecomp.ne.0) THEN !Rebuild the base grid, recenter stuff
        LADO0=MAX(X2-X1,Y2-Y1,Z2-Z1)
        LADO=LADO0-(LADO0/NX)
        CALL MALLA(NX,NY,NZ,LADO)

        WRITE(*,*)
        WRITE(*,*) '************************************************'
        WRITE(*,*) '                   NEW GRID                     '
        WRITE(*,*) '************************************************'
        WRITE(*,*) 'SIDE LENGTH=',LADO
        WRITE(*,*) 'NX,DX,RADX(1),RADX(NX)=',NX,DX,RADX(1),RADX(NX)
        WRITE(*,*)

        DDXC=(X1+X2)/2.
        DDYC=(Y1+Y2)/2.
        DDZC=(Z1+Z2)/2.
        SHIFT_X=DDXC
        SHIFT_Y=DDYC
        SHIFT_Z=DDZC

!$OMP PARALLEL DO SHARED(RXPA,RYPA,RZPA,N_PARTICLES,SHIFT_X,SHIFT_Y,
!$OMP+                   SHIFT_Z,shift_x0,shift_y0,shift_z0),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
        DO I=1,N_PARTICLES
         RXPA(I)=RXPA(I)-SHIFT_X
         RYPA(I)=RYPA(I)-SHIFT_Y
         RZPA(I)=RZPA(I)-SHIFT_Z
        END DO

        CIO_XC=CIO_XC+SHIFT_X/CIO_LENGTH 
        if (period_x.eq.1) cio_xc = cio_xc + shift_x0/cio_length
        CIO_YC=CIO_YC+SHIFT_Y/CIO_LENGTH 
        if (period_y.eq.1) cio_yc = cio_yc + shift_y0/cio_length
        CIO_ZC=CIO_ZC+SHIFT_Z/CIO_LENGTH 
        if (period_z.eq.1) cio_zc = cio_zc + shift_z0/cio_length

        WRITE(*,*)
        WRITE(*,*) 'After domain decomposition...'
        WRITE(*,*) 'x positions (min,max), in Mpc:',
     &     MINVAL(RXPA(1:N_DM+N_ST)),MAXVAL(RXPA(1:N_DM+N_ST))
        WRITE(*,*) 'y positions (min,max), in Mpc:',
     &     MINVAL(RYPA(1:N_DM+N_ST)),MAXVAL(RYPA(1:N_DM+N_ST))
        WRITE(*,*) 'z positions (min,max), in Mpc:',
     &     MINVAL(RZPA(1:N_DM+N_ST)),MAXVAL(RZPA(1:N_DM+N_ST))
        WRITE(*,*)

       END IF !(DO_DOMDECOMP.EQ.1)

       RETURN
      END

************************************************************************
      SUBROUTINE DEALLOC
************************************************************************
       USE PARTICLES
       IMPLICIT NONE

       PARTI=0
       DEALLOCATE(RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,ORIPA)
       DEALLOCATE(PARTICLES_PER_HALO)

       RETURN
      END
