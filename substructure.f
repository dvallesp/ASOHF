********************************************************************
       SUBROUTINE SEARCH_SUBSTRUCTURE_GRID(IR,NL,NX,NY,NZ,NPATCH,
     &             PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &             PATCHRY,PATCHRZ,PARE,NCLUS,MASA,RADIO,CLUSRX,CLUSRY,
     &             CLUSRZ,REALCLUS,LEVHAL,NHALLEV,BOUND,CONTRASTEC,RODO,
     &             SOLAP,VECINO,NVECI,CR0AMR,CR0AMR11,PATCHCLUS,
     &             VOL_SOLAP_LOW,CLUSRXCM,CLUSRYCM,CLUSRZCM,RSUB,MSUB,
     &             SUBS_LEV,UM,PROFILES)
********************************************************************
*      Pipeline for tentative halo finding over the grid
********************************************************************

       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

*      I/O DATA
       INTEGER NL,NX,NY,NZ
       INTEGER NPATCH(0:NLEVELS)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
       REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
       INTEGER PARE(NPALEV)
       INTEGER NCLUS
       REAL MASA(MAXNCLUS),RADIO(MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
       REAL*4 CLUSRXCM(MAXNCLUS),CLUSRYCM(MAXNCLUS),CLUSRZCM(MAXNCLUS)
       INTEGER REALCLUS(MAXNCLUS),LEVHAL(MAXNCLUS)
       INTEGER PATCHCLUS(MAXNCLUS)
       INTEGER NHALLEV(0:NLEVELS)
       REAL BOUND,CONTRASTEC,RODO
       INTEGER SOLAP(NAMRX,NAMRY,NAMRZ,NPALEV)
       INTEGER VECINO(NPALEV,NPALEV),NVECI(NPALEV)
       INTEGER CR0AMR(NMAX,NMAY,NMAZ)
       INTEGER CR0AMR11(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL VOL_SOLAP_LOW
       REAL RSUB(MAXNCLUS),MSUB(MAXNCLUS)
       INTEGER SUBS_LEV(0:NLEVELS)
       REAL UM
       REAL*4 PROFILES(NBINS,2,NMAXNCLUS)

*      GLOBAL VARIABLES
       REAL*4 DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       REAL*4  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
       COMMON /GRID/ RADX,RADY,RADZ

       REAL*4  RX(0:NAMRX+1,NPALEV),RY(0:NAMRX+1,NPALEV),
     &         RZ(0:NAMRX+1,NPALEV)
       COMMON /GRIDAMR/ RX,RY,RZ

       REAL*4 RETE,HTE,ROTE
       COMMON /BACK/ RETE,HTE,ROTE

       REAL*4 U1(NMAX,NMAY,NMAZ)
       REAL*4 U1G(NMAX,NMAY,NMAZ)
       REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL*4 U11G(NAMRX,NAMRY,NAMRZ,NPALEV)
       COMMON /VARIA/ U1,U11,U1G,U11G

       REAL*4 ACHE,T0,RE0
       COMMON /DOS/ ACHE,T0,RE0

*      LOCAL VARIABLES
       INTEGER CONTA(NMAX,NMAY,NMAZ)
       INTEGER CONTA1(NAMRX,NAMRY,NAMRZ,NPALEV)

       REAL MAXIMO(NPALEV)
       INTEGER VID(NLEVELS,NPALEV),NVID(NLEVELS)
       INTEGER RELEVANT_PATCHES(NPALEV),NRELEVANT_PATCHES(NLEVELS)

       INTEGER IR,IX,JY,KZ,I,J,K,II,JJ,KK,IPATCH,ICEN(3),NV_GOOD
       INTEGER L1,L2,L3,NX1,NX2,NY1,NY2,NZ1,NZ2,KK_ENTERO,ITER_GROW
       INTEGER N1,N2,N3,KONTA,LOW1,LOW2,I2,ICEN1(1),ICEN4(4),CEL,I1,J1
       INTEGER K1,BORAMR,IRR,BASINT,FLAG_ITER,BOR,FLAG,IHOSTHALO
       INTEGER FLAG_VIR,FLAG_JACOBI,NCLUS_INI,FLAG_DUP
       REAL PRUEBAX,PRUEBAY,PRUEBAZ,RMIN,BASMASS_SHELL,BASMASS,DELTA
       REAL ESP,ESP_LOG,BAS,KK_REAL,RSHELL,R_INT,R_EXT,RANT,LADO0
       REAL BASDELTA,AA,PI,VOLCELL,BASX,BASY,BASZ,BASVOL,DISTA,MINDISTA
       REAL X1,X2,Y1,Y2,Z1,Z2,DXPA,DYPA,DZPA,BASXX,BASYY,BASZZ
       REAL XCEN,YCEN,ZCEN,BOUNDIR,X3,Y3,Z3,X4,Y4,Z4,MINDERIV
       REAL VECDENS(1000),VECRAD(1000),DERIVATIVE(1000),BASVOL_SHELL
       REAL MHOST,DISTHOST,EQ_JACOBI_R,DXPAPA,DYPAPA,DZPAPA
       REAL CONTRASTECPEAK

       REAL*4, ALLOCATABLE::DDD(:)
       INTEGER, ALLOCATABLE::DDDX(:),DDDY(:),DDDZ(:),DDDP(:)

       CONTRASTECPEAK=CONTRASTEC/6.0

       NCLUS_INI=NCLUS
       WRITE(*,*) 'Looking for substructure at level', IR
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)

*      Find the cells covered by preexisting haloes
!$OMP PARALLEL DO SHARED(NPATCH,PATCHNX,PATCHNY,PATCHNZ,CONTA1),
!$OMP+            PRIVATE(I,N1,N2,N3,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        DO KZ=1,N3
        DO JY=1,N2
        DO IX=1,N1
         CONTA1(IX,JY,KZ,I)=0
        END DO
        END DO
        END DO
       END DO

!$OMP PARALLEL DO SHARED(NPATCH,PATCHNX,PATCHNY,PATCHNZ,PATCHRX,PATCHRY,
!$OMP+                   PATCHRZ,DXPA,DYPA,DZPA,CLUSRX,CLUSRY,CLUSRZ,
!$OMP+                   RADIO,RX,RY,RZ,CONTA1,NCLUS,REALCLUS,RSUB),
!$OMP+            PRIVATE(I,N1,N2,N3,X1,X2,Y1,Y2,Z1,Z2,BASX,BASY,BASZ,
!$OMP+                    BAS,IX,JY,KZ,AA,X3,Y3,Z3,X4,Y4,Z4),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        X1=PATCHRX(I)-DXPA
        Y1=PATCHRY(I)-DYPA
        Z1=PATCHRZ(I)-DZPA
        X2=X1+N1*DXPA
        Y2=Y1+N2*DYPA
        Z2=Z1+N3*DZPA
        DO II=1,NCLUS
         BASX=CLUSRX(II)
         BASY=CLUSRY(II)
         BASZ=CLUSRZ(II)
         IF (REALCLUS(II).EQ.0) THEN
          CYCLE
         ELSE IF (REALCLUS(II).EQ.-1) THEN
          BAS=RADIO(II)
         ELSE
          BAS=RSUB(II)
         END IF
         X3=BASX-BAS
         Y3=BASY-BAS
         Z3=BASZ-BAS
         X4=BASX+BAS
         Y4=BASY+BAS
         Z4=BASZ+BAS
         IF (X1.LE.X4.AND.X3.LE.X2.AND.
     &       Y1.LE.Y4.AND.Y3.LE.Y2.AND.
     &       Z1.LE.Z4.AND.Z3.LE.Z2) THEN
c           WRITE(*,*) BASX,BASY,BASZ,BAS
          DO KZ=1,N3
          DO JY=1,N2
          DO IX=1,N1
           AA=SQRT((RX(IX,I)-BASX)**2 +
     &             (RY(JY,I)-BASY)**2 +
     &             (RZ(KZ,I)-BASZ)**2)
           IF (AA.LE.BAS) CONTA1(IX,JY,KZ,I)=1
          END DO
          END DO
          END DO
         END IF
        END DO
       END DO

!      Unflag centers of haloes to avoid double-counting
!$OMP PARALLEL DO SHARED(NPATCH,PATCHNX,PATCHNY,PATCHNZ,PATCHRX,PATCHRY,
!$OMP+                   PATCHRZ,DXPA,DYPA,DZPA,CLUSRX,CLUSRY,CLUSRZ,
!$OMP+                   RADIO,RX,RY,RZ,CONTA1,NCLUS,REALCLUS),
!$OMP+            PRIVATE(I,N1,N2,N3,X1,X2,Y1,Y2,Z1,Z2,BASX,BASY,BASZ,
!$OMP+                    BAS,IX,JY,KZ,AA,X3,Y3,Z3,X4,Y4,Z4),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        X1=PATCHRX(I)-DXPA
        Y1=PATCHRY(I)-DYPA
        Z1=PATCHRZ(I)-DZPA
        X2=X1+N1*DXPA
        Y2=Y1+N2*DYPA
        Z2=Z1+N3*DZPA
        DO II=1,NCLUS
         BASX=CLUSRX(II)
         BASY=CLUSRY(II)
         BASZ=CLUSRZ(II)
         IF (REALCLUS(II).EQ.0) CYCLE
         BAS=2.0*DXPA !0.1*RADIO(II)
         X3=BASX-BAS
         Y3=BASY-BAS
         Z3=BASZ-BAS
         X4=BASX+BAS
         Y4=BASY+BAS
         Z4=BASZ+BAS
         IF (X1.LE.X4.AND.X3.LE.X2.AND.
     &       Y1.LE.Y4.AND.Y3.LE.Y2.AND.
     &       Z1.LE.Z4.AND.Z3.LE.Z2) THEN
c           WRITE(*,*) BASX,BASY,BASZ,BAS
          DO KZ=1,N3
          DO JY=1,N2
          DO IX=1,N1
           AA=SQRT((RX(IX,I)-BASX)**2 +
     &             (RY(JY,I)-BASY)**2 +
     &             (RZ(KZ,I)-BASZ)**2)
           IF (AA.LE.BAS) CONTA1(IX,JY,KZ,I)=0
          END DO
          END DO
          END DO
         END IF
        END DO
       END DO

        ! clean conta1 from overlaps
       CALL CLEAN_OVERLAPS_INT(NL,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                         SOLAP,CONTA1)

*       CONTA1=0 --> the cell is outside any previous halo, or in the
*                    core of a previous halo, or overlaps another
*                    cell and is slave
*       CONTA1=1 --> outside any halo, and not overlapped (or master)

       ESP=0.2*DXPA
       ESP_LOG=1.05
       BORAMR=1
       BOUNDIR=MAX(BOUND/1.5**IR,2.0)

*        estimation to allocate
       KK_ENTERO=0
!$OMP  PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,U11,CONTRASTECPEAK,
!$OMP+                    LOW1,LOW2,CONTA1,BORAMR),
!$OMP+             PRIVATE(N1,N2,N3,I,IX,JY,KZ,BASX,BASY,BASZ,BASXX,
!$OMP+                     BASYY,BASZZ),
!$OMP+             REDUCTION(+:KK_ENTERO),
!$OMP+             DEFAULT(NONE)
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        DO KZ=1+BORAMR,N3-BORAMR
        DO JY=1+BORAMR,N2-BORAMR
        DO IX=1+BORAMR,N1-BORAMR
         IF (U11(IX,JY,KZ,I).GE.CONTRASTECPEAK.AND.
     &       CONTA1(IX,JY,KZ,I).EQ.1) THEN
          BASX =U11(IX+1,JY,KZ,I)-U11(IX,JY,KZ,I)
          BASY =U11(IX,JY+1,KZ,I)-U11(IX,JY,KZ,I)
          BASZ =U11(IX,JY,KZ+1,I)-U11(IX,JY,KZ,I)
          BASXX=U11(IX,JY,KZ,I)  -U11(IX-1,JY,KZ,I)
          BASYY=U11(IX,JY,KZ,I)  -U11(IX,JY-1,KZ,I)
          BASZZ=U11(IX,JY,KZ,I)  -U11(IX,JY,KZ-1,I)
          IF (BASX.LT.0) THEN
          IF (BASY.LT.0) THEN
          IF (BASZ.LT.0) THEN
          IF (BASXX.GT.0) THEN
          IF (BASYY.GT.0) THEN
          IF (BASZZ.GT.0) THEN ! then it's a local maximum
           KK_ENTERO=KK_ENTERO+1
          END IF
          END IF
          END IF
          END IF
          END IF
          END IF
         END IF
        END DO
        END DO
        END DO
       END DO
       WRITE(*,*) '... Clean cells above virial contrast',KK_ENTERO

       ALLOCATE(DDD(KK_ENTERO))
       ALLOCATE(DDDX(KK_ENTERO))
       ALLOCATE(DDDY(KK_ENTERO))
       ALLOCATE(DDDZ(KK_ENTERO))
       ALLOCATE(DDDP(KK_ENTERO))

!$OMP PARALLEL DO SHARED(DDD,DDDX,DDDY,DDDZ,DDDP,KK_ENTERO),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
       DO I=1,KK_ENTERO
        DDD(I)=0.0
        DDDX(I)=0
        DDDY(I)=0
        DDDZ(I)=0
        DDDP(I)=0
       END DO

       II=0
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        DO KZ=1+BORAMR,N3-BORAMR
        DO JY=1+BORAMR,N2-BORAMR
        DO IX=1+BORAMR,N1-BORAMR
         IF (U11(IX,JY,KZ,I).GE.CONTRASTECPEAK.AND.
     &       CONTA1(IX,JY,KZ,I).EQ.1) THEN
          BASX =U11(IX+1,JY,KZ,I)-U11(IX,JY,KZ,I)
          BASY =U11(IX,JY+1,KZ,I)-U11(IX,JY,KZ,I)
          BASZ =U11(IX,JY,KZ+1,I)-U11(IX,JY,KZ,I)
          BASXX=U11(IX,JY,KZ,I)  -U11(IX-1,JY,KZ,I)
          BASYY=U11(IX,JY,KZ,I)  -U11(IX,JY-1,KZ,I)
          BASZZ=U11(IX,JY,KZ,I)  -U11(IX,JY,KZ-1,I)
          IF (BASX.LT.0) THEN
          IF (BASY.LT.0) THEN
          IF (BASZ.LT.0) THEN
          IF (BASXX.GT.0) THEN
          IF (BASYY.GT.0) THEN
          IF (BASZZ.GT.0) THEN ! then it's a local maximum
           II=II+1
           DDD(II)=U11(IX,JY,KZ,I)
           DDDX(II)=IX
           DDDY(II)=JY
           DDDZ(II)=KZ
           DDDP(II)=I
          END IF
          END IF
          END IF
          END IF
          END IF
          END IF
         END IF
        END DO
        END DO
        END DO
       END DO

       CALL SORT_CELLS_AMR(KK_ENTERO,DDD,DDDX,DDDY,DDDZ,DDDP)

       NV_GOOD=KK_ENTERO
       KK_ENTERO=0

*        Go through all the center candidates
       DO L1=1,NV_GOOD
C        write(*,*)
C        write(*,*) 'new candidate!!!!'
        ICEN4(1)=DDDX(L1)
        ICEN4(2)=DDDY(L1)
        ICEN4(3)=DDDZ(L1)
        ICEN4(4)=DDDP(L1)
        IX=ICEN4(1)
        JY=ICEN4(2)
        KZ=ICEN4(3)
        IPATCH=ICEN4(4)

        XCEN=RX(IX,IPATCH)
        YCEN=RY(JY,IPATCH)
        ZCEN=RZ(KZ,IPATCH)

        IF (IR.NE.NL) THEN
         CALL RECENTER_DENSITY_PEAK(XCEN,YCEN,ZCEN,IR,U1,U11,PATCHRX,
     &                              PATCHRY,PATCHRZ,PATCHNX,PATCHNY,
     &                              PATCHNZ,NPATCH,LADO0,NL,IX,JY,KZ,
     &                              IPATCH)

*        Assert we have not yet identified this same halo (dupplicated
*         due to recentering)
         FLAG_DUP=0
         DO II=NCLUS_INI+1,NCLUS
          IF ((XCEN-CLUSRX(II))**2+(YCEN-CLUSRY(II))**2+
     &        (ZCEN-CLUSRZ(II))**2.LT.RSUB(II)**2) THEN
           FLAG_DUP=1
           EXIT
          END IF
         END DO
         IF (FLAG_DUP.EQ.1) CYCLE
        END IF ! (IR.NE.NL)

c        KK_ENTERO=CONTA1(IX,JY,KZ,IPATCH)
c        if (kk_entero.eq.0) then
c         write(*,*) 'kk_entero=0 at',xcen,ycen,zcen
c         CALL FIND_HOST_HALO(XCEN,YCEN,ZCEN,DXPA,IR,SUBS_LEV,REALCLUS,
c     &                       CLUSRX,CLUSRY,CLUSRZ,RADIO,IHOSTHALO)
c         BASXX=CLUSRX(IHOSTHALO)
c         BASYY=CLUSRY(IHOSTHALO)
c         BASZZ=CLUSRZ(IHOSTHALO)
c         MHOST=MASA(IHOSTHALO)
c
c         DISTHOST=SQRT((BASXX-XCEN)**2+(BASYY-YCEN)**2+(BASZZ-ZCEN)**2)
c         write(*,*) 'indeed, dista=',disthost
c        end if

ccc        IF (KK_ENTERO.EQ.1) THEN ! this means this peak is not inside a halo yet

        CALL FIND_HOST_HALO(XCEN,YCEN,ZCEN,DXPA,IR,SUBS_LEV,REALCLUS,
     &                      CLUSRX,CLUSRY,CLUSRZ,RADIO,RSUB,IHOSTHALO)

        IF (IHOSTHALO.EQ.-1) THEN
C         write(*,*) 'not found parent',xcen,ycen,zcen
         CYCLE
        END IF

        BASXX=CLUSRX(IHOSTHALO)
        BASYY=CLUSRY(IHOSTHALO)
        BASZZ=CLUSRZ(IHOSTHALO)
        IF (REALCLUS(IHOSTHALO).EQ.-1) THEN
         MHOST=MASA(IHOSTHALO)
        ELSE
         MHOST=MSUB(IHOSTHALO)
        END IF
        DISTHOST=SQRT((BASXX-XCEN)**2+(BASYY-YCEN)**2+(BASZZ-ZCEN)**2)

        IF (DISTHOST.LT.2.0*DXPA) CYCLE

        ! MHOST is the mass inside DISTHOST. Let's interpolate
        DO JJ=1,NBINS
         IF (PROFILES(JJ,1,IHOSTHALO).GT.DISTHOST) EXIT
        END DO
        IF (JJ.EQ.1) THEN
         MHOST=PROFILES(JJ,2,IHOSTHALO)*
     &                           (DISTHOST/PROFILES(JJ,1,IHOSTHALO))**3
        ELSE IF (JJ.LE.NBINS) THEN
         MHOST=PROFILES(JJ-1,2,IHOSTHALO)+
     & (PROFILES(JJ,2,IHOSTHALO)-PROFILES(JJ-1,2,IHOSTHALO)) *
     & (DISTHOST**3-PROFILES(JJ-1,1,IHOSTHALO)**3) /
     & (PROFILES(JJ,1,IHOSTHALO)**3-PROFILES(JJ-1,1,IHOSTHALO)**3)
        ELSE
         MHOST=PROFILES(NBINS,2,IHOSTHALO)
        END IF

c        WRITE(*,*) 'MHOST INTERP', PROFILES(NBINS,2,IHOSTHALO),
c     &             DISTHOST/PROFILES(NBINS,1,IHOSTHALO),MHOST

        NCLUS=NCLUS+1
        REALCLUS(NCLUS)=IHOSTHALO
        LEVHAL(NCLUS)=IR
        !NHALLEV(IR)=NHALLEV(IR)+1
        PATCHCLUS(NCLUS)=IPATCH

        CLUSRX(NCLUS)=XCEN
        CLUSRY(NCLUS)=YCEN
        CLUSRZ(NCLUS)=ZCEN

        IF(NCLUS.GT.MAXNCLUS) THEN
         WRITE(*,*) 'WARNING: NCLUS>MAXNCLUS!!!',NCLUS,MAXNCLUS
         STOP
        END IF

*       tentative reach of the base grid
        BASX=XCEN-BOUNDIR
        NX1=INT(((BASX-RADX(1))/DX)+0.5)+1
        IF (NX1.LT.1) NX1=1

        BASX=XCEN+BOUNDIR
        NX2=INT(((BASX-RADX(1))/DX)+0.5)+1
        IF (NX2.GT.NX) NX2=NX

        BASY=YCEN-BOUNDIR
        NY1=INT(((BASY-RADY(1))/DY)+0.5)+1
        IF (NY1.LT.1) NY1=1

        BASY=YCEN+BOUNDIR
        NY2=INT(((BASY-RADY(1))/DY)+0.5)+1
        IF (NY2.GT.NY) NY2=NY

        BASZ=ZCEN-BOUNDIR
        NZ1=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
        IF (NZ1.LT.1) NZ1=1

        BASZ=ZCEN+BOUNDIR
        NZ2=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
        IF (NZ2.GT.NZ) NZ2=NZ

*          patches that could contain this halo
        CALL PATCHES_SPHERE(NPATCH,PATCHRX,PATCHRY,PATCHRZ,PATCHNX,
     &                      PATCHNY,PATCHNZ,XCEN,YCEN,ZCEN,BOUNDIR,
     &                      NL,NL,RELEVANT_PATCHES,NRELEVANT_PATCHES)

c           WRITE(*,*) CLUSRX(NCLUS),NX1,NX2,RADX(NX1),RADX(NX2)
c           DO I=1,NRELEVANT_PATCHES(1)
c            WRITE(*,*) I,PATCHRX(RELEVANT_PATCHES(I))-DX/2.0,
c     &                 PATCHRX(RELEVANT_PATCHES(I))+
c     &                 (PATCHNX(RELEVANT_PATCHES(I))-1)*DX/2.0
c           END DO

*          Now, we extend the cluster radially from its center
        BASX=0.0
        BASY=0.0
        BASZ=0.0
        BASDELTA=0.0

        BASMASS_SHELL=0.0
        BASMASS=0.0   !TOTAL MASS OF THE CLUSTER
        DELTA=0.0    !TOTAL CONTRAST OF THE CLUSTER
        BASVOL=0.0     !TOTAL VOLUME OF THE CLUSTER (sphere)

        R_INT=0.0
        R_EXT=0.2*DXPA ! major semidiagonal of a cube (max distance to a cell)

*          increase the radius until density falls below the virial value
        DELTA=10.0*CONTRASTEC*ROTE ! this is to ensure we enter the loop
        ITER_GROW=0
        II=0
        FLAG_VIR=0
        FLAG_JACOBI=0
        JJ=0

         !WRITE(*,*) '**',NX1,NX2,NY1,NY2,NZ1,NZ2,NRELEVANT_PATCHES(1:NL)

        DO WHILE(FLAG_JACOBI.EQ.0)
         ITER_GROW=ITER_GROW+1

         IF (ITER_GROW.GT.1) THEN
          IF (R_EXT.LE.BOUNDIR) THEN
           R_INT=R_EXT
           R_EXT=MAX(R_EXT+ESP, R_EXT*ESP_LOG)
          ELSE
           WRITE(*,*) 'WARNING: growing not converged', r_int, r_ext,
     &                 boundir,iter_grow,delta/rote,kk_entero,
     &                 xcen,ycen,zcen
           STOP
          END IF
         END IF

         BASMASS_SHELL=0.0

         VOLCELL=DX*DY*DZ
         DO K=NZ1,NZ2
         DO J=NY1,NY2
         DO I=NX1,NX2
          IF (CR0AMR(I,J,K).EQ.1) THEN
           AA=SQRT((RADX(I)-XCEN)**2 + (RADY(J)-YCEN)**2 +
     &             (RADZ(K)-ZCEN)**2)

           IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
            !CONTA(I,J,K)=1 ! do not try to find an additional halo here (at this level)
            II=II+1
            BASVOL=BASVOL+VOLCELL

            BAS=U1(I,J,K)*VOLCELL !U1 is not density contrast, but 1+delta = rho/rho_B!!!
            BASMASS_SHELL=BASMASS_SHELL+BAS

            BASX=BASX+RADX(I)*BAS
            BASY=BASY+RADY(J)*BAS
            BASZ=BASZ+RADZ(K)*BAS
            BASDELTA=BASDELTA+BAS
           END IF
          END IF
         END DO
         END DO
         END DO

         ! AND NOW LOOP UP THROUGH REFINEMENT LEVELS
         DO IRR=1,NL
          DXPAPA=DX/2.0**IRR
          DYPAPA=DY/2.0**IRR
          DZPAPA=DZ/2.0**IRR
          LOW1=SUM(NRELEVANT_PATCHES(1:IRR-1))+1
          LOW2=SUM(NRELEVANT_PATCHES(1:IRR))
          VOLCELL=DXPAPA*DYPAPA*DZPAPA
          DO BASINT=LOW1,LOW2
           IPATCH=RELEVANT_PATCHES(BASINT)
           N1=PATCHNX(IPATCH)
           N2=PATCHNY(IPATCH)
           N3=PATCHNZ(IPATCH)
           DO KZ=1,N3
           DO JY=1,N2
           DO IX=1,N1
            IF (CR0AMR11(IX,JY,KZ,IPATCH).EQ.1) THEN
             IF (SOLAP(IX,JY,KZ,IPATCH).EQ.1) THEN
              ! DO STUFF
              AA=SQRT((RX(IX,IPATCH)-XCEN)**2 +
     &                (RY(JY,IPATCH)-YCEN)**2 +
     &                (RZ(KZ,IPATCH)-ZCEN)**2)

              IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
               II=II+1
               BASVOL=BASVOL+VOLCELL

               BAS=U11(IX,JY,KZ,IPATCH)*VOLCELL !U1 is not density contrast, but 1+delta = rho/rho_B!!!
               BASMASS_SHELL=BASMASS_SHELL+BAS

               BASX=BASX+RX(IX,IPATCH)*BAS
               BASY=BASY+RY(JY,IPATCH)*BAS
               BASZ=BASZ+RZ(KZ,IPATCH)*BAS
               BASDELTA=BASDELTA+BAS

               CONTA1(IX,JY,KZ,IPATCH)=0
              END IF
             END IF
            END IF
           END DO
           END DO
           END DO
          END DO
         END DO

         BASMASS=BASMASS+BASMASS_SHELL*RODO*RE0**3
         IF (BASVOL.GT.0.0) DELTA=BASMASS/(BASVOL*RETE**3)

         IF (FLAG_VIR.EQ.0.AND.DELTA.LT.CONTRASTEC*ROTE) THEN
          FLAG_VIR=1
          MASA(NCLUS)=BASMASS*UM
          RADIO(NCLUS)=R_EXT*UM
         END IF

         EQ_JACOBI_R=(R_EXT/DISTHOST)-
     &                  ((BASMASS*UM)/(BASMASS*UM+3*MHOST))**(1.0/3.0)
C         WRITE(*,*) L1,ITER_GROW,R_EXT,BASMASS*UM,DELTA/ROTE,
C     &               II
C         WRITE(*,*) '... ',EQ_JACOBI_R,R_EXT/DISTHOST,
C     &          ((BASMASS*UM)/(BASMASS*UM+3*MHOST))**(1.0/3.0)

         IF (EQ_JACOBI_R.GT.0.0.AND.II.GT.15.AND.ITER_GROW.GE.4) THEN
          FLAG_JACOBI=1
          RSUB(NCLUS)=R_EXT
          MSUB(NCLUS)=BASMASS*UM
         END IF

        END DO   ! do while (DELTA)

        CLUSRXCM(NCLUS)=BASX/BASDELTA
        CLUSRYCM(NCLUS)=BASY/BASDELTA
        CLUSRZCM(NCLUS)=BASZ/BASDELTA

       END DO ! L1=1,NV_GOOD

       DEALLOCATE(DDD,DDDX,DDDY,DDDZ,DDDP)

       SUBS_LEV(IR)=NCLUS-NCLUS_INI
       WRITE(*,*) '... Halos found at level:',IR,SUBS_LEV(IR)

       RETURN
       END

********************************************************************
       SUBROUTINE FIND_HOST_HALO(BASX,BASY,BASZ,DXPA,IR,SUBS_LEV,
     &                           REALCLUS,CLUSRX,CLUSRY,CLUSRZ,RADIO,
     &                           RSUB,IHOSTHALO)
********************************************************************
*      Find the host halo (from deeper levels to coarser levels of
*       substructure)
********************************************************************

       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       REAL BASX,BASY,BASZ,DXPA
       INTEGER IR
       INTEGER SUBS_LEV(0:NLEVELS)
       INTEGER REALCLUS(MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
       REAL RADIO(MAXNCLUS),RSUB(MAXNCLUS)
       INTEGER IHOSTHALO

       INTEGER FLAG,IRR,LOW1,LOW2,II
       REAL MINDISTA,BASXX,BASYY,BASZZ,BAS,DISTA

       FLAG=0
       MINDISTA=100000.0

       DO IRR=IR-1,0,-1

        LOW1=SUM(SUBS_LEV(0:IRR-1))+1
        LOW2=SUM(SUBS_LEV(0:IRR))

        DO II=LOW1,LOW2

         IF (REALCLUS(II).EQ.0) CYCLE

         BASXX=CLUSRX(II)
         BASYY=CLUSRY(II)
         BASZZ=CLUSRZ(II)
         IF (REALCLUS(II).EQ.-1) THEN
          BAS=RADIO(II)
         ELSE
          BAS=RSUB(II)
         END IF

         DISTA=SQRT((BASX-BASXX)**2+(BASY-BASYY)**2+(BASZ-BASZZ)**2)
     &           / BAS

         IF (DISTA.LT.1.0) THEN
C          WRITE(*,*) 'PROG',II,BASX,BASY,BASZ,BASXX,BASYY,BASZZ,BAS,
C     &               DISTA
          IF (FLAG.EQ.0) FLAG=1
          IF (DISTA.LT.MINDISTA) THEN
           MINDISTA=DISTA
           IHOSTHALO=II
          END IF
         END IF !(DISTA.LT.1.0)
        END DO !II=LOW1,LOW2
        IF (FLAG.EQ.1) EXIT
       END DO !IRR=IR-1,0,-1

       IF (FLAG.NE.1) THEN
C        WRITE(*,*) 'Not found progenitor'
        IHOSTHALO=-1
       END IF

       RETURN
       END


**********************************************************************
       SUBROUTINE SUBSTRUCTURE_PARTICLES(IR,NL,NCLUS,MASA,RADIO,CLUSRX,
     &      CLUSRY,CLUSRZ,REALCLUS,CONCENTRA,ANGULARM,VMAXCLUS,IPLIP,VX,
     &      VY,VZ,VCMAX,MCMAX,RCMAX,M200C,M500C,M2500C,M200M,M500M,
     &      M2500M,MSUB,R200C,R500C,R2500C,R200M,R500M,R2500M,RSUB,
     &      DMPCLUS,LEVHAL,EIGENVAL,N_DM,RXPA,RYPA,RZPA,MASAP,U2DM,U3DM,
     &      U4DM,ORIPA2,CONTRASTEC,OMEGAZ,UM,UV,LADO0,CLUSRXCM,CLUSRYCM,
     &      CLUSRZCM,MEAN_VR,INERTIA_TENSOR,SUBS_LEV,PATCHCLUS,NPATCH,
     &      PROFILES)
**********************************************************************
*      Refines halo identification with DM particles
**********************************************************************
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER IR,NL,NCLUS
       REAL*4 MASA(MAXNCLUS),RADIO(MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
       INTEGER REALCLUS(MAXNCLUS)
       REAL*4 CONCENTRA(NMAXNCLUS),ANGULARM(3,NMAXNCLUS)
       REAL*4 VMAXCLUS(NMAXNCLUS)
       INTEGER IPLIP(NMAXNCLUS)
       REAL*4 VX(NMAXNCLUS),VY(NMAXNCLUS),VZ(NMAXNCLUS)
       REAL*4 VCMAX(NMAXNCLUS),MCMAX(NMAXNCLUS),RCMAX(NMAXNCLUS)
       REAL*4 M200C(NMAXNCLUS),R200C(NMAXNCLUS)
       REAL*4 M500C(NMAXNCLUS),R500C(NMAXNCLUS)
       REAL*4 M2500C(NMAXNCLUS),R2500C(NMAXNCLUS)
       REAL*4 M200M(NMAXNCLUS),R200M(NMAXNCLUS)
       REAL*4 M500M(NMAXNCLUS),R500M(NMAXNCLUS)
       REAL*4 M2500M(NMAXNCLUS),R2500M(NMAXNCLUS)
       REAL*4 MSUB(MAXNCLUS),RSUB(MAXNCLUS)
       INTEGER DMPCLUS(NMAXNCLUS),LEVHAL(MAXNCLUS)
       REAL*4 EIGENVAL(3,NMAXNCLUS)
       INTEGER N_DM
       REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED)
       REAL*4 U2DM(PARTIRED),U3DM(PARTIRED),U4DM(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       INTEGER ORIPA2(PARTIRED)
       REAL*4 CONTRASTEC,OMEGAZ,UM,UV,LADO0
       REAL*4 CLUSRXCM(MAXNCLUS),CLUSRYCM(MAXNCLUS),CLUSRZCM(MAXNCLUS)
       REAL*4 MEAN_VR(NMAXNCLUS),INERTIA_TENSOR(6,NMAXNCLUS)
       INTEGER SUBS_LEV(0:NLEVELS),PATCHCLUS(MAXNCLUS),NPATCH(0:NLEVELS)
       REAL*4 PROFILES(NBINS,2,NMAXNCLUS)

       REAL*4 PI,ACHE,T0,RE0,PI4ROD
       COMMON /DOS/ACHE,T0,RE0
       REAL*4 UNTERCIO,CGR,CGR2,ZI,RODO,ROI,REI
       COMMON /CONS/PI4ROD,REI,CGR,PI
       REAL*4 OMEGA0

       REAL*4 RETE,HTE,ROTE
       COMMON /BACK/ RETE,HTE,ROTE

       REAL*4 DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

*      Local variables
       INTEGER DIMEN,KONTA1,KONTA2,I,PABAS,NUMPARTBAS,KK_ENTERO,NROT,II
       INTEGER KONTA,J,CONTAERR,JJ,SALIDA,KONTA3,NSHELL_2,KONTA2PREV
       INTEGER IX,JY,KK1,KK2,FAC,ITER_SHRINK,COUNT_1,COUNT_2,JJCORE
       INTEGER FLAG200C,FLAG500C,FLAG2500C,FLAG200M,FLAG500M,FLAG2500M
       INTEGER FLAGVIR,JMINPROF,LOWH1,LOWH2,IHOSTHALO,FLAG_JACOBI
       INTEGER NCAPAS(NMAXNCLUS),IPATCH,IRR,MINJ,DIR,J_JACOBI,EACH_PROF
       INTEGER LIP(PARTIRED),CONTADM(PARTIRED)
       REAL ESP,REF,REF_MIN,REF_MAX,MASADM,BASMAS,DIS,VCM,MINOVERDENS
       REAL VVV2,VR,CONCEN,RS,BAS,AADM,CMX,CMY,CMZ,VCMX,VCMY,CX,CY,CZ
       REAL VCMZ,MASA2,NORMA,BAS1,BAS2,VOL,DELTA2,RSHELL,RCLUS,BASVEC(3)
       REAL DENSA,DENSB,DENSC,VKK,AA,BASX,BASY,BASZ,XP,YP,ZP,MP
       REAL INERTIA(3,3),BASEIGENVAL(3),RADII_ITER(7),BASVECCM(3)
       REAL BASVCM(3),MHOST,DISTHOST,EQ_JACOBI_R
       REAL DENSITOT(0:1000),RADIAL(0:1000),DENSR(0:1000)
       REAL LOGDERIV(0:1000)
       REAL DISTA(0:PARTIRED)

       PI=DACOS(-1.D0)

       DIMEN=3   !DIMENSION DE LOS HALOS
       NCAPAS=0
       ESP=0.0
       REF=0.0
       KONTA1=0
       KONTA2=0

       MINOVERDENS=200.0

       LOWH1=SUM(SUBS_LEV(0:IR-1))+1
       LOWH2=SUM(SUBS_LEV(0:IR))
       WRITE(*,*) 'Refining with particles, haloes:',LOWH1,LOWH2,
     &                                               LOWH2-LOWH1+1

       PABAS=PARTIRED_PLOT
       NUMPARTBAS=NUMPART
       NORMA=MAXVAL(MASAP)

!$OMP  PARALLEL DO SHARED(NCLUS,REALCLUS,
!$OMP+           LEVHAL,RXPA,RYPA,RZPA,CLUSRX,CLUSRY,CLUSRZ,NL,MASAP,
!$OMP+           U2DM,U3DM,U4DM,VX,VY,VZ,ACHE,PI,RETE,ROTE,VCMAX,
!$OMP+           MCMAX,RCMAX,CONTRASTEC,OMEGAZ,CGR,UM,UV,DMPCLUS,
!$OMP+           CONCENTRA,ORIPA2,ANGULARM,PABAS,IPLIP,DIMEN,EIGENVAL,
!$OMP+           NUMPARTBAS,RADIO,MASA,VMAXCLUS,N_DM,NORMA,R200M,
!$OMP+           R500M,R2500M,R200C,R500C,R2500C,M200M,M500M,M2500M,
!$OMP+           M200C,M500C,M2500C,RSUB,MSUB,MINOVERDENS,CLUSRXCM,
!$OMP+           CLUSRYCM,CLUSRZCM,DX,MEAN_VR,INERTIA_TENSOR,LOWH1,
!$OMP+           LOWH2,PATCHCLUS,NPATCH,PROFILES,NCAPAS),
!$OMP+   PRIVATE(I,INERTIA,REF_MIN,REF_MAX,KK_ENTERO,MASADM,KONTA,
!$OMP+           BASMAS,DIS,VCM,VVV2,VR,LIP,CONCEN,RS,KONTA2,BAS,J,
!$OMP+           AADM,KK1,KK2,CONTADM,CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MASA2,
!$OMP+           DISTA,FAC,CONTAERR,JJ,DENSITOT,RADIAL,SALIDA,BAS1,BAS2,
!$OMP+           VOL,DELTA2,RSHELL,KONTA3,NSHELL_2,KONTA1,
!$OMP+           DENSA,DENSB,DENSC,BASVEC,BASVECCM,VKK,AA,NROT,
!$OMP+           BASEIGENVAL,BASX,BASY,BASZ,XP,YP,ZP,MP,RCLUS,COUNT_1,
!$OMP+           COUNT_2,KONTA2PREV,FLAG200C,FLAG200M,FLAG500C,FLAG500M,
!$OMP+           FLAG2500C,FLAG2500M,FLAGVIR,JMINPROF,DENSR,LOGDERIV,
!$OMP+           CX,CY,CZ,JJCORE,RADII_ITER,BASVCM,IHOSTHALO,MINJ,DIR,
!$OMP+           FLAG_JACOBI,MHOST,DISTHOST,EQ_JACOBI_R,IPATCH,IRR,
!$OMP+           J_JACOBI,EACH_PROF),
!$OMP+   SCHEDULE(DYNAMIC), DEFAULT(NONE)
*****************************
       DO I=LOWH1,LOWH2
****************************
       KK_ENTERO=REALCLUS(I)
       IF (KK_ENTERO.NE.0) THEN

        INERTIA=0.0
        REF_MIN=10.0e+10
        REF_MAX=-1.0
        MASADM=0.0
        KONTA=0
        BASMAS=0.0
        DIS=1.0E+10    !Distance (to the center) of the most central particle
        VCM=0.0
        VVV2=0.0    ! v propia de las particulas respecto al CM
        VR=0.0      ! v radial
        CMX=0.0
        CMY=0.0
        CMZ=0.0
        CX=0.0
        CY=0.0
        CZ=0.0
        VCMX=0.0
        VCMY=0.0
        VCMZ=0.0
        LIP=0
        CONCEN=0.0       !concentration NFW profile
        RS=0.0           !scale radius NFW profile
        KONTA2=0
        BASX=0.0
        BASY=0.0
        BASZ=0.0

*********************************************************************
*       RECENTERING AND COMPUTING VCM OF HALO I (SHRINKING SPHERE)
*********************************************************************

        CX=CLUSRX(I)
        CY=CLUSRY(I)
        CZ=CLUSRZ(I)
        BAS=RSUB(I)
        IHOSTHALO=REALCLUS(I)
        IF (REALCLUS(IHOSTHALO).EQ.-1) THEN
         MHOST=MASA(IHOSTHALO)
        ELSE
         MHOST=MSUB(IHOSTHALO)
        END IF
        DISTHOST=SQRT((CLUSRX(IHOSTHALO)-CX)**2 +
     &                (CLUSRY(IHOSTHALO)-CY)**2 +
     &                (CLUSRZ(IHOSTHALO)-CZ)**2)


        ! MHOST is the mass inside DISTHOST. Let's interpolate
        DO JJ=1,NBINS
         IF (PROFILES(JJ,1,IHOSTHALO).GT.DISTHOST) EXIT
        END DO
        IF (JJ.EQ.1) THEN
         MHOST=PROFILES(JJ,2,IHOSTHALO)*
     &                           (DISTHOST/PROFILES(JJ,1,IHOSTHALO))**3
        ELSE IF (JJ.LE.NBINS) THEN
         MHOST=PROFILES(JJ-1,2,IHOSTHALO)+
     & (PROFILES(JJ,2,IHOSTHALO)-PROFILES(JJ-1,2,IHOSTHALO)) *
     & (DISTHOST**3-PROFILES(JJ-1,1,IHOSTHALO)**3) /
     & (PROFILES(JJ,1,IHOSTHALO)**3-PROFILES(JJ-1,1,IHOSTHALO)**3)
        ELSE
         MHOST=PROFILES(NBINS,2,IHOSTHALO)
        END IF

c        WRITE(*,*) 'MHOST INTERP', PROFILES(NBINS,2,IHOSTHALO),
c     &             DISTHOST/PROFILES(NBINS,1,IHOSTHALO),MHOST

        ! Find level used for recentering
        IPATCH=PATCHCLUS(I)
        DO IRR=1,NL
         IF (SUM(NPATCH(0:IRR-1))+1.LE.IPATCH.AND.
     &       IPATCH.LE.SUM(NPATCH(0:IRR))) EXIT
        END DO
c        WRITE(*,*) CX,CY,CZ,IRR
        CALL RECENTER_DENSITY_PEAK_PARTICLES(CX,CY,CZ,BAS,RXPA,RYPA,
     &                                       RZPA,MASAP,N_DM,
     &                                       DX/2.0**IRR)

        BAS=(CLUSRX(I)-CX)**2+(CLUSRY(I)-CY)**2+(CLUSRZ(I)-CZ)**2
        BAS=SQRT(BAS)
c        WRITE(*,*) 'Recentering shift', i, bas, bas/rsub(i)

        DISTHOST=SQRT((CLUSRX(IHOSTHALO)-CX)**2 +
     &                (CLUSRY(IHOSTHALO)-CY)**2 +
     &                (CLUSRZ(IHOSTHALO)-CZ)**2)
c        write(*,*) 'mhost,dishost',mhost,disthost

        CLUSRX(I)=CX
        CLUSRY(I)=CY
        CLUSRZ(I)=CZ

*       Find jacobi radius (substructure radius; grow if necessary)
        DELTA2=100.0*MINOVERDENS ! just to ensure it enters the loop
        FLAG_JACOBI=0
        RCLUS=RSUB(I)
        DO WHILE (FLAG_JACOBI.EQ.0)
         KONTA=0
         MASADM=0.0
         DO J=1,N_DM
          AADM=SQRT((RXPA(J)-CX)**2+(RYPA(J)-CY)**2+(RZPA(J)-CZ)**2)
          IF (AADM.LT.RCLUS) THEN
           KONTA=KONTA+1
           LIP(KONTA)=J
           MASADM=MASADM+MASAP(J)
          END IF
         END DO
         CONTADM=1
         CONTADM(1:KONTA)=0

         CALL REORDENAR(KONTA,CX,CY,CZ,RXPA,RYPA,RZPA,CONTADM,LIP,
     &                  DISTA,KONTA2,1)

c         WRITE(*,*) 'CHECK I,KONTA,KONTA2=',I,KONTA,KONTA2

         DELTA2=MASADM/(ROTE*RETE**3*(4*PI/3)*RCLUS**3)

*        Simplified equation (estimation)
         MASADM=0.0
         MINJ=MAX(10,KONTA/50)
         DO J=1,KONTA
          MASADM=MASADM+MASAP(LIP(J))
          EQ_JACOBI_R=(DISTA(J)/DISTHOST)-
     &                  ((MASADM*UM)/(MASADM*UM+3*MHOST))**(1.0/3.0)
c          WRITE(*,*) J,DISTA(J),MASADM*UM,EQ_JACOBI_R
          IF (EQ_JACOBI_R.GT.0.0.AND.J.GT.MINJ) THEN
           J_JACOBI=J
           EXIT
          END IF
         END DO

*        Exact equation (refining the solution)
         BASX=DISTA(J)/DISTHOST
         BASY=MASADM*UM/MHOST
         EQ_JACOBI_R=1/(1-BASX)**2-1-BASY/BASX**2+(1+BASY)*BASX
         IF (ABS(EQ_JACOBI_R).LT.0.01) THEN
          FLAG_JACOBI=1
         ELSE
          DIR=1
          IF (EQ_JACOBI_R.LT.0.0) DIR=-1

          IF (DIR.EQ.1) THEN
           DO J=J_JACOBI+1,KONTA
            BASX=DISTA(J)/DISTHOST
            MASADM=MASADM+MASAP(LIP(J))
            BASY=MASADM*UM/MHOST
            EQ_JACOBI_R=1/(1-BASX)**2-1-BASY/BASX**2+(1+BASY)*BASX
            IF (EQ_JACOBI_R.GT.0.0) THEN
             J_JACOBI=J
             FLAG_JACOBI=1
             EXIT
            END IF
           END DO
          ELSE
           DO J=J_JACOBI-1,MINJ,-1
            BASX=DISTA(J)/DISTHOST
            MASADM=MASADM-MASAP(LIP(J+1))
            BASY=MASADM*UM/MHOST
            EQ_JACOBI_R=1/(1-BASX)**2-1-BASY/BASX**2+(1+BASY)*BASX
c            WRITE(*,*) J,DISTA(J),MASADM*UM,EQ_JACOBI_R
            IF (EQ_JACOBI_R.LT.0.0) THEN
             J_JACOBI=J
             FLAG_JACOBI=1
             EXIT
            END IF
           END DO
          END IF
         END IF !(ABS(EQ_JACOBI_R).LT.0.01) / ELSE

         IF (FLAG_JACOBI.EQ.1) THEN
          RCLUS=DISTA(J_JACOBI)
          RSUB(I)=RCLUS
          MSUB(I)=MASADM*UM
          CONTADM(1:KONTA)=1
          CONTADM(1:J_JACOBI)=0
          KONTA=J_JACOBI
         END IF

         IF (FLAG_JACOBI.EQ.0) THEN
          RCLUS=1.1*RCLUS
         END IF
        END DO
c        WRITE(*,*) DELTA2,MASADM*UM,RCLUS,KONTA
c        WRITE(*,*) 'THUS, RJ,MJ,KONTA=',RCLUS,MASADM*UM,KONTA

        IF (MASADM.LE.0.0.OR.KONTA.EQ.0) THEN
         REALCLUS(I)=0
         CYCLE
        END IF

        CONTADM(1:KONTA)=0     !en principio todas estas ligadas
        CALL CENTROMASAS_PART(KONTA,CONTADM,LIP,U2DM,U3DM,U4DM,MASAP,
     &                        RXPA,RYPA,RZPA,CMX,CMY,CMZ,VCMX,VCMY,VCMZ,
     &                        MASA2)

        VCM=SQRT(VCMX**2+VCMY**2+VCMZ**2)

        CLUSRXCM(I)=CMX
        CLUSRYCM(I)=CMY
        CLUSRZCM(I)=CMZ
        VX(I)=VCMX
        VY(I)=VCMY
        VZ(I)=VCMZ
        RSUB(I)=RCLUS
        MSUB(I)=MASADM*UM
        DMPCLUS(I)=KONTA
c        WRITE(*,*) '*',sqrt((cmx-cx)**2+(cmy-cy)**2+(cmz-cz)**2)
c        write(*,*) '**',vcmx,vcmy,vcmz,vcm

*********************************************************************
*       END RECENTERING AND COMPUTING VCM OF HALO I (SHRINKING SPHERE)
*********************************************************************

********************************************************************
*      UNBINDING:SCAPE VELOCITY
********************************************************************

        CONTAERR=KONTA
        DISTA=0.0
        KONTA2=0
        CALL REORDENAR(KONTA,CX,CY,CZ,RXPA,RYPA,RZPA,CONTADM,LIP,
     &                 DISTA,KONTA2,1)
        REF_MAX=DISTA(KONTA2)
        REF_MIN=DISTA(1)

        FAC=0
        DO WHILE (CONTAERR.GT.0.OR.FAC.LT.3)
         FAC=FAC+1
         KONTA2PREV=KONTA2
         CALL UNBINDING8(FAC,I,REF_MIN,REF_MAX,DISTA,U2DM,U3DM,U4DM,
     &                   MASAP,RXPA,RYPA,RZPA,RSUB,MSUB,CLUSRXCM,
     &                   CLUSRYCM,CLUSRZCM,LIP,KONTA,CONTADM,VX,VY,VZ,
     &                   REALCLUS,KONTA2)
         CALL REORDENAR(KONTA,CX,CY,CZ,RXPA,RYPA,RZPA,CONTADM,LIP,
     &                  DISTA,KONTA2,0)
         REF_MAX=DISTA(KONTA2)
         REF_MIN=DISTA(1)
         CONTAERR=KONTA2PREV-KONTA2
        END DO

        count_1=konta-konta2
        count_2=konta2 !backup
c        write(*,*) 'Unbinding V_ESC',i,'. ',konta,'-->',konta2,
c     &             '. Pruned:',count_1,'. Iters:', FAC

********************************************************************
*      UNBINDING: PHASE SPACE
********************************************************************

        FAC=0
        CONTAERR=KONTA2
        DO WHILE (CONTAERR.GT.0.OR.FAC.LT.4)
         FAC=FAC+1
         KONTA2PREV=KONTA2
         CALL UNBINDING_SIGMA(FAC,I,REF_MIN,REF_MAX,U2DM,U3DM,U4DM,RXPA,
     &                        RYPA,RZPA,MASAP,RSUB,MSUB,CLUSRXCM,
     &                        CLUSRYCM,CLUSRZCM,LIP,KONTA,CONTADM,VX,
     &                        VY,VZ,REALCLUS,KONTA2)
         CALL REORDENAR(KONTA,CX,CY,CZ,RXPA,RYPA,RZPA,CONTADM,LIP,
     &                  DISTA,KONTA2,0)
         REF_MAX=DISTA(KONTA2)
         REF_MIN=DISTA(1)
         CONTAERR=KONTA2PREV-KONTA2
         !write(*,*) 'sigma unbinding: iter,unbound',fac,contaerr
        END DO

        count_2=count_2-konta2
c        write(*,*) 'Unbinding SIGMA',i,'. ',konta,'-->',konta2,
c     &             '. Pruned:',count_2,'. Iters:', FAC
c        write(*,*) '--'

********************************************************************
*      CORRECT RSUB, MSUB AFTER UNBINDING
********************************************************************
        MASADM=MSUB(I)/UM
        BASX=DISTA(KONTA2)/DISTHOST
        BASY=MASADM*UM/MHOST
        EQ_JACOBI_R=1/(1-BASX)**2-1-BASY/BASX**2+(1+BASY)*BASX
c        WRITE(*,*) 'after unbinding, eq_jacobi:',eq_jacobi_r
        IF (EQ_JACOBI_R.GT.0.0) THEN
         DO J=KONTA2,1,-1
          MASADM=MASADM-MASAP(LIP(J))
          BASX=DISTA(J)/DISTHOST
          BASY=MASADM*UM/MHOST
          EQ_JACOBI_R=1/(1-BASX)**2-1-BASY/BASX**2+(1+BASY)*BASX
          IF (EQ_JACOBI_R.LT.0.0) EXIT
         END DO
         RSUB(I)=DISTA(J+1)
         MSUB(I)=MASADM*UM
         KONTA2=J
         NCAPAS(I)=J
c         WRITE(*,*) 'FINALLY, RJ, MJ=',rsub(i),msub(i),konta2
         IF (J.EQ.KONTA2) THEN
          RSHELL=DISTA(J)
         ELSE
          RSHELL=DISTA(J+1)
         END IF
        END IF

********************************************************************
*      DISCARD POOR HALOES
********************************************************************
        IF (KONTA2.LT.NUMPARTBAS) THEN
         REALCLUS(I)=0
        ELSE
**************************************************************
*      DENSITY PROFILE
**************************************************************
         REF_MIN=DISTA(1)
         REF_MAX=DISTA(KONTA2) ! DISTANCE TO THE FURTHERST BOUND PARTICLE
         JJ=0
         DENSITOT=0.0
         RADIAL=0.0
****************** stopping condition
         SALIDA=0
******************
         BAS=0.0 ! to store cummulative mass profile
         BAS1=-1.0 ! to store vcmax
         BAS2=0.0 ! to store vc(r)
*        VCMAX=-1.0
         ! Initialise flags (whether each overdensity has been found)
         FLAG200C=0
         FLAG500C=0
         FLAG2500C=0
         FLAG200M=0
         FLAG500M=0
         FLAG2500M=0
         FLAGVIR=0

         FAC=MAX(100,INT(0.05*KONTA2))
         KONTA3=INT(REAL(KONTA2)/FAC)*FAC
         NSHELL_2=0
         JJCORE=INT(MAX(10.0,KONTA2/100.0))
         DO J=1,JJCORE
          BAS=BAS+(MASAP(LIP(J))/NORMA)
         END DO
         DO J=JJCORE+1,KONTA2      !!!!! DEJO 80 por ciento de BINS DE SEGURIDAD
          VOL=PI*(4.0/3.0)*(DISTA(J)*RETE)**3
          BAS=BAS+(MASAP(LIP(J))/NORMA)

          DELTA2=NORMA*BAS/VOL/ROTE ! overdensity

          BAS2=BAS/DISTA(J) ! vc(r)
          IF (BAS2.GT.BAS1) THEN
           VCMAX(I)=BAS2
           MCMAX(I)=BAS
           RCMAX(I)=DISTA(J)
           BAS1=BAS2
          END IF

          IF (FLAG2500C.EQ.0) THEN
           IF (DELTA2.LE.2500.0/OMEGAZ) THEN
            M2500C(I)=DELTA2*VOL*ROTE*UM
            R2500C(I)=DISTA(J)
            FLAG2500C=1
            IF (J.EQ.JJCORE+1) THEN
             M2500C(I)=-1.0
             R2500C(I)=-1.0
            END IF
           END IF
          END IF

          IF (FLAG500C.EQ.0) THEN
           IF (DELTA2.LE.500.0/OMEGAZ) THEN
            M500C(I)=DELTA2*VOL*ROTE*UM
            R500C(I)=DISTA(J)
            FLAG500C=1
            IF (J.EQ.JJCORE+1) THEN
             M500C(I)=-1.0
             R500C(I)=-1.0
            END IF
           END IF
          END IF

          IF (FLAG200C.EQ.0) THEN
           IF (DELTA2.LE.200.0/OMEGAZ) THEN
            M200C(I)=DELTA2*VOL*ROTE*UM
            R200C(I)=DISTA(J)
            FLAG200C=1
            IF (J.EQ.JJCORE+1) THEN
             M200C(I)=-1.0
             R200C(I)=-1.0
            END IF
           END IF
          END IF

          IF (FLAG2500M.EQ.0) THEN
           IF (DELTA2.LE.2500.0) THEN
            M2500M(I)=DELTA2*VOL*ROTE*UM
            R2500M(I)=DISTA(J)
            FLAG2500M=1
            IF (J.EQ.JJCORE+1) THEN
             M2500M(I)=-1.0
             R2500M(I)=-1.0
            END IF
           END IF
          END IF

          IF (FLAG500M.EQ.0) THEN
           IF (DELTA2.LE.500.0) THEN
            M500M(I)=DELTA2*VOL*ROTE*UM
            R500M(I)=DISTA(J)
            FLAG500M=1
            IF (J.EQ.JJCORE+1) THEN
             M500M(I)=-1.0
             R500M(I)=-1.0
            END IF
           END IF
          END IF

          IF (FLAG200M.EQ.0) THEN
           IF (DELTA2.LE.200.0) THEN
            M200M(I)=DELTA2*VOL*ROTE*UM
            R200M(I)=DISTA(J)
            FLAG200M=1
            IF (J.EQ.JJCORE+1) THEN
             M200M(I)=-1.0
             R200M(I)=-1.0
            END IF
           END IF
          END IF

          ! STOPPING CONDITIONS
          ! 1. Below virial overdensity
          IF (FLAGVIR.EQ.0) THEN
           IF (DELTA2.LE.CONTRASTEC.AND.J.GT.INT(0.1*KONTA2)) THEN
            SALIDA=1 ! this means we've found the virial radius
                     ! but we don't exit straightaway because we still
                     ! want to find 200m
            FLAGVIR=1
            MASA(I)=DELTA2*VOL*ROTE*UM
            RADIO(I)=DISTA(J)
           END IF
          END IF

          ! profile
          IF (MOD(J,FAC).EQ.0) THEN
           NSHELL_2=NSHELL_2+1
           DENSITOT(NSHELL_2)=NORMA*BAS*UM
           RADIAL(NSHELL_2)=DISTA(J)
          END IF

         END DO ! J=1,KONTA2

         !WRITE(*,*) RADIAL(1:NSHELL_2)
         !WRITE(*,*) DENSITOT(1:NSHELL_2)

         IF (MOD(J,FAC).NE.0) THEN
          NSHELL_2=NSHELL_2+1
          DENSITOT(NSHELL_2)=NORMA*BAS*UM
          RADIAL(NSHELL_2)=DISTA(J)
         END IF

c         WRITE(*,*) 'HALO I,KONTA2,NSHELL_2,KK_ENTERO=',
c     &               I,KONTA2,NSHELL_2,KK_ENTERO
c         WRITE(*,*) CLUSRX(I),CLUSRY(I),CLUSRZ(I)
c         WRITE(*,*) R2500C(I),R500C(I),R200C(I),R2500M(I),R500M(I),
c     &              R200M(I),RADIO(I)
c         WRITE(*,*) M2500C(I),M500C(I),M200C(I),M2500M(I),M500M(I),
c     &              M200M(I),MASA(I)
c         WRITE(*,*) '---'

         BAS=VCMAX(I)*NORMA*CGR/RETE
         VCMAX(I)=SQRT(BAS)*UV
         MCMAX(I)=MCMAX(I)*NORMA*UM
         RCMAX(I)=RCMAX(I)   !*RETE

***********************************************************
*      GUARDAMOS LAS PARTICULAS LIGADAS  DEL HALO I
***********************************************************
         KONTA=NCAPAS(I) ! DISTANCE TO PARTICLE AT DISTANCE RJ
         KONTA2=0
         BASMAS=0.0
         DMPCLUS(I)=0

         VCMX=VX(I)
         VCMY=VY(I)
         VCMZ=VZ(I)

         BASX=0.0
         BASY=0.0
         BASZ=0.0
         INERTIA=0.0

         DIS=1000000.0
         VMAXCLUS(I)=-1.0
         EACH_PROF=KONTA/NBINS
         DO J=1,KONTA
          IF (CONTADM(J).EQ.0) THEN
           JJ=LIP(J)
           BASVEC(1)=RXPA(JJ)-CX
           BASVEC(2)=RYPA(JJ)-CY
           BASVEC(3)=RZPA(JJ)-CZ
           AADM=SQRT(BASVEC(1)**2+BASVEC(2)**2+BASVEC(3)**2)
           IF (AADM.LE.RSHELL) THEN
            KONTA2=KONTA2+1
            BASMAS=BASMAS+(MASAP(JJ)/NORMA)

            IF (MOD(J,EACH_PROF).EQ.0) THEN
             IF (J/EACH_PROF.LE.NBINS) THEN
              PROFILES(J/EACH_PROF,1,I)=AADM
              PROFILES(J/EACH_PROF,2,I)=BASMAS*NORMA*UM
             END IF
            END IF

            BASVECCM(1)=RXPA(JJ)-CMX
            BASVECCM(2)=RYPA(JJ)-CMY
            BASVECCM(3)=RZPA(JJ)-CMZ

            BASVCM(1)=U2DM(JJ)-VCMX
            BASVCM(2)=U3DM(JJ)-VCMY
            BASVCM(3)=U4DM(JJ)-VCMZ

            VVV2=BASVCM(1)**2+BASVCM(2)**2+BASVCM(3)**2

**          ANGULAR MOMENTUM
            BASX=BASX+MASAP(JJ)*(BASVECCM(2)*BASVCM(3)
     &                          -BASVECCM(3)*BASVCM(2))
            BASY=BASY+MASAP(JJ)*(BASVECCM(3)*BASVCM(1)
     &                          -BASVECCM(1)*BASVCM(3))
            BASZ=BASZ+MASAP(JJ)*(BASVECCM(1)*BASVCM(2)
     &                          -BASVECCM(2)*BASVCM(1))

**          INERTIA TENSOR
            DO JY=1,3
            DO IX=1,3
              INERTIA(IX,JY)=INERTIA(IX,JY)
     &                       +MASAP(JJ)*BASVECCM(IX)*BASVECCM(JY)
            END DO
            END DO

**        VELOCITY OF THE FASTEST PARTICLE IN THE HALO
            VKK=0.0
            VKK=SQRT(U2DM(JJ)**2+U3DM(JJ)**2+U4DM(J)**2)

            IF (VKK.GT.VMAXCLUS(I)) THEN
             VMAXCLUS(I)=VKK
            END IF

**        CLOSEST PARTICLE TO THE CENTER OF THE HALO
            IF (AADM.LT.DIS) THEN
             DIS=AADM
             IPLIP(I)=ORIPA2(JJ)
            END IF

            IF(AADM.NE.0.0) THEN
             AA=(BASVEC(1)/AADM)*BASVCM(1)+
     &          (BASVEC(2)/AADM)*BASVCM(2)+
     &          (BASVEC(3)/AADM)*BASVCM(3)
             VR=VR+AA*MASAP(JJ)
            END IF
           END IF    !AADM.LT.RSHELL
          END IF     !CONTADM
         END DO      !KONTA

*******************************************************
*      SAVING MASSES, RADII, PROFILES AND SHAPES...
*******************************************************
         DMPCLUS(I)=KONTA2
         IF (KONTA2.LT.NUMPARTBAS) REALCLUS(I)=0
         !MASA(I)=BASMAS*NORMA*UM
         !RADIO(I)=RSHELL
         ANGULARM(1,I)=BASX*UV / (BASMAS*NORMA)
         ANGULARM(2,I)=BASY*UV / (BASMAS*NORMA)
         ANGULARM(3,I)=BASZ*UV / (BASMAS*NORMA)
         MEAN_VR(I)=VR/(BASMAS*NORMA)

         INERTIA(1:3,1:3)=INERTIA(1:3,1:3)/(BASMAS*NORMA)
         INERTIA_TENSOR(1,I)=INERTIA(1,1)
         INERTIA_TENSOR(2,I)=INERTIA(1,2)
         INERTIA_TENSOR(3,I)=INERTIA(1,3)
         INERTIA_TENSOR(4,I)=INERTIA(2,2)
         INERTIA_TENSOR(5,I)=INERTIA(2,3)
         INERTIA_TENSOR(6,I)=INERTIA(3,3)

         BASEIGENVAL(1:3)=0.0
         IF (DMPCLUS(I).GE.NUMPARTBAS) THEN
          CALL JACOBI(INERTIA,DIMEN,BASEIGENVAL,NROT)
          CALL SORT(BASEIGENVAL,DIMEN,DIMEN)
         END IF

         DO II=1,DIMEN
          EIGENVAL(II,I)=BASEIGENVAL(II)
          EIGENVAL(II,I)=SQRT(EIGENVAL(II,I))
         END DO

C         WRITE(*,*)
C         WRITE(*,*)I,CLUSRX(I),CLUSRY(I),CLUSRZ(I),
C     &         MASA(I),RADIO(I),DMPCLUS(I),
C     &         REALCLUS(I), LEVHAL(I),(INERTIA_TENSOR(IX,I),IX=1,6),
C     &         EIGENVAL(1,I),EIGENVAL(2,I),EIGENVAL(3,I),MEAN_VR(I)*UV,
C     &         CONCENTRA(I),(ANGULARM(J,I),J=1,3),
C     &         VCMAX(I),MCMAX(I),RCMAX(I),
C     &         R200M(I),M200M(I),R200C(I),M200C(I),
C     &         R500M(I),M500M(I),R500C(I),M500C(I),
C     &         R2500M(I),M2500M(I),R2500C(I),M2500C(I),
C     &         VX(I)*UV,VY(I)*UV,VZ(I)*UV

******************************************
************ FIN HALO I ******************
******************************************
        END IF ! (ELSE of KONTA2.LT.NUMPARTBAS) (poor haloes after unbinding)
       END IF ! (realclus(i).ne.0)

*****************
       END DO   !I=IP,IP2
****************

       WRITE(*,*) '... Finally, haloes found at level:',IR,
     &            COUNT(REALCLUS(LOWH1:LOWH2).NE.0)
       WRITE(*,*)

       RETURN
       END
