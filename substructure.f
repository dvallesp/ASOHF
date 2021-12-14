********************************************************************
       SUBROUTINE SEARCH_SUBSTRUCTURE_GRID(IR,NL,NX,NY,NZ,NPATCH,
     &             PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &             PATCHRY,PATCHRZ,PARE,NCLUS,MASA,RADIO,CLUSRX,CLUSRY,
     &             CLUSRZ,REALCLUS,LEVHAL,NHALLEV,BOUND,CONTRASTEC,RODO,
     &             SOLAP,VECINO,NVECI,CR0AMR,CR0AMR11,PATCHCLUS,
     &             VOL_SOLAP_LOW,CLUSRXCM,CLUSRYCM,CLUSRZCM,RSUB,MSUB,
     &             SUBS_LEV,UM)
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

       REAL*4, ALLOCATABLE::DDD(:)
       INTEGER, ALLOCATABLE::DDDX(:),DDDY(:),DDDZ(:),DDDP(:)

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
!$OMP+                   RADIO,RX,RY,RZ,CONTA1,NCLUS),
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
         BAS=RADIO(II)
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
!$OMP+                   RADIO,RX,RY,RZ,CONTA1,NCLUS),
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
!$OMP  PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,U11,CONTRASTEC,LOW1,
!$OMP+                    LOW2,CONTA1,BORAMR),
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
         IF (U11(IX,JY,KZ,I).GE.CONTRASTEC.AND.
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
       WRITE(*,*) 'Clean cells above virial contrast',KK_ENTERO

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
         IF (U11(IX,JY,KZ,I).GE.CONTRASTEC.AND.
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

        IF (IR.NE.NL) THEN
         XCEN=RX(IX,IPATCH)
         YCEN=RY(JY,IPATCH)
         ZCEN=RZ(KZ,IPATCH)
         CALL RECENTER_DENSITY_PEAK_SUBS(XCEN,YCEN,ZCEN,IR,U1,U11,
     &                               PATCHRX,PATCHRY,PATCHRZ,PATCHNX,
     &                               PATCHNY,PATCHNZ,NPATCH,LADO0,NL,
     &                               IX,JY,KZ,IPATCH)
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

*        Assert we have not yet identified this same halo (dupplicated
*         due to recentering)
        FLAG_DUP=0
        DO II=NCLUS_INI+1,NCLUS
         IF ((XCEN-CLUSRX(II))**2+(YCEN-CLUSRY(II))**2+
     &       (ZCEN-CLUSRZ(II))**2.LT.RSUB(II)**2) THEN
          FLAG_DUP=1
          EXIT
         END IF
        END DO
        IF (FLAG_DUP.EQ.1) CYCLE

        CALL FIND_HOST_HALO(XCEN,YCEN,ZCEN,DXPA,IR,SUBS_LEV,REALCLUS,
     &                      CLUSRX,CLUSRY,CLUSRZ,RADIO,IHOSTHALO)

        IF (IHOSTHALO.EQ.-1) THEN
         write(*,*) 'not found parent',xcen,ycen,zcen
         CYCLE
        END IF

        BASXX=CLUSRX(IHOSTHALO)
        BASYY=CLUSRY(IHOSTHALO)
        BASZZ=CLUSRZ(IHOSTHALO)
        MHOST=MASA(IHOSTHALO)
        DISTHOST=SQRT((BASXX-XCEN)**2+(BASYY-YCEN)**2+(BASZZ-ZCEN)**2)

        IF (DISTHOST.LT.2.0*DXPA) CYCLE

        NCLUS=NCLUS+1
        REALCLUS(NCLUS)=IHOSTHALO
        LEVHAL(NCLUS)=IR
        NHALLEV(IR)=NHALLEV(IR)+1
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
                   (RADZ(K)-ZCEN)**2)

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

C       WRITE(*,*)

       RETURN
       END

********************************************************************
       SUBROUTINE RECENTER_DENSITY_PEAK_SUBS(BASX,BASY,BASZ,HLEV,U1,U11,
     &                                       PATCHRX,PATCHRY,PATCHRZ,
     &                                       PATCHNX,PATCHNY,PATCHNZ,
     &                                       NPATCH,LADO0,NL,
     &                                       PIX,PJY,PKZ,PIPATCH)
********************************************************************
*      Refines the location of a density peak using finer AMR levels
*      Case for substructure (impose it is a relative maximum!)
********************************************************************
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       REAL BASX,BASY,BASZ
       INTEGER HLEV
       REAL U1(NMAX,NMAY,NMAZ),U11(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER NPATCH(0:NLEVELS)
       REAL LADO0
       INTEGER NL
       INTEGER PIX,PJY,PKZ,PIPATCH

       REAL*4  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
       COMMON /GRID/ RADX,RADY,RADZ

       REAL*4  RX(0:NAMRX+1,NPALEV),RY(0:NAMRX+1,NPALEV),
     &         RZ(0:NAMRX+1,NPALEV)
       COMMON /GRIDAMR/ RX,RY,RZ

       REAL*4 DX,DY,DZ,H2
       COMMON /ESPACIADO/ DX,DY,DZ

       REAL X1,X2,X3,X4,Y1,Y2,Y3,Y4,Z1,Z2,Z3,Z4,DXPA,DYPA,DZPA
       REAL DXP,DYP,DZP,DXM,DYM,DZM,CURRENTMAX
       REAL,ALLOCATABLE::UBAS(:,:,:)
       INTEGER IX,JY,KZ,II,JJ,KK,IR,LOW1,LOW2,FLAG_FOUND,I,I1,I2,J1,J2
       INTEGER K1,K2,INMAX(3),BOR

       BOR=2
       ALLOCATE(UBAS(2+2*BOR,2+2*BOR,2+2*BOR))

c       WRITE(*,*) hlev,BASX,BASY,BASZ

C       WRITE(*,*) 'RECENTER',HLEV,BASX,BASY,BASZ

       loop_levels: DO IR=HLEV+1,NL
        DXPA=DX/2.0**IR
        DYPA=DY/2.0**IR
        DZPA=DZ/2.0**IR

        X1=BASX-(1.0+BOR)*DXPA
        X2=BASX+(1.0+BOR)*DXPA
        Y1=BASY-(1.0+BOR)*DYPA
        Y2=BASY+(1.0+BOR)*DYPA
        Z1=BASZ-(1.0+BOR)*DZPA
        Z2=BASZ+(1.0+BOR)*DZPA

        FLAG_FOUND=0

        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))

        loop_patches: DO I=LOW1,LOW2
         X3=PATCHRX(I)-DXPA
         X4=X3+PATCHNX(I)*DXPA
         Y3=PATCHRY(I)-DYPA
         Y4=Y3+PATCHNY(I)*DYPA
         Z3=PATCHRZ(I)-DZPA
         Z4=Z3+PATCHNZ(I)*DZPA

         IF (X3.LT.X1.AND.X2.LT.X4.AND.   ! what we are looking is
     &       Y3.LT.Y1.AND.Y2.LT.Y4.AND.   ! completely enclosed in the
     &       Z3.LT.Z1.AND.Z2.LT.Z4) THEN  ! patch
          FLAG_FOUND=1
          EXIT loop_patches
         END IF

        END DO loop_patches

        IF (FLAG_FOUND.EQ.0) EXIT loop_levels

        I1=INT((X1-X3)/DXPA+0.5)+1
        J1=INT((Y1-Y3)/DYPA+0.5)+1
        K1=INT((Z1-Z3)/DZPA+0.5)+1
        I2=I1+1+2*BOR
        J2=J1+1+2*BOR
        K2=K1+1+2*BOR
c        WRITE(*,*) '--------------'
c        WRITE(*,*) IR
c        WRITE(*,*) X1,X2,X3,X4,I1,I2
c        WRITE(*,*) Y1,Y2,Y3,Y4,J1,J2
c        WRITE(*,*) Z1,Z2,Z3,Z4,K1,K2

        UBAS=U11(I1:I2,J1:J2,K1:K2,I)

        CURRENTMAX=-100000.0
        DO KK=3,4
        DO JJ=3,4
        DO II=3,4
         DXP=UBAS(II+1,JJ,KK)-UBAS(II,JJ,KK)
         DXM=UBAS(II,JJ,KK)-UBAS(II-1,JJ,KK)
         DYP=UBAS(II,JJ+1,KK)-UBAS(II,JJ,KK)
         DYM=UBAS(II,JJ,KK)-UBAS(II,JJ-1,KK)
         DZP=UBAS(II,JJ,KK+1)-UBAS(II,JJ,KK)
         DZM=UBAS(II,JJ,KK)-UBAS(II,JJ,KK-1)
         IF (DXP.LT.0.0.AND.DYP.LT.0.0.AND.DZP.LT.0.0.AND.
     &       DXM.GT.0.0.AND.DYM.GT.0.0.AND.DZM.GT.0.0) THEN
          IF (UBAS(II,JJ,KK).GT.CURRENTMAX) THEN
           CURRENTMAX=UBAS(II,JJ,KK)
           INMAX(1)=II
           INMAX(2)=JJ
           INMAX(3)=KK
          END IF
         END IF
        END DO
        END DO
        END DO

        IF (CURRENTMAX.LT.0.0) THEN
         DO KK=2,5
         DO JJ=2,5
         DO II=2,5
          DXP=UBAS(II+1,JJ,KK)-UBAS(II,JJ,KK)
          DXM=UBAS(II,JJ,KK)-UBAS(II-1,JJ,KK)
          DYP=UBAS(II,JJ+1,KK)-UBAS(II,JJ,KK)
          DYM=UBAS(II,JJ,KK)-UBAS(II,JJ-1,KK)
          DZP=UBAS(II,JJ,KK+1)-UBAS(II,JJ,KK)
          DZM=UBAS(II,JJ,KK)-UBAS(II,JJ,KK-1)
          IF (DXP.LT.0.0.AND.DYP.LT.0.0.AND.DZP.LT.0.0.AND.
     &        DXM.GT.0.0.AND.DYM.GT.0.0.AND.DZM.GT.0.0) THEN
           IF (UBAS(II,JJ,KK).GT.CURRENTMAX) THEN
            CURRENTMAX=UBAS(II,JJ,KK)
            INMAX(1)=II
            INMAX(2)=JJ
            INMAX(3)=KK
           END IF
          END IF
         END DO
         END DO
         END DO
        END IF

        IF (CURRENTMAX.LT.0.0) EXIT loop_levels

        IX=I1+INMAX(1)-1
        JY=J1+INMAX(2)-1
        KZ=K1+INMAX(3)-1

        BASX=RX(IX,I)
        BASY=RY(JY,I)
        BASZ=RZ(KZ,I)

        PIX=IX
        PJY=JY
        PKZ=KZ
        PIPATCH=I
c       WRITE(*,*) IR,BASX,BASY,BASZ,U11(IX,JY,KZ,I),I,inmax
c       write(*,*)
C       WRITE(*,*) 'RECENTER',IR,BASX,BASY,BASZ,INMAX(1:3)
       END DO loop_levels
c       WRITE(*,*) '-----'

       RETURN
       END

********************************************************************
       SUBROUTINE FIND_HOST_HALO(BASX,BASY,BASZ,DXPA,IR,SUBS_LEV,
     &                           REALCLUS,CLUSRX,CLUSRY,CLUSRZ,RADIO,
     &                           IHOSTHALO)
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
       REAL RADIO(MAXNCLUS)
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
         BAS=RADIO(II)

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
        WRITE(*,*) 'Not found progenitor'
        IHOSTHALO=-1
       END IF

       RETURN
       END
