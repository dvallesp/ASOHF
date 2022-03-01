********************************************************************
       SUBROUTINE BRYAN_NORMAN_98(CONTRASTEC,OMEGAZ,OMEGA0,ZETA)
****************************************************************
*      VIRIAL CONTRAST (in terms of the BACKGROUND MATTER density)
*      (Bryan & Norman ApJ, 1998)
*      Delta_vir,c = 18*pi^2 + 82 x - 39 x^2; x=Omega_m(z)-1
****************************************************************
        IMPLICIT NONE

        REAL CONTRASTEC,OMEGAZ ! intent:out
        REAL OMEGA0,ZETA ! intent in
        REAL CONTRASTEX,PI,OMEGALAMBDA0,BAS

        PI=DACOS(-1.D0)
        OMEGALAMBDA0=1.0-OMEGA0
        BAS=OMEGA0*(1+ZETA)**3
        OMEGAZ=BAS/(BAS+OMEGALAMBDA0)
        CONTRASTEX=OMEGAZ - 1.0
        CONTRASTEC= 18.0*PI**2 + 82.0*CONTRASTEX - 39.0*CONTRASTEX**2
        CONTRASTEC=CONTRASTEC/OMEGAZ

        RETURN

       END

********************************************************************
       SUBROUTINE OVERLAPPING(IR,NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                        PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                        PATCHRY,PATCHRZ,PARE,NCLUS,MASA,RADIO,
     &                        CLUSRX,CLUSRY,CLUSRZ,REALCLUS,LEVHAL,
     &                        NHALLEV,BOUND,CONTRASTEC,RODO,SOLAP,
     &                        CR0AMR,CR0AMR11,VOL_SOLAP_LOW)
********************************************************************
*      Detect and correct overlaps on the cluster catalogue
********************************************************************

       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

*      I/O DATA
       INTEGER IR,NL,NX,NY,NZ
       INTEGER NPATCH(0:NLEVELS)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
       REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
       INTEGER PARE(NPALEV)
       INTEGER NCLUS
       REAL MASA(MAXNCLUS),RADIO(MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
       INTEGER REALCLUS(MAXNCLUS),LEVHAL(MAXNCLUS)
       INTEGER NHALLEV(0:NLEVELS)
       REAL BOUND,CONTRASTEC,RODO
       INTEGER SOLAP(NAMRX,NAMRY,NAMRZ,NPALEV)
       INTEGER CR0AMR(NMAX,NMAY,NMAZ)
       INTEGER CR0AMR11(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL VOL_SOLAP_LOW

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
       REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
       COMMON /VARIA/ U1,U11

       REAL*4 ACHE,T0,RE0
       COMMON /DOS/ ACHE,T0,RE0

*      LOCAL VARIABLES
       INTEGER,ALLOCATABLE::SOLAPA(:,:),NSOLAP(:),NMAXSOLAP
       INTEGER IRR,I,J,K,IX,JY,KZ,IPATCH,LOWH1,LOWH2,LOWH3,BASINT
       INTEGER IMAXCLUS,IMINCLUS,CONTA,FLAG_ITER,NUM_ITERS,BASINT2
       REAL BAS,BASX,BASY,BASZ,R1,R2,DIST,X1,Y1,Z1,X2,Y2,Z2,VOL1,VOL2
       REAL VINT,XI,PI,SOLAP_LOWER_THR,M1,M2

       SOLAP_LOWER_THR=VOL_SOLAP_LOW

       PI=DACOS(-1.D0)
       NMAXSOLAP=200

       WRITE(*,*) '== HALOES OVERLAPPING IN IR =', IR
       WRITE(*,*) 'NCLUS,NHALLEV(IR)=', NCLUS, NHALLEV(IR)

       LOWH1=1 !SUM(NHALLEV(0:IR-1))+1
       LOWH2=SUM(NHALLEV(0:IR))
       LOWH3=SUM(NHALLEV(0:IR-1))+1

       ALLOCATE(NSOLAP(LOWH1:LOWH2), SOLAPA(LOWH1:LOWH2,NMAXSOLAP))

!$OMP PARALLEL DO SHARED(NSOLAP,SOLAPA,NMAXSOLAP,LOWH1,LOWH2),
!$OMP+            PRIVATE(I,J),
!$OMP+            DEFAULT(NONE)
       DO I=LOWH1,LOWH2
        NSOLAP(I)=0
        DO J=1,NMAXSOLAP
         SOLAPA(I,J)=0
        END DO
       END DO

       NUM_ITERS=0
       FLAG_ITER=1

       CONTA=0

*       1. Look for clusters completely included in other ones, and delete
*        the smallest ones.
       DO I=LOWH1,LOWH2
        BASINT=REALCLUS(I)
        IF (BASINT.NE.0) THEN
         R1=RADIO(I)
         X1=CLUSRX(I)
         Y1=CLUSRY(I)
         Z1=CLUSRZ(I)
         DO J=MAX(I+1,LOWH3),LOWH2
          BASINT2=REALCLUS(J)
          IF (BASINT2.NE.0) THEN
           R2=RADIO(J)
           X2=CLUSRX(J)
           Y2=CLUSRY(J)
           Z2=CLUSRZ(J)
           DIST=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
           IF (ABS(R1-R2).GE.DIST) THEN
            ! A sphere is completely included in the other one.
            ! We shall remove the smallest one.
            IMAXCLUS=I
            IMINCLUS=J
            IF (R2.GT.R1) THEN
             IMAXCLUS=J
             IMINCLUS=I
            END IF !(R2.GT.R1)

            REALCLUS(IMINCLUS)=0 !The smallest halo is removed
            CONTA=CONTA+1
           END IF !(ABS(R1-R2).GE.DIST)
          END IF !(BASINT2.NE.0)
         END DO !J=I+1,LOWH2
        END IF !(BASINT.NE.0) THEN
       END DO !I=LOWH1,LOWH2

*       2. Look for clusters which overlap a large amount of their
*        volume, and merge them.
       DO I=LOWH1,LOWH2
        BASINT=REALCLUS(I)
        IF (BASINT.NE.0) THEN
         R1=RADIO(I)
         X1=CLUSRX(I)
         Y1=CLUSRY(I)
         Z1=CLUSRZ(I)
         M1=MASA(I)
         DO J=MAX(I+1,LOWH3),LOWH2
          BASINT2=REALCLUS(J)
          IF (BASINT2.NE.0) THEN
           R2=RADIO(J)
           X2=CLUSRX(J)
           Y2=CLUSRY(J)
           Z2=CLUSRZ(J)
           M2=MASA(J)

           ! Compute the volume overlap
           DIST=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
           BAS=DIST/(R1+R2)
           IF (BAS.LT.1.0) THEN
            XI=BAS
            VINT=(PI/12.0) * (R1+R2)**3 * ((1-XI)**2/XI) *
     &               (XI**2 + 2*XI - 3*(R1-R2)**2/(R1+R2)**2)
            VINT=VINT/(4*PI/3*MIN(R1,R2)**3) ! fraction of the smallest cluster volume overlapped

            IF (VINT.GT.SOLAP_LOWER_THR) THEN
             IMAXCLUS=I
             IMINCLUS=J
             IF (R2.GT.1.5*R1) THEN
              IMAXCLUS=J
              IMINCLUS=I
             END IF !(R2.GT.R1)

             REALCLUS(IMINCLUS)=0 !The smallest halo is removed
             CONTA=CONTA+1

C              RADIO(IMAXCLUS)=(R1**3+R2**3)**(1.0/3.0)
C              CLUSRX(IMAXCLUS)=(M1*X1+M2*X2)/(M1+M2)
C              CLUSRY(IMAXCLUS)=(M1*Y1+M2*Y2)/(M1+M2)
C              CLUSRZ(IMAXCLUS)=(M1*Z1+M2*Z2)/(M1+M2)

C              write(*,*) '--------------'
C              write(*,*) i,j,dist,bas
C              write(*,*) x1,y1,z1,r1,m1
C              write(*,*) x2,y2,z2,r2,m2
C              write(*,*) vint
C              write(*,*) '-->',CLUSRX(IMAXCLUS),CLUSRY(IMAXCLUS),
C     &                         CLUSRZ(IMAXCLUS),RADIO(IMAXCLUS)

            END IF !(VINT.GT.SOLAP_LOWER_THR)
           END IF !(BAS.LT.1.0)
          END IF !(BASINT2.NE.0)
         END DO !J=I+1,LOWH2
        END IF !(BASINT.NE.0) THEN
       END DO !I=LOWH1,LOWH2


       BASINT=COUNT(REALCLUS(LOWH1:LOWH2).EQ.0)
       WRITE(*,*) 'REMOVED HALOS_0----->', BASINT
       BASINT=COUNT(REALCLUS(LOWH1:LOWH2).NE.0)
       WRITE(*,*) 'POSSIBLE HALOS_0----->', BASINT
       BASINT=COUNT(REALCLUS(LOWH1:LOWH2).EQ.-1)
       BASINT2=COUNT(REALCLUS(LOWH1:LOWH2).GT.0)
c       WRITE(*,*) '--> Of which free,substructure:',BASINT,BASINT2
c       WRITE(*,*)

       RETURN
       END

********************************************************************
       SUBROUTINE HALOFIND_GRID(NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                          PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                          PATCHRY,PATCHRZ,PARE,NCLUS,MASA,RADIO,
     &                          CLUSRX,CLUSRY,CLUSRZ,REALCLUS,LEVHAL,
     &                          NHALLEV,BOUND,CONTRASTEC,RODO,SOLAP,
     &                          CR0AMR,CR0AMR11,PATCHCLUS,VOL_SOLAP_LOW,
     &                          CLUSRXCM,CLUSRYCM,CLUSRZCM)
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
       INTEGER CR0AMR(NMAX,NMAY,NMAZ)
       INTEGER CR0AMR11(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL VOL_SOLAP_LOW

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
       REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
       COMMON /VARIA/ U1,U11

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
       INTEGER K1,BORAMR,IRR,BASINT,FLAG_ITER,BOR,FLAG,FLAG_DUP
       REAL PRUEBAX,PRUEBAY,PRUEBAZ,RMIN,BASMASS_SHELL,BASMASS,DELTA
       REAL ESP,ESP_LOG,BAS,KK_REAL,RSHELL,R_INT,R_EXT,RANT,LADO0
       REAL BASDELTA,AA,PI,VOLCELL,BASX,BASY,BASZ,BASVOL,BOUNDIR2
       REAL X1,X2,Y1,Y2,Z1,Z2,DXPA,DYPA,DZPA,BASXX,BASYY,BASZZ
       REAL XCEN,YCEN,ZCEN,BOUNDIR,X3,Y3,Z3,X4,Y4,Z4,MINDERIV
       REAL VECDENS(1000),VECRAD(1000),DERIVATIVE(1000),BASVOL_SHELL
       REAL DXPAPA,DYPAPA,DZPAPA,CONTRASTECPEAK

*      Additional variables to perform reductions on sums
c       REAL R_BASVOL,R_BASMASS_SHELL,R_BASX,R_BASY,R_BASZ,R_BASDELTA
c       INTEGER R_II

       REAL*4, ALLOCATABLE::DDD(:)
       INTEGER, ALLOCATABLE::DDDX(:),DDDY(:),DDDZ(:),DDDP(:)

       !This is to be able to find peaks which are not dense enough over
       ! the grid due to the smooth density interpolation
       CONTRASTECPEAK=CONTRASTEC/6.0

**************************************************************
*      NIVEL BASE!!
**************************************************************

       PI=DACOS(-1.D0)
       LADO0=NX*DX
       NCLUS=0

       IR=0
       ESP=0.2*DX
       ESP_LOG=1.05
       BOR=1
       WRITE(*,*) '<-------- BASE GRID -------->'
       WRITE(*,*) '--> IR,DX:',IR,DX

!$OMP PARALLEL DO SHARED(NX,NY,NZ,CONTA,BOR),
!$OMP+            PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO KZ=1+BOR,NZ-BOR
       DO JY=1+BOR,NY-BOR
       DO IX=1+BOR,NX-BOR
        CONTA(IX,JY,KZ)=1
       END DO
       END DO
       END DO

       KK_ENTERO=0
!$OMP PARALLEL DO SHARED(NX,NY,NZ,U1,CONTRASTECPEAK,BOR),
!$OMP+            PRIVATE(IX,JY,KZ,BASX,BASY,BASZ,BASXX,BASYY,BASZZ),
!$OMP+            REDUCTION(+:KK_ENTERO),
!$OMP+            DEFAULT(NONE)
       DO KZ=1+BOR,NZ-BOR
       DO JY=1+BOR,NY-BOR
       DO IX=1+BOR,NX-BOR
        IF (U1(IX,JY,KZ).GE.CONTRASTECPEAK) THEN
         BASX =U1(IX+1,JY,KZ)-U1(IX,JY,KZ)
         IF (BASX.LT.0.0) THEN
         BASY =U1(IX,JY+1,KZ)-U1(IX,JY,KZ)
         IF (BASY.LT.0.0) THEN
         BASZ =U1(IX,JY,KZ+1)-U1(IX,JY,KZ)
         IF (BASZ.LT.0.0) THEN
         BASXX=U1(IX,JY,KZ)  -U1(IX-1,JY,KZ)
         IF (BASXX.GT.0.0) THEN
         BASYY=U1(IX,JY,KZ)  -U1(IX,JY-1,KZ)
         IF (BASYY.GT.0.0) THEN
         BASZZ=U1(IX,JY,KZ)  -U1(IX,JY,KZ-1)
         IF (BASZZ.GT.0.0) THEN ! then it's a local maximum
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
       WRITE(*,*) 'Clean cells above virial contrast',KK_ENTERO
*      Ordenamos todas las celdas de los vecinos
       ALLOCATE(DDD(KK_ENTERO))
       ALLOCATE(DDDX(KK_ENTERO))
       ALLOCATE(DDDY(KK_ENTERO))
       ALLOCATE(DDDZ(KK_ENTERO))

!$OMP PARALLEL DO SHARED(DDD,DDDX,DDDY,DDDZ,KK_ENTERO),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
       DO I=1,KK_ENTERO
        DDD(I)=0.0
        DDDX(I)=0
        DDDY(I)=0
        DDDZ(I)=0
       END DO

       II=0
       DO KZ=1+BOR,NZ-BOR
       DO JY=1+BOR,NY-BOR
       DO IX=1+BOR,NX-BOR
        BAS=U1(IX,JY,KZ)
        IF (BAS.GE.CONTRASTECPEAK) THEN
         BASX =U1(IX+1,JY,KZ)-U1(IX,JY,KZ)
         IF (BASX.LT.0.0) THEN
         BASY =U1(IX,JY+1,KZ)-U1(IX,JY,KZ)
         IF (BASY.LT.0.0) THEN
         BASZ =U1(IX,JY,KZ+1)-U1(IX,JY,KZ)
         IF (BASZ.LT.0.0) THEN
         BASXX=U1(IX,JY,KZ)  -U1(IX-1,JY,KZ)
         IF (BASXX.GT.0.0) THEN
         BASYY=U1(IX,JY,KZ)  -U1(IX,JY-1,KZ)
         IF (BASYY.GT.0.0) THEN
         BASZZ=U1(IX,JY,KZ)  -U1(IX,JY,KZ-1)
         IF (BASZZ.GT.0.0) THEN ! then it's a local maximum
          II=II+1
          DDD(II)=BAS
          DDDX(II)=IX
          DDDY(II)=JY
          DDDZ(II)=KZ
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

       IF (II.NE.KK_ENTERO) THEN
         WRITE(*,*) 'WARNING: bad allocation of DDD', II, KK_ENTERO
         STOP
       END IF

       CALL SORT_CELLS(KK_ENTERO,DDD,DDDX,DDDY,DDDZ)

       NV_GOOD=KK_ENTERO
       KK_ENTERO=0

*      Go through all the center candidates
       DO L1=1,NV_GOOD
        ICEN(1)=DDDX(L1)
        ICEN(2)=DDDY(L1)
        ICEN(3)=DDDZ(L1)

        KK_ENTERO=CONTA(ICEN(1),ICEN(2),ICEN(3))
        IF(KK_ENTERO.EQ.1) THEN ! this means this peak is not inside a halo yet
c         WRITE(*,*) U1(ICEN(1),ICEN(2),ICEN(3))
         IF(NCLUS.GT.MAXNCLUS) THEN
          WRITE(*,*) 'WARNING: NCLUS>MAXNCLUS!!!',NCLUS,MAXNCLUS
          STOP
         END IF

         XCEN=RADX(ICEN(1))
         YCEN=RADY(ICEN(2))
         ZCEN=RADZ(ICEN(3))

         IF (NL.GT.0) THEN
          IPATCH=0
          CALL RECENTER_DENSITY_PEAK(XCEN,YCEN,ZCEN,0,U1,U11,PATCHRX,
     &                               PATCHRY,PATCHRZ,PATCHNX,PATCHNY,
     &                               PATCHNZ,NPATCH,LADO0,NL,IX,JY,KZ,
     &                               IPATCH)
*         Assert we have not yet identified this same halo (dupplicated
*          due to recentering)
          FLAG_DUP=0
          DO II=1,NCLUS
           IF ((XCEN-CLUSRX(II))**2+(YCEN-CLUSRY(II))**2+
     &         (ZCEN-CLUSRZ(II))**2.LT.RADIO(II)**2) THEN
            FLAG_DUP=1
            EXIT
           END IF
          END DO
          IF (FLAG_DUP.EQ.1) CYCLE
         END IF

         NCLUS=NCLUS+1
         REALCLUS(NCLUS)=-1
         LEVHAL(NCLUS)=IR
         NHALLEV(IR)=NHALLEV(IR)+1
         PATCHCLUS(NCLUS)=IPATCH

         CLUSRX(NCLUS)=XCEN
         CLUSRY(NCLUS)=YCEN
         CLUSRZ(NCLUS)=ZCEN

*        tentative reach of the halo ---> build a mini box around it
*        (by excess, set as a parameter in asohf.dat)

*        Now, we extend the cluster radially from its center
         BASX=0.0
         BASY=0.0
         BASZ=0.0
         BASDELTA=0.0

         BASMASS_SHELL=0.0
         BASMASS=0.0   !TOTAL MASS OF THE CLUSTER
         DELTA=0.0    !TOTAL CONTRAST OF THE CLUSTER
         BASVOL=0.0     !TOTAL VOLUME OF THE CLUSTER (sphere)

         R_INT=0.0
         R_EXT=0.8*DX

*        increase the radius until density falls below the virial value
         DELTA=10.0*CONTRASTEC*ROTE ! this is to ensure we enter the loop
         ITER_GROW=0
         DO WHILE(DELTA.GT.CONTRASTEC*ROTE)
          ITER_GROW=ITER_GROW+1

          IF (ITER_GROW.GT.1) THEN
           IF (R_EXT.LE.BOUND) THEN
            R_INT=R_EXT
            R_EXT=MAX(R_EXT+ESP, R_EXT*ESP_LOG)
           ELSE
            WRITE(*,*) 'WARNING: growing not converged', r_int, r_ext,
     &                                                   bound
            STOP
           END IF
          END IF

          BASX=XCEN-R_EXT
          NX1=INT(((BASX-RADX(1))/DX)+0.5)+1
          IF (NX1.LT.1) NX1=1

          BASX=XCEN+R_EXT
          NX2=INT(((BASX-RADX(1))/DX)+0.5)+1
          IF (NX2.GT.NX) NX2=NX

          BASY=YCEN-R_EXT
          NY1=INT(((BASY-RADY(1))/DY)+0.5)+1
          IF (NY1.LT.1) NY1=1

          BASY=YCEN+R_EXT
          NY2=INT(((BASY-RADY(1))/DY)+0.5)+1
          IF (NY2.GT.NY) NY2=NY

          BASZ=ZCEN-R_EXT
          NZ1=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
          IF (NZ1.LT.1) NZ1=1

          BASZ=ZCEN+R_EXT
          NZ2=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
          IF (NZ2.GT.NZ) NZ2=NZ

          VOLCELL=DX*DY*DZ
          BASMASS_SHELL=0.0
          II=0
!$OMP PARALLEL DO SHARED(NZ1,NZ2,NY1,NY2,NX1,NX2,RADX,RADY,RADZ,XCEN,
!$OMP+                   YCEN,ZCEN,R_INT,R_EXT,CONTA,VOLCELL,U1),
!$OMP+            PRIVATE(I,J,K,AA,BAS),
!$OMP+            REDUCTION(+:BASVOL,BASMASS_SHELL,BASX,BASY,BASZ,
!$OMP+                        BASDELTA,II),
!$OMP+            DEFAULT(NONE)
          DO K=NZ1,NZ2
          DO J=NY1,NY2
          DO I=NX1,NX2
           AA=SQRT((RADX(I)-XCEN)**2 + (RADY(J)-YCEN)**2 +
     &             (RADZ(K)-ZCEN)**2)

           IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
            CONTA(I,J,K)=0 ! do not try to find an additional halo here (at this level)

            BASVOL=BASVOL+VOLCELL

            BAS=U1(I,J,K)*VOLCELL !U1 is not density contrast, but 1+delta = rho/rho_B!!!
            BASMASS_SHELL=BASMASS_SHELL+BAS

            BASX=BASX+RADX(I)*BAS
            BASY=BASY+RADY(J)*BAS
            BASZ=BASZ+RADZ(K)*BAS
            BASDELTA=BASDELTA+BAS

            II=II+1
           END IF

          END DO
          END DO
          END DO

          BASMASS=BASMASS+BASMASS_SHELL*RODO*RE0**3
          IF (BASVOL.GT.0.0) DELTA=BASMASS/(BASVOL*RETE**3)
c         WRITE(*,*) DELTA/ROTE, II
         END DO   ! do while (DELTA)

         RADIO(NCLUS)=R_EXT
         MASA(NCLUS)=BASMASS

         CLUSRXCM(NCLUS)=BASX/BASDELTA
         CLUSRYCM(NCLUS)=BASY/BASDELTA
         CLUSRZCM(NCLUS)=BASZ/BASDELTA

C         WRITE(*,*) '*** CLUSTER',NCLUS,'FINAL ********'
C         WRITE(*,*) CLUSRX(NCLUS),CLUSRY(NCLUS),CLUSRZ(NCLUS),
C     &              RADIO(NCLUS),MASA(NCLUS)*9.1717E18,IR
C         write(*,*) 'offset',SQRT((CLUSRX(NCLUS)-CLUSRXCM(NCLUS))**2 +
C     &                            (CLUSRY(NCLUS)-CLUSRYCM(NCLUS))**2 +
C     &                            (CLUSRZ(NCLUS)-CLUSRZCM(NCLUS))**2)

        END IF ! KK_ENTERO.EQ.0
       END DO ! I=1,NV_GOOD

       DEALLOCATE(DDD,DDDX,DDDY,DDDZ)

       WRITE(*,*) 'At level',IR,', num. haloes:', NHALLEV(IR)

       CALL OVERLAPPING(IR,NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                  PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &                  PARE,NCLUS,MASA,RADIO,CLUSRX,CLUSRY,CLUSRZ,
     &                  REALCLUS,LEVHAL,NHALLEV,BOUND,CONTRASTEC,RODO,
     &                  SOLAP,CR0AMR,CR0AMR11,
     &                  VOL_SOLAP_LOW)

       WRITE(*,*) 'End of base level', 0



       IF (NL.GT.0) THEN
        WRITE(*,*)
        WRITE(*,*) '<-------- Now proceeding with the',NL,
     &             'AMR levels -------->'
*      Find the l=1 cells covered by l=0 haloes (they will
*      potentially host substructures)
!$OMP PARALLEL DO SHARED(NPATCH,PATCHNX,PATCHNY,PATCHNZ,CONTA1),
!$OMP+            PRIVATE(I,N1,N2,N3,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
         DO I=1,NPATCH(1)
          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)
          DO KZ=1,N3
          DO JY=1,N2
          DO IX=1,N1
           CONTA1(IX,JY,KZ,I)=1
          END DO
          END DO
          END DO
         END DO

!$OMP PARALLEL DO SHARED(NPATCH,PATCHNX,PATCHNY,PATCHNZ,PATCHRX,PATCHRY,
!$OMP+                   PATCHRZ,DX,DY,DZ,CLUSRX,CLUSRY,CLUSRZ,RADIO,
!$OMP+                   RX,RY,RZ,CONTA1,NCLUS),
!$OMP+            PRIVATE(I,N1,N2,N3,X1,X2,Y1,Y2,Z1,Z2,BASX,BASY,BASZ,
!$OMP+                    BAS,IX,JY,KZ,AA,X3,Y3,Z3,X4,Y4,Z4),
!$OMP+            DEFAULT(NONE)
         DO I=1,NPATCH(1)
          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)
          X1=PATCHRX(I)-DX/2.0
          Y1=PATCHRY(I)-DY/2.0
          Z1=PATCHRZ(I)-DZ/2.0
          X2=X1+N1*DX/2.0
          Y2=Y1+N2*DY/2.0
          Z2=Z1+N3*DZ/2.0
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
     &         Y1.LE.Y4.AND.Y3.LE.Y2.AND.
     &         Z1.LE.Z4.AND.Z3.LE.Z2) THEN
c           WRITE(*,*) BASX,BASY,BASZ,BAS
            DO KZ=1,N3
            DO JY=1,N2
            DO IX=1,N1
             AA=SQRT((RX(IX,I)-BASX)**2 +
     &               (RY(JY,I)-BASY)**2 +
     &               (RZ(KZ,I)-BASZ)**2)
             IF (AA.LE.BAS) CONTA1(IX,JY,KZ,I)=0
            END DO
            END DO
            END DO
           END IF
          END DO
         END DO

        ! clean conta1 from overlaps
        CALL CLEAN_OVERLAPS_INT(NL,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                          SOLAP,CONTA1)

*       CONTA1=0 --> the cell is overlapping another IR cell, and slave (not to be counted),
*                    OR already belongs to a halo (not to locate a center here)
*       CONTA1=1 --> outside any halo, and not overlapped (or master)

        DO IR=1,NL
        ! then proceed with the halo finding as written in the asohf.f now, but with the new modifications
        ! especially, cycling not only around the vecinos, but also around the pares....
         DXPA=DX/(2.0**IR)
         DYPA=DY/(2.0**IR)
         DZPA=DZ/(2.0**IR)

         WRITE(*,*) '----> NEW LEVEL: IR,DX=',IR,DXPA,'<----'

         LOW1=SUM(NPATCH(0:IR-1))+1
         LOW2=SUM(NPATCH(0:IR))

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
     &         CONTA1(IX,JY,KZ,I).EQ.1) THEN
            BASX =U11(IX+1,JY,KZ,I)-U11(IX,JY,KZ,I)
            IF (BASX.LT.0.0) THEN
            BASY =U11(IX,JY+1,KZ,I)-U11(IX,JY,KZ,I)
            IF (BASY.LT.0.0) THEN
            BASZ =U11(IX,JY,KZ+1,I)-U11(IX,JY,KZ,I)
            IF (BASZ.LT.0.0) THEN
            BASXX=U11(IX,JY,KZ,I)  -U11(IX-1,JY,KZ,I)
            IF (BASXX.GT.0.0) THEN
            BASYY=U11(IX,JY,KZ,I)  -U11(IX,JY-1,KZ,I)
            IF (BASYY.GT.0.0) THEN
            BASZZ=U11(IX,JY,KZ,I)  -U11(IX,JY,KZ-1,I)
            IF (BASZZ.GT.0.0) THEN ! then it's a local maximum
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
           IF (U11(IX,JY,KZ,I).GE.CONTRASTECPEAK.AND.
     &         CONTA1(IX,JY,KZ,I).EQ.1) THEN
            BASX =U11(IX+1,JY,KZ,I)-U11(IX,JY,KZ,I)
            IF (BASX.LT.0.0) THEN
            BASY =U11(IX,JY+1,KZ,I)-U11(IX,JY,KZ,I)
            IF (BASY.LT.0.0) THEN
            BASZ =U11(IX,JY,KZ+1,I)-U11(IX,JY,KZ,I)
            IF (BASZ.LT.0.0) THEN
            BASXX=U11(IX,JY,KZ,I)  -U11(IX-1,JY,KZ,I)
            IF (BASXX.GT.0.0) THEN
            BASYY=U11(IX,JY,KZ,I)  -U11(IX,JY-1,KZ,I)
            IF (BASYY.GT.0.0) THEN
            BASZZ=U11(IX,JY,KZ,I)  -U11(IX,JY,KZ-1,I)
            IF (BASZZ.GT.0.0) THEN ! then it's a local maximum
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

         !WRITE(*,*) 'CHECK KK_ENTERO,II=',KK_ENTERO,II

         CALL SORT_CELLS_AMR(KK_ENTERO,DDD,DDDX,DDDY,DDDZ,DDDP)

c         DO IX=1,KK_ENTERO
c          WRITE(*,*) DDD(IX),DDDP(IX),DDDX(IX),DDDY(IX),DDDZ(IX),
c     &               CONTA1(DDDX(IX),DDDY(IX),DDDZ(IX),DDDP(IX))
c         END DO

         NV_GOOD=KK_ENTERO
         KK_ENTERO=0

*        Go through all the center candidates
         DO L1=1,NV_GOOD
          ICEN4(1)=DDDX(L1)
          ICEN4(2)=DDDY(L1)
          ICEN4(3)=DDDZ(L1)
          ICEN4(4)=DDDP(L1)
          IX=ICEN4(1)
          JY=ICEN4(2)
          KZ=ICEN4(3)
          IPATCH=ICEN4(4)

          KK_ENTERO=CONTA1(IX,JY,KZ,IPATCH)
          IF (KK_ENTERO.EQ.1) THEN ! this means this peak is not inside a halo yet
           XCEN=RX(IX,IPATCH)
           YCEN=RY(JY,IPATCH)
           ZCEN=RZ(KZ,IPATCH)

           IF (IR.NE.NL) THEN
            CALL RECENTER_DENSITY_PEAK(XCEN,YCEN,ZCEN,IR,U1,U11,PATCHRX,
     &                                 PATCHRY,PATCHRZ,PATCHNX,PATCHNY,
     &                                 PATCHNZ,NPATCH,LADO0,NL,IX,JY,KZ,
     &                                 IPATCH)
*         Assert we have not yet identified this same halo (dupplicated
*          due to recentering)
            FLAG_DUP=0
!$OMP PARALLEL DO SHARED(NCLUS,CLUSRX,CLUSRY,CLUSRZ,RADIO,XCEN,YCEN,
!$OMP+                   ZCEN),
!$OMP+            PRIVATE(II),
!$OMP+            REDUCTION(MAX:FLAG_DUP),
!$OMP+            DEFAULT(NONE)
            DO II=1,NCLUS
             IF ((XCEN-CLUSRX(II))**2+(YCEN-CLUSRY(II))**2+
     &           (ZCEN-CLUSRZ(II))**2.LT.RADIO(II)**2) THEN
              FLAG_DUP=1
             END IF
            END DO
            IF (FLAG_DUP.EQ.1) CYCLE
           END IF ! (IR.NE.NL)

           NCLUS=NCLUS+1
           REALCLUS(NCLUS)=-1
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
           R_EXT=0.87*DXPA ! major semidiagonal of a cube (max distance to a cell)

*          increase the radius until density falls below the virial value
           DELTA=10.0*CONTRASTEC*ROTE ! this is to ensure we enter the loop
           ITER_GROW=0
           II=0
           JJ=0

           DO WHILE(DELTA.GT.CONTRASTEC*ROTE)
            ITER_GROW=ITER_GROW+1

            IF (ITER_GROW.GT.1) THEN
             IF (R_EXT.LE.BOUNDIR) THEN
              R_INT=R_EXT
              R_EXT=MAX(R_EXT+ESP, R_EXT*ESP_LOG)
             ELSE
              WRITE(*,*) 'WARNING: growing not converged', r_int, r_ext,
     &                    boundir,iter_grow,delta/rote,kk_entero,
     &                    xcen,ycen,zcen
              STOP
             END IF
            END IF

            BOUNDIR2=R_EXT
*           tentative reach of the base grid
            BASX=XCEN-BOUNDIR2
            NX1=INT(((BASX-RADX(1))/DX)+0.5)+1
            IF (NX1.LT.1) NX1=1

            BASX=XCEN+BOUNDIR2
            NX2=INT(((BASX-RADX(1))/DX)+0.5)+1
            IF (NX2.GT.NX) NX2=NX

            BASY=YCEN-BOUNDIR2
            NY1=INT(((BASY-RADY(1))/DY)+0.5)+1
            IF (NY1.LT.1) NY1=1

            BASY=YCEN+BOUNDIR2
            NY2=INT(((BASY-RADY(1))/DY)+0.5)+1
            IF (NY2.GT.NY) NY2=NY

            BASZ=ZCEN-BOUNDIR2
            NZ1=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
            IF (NZ1.LT.1) NZ1=1

            BASZ=ZCEN+BOUNDIR2
            NZ2=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
            IF (NZ2.GT.NZ) NZ2=NZ

*          patches that could contain this halo
            CALL PATCHES_SPHERE(NPATCH,PATCHRX,PATCHRY,PATCHRZ,PATCHNX,
     &                          PATCHNY,PATCHNZ,XCEN,YCEN,ZCEN,BOUNDIR2,
     &                         IR,NL,RELEVANT_PATCHES,NRELEVANT_PATCHES)

            BASMASS_SHELL=0.0
            VOLCELL=DX*DY*DZ
!$OMP PARALLEL DO SHARED(NZ1,NZ2,NY1,NY2,NX1,NX2,RADX,RADY,RADZ,XCEN,
!$OMP+                   YCEN,ZCEN,R_INT,R_EXT,VOLCELL,U1,CR0AMR),
!$OMP+            PRIVATE(I,J,K,AA,BAS),
!$OMP+            REDUCTION(+:BASVOL,BASMASS_SHELL,BASX,BASY,BASZ,
!$OMP+                        BASDELTA,II),
!$OMP+            DEFAULT(NONE)
            DO K=NZ1,NZ2
            DO J=NY1,NY2
            DO I=NX1,NX2
             IF (CR0AMR(I,J,K).EQ.1) THEN
              AA=SQRT((RADX(I)-XCEN)**2 + (RADY(J)-YCEN)**2 +
     &                (RADZ(K)-ZCEN)**2)

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

            ! AND NOW LOOP UP TO IR-1, AND THEN ANOTHER LOOP ON IR, OVER RELEVANT_PATCHES
            DO IRR=1,IR-1
             DXPAPA=DX/2.0**IRR
             DYPAPA=DY/2.0**IRR
             DZPAPA=DZ/2.0**IRR
             LOW1=SUM(NRELEVANT_PATCHES(1:IRR-1))+1
             LOW2=SUM(NRELEVANT_PATCHES(1:IRR))
             VOLCELL=DXPAPA*DYPAPA*DZPAPA
!$OMP PARALLEL DO SHARED(RX,RY,RZ,XCEN,YCEN,ZCEN,R_INT,R_EXT,CONTA,
!$OMP+                   VOLCELL,U11,CR0AMR11,SOLAP,RELEVANT_PATCHES,
!$OMP+	                  PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2),
!$OMP+            PRIVATE(IX,JY,KZ,AA,BAS,BASINT,IPATCH,N1,N2,N3),
!$OMP+            REDUCTION(+:BASVOL,BASMASS_SHELL,BASX,BASY,BASZ,
!$OMP+                        BASDELTA,II),
!$OMP+            DEFAULT(NONE)
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
     &                   (RY(JY,IPATCH)-YCEN)**2 +
     &                   (RZ(KZ,IPATCH)-ZCEN)**2)

                 IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
                  II=II+1
                  BASVOL=BASVOL+VOLCELL

                  BAS=U11(IX,JY,KZ,IPATCH)*VOLCELL !U1 is not density contrast, but 1+delta = rho/rho_B!!!
                  BASMASS_SHELL=BASMASS_SHELL+BAS

                  BASX=BASX+RX(IX,IPATCH)*BAS
                  BASY=BASY+RY(JY,IPATCH)*BAS
                  BASZ=BASZ+RZ(KZ,IPATCH)*BAS
                  BASDELTA=BASDELTA+BAS
                 END IF
                END IF
               END IF
              END DO
              END DO
              END DO
             END DO
            END DO

            IRR=IR
            DXPAPA=DX/(2.0**IRR)
            DYPAPA=DY/(2.0**IRR)
            DZPAPA=DZ/(2.0**IRR)
            LOW1=SUM(NRELEVANT_PATCHES(1:IRR-1))+1
            LOW2=SUM(NRELEVANT_PATCHES(1:IRR))
            VOLCELL=DXPAPA*DYPAPA*DZPAPA
!$OMP PARALLEL DO SHARED(RX,RY,RZ,XCEN,YCEN,ZCEN,R_INT,R_EXT,
!$OMP+                   CONTA1,VOLCELL,U11,SOLAP,LOW1,LOW2,
!$OMP+		                 RELEVANT_PATCHES,PATCHNX,PATCHNY,PATCHNZ),
!$OMP+            PRIVATE(IX,JY,KZ,AA,BAS,BASINT,IPATCH,N1,N2,N3),
!$OMP+            REDUCTION(+:BASVOL,BASMASS_SHELL,BASX,BASY,BASZ,
!$OMP+                        BASDELTA,JJ),
!$OMP+            DEFAULT(NONE)
            DO BASINT=LOW1,LOW2
             IPATCH=RELEVANT_PATCHES(BASINT)
             !if (ir.gt.1) write(*,*) 'ipatch=',ipatch,low1,low2,basint
             N1=PATCHNX(IPATCH)
             N2=PATCHNY(IPATCH)
             N3=PATCHNZ(IPATCH)
             DO KZ=1,N3
             DO JY=1,N2
             DO IX=1,N1
              IF (SOLAP(IX,JY,KZ,IPATCH).EQ.1) THEN

               AA=SQRT((RX(IX,IPATCH)-XCEN)**2 +
     &                 (RY(JY,IPATCH)-YCEN)**2 +
     &                 (RZ(KZ,IPATCH)-ZCEN)**2)

               IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
                JJ=JJ+1
                CONTA1(IX,JY,KZ,IPATCH)=0
                BASVOL=BASVOL+VOLCELL

                BAS=U11(IX,JY,KZ,IPATCH)*VOLCELL !U1 is not density contrast, but 1+delta = rho/rho_B!!!
                BASMASS_SHELL=BASMASS_SHELL+BAS

                BASX=BASX+RX(IX,IPATCH)*BAS
                BASY=BASY+RY(JY,IPATCH)*BAS
                BASZ=BASZ+RZ(KZ,IPATCH)*BAS
                BASDELTA=BASDELTA+BAS
               END IF
              END IF
             END DO
             END DO
             END DO
            END DO

            BASMASS=BASMASS+BASMASS_SHELL*RODO*RE0**3
            IF (BASVOL.GT.0.0) DELTA=BASMASS/(BASVOL*RETE**3)

c           WRITE(*,*) DELTA/ROTE,II,JJ

c            write(*,*) iter_grow,r_ext,
c     &                 clusrx(nclus),clusry(nclus),clusrz(nclus),
c     &                 delta/rote
           END DO   ! do while (DELTA)

           RADIO(NCLUS)=R_EXT
           MASA(NCLUS)=BASMASS

           CLUSRXCM(NCLUS)=BASX/BASDELTA
           CLUSRYCM(NCLUS)=BASY/BASDELTA
           CLUSRZCM(NCLUS)=BASZ/BASDELTA

c           WRITE(*,*) '*** CLUSTER',NCLUS,'FINAL ********'
c           WRITE(*,*) CLUSRX(NCLUS),CLUSRY(NCLUS),CLUSRZ(NCLUS),
c     &                RADIO(NCLUS),MASA(NCLUS)*9.1717E18,IR
          END IF ! KK_ENTERO.EQ.1
         END DO ! L1=1,NV_GOOD

         DEALLOCATE(DDD,DDDX,DDDY,DDDZ,DDDP)

         WRITE(*,*) 'At level', IR,', no. haloes:', NHALLEV(IR)

         CALL OVERLAPPING(IR,NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &                    PARE,NCLUS,MASA,RADIO,CLUSRX,CLUSRY,CLUSRZ,
     &                    REALCLUS,LEVHAL,NHALLEV,BOUND,CONTRASTEC,RODO,
     &                    SOLAP,CR0AMR,CR0AMR11,VOL_SOLAP_LOW)

         WRITE(*,*)

         IF (IR.LT.NL) THEN
          LOW1=SUM(NPATCH(0:IR))+1
          LOW2=SUM(NPATCH(0:IR+1))

!$OMP PARALLEL DO SHARED(NPATCH,PATCHNX,PATCHNY,PATCHNZ,CONTA1,LOW1,
!$OMP+                   LOW2),
!$OMP+            PRIVATE(I,N1,N2,N3,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
          DO I=LOW1,LOW2
           N1=PATCHNX(I)
           N2=PATCHNY(I)
           N3=PATCHNZ(I)
           DO KZ=1,N3
           DO JY=1,N2
           DO IX=1,N1
            CONTA1(IX,JY,KZ,I)=1
           END DO
           END DO
           END DO
          END DO

          DXPAPA=DX/(2.0**(IR+1))
          DYPAPA=DY/(2.0**(IR+1))
          DZPAPA=DZ/(2.0**(IR+1))

!$OMP PARALLEL DO SHARED(NPATCH,PATCHNX,PATCHNY,PATCHNZ,PATCHRX,PATCHRY,
!$OMP+                   PATCHRZ,DX,DY,DZ,CLUSRX,CLUSRY,CLUSRZ,RADIO,
!$OMP+                   RX,RY,RZ,CONTA1,NCLUS,IR,LOW1,LOW2,DXPAPA,
!$OMP+                   DYPAPA,DZPAPA),
!$OMP+            PRIVATE(I,N1,N2,N3,X1,X2,Y1,Y2,Z1,Z2,BASX,BASY,BASZ,
!$OMP+                    BAS,IX,JY,KZ,AA,X3,Y3,Z3,X4,Y4,Z4),
!$OMP+            DEFAULT(NONE)
          DO I=LOW1,LOW2
           N1=PATCHNX(I)
           N2=PATCHNY(I)
           N3=PATCHNZ(I)
           X1=PATCHRX(I)-DXPAPA
           Y1=PATCHRY(I)-DYPAPA
           Z1=PATCHRZ(I)-DZPAPA
           X2=X1+N1*DXPAPA
           Y2=Y1+N2*DYPAPA
           Z2=Z1+N3*DZPAPA
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
     &          Y1.LE.Y4.AND.Y3.LE.Y2.AND.
     &          Z1.LE.Z4.AND.Z3.LE.Z2) THEN
c           WRITE(*,*) BASX,BASY,BASZ,BAS
             DO KZ=1,N3
             DO JY=1,N2
             DO IX=1,N1
              AA=SQRT((RX(IX,I)-BASX)**2 +
     &                (RY(JY,I)-BASY)**2 +
     &                (RZ(KZ,I)-BASZ)**2)
              IF (AA.LE.BAS) CONTA1(IX,JY,KZ,I)=0
             END DO
             END DO
             END DO
            END IF
           END DO
          END DO

          ! clean conta1 from overlaps
          CALL CLEAN_OVERLAPS_INT(NL,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                            SOLAP,CONTA1)

         END IF
        END DO !(IR=1,NL)

       END IF !(NL.GT.0)

       BASINT=COUNT(REALCLUS(1:NCLUS).EQ.-1)
       WRITE(*,*) '====> TOTAL NUMBER OF TENTATIVE FREE HALOES',BASINT
c       BASINT=COUNT(REALCLUS(1:NCLUS).GT.0)
c       WRITE(*,*) '====> TOTAL NUMBER OF TENTATIVE SUBSTRUCTURES',BASINT

       RETURN
       END

********************************************************************
       SUBROUTINE RECENTER_DENSITY_PEAK(BASX,BASY,BASZ,HLEV,U1,U11,
     &                                  PATCHRX,PATCHRY,PATCHRZ,PATCHNX,
     &                                  PATCHNY,PATCHNZ,NPATCH,LADO0,NL,
     &                                  PIX,PJY,PKZ,PIPATCH)
********************************************************************
*      Refines the location of a density peak using finer AMR levels
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
       SUBROUTINE PATCHES_SPHERE(NPATCH,PATCHRX,PATCHRY,PATCHRZ,
     &                           PATCHNX,PATCHNY,PATCHNZ,
     &                           XCEN,YCEN,ZCEN,REACH,MAXLEVEL,NL,
     &                           RELEVANT_PATCHES,NRELEVANT_PATCHES)
********************************************************************
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NPATCH(0:NLEVELS)
       REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       REAL XCEN,YCEN,ZCEN,REACH
       INTEGER MAXLEVEL,NL
       INTEGER RELEVANT_PATCHES(NPALEV),NRELEVANT_PATCHES(NLEVELS)

       REAL*4 DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       REAL X1,X2,X3,X4,Y1,Y2,Y3,Y4,Z1,Z2,Z3,Z4,DXPA,DYPA,DZPA
       INTEGER I,J,K,IX,JY,KZ,II,JJ,KK,LOW1,LOW2,IR,N1,N2,N3,LOW3,CONTA

       DO IR=1,NL
        NRELEVANT_PATCHES(IR)=0
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2
         RELEVANT_PATCHES(I)=0
        END DO
       END DO

       X1=XCEN-REACH
       X2=XCEN+REACH
       Y1=YCEN-REACH
       Y2=YCEN+REACH
       Z1=ZCEN-REACH
       Z2=ZCEN+REACH

       LOW3=0
       DO IR=1,MAXLEVEL
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))

        DXPA=DX/2.0**IR
        DYPA=DY/2.0**IR
        DZPA=DZ/2.0**IR

        CONTA=0
        DO I=LOW1,LOW2
         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         X3=PATCHRX(I)-DXPA
         Y3=PATCHRY(I)-DYPA
         Z3=PATCHRZ(I)-DZPA

         X4=X3+N1*DXPA
         Y4=Y3+N2*DYPA
         Z4=Z3+N3*DZPA

         IF (X1.LE.X4.AND.X3.LE.X2) THEN
         IF (Y1.LE.Y4.AND.Y3.LE.Y2) THEN
         IF (Z1.LE.Z4.AND.Z3.LE.Z2) THEN
          CONTA=CONTA+1
          RELEVANT_PATCHES(LOW3+CONTA)=I
         END IF
         END IF
         END IF
        END DO
        NRELEVANT_PATCHES(IR)=CONTA
        LOW3=LOW3+NRELEVANT_PATCHES(IR)
       END DO


       RETURN
       END

********************************************************************
       SUBROUTINE SORT_CELLS(KK_ENTERO,DDD,DDDX,DDDY,DDDZ)
********************************************************************
*      Sorts cells decreasingly in density (DDD)
********************************************************************

       IMPLICIT NONE

       INTEGER KK_ENTERO
       REAL DDD(KK_ENTERO)
       INTEGER DDDX(KK_ENTERO),DDDY(KK_ENTERO),DDDZ(KK_ENTERO)

       REAL DDD2(KK_ENTERO)
       INTEGER DDDX2(KK_ENTERO),DDDY2(KK_ENTERO),DDDZ2(KK_ENTERO)
       INTEGER INDICE2(KK_ENTERO)
       INTEGER I

!$OMP PARALLEL DO SHARED(KK_ENTERO,DDD,DDDX,DDDY,DDDZ,DDD2,DDDX2,
!$OMP+                   DDDY2,DDDZ2),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
       DO I=1,KK_ENTERO
        DDD2(I)=1.0/DDD(I)
        DDDX2(I)=DDDX(I)
        DDDY2(I)=DDDY(I)
        DDDZ2(I)=DDDZ(I)
       END DO

       CALL INDEXX(KK_ENTERO,DDD2,INDICE2)

!$OMP PARALLEL DO SHARED(KK_ENTERO,DDD,DDDX,DDDY,DDDZ,DDD2,DDDX2,
!$OMP+                   DDDY2,DDDZ2,INDICE2),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
       DO I=1,KK_ENTERO
        DDD(I)=1.0/DDD2(INDICE2(I))
        DDDX(I)=DDDX2(INDICE2(I))
        DDDY(I)=DDDY2(INDICE2(I))
        DDDZ(I)=DDDZ2(INDICE2(I))
       END DO

       RETURN
       END

********************************************************************
       SUBROUTINE SORT_CELLS_AMR(KK_ENTERO,DDD,DDDX,DDDY,DDDZ,DDDP)
********************************************************************
*      Sorts cells decreasingly in density (DDD) (including patch number)
********************************************************************

       IMPLICIT NONE

       INTEGER KK_ENTERO
       REAL DDD(KK_ENTERO)
       INTEGER DDDX(KK_ENTERO),DDDY(KK_ENTERO),DDDZ(KK_ENTERO)
       INTEGER DDDP(KK_ENTERO)

       REAL DDD2(KK_ENTERO)
       INTEGER DDDX2(KK_ENTERO),DDDY2(KK_ENTERO),DDDZ2(KK_ENTERO)
       INTEGER INDICE2(KK_ENTERO),DDDP2(KK_ENTERO)
       INTEGER I

!$OMP PARALLEL DO SHARED(KK_ENTERO,DDD,DDDX,DDDY,DDDZ,DDD2,DDDX2,
!$OMP+                   DDDY2,DDDZ2,DDDP,DDDP2),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
       DO I=1,KK_ENTERO
        DDD2(I)=1.0/DDD(I)
        DDDX2(I)=DDDX(I)
        DDDY2(I)=DDDY(I)
        DDDZ2(I)=DDDZ(I)
        DDDP2(I)=DDDP(I)
       END DO

       CALL INDEXX(KK_ENTERO,DDD2,INDICE2)

!$OMP PARALLEL DO SHARED(KK_ENTERO,DDD,DDDX,DDDY,DDDZ,DDD2,DDDX2,
!$OMP+                   DDDY2,DDDZ2,INDICE2,DDDP,DDDP2),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
       DO I=1,KK_ENTERO
        DDD(I)=1.0/DDD2(INDICE2(I))
        DDDX(I)=DDDX2(INDICE2(I))
        DDDY(I)=DDDY2(INDICE2(I))
        DDDZ(I)=DDDZ2(INDICE2(I))
        DDDP(I)=DDDP2(INDICE2(I))
       END DO

       RETURN
       END
