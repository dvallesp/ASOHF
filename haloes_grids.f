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
     &                        VECINO,NVECI,CR0AMR,CR0AMR11,
     &                        VOL_SOLAP_LOW)
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
       INTEGER VECINO(NPALEV,NPALEV),NVECI(NPALEV)
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
       REAL*4 U1G(NMAX,NMAY,NMAZ)
       REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL*4 U11G(NAMRX,NAMRY,NAMRZ,NPALEV)
       COMMON /VARIA/ U1,U11,U1G,U11G

       REAL*4 ACHE,T0,RE0
       COMMON /DOS/ ACHE,T0,RE0

*      LOCAL VARIABLES
       INTEGER,ALLOCATABLE::SOLAPA(:,:),NSOLAP(:)
       INTEGER IRR,I,J,K,IX,JY,KZ,IPATCH,LOWH1,LOWH2,LOWH3,BASINT
       INTEGER IMAXCLUS,IMINCLUS,CONTA,FLAG_ITER,NUM_ITERS,BASINT2
       REAL BAS,BASX,BASY,BASZ,R1,R2,DIST,X1,Y1,Z1,X2,Y2,Z2,VOL1,VOL2
       REAL VINT,XI,PI,SOLAP_LOWER_THR,M1,M2

       SOLAP_LOWER_THR=VOL_SOLAP_LOW

       PI=DACOS(-1.D0)

       WRITE(*,*) '== HALOES OVERLAPPING IN IR =', IR
       WRITE(*,*) 'NCLUS,NHALLEV(IR)=', NCLUS, NHALLEV(IR)

       LOWH1=SUM(NHALLEV(0:IR-1))+1
       LOWH2=SUM(NHALLEV(0:IR))

       ALLOCATE(NSOLAP(LOWH1:LOWH2), SOLAPA(LOWH1:LOWH2,NMAXSOLAP))

!$OMP PARALLEL DO SHARED(NSOLAP,SOLAPA), PRIVATE(I,J), DEFAULT(NONE)
       DO I=LOWH1,LOWH2
        NSOLAP(I)=0
        DO J=1,NMAXSOLAP
         SOLAPA(I,J)=0
        END DO
       END DO

       NUM_ITERS=0
       FLAG_ITER=1

       DO WHILE (FLAG_ITER.EQ.1)
        NUM_ITERS=NUM_ITERS+1
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
          DO J=I+1,LOWH2
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
          DO J=I+1,LOWH2
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
              IF (R2.GT.R1) THEN
               IMAXCLUS=J
               IMINCLUS=I
              END IF !(R2.GT.R1)

              REALCLUS(IMINCLUS)=0 !The smallest halo is removed
              CONTA=CONTA+1

              RADIO(IMAXCLUS)=(R1**3+R2**3)**(1.0/3.0)
              CLUSRX(IMAXCLUS)=(M1*X1+M2*X2)/(M1+M2)
              CLUSRY(IMAXCLUS)=(M1*Y1+M2*Y2)/(M1+M2)
              CLUSRZ(IMAXCLUS)=(M1*Z1+M2*Z2)/(M1+M2)

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

        IF (CONTA.EQ.0) FLAG_ITER=0
        !WRITE(*,*) 'OVERLAPPING',IR,NUM_ITERS,CONTA
       END DO !WHILE (FLAG_ITER.EQ.1)

       BASINT=COUNT(REALCLUS(LOWH1:LOWH2).EQ.0)
       WRITE(*,*) 'REMOVED HALOS_0----->', BASINT
       BASINT=COUNT(REALCLUS(LOWH1:LOWH2).NE.0)
       WRITE(*,*) 'POSSIBLE HALOS_0----->', BASINT
       BASINT=COUNT(REALCLUS(LOWH1:LOWH2).EQ.-1)
       BASINT2=COUNT(REALCLUS(LOWH1:LOWH2).GT.0)
       WRITE(*,*) '--> Of which free,substructure:',BASINT,BASINT2
       WRITE(*,*)

       RETURN
       END

********************************************************************
       SUBROUTINE HALOFIND_GRID(NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                          PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                          PATCHRY,PATCHRZ,PARE,NCLUS,MASA,RADIO,
     &                          CLUSRX,CLUSRY,CLUSRZ,REALCLUS,LEVHAL,
     &                          NHALLEV,BOUND,CONTRASTEC,RODO,SOLAP,
     &                          VECINO,NVECI,CR0AMR,CR0AMR11,PATCHCLUS,
     &                          VOL_SOLAP_LOW)
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
       INTEGER REALCLUS(MAXNCLUS),LEVHAL(MAXNCLUS)
       INTEGER PATCHCLUS(MAXNCLUS)
       INTEGER NHALLEV(0:NLEVELS)
       REAL BOUND,CONTRASTEC,RODO
       INTEGER SOLAP(NAMRX,NAMRY,NAMRZ,NPALEV)
       INTEGER VECINO(NPALEV,NPALEV),NVECI(NPALEV)
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
       REAL*4 U1G(NMAX,NMAY,NMAZ)
       REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL*4 U11G(NAMRX,NAMRY,NAMRZ,NPALEV)
       COMMON /VARIA/ U1,U11,U1G,U11G

       REAL*4 ACHE,T0,RE0
       COMMON /DOS/ ACHE,T0,RE0

*      LOCAL VARIABLES
       INTEGER CONTA(NMAX,NMAY,NMAZ)
       INTEGER CONTA1(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL UBAS1(NMAX,NMAY,NMAZ)

       REAL MAXIMO(NPALEV)
       INTEGER VID(NLEVELS,NPALEV),NVID(NLEVELS)
       INTEGER RELEVANT_PATCHES(NPALEV),NRELEVANT_PATCHES(NLEVELS)

       INTEGER IR,IX,JY,KZ,I,J,K,II,JJ,KK,IPATCH,ICEN(3),NV_GOOD
       INTEGER L1,L2,L3,NX1,NX2,NY1,NY2,NZ1,NZ2,KK_ENTERO,ITER_GROW
       INTEGER N1,N2,N3,KONTA,LOW1,LOW2,I2,ICEN1(1),ICEN4(4),CEL,I1,J1
       INTEGER K1,BORAMR,IRR,BASINT,FLAG_ITER
       REAL PRUEBAX,PRUEBAY,PRUEBAZ,RMIN,BASMASS_SHELL,BASMASS,DELTA
       REAL ESP,ESP_LOG,BAS,KK_REAL,RSHELL,R_INT,R_EXT,RANT
       REAL BASDELTA,AA,PI,VOLCELL,BASX,BASY,BASZ,BASVOL
       REAL X1,X2,Y1,Y2,Z1,Z2,DXPA,DYPA,DZPA,BASXX,BASYY,BASZZ
       REAL XCEN,YCEN,ZCEN,BOUNDIR,X3,Y3,Z3,X4,Y4,Z4,MINDERIV
       REAL VECDENS(1000),VECRAD(1000),DERIVATIVE(1000),BASVOL_SHELL

       REAL*4, ALLOCATABLE::DDD(:)
       INTEGER, ALLOCATABLE::DDDX(:),DDDY(:),DDDZ(:),DDDP(:)

**************************************************************
*      NIVEL BASE!!
**************************************************************

       PI=DACOS(-1.D0)

       IR=0
       ESP=0.2*DX
       ESP_LOG=1.05
       WRITE(*,*) '<-------- BASE GRID -------->'
       WRITE(*,*) '--> IR,DX:',IR,DX

!$OMP PARALLEL DO SHARED(NX,NY,NZ,CONTA,UBAS1,U1),
!$OMP+            PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO KZ=1,NZ
       DO JY=1,NY
       DO IX=1,NX
        CONTA(IX,JY,KZ)=0
        UBAS1(IX,JY,KZ)=U1(IX,JY,KZ)
       END DO
       END DO
       END DO

       KK_ENTERO=0
!$OMP PARALLEL DO SHARED(NX,NY,NZ,UBAS1,CONTRASTEC),
!$OMP+            PRIVATE(IX,JY,KZ),
!$OMP+            REDUCTION(+:KK_ENTERO),
!$OMP+            DEFAULT(NONE)
       DO KZ=1,NZ
       DO JY=1,NY
       DO IX=1,NX
        IF (UBAS1(IX,JY,KZ).GE.CONTRASTEC) KK_ENTERO=KK_ENTERO+1
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
       DO KZ=1,NZ
       DO JY=1,NY
       DO IX=1,NX
        BAS=UBAS1(IX,JY,KZ)
        IF (BAS.GE.CONTRASTEC) THEN
         II=II+1
         DDD(II)=BAS
         DDDX(II)=IX
         DDDY(II)=JY
         DDDZ(II)=KZ
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
        IF(KK_ENTERO.EQ.0) THEN ! this means this peak is not inside a halo yet
c         WRITE(*,*) U1(ICEN(1),ICEN(2),ICEN(3))
         NCLUS=NCLUS+1
         REALCLUS(NCLUS)=-1
         LEVHAL(NCLUS)=IR
         NHALLEV(IR)=NHALLEV(IR)+1
         PATCHCLUS(NCLUS)=0

         IF(NCLUS.GT.MAXNCLUS) THEN
          WRITE(*,*) 'WARNING: NCLUS>MAXNCLUS!!!',NCLUS,MAXNCLUS
          STOP
         END IF

         CLUSRX(NCLUS)=RADX(ICEN(1))
         CLUSRY(NCLUS)=RADY(ICEN(2))
         CLUSRZ(NCLUS)=RADZ(ICEN(3))

*        tentative reach of the halo ---> build a mini box around it
*        (by excess, set as a parameter in asohf.dat)
         BASX=CLUSRX(NCLUS)-BOUND
         NX1=INT(((BASX-RADX(1))/DX)+0.5)+1
         IF (NX1.LT.1) NX1=1

         BASX=CLUSRX(NCLUS)+BOUND
         NX2=INT(((BASX-RADX(1))/DX)+0.5)+1
         IF (NX2.GT.NX) NX2=NX

         BASY=CLUSRY(NCLUS)-BOUND
         NY1=INT(((BASY-RADY(1))/DY)+0.5)+1
         IF (NY1.LT.1) NY1=1

         BASY=CLUSRY(NCLUS)+BOUND
         NY2=INT(((BASY-RADY(1))/DY)+0.5)+1
         IF (NY2.GT.NY) NY2=NY

         BASZ=CLUSRZ(NCLUS)-BOUND
         NZ1=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
         IF (NZ1.LT.1) NZ1=1

         BASZ=CLUSRZ(NCLUS)+BOUND
         NZ2=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
         IF (NZ2.GT.NZ) NZ2=NZ

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
         R_EXT=0.5*DX

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

          VOLCELL=DX*DY*DZ
          BASMASS_SHELL=0.0
          II=0
          DO K=NZ1,NZ2
          DO J=NY1,NY2
          DO I=NX1,NX2
           AA=SQRT((RADX(I)-CLUSRX(NCLUS))**2 +
     &             (RADY(J)-CLUSRY(NCLUS))**2 +
     &             (RADZ(K)-CLUSRZ(NCLUS))**2)

           IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
            CONTA(I,J,K)=1 ! do not try to find an additional halo here (at this level)

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
          DELTA=BASMASS/(BASVOL*RETE**3)
c          WRITE(*,*) DELTA/ROTE, II
         END DO   ! do while (DELTA)

         RADIO(NCLUS)=R_EXT
         MASA(NCLUS)=BASMASS

         CLUSRX(NCLUS)=BASX/BASDELTA
         CLUSRY(NCLUS)=BASY/BASDELTA
         CLUSRZ(NCLUS)=BASZ/BASDELTA

C         WRITE(*,*) CLUSRX(NCLUS),CLUSRY(NCLUS),CLUSRZ(NCLUS),
C     &             RADIO(NCLUS),MASA(NCLUS)*9.1717E18,IR

        END IF ! KK_ENTERO.EQ.0
       END DO ! I=1,NV_GOOD

       DEALLOCATE(DDD,DDDX,DDDY,DDDZ)

       WRITE(*,*) 'At level',IR,', num. haloes:', NHALLEV(IR)
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
          DO KZ=1,N1
          DO JY=1,N2
          DO IX=1,N3
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
            DO KZ=1,N1
            DO JY=1,N2
            DO IX=1,N3
             AA=SQRT((RX(IX,I)-BASX)**2 +
     &               (RY(JY,I)-BASY)**2 +
     &               (RZ(KZ,I)-BASZ)**2)
             IF (AA.LE.BAS) CONTA1(IX,JY,KZ,I)=-1
            END DO
            END DO
            END DO
           END IF
          END DO
         END DO

*       And mark the centers of haloes (to avoid identifying haloes at same levels as substructure)
!$OMP PARALLEL DO SHARED(NPATCH,DX,DY,DZ,PATCHNX,PATCHNY,PATCHNZ,
!$OMP+                   PATCHRX,PATCHRY,PATCHRZ,NCLUS,CLUSRX,CLUSRY,
!$OMP+                   CLUSRZ,CONTA1),
!$OMP+            PRIVATE(I,N1,N2,N3,X1,X2,Y1,Y2,Z1,Z2,II,BASX,BASY,
!$OMP+                    BASZ,IX,JY,KZ),
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
          BASX=(BASX-X1)*(X2-BASX)
          BASY=(BASY-Y1)*(Y2-BASY)
          BASZ=(BASZ-Z1)*(Z2-BASZ)
          IF (BASX.GT.0.AND.BASY.GT.0.AND.BASZ.GT.0) THEN
           BASX=CLUSRX(II)-X1
           BASY=CLUSRY(II)-Y1
           BASZ=CLUSRZ(II)-Z1
           IX=INT(BASX/(DX/2.0))+1
           JY=INT(BASY/(DY/2.0))+1
           KZ=INT(BASZ/(DZ/2.0))+1
           DO I1=IX-1,IX+1
           DO J1=JY-1,JY+1
           DO K1=KZ-1,KZ+1
            IF (I1.GE.1.AND.I1.LE.N1.AND.
     &          J1.GE.1.AND.J1.LE.N2.AND.
     &          K1.GE.1.AND.K1.LE.N3) THEN
             CONTA1(I1,J1,K1,I)=-2
            END IF
           END DO
           END DO
           END DO
          END IF
         END DO
        END DO

        ! clean conta1 from overlaps
        CALL CLEAN_OVERLAPS_INT(NL,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                          SOLAP,CONTA1)

*       CONTA1=0 --> the cell is overlapping another IR cell, and slave (not to be counted)
*       CONTA1=1 --> outside any halo
*       CONTA1=2 --> inside a halo at this same level IR
*       CONTA1=-1 --> outside haloes at this same level, but inside an IR-1 halo
*       CONTA1=-2 --> outside haloes at this same level, but center of a IR-1 halo

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
         BOUNDIR=BOUND/1.75**IR

*        estimation to allocate
         KK_ENTERO=0
!$OMP  PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,U11,CONTRASTEC,LOW1,
!$OMP+                    LOW2,CONTA1,BORAMR),
!$OMP+             PRIVATE(N1,N2,N3,I,IX,JY,KZ),
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
     &         CONTA1(IX,JY,KZ,I).NE.0) KK_ENTERO=KK_ENTERO+1
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
     &         CONTA1(IX,JY,KZ,I).NE.0) THEN
            II=II+1
            DDD(II)=U11(IX,JY,KZ,I)
            DDDX(II)=IX
            DDDY(II)=JY
            DDDZ(II)=KZ
            DDDP(II)=I
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

         ! cycle through the cells, and see how to grow the halo, depending
         ! on whether it is virialised halo or substructure
         ! problem: we might need to check whether candidates to substructure
         ! are bonafide density maxima, since excluding the center cell
         ! does not ensure the code will not pick a dense cell just next to it.
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
          XCEN=RX(IX,IPATCH)
          YCEN=RY(JY,IPATCH)
          ZCEN=RZ(KZ,IPATCH)

          KK_ENTERO=CONTA1(IX,JY,KZ,IPATCH)
          IF (KK_ENTERO.EQ.1) THEN ! this means this peak is not inside a halo yet
           NCLUS=NCLUS+1
           REALCLUS(NCLUS)=-1
           LEVHAL(NCLUS)=IR
           NHALLEV(IR)=NHALLEV(IR)+1
           PATCHCLUS(NCLUS)=IPATCH

           IF(NCLUS.GT.MAXNCLUS) THEN
            WRITE(*,*) 'WARNING: NCLUS>MAXNCLUS!!!',NCLUS,MAXNCLUS
            STOP
           END IF

           CLUSRX(NCLUS)=XCEN
           CLUSRY(NCLUS)=YCEN
           CLUSRZ(NCLUS)=ZCEN

*          tentative reach of the base grid
           BASX=CLUSRX(NCLUS)-BOUNDIR
           NX1=INT(((BASX-RADX(1))/DX)+0.5)+1
           IF (NX1.LT.1) NX1=1

           BASX=CLUSRX(NCLUS)+BOUNDIR
           NX2=INT(((BASX-RADX(1))/DX)+0.5)+1
           IF (NX2.GT.NX) NX2=NX

           BASY=CLUSRY(NCLUS)-BOUNDIR
           NY1=INT(((BASY-RADY(1))/DY)+0.5)+1
           IF (NY1.LT.1) NY1=1

           BASY=CLUSRY(NCLUS)+BOUNDIR
           NY2=INT(((BASY-RADY(1))/DY)+0.5)+1
           IF (NY2.GT.NY) NY2=NY

           BASZ=CLUSRZ(NCLUS)-BOUNDIR
           NZ1=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
           IF (NZ1.LT.1) NZ1=1

           BASZ=CLUSRZ(NCLUS)+BOUNDIR
           NZ2=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
           IF (NZ2.GT.NZ) NZ2=NZ

*          patches that could contain this halo
           CALL PATCHES_SPHERE(NPATCH,PATCHRX,PATCHRY,PATCHRZ,PATCHNX,
     &                         PATCHNY,PATCHNZ,XCEN,YCEN,ZCEN,BOUNDIR,
     &                         IR,NL,RELEVANT_PATCHES,NRELEVANT_PATCHES)

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
           R_EXT=0.5*DXPA

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
     &                    boundir,iter_grow,delta/rote,kk_entero
              STOP
             END IF
            END IF

            BASMASS_SHELL=0.0

            VOLCELL=DX*DY*DZ
            DO K=NZ1,NZ2
            DO J=NY1,NY2
            DO I=NX1,NX2
             IF (CR0AMR(I,J,K).EQ.1) THEN
              AA=SQRT((RADX(I)-CLUSRX(NCLUS))**2 +
     &                (RADY(J)-CLUSRY(NCLUS))**2 +
     &                (RADZ(K)-CLUSRZ(NCLUS))**2)

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
             DXPA=DX/2.0**IRR
             DYPA=DY/2.0**IRR
             DZPA=DZ/2.0**IRR
             LOW1=SUM(NRELEVANT_PATCHES(1:IRR-1))+1
             LOW2=SUM(NRELEVANT_PATCHES(1:IRR))
             VOLCELL=DXPA*DYPA*DZPA
             DO BASINT=LOW1,LOW2
              IPATCH=RELEVANT_PATCHES(BASINT)
              N1=PATCHNX(IPATCH)
              N2=PATCHNY(IPATCH)
              N3=PATCHNZ(IPATCH)
              DO KZ=1,N3
              DO JY=1,N2
              DO IX=1,N1
               IF (CR0AMR11(IX,JY,KZ,IPATCH).EQ.1) THEN
                IF (CONTA1(IX,JY,KZ,IPATCH).NE.0) THEN
                 ! DO STUFF
                 AA=SQRT((RX(IX,IPATCH)-CLUSRX(NCLUS))**2 +
     &                   (RY(JY,IPATCH)-CLUSRY(NCLUS))**2 +
     &                   (RZ(KZ,IPATCH)-CLUSRZ(NCLUS))**2)

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
            DXPA=DX/(2.0**IRR)
            DYPA=DY/(2.0**IRR)
            DZPA=DZ/(2.0**IRR)
            LOW1=SUM(NRELEVANT_PATCHES(1:IRR-1))+1
            LOW2=SUM(NRELEVANT_PATCHES(1:IRR))
            VOLCELL=DXPA*DYPA*DZPA
            DO BASINT=LOW1,LOW2
             IPATCH=RELEVANT_PATCHES(BASINT)
             !if (ir.gt.1) write(*,*) 'ipatch=',ipatch,low1,low2,basint
             N1=PATCHNX(IPATCH)
             N2=PATCHNY(IPATCH)
             N3=PATCHNZ(IPATCH)
             DO KZ=1,N3
             DO JY=1,N2
             DO IX=1,N1
              IF (CONTA1(IX,JY,KZ,IPATCH).NE.0) THEN

               AA=SQRT((RX(IX,IPATCH)-CLUSRX(NCLUS))**2 +
     &                 (RY(JY,IPATCH)-CLUSRY(NCLUS))**2 +
     &                 (RZ(KZ,IPATCH)-CLUSRZ(NCLUS))**2)

               IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
                JJ=JJ+1
                CONTA1(IX,JY,KZ,IPATCH)=2 ! do not try to find an additional halo here (at this level)
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
            DELTA=BASMASS/(BASVOL*RETE**3)

c            write(*,*) iter_grow,r_ext,
c     &                 clusrx(nclus),clusry(nclus),clusrz(nclus),
c     &                 delta/rote
           END DO   ! do while (DELTA)

           RADIO(NCLUS)=R_EXT
           MASA(NCLUS)=BASMASS

           CLUSRX(NCLUS)=BASX/BASDELTA
           CLUSRY(NCLUS)=BASY/BASDELTA
           CLUSRZ(NCLUS)=BASZ/BASDELTA

c           WRITE(*,*) IR,CLUSRX(NCLUS),CLUSRY(NCLUS),CLUSRZ(NCLUS),
c     &             RADIO(NCLUS),MASA(NCLUS)*9.1717E18,II,JJ,DELTA/ROTE,
c     &             NCLUS
          ELSE IF(KK_ENTERO.EQ.-1) THEN ! this mean this peak must be a substructure
           BASX =U11(IX+1,JY,KZ,IPATCH)-U11(IX,JY,KZ,IPATCH)
           BASY =U11(IX,JY+1,KZ,IPATCH)-U11(IX,JY,KZ,IPATCH)
           BASZ =U11(IX,JY,KZ+1,IPATCH)-U11(IX,JY,KZ,IPATCH)
           BASXX=U11(IX,JY,KZ,IPATCH)  -U11(IX-1,JY,KZ,IPATCH)
           BASYY=U11(IX,JY,KZ,IPATCH)  -U11(IX,JY-1,KZ,IPATCH)
           BASZZ=U11(IX,JY,KZ,IPATCH)  -U11(IX,JY,KZ-1,IPATCH)
           IF (BASX.LT.0) THEN
           IF (BASY.LT.0) THEN
           IF (BASZ.LT.0) THEN
           IF (BASXX.GT.0) THEN
           IF (BASYY.GT.0) THEN
           IF (BASZZ.GT.0) THEN ! then it's a local maximum
C            WRITE(*,*) XCEN,YCEN,ZCEN,'substructure',ir
C            WRITE(*,*) 'CHECK:',
C     &                 MINVAL(U11(IX-1:IX+1,JY-1:JY+1,KZ-1:KZ+1,IPATCH))

            NCLUS=NCLUS+1
            REALCLUS(NCLUS)=-1
            LEVHAL(NCLUS)=IR
            NHALLEV(IR)=NHALLEV(IR)+1
            PATCHCLUS(NCLUS)=IPATCH

            IF(NCLUS.GT.MAXNCLUS) THEN
             WRITE(*,*) 'WARNING: NCLUS>MAXNCLUS!!!',NCLUS,MAXNCLUS
             STOP
            END IF

            CLUSRX(NCLUS)=XCEN
            CLUSRY(NCLUS)=YCEN
            CLUSRZ(NCLUS)=ZCEN

*           tentative reach of the base grid
            BASX=CLUSRX(NCLUS)-BOUNDIR
            NX1=INT(((BASX-RADX(1))/DX)+0.5)+1
            IF (NX1.LT.1) NX1=1

            BASX=CLUSRX(NCLUS)+BOUNDIR
            NX2=INT(((BASX-RADX(1))/DX)+0.5)+1
            IF (NX2.GT.NX) NX2=NX

            BASY=CLUSRY(NCLUS)-BOUNDIR
            NY1=INT(((BASY-RADY(1))/DY)+0.5)+1
            IF (NY1.LT.1) NY1=1

            BASY=CLUSRY(NCLUS)+BOUNDIR
            NY2=INT(((BASY-RADY(1))/DY)+0.5)+1
            IF (NY2.GT.NY) NY2=NY

            BASZ=CLUSRZ(NCLUS)-BOUNDIR
            NZ1=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
            IF (NZ1.LT.1) NZ1=1

            BASZ=CLUSRZ(NCLUS)+BOUNDIR
            NZ2=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
            IF (NZ2.GT.NZ) NZ2=NZ

*           patches that could contain this halo
            CALL PATCHES_SPHERE(NPATCH,PATCHRX,PATCHRY,PATCHRZ,PATCHNX,
     &                          PATCHNY,PATCHNZ,XCEN,YCEN,ZCEN,BOUNDIR,
     &                          IR,NL,RELEVANT_PATCHES,
     &                          NRELEVANT_PATCHES)

*           Now, we extend the cluster radially from its center
            BASX=0.0
            BASY=0.0
            BASZ=0.0
            BASDELTA=0.0

            BASMASS_SHELL=0.0
            BASMASS=0.0   !TOTAL MASS OF THE CLUSTER
            DELTA=0.0    !TOTAL CONTRAST OF THE CLUSTER
            BASVOL=0.0     !TOTAL VOLUME OF THE CLUSTER (sphere)

            R_INT=0.0
            R_EXT=0.5*DXPA

*           increase the radius until density falls below the virial value
            DELTA=10.0*CONTRASTEC*ROTE ! this is to ensure we enter the loop
            ITER_GROW=0
            FLAG_ITER=1
            JJ=0
            MINDERIV=1000.0

            DO WHILE(DELTA.GT.CONTRASTEC*ROTE.AND.FLAG_ITER.EQ.1)
             ITER_GROW=ITER_GROW+1
             II=0

             IF (ITER_GROW.GT.1) THEN
              IF (R_EXT.LE.BOUNDIR) THEN
               R_INT=R_EXT
               R_EXT=MAX(R_EXT+ESP, R_EXT*ESP_LOG)
              ELSE
               WRITE(*,*) 'WARNING: growing not converged', r_int,
     &                    r_ext,boundir,iter_grow,delta/rote,kk_entero
               STOP
              END IF
             END IF

             BASMASS_SHELL=0.0
             BASVOL_SHELL=0.0

             VOLCELL=DX*DY*DZ
             DO K=NZ1,NZ2
             DO J=NY1,NY2
             DO I=NX1,NX2
              IF (CR0AMR(I,J,K).EQ.1) THEN
               AA=SQRT((RADX(I)-CLUSRX(NCLUS))**2 +
     &                 (RADY(J)-CLUSRY(NCLUS))**2 +
     &                 (RADZ(K)-CLUSRZ(NCLUS))**2)

               IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
                !CONTA(I,J,K)=1 ! do not try to find an additional halo here (at this level)
                II=II+1
                BASVOL=BASVOL+VOLCELL

                BAS=U1(I,J,K)*VOLCELL !U1 is not density contrast, but 1+delta = rho/rho_B!!!
                BASMASS_SHELL=BASMASS_SHELL+BAS
                BASVOL_SHELL=BASVOL_SHELL+VOLCELL

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
              DXPA=DX/2.0**IRR
              DYPA=DY/2.0**IRR
              DZPA=DZ/2.0**IRR
              LOW1=SUM(NRELEVANT_PATCHES(1:IRR-1))+1
              LOW2=SUM(NRELEVANT_PATCHES(1:IRR))
              VOLCELL=DXPA*DYPA*DZPA
              DO BASINT=LOW1,LOW2
               IPATCH=RELEVANT_PATCHES(BASINT)
               N1=PATCHNX(IPATCH)
               N2=PATCHNY(IPATCH)
               N3=PATCHNZ(IPATCH)
               DO KZ=1,N3
               DO JY=1,N2
               DO IX=1,N1
                IF (CR0AMR11(IX,JY,KZ,IPATCH).EQ.1) THEN
                 IF (CONTA1(IX,JY,KZ,IPATCH).NE.0) THEN
                  ! DO STUFF
                  AA=SQRT((RX(IX,IPATCH)-CLUSRX(NCLUS))**2 +
     &                    (RY(JY,IPATCH)-CLUSRY(NCLUS))**2 +
     &                    (RZ(KZ,IPATCH)-CLUSRZ(NCLUS))**2)

                  IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
                   II=II+1
                   BASVOL=BASVOL+VOLCELL

                   BAS=U11(IX,JY,KZ,IPATCH)*VOLCELL !U1 is not density contrast, but 1+delta = rho/rho_B!!!
                   BASMASS_SHELL=BASMASS_SHELL+BAS
                   BASVOL_SHELL=BASVOL_SHELL+VOLCELL

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
             DXPA=DX/(2.0**IRR)
             DYPA=DY/(2.0**IRR)
             DZPA=DZ/(2.0**IRR)
             LOW1=SUM(NRELEVANT_PATCHES(1:IRR-1))+1
             LOW2=SUM(NRELEVANT_PATCHES(1:IRR))
             VOLCELL=DXPA*DYPA*DZPA
             DO BASINT=LOW1,LOW2
              IPATCH=RELEVANT_PATCHES(BASINT)
              !if (ir.gt.1) write(*,*) 'ipatch=',ipatch,low1,low2,basint
              N1=PATCHNX(IPATCH)
              N2=PATCHNY(IPATCH)
              N3=PATCHNZ(IPATCH)
              DO KZ=1,N3
              DO JY=1,N2
              DO IX=1,N1
               IF (CONTA1(IX,JY,KZ,IPATCH).NE.0) THEN

                AA=SQRT((RX(IX,IPATCH)-CLUSRX(NCLUS))**2 +
     &                  (RY(JY,IPATCH)-CLUSRY(NCLUS))**2 +
     &                  (RZ(KZ,IPATCH)-CLUSRZ(NCLUS))**2)

                IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
                 II=II+1
                 CONTA1(IX,JY,KZ,IPATCH)=2 ! do not try to find an additional halo here (at this level)
                 BASVOL=BASVOL+VOLCELL

                 BAS=U11(IX,JY,KZ,IPATCH)*VOLCELL !U1 is not density contrast, but 1+delta = rho/rho_B!!!
                 BASMASS_SHELL=BASMASS_SHELL+BAS
                 BASVOL_SHELL=BASVOL_SHELL+VOLCELL

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
             DELTA=BASMASS/(BASVOL*RETE**3)

             IF (II.GT.0) THEN
              JJ=JJ+1
              VECRAD(JJ)=0.5*(R_INT+R_EXT)
              VECDENS(JJ)=BASMASS_SHELL*RODO*RE0**3 /
     &                     (BASVOL_SHELL*RETE**3)
              IF (JJ.EQ.1) THEN
               DERIVATIVE(JJ)=1000.0
              ELSE
               BAS=LOG(VECDENS(JJ))-LOG(VECDENS(JJ-1))
               DERIVATIVE(JJ)=BAS/(LOG(VECRAD(JJ))-LOG(VECRAD(JJ-1)))
               IF (DERIVATIVE(JJ).GT.MINDERIV) FLAG_ITER=0
              END IF
              BAS=DERIVATIVE(JJ)
              MINDERIV=MIN(MINDERIV,BAS+0.1*ABS(BAS))
C              WRITE(*,*) JJ,DERIVATIVE(JJ),II,DELTA/ROTE
             END IF

c            write(*,*) iter_grow,r_ext,
c     &                 clusrx(nclus),clusry(nclus),clusrz(nclus),
c     &                 delta/rote
            END DO   ! do while (DELTA.AND.FLAG_ITER)

            RADIO(NCLUS)=R_EXT
            MASA(NCLUS)=BASMASS

            CLUSRX(NCLUS)=BASX/BASDELTA
            CLUSRY(NCLUS)=BASY/BASDELTA
            CLUSRZ(NCLUS)=BASZ/BASDELTA

*           Find its parent (from higher to lower level)
            busca_pare: DO I1=IR-1,0,-1
             LOW1=SUM(NHALLEV(0:I1-1))+1
             LOW2=SUM(NHALLEV(0:I1))
             X1=CLUSRX(NCLUS)
             Y1=CLUSRY(NCLUS)
             Z1=CLUSRZ(NCLUS)
             DO II=LOW1,LOW2
              X2=CLUSRX(II)
              Y2=CLUSRY(II)
              Z2=CLUSRZ(II)
              BAS=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
              IF (BAS/(1.05*RADIO(II)).LT.1) THEN
               REALCLUS(NCLUS)=II
               EXIT busca_pare
              END IF
             END DO
            END DO busca_pare

c            WRITE(*,*) IR,CLUSRX(NCLUS),CLUSRY(NCLUS),CLUSRZ(NCLUS),
c     &             RADIO(NCLUS),MASA(NCLUS)*9.1717E18,DELTA/ROTE,
c     &             NCLUS,REALCLUS(NCLUS)

           END IF
           END IF
           END IF
           END IF
           END IF
           END IF
          END IF ! KK_ENTERO.EQ.1 (FREE HALO) OR .EQ.-1 (SUBSTRUCTURE)
         END DO ! L1=1,NV_GOOD

         DEALLOCATE(DDD,DDDX,DDDY,DDDZ,DDDP)

         WRITE(*,*) 'At level', IR,', no. haloes:', NHALLEV(IR)
         LOW1=SUM(NHALLEV(0:IR-1))+1
         LOW2=SUM(NHALLEV(0:IR))
         WRITE(*,*) '--> Of which, potential substructure:',
     &               COUNT(REALCLUS(LOW1:LOW2).GT.0)
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
           DO KZ=1,N1
           DO JY=1,N2
           DO IX=1,N3
            CONTA1(IX,JY,KZ,I)=1
           END DO
           END DO
           END DO
          END DO

          DXPA=DX/(2.0**(IR+1))
          DYPA=DY/(2.0**(IR+1))
          DZPA=DZ/(2.0**(IR+1))

!$OMP PARALLEL DO SHARED(NPATCH,PATCHNX,PATCHNY,PATCHNZ,PATCHRX,PATCHRY,
!$OMP+                   PATCHRZ,DX,DY,DZ,CLUSRX,CLUSRY,CLUSRZ,RADIO,
!$OMP+                   RX,RY,RZ,CONTA1,NCLUS,IR,LOW1,LOW2,DXPA,DYPA,
!$OMP+                   DZPA),
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
     &          Y1.LE.Y4.AND.Y3.LE.Y2.AND.
     &          Z1.LE.Z4.AND.Z3.LE.Z2) THEN
c           WRITE(*,*) BASX,BASY,BASZ,BAS
             DO KZ=1,N1
             DO JY=1,N2
             DO IX=1,N3
              AA=SQRT((RX(IX,I)-BASX)**2 +
     &                (RY(JY,I)-BASY)**2 +
     &                (RZ(KZ,I)-BASZ)**2)
              IF (AA.LE.BAS) CONTA1(IX,JY,KZ,I)=-1
             END DO
             END DO
             END DO
            END IF
           END DO
          END DO

*       And mark the centers of haloes (to avoid identifying haloes at same levels as substructure)
!$OMP PARALLEL DO SHARED(NPATCH,DX,DY,DZ,PATCHNX,PATCHNY,PATCHNZ,
!$OMP+                   PATCHRX,PATCHRY,PATCHRZ,NCLUS,CLUSRX,CLUSRY,
!$OMP+                   CLUSRZ,CONTA1,LOW1,LOW2,DXPA,DYPA,DZPA),
!$OMP+            PRIVATE(I,N1,N2,N3,X1,X2,Y1,Y2,Z1,Z2,II,BASX,BASY,
!$OMP+                    BASZ,IX,JY,KZ),
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
            BASX=(BASX-X1)*(X2-BASX)
            BASY=(BASY-Y1)*(Y2-BASY)
            BASZ=(BASZ-Z1)*(Z2-BASZ)
            IF (BASX.GT.0.AND.BASY.GT.0.AND.BASZ.GT.0) THEN
             BASX=CLUSRX(II)-X1
             BASY=CLUSRY(II)-Y1
             BASZ=CLUSRZ(II)-Z1
             IX=INT(BASX/DXPA)+1
             JY=INT(BASY/DYPA)+1
             KZ=INT(BASZ/DZPA)+1
             DO I1=IX-1,IX+1
             DO J1=JY-1,JY+1
             DO K1=KZ-1,KZ+1
              IF (I1.GE.1.AND.I1.LE.N1.AND.
     &            J1.GE.1.AND.J1.LE.N2.AND.
     &            K1.GE.1.AND.K1.LE.N3) THEN
               CONTA1(I1,J1,K1,I)=-2
              END IF
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

       DO IR=0,NL
        IF (NHALLEV(IR).GT.0) THEN
         CALL OVERLAPPING(IR,NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                        PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                        PATCHRY,PATCHRZ,PARE,NCLUS,MASA,RADIO,
     &                        CLUSRX,CLUSRY,CLUSRZ,REALCLUS,LEVHAL,
     &                        NHALLEV,BOUND,CONTRASTEC,RODO,SOLAP,
     &                        VECINO,NVECI,CR0AMR,CR0AMR11,
     &                        VOL_SOLAP_LOW)
        END IF
       END DO

       BASINT=COUNT(REALCLUS(1:NCLUS).EQ.-1)
       WRITE(*,*) '====> TOTAL NUMBER OF TENTATIVE FREE HALOES',BASINT
       BASINT=COUNT(REALCLUS(1:NCLUS).GT.0)
       WRITE(*,*) '====> TOTAL NUMBER OF TENTATIVE SUBSTRUCTURES',BASINT

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
       SUBROUTINE FIND_NEIGHBOURING_PATCHES(IR,NL,CEL,PARE,VECINO,
     &                                      NVECI,NPATCH,VID,NVID)
********************************************************************
*      Finds recursively all the patches touching a given one,
*      and does the same with their parents up to level 1
*      UNUSED NOW / NOT WORKING
********************************************************************
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER IR,NL,CEL
       INTEGER VECINO(NPALEV,NPALEV),NVECI(NPALEV),NPATCH(0:NLEVELS)
       INTEGER PARE(NPALEV),VID(NLEVELS,NPALEV),NVID(NLEVELS)

       INTEGER IRR,FLAG,NV,NV2,NV3,II,JJ,LOW1,LOW2,KK_ENTERO,IPATCH
       INTEGER CONTAP(NPALEV)

       DO II=1,NL
        NVID(II)=0
        DO JJ=1,NPALEV
         VID(II,JJ)=0
        END DO
       END DO

       IPATCH=CEL

       DO IRR=IR,1,-1

        FLAG=0

        NV=1        !first neighbour: the patch itself
        VID(IRR,NV)=IPATCH

        NV2=0
        NV3=NV

        LOW1=SUM(NPATCH(0:IRR-1))+1
        LOW2=SUM(NPATCH(0:IRR))
        write(*,*) irr,low1,low2
        CONTAP(LOW1:LOW2)=1    !(1 vale, 0 no)
        CONTAP(VID(IRR,NV))=0      !asi este no se puede autocontar
        write(*,*) irr,ipatch,nveci(vid(irr,1))
        DO WHILE (FLAG.EQ.0)
         DO II=NV2+1, NV3
          DO JJ=1, NVECI(VID(IRR,II))
           KK_ENTERO=0
           KK_ENTERO=CONTAP(VECINO(JJ,VID(IRR,II)))
           if (irr.eq.1) write(*,*) ii,jj,VECINO(JJ,VID(IRR,II)),
     &                              VID(IRR,II),kk_entero
           IF (KK_ENTERO.EQ.1) THEN
            NV=NV+1
            VID(IRR,NV)=VECINO(JJ,VID(IRR,II))
            CONTAP(VECINO(JJ,VID(IRR,II)))=0
           END IF
          END DO
         END DO

         IF (NV.EQ.NV3) THEN
          FLAG=1
         ELSE
          NV2=NV3
          NV3=NV
         END IF

        END DO  ! WHILE (FLAG)

        NVID(IRR)=NV

        IF (IRR.GT.1) IPATCH=PARE(IPATCH)

       END DO !IRR=IR,1,-1

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
