************************************************************************
      SUBROUTINE CREATE_MESH(ITER,NX,NY,NZ,NL_MESH,NPATCH,PARE,
     &           PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &           PATCHRX,PATCHRY,PATCHRZ,RXPA,RYPA,RZPA,U2DM,U3DM,
     &           U4DM,MASAP,N_PARTICLES,N_DM,N_GAS,LADO0,T,ZETA,
     &           REFINE_THR,MIN_PATCHSIZE,FRAC_REFINABLE,BOR,BORAMR,
     &           BOR_OVLP,NPART_ESP,FW1)
************************************************************************
*     Creats a mesh hierarchy for the given particle distribution
************************************************************************

      IMPLICIT NONE

      INCLUDE 'input_files/asohf_parameters.dat'

*     function parameters
      INTEGER ITER,NX,NY,NZ,NL_MESH,N_PARTICLES,N_DM,N_GAS
      INTEGER NPATCH(0:NLEVELS),PARE(NPALEV)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
      REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED),
     &       U2DM(PARTIRED),U3DM(PARTIRED),U4DM(PARTIRED),
     &       MASAP(PARTIRED)
      REAL LADO0,T,ZETA
      INTEGER REFINE_THR,MIN_PATCHSIZE
      REAL FRAC_REFINABLE
      INTEGER BOR,BORAMR
      INTEGER NPART_ESP(0:N_ESP-1),FW1

*     COMMON VARIABLES
      REAL DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      REAL  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
      COMMON /GRID/  RADX,RADY,RADZ

      REAL*4  RX(0:NAMRX+1,NPALEV),RY(0:NAMRX+1,NPALEV),
     &        RZ(0:NAMRX+1,NPALEV)
      COMMON /GRIDAMR/ RX,RY,RZ

*     LOCAL VARIABLES
      INTEGER PLEV(PARTIRED)
      INTEGER,ALLOCATABLE::CR0(:,:,:)
      INTEGER,ALLOCATABLE::CR01(:,:,:,:)
      INTEGER,ALLOCATABLE::CONTA1(:,:,:)
      INTEGER,ALLOCATABLE::CONTA11(:,:,:,:)
      REAL MAP,XL,YL,ZL,DXPA,DYPA,DZPA,BAS
      INTEGER I,IX,JY,KZ,REFINE_COUNT
      INTEGER INI_EXTENSION,NBIS,IRPA,LOW1,LOW2,IPATCH,IPARE
      INTEGER INMAX(3),INMAX2(2),I1,I2,J1,J2,K1,K2,N1,N2,N3,IR,MARCA
      INTEGER NP1,NP2,NP3,BASINT,NPALEV3,II,JJ,KK,BOR_OVLP
      INTEGER I1BIS,I2BIS,J1BIS,J2BIS,K1BIS,K2BIS

      INTEGER,ALLOCATABLE::LNPATCH(:)
      INTEGER,ALLOCATABLE::LPATCHNX(:,:),LPATCHNY(:,:),LPATCHNZ(:,:)
      INTEGER,ALLOCATABLE::LPATCHX(:,:),LPATCHY(:,:),LPATCHZ(:,:)
      REAL,ALLOCATABLE::LPATCHRX(:,:),LPATCHRY(:,:),LPATCHRZ(:,:)
      INTEGER,ALLOCATABLE::LVAL(:,:)

      REAL,ALLOCATABLE::DDD(:)
      INTEGER,ALLOCATABLE::DDDX(:),DDDY(:),DDDZ(:)

      CHARACTER*5 ITER_STRING
      WRITE(ITER_STRING, '(I5.5)') ITER

!     hard-coded parameters (for now, at least)
      INI_EXTENSION=2 !initial extension of a patch around a cell (on each direction)
      NPALEV3=(INT(NAMRX/5)**3)+1
      write(*,*) 'NPALEV3=',NPALEV3

      MAP=MAXVAL(MASAP(1:N_PARTICLES))

      PLEV=0
      DO IR=0,N_ESP-1
       LOW1=SUM(NPART_ESP(0:IR-1))+1
       LOW2=SUM(NPART_ESP(0:IR))
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PLEV,IR), PRIVATE(I), DEFAULT(NONE)
       DO I=LOW1,LOW2
        PLEV(I)=IR
       END DO
      END DO

      WRITE(*,*) 'Particle levels: min and max values:', MINVAL(PLEV),
     &           MAXVAL(PLEV)

      XL=-FLOAT(NMAX)*DX/2.
      YL=-FLOAT(NMAY)*DY/2.
      ZL=-FLOAT(NMAZ)*DZ/2.

*     FIRST LEVEL OF REFINEMENT ========================================
      IR=1
      DXPA=DX/(2.0**IR)
      DYPA=DY/(2.0**IR)
      DZPA=DZ/(2.0**IR)
      !ALLOCATE(U1(NMAX,NMAY,NMAZ))
      ALLOCATE(CONTA1(NMAX,NMAY,NMAZ))
      ALLOCATE(CR0(NMAX,NMAY,NMAZ))

!$OMP PARALLEL DO SHARED(CONTA1,CR0,NX,NY,NZ),PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
      DO KZ=1,NX
      DO JY=1,NY
      DO IX=1,NZ
       CONTA1(IX,JY,KZ)=0
       CR0(IX,JY,KZ)=0
      END DO
      END DO
      END DO

!$OMP PARALLEL DO SHARED(N_PARTICLES,RXPA,RYPA,RZPA,XL,YL,ZL,DX,DY,DZ,
!$OMP+                   NX,NY,NZ,PLEV,REFINE_THR),
!$OMP+            PRIVATE(I,IX,JY,KZ), DEFAULT(NONE)
!$OMP+            REDUCTION(+: CONTA1)
      DO I=1,N_PARTICLES
       IX=INT((RXPA(I)-XL)/DX)+1
       JY=INT((RYPA(I)-YL)/DY)+1
       KZ=INT((RZPA(I)-ZL)/DZ)+1
       IF (IX.LT.1) IX=1
       IF (IX.GT.NX) IX=NX
       IF (JY.LT.1) JY=1
       IF (JY.GT.NY) JY=NY
       IF (KZ.LT.1) KZ=1
       IF (KZ.GT.NZ) KZ=NZ

       !U1(IX,JY,KZ)=U1(IX,JY,KZ)+MASAP(I)

       IF (PLEV(I).EQ.0) THEN
        CONTA1(IX,JY,KZ)=CONTA1(IX,JY,KZ)+1
       ELSE
        CONTA1(IX,JY,KZ)=CONTA1(IX,JY,KZ)+REFINE_THR
       END IF
      END DO

!$OMP PARALLEL DO SHARED(NX,NY,NZ,BOR,CONTA1,CR0),
!$OMP+            PRIVATE(IX,JY,KZ), DEFAULT(NONE)
      DO KZ=1,NZ
      DO JY=1,NY
      DO IX=1,NX
       IF(IX.LE.BOR.OR.IX.GE.NX-BOR+1.OR.
     &    JY.LE.BOR.OR.JY.GE.NY-BOR+1.OR.
     &    KZ.LE.BOR.OR.KZ.GE.NZ-BOR+1) THEN
         CONTA1(IX,JY,KZ)=0
       END IF
       CR0(IX,JY,KZ)=CONTA1(IX,JY,KZ)
      END DO
      END DO
      END DO

      REFINE_COUNT=COUNT(CR0.GE.REFINE_THR)
      WRITE(*,*) 'REFINABLE CELLS:', REFINE_COUNT

      ALLOCATE(DDD(REFINE_COUNT),DDDX(REFINE_COUNT),
     &         DDDY(REFINE_COUNT),DDDZ(REFINE_COUNT))

      I=0
      DO KZ=1,NZ
      DO JY=1,NY
      DO IX=1,NX
       IF (CR0(IX,JY,KZ).GE.REFINE_THR) THEN
        I=I+1
        DDD(I)=FLOAT(CR0(IX,JY,KZ))
        DDDX(I)=IX
        DDDY(I)=JY
        DDDZ(I)=KZ
       END IF
      END DO
      END DO
      END DO

      CALL SORT_CELLS(REFINE_COUNT,DDD,DDDX,DDDY,DDDZ)
      write(*,*) !ddd(1:10000)

      IPATCH=0

      DO I=1,REFINE_COUNT
       IX=DDDX(I)
       JY=DDDY(I)
       KZ=DDDZ(I)
       IF (CR0(IX,JY,KZ).LT.REFINE_THR) CYCLE
       !IF (CONTA1(IX,JY,KZ).LT.REFINE_THR) EXIT

       I1=MAX(IX-INI_EXTENSION,BOR+1)
       I2=MIN(IX+INI_EXTENSION,NX-BOR)
       J1=MAX(JY-INI_EXTENSION,BOR+1)
       J2=MIN(JY+INI_EXTENSION,NY-BOR)
       K1=MAX(KZ-INI_EXTENSION,BOR+1)
       K2=MIN(KZ+INI_EXTENSION,NZ-BOR)

       N1=2*(I2-I1+1)
       N2=2*(J2-J1+1)
       N3=2*(K2-K1+1)
       !NBAS=MAXVAL(N1,N2,N3)
       !NBIS=MINVAL(N1,N2,N3)

       MARCA = 1
       DO WHILE (MARCA.EQ.1) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MARCA=0
        IF (N1.LE.NAMRX-2.AND.I1.GT.BOR+1) THEN
         BAS=FLOAT(4*COUNT(CONTA1(I1-1,J1:J2,K1:K2).GE.REFINE_THR))
     &       /(N2*N3)
         IF (BAS.GT.FRAC_REFINABLE) THEN
          I1=I1-1
          N1=2*(I2-I1+1)
          MARCA=1
         END IF
        END IF

        IF (N1.LE.NAMRX-2.AND.I2.LT.NX-BOR) THEN
         BAS=FLOAT(4*COUNT(CONTA1(I2+1,J1:J2,K1:K2).GE.REFINE_THR))
     &       /(N2*N3)
         IF (BAS.GT.FRAC_REFINABLE) THEN
          I2=I2+1
          N1=2*(I2-I1+1)
          MARCA=1
         END IF
        END IF

        IF (N2.LE.NAMRY-2.AND.J1.GT.BOR+1) THEN
         BAS=FLOAT(4*COUNT(CONTA1(I1:I2,J1-1,K1:K2).GE.REFINE_THR))
     &       /(N1*N3)
         IF (BAS.GT.FRAC_REFINABLE) THEN
          J1=J1-1
          N2=2*(J2-J1+1)
          MARCA=1
         END IF
        END IF

        IF (N2.LE.NAMRY-2.AND.J2.LT.NY-BOR) THEN
         BAS=FLOAT(4*COUNT(CONTA1(I1:I2,J2+1,K1:K2).GE.REFINE_THR))
     &       /(N1*N3)
         IF (BAS.GT.FRAC_REFINABLE) THEN
          J2=J2+1
          N2=2*(J2-J1+1)
          MARCA=1
         END IF
        END IF

        IF (N3.LE.NAMRZ-2.AND.K1.GT.BOR+1) THEN
         BAS=FLOAT(4*COUNT(CONTA1(I1:I2,J1:J2,K1-1).GE.REFINE_THR))
     &       /(N1*N2)
         IF (BAS.GT.FRAC_REFINABLE) THEN
          K1=K1-1
          N3=2*(K2-K1+1)
          MARCA=1
         END IF
        END IF

        IF (N3.LE.NAMRZ-2.AND.K2.LT.NZ-BOR) THEN
         BAS=FLOAT(4*COUNT(CONTA1(I1:I2,J1:J2,K2+1).GE.REFINE_THR))
     &       /(N1*N2)
         IF (BAS.GT.FRAC_REFINABLE) THEN
          K2=K2+1
          N3=2*(K2-K1+1)
          MARCA=1
         END IF
        END IF

       END DO !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       NBIS=MIN(N1,N2,N3)
       IF (NBIS.LE.MIN_PATCHSIZE) THEN
        DO II=I1,I2
        DO JJ=J1,J2
        DO KK=K1,K2
         IF (CR0(II,JJ,KK).GT.0) CR0(II,JJ,KK)=0
        END DO
        END DO
        END DO
        CONTA1(I1:I2,J1:J2,K1:K2)=0
       ELSE
        IPATCH=IPATCH+1
c       WRITE(*,*) IPATCH,N1,N2,N3,I1,I2,J1,J2,K1,K2,
c     &             COUNT(CONTA1(I1:I2,J1:J2,K1:K2).GE.REFINE_THR)
c       WRITE(*,*) '*',I1+BOR_OVLP,I2-BOR_OVLP,J1+BOR_OVLP,J2-BOR_OVLP,
c     &         K1+BOR_OVLP,K2-BOR_OVLP
c       write(*,*) '**',IX,JY,KZ
c       write(*,*) '***',REFINE_COUNT

        I1BIS=I1+BOR_OVLP
        I2BIS=I2-BOR_OVLP
        J1BIS=J1+BOR_OVLP
        J2BIS=J2-BOR_OVLP
        K1BIS=K1+BOR_OVLP
        K2BIS=K2-BOR_OVLP
        IF (I1.EQ.IX) I1BIS=I1
        IF (I2.EQ.IX) I2BIS=I2
        IF (J1.EQ.JY) J1BIS=J1
        IF (J2.EQ.JY) J2BIS=J2
        IF (K1.EQ.KZ) K1BIS=K1
        IF (K2.EQ.KZ) K2BIS=K2

        CONTA1(I1BIS:I2BIS,J1BIS:J2BIS,K1BIS:K2BIS)=0
        CR0(I1BIS:I2BIS,J1BIS:J2BIS,K1BIS:K2BIS)=-1

        PATCHNX(IPATCH)=N1
        PATCHNY(IPATCH)=N2
        PATCHNZ(IPATCH)=N3

        PATCHX(IPATCH)=I1
        PATCHY(IPATCH)=J1
        PATCHZ(IPATCH)=K1

        PATCHRX(IPATCH)=RADX(I1)
        PATCHRY(IPATCH)=RADY(J1)
        PATCHRZ(IPATCH)=RADZ(K1)

        PARE(IPATCH)=0
       END IF

       REFINE_COUNT=COUNT(CR0.GE.REFINE_THR)
       !WRITE(*,*) REFINE_COUNT


      END DO  !-----------------------------------------------------

      DEALLOCATE(DDD,DDDX,DDDY,DDDZ)

      NPATCH(IR)=IPATCH

      WRITE(*,*) 'At l=',1,', patches:', NPATCH(IR)
      WRITE(*,*) '  --> l=',0,' cells refined:', COUNT(CR0.EQ.-1)
      DEALLOCATE(CONTA1,CR0)
*     END FIRST LEVEL OF REFINEMENT ====================================

*     START SUBSEQUENT LEVELS OF REFINEMENT ============================
      pare_levels: DO IRPA=1,NL_MESH-1 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       IF (NPATCH(IRPA).EQ.0) THEN
         WRITE(*,*) 'Mesh building stops at level: ', IRPA
         WRITE(*,*) 'There are no more candidate patches'
         EXIT pare_levels
       END IF

       DXPA=DX/(2.0**IRPA)
       DYPA=DY/(2.0**IRPA)
       DZPA=DZ/(2.0**IRPA)

       LOW1=SUM(NPATCH(0:IRPA-1))+1
       LOW2=SUM(NPATCH(0:IRPA))
       !WRITE(*,*) IRPA, LOW1,LOW2

       ALLOCATE(CR01(1:NAMRX,1:NAMRY,1:NAMRZ,LOW1:LOW2))
       ALLOCATE(CONTA11(1:NAMRX,1:NAMRY,1:NAMRZ,LOW1:LOW2))

!$OMP PARALLEL DO SHARED(LOW1,LOW2,CR01,CONTA11),
!$OMP+            PRIVATE(IPATCH,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2
        DO KZ=1,NAMRZ
        DO JY=1,NAMRY
        DO IX=1,NAMRX
         CR01(IX,JY,KZ,IPATCH)=0
         CONTA11(IX,JY,KZ,IPATCH)=0
        END DO
        END DO
        END DO
       END DO

!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHRX,PATCHRY,PATCHRZ,DXPA,DYPA,
!$OMP+                   DZPA,PATCHNX,PATCHNY,PATCHNZ,N_PARTICLES,RXPA,
!$OMP+                   RYPA,RZPA,CONTA11,REFINE_THR,IRPA,PLEV,BORAMR,
!$OMP+                   CR01),
!$OMP+            PRIVATE(IPATCH,XL,YL,ZL,N1,N2,N3,I,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2 != = = = = = = = = = = = = = = = = = = = = =
        !WRITE(*,*) IPATCH, LOW2
        XL=PATCHRX(IPATCH)-DXPA
        YL=PATCHRY(IPATCH)-DYPA
        ZL=PATCHRZ(IPATCH)-DZPA

        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)

        DO I=1,N_PARTICLES
         IX=INT((RXPA(I)-XL)/DXPA)+1
         JY=INT((RYPA(I)-YL)/DYPA)+1
         KZ=INT((RZPA(I)-ZL)/DZPA)+1
         IF (IX.GE.1.AND.IX.LE.N1.AND.
     &       JY.GE.1.AND.JY.LE.N2.AND.
     &       KZ.GE.1.AND.KZ.LE.N3) THEN !*****************************
          !U1(IX,JY,KZ)=U1(IX,JY,KZ)+MASAP(I)
          IF (PLEV(I).LE.IRPA) THEN
           CONTA11(IX,JY,KZ,IPATCH)=CONTA11(IX,JY,KZ,IPATCH)+1
          ELSE
           CONTA11(IX,JY,KZ,IPATCH)=CONTA11(IX,JY,KZ,IPATCH)+REFINE_THR
          END IF
         END IF !*****************************************************
        END DO

        DO KZ=1,N3
        DO JY=1,N2
        DO IX=1,N1
         IF(IX.LE.BORAMR.OR.IX.GE.N1-BORAMR+1.OR.
     &      JY.LE.BORAMR.OR.JY.GE.N2-BORAMR+1.OR.
     &      KZ.LE.BORAMR.OR.KZ.GE.N3-BORAMR+1) THEN
           CONTA11(IX,JY,KZ,IPATCH)=0
         END IF
         CR01(IX,JY,KZ,IPATCH)=CONTA11(IX,JY,KZ,IPATCH)
        END DO
        END DO
        END DO

       END DO != = = = = = = = = = = = = = = = = = = = = = = = = = = = =

       WRITE(*,*) '  --> Max particles at a cell:',
     &            MAXVAL(CR01(:,:,:,LOW1:LOW2))

c       REFINE_COUNT=COUNT(CR01(:,:,:,LOW1:LOW2).GE.REFINE_THR)
c       WRITE(*,*) 'Refinable cells BEFORE cleaning at l=',IRPA,
c     &            REFINE_COUNT

       CALL VEINSGRID_REDUCED(IRPA,NPATCH,PARE,PATCHNX,PATCHNY,
     &      PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,CR01,
     &      CONTA11,LOW1,LOW2)

       REFINE_COUNT=COUNT(CR01(:,:,:,LOW1:LOW2).GE.REFINE_THR)
       WRITE(*,*) '  --> Refinable cells AFTER cleaning:',REFINE_COUNT


       ! mesh creation at the next level
       IR=IRPA+1
       !DXPA=DX/(2.0**IR)
       !DYPA=DY/(2.0**IR)
       !DZPA=DZ/(2.0**IR)

       ALLOCATE(LNPATCH(LOW1:LOW2))
       ALLOCATE(LPATCHNX(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHNY(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHNZ(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHX(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHY(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHZ(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHRX(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHRY(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHRZ(NPALEV3,LOW1:LOW2))
       ALLOCATE(LVAL(NPALEV3,LOW1:LOW2))

       LNPATCH(:)=0
       LVAL(:,:)=0

c       WRITE(*,*) 'REFINABLE CELLS:', REFINE_COUNT

!$OMP PARALLEL DO SHARED(LOW1,LOW2,CR01,REFINE_THR,NPALEV3,PATCHNX,
!$OMP+                   PATCHNY,PATCHNZ,INI_EXTENSION,BORAMR,CONTA11,
!$OMP+                   MIN_PATCHSIZE,DXPA,DYPA,DZPA,LPATCHNX,
!$OMP+                   LPATCHNY,LPATCHNZ,LPATCHX,LPATCHY,LPATCHZ,
!$OMP+                   LPATCHRX,LPATCHRY,LPATCHRZ,LVAL,PATCHRX,
!$OMP+                   PATCHRY,PATCHRZ,FRAC_REFINABLE,BOR_OVLP),
!$OMP+            PRIVATE(IPARE,REFINE_COUNT,IPATCH,INMAX,IX,JY,KZ,
!$OMP+                    BASINT,NP1,NP2,NP3,I1,I2,J1,J2,K1,K2,N1,N2,N3,
!$OMP+                    MARCA,NBIS,BAS,I1BIS,I2BIS,J1BIS,J2BIS,K1BIS,
!$OMP+                    K2BIS),
!$OMP+            DEFAULT(NONE)
       DO IPARE=LOW1,LOW2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REFINE_COUNT=COUNT(CR01(:,:,:,IPARE).GE.REFINE_THR)
        IPATCH=0
        DO WHILE (REFINE_COUNT.GT.0.AND.IPATCH.LT.NPALEV3) !------------
         INMAX=MAXLOC(CR01(:,:,:,IPARE))
         IX=INMAX(1)
         JY=INMAX(2)
         KZ=INMAX(3)
         BASINT=CR01(IX,JY,KZ,IPARE)
         !IF (CONTA1(IX,JY,KZ).LT.REFINE_THR) EXIT

         NP1=PATCHNX(IPARE)
         NP2=PATCHNY(IPARE)
         NP3=PATCHNZ(IPARE)

         I1=MAX(IX-INI_EXTENSION,BORAMR+1)
         I2=MIN(IX+INI_EXTENSION,NP1-BORAMR)
         J1=MAX(JY-INI_EXTENSION,BORAMR+1)
         J2=MIN(JY+INI_EXTENSION,NP2-BORAMR)
         K1=MAX(KZ-INI_EXTENSION,BORAMR+1)
         K2=MIN(KZ+INI_EXTENSION,NP3-BORAMR)

         N1=2*(I2-I1+1)
         N2=2*(J2-J1+1)
         N3=2*(K2-K1+1)
         !NBAS=MAXVAL(N1,N2,N3)
         !NBIS=MINVAL(N1,N2,N3)

         MARCA = 1
         DO WHILE (MARCA.EQ.1) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          MARCA=0
          IF (N1.LE.NAMRX-2.AND.I1.GT.BORAMR+1) THEN
           BAS=FLOAT(4*COUNT(CONTA11(I1-1,J1:J2,K1:K2,IPARE).GE.
     &               REFINE_THR))/(N2*N3)
           IF (BAS.GT.FRAC_REFINABLE) THEN
            I1=I1-1
            N1=2*(I2-I1+1)
            MARCA=1
           END IF
          END IF

          IF (N1.LE.NAMRX-2.AND.I2.LT.NP1-BORAMR) THEN
           BAS=FLOAT(4*COUNT(CONTA11(I2+1,J1:J2,K1:K2,IPARE).GE.
     &               REFINE_THR))/(N2*N3)
           IF (BAS.GT.FRAC_REFINABLE) THEN
            I2=I2+1
            N1=2*(I2-I1+1)
            MARCA=1
           END IF
          END IF

          IF (N2.LE.NAMRY-2.AND.J1.GT.BORAMR+1) THEN
           BAS=FLOAT(4*COUNT(CONTA11(I1:I2,J1-1,K1:K2,IPARE).GE.
     &               REFINE_THR))/(N1*N3)
           IF (BAS.GT.FRAC_REFINABLE) THEN
            J1=J1-1
            N2=2*(J2-J1+1)
            MARCA=1
           END IF
          END IF

          IF (N2.LE.NAMRY-2.AND.J2.LT.NP2-BORAMR) THEN
           BAS=FLOAT(4*COUNT(CONTA11(I1:I2,J2+1,K1:K2,IPARE).GE.
     &               REFINE_THR))/(N1*N3)
           IF (BAS.GT.FRAC_REFINABLE) THEN
            J2=J2+1
            N2=2*(J2-J1+1)
            MARCA=1
           END IF
          END IF

          IF (N3.LE.NAMRZ-2.AND.K1.GT.BORAMR+1) THEN
           BAS=FLOAT(4*COUNT(CONTA11(I1:I2,J1:J2,K1-1,IPARE).GE.
     &               REFINE_THR))/(N1*N2)
           IF (BAS.GT.FRAC_REFINABLE) THEN
            K1=K1-1
            N3=2*(K2-K1+1)
            MARCA=1
           END IF
          END IF

          IF (N3.LE.NAMRZ-2.AND.K2.LT.NP3-BORAMR) THEN
           BAS=FLOAT(4*COUNT(CONTA11(I1:I2,J1:J2,K2+1,IPARE).GE.
     &               REFINE_THR))/(N1*N2)
           IF (BAS.GT.FRAC_REFINABLE) THEN
            K2=K2+1
            N3=2*(K2-K1+1)
            MARCA=1
           END IF
          END IF

         END DO !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         NBIS=MIN(N1,N2,N3)
         IF (NBIS.LE.MIN_PATCHSIZE) THEN
          CR01(I1:I2,J1:J2,K1:K2,IPARE)=0
          CONTA11(I1:I2,J1:J2,K1:K2,IPARE)=0
         ELSE
          IPATCH=IPATCH+1
c          WRITE(*,*) 'new,pare:',IPATCH,IPARE
c          WRITE(*,*) 'N1,N2,N3,refinable:',N1,N2,N3,
c     &             COUNT(CONTA11(I1:I2,J1:J2,K1:K2,IPARE).GE.REFINE_THR)
c          write(*,*) 'x,y,z',i1,j1,k1

          I1BIS=I1+BOR_OVLP
          I2BIS=I2-BOR_OVLP
          J1BIS=J1+BOR_OVLP
          J2BIS=J2-BOR_OVLP
          K1BIS=K1+BOR_OVLP
          K2BIS=K2-BOR_OVLP
          IF (I1.EQ.IX) I1BIS=I1
          IF (I2.EQ.IX) I2BIS=I2
          IF (J1.EQ.JY) J1BIS=J1
          IF (J2.EQ.JY) J2BIS=J2
          IF (K1.EQ.KZ) K1BIS=K1
          IF (K2.EQ.KZ) K2BIS=K2

          CONTA11(I1BIS:I2BIS,J1BIS:J2BIS,K1BIS:K2BIS,IPARE)=0
          CR01(I1BIS:I2BIS,J1BIS:J2BIS,K1BIS:K2BIS,IPARE)=-1

          LPATCHNX(IPATCH,IPARE)=N1
          LPATCHNY(IPATCH,IPARE)=N2
          LPATCHNZ(IPATCH,IPARE)=N3

          LPATCHX(IPATCH,IPARE)=I1
          LPATCHY(IPATCH,IPARE)=J1
          LPATCHZ(IPATCH,IPARE)=K1

          ! remember that dxpa is the cellsize of the parent!!!
          LPATCHRX(IPATCH,IPARE)=PATCHRX(IPARE)+(I1-1.5)*DXPA
          LPATCHRY(IPATCH,IPARE)=PATCHRY(IPARE)+(J1-1.5)*DYPA
          LPATCHRZ(IPATCH,IPARE)=PATCHRZ(IPARE)+(K1-1.5)*DZPA

          LVAL(IPATCH,IPARE)=BASINT
         END IF

         REFINE_COUNT=COUNT(CR01(:,:,:,IPARE).GE.REFINE_THR)
         !WRITE(*,*) REFINE_COUNT
        END DO  !-------------------------------------------------------
       END DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       IPATCH=LOW2
       DO WHILE(COUNT(LVAL(:,:).GT.0).GT.0.AND.IPATCH.LT.NPALEV)
        INMAX2=MAXLOC(LVAL)
        I=INMAX2(1)
        IPARE=LOW1-1+INMAX2(2)
C        WRITE(*,*) LVAL(I,IPARE)
        LVAL(I,IPARE)=0

        IPATCH=IPATCH+1
        IF (IPATCH.GT.NPALEV) EXIT

        PATCHNX(IPATCH)=LPATCHNX(I,IPARE)
        PATCHNY(IPATCH)=LPATCHNY(I,IPARE)
        PATCHNZ(IPATCH)=LPATCHNZ(I,IPARE)

        PATCHX(IPATCH)=LPATCHX(I,IPARE)
        PATCHY(IPATCH)=LPATCHY(I,IPARE)
        PATCHZ(IPATCH)=LPATCHZ(I,IPARE)

        PATCHRX(IPATCH)=LPATCHRX(I,IPARE)
        PATCHRY(IPATCH)=LPATCHRY(I,IPARE)
        PATCHRZ(IPATCH)=LPATCHRZ(I,IPARE)

        PARE(IPATCH)=IPARE
       END DO

       NPATCH(IR)=IPATCH-SUM(NPATCH(0:IR-1))
       IF (SUM(NPATCH).GE.NPALEV) STOP 'NPALEV too small'

       DEALLOCATE(LNPATCH,LPATCHNX,LPATCHNY,LPATCHNZ,LPATCHRX,LPATCHRY,
     &            LPATCHRZ,LPATCHX,LPATCHY,LPATCHZ,LVAL)

       LOW1=SUM(NPATCH(0:IRPA-1))+1
       LOW2=SUM(NPATCH(0:IRPA))
       WRITE(*,*) 'At l=',IR,', patches:', NPATCH(IR)
       WRITE(*,*) '  --> l=',IRPA,' cells refined:',
     &            COUNT(CR01(:,:,:,LOW1:LOW2).EQ.-1)

       DEALLOCATE(CR01, CONTA11)

      END DO pare_levels !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
*     END SUBSEQUENT LEVELS OF REFINEMENT

      DO IR=1,NL_MESH
       IF (NPATCH(IR).EQ.0) THEN
        NL_MESH=IR-1
        EXIT
       END IF
      END DO

*     WRITING GRID DATA ON A FILE
      IF (FW1.EQ.1) THEN
       OPEN(33,FILE='./output_files/grids_asohf'//ITER_STRING//'.res',
     &      STATUS='UNKNOWN')

       WRITE(33,*) ITER,T,NL_MESH,MAXVAL(MASAP),0.0
       WRITE(33,*) ZETA
       WRITE(33,*) 0,0,0,NX,NY,NZ
       DO IR=1,NL_MESH
        WRITE(33,*) IR,NPATCH(IR),0,0,0
        WRITE(33,*) '------ within level l=',IR,'-------'
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2
         WRITE(33,*) PATCHNX(I),PATCHNY(I),PATCHNZ(I)
         WRITE(33,*) PATCHX(I),PATCHY(I),PATCHZ(I)
         WRITE(33,*) PATCHRX(I),PATCHRY(I),PATCHRZ(I)
         WRITE(33,*) PARE(I)
        END DO
       END DO
       CLOSE(33)
      END IF

*     Build the AMR mesh
      DO IR=1,NL_MESH
       DXPA=DX/2.0**IR
       DYPA=DY/2.0**IR
       DZPA=DZ/2.0**IR
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,PATCHRX,
!$OMP+                   PATCHRY,PATCHRZ,DXPA,DYPA,DZPA,RX,RY,RZ),
!$OMP+            PRIVATE(I,N1,N2,N3,XL,YL,ZL,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        XL=PATCHRX(I)-0.5*DXPA
        YL=PATCHRY(I)-0.5*DYPA
        ZL=PATCHRZ(I)-0.5*DZPA
        DO IX=0,N1+1
         RX(IX,I)=XL+(IX-1)*DXPA
        END DO
        DO JY=0,N2+1
         RY(JY,I)=YL+(JY-1)*DYPA
        END DO
        DO KZ=0,N3+1
         RZ(KZ,I)=ZL+(KZ-1)*DZPA
        END DO
       END DO
      END DO

      WRITE(*,*) '     ===> TOTAL NUMBER OF PATCHES:',
     &           SUM(NPATCH(0:NL_MESH)),'<==='

      RETURN
      END

************************************************************************
      SUBROUTINE DENSITY(ITER,NX,NY,NZ,NL,NPATCH,PARE,PATCHNX,PATCHNY,
     &              PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,
     &              PATCHRZ,RXPA,RYPA,RZPA,MASAP,N_PARTICLES,N_DM,
     &              N_GAS,LADO0,T,ZETA,NPART_ESP,INTERP_DEGREE)
************************************************************************
*      Computes the density field on the AMR hierarchy (including base)
*       grid) by assigning each particle a cloud its size in the
*       initial conditions
************************************************************************

      IMPLICIT NONE
      INCLUDE 'input_files/asohf_parameters.dat'

*     function parameters
      INTEGER ITER,NX,NY,NZ,NL,N_PARTICLES,N_DM,N_GAS
      INTEGER NPATCH(0:NLEVELS),PARE(NPALEV)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
      REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED),
     &       MASAP(PARTIRED)
      REAL LADO0,T,ZETA
      INTEGER NPART_ESP(0:N_ESP-1),INTERP_DEGREE

*     VARIABLES
      REAL*4 U1(NMAX,NMAY,NMAZ)
      REAL*4 U1G(NMAX,NMAY,NMAZ)
      REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
      REAL*4 U11G(NAMRX,NAMRY,NAMRZ,NPALEV)
      COMMON /VARIA/ U1,U11,U1G,U11G

      REAL DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      REAL*4 RETE,HTE,ROTE
      COMMON /BACK/ RETE,HTE,ROTE

*     KERNELS
      INTEGER NMAX_KERNELS
      PARAMETER (NMAX_KERNELS=2**NLEVELS)
      REAL KERNELS(-NMAX_KERNELS:NMAX_KERNELS,0:NLEVELS)

*     other local variables
      INTEGER IPATCH,I,IX,JY,KZ,LOWP1,LOWP2,N,IP,II,JJ,KK,I1,J1,K1,IR
      INTEGER LOW1,LOW2,IRKERN,IRPART,N1,N2,N3,I2,J2,K2,I3,J3,K3,I4,J4
      INTEGER K4
      REAL BAS,BASX,BASY,BASZ,XL,YL,ZL,XP,YP,ZP,DENBAS,MAXKERNELEXT
      REAL DXPA,DYPA,DZPA,XR,YR,ZR,XLFIX,YLFIX,ZLFIX
      REAL,ALLOCATABLE::KERN1D(:)

      WRITE(*,*) '==== Density interpolation'

      WRITE(*,*) 'Computing 1D kernels of degree',INTERP_DEGREE,'...'
      CALL COMPUTE_KERNELS(KERNELS,NL,INTERP_DEGREE)

*     BASE GRID
      IR=0

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U1), PRIVATE(IX,JY,KZ), DEFAULT(NONE)
      DO KZ=1,NZ
      DO JY=1,NY
      DO IX=1,NX
       U1(IX,JY,KZ)=0.0
      END DO
      END DO
      END DO

      IRKERN=0
      N=2**IRKERN !1
      ALLOCATE(KERN1D(-N:N))
      DO I=-N,N
       KERN1D(I)=KERNELS(I,IR)
      END DO

      XL=-LADO0/2.0
      YL=-LADO0/2.0
      ZL=-LADO0/2.0

      LOWP1=1
      LOWP2=SUM(NPART_ESP)

!$OMP PARALLEL DO SHARED(LOWP1,LOWP2,MASAP,XL,YL,ZL,DX,DY,DZ,RXPA,RYPA,
!$OMP+                   RZPA,N,NX,NY,NZ,KERN1D),
!$OMP+            PRIVATE(IP,BAS,IX,JY,KZ,XP,YP,ZP,II,JJ,KK,I1,J1,K1),
!$OMP+            REDUCTION(+:U1)
!$OMP+            DEFAULT(NONE)
      DO IP=LOWP1,LOWP2
       BAS=MASAP(IP)

       XP=RXPA(IP)
       YP=RYPA(IP)
       ZP=RZPA(IP)
       IX=INT((XP-XL)/DX)+1
       JY=INT((YP-YL)/DY)+1
       KZ=INT((ZP-ZL)/DZ)+1

       DO KK=-N,N
       DO JJ=-N,N
       DO II=-N,N
        I1=IX+II
        J1=JY+JJ
        K1=KZ+KK
        IF (I1.GT.NX) I1=I1-NX
        IF (I1.LT.1) I1=I1+NX
        IF (J1.GT.NY) J1=J1-NY
        IF (J1.LT.1) J1=J1+NY
        IF (K1.GT.NZ) K1=K1-NZ
        IF (K1.LT.1) K1=K1+NZ

        U1(I1,J1,K1)=U1(I1,J1,K1)+BAS*KERN1D(II)*KERN1D(JJ)*KERN1D(KK)
       END DO
       END DO
       END DO

      END DO !IP=LOWP1,LOWP2

      DEALLOCATE(KERN1D)

      DENBAS=DX*DY*DZ*ROTE*RETE**3
!$OMP PARALLEL DO SHARED(NX,NY,NZ,U1,DENBAS),
!$OMP+            PRIVATE(IX,JY,KZ), DEFAULT(NONE)
      DO KZ=1,NZ
      DO JY=1,NY
      DO IX=1,NX
       U1(IX,JY,KZ)=U1(IX,JY,KZ)/DENBAS
      END DO
      END DO
      END DO

      WRITE(*,*) 'At level',0,MINVAL(U1),MAXVAL(U1)

      MAXKERNELEXT=DX

      DO IR=1,NL
       DXPA=DX/2.0**IR
       DYPA=DY/2.0**IR
       DZPA=DZ/2.0**IR

       DENBAS=DXPA*DYPA*DZPA*ROTE*RETE**3

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))

!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,U11,PATCHRX,
!$OMP+                   PATCHRY,PATCHRZ,DXPA,DYPA,DZPA,MAXKERNELEXT,
!$OMP+                   IR,KERNELS,NPART_ESP,RXPA,RYPA,RZPA,MASAP,
!$OMP+                   DENBAS),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,IX,JY,KZ,XLFIX,YLFIX,ZLFIX,
!$OMP+                    XL,YL,ZL,XR,YR,ZR,IRPART,IRKERN,N,KERN1D,
!$OMP+                    LOWP1,LOWP2,IP,XP,YP,ZP,BAS,I1,I2,J1,J2,K1,K2,
!$OMP+                    II,JJ,KK,I3,J3,K3),
!$OMP+            DEFAULT(NONE), SCHEDULE(DYNAMIC)
       DO IPATCH=LOW1,LOW2
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)

        DO KZ=1,N3
        DO JY=1,N2
        DO IX=1,N1
         U11(IX,JY,KZ,IPATCH)=0.0
        END DO
        END DO
        END DO

        XLFIX=PATCHRX(IPATCH)-DXPA
        YLFIX=PATCHRY(IPATCH)-DYPA
        ZLFIX=PATCHRZ(IPATCH)-DZPA

        XL=XLFIX
        YL=YLFIX
        ZL=ZLFIX
        XR=XL+N1*DXPA
        YR=YL+N2*DYPA
        ZR=ZL+N3*DZPA
        XL=XL-MAXKERNELEXT
        XR=XR+MAXKERNELEXT
        YL=YL-MAXKERNELEXT
        YR=YR+MAXKERNELEXT
        ZL=ZL-MAXKERNELEXT
        ZR=ZR+MAXKERNELEXT

        DO IRPART=0,N_ESP-1
         IRKERN=IR-IRPART
         IF (IRKERN.LT.0) IRKERN=0 ! SIZE (level) OF THE KERNEL TO USE
         N=2**IRKERN
         ALLOCATE(KERN1D(-N:N))
         DO I=-N,N
          KERN1D(I)=KERNELS(I,IRKERN)
         END DO
         !write(*,*) ir,ipatch,irpart,irkern,n
         LOWP1=SUM(NPART_ESP(0:IRPART-1))+1
         LOWP2=SUM(NPART_ESP(0:IRPART))

         DO IP=LOWP1,LOWP2
          XP=RXPA(IP)
          IF (XL.LT.XP) THEN
          IF (XP.LT.XR) THEN
          YP=RYPA(IP)
          IF (YL.LT.YP) THEN
          IF (YP.LT.YR) THEN
          ZP=RZPA(IP)
          IF (ZL.LT.ZP) THEN
          IF (ZP.LT.ZR) THEN
           BAS=MASAP(IP)

           IX=FLOOR((XP-XLFIX)/DXPA)+1
           JY=FLOOR((YP-YLFIX)/DYPA)+1
           KZ=FLOOR((ZP-ZLFIX)/DZPA)+1

!          Cloud bounds over the patch grid
           I1=IX-N
           I2=IX+N
           J1=JY-N
           J2=JY+N
           K1=KZ-N
           K2=KZ+N

!          Cycle if the cloud does not overlap the patch
           IF (I1.GT.N1.OR.I2.LT.1.OR.
     &         J1.GT.N2.OR.J2.LT.1.OR.
     &         K1.GT.N3.OR.J2.LT.1) CYCLE

!          Prune parts of the cloud outside the patch
           IF (I1.LT.1) I1=1
           IF (I2.GT.N1) I2=N1
           IF (J1.LT.1) J1=1
           IF (J2.GT.N2) J2=N2
           IF (K1.LT.1) K1=1
           IF (K2.GT.N3) K2=N3

           DO KK=K1,K2
           DO JJ=J1,J2
           DO II=I1,I2
            ! II,JJ,KK: patch cell indices
            ! I3,J3,K3: kernel indices
            I3=II-IX
            J3=JJ-JY
            K3=KK-KZ

            U11(II,JJ,KK,IPATCH)=U11(II,JJ,KK,IPATCH)+
     &                         BAS*KERN1D(I3)*KERN1D(J3)*KERN1D(K3)

           END DO
           END DO
           END DO

          END IF
          END IF
          END IF
          END IF
          END IF
          END IF

         END DO !IP=LOWP1,LOWP2
         DEALLOCATE(KERN1D)
        END DO !IRPART=0,N_ESP-1

        DO KZ=1,N3
        DO JY=1,N2
        DO IX=1,N1
         U11(IX,JY,KZ,IPATCH)=U11(IX,JY,KZ,IPATCH)/DENBAS
        END DO
        END DO
        END DO

       END DO !IPATCH=LOW1,LOW2

       bas=1000000.0
       do i=low1,low2
        n1=patchnx(i)
        n2=patchny(i)
        n3=patchnz(i)
        bas=min(bas,minval(u11(1:n1,1:n2,1:n3,i)))
       end do
       WRITE(*,*) 'At level',IR,bas,
     &                          maxval(u11(:,:,:,low1:low2))

      END DO

      WRITE(*,*) '==== End density interpolation'

      RETURN
      END

************************************************************************
      SUBROUTINE COMPUTE_KERNELS(KERNELS,NL,INTERP_DEGREE)
************************************************************************
*     Computes the kernel for density interpolation
************************************************************************

      IMPLICIT NONE
      INCLUDE 'input_files/asohf_parameters.dat'

      INTEGER NMAX_KERNELS
      PARAMETER (NMAX_KERNELS=2**NLEVELS)

*     function parameters
      REAL KERNELS(-NMAX_KERNELS:NMAX_KERNELS,0:NLEVELS)
      INTEGER NL,INTERP_DEGREE

      INTEGER IR,N,I,J,K
      REAL H,HDIV2,H32,BAS

      DO IR=0,NL
       N=2**IR

       H32=FLOAT(N)+0.5
       H=2.0*H32/3.0
       HDIV2=H/2.0

       IF (INTERP_DEGREE.EQ.2) THEN ! TSC-like KERNEL
!$OMP PARALLEL DO SHARED(N,H,KERNELS,IR),
!$OMP+            PRIVATE(I,BAS), DEFAULT(NONE)
        DO I=-N,N
         BAS=FLOAT(ABS(I))/H
         IF (BAS.LE.0.5) THEN
          KERNELS(I,IR)=0.75-BAS**2
         ELSE
          KERNELS(I,IR)=0.5*(1.5-BAS)**2
         END IF
        END DO
       ELSE IF (INTERP_DEGREE.EQ.1) THEN ! CIC-like KERNEL
!$OMP PARALLEL DO SHARED(N,H32,KERNELS,IR),
!$OMP+            PRIVATE(I,BAS), DEFAULT(NONE)
        DO I=-N,N
         BAS=FLOAT(ABS(I))/H32
         KERNELS(I,IR)=1-BAS
        END DO
       END IF ! INTERP_DEGREE

       BAS=0.0
!$OMP PARALLEL DO SHARED(N,KERNELS,IR), PRIVATE(I,J,K),
!$OMP+            REDUCTION(+:BAS), DEFAULT(NONE)
       DO I=-N,N
       DO J=-N,N
       DO K=-N,N
        BAS=BAS+KERNELS(I,IR)*KERNELS(J,IR)*KERNELS(K,IR)
       END DO
       END DO
       END DO

       !WRITE(*,*) IR,BAS,BAS**(1.0/3)

       BAS=1/BAS**(1.0/3.0)
!$OMP PARALLEL DO SHARED(N,KERNELS,BAS,IR), PRIVATE(I), DEFAULT(NONE)
       DO I=-N,N
        KERNELS(I,IR)=BAS*KERNELS(I,IR)
       END DO

       !WRITE(*,*) KERNELS(-N:N,IR)
       WRITE(*,*) IR,N,MINVAL(KERNELS(-N:N,IR)),MAXVAL(KERNELS(-N:N,IR))

      END DO !IR=0,NL

      RETURN
      END

************************************************************************
      SUBROUTINE INTERPOLATE_DENSITY_TSC(NX,NY,NZ,NL_MESH,NPATCH,
     &           PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &           PATCHRX,PATCHRY,PATCHRZ,RXPA,RYPA,RZPA,MASAP,
     &           N_PARTICLES,N_DM,LADO0,U1,U11)
************************************************************************
*     Interpolates density field (TSC)
*     Used only for solving Poisson equation
************************************************************************

      IMPLICIT NONE

      INCLUDE 'input_files/asohf_parameters.dat'

*     function parameters
      INTEGER ITER,NX,NY,NZ,NL_MESH,N_PARTICLES,N_DM
      INTEGER NPATCH(0:NLEVELS),PARE(NPALEV)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
      REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED),
     &       MASAP(PARTIRED)
      REAL LADO0
!     Not the COMMON one !!!!!!!!
      REAL*4 U1(NMAX,NMAY,NMAZ)
      REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!

*     COMMON VARIABLES
      REAL DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      REAL  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
      COMMON /GRID/  RADX,RADY,RADZ

      REAL*4 RETE,HTE,ROTE
      COMMON /BACK/ RETE,HTE,ROTE

*     LOCAL VARIABLES
      INTEGER PLEV(PARTIRED) ! Maximum mesh level of a particle.
      REAL XL,YL,ZL,DXPA,DYPA,DZPA,XR,YR,ZR,XP,YP,ZP,BAS,DENBAS
      INTEGER I,IX,JY,KZ,II,JJ,KK,N1,N2,N3,I1,I2,J1,J2,K1,K2,IR,IPATCH
      INTEGER BUF,LOW1,LOW2,NBAS,LOWP1,LOWP2
      REAL,ALLOCATABLE::RXPA2(:),RYPA2(:),RZPA2(:),MASAP2(:)
      REAL MINIRADX(-2:NAMRX+3),MINIRADY(-2:NAMRY+3),
     &     MINIRADZ(-2:NAMRZ+3)
      REAL VX(-1:1),VY(-1:1),VZ(-1:1)

      BUF=1 ! (1 cell extra buffer (TSC))

!$OMP PARALLEL DO SHARED(N_PARTICLES,PLEV), PRIVATE(I), DEFAULT(NONE)
      DO I=1,N_PARTICLES
       PLEV(I)=0
      END DO

      DO IR=NL_MESH,1,-1
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DXPA=DX/2.0**IR
       DYPA=DY/2.0**IR
       DZPA=DZ/2.0**IR
!$OMP PARALLEL DO SHARED(N_PARTICLES,PLEV,LOW1,LOW2,BUF,DXPA,DYPA,DZPA,
!$OMP+                   PATCHRX,PATCHRY,PATCHRZ,PATCHNX,PATCHNY,
!$OMP+                   PATCHNZ,RXPA,RYPA,RZPA,IR),
!$OMP+            PRIVATE(I,IPATCH,XL,YL,ZL,XR,YR,ZR),
!$OMP+            DEFAULT(NONE)
       DO I=1,N_PARTICLES
        IF (PLEV(I).EQ.0) THEN
         loop_patches: DO IPATCH=LOW1,LOW2
          XL=PATCHRX(IPATCH)-(1+BUF)*DXPA
          YL=PATCHRY(IPATCH)-(1+BUF)*DYPA
          ZL=PATCHRZ(IPATCH)-(1+BUF)*DZPA
          XR=XL+(PATCHNX(IPATCH)+2*BUF)*DXPA
          YR=YL+(PATCHNY(IPATCH)+2*BUF)*DYPA
          ZR=ZL+(PATCHNZ(IPATCH)+2*BUF)*DZPA
          IF (XL.LT.RXPA(I).AND.RXPA(I).LT.XR) THEN
           IF (YL.LT.RYPA(I).AND.RYPA(I).LT.YR) THEN
            IF (ZL.LT.RZPA(I).AND.RZPA(I).LT.ZR) THEN
             PLEV(I)=IR
             EXIT loop_patches
            END IF
           END IF
          END IF
         END DO loop_patches
        END IF
       END DO
      END DO

      DO IR=0,NL_MESH
       WRITE(*,*) 'Particles at level',IR,
     &             COUNT(PLEV(1:N_PARTICLES).EQ.IR)
      END DO

      ! NOW, BASE GRID.
      XL=-LADO0/2.0
      YL=-LADO0/2.0
      ZL=-LADO0/2.0

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U1), PRIVATE(IX,JY,KZ), DEFAULT(NONE)
      DO KZ=1,NZ
      DO JY=1,NY
      DO IX=1,NX
       U1(IX,JY,KZ)=0.0
      END DO
      END DO
      END DO

!$OMP PARALLEL DO SHARED(N_PARTICLES,RXPA,RYPA,RZPA,DX,DY,DZ,XL,YL,ZL,
!$OMP+                   RADX,RADY,RADZ,NX,NY,NZ,MASAP),
!$OMP+            PRIVATE(I,XP,YP,ZP,IX,JY,KZ,BAS,VX,VY,VZ,II,JJ,KK,
!$OMP+                    I1,J1,K1),
!$OMP+            REDUCTION(+:U1)
!$OMP+            DEFAULT(NONE)
      DO I=1,N_PARTICLES
       XP=RXPA(I)
       YP=RYPA(I)
       ZP=RZPA(I)

       IX=INT((XP-XL)/DX)+1
       JY=INT((YP-YL)/DY)+1
       KZ=INT((ZP-ZL)/DZ)+1
       IF (IX.LT.1) IX=1
       IF (IX.GT.NX) IX=NX
       IF (JY.LT.1) JY=1
       IF (JY.GT.NY) JY=NY
       IF (KZ.LT.1) KZ=1
       IF (KZ.GT.NZ) KZ=NZ

       BAS=ABS(XP-RADX(IX-1))/DX
       VX(-1)=0.5*(1.5-BAS)**2
       BAS=ABS(XP-RADX(IX))/DX
       VX(0)=0.75-BAS**2
       BAS=ABS(XP-RADX(IX+1))/DX
       VX(1)=0.5*(1.5-BAS)**2

       BAS=ABS(YP-RADY(JY-1))/DY
       VY(-1)=0.5*(1.5-BAS)**2
       BAS=ABS(YP-RADY(JY))/DY
       VY(0)=0.75-BAS**2
       BAS=ABS(YP-RADY(JY+1))/DY
       VY(1)=0.5*(1.5-BAS)**2

       BAS=ABS(ZP-RADZ(KZ-1))/DZ
       VZ(-1)=0.5*(1.5-BAS)**2
       BAS=ABS(ZP-RADZ(KZ))/DZ
       VZ(0)=0.75-BAS**2
       BAS=ABS(ZP-RADZ(KZ+1))/DZ
       VZ(1)=0.5*(1.5-BAS)**2

       DO KK=-1,1
       DO JJ=-1,1
       DO II=-1,1
        I1=IX+II
        J1=JY+JJ
        K1=KZ+KK
        IF (I1.EQ.NX+1) I1=1
        IF (I1.EQ.0) I1=NX
        IF (J1.EQ.NY+1) J1=1
        IF (J1.EQ.0) J1=NY
        IF (K1.EQ.NZ+1) K1=1
        IF (K1.EQ.0) K1=NZ

        U1(I1,J1,K1)= U1(I1,J1,K1) + MASAP(I)*VX(II)*VY(JJ)*VZ(KK)
       END DO
       END DO
       END DO

      END DO

      DENBAS=DX*DY*DZ*ROTE*RETE**3
!$OMP PARALLEL DO SHARED(NX,NY,NZ,DENBAS,U1), PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
      DO KZ=1,NZ
      DO JY=1,NY
      DO IX=1,NX
       U1(IX,JY,KZ)=U1(IX,JY,KZ)/DENBAS
      END DO
      END DO
      END DO

      WRITE(*,*) 'At level',0,minval(u1(:,:,:)),maxval(u1(:,:,:))

      DO IR=1,NL_MESH
       NBAS=COUNT(PLEV(1:N_PARTICLES).GE.IR)
       ALLOCATE(RXPA2(NBAS),RYPA2(NBAS),RZPA2(NBAS),MASAP2(NBAS))

       II=0
       DO I=1,N_PARTICLES
        IF (PLEV(I).GE.IR) THEN
         II=II+1
         RXPA2(II)=RXPA(I)
         RYPA2(II)=RYPA(I)
         RZPA2(II)=RZPA(I)
         MASAP2(II)=MASAP(I)
        END IF
       END DO

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DXPA=DX/2.0**IR
       DYPA=DY/2.0**IR
       DZPA=DZ/2.0**IR
       DENBAS=DXPA*DYPA*DZPA*ROTE*RETE**3

!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,U11),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)
        DO KZ=1,N3
        DO JY=1,N2
        DO IX=1,N1
         U11(IX,JY,KZ,IPATCH)=0.0
        END DO
        END DO
        END DO
       END DO

!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHRX,PATCHRY,PATCHRZ,BUF,DXPA,
!$OMP+                   DYPA,DZPA,PATCHNX,PATCHNY,PATCHNZ,NBAS,
!$OMP+                   RXPA2,RYPA2,RZPA2,MASAP2,U11,IR),
!$OMP+            PRIVATE(IPATCH,XL,YL,ZL,XR,YR,ZR,N1,N2,N3,II,JJ,KK,
!$OMP+                    MINIRADX,MINIRADY,MINIRADZ,I,XP,YP,ZP,IX,
!$OMP+                    JY,KZ,BAS,VX,VY,VZ,I1,J1,K1),
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2
        XL=PATCHRX(IPATCH)-(1+BUF)*DXPA
        YL=PATCHRY(IPATCH)-(1+BUF)*DYPA
        ZL=PATCHRZ(IPATCH)-(1+BUF)*DZPA
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)
        XR=XL+(N1+2*BUF)*DXPA
        YR=YL+(N2+2*BUF)*DYPA
        ZR=ZL+(N3+2*BUF)*DZPA

        DO II=-2,N1+3
         MINIRADX(II)=PATCHRX(IPATCH)+(FLOAT(II)-1.5)*DXPA
        END DO
        DO JJ=-2,N2+3
         MINIRADY(JJ)=PATCHRY(IPATCH)+(FLOAT(JJ)-1.5)*DYPA
        END DO
        DO KK=-2,N3+3
         MINIRADZ(KK)=PATCHRZ(IPATCH)+(FLOAT(KK)-1.5)*DZPA
        END DO

        DO I=1,NBAS
         XP=RXPA2(I)
         YP=RYPA2(I)
         ZP=RZPA2(I)
         IF (XL.LT.XP.AND.XP.LT.XR) THEN
          IF (YL.LT.YP.AND.YP.LT.YR) THEN
           IF (ZL.LT.ZP.AND.ZP.LT.ZR) THEN
            IX=INT((XP-XL)/DXPA)+1-BUF ! -BUF to account for the buffer
            JY=INT((YP-YL)/DYPA)+1-BUF
            KZ=INT((ZP-ZL)/DZPA)+1-BUF

            BAS=ABS(XP-MINIRADX(IX-1))/DXPA
            VX(-1)=0.5*(1.5-BAS)**2
            BAS=ABS(XP-MINIRADX(IX))/DXPA
            VX(0)=0.75-BAS**2
            BAS=ABS(XP-MINIRADX(IX+1))/DXPA
            VX(1)=0.5*(1.5-BAS)**2

            BAS=ABS(YP-MINIRADY(JY-1))/DYPA
            VY(-1)=0.5*(1.5-BAS)**2
            BAS=ABS(YP-MINIRADY(JY))/DYPA
            VY(0)=0.75-BAS**2
            BAS=ABS(YP-MINIRADY(JY+1))/DYPA
            VY(1)=0.5*(1.5-BAS)**2

            BAS=ABS(ZP-MINIRADZ(KZ-1))/DZPA
            VZ(-1)=0.5*(1.5-BAS)**2
            BAS=ABS(ZP-MINIRADZ(KZ))/DZPA
            VZ(0)=0.75-BAS**2
            BAS=ABS(ZP-MINIRADZ(KZ+1))/DZPA
            VZ(1)=0.5*(1.5-BAS)**2

            DO KK=-1,1
            DO JJ=-1,1
            DO II=-1,1
             I1=IX+II
             J1=JY+JJ
             K1=KZ+KK
             IF (I1.GT.0.AND.I1.LE.N1.AND.
     &           J1.GT.0.AND.J1.LE.N2.AND.
     &           K1.GT.0.AND.K1.LE.N3) THEN
              U11(I1,J1,K1,IPATCH)=U11(I1,J1,K1,IPATCH)
     &                             +MASAP2(I)*VX(II)*VY(JJ)*VZ(KK)
             END IF
            END DO
            END DO
            END DO

           END IF
          END IF
         END IF
        END DO
       END DO

!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,DENBAS,U11),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)
        DO KZ=1,N3
        DO JY=1,N2
        DO IX=1,N1
         U11(IX,JY,KZ,IPATCH)=U11(IX,JY,KZ,IPATCH)/DENBAS
        END DO
        END DO
        END DO
       END DO

       DEALLOCATE(RXPA2,RYPA2,RZPA2,MASAP2)

       bas=1000.0
       do i=low1,low2
        n1=patchnx(i)
        n2=patchny(i)
        n3=patchnz(i)
        bas=min(bas,minval(u11(1:n1,1:n2,1:n3,i)))
       end do
       WRITE(*,*) 'At level',IR,bas,
     &                          maxval(u11(:,:,:,low1:low2))
       END DO !ir=1,NL_MESH

      RETURN
      END


************************************************************************
      SUBROUTINE INTERPOLATE_DENSITY_KERNEL(ITER,NX,NY,NZ,NL_TSC,
     &           NL_MESH,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,
     &           PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,RXPA,RYPA,RZPA,
     &           MASAP,N_PARTICLES,N_DM,N_GAS,LADO0,T,ZETA,NPART_ESP)
************************************************************************
*     Interpolates density field (assuring the field will be continuous)
*     Old (deprecated now)
************************************************************************

      IMPLICIT NONE

      INCLUDE 'input_files/asohf_parameters.dat'

*     function parameters
      INTEGER ITER,NX,NY,NZ,NL_TSC,NL_MESH,N_PARTICLES,N_DM,N_GAS
      INTEGER NPATCH(0:NLEVELS),PARE(NPALEV)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
      REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED),
     &       MASAP(PARTIRED)
      REAL LADO0,T,ZETA
      INTEGER NPART_ESP(0:N_ESP-1)

*     COMMON VARIABLES
      REAL DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      REAL  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
      COMMON /GRID/  RADX,RADY,RADZ

      REAL*4 RETE,HTE,ROTE
      COMMON /BACK/ RETE,HTE,ROTE

      REAL*4 U1(NMAX,NMAY,NMAZ)
      REAL*4 U1G(NMAX,NMAY,NMAZ)
      REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
      REAL*4 U11G(NAMRX,NAMRY,NAMRZ,NPALEV)
      COMMON /VARIA/ U1,U11,U1G,U11G

*     LOCAL VARIABLES
      INTEGER PLEV(PARTIRED)
      REAL XL,YL,ZL,DXPA,DYPA,DZPA,XR,YR,ZR,XP,YP,ZP,BAS,DENBAS,PI
      REAL RRR,MMM,VOLCORRECT
      INTEGER I,IX,JY,KZ,II,JJ,KK,N1,N2,N3,I1,I2,J1,J2,K1,K2,IR,IPATCH
      INTEGER BUF,LOW1,LOW2,NBAS,NBAS2,LINT,NBASPART,NMIN_PART,IBASPART
      INTEGER IXX,JYY,KZZ,L1,L2,L3,LL1,LL2,LL3,MAXBUF
      INTEGER,ALLOCATABLE::SCRPAINT(:),SCRIX(:),SCRJY(:),SCRKZ(:)
      INTEGER,ALLOCATABLE::INDICE(:)
      INTEGER,ALLOCATABLE::SCRINT(:,:,:)
      REAL,ALLOCATABLE::MASAP2(:),RADPART(:),SCR41(:),SCR42(:)
      REAL MINIRADX(-2:NAMRX+3),MINIRADY(-2:NAMRY+3),
     &     MINIRADZ(-2:NAMRZ+3)

      REAL CPUTIME1,CPUTIME2
      INTEGER WALLTIME1,WALLTIME2,TIME

      BUF=12 ! extra buffer around each patch, to find radius around each cell
      NMIN_PART=8
      PI=DACOS(-1.D0)
      VOLCORRECT=6.0/PI !volume relation from cube to sphere (same side than diameter)

*     Find the (maximum) level a particle belongs to
*      (this is done to reduce computational burden)
!$OMP PARALLEL DO SHARED(N_PARTICLES,PLEV), PRIVATE(I), DEFAULT(NONE)
      DO I=1,N_PARTICLES
       PLEV(I)=0
      END DO

      DO IR=NL_MESH,NL_TSC+1,-1
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DXPA=DX/2.0**IR
       DYPA=DY/2.0**IR
       DZPA=DZ/2.0**IR
!$OMP PARALLEL DO SHARED(N_PARTICLES,PLEV,LOW1,LOW2,BUF,DXPA,DYPA,DZPA,
!$OMP+                   PATCHRX,PATCHRY,PATCHRZ,PATCHNX,PATCHNY,
!$OMP+                   PATCHNZ,RXPA,RYPA,RZPA,IR),
!$OMP+            PRIVATE(I,IPATCH,XL,YL,ZL,XR,YR,ZR),
!$OMP+            DEFAULT(NONE)
       DO I=1,N_PARTICLES
        IF (PLEV(I).EQ.0) THEN
         loop_patches: DO IPATCH=LOW1,LOW2
          XL=PATCHRX(IPATCH)-(1+BUF)*DXPA
          YL=PATCHRY(IPATCH)-(1+BUF)*DYPA
          ZL=PATCHRZ(IPATCH)-(1+BUF)*DZPA
          XR=XL+(PATCHNX(IPATCH)+2*BUF)*DXPA
          YR=YL+(PATCHNY(IPATCH)+2*BUF)*DYPA
          ZR=ZL+(PATCHNZ(IPATCH)+2*BUF)*DZPA
          IF (XL.LT.RXPA(I).AND.RXPA(I).LT.XR) THEN
           IF (YL.LT.RYPA(I).AND.RYPA(I).LT.YR) THEN
            IF (ZL.LT.RZPA(I).AND.RZPA(I).LT.ZR) THEN
             PLEV(I)=IR
             EXIT loop_patches
            END IF
           END IF
          END IF
         END DO loop_patches
        END IF
       END DO
      END DO

      DO IR=NL_TSC+1,NL_MESH
       WRITE(*,*) 'Particles at level',IR,
     &             COUNT(PLEV(1:N_PARTICLES).EQ.IR)
      END DO

*     Go from the finest to the coarsest level
      DO IR=NL_MESH,NL_TSC+1,-1

       CALL CPU_TIME(CPUTIME1)
       WALLTIME1=TIME()

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DXPA=DX/2.0**IR
       DYPA=DY/2.0**IR
       DZPA=DZ/2.0**IR
       DENBAS=(4.0*PI/3.0)*ROTE*RETE**3

       NBAS=COUNT(PLEV(1:N_PARTICLES).GE.IR)
       !WRITE(*,*) 'NBAS IN IR',IR,NBAS

!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,U11),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)
        DO KZ=1,N3
        DO JY=1,N2
        DO IX=1,N1
         U11(IX,JY,KZ,IPATCH)=0.0
        END DO
        END DO
        END DO
       END DO

!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHRX,PATCHRY,PATCHRZ,DXPA,DYPA,
!$OMP+                   DZPA,BUF,PATCHNX,PATCHNY,PATCHNZ,NBAS,
!$OMP+                   N_PARTICLES,PLEV,IR,RXPA,RYPA,RZPA,
!$OMP+                   NMIN_PART,DENBAS,U11,MASAP,VOLCORRECT),
!$OMP+            PRIVATE(IPATCH,XL,YL,ZL,N1,N2,N3,XR,YR,ZR,MINIRADX,
!$OMP+                    MINIRADY,MINIRADZ,II,JJ,KK,SCRPAINT,I,XP,YP,
!$OMP+                    ZP,SCRINT,IX,JY,KZ,NBAS2,SCRIX,SCRJY,SCRKZ,
!$OMP+                    LINT,NBASPART,I1,I2,J1,J2,K1,K2,MASAP2,
!$OMP+                    RADPART,IBASPART,IXX,JYY,KZZ,BAS,INDICE,
!$OMP+                    SCR41,SCR42,RRR,MMM,L1,L2,L3,LL1,LL2,LL3),
!$OMP+            DEFAULT(NONE),
!$OMP+            REDUCTION(MAX:MAXBUF),
!$OMP+            SCHEDULE(DYNAMIC)
       DO IPATCH=LOW1,LOW2
        XL=PATCHRX(IPATCH)-(1+BUF)*DXPA
        YL=PATCHRY(IPATCH)-(1+BUF)*DYPA
        ZL=PATCHRZ(IPATCH)-(1+BUF)*DZPA
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)
        XR=XL+(N1+2*BUF)*DXPA
        YR=YL+(N2+2*BUF)*DYPA
        ZR=ZL+(N3+2*BUF)*DZPA

        DO II=-2,N1+3
         MINIRADX(II)=PATCHRX(IPATCH)+(FLOAT(II)-1.5)*DXPA
        END DO
        DO JJ=-2,N2+3
         MINIRADY(JJ)=PATCHRY(IPATCH)+(FLOAT(JJ)-1.5)*DYPA
        END DO
        DO KK=-2,N3+3
         MINIRADZ(KK)=PATCHRZ(IPATCH)+(FLOAT(KK)-1.5)*DZPA
        END DO

        ALLOCATE(SCRPAINT(NBAS))
        II=0
        DO I=1,N_PARTICLES
         IF (PLEV(I).GE.IR) THEN
          XP=RXPA(I)
          YP=RYPA(I)
          ZP=RZPA(I)
          IF (XL.LT.XP) THEN
          IF (XP.LT.XR) THEN
          IF (YL.LT.YP) THEN
          IF (YP.LT.YR) THEN
          IF (ZL.LT.ZP) THEN
          IF (ZP.LT.ZR) THEN
           II=II+1
           SCRPAINT(II)=I
          END IF
          END IF
          END IF
          END IF
          END IF
          END IF
         END IF
        END DO

        NBAS2=II
        !WRITE(*,*) 'IN PATCH',IPATCH,'NBAS2=',NBAS2,xl,xr,yl,yr,zl,zr
        ALLOCATE(SCRIX(NBAS2),SCRJY(NBAS2),SCRKZ(NBAS2))
        ALLOCATE(SCRINT(1-BUF:N1+BUF,1-BUF:N2+BUF,1-BUF:N3+BUF))

        DO KZ=1-BUF,N3+BUF
        DO JY=1-BUF,N2+BUF
        DO IX=1-BUF,N1+BUF
         SCRINT(IX,JY,KZ)=0
        END DO
        END DO
        END DO

        L1=1-BUF
        L2=1-BUF
        L3=1-BUF
        LL1=N1+BUF
        LL2=N2+BUF
        LL3=N3+BUF
        DO II=1,NBAS2
         I=SCRPAINT(II)

         XP=RXPA(I)
         YP=RYPA(I)
         ZP=RZPA(I)
         IX=INT((XP-XL)/DXPA)+L1
         JY=INT((YP-YL)/DYPA)+L2
         KZ=INT((ZP-ZL)/DZPA)+L3
         IF (IX.EQ.L1-1) THEN
          IX=L1
         ELSE IF (IX.EQ.LL1+1) THEN
          IX=LL1
         END IF
         IF (JY.EQ.L2-1) THEN
          JY=L2
         ELSE IF (JY.EQ.LL2+1) THEN
          JY=LL2
         END IF
         IF (KZ.EQ.L3-1) THEN
          KZ=L3
         ELSE IF (KZ.EQ.LL3+1) THEN
          KZ=LL3
         END IF

         SCRIX(II)=IX
         SCRJY(II)=JY
         SCRKZ(II)=KZ

         SCRINT(IX,JY,KZ)=SCRINT(IX,JY,KZ)+1
        END DO
        !write(*,*) ipatch,nbas2,sum(scrint)

        DO KZ=1,N3
        DO JY=1,N2
        DO IX=1,N1
          LINT=0
          NBASPART=0

          XL=MINIRADX(IX)
          YL=MINIRADY(JY)
          ZL=MINIRADZ(KZ)

          DO WHILE (NBASPART.LT.NMIN_PART)
           LINT=LINT+1
           I1=IX-LINT
           I2=IX+LINT
           J1=JY-LINT
           J2=JY+LINT
           K1=KZ-LINT
           K2=KZ+LINT
           IF (I1.LE.L1.OR.I2.GE.LL1.OR.
     &         J1.LE.L2.OR.J2.GE.LL2.OR.
     &         K1.LE.L3.OR.K2.GE.LL3) THEN
             WRITE(*,*) 'increase BUF in INTERPOLATE_DENSITY_KERNEL',
     &                  ipatch,nbaspart,ix,jy,kz,i1,i2,j1,j2,k1,k2,
     &                  n1,n2,n3
             STOP
            END IF
            NBASPART=SUM(SCRINT(I1:I2,J1:J2,K1:K2))
          END DO

          IF (LINT.GT.MAXBUF) MAXBUF=LINT

          ALLOCATE(MASAP2(NBASPART),RADPART(NBASPART))

          IBASPART=0
          DO II=1,NBAS2
           IXX=SCRIX(II)
           IF (IXX.GE.I1) THEN
           IF (IXX.LE.I2) THEN
           JYY=SCRJY(II)
           IF (JYY.GE.J1) THEN
           IF (JYY.LE.J2) THEN
           KZZ=SCRKZ(II)
           IF (KZZ.GE.K1) THEN
           IF (KZZ.LE.K2) THEN
            IBASPART=IBASPART+1
            I=SCRPAINT(II)
            XP=RXPA(I)
            YP=RYPA(I)
            ZP=RZPA(I)
            MASAP2(IBASPART)=MASAP(I)
            BAS=(XP-XL)**2+(YP-YL)**2+(ZP-ZL)**2
            RADPART(IBASPART)=SQRT(BAS)
           END IF
           END IF
           END IF
           END IF
           END IF
           END IF
          END DO

          IF (LINT.EQ.1) THEN
           RRR=1.5*DXPA
           ! to avoid increasing buffers:
           !  we account for the mass in the cube (instead of sphere)
           !  and correct the volumes
           MMM=SUM(MASAP2(1:NBASPART))/VOLCORRECT
          ELSE
           ALLOCATE(INDICE(NBASPART),SCR41(NBASPART),SCR42(NBASPART))
           CALL INDEXX(NBASPART,RADPART,INDICE)
           DO II=1,NBASPART
            SCR41(II)=RADPART(INDICE(II))
            SCR42(II)=MASAP2(INDICE(II))
           END DO
           DO II=1,NBASPART
            RADPART(II)=SCR41(II)
            MASAP2(II)=SCR42(II)
           END DO
           DEALLOCATE(INDICE,SCR41,SCR42)

           IF (NBASPART.GT.NMIN_PART) THEN
            RRR=0.5*(RADPART(NMIN_PART)+RADPART(NMIN_PART+1))
            MMM=SUM(MASAP2(1:NMIN_PART))
           ELSE
            RRR=RADPART(NBASPART)
            MMM=SUM(MASAP2(1:NBASPART))
           END IF
         END IF

         DEALLOCATE(MASAP2,RADPART)

         U11(IX,JY,KZ,IPATCH)=MMM/(DENBAS*RRR**3)
        END DO
        END DO
        END DO

        DEALLOCATE(SCRIX,SCRJY,SCRKZ,SCRINT,SCRPAINT)

       END DO !IPATCH

       bas=1000.0
       do i=low1,low2
        n1=patchnx(i)
        n2=patchny(i)
        n3=patchnz(i)
        bas=min(bas,minval(u11(1:n1,1:n2,1:n3,i)))
       end do

       CALL CPU_TIME(CPUTIME2)
       WALLTIME2=TIME()
       CPUTIME2=CPUTIME2-CPUTIME1
       WALLTIME2=WALLTIME2-WALLTIME1

       WRITE(*,*) 'At level',IR,bas,
     &                          maxval(u11(:,:,:,low1:low2))
       WRITE(*,*) '--> CPU, Wall time (s):', CPUTIME2,FLOAT(WALLTIME2)
       WRITE(*,*) '--> Max kernel size (cells):',MAXBUF

      END DO !IR=NL_MESH,NL_TSC+1,-1


      RETURN
      END

************************************************************************
      SUBROUTINE VEINSGRID_REDUCED(IR,NPATCH,PARE,
     &      PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &      PATCHRY,PATCHRZ,CR01,CONTA11,LOW1,LOW2)
************************************************************************
*     small fraction of patches with rare geometry were not detected
*     when overlapping

      IMPLICIT NONE

      INCLUDE 'input_files/asohf_parameters.dat'

      INTEGER NPALEV2

*     U11(PATCHNX,PATCHNY,PATCHNZ,NLEVEL,NPALEV)
*     PATCHNX,PATCHNY,PATCHNZ patches dimensions
*     IPATCH number of patches per level
*     NLEVELS total number of levels

      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV)
      INTEGER PATCHNY(NPALEV)
      INTEGER PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV)
      INTEGER PATCHY(NPALEV)
      INTEGER PATCHZ(NPALEV)
      REAL*4  PATCHRX(NPALEV)
      REAL*4  PATCHRY(NPALEV)
      REAL*4  PATCHRZ(NPALEV)
      INTEGER PARE(NPALEV)

      INTEGER CR1,CR2,CR3,CR4,CR5,CR6
      INTEGER IR,I,J,IX,JY,KZ,II,JJ,KK,I2,J2
      INTEGER N1,N2,N3,L1,L2,L3,NN1,NN2,NN3,LL1,LL2,LL3
      INTEGER KK2,JJ2,II2,KZ2,JY2,IX2
      INTEGER NV,A2,B2,C2,K

      INTEGER LOW1,LOW2
      INTEGER SOLAP_PATCH(NAMRX,NAMRY,NAMRZ,LOW1:LOW2)
      INTEGER CR01(NAMRX,NAMRY,NAMRZ,LOW1:LOW2)
      INTEGER CONTA11(NAMRX,NAMRY,NAMRZ,LOW1:LOW2)

      REAL*4 A1,B1,C1,RIV1,RIV2,RIV3
      INTEGER CONTROL
      INTEGER CORNX1,CORNXX1,CORNX2,CORNXX2
      INTEGER CORNY1,CORNYY1,CORNY2,CORNYY2
      INTEGER CORNZ1,CORNZZ1,CORNZ2,CORNZZ2
      REAL*4 RX1,RXX1,RX2,RXX2,RY1,RYY1,RY2,RYY2
      REAL*4 RZ1,RZZ1,RZ2,RZZ2,ORXX1,ORYY1,ORZZ1

      REAL*4 DXPA,DYPA,DZPA
      REAL*4 DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      REAL  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
      COMMON /GRID/  RADX,RADY,RADZ

      INTEGER,ALLOCATABLE::VECINO(:,:)
      INTEGER,ALLOCATABLE::NVECI(:)

      INTEGER IG1,IG2,JG1,JG2,KG1,KG2,IG3,JG3,KG3,IG4,JG4,KG4
      REAL*4 RXFIX,RYFIX,RZFIX

      REAL*4 OVERLAP, OVERLAP_THR
*
       NPALEV2=MAX(100,INT(NPALEV/5))
       OVERLAP_THR=0.1

       ALLOCATE(VECINO(NPALEV2,NPATCH(IR)))
       ALLOCATE(NVECI(NPATCH(IR)))

       DXPA=DX/(2.**IR)
       DYPA=DY/(2.**IR)
       DZPA=DZ/(2.**IR)

*      built auxiliar grid for comparison
       RXFIX=RADX(1) - DX*0.5 + 0.5*DXPA
       RYFIX=RADY(1) - DY*0.5 + 0.5*DYPA
       RZFIX=RADZ(1) - DZ*0.5 + 0.5*DZPA

!$OMP   PARALLEL DO SHARED(IR,NPATCH,PARE,PATCHX,PATCHY,PATCHZ,
!$OMP+        PATCHNX,PATCHNY,PATCHNZ,VECINO,NVECI,
!$OMP+        DXPA,DYPA,DZPA,PATCHRX,PATCHRY,PATCHRZ,
!$OMP+        SOLAP_PATCH,LOW1,LOW2,RXFIX,RYFIX,RZFIX),
!$OMP+  PRIVATE(I,N1,N2,N3,NV,J,NN1,NN2,NN3,
!$OMP+          RX1,RY1,RZ1,RXX1,RYY1,RZZ1,I2,
!$OMP+          IG1,IG2,JG1,JG2,KG1,KG2,IG3,JG3,KG3,IG4,JG4,KG4),
!$OMP+  DEFAULT(NONE)
       DO I=LOW1,LOW2

         I2=I-LOW1+1

         NVECI(I2)=0
         VECINO(:,I2)=0

         SOLAP_PATCH(:,:,:,I)=0

         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         NV=0

         RX1=PATCHRX(I)-0.5*DXPA
         RY1=PATCHRY(I)-0.5*DYPA
         RZ1=PATCHRZ(I)-0.5*DZPA

         IG1=INT(((RX1-RXFIX)/DXPA)+0.5) + 1
         JG1=INT(((RY1-RYFIX)/DYPA)+0.5) + 1
         KG1=INT(((RZ1-RZFIX)/DZPA)+0.5) + 1

         IG2=IG1 + N1 - 1
         JG2=JG1 + N2 - 1
         KG2=KG1 + N3 - 1

         DO J=LOW1,LOW2
          IF (J.NE.I) THEN

          NN1=PATCHNX(J)
          NN2=PATCHNY(J)
          NN3=PATCHNZ(J)

          RXX1=PATCHRX(J)-0.5*DXPA
          RYY1=PATCHRY(J)-0.5*DYPA
          RZZ1=PATCHRZ(J)-0.5*DZPA

          IG3=INT(((RXX1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RYY1-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RZZ1-RZFIX)/DZPA)+0.5) + 1

          IG4=IG3 + NN1 - 1
          JG4=JG3 + NN2 - 1
          KG4=KG3 + NN3 - 1

          IF (IG1.LE.IG4.AND.IG3.LE.IG2.AND.
     &        JG1.LE.JG4.AND.JG3.LE.JG2.AND.
     &        KG1.LE.KG4.AND.KG3.LE.KG2) THEN
           NV=NV+1
           VECINO(NV,I2)=J
          END IF

          END IF

         END DO
         NVECI(I2)=NV
       END DO


       IF (MAXVAL(NVECI(1:NPATCH(IR))).GT.NPALEV2) THEN
         WRITE(*,*) 'ERROR: gvecino ST second dimension too large',
     &     MAXVAL(NVECI(1:NPATCH(IR)))
         STOP
       END IF


       DO I=LOW1,LOW2

         L1=PATCHX(I)
         L2=PATCHY(I)
         L3=PATCHZ(I)

         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         RX1=PATCHRX(I)-0.5*DXPA
         RY1=PATCHRY(I)-0.5*DYPA
         RZ1=PATCHRZ(I)-0.5*DZPA
         RX2=PATCHRX(I)-0.5*DXPA+(N1-1)*DXPA
         RY2=PATCHRY(I)-0.5*DYPA+(N2-1)*DYPA
         RZ2=PATCHRZ(I)-0.5*DZPA+(N3-1)*DZPA

         I2=I-LOW1+1

         DO K=1,NVECI(I2)
         J=VECINO(K,I2)
         J2=J-LOW1+1

         LL1=PATCHX(J)
         LL2=PATCHY(J)
         LL3=PATCHZ(J)

         NN1=PATCHNX(J)
         NN2=PATCHNY(J)
         NN3=PATCHNZ(J)

         RXX1=PATCHRX(J)-0.5*DXPA
         RYY1=PATCHRY(J)-0.5*DYPA
         RZZ1=PATCHRZ(J)-0.5*DZPA
         RXX2=PATCHRX(J)-0.5*DXPA+(NN1-1)*DXPA
         RYY2=PATCHRY(J)-0.5*DYPA+(NN2-1)*DYPA
         RZZ2=PATCHRZ(J)-0.5*DZPA+(NN3-1)*DZPA

*        X
         IF (RXX1.GE.RX1.AND.RXX2.LE.RX2) THEN
            CORNX1=INT(((RXX1-RX1)/DXPA)+0.5) + 1
            CORNX2=INT(((RXX2-RX1)/DXPA)+0.5) + 1
            CORNXX1=1
            CORNXX2=NN1
         END IF
         IF (RXX1.GE.RX1.AND.RXX2.GT.RX2) THEN
            CORNX1=INT(((RXX1-RX1)/DXPA)+0.5) + 1
            CORNX2=N1
            CORNXX1=1
            CORNXX2=INT(((RX2-RXX1)/DXPA)+0.5) +1
         END IF
         IF (RXX2.LE.RX2.AND.RXX1.LT.RX1) THEN
            CORNX1=1
            CORNX2=INT(((RXX2-RX1)/DXPA)+0.5) + 1
            CORNXX1=INT(((RX1-RXX1)/DXPA)+0.5) + 1
            CORNXX2=NN1
         END IF
         IF (RXX1.LT.RX1.AND.RXX2.GT.RX2) THEN
            CORNX1=1
            CORNX2=N1
            CORNXX1=INT(((RX1-RXX1)/DXPA)+0.5) + 1
            CORNXX2=INT(((RX2-RXX1)/DXPA)+0.5) + 1
         END IF

*        Y
         IF (RYY1.GE.RY1.AND.RYY2.LE.RY2) THEN
            CORNY1=INT(((RYY1-RY1)/DYPA)+0.5) + 1
            CORNY2=INT(((RYY2-RY1)/DYPA)+0.5) + 1
            CORNYY1=1
            CORNYY2=NN2
         END IF
         IF (RYY1.GE.RY1.AND.RYY2.GT.RY2) THEN
            CORNY1=INT(((RYY1-RY1)/DYPA)+0.5) + 1
            CORNY2=N2
            CORNYY1=1
            CORNYY2=INT(((RY2-RYY1)/DYPA)+0.5) +1
         END IF
         IF (RYY2.LE.RY2.AND.RYY1.LT.RY1) THEN
            CORNY1=1
            CORNY2=INT(((RYY2-RY1)/DYPA)+0.5) + 1
            CORNYY1=INT(((RY1-RYY1)/DYPA)+0.5) + 1
            CORNYY2=NN2
         END IF
         IF (RYY1.LT.RY1.AND.RYY2.GT.RY2) THEN
            CORNY1=1
            CORNY2=N2
            CORNYY1=INT(((RY1-RYY1)/DYPA)+0.5) + 1
            CORNYY2=INT(((RY2-RYY1)/DYPA)+0.5) + 1
         END IF

*        Z
         IF (RZZ1.GE.RZ1.AND.RZZ2.LE.RZ2) THEN
            CORNZ1=INT(((RZZ1-RZ1)/DZPA)+0.5) + 1
            CORNZ2=INT(((RZZ2-RZ1)/DZPA)+0.5) + 1
            CORNZZ1=1
            CORNZZ2=NN3
         END IF
         IF (RZZ1.GE.RZ1.AND.RZZ2.GT.RZ2) THEN
            CORNZ1=INT(((RZZ1-RZ1)/DZPA)+0.5) + 1
            CORNZ2=N3
            CORNZZ1=1
            CORNZZ2=INT(((RZ2-RZZ1)/DZPA)+0.5) +1
         END IF
         IF (RZZ2.LE.RZ2.AND.RZZ1.LT.RZ1) THEN
            CORNZ1=1
            CORNZ2=INT(((RZZ2-RZ1)/DZPA)+0.5) + 1
            CORNZZ1=INT(((RZ1-RZZ1)/DZPA)+0.5) + 1
            CORNZZ2=NN3
         END IF
         IF (RZZ1.LT.RZ1.AND.RZZ2.GT.RZ2) THEN
            CORNZ1=1
            CORNZ2=N3
            CORNZZ1=INT(((RZ1-RZZ1)/DZPA)+0.5) + 1
            CORNZZ2=INT(((RZ2-RZZ1)/DZPA)+0.5) + 1
         END IF

*        overlap is the fraction of volumen overlapping
*        the idea is to fix a thershold of overlapping that below
*        this therhold for instance 10%, the cells are not marked as overlap
         OVERLAP=(CORNZZ2-CORNZZ1+1)*(CORNYY2-CORNYY1+1)*
     &           (CORNXX2-CORNXX1+1)
         OVERLAP=ABS(OVERLAP)/PATCHNX(J)/PATCHNY(J)/PATCHNZ(J)

**       celdas madre del nivel inferior
         IF (OVERLAP.GT.OVERLAP_THR) THEN
           DO KK=CORNZZ1,CORNZZ2
           DO JJ=CORNYY1,CORNYY2
           DO II=CORNXX1,CORNXX2
            IX=II-CORNXX1+CORNX1
            JY=JJ-CORNYY1+CORNY1
            KZ=KK-CORNZZ1+CORNZ1
            IF (SOLAP_PATCH(IX,JY,KZ,I).EQ.0) THEN
*            the cell ii,jj,kk,ir,j overlaps ix,jy,kz,ir,i
             SOLAP_PATCH(II,JJ,KK,J)=1
             CR01(II,JJ,KK,J)=0
             CONTA11(II,JJ,KK,J)=0
            END IF
           END DO
           END DO
           END DO
          END IF
       END DO
       END DO

      DEALLOCATE(VECINO)
      DEALLOCATE(NVECI)

      RETURN
      END

***********************************************************************
      SUBROUTINE COMPUTE_CR0AMR(NL,NX,NY,NZ,NPATCH,PARE,PATCHNX,PATCHNY,
     &                     PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                     PATCHRY,PATCHRZ,CR0AMR,CR0AMR1,LADO0)
***********************************************************************

       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NL,NX,NY,NZ
       INTEGER NPATCH(0:NLEVELS)
       INTEGER PARE(NPALEV)
       INTEGER PATCHNX(NPALEV)
       INTEGER PATCHNY(NPALEV)
       INTEGER PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV)
       INTEGER PATCHY(NPALEV)
       INTEGER PATCHZ(NPALEV)
       REAL*4  PATCHRX(NPALEV)
       REAL*4  PATCHRY(NPALEV)
       REAL*4  PATCHRZ(NPALEV)
       INTEGER CR0AMR(NMAX,NMAY,NMAZ)
       INTEGER CR0AMR1(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL*4 LADO0

       INTEGER IX,JY,KZ,IR,LOW1,LOW2,I,L1,L2,L3,CR1,CR2,CR3,MARCA,II
       INTEGER LOW3,LOW4
       REAL DXPA,DYPA,DZPA,DX,DY,DZ,XX,XX1,XX2,YY,YY1,YY2,ZZ,ZZ1,ZZ2
       REAL BAS1,BAS2,BAS3

!$OMP PARALLEL DO SHARED(NX,NY,NZ,CR0AMR),PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO KZ=1,NZ
       DO JY=1,NY
       DO IX=1,NX
        CR0AMR(IX,JY,KZ)=1
       END DO
       END DO
       END DO

       DO IR=1,NL
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(LOW1,LOW2,CR0AMR1),PRIVATE(I),DEFAULT(NONE)
        DO I=LOW1,LOW2
         CR0AMR1(:,:,:,i)=1
        END DO
       END DO

       IR=1
!$OMP PARALLEL DO SHARED(IR,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
!$OMP+                   PATCHX,PATCHY,PATCHZ),
!$OMP+            PRIVATE(I,L1,L2,L3,IX,JY,KZ,CR1,CR2,CR3),
!$OMP+            REDUCTION(MIN:CR0AMR),DEFAULT(NONE)
       DO I=1,NPATCH(IR)
        L1=PATCHX(I)
        L2=PATCHY(I)
        L3=PATCHZ(I)
        DO KZ=1,PATCHNZ(I),2
        DO JY=1,PATCHNY(I),2
        DO IX=1,PATCHNX(I),2
*        celdas madre del nivel inferior
         CR1=L1-1+INT((IX+1)/2)
         CR2=L2-1+INT((JY+1)/2)
         CR3=L3-1+INT((KZ+1)/2)
         CR0AMR(CR1,CR2,CR3)=0
        END DO
        END DO
        END DO
       END DO

       DO IR=2,NL
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
!$OMP  PARALLEL DO SHARED(IR,PATCHX,PATCHY,PATCHZ,PARE,
!$OMP+         PATCHNX,PATCHNY,PATCHNZ,LOW1,LOW2,CR0AMR1),
!$OMP+         PRIVATE(I,L1,L2,L3,IX,JY,KZ,CR1,CR2,CR3),DEFAULT(NONE)
        DO I=LOW1,LOW2
         L1=PATCHX(I)
         L2=PATCHY(I)
         L3=PATCHZ(I)
         DO KZ=1,PATCHNZ(I)
         DO JY=1,PATCHNY(I)
         DO IX=1,PATCHNX(I)
*        celdas madre del nivel inferior
          CR1=L1-1+INT((IX+1)/2)
          CR2=L2-1+INT((JY+1)/2)
          CR3=L3-1+INT((KZ+1)/2)
          CR0AMR1(CR1,CR2,CR3,PARE(I))=0
         END DO
         END DO
         END DO
        END DO
       END DO

       DX=LADO0/NX
       DY=LADO0/NY
       DZ=LADO0/NZ
       DO IR=1,NL-1
        DXPA=DX/(2.**IR)
        DYPA=DY/(2.**IR)
        DZPA=DZ/(2.**IR)
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(IR,NL,DXPA,DYPA,DZPA,DX,DY,DZ,LOW1,LOW2,
!$OMP+            NPATCH,PATCHNZ,PATCHNY,PATCHNX,CR0AMR1,
!$OMP+            PATCHRX,PATCHRY,PATCHRZ,PARE),
!$OMP+      PRIVATE(I,KZ,JY,IX,XX,YY,ZZ,LOW3,LOW4,II,
!$OMP+            XX1,YY1,ZZ1,XX2,YY2,ZZ2,BAS1,BAS2,BAS3,MARCA),
!$OMP+      DEFAULT(NONE)
        DO I=LOW1,LOW2
         DO KZ=1,PATCHNZ(I)
         DO JY=1,PATCHNY(I)
         DO IX=1,PATCHNX(I)
          IF (CR0AMR1(IX,JY,KZ,I).EQ.1) THEN
*         celdas madre del nivel inferior
           XX=PATCHRX(I)-0.5*DXPA+(IX-1)*DXPA
           YY=PATCHRY(I)-0.5*DYPA+(JY-1)*DYPA
           ZZ=PATCHRZ(I)-0.5*DZPA+(KZ-1)*DZPA
           MARCA=0

*      FILLS
           LOW3=SUM(NPATCH(0:IR))+1
           LOW4=SUM(NPATCH(0:IR+1))
           DO II=LOW3,LOW4
            IF (PARE(II).NE.I) THEN
             XX1=PATCHRX(II)-0.5*(DXPA*0.5)
             YY1=PATCHRY(II)-0.5*(DYPA*0.5)
             ZZ1=PATCHRZ(II)-0.5*(DZPA*0.5)
             XX2=XX1+(PATCHNX(II)-1)*DXPA*0.5
             YY2=YY1+(PATCHNY(II)-1)*DYPA*0.5
             ZZ2=ZZ1+(PATCHNZ(II)-1)*DZPA*0.5

             BAS1=(XX-XX1)*(XX2-XX)
             BAS2=(YY-YY1)*(YY2-YY)
             BAS3=(ZZ-ZZ1)*(ZZ2-ZZ)
             IF (BAS1.GE.0.0.AND.BAS2.GE.0.0.
     &           AND.BAS3.GE.0.0) THEN
              CR0AMR1(IX,JY,KZ,I)=0   !quiere decir que esta refinada esta celda
              MARCA=1
             END IF
            END IF
            IF (MARCA.EQ.1) EXIT
           END DO
          END IF
         END DO
         END DO
         END DO
        END DO
       END DO

      END

************************************************************************
        SUBROUTINE VEINSGRID(IR,NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &             PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,SOLAP,
     &             VECINO,NVECI)
************************************************************************
*       Finds the overlaps between patches in a given level of the
*       grid hierarchy
************************************************************************

       IMPLICIT NONE

       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NPATCH(0:NLEVELS)
       INTEGER PATCHNX(NPALEV)
       INTEGER PATCHNY(NPALEV)
       INTEGER PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV)
       INTEGER PATCHY(NPALEV)
       INTEGER PATCHZ(NPALEV)
       REAL*4  PATCHRX(NPALEV)
       REAL*4  PATCHRY(NPALEV)
       REAL*4  PATCHRZ(NPALEV)
       INTEGER PARE(NPALEV)

       INTEGER IR,I,J,IX,JY,KZ,II,JJ,KK
       INTEGER N1,N2,N3,L1,L2,L3,NL
       INTEGER NN1,NN2,NN3,LL1,LL2,LL3
       INTEGER KZ2,JY2,IX2,I2
       INTEGER NV,A2,B2,C2,K,LOW1,LOW2

       INTEGER VECINO(NPALEV,NPALEV),NVECI(NPALEV)
       INTEGER SOLAP(NAMRX,NAMRY,NAMRZ,NPALEV)

       INTEGER CORNX1,CORNXX1,CORNX2,CORNXX2
       INTEGER CORNY1,CORNYY1,CORNY2,CORNYY2
       INTEGER CORNZ1,CORNZZ1,CORNZ2,CORNZZ2
       REAL*4 RX1,RXX1,RX2,RXX2,RY1,RYY1,RY2,RYY2
       REAL*4 RZ1,RZZ1,RZ2,RZZ2

       REAL*4 DXPA,DYPA,DZPA
       REAL*4 DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       REAL*4  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
       COMMON /GRID/ RADX,RADY,RADZ

       INTEGER IG1,IG2,JG1,JG2,KG1,KG2,IG3,JG3,KG3,IG4,JG4,KG4
       REAL*4 RXFIX,RYFIX,RZFIX
       INTEGER NPALEV2

       NPALEV2=MAX(100,INT(NPALEV/10))

       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)

*      auxiliar grid for comparison
       RXFIX=RADX(1) - DX*0.5 + 0.5*DXPA
       RYFIX=RADY(1) - DY*0.5 + 0.5*DYPA
       RZFIX=RADZ(1) - DZ*0.5 + 0.5*DZPA

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))

*      identify neighbouring patches
       DO I=LOW1,LOW2
         I2=I-LOW1+1

         NVECI(I2)=0
         VECINO(:,I2)=0

         L1=PATCHX(I)
         L2=PATCHY(I)
         L3=PATCHZ(I)

         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         SOLAP(:,:,:,I)=1

         NV=0

         RX1=PATCHRX(I)-0.5*DXPA
         RY1=PATCHRY(I)-0.5*DYPA
         RZ1=PATCHRZ(I)-0.5*DZPA
         IG1=INT(((RX1-RXFIX)/DXPA)+0.5) + 1
         JG1=INT(((RY1-RYFIX)/DYPA)+0.5) + 1
         KG1=INT(((RZ1-RZFIX)/DZPA)+0.5) + 1

         IG2=IG1 + N1 - 1
         JG2=JG1 + N2 - 1
         KG2=KG1 + N3 - 1

         DO J=LOW1,LOW2
          IF (J.NE.I) THEN

           NN1=PATCHNX(J)
           NN2=PATCHNY(J)
           NN3=PATCHNZ(J)

           RXX1=PATCHRX(J)-0.5*DXPA
           RYY1=PATCHRY(J)-0.5*DYPA
           RZZ1=PATCHRZ(J)-0.5*DZPA

           IG3=INT(((RXX1-RXFIX)/DXPA)+0.5) + 1
           JG3=INT(((RYY1-RYFIX)/DYPA)+0.5) + 1
           KG3=INT(((RZZ1-RZFIX)/DZPA)+0.5) + 1

           IG4=IG3 + NN1 - 1
           JG4=JG3 + NN2 - 1
           KG4=KG3 + NN3 - 1

           IF (IG1.LE.IG4.AND.IG3.LE.IG2.AND.
     &         JG1.LE.JG4.AND.JG3.LE.JG2.AND.
     &         KG1.LE.KG4.AND.KG3.LE.KG2) THEN
            NV=NV+1
            VECINO(NV,I2)=J
           END IF

          END IF
         END DO
         NVECI(I2)=NV
       END DO

       IF (MAXVAL(NVECI(1:NPATCH(IR))).GT.NPALEV2) WRITE(*,*)
     &    'ERROR: gvecino second dimension too large',
     &     MAXVAL(NVECI(1:NPATCH(IR)))

*      Identify overlapping cells
       DO I=LOW1,LOW2

         L1=PATCHX(I)
         L2=PATCHY(I)
         L3=PATCHZ(I)

         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         RX1=PATCHRX(I)-0.5*DXPA
         RY1=PATCHRY(I)-0.5*DYPA
         RZ1=PATCHRZ(I)-0.5*DZPA
         RX2=PATCHRX(I)-0.5*DXPA+(N1-1)*DXPA
         RY2=PATCHRY(I)-0.5*DYPA+(N2-1)*DYPA
         RZ2=PATCHRZ(I)-0.5*DZPA+(N3-1)*DZPA

         I2=I-LOW1+1

         DO K=1,NVECI(I2)
          J=VECINO(K,I2)

          LL1=PATCHX(J)
          LL2=PATCHY(J)
          LL3=PATCHZ(J)

          NN1=PATCHNX(J)
          NN2=PATCHNY(J)
          NN3=PATCHNZ(J)

          RXX1=PATCHRX(J)-0.5*DXPA
          RYY1=PATCHRY(J)-0.5*DYPA
          RZZ1=PATCHRZ(J)-0.5*DZPA
          RXX2=PATCHRX(J)-0.5*DXPA+(NN1-1)*DXPA
          RYY2=PATCHRY(J)-0.5*DYPA+(NN2-1)*DYPA
          RZZ2=PATCHRZ(J)-0.5*DZPA+(NN3-1)*DZPA

*         X
          IF (RXX1.GE.RX1.AND.RXX2.LE.RX2) THEN
             CORNX1=INT(((RXX1-RX1)/DXPA)+0.5) + 1
             CORNX2=INT(((RXX2-RX1)/DXPA)+0.5) + 1
             CORNXX1=1
             CORNXX2=NN1
          END IF
          IF (RXX1.GE.RX1.AND.RXX2.GT.RX2) THEN
             CORNX1=INT(((RXX1-RX1)/DXPA)+0.5) + 1
             CORNX2=N1
             CORNXX1=1
             CORNXX2=INT(((RX2-RXX1)/DXPA)+0.5) +1
          END IF
          IF (RXX2.LE.RX2.AND.RXX1.LT.RX1) THEN
             CORNX1=1
             CORNX2=INT(((RXX2-RX1)/DXPA)+0.5) + 1
             CORNXX1=INT(((RX1-RXX1)/DXPA)+0.5) + 1
             CORNXX2=NN1
          END IF
          IF (RXX1.LT.RX1.AND.RXX2.GT.RX2) THEN
             CORNX1=1
             CORNX2=N1
             CORNXX1=INT(((RX1-RXX1)/DXPA)+0.5) + 1
             CORNXX2=INT(((RX2-RXX1)/DXPA)+0.5) + 1
          END IF

*         Y
          IF (RYY1.GE.RY1.AND.RYY2.LE.RY2) THEN
             CORNY1=INT(((RYY1-RY1)/DYPA)+0.5) + 1
             CORNY2=INT(((RYY2-RY1)/DYPA)+0.5) + 1
             CORNYY1=1
             CORNYY2=NN2
          END IF
          IF (RYY1.GE.RY1.AND.RYY2.GT.RY2) THEN
             CORNY1=INT(((RYY1-RY1)/DYPA)+0.5) + 1
             CORNY2=N2
             CORNYY1=1
             CORNYY2=INT(((RY2-RYY1)/DYPA)+0.5) +1
          END IF
          IF (RYY2.LE.RY2.AND.RYY1.LT.RY1) THEN
             CORNY1=1
             CORNY2=INT(((RYY2-RY1)/DYPA)+0.5) + 1
             CORNYY1=INT(((RY1-RYY1)/DYPA)+0.5) + 1
             CORNYY2=NN2
          END IF
          IF (RYY1.LT.RY1.AND.RYY2.GT.RY2) THEN
             CORNY1=1
             CORNY2=N2
             CORNYY1=INT(((RY1-RYY1)/DYPA)+0.5) + 1
             CORNYY2=INT(((RY2-RYY1)/DYPA)+0.5) + 1
          END IF

*         Z
          IF (RZZ1.GE.RZ1.AND.RZZ2.LE.RZ2) THEN
             CORNZ1=INT(((RZZ1-RZ1)/DZPA)+0.5) + 1
             CORNZ2=INT(((RZZ2-RZ1)/DZPA)+0.5) + 1
             CORNZZ1=1
             CORNZZ2=NN3
          END IF
          IF (RZZ1.GE.RZ1.AND.RZZ2.GT.RZ2) THEN
             CORNZ1=INT(((RZZ1-RZ1)/DZPA)+0.5) + 1
             CORNZ2=N3
             CORNZZ1=1
             CORNZZ2=INT(((RZ2-RZZ1)/DZPA)+0.5) +1
          END IF
          IF (RZZ2.LE.RZ2.AND.RZZ1.LT.RZ1) THEN
             CORNZ1=1
             CORNZ2=INT(((RZZ2-RZ1)/DZPA)+0.5) + 1
             CORNZZ1=INT(((RZ1-RZZ1)/DZPA)+0.5) + 1
             CORNZZ2=NN3
          END IF
          IF (RZZ1.LT.RZ1.AND.RZZ2.GT.RZ2) THEN
             CORNZ1=1
             CORNZ2=N3
             CORNZZ1=INT(((RZ1-RZZ1)/DZPA)+0.5) + 1
             CORNZZ2=INT(((RZ2-RZZ1)/DZPA)+0.5) + 1
          END IF

***        celdas madre del nivel inferior
           DO KK=CORNZZ1,CORNZZ2
           DO JJ=CORNYY1,CORNYY2
           DO II=CORNXX1,CORNXX2
            IX=II-CORNXX1+CORNX1
            JY=JJ-CORNYY1+CORNY1
            KZ=KK-CORNZZ1+CORNZ1
              ! MODIFIED DV 02-11-2021: AVOID BORDERS
            IF (SOLAP(IX,JY,KZ,I).EQ.1.AND.
     &          SOLAP(II,JJ,KK,J).EQ.1) THEN !there's an overlap we need to take care of
             IF (IX.EQ.1.OR.IX.EQ.N1.OR.
     &           JY.EQ.1.OR.JY.EQ.N2.OR.
     &           KZ.EQ.1.OR.KZ.EQ.N3) THEN ! if the master candidate is close to the patch border
              SOLAP(IX,JY,KZ,I)=0 ! then make master the other one
             ELSE
              SOLAP(II,JJ,KK,J)=0 ! default veinsgrid behaviour
             END IF
            END IF
           END DO
           END DO
           END DO

       END DO
       END DO

      RETURN
      END

*************************************************************
       SUBROUTINE MALLA(NX,NY,NZ,LADO)
*************************************************************

       IMPLICIT NONE
       INTEGER NX,I,NY,J,NZ,K

       INCLUDE 'input_files/asohf_parameters.dat'

       REAL*4 A,B,C,LADO

       REAL*4  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
       COMMON /GRID/   RADX,RADY,RADZ

       REAL*4 DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ
*      GENERAL INITIAL CONDITIONS

*      GRID LIMITS
       A=-LADO/2.0
       B=LADO/2.0


*      GRID

*      X-AXIS
       C=(B-A)/(NX-1)
       RADX(1)=A
       DO I=2,NX
        RADX(I)=RADX(1)+(I-1)*C
       END DO

*      FICTICIUS CELLS
       RADX(0)=RADX(1)-C
       RADX(NX+1)=RADX(NX)+C

*      Y-AXIS
       C=(B-A)/(NY-1)
       RADY(1)=A
       DO J=2,NY
        RADY(J)=RADY(1)+(J-1)*C
       END DO

*     FICTICIUS CELLS
       RADY(0)=RADY(1)-C
       RADY(NY+1)=RADY(NY)+C

*      Z-AXIS
       C=(B-A)/(NZ-1)
       RADZ(1)=A
       DO K=2,NZ
        RADZ(K)=RADZ(1)+(K-1)*C
       END DO

*      FICTICIUS CELLS
       RADZ(0)=RADZ(1)-C
       RADZ(NZ+1)=RADZ(NZ)+C


       DX=RADX(2)-RADX(1)
       DY=RADY(2)-RADY(1)
       DZ=RADZ(2)-RADZ(1)


       RETURN
       END

*************************************************************
       SUBROUTINE CLEAN_OVERLAPS(NL,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                           SOLAP,FIELD)
*************************************************************
*      Sets to 0 the values in FIELD which are overlapped
*      by other AMR cells (and not the master)
*************************************************************

       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NL
       INTEGER NPATCH(0:NL)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER SOLAP(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL FIELD(NAMRX,NAMRY,NAMRZ,NPALEV)

       INTEGER IR,LOW1,LOW2,N1,N2,N3,IX,JY,KZ,I

       DO IR=1,NL
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))

!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,FIELD,SOLAP,LOW1,
!$OMP+                   LOW2),
!$OMP+            PRIVATE(N1,N2,N3,I,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
        DO I=LOW1,LOW2
         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         DO KZ=1,N3
         DO JY=1,N2
         DO IX=1,N1
          FIELD(IX,JY,KZ,I)=FIELD(IX,JY,KZ,I)*SOLAP(IX,JY,KZ,I)
         END DO
         END DO
         END DO

        END DO
       END DO

       RETURN
       END

*************************************************************
       SUBROUTINE CLEAN_OVERLAPS_INT(NL,NPATCH,PATCHNX,PATCHNY,
     &                               PATCHNZ,SOLAP,FIELD)
*************************************************************
*      Sets to 0 the values in FIELD which are overlapped
*      by other AMR cells (and not the master)
*************************************************************

       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NL
       INTEGER NPATCH(0:NL)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER SOLAP(NAMRX,NAMRY,NAMRZ,NPALEV)
       INTEGER FIELD(NAMRX,NAMRY,NAMRZ,NPALEV)

       INTEGER IR,LOW1,LOW2,N1,N2,N3,IX,JY,KZ,I

       DO IR=1,NL
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))

!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,FIELD,SOLAP,LOW1,
!$OMP+                   LOW2),
!$OMP+            PRIVATE(N1,N2,N3,I,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
        DO I=LOW1,LOW2
         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         DO KZ=1,N3
         DO JY=1,N2
         DO IX=1,N1
          FIELD(IX,JY,KZ,I)=FIELD(IX,JY,KZ,I)*SOLAP(IX,JY,KZ,I)
         END DO
         END DO
         END DO

        END DO
       END DO

       RETURN
       END

*************************************************************
       SUBROUTINE RENORM_DENSITY(NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                           PATCHNZ,CR0AMR,CR0AMR11,SOLAP,U1,U11,
     &                           LADO0,RODO,RE0)
*************************************************************
*      Normalizes density so that the mean density contrast over the
*      box is 1 (to account for unaccounted species, e.g. gas)
*************************************************************
        IMPLICIT NONE
        INCLUDE 'input_files/asohf_parameters.dat'

        INTEGER NL,NX,NY,NZ
        INTEGER NPATCH(0:NL)
        INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
        INTEGER CR0AMR(NMAX,NMAY,NMAZ)
        INTEGER CR0AMR11(NAMRX,NAMRY,NAMRZ,NPALEV)
        INTEGER SOLAP(NAMRX,NAMRY,NAMRZ,NPALEV)
        REAL*4 U1(NMAX,NMAY,NMAZ)
        REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
        REAL LADO0,RODO,RE0

        REAL*4 DX,DY,DZ
        COMMON /ESPACIADO/ DX,DY,DZ

        REAL BASVOL,BASMASS,BASMASS2,DXPA,DYPA,DZPA
        INTEGER IX,JY,KZ,I,IR,LOW1,LOW2,N1,N2,N3

        BASMASS=0.0
        BASVOL=DX*DY*DZ*RE0**3
!$OMP PARALLEL DO SHARED(NX,NY,NZ,CR0AMR,BASVOL,U1),
!$OMP+            PRIVATE(IX,JY,KZ),
!$OMP+            REDUCTION(+:BASMASS),
!$OMP+            DEFAULT(NONE)
        DO KZ=1,NZ
        DO JY=1,NY
        DO IX=1,NX
         IF (CR0AMR(IX,JY,KZ).EQ.1) BASMASS=BASMASS+U1(IX,JY,KZ)*BASVOL
        END DO
        END DO
        END DO

        DO IR=1,NL
         DXPA=DX/2.0**IR
         DYPA=DY/2.0**IR
         DZPA=DZ/2.0**IR
         BASVOL=DXPA*DYPA*DZPA*RE0**3

         LOW1=SUM(NPATCH(0:IR-1))+1
         LOW2=SUM(NPATCH(0:IR))

         BASMASS2=0.0
!$OMP PARALLEL DO SHARED(CR0AMR11,SOLAP,BASVOL,U11,PATCHNX,PATCHNY,
!$OMP+                   PATCHNZ,LOW1,LOW2),
!$OMP+            PRIVATE(IX,JY,KZ,I,N1,N2,N3),
!$OMP+            REDUCTION(+:BASMASS2),
!$OMP+            DEFAULT(NONE)
         DO I=LOW1,LOW2
          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)

          DO KZ=1,N3
          DO JY=1,N2
          DO IX=1,N1
           IF (CR0AMR11(IX,JY,KZ,I).EQ.1.AND.
     &         SOLAP(IX,JY,KZ,I).EQ.1) THEN
            BASMASS2=BASMASS2+U11(IX,JY,KZ,I)*BASVOL
           END IF
          END DO
          END DO
          END DO
         END DO

         BASMASS=BASMASS+BASMASS2
        END DO

        ! TOTAL MASS IN CODE UNITS
        BASMASS=BASMASS*RODO
        ! OVERDENSITY IN THE BOX
        BASMASS=BASMASS/(RODO*RE0**3*LADO0**3)

!$OMP PARALLEL DO SHARED(NX,NY,NZ,BASMASS,U1),
!$OMP+            PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
        DO KZ=1,NZ
        DO JY=1,NY
        DO IX=1,NX
         U1(IX,JY,KZ)=U1(IX,JY,KZ)/BASMASS
        END DO
        END DO
        END DO

        DO IR=1,NL
         LOW1=SUM(NPATCH(0:IR-1))+1
         LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,U11,BASMASS),
!$OMP+            PRIVATE(I,N1,N2,N3,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
         DO I=LOW1,LOW2
          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)
          DO KZ=1,N3
          DO JY=1,N2
          DO IX=1,N1
           U11(IX,JY,KZ,I)=U11(IX,JY,KZ,I)/BASMASS
          END DO
          END DO
          END DO
         END DO
        END DO



        RETURN
        END
