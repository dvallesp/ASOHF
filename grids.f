************************************************************************
      SUBROUTINE CREATE_MESH(ITER,NX,NY,NZ,NL_MESH,NPATCH,PARE,
     &           PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &           PATCHRX,PATCHRY,PATCHRZ,RXPA,RYPA,RZPA,U2DM,U3DM,
     &           U4DM,MASAP,N_PARTICLES,N_DM,N_GAS,LADO0,T,ZETA)
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

*     COMMON VARIABLES
      REAL DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      REAL  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
      COMMON /GRID/  RADX,RADY,RADZ

*     LOCAL VARIABLES
      INTEGER PLEV(PARTIRED)
      INTEGER,ALLOCATABLE::CR0(:,:,:)
      INTEGER,ALLOCATABLE::CR01(:,:,:,:)
      INTEGER,ALLOCATABLE::CONTA1(:,:,:)
      INTEGER,ALLOCATABLE::CONTA11(:,:,:,:)
      REAL MAP,XL,YL,ZL,DXPA,DYPA,DZPA
      INTEGER I,IX,JY,KZ,REFINE_THR,REFINE_COUNT,BOR,MIN_PATCHSIZE
      INTEGER INI_EXTENSION,NBIS,IRPA,BORAMR,LOW1,LOW2,IPATCH,IPARE
      INTEGER INMAX(3),INMAX2(2),I1,I2,J1,J2,K1,K2,N1,N2,N3,IR,MARCA
      INTEGER NP1,NP2,NP3,BASINT,NPALEV3,II,JJ,KK

      INTEGER,ALLOCATABLE::LNPATCH(:)
      INTEGER,ALLOCATABLE::LPATCHNX(:,:),LPATCHNY(:,:),LPATCHNZ(:,:)
      INTEGER,ALLOCATABLE::LPATCHX(:,:),LPATCHY(:,:),LPATCHZ(:,:)
      REAL,ALLOCATABLE::LPATCHRX(:,:),LPATCHRY(:,:),LPATCHRZ(:,:)
      INTEGER,ALLOCATABLE::LVAL(:,:)

      CHARACTER*15 FILE6
      CHARACTER*30 FILERR

!     hard-coded parameters (for now, at least)
      REFINE_THR=3
      BOR=8
      BORAMR=0
      INI_EXTENSION=2 !initial extension of a patch around a cell (on each direction)
      MIN_PATCHSIZE=14 !minimum size (child cells) to be accepted
      NPALEV3=(INT(NAMRX/5)**3)+1
      write(*,*) 'NPALEV3=',NPALEV3

      MAP=MAXVAL(MASAP(1:N_PARTICLES))

      PLEV=0
!$OMP PARALLEL DO SHARED(N_PARTICLES,PLEV,MAP,MASAP),PRIVATE(I),
!$OMP+            DEFAULT(NONE)
      DO I=1,N_PARTICLES
       PLEV(I)=LOG(MAP/MASAP(I)+.5)/LOG(8.)
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
      DO IX=1,NX
      DO JY=1,NY
      DO KZ=1,NZ
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
      DO IX=1,NX
      DO JY=1,NY
      DO KZ=1,NZ
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

      IPATCH=0

      DO WHILE (REFINE_COUNT.GT.0.AND.IPATCH.LT.NPALEV) !--------------
       INMAX=MAXLOC(CR0)
       IX=INMAX(1)
       JY=INMAX(2)
       KZ=INMAX(3)
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
         IF (COUNT(CONTA1(I1-1,J1:J2,K1:K2).GE.REFINE_THR).GT.0) THEN
          I1=I1-1
          N1=2*(I2-I1+1)
          MARCA=1
         END IF
        END IF

        IF (N1.LE.NAMRX-2.AND.I2.LT.NX-BOR) THEN
         !IF (IPATCH.EQ.153) WRITE(*,*) IX,JY,KZ,I1,I2,J1,J2,K1,K2
         IF (COUNT(CONTA1(I2+1,J1:J2,K1:K2).GE.REFINE_THR).GT.0) THEN
          I2=I2+1
          N1=2*(I2-I1+1)
          MARCA=1
         END IF
        END IF

        IF (N2.LE.NAMRY-2.AND.J1.GT.BOR+1) THEN
         IF (COUNT(CONTA1(I1:I2,J1-1,K1:K2).GE.REFINE_THR).GT.0) THEN
          J1=J1-1
          N2=2*(J2-J1+1)
          MARCA=1
         END IF
        END IF

        IF (N2.LE.NAMRY-2.AND.J2.LT.NY-BOR) THEN
         IF (COUNT(CONTA1(I1:I2,J2+1,K1:K2).GE.REFINE_THR).GT.0) THEN
          J2=J2+1
          N2=2*(J2-J1+1)
          MARCA=1
         END IF
        END IF

        IF (N3.LE.NAMRZ-2.AND.K1.GT.BOR+1) THEN
         IF (COUNT(CONTA1(I1:I2,J1:J2,K1-1).GE.REFINE_THR).GT.0) THEN
          K1=K1-1
          N3=2*(K2-K1+1)
          MARCA=1
         END IF
        END IF

        IF (N3.LE.NAMRZ-2.AND.K2.LT.NZ-BOR) THEN
         IF (COUNT(CONTA1(I1:I2,J1:J2,K2+1).GE.REFINE_THR).GT.0) THEN
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
*       WRITE(*,*) IPATCH,N1,N2,N3,
*     &             COUNT(CONTA1(I1:I2,J1:J2,K1:K2).GE.REFINE_THR)
        CONTA1(I1:I2,J1:J2,K1:K2)=0
        CR0(I1:I2,J1:J2,K1:K2)=-1

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
        DO IX=1,NAMRX
        DO JY=1,NAMRY
        DO KZ=1,NAMRZ
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

        DO IX=1,N1
        DO JY=1,N2
        DO KZ=1,N3
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
!$OMP+                   PATCHRY,PATCHRZ),
!$OMP+            PRIVATE(IPARE,REFINE_COUNT,IPATCH,INMAX,IX,JY,KZ,
!$OMP+                    BASINT,NP1,NP2,NP3,I1,I2,J1,J2,K1,K2,N1,N2,N3,
!$OMP+                    MARCA,NBIS),
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
           IF (COUNT(CONTA11(I1-1,J1:J2,K1:K2,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            I1=I1-1
            N1=2*(I2-I1+1)
            MARCA=1
           END IF
          END IF

          IF (N1.LE.NAMRX-2.AND.I2.LT.NP1-BORAMR) THEN
           IF (COUNT(CONTA11(I2+1,J1:J2,K1:K2,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            I2=I2+1
            N1=2*(I2-I1+1)
            MARCA=1
           END IF
          END IF

          IF (N2.LE.NAMRY-2.AND.J1.GT.BORAMR+1) THEN
           IF (COUNT(CONTA11(I1:I2,J1-1,K1:K2,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            J1=J1-1
            N2=2*(J2-J1+1)
            MARCA=1
           END IF
          END IF

          IF (N2.LE.NAMRY-2.AND.J2.LT.NP2-BORAMR) THEN
           IF (COUNT(CONTA11(I1:I2,J2+1,K1:K2,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            J2=J2+1
            N2=2*(J2-J1+1)
            MARCA=1
           END IF
          END IF

          IF (N3.LE.NAMRZ-2.AND.K1.GT.BORAMR+1) THEN
           IF (COUNT(CONTA11(I1:I2,J1:J2,K1-1,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            K1=K1-1
            N3=2*(K2-K1+1)
            MARCA=1
           END IF
          END IF

          IF (N3.LE.NAMRZ-2.AND.K2.LT.NP3-BORAMR) THEN
           IF (COUNT(CONTA11(I1:I2,J1:J2,K2+1,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
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

          CONTA11(I1:I2,J1:J2,K1:K2,IPARE)=0
          CR01(I1:I2,J1:J2,K1:K2,IPARE)=-1

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

*     WRITING GRID DATA ON A FILE
      CALL NOMFILE6(ITER,FILE6)
      FILERR='./output_files/'//FILE6
      OPEN(33,FILE=FILERR,STATUS='UNKNOWN')

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

      RETURN
      END

************************************************************************
      SUBROUTINE INTERPOLATE_DENSITY(ITER,NX,NY,NZ,NL_MESH,NPATCH,PARE,
     &           PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &           PATCHRX,PATCHRY,PATCHRZ,RXPA,RYPA,RZPA,MASAP,
     &           N_PARTICLES,N_DM,N_GAS,LADO0,T,ZETA)
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
     &       MASAP(PARTIRED)
      REAL LADO0,T,ZETA

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
      INTEGER PLEV(PARTIRED) ! This variable has now a different content:
      ! Now it will store the maximum mesh level of a particle.
      REAL XL,YL,ZL,DXPA,DYPA,DZPA,XR,YR,ZR,XP,YP,ZP,BAS,DENBAS
      INTEGER I,IX,JY,KZ,II,JJ,KK,N1,N2,N3,I1,I2,J1,J2,K1,K2,IR,IPATCH
      INTEGER BUF,LOW1,LOW2,NBAS
      REAL,ALLOCATABLE::RXPA2(:),RYPA2(:),RZPA2(:),MASAP2(:)
      REAL MINIRADX(-2:NAMRX+3),MINIRADY(-2:NAMRY+3),
     &     MINIRADZ(-2:NAMRZ+3)
      REAL VX(-1:1),VY(-1:1),VZ(-1:1)

      BUF=1 ! (1 cell extra buffer (TSC))

!$OMP PARALLEL DO SHARED(N_PARTICLES,PLEV), PRIVATE(I), DEFAULT(NONE)
      DO I=1,N_PARTICLES
       PLEV(I)=0
      END DO

      DO IR=9,1,-1
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
          XR=XL+(PATCHNX(IPATCH)+2)*DXPA
          YR=YL+(PATCHNY(IPATCH)+2)*DYPA
          ZR=ZL+(PATCHNZ(IPATCH)+2)*DZPA
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

      DO IR=0,9
       WRITE(*,*) 'Particles at level',IR,
     &             COUNT(PLEV(1:N_PARTICLES).EQ.IR)
      END DO

      DO IR=9,1,-1
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
        DO IX=1,N1
        DO JY=1,N2
        DO KZ=1,N3
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
        DO IX=1,N1
        DO JY=1,N2
        DO KZ=1,N3
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
      END DO !ir=9,1,-1

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
!$OMP+                   RADX,RADY,RADZ,NX,NY,NZ,U1,MASAP),
!$OMP+            PRIVATE(I,XP,YP,ZP,IX,JY,KZ,BAS,VX,VY,VZ,II,JJ,KK,
!$OMP+                    I1,J1,K1),
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
        if (MASAP(I)*VX(II)*VY(JJ)*VZ(KK).lt.0) then
         write(*,*) '-------'
         write(*,*) xp,yp,zp
         write(*,*) ix,jy,kz
         write(*,*) ii,jj,kk
         write(*,*) vx(ii),vy(jj),vz(kk)
        end if
        U1(I1,J1,K1)= U1(I1,J1,K1) + MASAP(I)*VX(II)*VY(JJ)*VZ(KK)
       END DO
       END DO
       END DO

      END DO

      DENBAS=DX*DY*DZ*ROTE*RETE**3
!$OMP PARALLEL DO SHARED(NX,NY,NZ,DENBAS,U1), PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
      DO IX=1,NX
      DO JY=1,NY
      DO KZ=1,NZ
       U1(IX,JY,KZ)=U1(IX,JY,KZ)/DENBAS
      END DO
      END DO
      END DO

      WRITE(*,*) 'At level',0,minval(u1(:,:,:)),maxval(u1(:,:,:))

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

************************************************************************
       SUBROUTINE MESHRENOEF_OLD(ITER,NX,NY,NZ,NL,COTA,NPATCH,
     &                   NPART,PATCHNX,PATCHNY,PATCHNZ,
     &                   PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,
     &                   PATCHRZ,PARE,U2DM,U3DM,U4DM,MASAP,MAP,
     &                   RXPA,RYPA,RZPA,ZETA,T,LADO0,FLAG_MASCLET,PLOT)
************************************************************************
*      Builds the grid from a set of particles
************************************************************************
*      OJO EN ESTA SRUTINA CON: 1) BOR E INCRE

       IMPLICIT NONE

       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NPALEV3       !max num de parches q pueden salir de un parche
*       PARAMETER (NPALEV3=(N_ESP*(INT(NAMRX/3)+1)**3))
**       PARAMETER (NPALEV3=(N_ESP*(INT(NAMRX/9)**3)+1))

       REAL*4 RETE,HTE,ROTE
       COMMON /BACK/ RETE,HTE,ROTE

       INTEGER NX,NY,NZ,NL1,NL10,NL,ITER

       REAL*4  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
       COMMON /GRID/   RADX,RADY,RADZ

       REAL*4  UP1(NMAX,NMAY,NMAZ,N_ESP)
       REAL*4  UP11(NAMRX,NAMRY,NAMRZ,NPALEV,N_ESP)

       REAL*4 DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       REAL*4 COTA(NCOTAS,0:NLEVELS)
       REAL*4 ZETA,T,MAP,MAXMAP
       REAL*4 MASA_TEMP(NAMRX,NAMRY,NAMRZ)
       REAL*4 MASA_TEMP0(NMAX,NMAY,NMAZ)
       INTEGER N_GAS,N_DM,FLAG_MASCLET

       CHARACTER*9 FILNOM1,FILNOM2
       CHARACTER*10 FILNOM3
       CHARACTER*25 FIL1,FIL2
       CHARACTER*26 FIL3
       CHARACTER*15 FILE6
       CHARACTER*30 FILERR

       INTEGER NUM
       COMMON /PROCESADORES/ NUM

       INTEGER,ALLOCATABLE:: LNPATCH(:)        ! variables auxiliares paralelizacion
       INTEGER,ALLOCATABLE:: LPARE(:,:)
       INTEGER,ALLOCATABLE:: LPATCHNX(:,:)
       INTEGER,ALLOCATABLE:: LPATCHNY(:,:)
       INTEGER,ALLOCATABLE:: LPATCHNZ(:,:)
       INTEGER,ALLOCATABLE:: LPATCHX(:,:)
       INTEGER,ALLOCATABLE:: LPATCHY(:,:)
       INTEGER,ALLOCATABLE:: LPATCHZ(:,:)
       REAL*4,ALLOCATABLE :: LPATCHRX(:,:)
       REAL*4,ALLOCATABLE :: LPATCHRY(:,:)
       REAL*4,ALLOCATABLE :: LPATCHRZ(:,:)

       INTEGER CR0(-1:NMAX+2,-1:NMAY+2,-1:NMAZ+2,N_ESP)   !ANAYDIMOS 2 CELDAS FICTICIAS
       REAL*4 UBAS(NMAX,NMAY,NMAZ),UBAS2(NAMRX,NAMRY,NAMRZ)
       INTEGER CR02(-1:NAMRX+2,-1:NAMRY+2,-1:NAMRZ+2,N_ESP)

       INTEGER I,J,K,IX,JY,KZ,I1,J1,K1,IR,IPALE,BOR,IRR
       INTEGER NEF,NCELL,PAX1,PAX2,PAY1,PAY2,PAZ1,PAZ2
       INTEGER IPATCH,II,JJ,KK,N1,N2,N3,NEF1,NEF2
       INTEGER INMAX(3),KONTA,NPARTICULAS,MARCA
       INTEGER NBAS,IPA2,NPX,NPY,NPZ,INCRE,NBIS
       REAL*4 DXPA,DYPA,DZPA,LADO0,LADO,U1MIN
       INTEGER IP,NP1,NP2,NP3,L1,L2,L3, PLOT

       INTEGER IESP,NDM_ESP(N_ESP),NDXYZ
       INTEGER KK1, KK2,IPATCH_ESP, CONTA
       INTEGER NESP_BAS


*////OUTPUTS DE ESTA RUTINA//////////////
*      VARIABLES
       REAL*4 U1(NMAX,NMAY,NMAZ)
       REAL*4 U1G(NMAX,NMAY,NMAZ)
       REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL*4 U11G(NAMRX,NAMRY,NAMRZ,NPALEV)
       COMMON /VARIA/ U1,U11,U1G,U11G

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

*      DM PARTICLES
       INTEGER NPART(0:NLEVELS)
       REAL*4 U2DM(PARTIRED)
       REAL*4 U3DM(PARTIRED)
       REAL*4 U4DM(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       REAL*4 RXPA(PARTIRED)
       REAL*4 RYPA(PARTIRED)
       REAL*4 RZPA(PARTIRED)
       INTEGER ORIPA1(PARTIRED),ORIPA2(PARTIRED)
       COMMON /PUNTEROS/ ORIPA1, ORIPA2

*      GAS PARTICLES
       REAL*4, ALLOCATABLE::RXG(:)
       REAL*4, ALLOCATABLE::RYG(:)
       REAL*4, ALLOCATABLE::RZG(:)
       REAL*4, ALLOCATABLE::VXG(:)
       REAL*4, ALLOCATABLE::VYG(:)
       REAL*4, ALLOCATABLE::VZG(:)
       REAL*4, ALLOCATABLE::MASAG(:)
*//////////////////////////////////////////

       LOGICAL PARALEL

       INTEGER NL_2,PABAS
       INTEGER LOW1, LOW2
       INTEGER IPATCH2,IPALE2

       INTEGER PARCHLIM       ! =0, sin limite por nivel, =/ 0 limites
       INTEGER MPAPOLEV(100)  !maximo numero de parches por nivel, 100 niveles max.
       COMMON /LIMPARCH/ PARCHLIM,MPAPOLEV

       NPALEV3=(INT(NAMRX/3)+1)**3
       NPALEV3=NPALEV3*N_ESP

*////////////////////////////////////////////
*      Reading external file with particles
*////////////////////////////////////////////

       LADO=LADO0
       CR0=0
       CR02=0

       WRITE(*,*)'NPALEV3=', NPALEV3

       NPART=0
       NPARTICULAS=0
       N_GAS=0
       N_DM=0
       KONTA=0

****
       PABAS=PARTIRED
!$OMP PARALLEL DO SHARED(PABAS,U2DM,U3DM,U4DM,RXPA,RYPA,RZPA,
!$OMP+                   MASAP,ORIPA1,ORIPA2),
!$OMP+            PRIVATE(I)
       DO I=1,PABAS
        U2DM(I)=0.0
        U3DM(I)=0.0
        U4DM(I)=0.0
        RXPA(I)=0.0
        RYPA(I)=0.0
        RZPA(I)=0.0
        MASAP(I)=0.0
        ORIPA1(I)=0
        ORIPA2(I)=0
       END DO
****

*------------------------------------------*
       IF (FLAG_MASCLET.EQ.0) THEN   !!From a external file of particles
*------------------------------------------*

       OPEN (5,FILE='particle_list.dat',
     &               STATUS='UNKNOWN',ACTION='READ')

       READ(5,*) NPARTICULAS,ZETA,T,N_GAS,N_DM
       ORIPA1(1:N_DM)=0   !todas en nivel 0

       IF (NPARTICULAS.NE.N_GAS+N_DM) THEN
        WRITE(*,*) 'WARNING: NPARTICULAS.NE.N_GAS+N_DM!!',
     &              NPARTICULAS, N_GAS, N_DM
        STOP
       END IF
       WRITE(*,*)'N_DM,N_GAS=',N_DM,N_GAS

       IF (N_GAS.GT.0) THEN

        ALLOCATE(RXG(N_GAS))
        ALLOCATE(RYG(N_GAS))
        ALLOCATE(RZG(N_GAS))
        ALLOCATE(VXG(N_GAS))
        ALLOCATE(VYG(N_GAS))
        ALLOCATE(VZG(N_GAS))
        ALLOCATE(MASAG(N_GAS))

        RXG=0.0
        RYG=0.0
        RZG=0.0
        VXG=0.0
        VYG=0.0
        VZG=0.0
        MASAG=0.0

       END IF


**OJO!! Tal y como leemos a continuacion, estamos pasando del gas al
*        construir la malla aunque si que lo leemos

       DO I=1,NPARTICULAS

       IF(I.LE.N_GAS.AND.N_GAS.GT.0) THEN

       READ(5,*) JJ, RXG(I),RYG(I),RZG(I),
     &           VXG(I),VYG(I),VZG(I),
     &           MASAG(I)

       ELSE

       KONTA=KONTA+1

       READ(5,*) JJ, RXPA(KONTA),RYPA(KONTA),RZPA(KONTA),
     &           U2DM(KONTA),U3DM(KONTA),U4DM(KONTA),
     &           MASAP(KONTA)

       ORIPA2(KONTA)=JJ
       END IF

       END DO

       WRITE(*,*) 'ORIPA1=',MAXVAL(ORIPA1(1:KONTA)),
     &                      MINVAL(ORIPA1(1:KONTA))
       WRITE(*,*) 'ORIPA2=',MAXVAL(ORIPA2(1:KONTA)),
     &                      MINVAL(ORIPA2(1:KONTA))

       CLOSE(5)

*DM
       WRITE(*,*) '///DM///'
       WRITE(*,*) MAXVAL(RXPA(1:N_DM)),
     &            MINVAL(RXPA(1:N_DM))
       WRITE(*,*) MAXVAL(RYPA(1:N_DM)),
     &            MINVAL(RYPA(1:N_DM))
       WRITE(*,*) MAXVAL(RZPA(1:N_DM)),
     &            MINVAL(RZPA(1:N_DM))
       WRITE(*,*) MAXVAL(MASAP(1:N_DM)),
     &            MINVAL(MASAP(1:N_DM))

*!OJO CON ESTO Q PARA EL GAS SERIA DIFERENTE NO?
       MAXMAP=0.0
       MAXMAP=MAXVAL(MASAP(1:N_DM))
       MASAP=MASAP/MAXMAP

*      CORREGIMOS PARA TENER LAS COORDENADAS EN [-L/2,L/2]
       RXPA(1:N_DM)=RXPA(1:N_DM)-LADO0*0.5
       RYPA(1:N_DM)=RYPA(1:N_DM)-LADO0*0.5
       RZPA(1:N_DM)=RZPA(1:N_DM)-LADO0*0.5

       WRITE(*,*)'CHECKING BOX SIDE', LADO
       IX=0
       IX=COUNT(ABS(RXPA(1:N_DM)).GT.LADO*0.5001)
       JY=0
       JY=COUNT(ABS(RYPA(1:N_DM)).GT.LADO*0.5001)
       KZ=0
       KZ=COUNT(ABS(RZPA(1:N_DM)).GT.LADO*0.5001)
       IF(IX.GT.0.OR.JY.GT.0.OR.KZ.GT.0) THEN
        WRITE(*,*) IX,JY,KZ
        WRITE(*,*)'WARNING DM: LADO!!'
        STOP
       ENDIF

       WRITE(*,*) '///DM///'
       WRITE(*,*) MAXVAL(RXPA(1:N_DM)),
     &            MINVAL(RXPA(1:N_DM))
       WRITE(*,*) MAXVAL(RYPA(1:N_DM)),
     &            MINVAL(RYPA(1:N_DM))
       WRITE(*,*) MAXVAL(RZPA(1:N_DM)),
     &            MINVAL(RZPA(1:N_DM))


       IF(N_GAS.GT.0) THEN

       WRITE(*,*) '///GAS///'
       WRITE(*,*) MAXVAL(RXG(1:N_GAS)),
     &            MINVAL(RXG(1:N_GAS))
       WRITE(*,*) MAXVAL(RYG(1:N_GAS)),
     &            MINVAL(RYG(1:N_GAS))
       WRITE(*,*) MAXVAL(RZG(1:N_GAS)),
     &            MINVAL(RZG(1:N_GAS))

       RXG(1:N_GAS)=RXG(1:N_GAS)-LADO0*0.5
       RYG(1:N_GAS)=RYG(1:N_GAS)-LADO0*0.5
       RZG(1:N_GAS)=RZG(1:N_GAS)-LADO0*0.5

       IX=0
       IX=COUNT(ABS(RXG(1:N_GAS)).GT.LADO*0.5001)
       JY=0
       JY=COUNT(ABS(RYG(1:N_GAS)).GT.LADO*0.5001)
       KZ=0
       KZ=COUNT(ABS(RZG(1:N_GAS)).GT.LADO*0.5001)
       IF(IX.GT.0.OR.JY.GT.0.OR.KZ.GT.0) THEN
        WRITE(*,*) IX,JY,KZ
        WRITE(*,*)'WARNING GAS: LADO!!'
        STOP
       ENDIF

       WRITE(*,*) '///GAS///'
       WRITE(*,*) MAXVAL(RXG(1:N_GAS)),
     &            MINVAL(RXG(1:N_GAS))
       WRITE(*,*) MAXVAL(RYG(1:N_GAS)),
     &            MINVAL(RYG(1:N_GAS))
       WRITE(*,*) MAXVAL(RZG(1:N_GAS)),
     &            MINVAL(RZG(1:N_GAS))

       END IF !N_GAS



      END IF !FLAG_MASCLET


*------------------------------------------*
       IF (FLAG_MASCLET.EQ.1) THEN    !!!STAND-ALONE
*------------------------------------------*

       NPATCH=0
       NPART=0
       NL_2=0
       MAP=0.0

       CALL NOMFILE(ITER,FILNOM1,FILNOM2,FILNOM3)
       WRITE(*,*) 'leyendo',ITER,' ',FILNOM2,FILNOM3

       FIL2='simu_masclet/'//FILNOM2
       FIL3='simu_masclet/'//FILNOM3

       OPEN (33,FILE=FIL3,STATUS='UNKNOWN',ACTION='READ')
       OPEN (32,FILE=FIL2,
     &       STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')


**      GRID DATA
       READ(33,*) IRR,T,NL_2,MAP
       READ(33,*) ZETA
       READ(33,*) IR,NDXYZ
       WRITE(*,*) 'IR,NL_2,NDXYZ,MAP', IR,NL_2,NDXYZ,MAP
       DO IR=1,NL_2
       READ(33,*) IRR,NPATCH(IR), NPART(IR)
       READ(33,*)
       WRITE(*,*) 'IR,NPATCH(IR),NPART(IR)',IR,NPATCH(IR),NPART(IR)
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       IF (IR.NE.IRR) WRITE(*,*)'Warning: fail in restart'
       DO I=LOW1,LOW2
        READ(33,*) !PATCHNX(I),PATCHNY(I),PATCHNZ(I)
        READ(33,*) !PATCHX(I),PATCHY(I),PATCHZ(I)
        READ(33,*) !AAA,BBB,CCC
        READ(33,*) !PARE(I)
       END DO


       END DO
       CLOSE(33)

       NPART(0)=NDXYZ


**     DARK MATTER
        READ(32)
        READ(32) !(((U1(I,J,K),K=1,NZ),J=1,NY),I=1,NX)
        READ(32) (RXPA(I),I=1,NDXYZ)
        READ(32) (RYPA(I),I=1,NDXYZ)
        READ(32) (RZPA(I),I=1,NDXYZ)
        READ(32) (U2DM(I),I=1,NDXYZ)
        READ(32) (U3DM(I),I=1,NDXYZ)
        READ(32) (U4DM(I),I=1,NDXYZ)
C        READ(32) (ORIPA1(I),I=1,NDXYZ)      !OJO! las nuevas versioens de MASCLET no lo tienen
        READ(32) (ORIPA2(I),I=1,NDXYZ)       !partcile ID
        KONTA=NDXYZ
        MASAP(1:NDXYZ)=MAP

cx     !OJO! redefino ORIPA1 para poder usar el merger tree
       IF (PLOT.GT.1)  ORIPA1(1:NDXYZ)=0    !!!IR
cx

        CONTA=KONTA
       DO IR=1,NL_2
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2
        READ(32) !(((U11(IX,J,K,I),K=1,N3),J=1,N2),IX=1,N1)
       END DO

        READ(32) RXPA(CONTA+1:CONTA+NPART(IR))
        READ(32) RYPA(CONTA+1:CONTA+NPART(IR))
        READ(32) RZPA(CONTA+1:CONTA+NPART(IR))
        READ(32) U2DM(CONTA+1:CONTA+NPART(IR))
        READ(32) U3DM(CONTA+1:CONTA+NPART(IR))
        READ(32) U4DM(CONTA+1:CONTA+NPART(IR))
        READ(32) MASAP(CONTA+1:CONTA+NPART(IR))
C        READ(32) ORIPA1(CONTA+1:CONTA+NPART(IR))       !OJO! las nuevas versioens de MASCLET no lo tienen
        READ(32) ORIPA2(CONTA+1:CONTA+NPART(IR))

cx     !OJO! redefino ORIPA1 para poder usar el merger tree
       IF(PLOT.GT.1) THEN
        DO I=CONTA+1, CONTA+NPART(IR)
         IF (ORIPA2(I).LT.0) THEN
           ORIPA1(I)=0    !!!IR
           ORIPA2(I)=ABS(ORIPA2(I))
         ELSE
          ORIPA1(I)=IR
         END IF
        END DO
       END IF  !PLOT
cx

       CONTA=CONTA+NPART(IR)

       END DO   !IR
       CLOSE(32)

       KONTA=CONTA
       N_DM=SUM(NPART(0:NL_2))
       WRITE(*,*) 'PARTICULAS TOTALES EN ITER=',N_DM

       WRITE(*,*) 'ORIPA1=',MAXVAL(ORIPA1(1:KONTA)),
     &                      MINVAL(ORIPA1(1:KONTA))
       WRITE(*,*) 'ORIPA2=',MAXVAL(ORIPA2(1:KONTA)),
     &                      MINVAL(ORIPA2(1:KONTA))


**************************************
*      SEPARO LAS PARTI EN ESPECIES
**************************************
       NDM_ESP(1)=SUM(NPART(0:NL_2))
CX       NDM_ESP(1)=NPART(0)
CX       NDM_ESP(2)=SUM(NPART(1:NL_2))

       DO I=1, N_ESP
        WRITE(*,*) 'Particles of especie ', I,'= ',NDM_ESP(I)
       END DO

       IF (N_DM.NE.SUM(NDM_ESP(1:N_ESP))) THEN
        WRITE(*,*) 'WARNING: ESPECIES!'
        STOP
       END IF
***********************

       WRITE(*,*) '///DM///'
       WRITE(*,*) MAXVAL(RXPA(1:N_DM)),
     &            MINVAL(RXPA(1:N_DM))
       WRITE(*,*) MAXVAL(RYPA(1:N_DM)),
     &            MINVAL(RYPA(1:N_DM))
       WRITE(*,*) MAXVAL(RZPA(1:N_DM)),
     &            MINVAL(RZPA(1:N_DM))
       WRITE(*,*) MAXVAL(MASAP(1:N_DM)),
     &            MINVAL(MASAP(1:N_DM))

       MAXMAP=0.0
       MAXMAP=MAXVAL(MASAP(1:N_DM))
       MASAP=MASAP/MAXMAP

       NPATCH=0
       UBAS=0.0
       UBAS2=0.0
*------------------------------------------*
       END IF  !FLAG_MASCLET
*------------------------------------------*

*///////////////////////////////////
*      FIN LEEMOS FICHERO EXTERNO CON PARTI
*//////////////////////////////////


***//AQUI EMPIEZA LA CONSTRUCCION DE LA MALLA!!!

*////////////////////////////////////////////////////
*      CALCULO UP1=NUM.PART./CELDA DEL NIVEL BASE
*      !!IR=0!!
*////////////////////////////////////////////////////

       UP1=0.0
       UP11=0.0

       MASA_TEMP0=0.0
       KONTA=0

*-----
       DO IESP=N_ESP, 1, -1
*-----

       KK1=0
       KK2=0
       IF (IESP.EQ.1) THEN
        KK1=0
        KK2=SUM(NDM_ESP(1:IESP))
       ELSE
        KK1=SUM(NDM_ESP(1:IESP-1))
        KK2=SUM(NDM_ESP(1:IESP))
       END IF

       WRITE(*,*)'CHECKING SPECIES', KK1+1, KK2

       DO I=KK1+1,KK2

       IX=INT(((RXPA(I)-RADX(1))/DX)+0.5)+1
       JY=INT(((RYPA(I)-RADY(1))/DY)+0.5)+1
       KZ=INT(((RZPA(I)-RADZ(1))/DZ)+0.5)+1

       IF (IX.EQ.NX+1) IX=NX
       IF (JY.EQ.NY+1) JY=NY
       IF (KZ.EQ.NZ+1) KZ=NZ


       IF(IX.LT.1.OR.JY.LT.1.OR.KZ.LT.1)THEN
        WRITE(*,*) 'WARNING:'
        WRITE(*,*) IX,JY,KZ,IR
        WRITE(*,*) I,RXPA(I),RYPA(I),RZPA(I)
        STOP
       END IF

       IF(IX.GT.NX.OR.JY.GT.NY.OR.KZ.GT.NZ)THEN
        WRITE(*,*) 'WARNING_2:'
        WRITE(*,*) IX,JY,KZ,IR
        WRITE(*,*) I,RXPA(I),RYPA(I),RZPA(I)
        STOP
       END IF

       UP1(IX,JY,KZ,IESP)=UP1(IX,JY,KZ,IESP)+1.0

       MASA_TEMP0(IX,JY,KZ)=MASA_TEMP0(IX,JY,KZ)+MASAP(I)

       KONTA=KONTA+1


       END DO

       write(*,*) 'up1=',IESP, MAXVAL(UP1(1:NX,1:NY,1:NZ,IESP))
     &                        ,MINVAL(UP1(1:NX,1:NY,1:NZ,IESP))


*-----
       END DO   !ESPECIES!!!
*-----

       U1(1:NX,1:NY,1:NZ)=MASA_TEMP0(1:NX,1:NY,1:NZ)/(DX*DY*DZ)


       write(*,*) 'u1=', MAXVAL(U1(1:NX,1:NY,1:NZ))
     &                  ,MINVAL(U1(1:NX,1:NY,1:NZ))
       WRITE(*,*)'KONTA TRAS UP1!!', KONTA
       KONTA=0
*///////////////////////////


*/////////GAS//////////////////
*      CALCULO UP1=NUM.PART./CELDA DEL NIVEL BASE
*      !!IR=0!!
*///////////////////////////

       IF (N_GAS.GT.0) THEN

       UBAS=0.0

       DO I=1,N_GAS

       IX=INT(((RXG(I)-RADX(1))/DX)+0.5)+1
       JY=INT(((RYG(I)-RADY(1))/DY)+0.5)+1
       KZ=INT(((RZG(I)-RADZ(1))/DZ)+0.5)+1

       IF (IX.EQ.NX+1) IX=NX
       IF (JY.EQ.NY+1) JY=NY
       IF (KZ.EQ.NZ+1) KZ=NZ

       IF(IX.LT.1.OR.JY.LT.1.OR.KZ.LT.1)THEN
        WRITE(*,*) 'WARNING: GAS!'
        WRITE(*,*) IX,JY,KZ,IR
        WRITE(*,*) I,RXG(I),RYG(I),RZG(I)
        STOP
       END IF

       IF(IX.GT.NX.OR.JY.GT.NY.OR.KZ.GT.NZ)THEN
        WRITE(*,*) 'WARNING_2:GAS!'
        WRITE(*,*) IX,JY,KZ,IR
        WRITE(*,*) I,RXG(I),RYG(I),RZG(I)
        STOP
       END IF

       U1G(IX,JY,KZ)=U1G(IX,JY,KZ)+1.0     !MASAG(I)

       UBAS(IX,JY,KZ)=UBAS(IX,JY,KZ)+1.0
       END DO

       write(*,*) 'u1g=', MAXVAL(U1G(1:NX,1:NY,1:NZ))
     &                   ,MINVAL(U1G(1:NX,1:NY,1:NZ))

       U1G(1:NX,1:NY,1:NZ)=U1G(1:NX,1:NY,1:NZ)
     &                     *MASAG(1)/(DX*DY*DZ)

       write(*,*) 'u1GAS=', MAXVAL(U1G(1:NX,1:NY,1:NZ))
     &                     ,MINVAL(U1G(1:NX,1:NY,1:NZ))

       END IF    !N_GAS>0

*///////////////////////////


*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       BOR=1      !OJO!!!!!!!!!!
*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*      1 LEVEL , REFINMENT DX/2

*****************************
       IR=1
*****************************
       KONTA=0

       DXPA=0.0
       DYPA=0.0
       DZPA=0.0
       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)

       INCRE=1
       NPATCH(IR)=0
       IPATCH=0
*******************

       NPX=NAMRX
       NPY=NAMRY
       NPZ=NAMRZ

       N1=NX
       N2=NY
       N3=NZ


*-----
       DO IESP=N_ESP, 1, -1
*-----
       IPATCH_ESP=0

       U1MIN=0.0

!$OMP  PARALLEL DO SHARED(NX,NY,NZ,UBAS,U1MIN),PRIVATE(I,J,K)
       DO K=1,NZ
       DO J=1,NY
       DO I=1,NX
        UBAS(I,J,K)=U1MIN
       END DO
       END DO
       END DO


!$OMP  PARALLEL DO SHARED(NX,NY,NZ,UBAS,UP1,BOR,IESP),
!$OMP+         PRIVATE(I,J,K)
       DO K=BOR,NZ-BOR+1
       DO J=BOR,NY-BOR+1
       DO I=BOR,NX-BOR+1
        UBAS(I,J,K)=UP1(I,J,K,IESP)
       END DO
       END DO
       END DO

       NL1=0
*!$OMP  PARALLEL DO SHARED(CR0,BOR,UBAS,COTA,IR,NX,NY,NZ,IESP),
*!$OMP+        PRIVATE(I,J,K),REDUCTION(+:NL1)
       DO K=BOR,NZ-BOR+1
       DO J=BOR,NY-BOR+1
       DO I=BOR,NX-BOR+1
       IF (UBAS(I,J,K).GT.COTA(1,IR)) THEN
        NL1=NL1+1
        CR0(I,J,K,IESP)=NL1
       ENDIF
       END DO
       END DO
       END DO



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO WHILE (NL1.GT.0.AND.IPATCH.LT.NPALEV)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       NL10=NL1

       INMAX=MAXLOC(UBAS)

       IX=INMAX(1)
       JY=INMAX(2)
       KZ=INMAX(3)

       PAX1=IX-1   !!OJO!! AKI HABIA +/-2
       PAX2=IX+1
       PAY1=JY-1
       PAY2=JY+1
       PAZ1=KZ-1
       PAZ2=KZ+1

       IF (PAX1.LT.BOR) PAX1=BOR
       IF (PAY1.LT.BOR) PAY1=BOR
       IF (PAZ1.LT.BOR) PAZ1=BOR
       IF (PAX2.GT.(NX-BOR+1)) PAX2=NX-BOR+1
       IF (PAY2.GT.(NY-BOR+1)) PAY2=NY-BOR+1
       IF (PAZ2.GT.(NZ-BOR+1)) PAZ2=NZ-BOR+1


       NBAS=MAX(PAX2-PAX1,PAY2-PAY1,PAZ2-PAZ1)
       NBAS=2*(NBAS+1)

*      SHRINKING THE PATCH

       NEF=0
       NEF=COUNT(CR0(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
       NCELL=(PAX2-PAX1+1)*(PAY2-PAY1+1)*(PAZ2-PAZ1+1)

       MARCA=0
       DO WHILE (NBAS.LT.NPX.AND.MARCA.NE.1)

        NEF1=0
        NEF1=COUNT(CR0(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
        NEF2=0
        NEF2=COUNT(CR0(PAX1:PAX2,PAY1:PAY2,PAZ1-INCRE:PAZ2,IESP)
     &                                                        > 0 )
        NCELL=(PAX2-PAX1+1)*(PAY2-PAY1+1)*(PAZ2-(PAZ1-INCRE)+1)
        IF (NEF2.GT.NEF1) PAZ1=PAZ1-INCRE
        IF (2*(PAZ2-PAZ1+1).GT.NPZ)  PAZ1=PAZ1+INCRE
        PAZ1=MAX(PAZ1,BOR)

        NEF1=0
        NEF1=COUNT(CR0(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
        NEF2=0
        NEF2=COUNT(CR0(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2+INCRE,IESP)
     &                                                        > 0 )
        NCELL=(PAX2-PAX1+1)*(PAY2-PAY1+1)*((PAZ2+INCRE)-PAZ1+1)
        IF (NEF2.GT.NEF1) PAZ2=PAZ2+INCRE
        IF (2*(PAZ2-PAZ1+1).GT.NPZ)  PAZ2=PAZ2-INCRE
        PAZ2=MIN(PAZ2,N3-BOR+1)

        NEF1=0
        NEF1=COUNT(CR0(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
        NEF2=0
        NEF2=COUNT(CR0(PAX1:PAX2,PAY1-INCRE:PAY2,PAZ1:PAZ2,IESP)
     &                                                        > 0 )
        NCELL=(PAX2-PAX1+1)*(PAY2-(PAY1-INCRE)+1)*(PAZ2-PAZ1+1)
        IF (NEF2.GT.NEF1) PAY1=PAY1-INCRE
        IF (2*(PAY2-PAY1+1).GT.NPY)  PAY1=PAY1+INCRE
        PAY1=MAX(PAY1,BOR)

        NEF1=0
        NEF1=COUNT(CR0(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
        NEF2=0
        NEF2=COUNT(CR0(PAX1:PAX2,PAY1:PAY2+INCRE,PAZ1:PAZ2,IESP)
     &                                                        > 0 )
        NCELL=(PAX2-PAX1+1)*((PAY2+INCRE)-PAY1+1)*(PAZ2-PAZ1+1)
        IF (NEF2.GT.NEF1) PAY2=PAY2+INCRE
        IF (2*(PAY2-PAY1+1).GT.NPY) PAY2=PAY2-INCRE
        PAY2=MIN(PAY2,N2-BOR+1)

        NEF1=0
        NEF1=COUNT(CR0(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
        NEF2=0
        NEF2=COUNT(CR0(PAX1-INCRE:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP)
     &                                                        > 0 )
        NCELL=(PAX2-(PAX1-INCRE)+1)*(PAY2-PAY1+1)*(PAZ2-PAZ1+1)
        IF (NEF2.GT.NEF1) PAX1=PAX1-INCRE
        IF (2*(PAX2-PAX1+1).GT.NPX)   PAX1=PAX1+INCRE
        PAX1=MAX(PAX1,BOR)

        NEF1=0
        NEF1=COUNT(CR0(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
        NEF2=0
        NEF2=COUNT(CR0(PAX1:PAX2+INCRE,PAY1:PAY2,PAZ1:PAZ2,IESP)
     &                                                        > 0 )
        NCELL=((PAX2+INCRE)-PAX1+1)*(PAY2-PAY1+1)*(PAZ2-PAZ1+1)
        IF (NEF2.GT.NEF1) PAX2=PAX2+INCRE
        IF (2*(PAX2-PAX1+1).GT.NPX)  PAX2=PAX2-INCRE
        PAX2=MIN(PAX2,N1-BOR+1)


*       FINAL EFFICENCY

        NEF1=0
        NEF1=COUNT(CR0(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
        NCELL=(PAX2-PAX1+1)*(PAY2-PAY1+1)*(PAZ2-PAZ1+1)

        IF (NEF1.LE.NEF) THEN
         IF (MARCA.EQ.0) THEN
          MARCA=2
         ELSE
          MARCA=1
         END IF
        END IF

        NBAS=MAX(PAX2-PAX1,PAY2-PAY1,PAZ2-PAZ1)
        NBAS=2*(NBAS+1)

        NEF=NEF1
       END DO

*       OVERSIZE CONTROL
        IF (PAX1.LT.BOR) PAX1=BOR
        IF (PAX2.GT.(NX-BOR+1)) PAX2=NX-BOR+1
        IF (PAY1.LT.BOR) PAY1=BOR
        IF (PAY2.GT.(NY-BOR+1)) PAY2=NY-BOR+1
        IF (PAZ1.LT.BOR) PAZ1=BOR
        IF (PAZ2.GT.(NZ-BOR+1)) PAZ2=NZ-BOR+1

        NBAS=MAX(PAX2-PAX1,PAY2-PAY1,PAZ2-PAZ1)
        NBAS=2*(NBAS+1)
        NBIS=MIN(PAX2-PAX1,PAY2-PAY1,PAZ2-PAZ1)
        NBIS=2*(NBAS+1)


         IPATCH_ESP=IPATCH_ESP+1
         IPATCH=IPATCH+1

         IF (PARCHLIM.NE.0.AND.IPATCH.GT.MPAPOLEV(IR)) THEN
          IPATCH=IPATCH-1
          WRITE(*,*) 'exit por PARCHLIM.NE.0.AND.IPATCH.GT.MPAPOLEV(IR)'
          EXIT
         ENDIF


         PATCHNX(IPATCH)=2*(PAX2-PAX1+1)
         PATCHNY(IPATCH)=2*(PAY2-PAY1+1)
         PATCHNZ(IPATCH)=2*(PAZ2-PAZ1+1)
*        LEFT-BOTTOM LIMIT OF THE RECTANGLE
         PATCHX(IPATCH)=PAX1
         PATCHY(IPATCH)=PAY1
         PATCHZ(IPATCH)=PAZ1
         PATCHRX(IPATCH)=RADX(PAX1)
         PATCHRY(IPATCH)=RADY(PAY1)
         PATCHRZ(IPATCH)=RADZ(PAZ1)
         PARE(IPATCH)=0


         NEF=0
         DO K=PAZ1,PAZ2
         DO J=PAY1,PAY2
         DO I=PAX1,PAX2
           IF (CR0(I,J,K,IESP).NE.0) NEF=NEF+1
           CR0(I,J,K,IESP)=0
           UBAS(I,J,K)=U1MIN
         END DO
         END DO
         END DO
         NL1=NL1-NEF

       IF (NL10.EQ.NL1) EXIT

!!!!!!!!!!!!!!!!!!!
       END DO
!!!!!!!!!!!!!!!!!!!

CX       WRITE(*,*) 'IESP, IPATCH',IESP,IPATCH_ESP,IPATCH

*-----
       END DO     !IESP
*-----


*      TOTAL NUMBER OF PATCHES
       NPATCH(IR)=IPATCH
       WRITE(*,*) 'Total number of patches',ipatch
       WRITE(*,*) 'Starting interpolation...'



*      INTERPOLACION

       NESP_BAS=N_ESP
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))


!$OMP  PARALLEL DO SHARED(UP11,UP1,NESP_BAS,
!$OMP+     PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
!$OMP+     PARE,PATCHRX,PATCHRY,PATCHRZ,MASAP,LOW1,LOW2,
!$OMP+     DX,DY,DZ,NX,NY,NZ,NPARTICULAS,NPART,NDM_ESP,
!$OMP+     RXPA,RYPA,RZPA,U1,U11,DXPA,DYPA,DZPA),
!$OMP+   PRIVATE(I,N1,N2,N3,KZ,JY,IX,II,JJ,KK,I1,J1,K1,L1,L2,L3,
!$OMP+     IP,NP1,NP2,NP3,MASA_TEMP,KONTA,IESP,KK1,KK2)
       DO I=LOW1,LOW2

       NP1=NX
       NP2=NY
       NP3=NZ

       L1=PATCHX(I)
       L2=PATCHY(I)
       L3=PATCHZ(I)

       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       MASA_TEMP=0.0
       KONTA=0

*-----
       DO IESP=NESP_BAS, 1, -1
*-----

       KK1=0
       KK2=0
       IF (IESP.EQ.1) THEN
        KK1=0
        KK2=SUM(NDM_ESP(1:IESP))
       ELSE
        KK1=SUM(NDM_ESP(1:IESP-1))
        KK2=SUM(NDM_ESP(1:IESP))
       END IF

       DO IP=KK1+1,KK2

       IX=INT(((RXPA(IP)-(PATCHRX(I)-0.5*DXPA))/DXPA)+0.5)+1
       JY=INT(((RYPA(IP)-(PATCHRY(I)-0.5*DYPA))/DYPA)+0.5)+1
       KZ=INT(((RZPA(IP)-(PATCHRZ(I)-0.5*DZPA))/DZPA)+0.5)+1


       IF(IX.LE.N1.AND.IX.GE.1.AND.JY.LE.N2.AND.
     &    JY.GE.1.AND.KZ.LE.N3.AND.KZ.GE.1) THEN

       UP11(IX,JY,KZ,I,IESP)=UP11(IX,JY,KZ,I,IESP)+1.0

       MASA_TEMP(IX,JY,KZ)=MASA_TEMP(IX,JY,KZ)+MASAP(IP)
       KONTA=KONTA+1

       END IF

       END DO


*-----
       END DO !IESP
*-----


       U11(1:N1,1:N2,1:N3,I)=MASA_TEMP(1:N1,1:N2,1:N3)
     &                          /(DXPA*DYPA*DZPA)

*//////////////GAS/////////
CX       IF (N_GAS.GT.0) THEN

CX       MASA_TEMP=0.0
CX       KONTA=0

CX       DO IP=1,N_GAS

CX       IX=INT(((RXG(IP)-(PATCHRX(IR,I)-0.5*DXPA))/DXPA)+0.5)+1
CX       JY=INT(((RYG(IP)-(PATCHRY(IR,I)-0.5*DYPA))/DYPA)+0.5)+1
CX       KZ=INT(((RZG(IP)-(PATCHRZ(IR,I)-0.5*DZPA))/DZPA)+0.5)+1


CX       IF(IX.LE.N1.AND.IX.GE.1.AND.JY.LE.N2.AND.
CX     &    JY.GE.1.AND.KZ.LE.N3.AND.KZ.GE.1) THEN
CX
CX       U11G(IX,JY,KZ,IR,I)=U11G(IX,JY,KZ,IR,I)+1.0
CX       U12(IX,JY,KZ,IR,I)=U12(IX,JY,KZ,IR,I)+VXG(IP)*MASAG(IP)
CX       U13(IX,JY,KZ,IR,I)=U13(IX,JY,KZ,IR,I)+VYG(IP)*MASAG(IP)
CX       U14(IX,JY,KZ,IR,I)=U14(IX,JY,KZ,IR,I)+VZG(IP)*MASAG(IP)

CX       MASA_TEMP(IX,JY,KZ)=MASA_TEMP(IX,JY,KZ)+MASAG(IP)
CX       KONTA=KONTA+1

CX       END IF

CX       END DO


CX       U11G(1:N1,1:N2,1:N3,IR,I)=MASA_TEMP(1:N1,1:N2,1:N3)
CX     &                           /(DXPA*DYPA*DZPA)
CX       U12(1:N1,1:N2,1:N3,IR,I)=U12(1:N1,1:N2,1:N3,IR,I)
CX     &                          /MASA_TEMP(1:N1,1:N2,1:N3)
CX       U13(1:N1,1:N2,1:N3,IR,I)=U13(1:N1,1:N2,1:N3,IR,I)
CX     &                          /MASA_TEMP(1:N1,1:N2,1:N3)
CX       U14(1:N1,1:N2,1:N3,IR,I)=U14(1:N1,1:N2,1:N3,IR,I)
CX     &                          /MASA_TEMP(1:N1,1:N2,1:N3)

CX       END IF
*///////////////////////


       END DO


       WRITE(*,*) 'LEVEL=',IR,' patches=',NPATCH(IR)
*************************
*************************

*------------------------------*
********************************
*  fin nivel 1...otros niveles *
********************************
*------------------------------*

       WRITE(*,*) '*----------LEVELS!!! FROM 2 UP TO NL------------*'

       NPX=NAMRX
       NPY=NAMRY
       NPZ=NAMRZ

       BOR=1  !OJO

*****************************
       DO IR=2,NL
*****************************
       KONTA=0

       DXPA=0.0
       DYPA=0.0
       DZPA=0.0
       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)

       NPATCH(IR)=0
*********************************
       IF (NPATCH(IR-1).GT.0) THEN
*********************************

*! variables auxiliares paralelizacion
       ALLOCATE(LNPATCH(NPATCH(IR-1)))          !num de parches q saldran de cada uno de IR-1
       ALLOCATE(LPATCHNX(NPALEV3,NPATCH(IR-1)))
       ALLOCATE(LPATCHNY(NPALEV3,NPATCH(IR-1)))
       ALLOCATE(LPATCHNZ(NPALEV3,NPATCH(IR-1)))
       ALLOCATE(LPATCHX(NPALEV3,NPATCH(IR-1)))
       ALLOCATE(LPATCHY(NPALEV3,NPATCH(IR-1)))
       ALLOCATE(LPATCHZ(NPALEV3,NPATCH(IR-1)))
       ALLOCATE(LPATCHRX(NPALEV3,NPATCH(IR-1)))
       ALLOCATE(LPATCHRY(NPALEV3,NPATCH(IR-1)))
       ALLOCATE(LPATCHRZ(NPALEV3,NPATCH(IR-1)))
       ALLOCATE(LPARE(NPALEV3,NPATCH(IR-1)))



!$OMP  PARALLEL DO SHARED(IR,NPATCH,LNPATCH,LPATCHNX,LPATCHNY,
!$OMP+        LPATCHNZ,LPATCHX,LPATCHY,LPATCHZ,LPATCHRX,
!$OMP+        LPATCHRY,LPATCHRZ,LPARE),PRIVATE(I)
       DO I=1,NPATCH(IR-1)
        LNPATCH(I)=0
        LPATCHNX(:,I)=0
        LPATCHNY(:,I)=0
        LPATCHNZ(:,I)=0
        LPATCHX(:,I)=0
        LPATCHY(:,I)=0
        LPATCHZ(:,I)=0
        LPATCHRX(:,I)=0.0
        LPATCHRY(:,I)=0.0
        LPATCHRZ(:,I)=0.0
        LPARE(:,I)=0
       END DO

*       write(*,*) 'prueba 1'

*      INCRE, PATCH EXTENSION UNITS
        INCRE=1

       NESP_BAS=N_ESP

*      OJO A LA PARELIZACION
!$OMP  PARALLEL DO SHARED(NPATCH,IR,BOR,PATCHNZ,
!$OMP+     PATCHNY,PATCHNX,
!$OMP+     COTA,UP11,NL,PATCHRX,PATCHRY,PATCHRZ,DX,DY,DZ,
!$OMP+     NPX,NPY,NPZ,PARE,PATCHX,PATCHY,PATCHZ,
!$OMP+     LNPATCH,LPATCHNX,LPATCHNY,LPATCHNZ,LPATCHX,LPATCHY,
!$OMP+     LPATCHZ,LPATCHRX,LPATCHRY,LPATCHRZ,LPARE,INCRE,
!$OMP+     NESP_BAS),
!$OMP+ PRIVATE(NL1,NL10,K,I,J,N1,N2,N3,UBAS,U1MIN,IPALE,
!$OMP+     INMAX,IX,JY,KZ,PAX1,PAX2,PAY1,PAY2,PAZ1,PAZ2,
!$OMP+     NBAS,NEF,NEF1,NEF2,NCELL,
!$OMP+     IPATCH,MARCA,NBIS,
!$OMP+     IESP),
!$OMP+  FIRSTPRIVATE(CR02,UBAS2)

*SE PUEDE PAralleliza el do wwhile siguiente???!!

*----------------------------
       DO IPALE2=1,NPATCH(IR-1)
*----------------------------

       IPALE=SUM(NPATCH(0:IR-2)) + IPALE2     !ipale de (ir-1)

       IPATCH=0

C*-----
       DO IESP=NESP_BAS, 1, -1
C*-----
       IPATCH_ESP=0

       NBAS=MIN(PATCHNZ(IPALE),PATCHNY(IPALE),PATCHNX(IPALE))

       N1=PATCHNX(IPALE)
       N2=PATCHNY(IPALE)
       N3=PATCHNZ(IPALE)

       U1MIN=0.0
       UBAS2=U1MIN

       UBAS2(BOR:N1-BOR+1,BOR:N2-BOR+1,BOR:N3-BOR+1)=
     & UP11(BOR:N1-BOR+1,BOR:N2-BOR+1,BOR:N3-BOR+1,IPALE,IESP)

*      OVER avoids overlapping
       NL1=0
       DO K=BOR,PATCHNZ(IPALE)-BOR+1
       DO J=BOR,PATCHNY(IPALE)-BOR+1
       DO I=BOR,PATCHNX(IPALE)-BOR+1
       IF (UBAS2(I,J,K).GT.COTA(1,IR)) THEN
        NL1=NL1+1
        CR02(I,J,K,IESP)=NL1
       ENDIF
       END DO
       END DO
       END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO WHILE (NL1.GT.0)
!!!!!!!!!!!!!!!!!!!!!!!!!!

       NL10=NL1

       INMAX=MAXLOC(UBAS2)

       IX=INMAX(1)
       JY=INMAX(2)
       KZ=INMAX(3)

       PAX1=IX-1         !AKI HABIA 4
       PAX2=IX+1
       PAY1=JY-1
       PAY2=JY+1
       PAZ1=KZ-1
       PAZ2=KZ+1

       IF (PAX1.LT.BOR) PAX1=BOR
       IF (PAY1.LT.BOR) PAY1=BOR
       IF (PAZ1.LT.BOR) PAZ1=BOR
       IF (PAX2.GT.(N1-BOR+1)) PAX2=N1-BOR+1
       IF (PAY2.GT.(N2-BOR+1)) PAY2=N2-BOR+1
       IF (PAZ2.GT.(N3-BOR+1)) PAZ2=N3-BOR+1


       NBAS=MAX(PAX2-PAX1,PAY2-PAY1,PAZ2-PAZ1)
       NBAS=2*(NBAS+1)

*      SHRINKING THE PATCH

       NEF=0
       NEF=COUNT(CR02(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
       NCELL=(PAX2-PAX1+1)*(PAY2-PAY1+1)*(PAZ2-PAZ1+1)

       MARCA=0
       DO WHILE (NBAS.LT.NPX.AND.MARCA.NE.1)

        NEF1=0
        NEF1=COUNT(CR02(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
        NEF2=0
        NEF2=COUNT(CR02(PAX1:PAX2,PAY1:PAY2,PAZ1-INCRE:PAZ2,IESP)
     &                                                         > 0 )
        NCELL=(PAX2-PAX1+1)*(PAY2-PAY1+1)*(PAZ2-(PAZ1-INCRE)+1)
        IF (NEF2.GT.NEF1) PAZ1=PAZ1-INCRE
        IF (2*(PAZ2-PAZ1+1).GT.NPZ)  PAZ1=PAZ1+INCRE
        PAZ1=MAX(PAZ1,BOR)

        NEF1=0
        NEF1=COUNT(CR02(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
        NEF2=0
        NEF2=COUNT(CR02(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2+INCRE,IESP)
     &                                                         > 0 )
        NCELL=(PAX2-PAX1+1)*(PAY2-PAY1+1)*((PAZ2+INCRE)-PAZ1+1)
        IF (NEF2.GT.NEF1) PAZ2=PAZ2+INCRE
        IF (2*(PAZ2-PAZ1+1).GT.NPZ)  PAZ2=PAZ2-INCRE
        PAZ2=MIN(PAZ2,N3-BOR+1)


        NEF1=0
        NEF1=COUNT(CR02(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
        NEF2=0
        NEF2=COUNT(CR02(PAX1:PAX2,PAY1-INCRE:PAY2,PAZ1:PAZ2,IESP)
     &                                                         > 0 )
        NCELL=(PAX2-PAX1+1)*(PAY2-(PAY1-INCRE)+1)*(PAZ2-PAZ1+1)
        IF (NEF2.GT.NEF1) PAY1=PAY1-INCRE
        IF (2*(PAY2-PAY1+1).GT.NPY)  PAY1=PAY1+INCRE
        PAY1=MAX(PAY1,BOR)


        NEF1=0
        NEF1=COUNT(CR02(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
        NEF2=0
        NEF2=COUNT(CR02(PAX1:PAX2,PAY1:PAY2+INCRE,PAZ1:PAZ2,IESP)
     &                                                         > 0 )
        NCELL=(PAX2-PAX1+1)*((PAY2+INCRE)-PAY1+1)*(PAZ2-PAZ1+1)
        IF (NEF2.GT.NEF1) PAY2=PAY2+INCRE
        IF (2*(PAY2-PAY1+1).GT.NPY)  PAY2=PAY2-INCRE
        PAY2=MIN(PAY2,N2-BOR+1)


        NEF1=0
        NEF1=COUNT(CR02(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
        NEF2=0
        NEF2=COUNT(CR02(PAX1-INCRE:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP)
     &                                                         > 0 )
        NCELL=(PAX2-(PAX1-INCRE)+1)*(PAY2-PAY1+1)*(PAZ2-PAZ1+1)
        IF (NEF2.GT.NEF1) PAX1=PAX1-INCRE
        IF (2*(PAX2-PAX1+1).GT.NPX)  PAX1=PAX1+INCRE
        PAX1=MAX(PAX1,BOR)


        NEF1=0
        NEF1=COUNT(CR02(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
        NEF2=0
        NEF2=COUNT(CR02(PAX1:PAX2+INCRE,PAY1:PAY2,PAZ1:PAZ2,IESP)
     &                                                         > 0 )
        NCELL=((PAX2+INCRE)-PAX1+1)*(PAY2-PAY1+1)*(PAZ2-PAZ1+1)
        IF (NEF2.GT.NEF1) PAX2=PAX2+INCRE
        IF (2*(PAX2-PAX1+1).GT.NPX)  PAX2=PAX2-INCRE
        PAX2=MIN(PAX2,N1-BOR+1)


*       FINAL EFFICENCY

        NEF1=0
        NEF1=COUNT(CR02(PAX1:PAX2,PAY1:PAY2,PAZ1:PAZ2,IESP) > 0 )
        NCELL=(PAX2-PAX1+1)*(PAY2-PAY1+1)*(PAZ2-PAZ1+1)

        IF (NEF1.LE.NEF) THEN
         IF (MARCA.EQ.0) THEN
          MARCA=2
         ELSE
          MARCA=1
         END IF
        END IF

        NBAS=MAX(PAX2-PAX1,PAY2-PAY1,PAZ2-PAZ1)
        NBAS=2*(NBAS+1)

        NEF=NEF1
       END DO    !while (NBAS.LT.NPX.AND.MARCA.NE.1)


*       OVERSIZE CONTROL
        IF (PAX1.LT.BOR) PAX1=BOR
        IF (PAX2.GT.(N1-BOR+1)) PAX2=PATCHNX(IPALE)-BOR+1
        IF (PAY1.LT.BOR) PAY1=BOR
        IF (PAY2.GT.(N2-BOR+1)) PAY2=PATCHNY(IPALE)-BOR+1
        IF (PAZ1.LT.BOR) PAZ1=BOR
        IF (PAZ2.GT.(N3-BOR+1)) PAZ2=PATCHNZ(IPALE)-BOR+1

        NBAS=MAX(PAX2-PAX1,PAY2-PAY1,PAZ2-PAZ1)
        NBAS=2*(NBAS+1)
        NBIS=MIN(PAX2-PAX1,PAY2-PAY1,PAZ2-PAZ1)
        NBIS=2*NBIS


        IPATCH_ESP=IPATCH_ESP+1
        IPATCH=IPATCH+1                !parches dentro de este parche

        IF (IPATCH.GT.NPALEV3) THEN
         WRITE(*,*) 'WARNING: ipatch > npalev3', ipatch,npalev3,ipale
*        STOP
        END IF

        LPATCHNX(IPATCH,IPALE2)=2*(PAX2-PAX1+1)
        LPATCHNY(IPATCH,IPALE2)=2*(PAY2-PAY1+1)
        LPATCHNZ(IPATCH,IPALE2)=2*(PAZ2-PAZ1+1)

*       LEFT-BOTTOM LIMIT OF THE RECTANGLE
        LPATCHX(IPATCH,IPALE2)=PAX1
        LPATCHY(IPATCH,IPALE2)=PAY1
        LPATCHZ(IPATCH,IPALE2)=PAZ1

        LPATCHRX(IPATCH,IPALE2)=FLOAT((PAX1-1))*(DX/(2**(IR-1)))+
     &                     PATCHRX(IPALE)-0.5*(DX/(2**(IR-1)))
        LPATCHRY(IPATCH,IPALE2)=FLOAT((PAY1-1))*(DY/(2**(IR-1)))+
     &                     PATCHRY(IPALE)-0.5*(DY/(2**(IR-1)))
        LPATCHRZ(IPATCH,IPALE2)=FLOAT((PAZ1-1))*(DZ/(2**(IR-1)))+
     &                     PATCHRZ(IPALE)-0.5*(DZ/(2**(IR-1)))

*       PARENT PATCH
        LPARE(IPATCH,IPALE2)=IPALE
        LNPATCH(IPALE2)=IPATCH

        IF (LPATCHNX(IPATCH,IPALE2).GT.NPX)
     &     WRITE(*,*) 'WARNING: PARCHE X DEMASIADO GRANDE',IR
        IF (LPATCHNY(IPATCH,IPALE2).GT.NPY)
     &     WRITE(*,*) 'WARNING: PARCHE Y DEMASIADO GRANDE',IR
        IF (LPATCHNZ(IPATCH,IPALE2).GT.NPZ)
     &     WRITE(*,*) 'WARNING: PARCHE Z DEMASIADO GRANDE',IR



        NEF=0
        DO K=PAZ1,PAZ2
        DO J=PAY1,PAY2
        DO I=PAX1,PAX2
          IF (CR02(I,J,K,IESP).NE.0) NEF=NEF+1
          CR02(I,J,K,IESP)=0
          UBAS2(I,J,K)=U1MIN
        END DO
        END DO
        END DO
        NL1=NL1-NEF

!!!!!!!!!!!!!!!!!!!
       END DO     !while(nl1.gt.0)
!!!!!!!!!!!!!!!!!!!


*-----------------
       END DO    !IESP
*-----------------


*-----------------
       END DO    !IPALE
*-----------------

*****************************************************************
*      REUNIFICAMOS Y DIMENSIONAMOS TODOS LOS PARCHES DE IR
*****************************************************************

       IPATCH=SUM(NPATCH(0:IR-1))
       IPATCH2=0

       DO IPALE2=1,NPATCH(IR-1)
       DO IPA2=1,LNPATCH(IPALE2)         !num de parches en los q se divide IPALE2

        IPATCH=IPATCH+1
        IPATCH2=IPATCH2+1

*       si el numero total de parches es muy grande se detiene
        IF (IPATCH.GT.NPALEV) THEN
         IPATCH=IPATCH-1
         IPATCH2=IPATCH2-1
C         WRITE(*,*) 'WARNING IN PATCH NUMBER',IPATCH,NPALEV
         EXIT
        END IF

*       si el numero de parches en este nivel es muy grande se detiene
        IF (PARCHLIM.NE.0.AND.IPATCH2.GT.MPAPOLEV(IR)) THEN
         IPATCH=IPATCH-1
         IPATCH2=IPATCH2-1
C         WRITE(*,*) 'WARNING IN MPAPOLEV(IR)',IPATCH2
         EXIT
        END IF

         PATCHNX(IPATCH)=LPATCHNX(IPA2,IPALE2)
         PATCHNY(IPATCH)=LPATCHNY(IPA2,IPALE2)
         PATCHNZ(IPATCH)=LPATCHNZ(IPA2,IPALE2)

         PATCHX(IPATCH)=LPATCHX(IPA2,IPALE2)
         PATCHY(IPATCH)=LPATCHY(IPA2,IPALE2)
         PATCHZ(IPATCH)=LPATCHZ(IPA2,IPALE2)

         PATCHRX(IPATCH)=LPATCHRX(IPA2,IPALE2)
         PATCHRY(IPATCH)=LPATCHRY(IPA2,IPALE2)
         PATCHRZ(IPATCH)=LPATCHRZ(IPA2,IPALE2)

         PARE(IPATCH)=LPARE(IPA2,IPALE2)


       END DO
       END DO

*      TOTAL PATCHES NUMBER AT A FIXED LEVEL
       NPATCH(IR)=IPATCH2

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))

*      INTERPOLACION
       NESP_BAS=N_ESP
cx*!$OMP  PARALLEL DO SHARED(IR,NPATCH,UP11,NESP_BAS,
cx*!$OMP+     PATCHNX,PATCHNY,PATCHNZ,PARE,PATCHX,PATCHY,PATCHZ,
cx*!$OMP+     PATCHRX,PATCHRY,PATCHRZ,
cx*!$OMP+     DX,DY,DZ,NPARTICULAS,NPART,NDM_ESP,LOW1,LOW2,
cx*!$OMP+     RXPA,RYPA,RZPA,U1,U11,U12,U13,U14,DXPA,DYPA,DZPA),
cx*!$OMP+   PRIVATE(I,N1,N2,N3,KZ,JY,IX,II,JJ,KK,I1,J1,K1,L1,L2,L3,
cx*!$OMP+     IP,NP1,NP2,NP3,MASA_TEMP,KONTA,IESP,KK1,KK2)
*
*!$OMP  PARALLEL DO SHARED(IR,NPATCH,UP11,NESP_BAS,
*!$OMP+     PATCHNX,PATCHNY,PATCHNZ,PARE,PATCHX,PATCHY,PATCHZ,
*!$OMP+     PATCHRX,PATCHRY,PATCHRZ,
*!$OMP+     DX,DY,DZ,NPARTICULAS,NPART,NDM_ESP,LOW1,LOW2,
*!$OMP+     RXPA,RYPA,RZPA,U1,U11,DXPA,DYPA,DZPA),
*!$OMP+   PRIVATE(I,N1,N2,N3,KZ,JY,IX,II,JJ,KK,I1,J1,K1,L1,L2,L3,
*!$OMP+     IP,NP1,NP2,NP3,MASA_TEMP,KONTA,IESP,KK1,KK2)


       DO I=LOW1,LOW2

       NP1=PATCHNX(PARE(I))
       NP2=PATCHNY(PARE(I))
       NP3=PATCHNZ(PARE(I))

       L1=PATCHX(I)
       L2=PATCHY(I)
       L3=PATCHZ(I)
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       MASA_TEMP=0.0
       KONTA=0


*-----
       DO IESP=NESP_BAS, 1, -1
*-----

       KK1=0
       KK2=0

       IF (IESP.EQ.1) THEN
        KK1=0
        KK2=SUM(NDM_ESP(1:IESP))
       ELSE
        KK1=SUM(NDM_ESP(1:IESP-1))
        KK2=SUM(NDM_ESP(1:IESP))
       END IF


       DO IP=KK1+1,KK2

       IX=INT(((RXPA(IP)-(PATCHRX(I)-0.5*DXPA))/DXPA)+0.5)+1
       JY=INT(((RYPA(IP)-(PATCHRY(I)-0.5*DYPA))/DYPA)+0.5)+1
       KZ=INT(((RZPA(IP)-(PATCHRZ(I)-0.5*DZPA))/DZPA)+0.5)+1


       IF(IX.LE.N1.AND.IX.GE.1.AND.JY.LE.N2.AND.
     &    JY.GE.1.AND.KZ.LE.N3.AND.KZ.GE.1) THEN

       UP11(IX,JY,KZ,I,IESP)=UP11(IX,JY,KZ,I,IESP)+1.0

       MASA_TEMP(IX,JY,KZ)=MASA_TEMP(IX,JY,KZ)+MASAP(IP)
       KONTA=KONTA+1

       END IF

       END DO

*---
       END DO  !IESP!!!
*---


       U11(1:N1,1:N2,1:N3,I)=MASA_TEMP(1:N1,1:N2,1:N3)
     &                          /(DXPA*DYPA*DZPA)


*/////////////GAS
CX       IF (N_GAS.GT.0) THEN

CX       MASA_TEMP=0.0
CX       KONTA=0

CX       DO IP=1,N_GAS

CX       IX=INT(((RXG(IP)-(PATCHRX(IR,I)-0.5*DXPA))/DXPA)+0.5)+1
CX       JY=INT(((RYG(IP)-(PATCHRY(IR,I)-0.5*DYPA))/DYPA)+0.5)+1
CX       KZ=INT(((RZG(IP)-(PATCHRZ(IR,I)-0.5*DZPA))/DZPA)+0.5)+1


CX       IF(IX.LE.N1.AND.IX.GE.1.AND.JY.LE.N2.AND.
CX     &    JY.GE.1.AND.KZ.LE.N3.AND.KZ.GE.1) THEN

CX       U12(IX,JY,KZ,IR,I)=U12(IX,JY,KZ,IR,I)+VXG(IP)*MASAG(IP)
CX       U13(IX,JY,KZ,IR,I)=U13(IX,JY,KZ,IR,I)+VYG(IP)*MASAG(IP)
CX       U14(IX,JY,KZ,IR,I)=U14(IX,JY,KZ,IR,I)+VZG(IP)*MASAG(IP)
CX
CX       MASA_TEMP(IX,JY,KZ)=MASA_TEMP(IX,JY,KZ)+MASAG(IP)
CX       KONTA=KONTA+1

CX       END IF
CX
CX       END DO
CX

CX       U11G(1:N1,1:N2,1:N3,IR,I)=MASA_TEMP(1:N1,1:N2,1:N3)
CX     &                          /(DXPA*DYPA*DZPA)
CX       U12(1:N1,1:N2,1:N3,IR,I)=U12(1:N1,1:N2,1:N3,IR,I)
CX     &                           /MASA_TEMP(1:N1,1:N2,1:N3)
CX       U13(1:N1,1:N2,1:N3,IR,I)=U13(1:N1,1:N2,1:N3,IR,I)
CX     &                           /MASA_TEMP(1:N1,1:N2,1:N3)
CX       U14(1:N1,1:N2,1:N3,IR,I)=U14(1:N1,1:N2,1:N3,IR,I)
CX     &                           /MASA_TEMP(1:N1,1:N2,1:N3)


CX       END IF

       END DO

*////////////////////////////

       DEALLOCATE(LNPATCH)    ! variables auxiliares paralelizacion
       DEALLOCATE(LPATCHNX)
       DEALLOCATE(LPATCHNY)
       DEALLOCATE(LPATCHNZ)
       DEALLOCATE(LPATCHX)
       DEALLOCATE(LPATCHY)
       DEALLOCATE(LPATCHZ)
       DEALLOCATE(LPATCHRX)
       DEALLOCATE(LPATCHRY)
       DEALLOCATE(LPATCHRZ)
       DEALLOCATE(LPARE)


*********************
       END IF        !IF (NPATCH(IR-1).GT.0) THEN
*********************

       WRITE(*,*) 'LEVEL=',IR,' patch=',NPATCH(IR)

*********************
       END DO     !IR
*********************

       !////IF ONLY DM///
       IF (N_GAS.EQ.0) THEN
        U1G=0.0
        U11G=0.0
       END IF

       IF (N_GAS.GT.0) THEN
        WRITE(*,*)'MASAG,*UM',MASAG(1),MASAG(1)*9.1717E18
        DEALLOCATE(RXG)
        DEALLOCATE(RYG)
        DEALLOCATE(RZG)
        DEALLOCATE(VXG)
        DEALLOCATE(VYG)
        DEALLOCATE(VZG)
        DEALLOCATE(MASAG)
       END IF

       NPART(0)=N_DM
       NPART(1:NL)=0


***
       MASAP=MASAP*MAXMAP
       U1=U1*MAXMAP
       U11=U11*MAXMAP
**
       DO IR=0, NL
        WRITE(*,*) IR,NPART(IR)
       END DO

       WRITE(*,*)'MASAP=',MINVAL(MASAP(1:N_DM)),
     &                    MAXVAL(MASAP(1:N_DM))
       WRITE(*,*)'Total part.', SUM(NPART(0:NL))
       WRITE(*,*) 'ORIPA1=',MAXVAL(ORIPA1(1:N_DM)),
     &                      MINVAL(ORIPA1(1:N_DM))
       WRITE(*,*) 'ORIPA2=',MAXVAL(ORIPA2(1:N_DM)),
     &                      MINVAL(ORIPA2(1:N_DM))



*      WRITING GRID DATA ON A FILE
       CALL NOMFILE6(ITER,FILE6)
       FILERR='./output_files/'//FILE6
       OPEN(33,FILE=FILERR,STATUS='UNKNOWN')

       WRITE(33,*) NL,N_DM,N_GAS
       DO IR=1,NL
        WRITE(33,*) '------new level-------'
        WRITE(33,*) IR,NPATCH(IR), NPART(IR)

       DO I=LOW1,LOW2
        WRITE(33,*) PATCHNX(I),PATCHNY(I),PATCHNZ(I)
        WRITE(33,*) PATCHX(I),PATCHY(I),PATCHZ(I)
        WRITE(33,*) PATCHRX(I),PATCHRY(I),PATCHRZ(I)
        WRITE(33,*) PARE(I)
       END DO
       END DO
       CLOSE(33)

       RETURN
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
              IF (SOLAP(IX,JY,KZ,I).EQ.1) SOLAP(II,JJ,KK,J)=0
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
