*********************************************************************
       SUBROUTINE READ_MASCLET(VAR,ITER,NX,NY,NZ,NDXYZ,T,ZETA,NL,NPATCH,
     &           PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &           PATCHRX,PATCHRY,PATCHRZ,
     &           U2DM,U3DM,U4DM,MASAP,NPART,RXPA,RYPA,RZPA,ORIPA,N_DM)
*********************************************************************
*      Reads MASCLET data: grids, gas density (clus files) and
*      DM particles information.
*      Must be checked depending on the version/flavour of MASCLET
*      the simulation has been run with
*********************************************************************

       IMPLICIT NONE

       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NX,NY,NZ,ITER,NDXYZ,N_DM
       REAL*4 T,AAA,BBB,CCC,MAP,ZETA

       INTEGER I,J,K,IX,NL,IR,IRR,VAR,N1,N2,N3

*      VARIABLES
       REAL*4 U1(NMAX,NMAY,NMAZ)
       REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
       COMMON /VARIA/ U1,U11  !!,U12,U13,U14

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

       INTEGER NPART(0:NLEVELS)
       REAL*4 U2DM(PARTI_READ)
       REAL*4 U3DM(PARTI_READ)
       REAL*4 U4DM(PARTI_READ)
       REAL*4 MASAP(PARTI_READ)
       REAL*4 RXPA(PARTI_READ)
       REAL*4 RYPA(PARTI_READ)
       REAL*4 RZPA(PARTI_READ)

       INTEGER ORIPA(PARTI_READ)

       REAL*4, ALLOCATABLE::SCR(:,:,:)

       REAL*4 UBAS(0:PARTI_READ)
       INTEGER UBAS2(0:PARTI_READ),CONTA,LOW1,LOW2

       CHARACTER*5 ITER_STRING


*      READING DATA
       WRITE(ITER_STRING, '(I5.5)') ITER !For saving files to disk
       WRITE(*,*) 'Reading iter',ITER

       OPEN (33,FILE='./simu_masclet/grids'//ITER_STRING,
     &       STATUS='UNKNOWN',ACTION='READ')
       OPEN (31,FILE='./simu_masclet/clus'//ITER_STRING,
     &       STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')
       OPEN (32,FILE='./simu_masclet/cldm'//ITER_STRING,
     &       STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')


*      GRID DATA
       READ(33,*) IRR,T,NL,MAP
       READ(33,*) ZETA
       READ(33,*) IR,NDXYZ
       WRITE(*,*) 'IR,NL,NDXYZ,MAP', IR,NL,NDXYZ,MAP

       DO IR=1,NL
       READ(33,*) IRR,NPATCH(IR), NPART(IR)
       WRITE(*,*) 'NPATCH(IR), NPART(IR)',NPATCH(IR), NPART(IR)
       READ(33,*)

       IF (IR.NE.IRR) WRITE(*,*)'Warning: fail in restart'
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
*       DO I=1,NPATCH(IR)
       DO I=LOW1,LOW2
        READ(33,*) PATCHNX(I),PATCHNY(I),PATCHNZ(I)
        READ(33,*) PATCHX(I),PATCHY(I),PATCHZ(I)
        READ(33,*) AAA,BBB,CCC
        PATCHRX(I)=AAA
        PATCHRY(I)=BBB
        PATCHRZ(I)=CCC
        READ(33,*) PARE(I)
       END DO
       END DO
       CLOSE(33)
       NPART(0)=NDXYZ

C       IF (VAR.EQ.1) THEN
*      BARYONIC
       READ(31)
       IR=0
        N1=NX
        N2=NY
        N3=NZ
        READ(31) !(((U1G(I,J,K),I=1,N1),J=1,N2),K=1,N3)
        READ(31) !U2
        READ(31) !U3
        READ(31) !U4
        READ(31) !PRES
        READ(31) !POT
        READ(31) !!OPOT
        READ(31) !!CAUTION with this line!! depends on MASCLET version: T
        READ(31) !!new: metalicity!! depends on MASCLET version!: TRACER
        READ(31) !!ABS(CR0AMR(1:N1,1:N2,1:N3)-1)
        READ(31) !Bx
        READ(31) !By
        READ(31) !Bz

       ALLOCATE(SCR(NAMRX,NAMRY,NAMRZ))
       SCR=0.0
       DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        READ(31) !(((SCR(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
           !U11G(1:N1,1:N2,1:N3,I)=SCR(1:N1,1:N2,1:N3)
        READ(31) !!(((SCR(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
           !!U12(1:N1,1:N2,1:N3,I)=SCR(1:N1,1:N2,1:N3)
        READ(31) !!(((SCR(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
           !!U13(1:N1,1:N2,1:N3,I)=SCR(1:N1,1:N2,1:N3)
        READ(31) !!(((SCR(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
           !!U14(1:N1,1:N2,1:N3,I)=SCR(1:N1,1:N2,1:N3)
        READ(31) !PRES21
        READ(31) !POT1
        READ(31) !OPOT
        READ(31) !!CAUTION with this line!! depends on MASCLET version: T
        READ(31) !!new: metalicity!! depends on MASCLET version! TRACER
        READ(31) !!ABS(CR0AMR11(:,:,:,I)-1) (celdas=0 se eliminan del nivel por estar refinadas)
        READ(31) !!ABS(SOLAPST(:,:,:,I)-1): se guarda solapst para saber que celdas estan solapadas
        READ(31) !Bx
        READ(31) !By
        READ(31) !Bz


       END DO
       END DO
       DEALLOCATE(SCR)

C      END IF
      CLOSE(31)

***       IF (VAR.EQ.2) THEN
**     DARK MATTER
       READ(32)
       IR=0
        N1=NX
        N2=NY
        N3=NZ

        READ(32) (((U1(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
        READ(32) (RXPA(I),I=1,NDXYZ)
        READ(32) (RYPA(I),I=1,NDXYZ)
        READ(32) (RZPA(I),I=1,NDXYZ)
C        WRITE(*,*)'HOLA2', MAXVAL(RXPA(1:NDXYZ)), MAXVAL(RYPA(1:NDXYZ))

        READ(32) (U2DM(I),I=1,NDXYZ)
        READ(32) (U3DM(I),I=1,NDXYZ)
        READ(32) (U4DM(I),I=1,NDXYZ)
        READ(32) (ORIPA(I),I=1,NDXYZ)    !particle ID
        CONTA=NDXYZ
        MASAP(1:NDXYZ)=MAP
        WRITE(*,*) 'ORIPA=',MAXVAL(ORIPA(1:NDXYZ)),
     &                       MINVAL(ORIPA(1:NDXYZ))
        WRITE(*,*) 'NPART(0)=',IR, NPART(IR),CONTA


       ALLOCATE(SCR(NAMRX,NAMRY,NAMRZ))
       SCR=0.0
       DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)

        READ(32) (((SCR(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
           U11(1:N1,1:N2,1:N3,I)=SCR(1:N1,1:N2,1:N3)
       END DO

        UBAS=0.0
        UBAS2=0
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        RXPA(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        RYPA(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        RZPA(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        U2DM(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        U3DM(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        U4DM(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        MASAP(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))

        READ(32) (UBAS2(IX),IX=1,NPART(IR))
        IF (NPART(IR).GT.0)
     &   ORIPA(CONTA+1:CONTA+NPART(IR))=UBAS2(1:NPART(IR))

        IF (NPART(IR).GT.0) THEN
        WRITE(*,*) 'ORIPA=',MAXVAL(ORIPA(CONTA+1:CONTA+NPART(IR))),
     &                       MINVAL(ORIPA(CONTA+1:CONTA+NPART(IR)))
        END IF

        CONTA=CONTA+NPART(IR)
        WRITE(*,*) 'NPART(IR)=',IR,NPART(IR),CONTA


       END DO

       DEALLOCATE(SCR)

       CLOSE(32)
***       END IF

       WRITE(*,*) 'TOTAL PARTICLES IN ITER=',CONTA
       IF (CONTA.NE.N_DM) THEN
        WRITE(*,*) 'WARNING: CONTA != N_DM',CONTA,N_DM
        STOP
       END IF

       RETURN
       END

*********************************************************************
       SUBROUTINE READ_PARTICLES_MASCLET(ITER,NX,NY,NZ,T,ZETA,
     &             U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,ORIPA,
     &             N_DM,VAR,N_ST)
*********************************************************************
*      Reads MASCLET data: grids and DM particles information.
*      Must be checked depending on the version/flavour of MASCLET
*      the simulation has been run with
*********************************************************************

       IMPLICIT NONE

       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NX,NY,NZ,ITER,NDXYZ
       REAL*4 T,AAA,BBB,CCC,MAP,ZETA
       INTEGER VAR !(=1: only DM; =2: DM+stars)

       REAL BAS
       INTEGER I,J,K,IX,NL,IR,IRR,N1,N2,N3,N_DM,N_ST,NBAS,ARE_BH,NST0

       INTEGER,ALLOCATABLE::NPATCH(:)

       INTEGER,ALLOCATABLE::NPART(:),NPARTST(:),NPARTBH(:)
       REAL*4 U2DM(PARTI_READ)
       REAL*4 U3DM(PARTI_READ)
       REAL*4 U4DM(PARTI_READ)
       REAL*4 MASAP(PARTI_READ)
       REAL*4 RXPA(PARTI_READ)
       REAL*4 RYPA(PARTI_READ)
       REAL*4 RZPA(PARTI_READ)

       INTEGER ORIPA(PARTI_READ)

       REAL*4 UBAS(0:PARTI_READ)
       INTEGER UBAS2(0:PARTI_READ),CONTA,LOW1,LOW2

       CHARACTER*5 ITER_STRING

       ARE_BH=1 ! Depends on MASCLET version (are there BHs??)

*      READING DATA
       WRITE(ITER_STRING, '(I5.5)') ITER !For saving files to disk
       WRITE(*,*) 'Reading iter',ITER

*      GRID DATA
       OPEN (33,FILE='./simu_masclet/grids'//ITER_STRING,
     &       STATUS='UNKNOWN',ACTION='READ')
       READ(33,*) IRR,T,NL,MAP
       READ(33,*) ZETA
       READ(33,*) IR,NDXYZ,NST0
       !WRITE(*,*) 'IR,NL,NDXYZ,MAP', IR,NL,NDXYZ,MAP

       ALLOCATE(NPATCH(0:NL),NPART(0:NL),NPARTST(0:NL),NPARTBH(0:NL))
       NPATCH=0
       NPART=0
       NPARTST=0
       NPARTBH=0

       DO IR=1,NL
       IF (ARE_BH.EQ.0) THEN
        READ(33,*) IRR,NPATCH(IR),NPART(IR),NPARTST(IR)!,NPARTBH(IR)
       ELSE
        READ(33,*) IRR,NPATCH(IR),NPART(IR),NPARTST(IR),NPARTBH(IR)
       END IF
       !WRITE(*,*) 'NPATCH(IR), NPART(IR)',NPATCH(IR), NPART(IR)
       READ(33,*)

       IF (IR.NE.IRR) WRITE(*,*)'Warning: fail in restart'
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
*       DO I=1,NPATCH(IR)
       DO I=LOW1,LOW2
        READ(33,*) !PATCHNX(I),PATCHNY(I),PATCHNZ(I)
        READ(33,*) !PATCHX(I),PATCHY(I),PATCHZ(I)
        READ(33,*) !PATCHRX(I),PATCHRY(I),PATCHRZ(I)
        READ(33,*) !PARE(I)
       END DO
       END DO
       NPART(0)=NDXYZ
       NPARTST(0)=NST0

       CLOSE(33)

**     DARK MATTER
       OPEN (32,FILE='./simu_masclet/cldm'//ITER_STRING,
     &       STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')
       READ(32)
       !IR=0
       READ(32) !(((U1(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       READ(32) (RXPA(I),I=1,NDXYZ)
       READ(32) (RYPA(I),I=1,NDXYZ)
       READ(32) (RZPA(I),I=1,NDXYZ)
       READ(32) (U2DM(I),I=1,NDXYZ)
       READ(32) (U3DM(I),I=1,NDXYZ)
       READ(32) (U4DM(I),I=1,NDXYZ)
       READ(32) (ORIPA(I),I=1,NDXYZ)
       CONTA=NDXYZ
       MASAP(1:NDXYZ)=MAP

C       WRITE(*,*) 'ORIPA=',MAXVAL(ORIPA(1:NDXYZ)),
C     &                      MINVAL(ORIPA(1:NDXYZ))
       WRITE(*,*) 'NPART(IR)=',0, NPART(0),CONTA

       DO IR=1,NL
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))

        DO I=LOW1,LOW2
         READ(32) !(((SCR(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
        END DO

        UBAS=0.0
        UBAS2=0
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        RXPA(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        RYPA(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        RZPA(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        U2DM(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        U3DM(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        U4DM(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        MASAP(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))

        READ(32) (UBAS2(IX),IX=1,NPART(IR))
        ORIPA(CONTA+1:CONTA+NPART(IR))=UBAS2(1:NPART(IR))

        CONTA=CONTA+NPART(IR)
        WRITE(*,*) 'NPART(IR)=',IR,NPART(IR),CONTA
       END DO

       CLOSE(32)

       WRITE(*,*) 'TOTAL DM PARTICLES IN ITER=',CONTA
       N_DM=SUM(NPART(0:NL))
       IF (CONTA.NE.N_DM) THEN
        WRITE(*,*) 'WARNING: N_DM rewritten:',N_DM,'-->',CONTA
        N_DM=CONTA
       END IF

*      Fix ORIPAs: particles of the heavier species get negative
       BAS=MAXVAL(MASAP(1:N_DM))
       IF (BAS.GT.4.0*MINVAL(MASAP(1:N_DM))) THEN
        BAS=0.9*BAS
!$OMP PARALLEL DO SHARED(N_DM,BAS,ORIPA,MASAP),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
        DO I=1,N_DM
         IF (MASAP(I).GT.BAS) ORIPA(I)=-ABS(ORIPA(I))
        END DO
       END IF
*      END Fix ORIPAs

       IF (VAR.EQ.2) THEN
        N_ST=SUM(NPARTST(0:NL))+SUM(NPARTBH(0:NL))

        IF (N_DM+N_ST.GT.PARTI_READ) THEN
         WRITE(*,*) 'WARNING: bad dimensioning of PARTI_READ',
     &               N_DM+N_ST,'>',PARTI_READ
         STOP
        END IF

        OPEN (34,FILE='./simu_masclet/clst'//ITER_STRING,
     &        STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')

        READ(34) !ITER,T4,ZETA
        !IR=0
        NBAS=NPARTST(0)+NPARTBH(0)
        READ(34) !(((U1ST(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
        READ(34) (UBAS(I),I=1,NBAS)
        RXPA(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
        READ(34) (UBAS(I),I=1,NBAS)
        RYPA(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
        READ(34) (UBAS(I),I=1,NBAS)
        RZPA(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
        READ(34) (UBAS(I),I=1,NBAS)
        U2DM(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
        READ(34) (UBAS(I),I=1,NBAS)
        U3DM(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
        READ(34) (UBAS(I),I=1,NBAS)
        U4DM(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
        READ(34) (UBAS(I),I=1,NBAS)
        MASAP(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
        READ(34) !(UBAS(I),I=1,NDXYZ)
        READ(34) !(UBAS(I),I=1,NDXYZ)
        CONTA=CONTA+NBAS

        WRITE(*,*) 'NST(IR)=',0,NBAS,CONTA

        DO IR=1,NL
         LOW1=SUM(NPATCH(0:IR-1))+1
         LOW2=SUM(NPATCH(0:IR))
         DO I=LOW1,LOW2
          READ(34) !(((SCR(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
         END DO

         NBAS=NPARTST(IR)
         READ(34) (UBAS(IX),IX=1,NBAS)
         RXPA(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
         READ(34) (UBAS(IX),IX=1,NBAS)
         RYPA(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
         READ(34) (UBAS(IX),IX=1,NBAS)
         RZPA(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
         READ(34) (UBAS(IX),IX=1,NBAS)
         U2DM(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
         READ(34) (UBAS(IX),IX=1,NBAS)
         U3DM(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
         READ(34) (UBAS(IX),IX=1,NBAS)
         U4DM(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
         READ(34) (UBAS(IX),IX=1,NBAS)
         MASAP(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
         READ(34)
         READ(34)
         READ(34) ! ORIPAST
         CONTA=CONTA+NBAS

         NBAS=NPARTBH(IR)
         IF (ARE_BH.EQ.1) THEN
          READ(34) (UBAS(IX),IX=1,NBAS)
          RXPA(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
          READ(34) (UBAS(IX),IX=1,NBAS)
          RYPA(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
          READ(34) (UBAS(IX),IX=1,NBAS)
          RZPA(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
          READ(34) (UBAS(IX),IX=1,NBAS)
          U2DM(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
          READ(34) (UBAS(IX),IX=1,NBAS)
          U3DM(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
          READ(34) (UBAS(IX),IX=1,NBAS)
          U4DM(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
          READ(34) (UBAS(IX),IX=1,NBAS)
          MASAP(CONTA+1:CONTA+NBAS)=UBAS(1:NBAS)
          READ(34)
          READ(34)
          CONTA=CONTA+NBAS
         END IF

         WRITE(*,*) 'NST(IR)=',IR,NPARTST(IR)+NPARTBH(IR),CONTA
        END DO

        CLOSE(34)
       ELSE
        N_ST=0
       END IF !(VAR.EQ.2)

       DEALLOCATE(NPATCH,NPART,NPARTST,NPARTBH)

       RETURN
       END


*********************************************************************
      SUBROUTINE READ_PARTICLES_GENERAL(ITER,NX,NY,NZ,T,ZETA,
     &             U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,ORIPA,
     &             N_DM,VAR,N_ST,UV,UM,HUBBLE_LITTLEH)
*********************************************************************
*      Reads data from generic input format
*********************************************************************

      IMPLICIT NONE

      INCLUDE 'input_files/asohf_parameters.dat'

      INTEGER NX,NY,NZ,ITER,NDXYZ
      REAL*4 T,AAA,BBB,CCC,MAP,ZETA,UV,UM,HUBBLE_LITTLEH
      INTEGER VAR !(=1: only DM; =2: DM+stars)

      INTEGER I,J,K,IX,NL,IR,IRR,N_DM,N_ST,NBAS,NST0

      REAL*4 U2DM(PARTI_READ)
      REAL*4 U3DM(PARTI_READ)
      REAL*4 U4DM(PARTI_READ)
      REAL*4 MASAP(PARTI_READ)
      REAL*4 RXPA(PARTI_READ)
      REAL*4 RYPA(PARTI_READ)
      REAL*4 RZPA(PARTI_READ)

      INTEGER ORIPA(PARTI_READ)

      REAL*4 UBAS(0:PARTI_READ)
      INTEGER UBASINT(0:PARTI_READ)

      REAL CIO_MASS,CIO_SPEED,CIO_LENGTH,CIO_ALPHA,CIO_XC,CIO_YC,CIO_ZC
      COMMON /CONV_IO/ CIO_MASS,CIO_SPEED,CIO_LENGTH,CIO_ALPHA,CIO_XC,
     &                 CIO_YC,CIO_ZC
      REAL FACT_MASS,FACT_SPEED,FACT_LENGTH

      CHARACTER*5 ITER_STRING

*     READING DATA
      WRITE(ITER_STRING, '(I5.5)') ITER !For saving files to disk
      WRITE(*,*) 'Reading iter',ITER

**    DARK MATTER
      OPEN (32,FILE='./simulation/particles'//ITER_STRING,
     &         STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')
      READ(32) ZETA
      READ(32) N_DM

      READ(32) (UBAS(I),I=1,N_DM)
      RXPA(1:N_DM)=UBAS(1:N_DM)
      READ(32) (UBAS(I),I=1,N_DM)
      RYPA(1:N_DM)=UBAS(1:N_DM)
      READ(32) (UBAS(I),I=1,N_DM)
      RZPA(1:N_DM)=UBAS(1:N_DM)
      READ(32) (UBAS(I),I=1,N_DM)
      U2DM(1:N_DM)=UBAS(1:N_DM)
      READ(32) (UBAS(I),I=1,N_DM)
      U3DM(1:N_DM)=UBAS(1:N_DM)
      READ(32) (UBAS(I),I=1,N_DM)
      U4DM(1:N_DM)=UBAS(1:N_DM)
      READ(32) (UBAS(I),I=1,N_DM)
      MASAP(1:N_DM)=UBAS(1:N_DM)
      READ(32) (UBASINT(I),I=1,N_DM)
      ORIPA(1:N_DM)=UBASINT(1:N_DM)

      WRITE(*,*)
      WRITE(*,*) 'INPUT. DM x positions (min,max):',
     &            MINVAL(RXPA(1:N_DM)),MAXVAL(RXPA(1:N_DM))
      WRITE(*,*) 'INPUT. DM y positions (min,max):',
     &            MINVAL(RYPA(1:N_DM)),MAXVAL(RYPA(1:N_DM))
      WRITE(*,*) 'INPUT. DM z positions (min,max):',
     &            MINVAL(RZPA(1:N_DM)),MAXVAL(RZPA(1:N_DM))
      WRITE(*,*) 'INPUT. DM x velocities (min,max):',
     &            MINVAL(U2DM(1:N_DM)),MAXVAL(U2DM(1:N_DM))
      WRITE(*,*) 'INPUT. DM y velocities (min,max):',
     &            MINVAL(U3DM(1:N_DM)),MAXVAL(U3DM(1:N_DM))
      WRITE(*,*) 'INPUT. DM z velocities (min,max):',
     &            MINVAL(U4DM(1:N_DM)),MAXVAL(U4DM(1:N_DM))
      WRITE(*,*) 'INPUT. DM masses (min,max):',
     &            MINVAL(MASAP(1:N_DM)),MAXVAL(MASAP(1:N_DM))
      WRITE(*,*) 'INPUT. DM unique IDs (min,max):',
     &            MINVAL(ORIPA(1:N_DM)),MAXVAL(ORIPA(1:N_DM))
      WRITE(*,*)
      WRITE(*,*) 'TOTAL DM PARTICLES IN ITER=',N_DM


      IF (VAR.EQ.2) THEN
       READ(32) N_ST
       IF (N_DM+N_ST.GT.PARTI_READ) THEN
        WRITE(*,*) 'WARNING: bad dimensioning of PARTI_READ',
     &              N_DM+N_ST,'>',PARTI_READ
        STOP
       END IF

       READ(32) (UBAS(I),I=1,N_ST)
       RXPA(N_DM+1:N_DM+N_ST)=UBAS(1:N_ST)
       READ(32) (UBAS(I),I=1,N_ST)
       RYPA(N_DM+1:N_DM+N_ST)=UBAS(1:N_ST)
       READ(32) (UBAS(I),I=1,N_ST)
       RZPA(N_DM+1:N_DM+N_ST)=UBAS(1:N_ST)
       READ(32) (UBAS(I),I=1,N_ST)
       U2DM(N_DM+1:N_DM+N_ST)=UBAS(1:N_ST)
       READ(32) (UBAS(I),I=1,N_ST)
       U3DM(N_DM+1:N_DM+N_ST)=UBAS(1:N_ST)
       READ(32) (UBAS(I),I=1,N_ST)
       U4DM(N_DM+1:N_DM+N_ST)=UBAS(1:N_ST)
       READ(32) (UBAS(I),I=1,N_ST)
       MASAP(N_DM+1:N_DM+N_ST)=UBAS(1:N_ST)
       READ(32) (UBASINT(I),I=1,N_ST)
       ORIPA(N_DM+1:N_DM+N_ST)=UBASINT(1:N_ST)

       WRITE(*,*)
       WRITE(*,*) 'INPUT. ST x positions (min,max):',
     &     MINVAL(RXPA(N_DM+1:N_DM+N_ST)),MAXVAL(RXPA(N_DM+1:N_DM+N_ST))
       WRITE(*,*) 'INPUT. ST y positions (min,max):',
     &     MINVAL(RYPA(N_DM+1:N_DM+N_ST)),MAXVAL(RYPA(N_DM+1:N_DM+N_ST))
       WRITE(*,*) 'INPUT. ST z positions (min,max):',
     &     MINVAL(RZPA(N_DM+1:N_DM+N_ST)),MAXVAL(RZPA(N_DM+1:N_DM+N_ST))
       WRITE(*,*) 'INPUT. ST x velocities (min,max):',
     &     MINVAL(U2DM(N_DM+1:N_DM+N_ST)),MAXVAL(U2DM(N_DM+1:N_DM+N_ST))
       WRITE(*,*) 'INPUT. ST y velocities (min,max):',
     &     MINVAL(U3DM(N_DM+1:N_DM+N_ST)),MAXVAL(U3DM(N_DM+1:N_DM+N_ST))
       WRITE(*,*) 'INPUT. ST z velocities (min,max):',
     &     MINVAL(U4DM(N_DM+1:N_DM+N_ST)),MAXVAL(U4DM(N_DM+1:N_DM+N_ST))
       WRITE(*,*) 'INPUT. ST masses (min,max):',
     &   MINVAL(MASAP(N_DM+1:N_DM+N_ST)),MAXVAL(MASAP(N_DM+1:N_DM+N_ST))
       WRITE(*,*) 'INPUT. ST unique IDs (min,max):',
     &   MINVAL(ORIPA(N_DM+1:N_DM+N_ST)),MAXVAL(ORIPA(N_DM+1:N_DM+N_ST))
       WRITE(*,*)
       WRITE(*,*) 'TOTAL STELLAR PARTICLES IN ITER=',N_ST
      ELSE
       N_ST=0
      END IF

      CLOSE(32)

      FACT_MASS=CIO_MASS/UM
      FACT_SPEED=(CIO_SPEED/UV)*(1+ZETA)**(CIO_ALPHA-1.0)
      FACT_LENGTH=CIO_LENGTH

!$OMP PARALLEL DO SHARED(N_DM,N_ST,RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,
!$OMP+                   FACT_LENGTH,FACT_SPEED,FACT_MASS,CIO_XC,CIO_YC,
!$OMP+                   CIO_ZC),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
      DO I=1,N_DM+N_ST
       RXPA(I)=(RXPA(I)-CIO_XC)*FACT_LENGTH
       RYPA(I)=(RYPA(I)-CIO_YC)*FACT_LENGTH
       RZPA(I)=(RZPA(I)-CIO_ZC)*FACT_LENGTH

       U2DM(I)=U2DM(I)*FACT_SPEED
       U3DM(I)=U3DM(I)*FACT_SPEED
       U4DM(I)=U4DM(I)*FACT_SPEED

       MASAP(I)=MASAP(I)*FACT_MASS
      END DO

      WRITE(*,*)
      WRITE(*,*) 'After unit conversion...'
      WRITE(*,*) 'x positions (min,max), in Mpc:',
     &     MINVAL(RXPA(1:N_DM+N_ST)),MAXVAL(RXPA(1:N_DM+N_ST))
      WRITE(*,*) 'y positions (min,max), in Mpc:',
     &     MINVAL(RYPA(1:N_DM+N_ST)),MAXVAL(RYPA(1:N_DM+N_ST))
      WRITE(*,*) 'z positions (min,max), in Mpc:',
     &     MINVAL(RZPA(1:N_DM+N_ST)),MAXVAL(RZPA(1:N_DM+N_ST))
      WRITE(*,*) 'x velocities (min,max), in c:',
     &     MINVAL(U2DM(1:N_DM+N_ST)),MAXVAL(U2DM(1:N_DM+N_ST))
      WRITE(*,*) 'y velocities (min,max), in c:',
     &     MINVAL(U3DM(1:N_DM+N_ST)),MAXVAL(U3DM(1:N_DM+N_ST))
      WRITE(*,*) 'z velocities (min,max), in c:',
     &     MINVAL(U4DM(1:N_DM+N_ST)),MAXVAL(U4DM(1:N_DM+N_ST))
      WRITE(*,*) 'Masses (min,max), in internal units:',
     &   MINVAL(MASAP(1:N_DM+N_ST)),MAXVAL(MASAP(1:N_DM+N_ST))
      WRITE(*,*) 'Unique IDs (min,max):',
     &   MINVAL(ORIPA(1:N_DM+N_ST)),MAXVAL(ORIPA(1:N_DM+N_ST))
      WRITE(*,*)

      RETURN
      END


*********************************************************************
      SUBROUTINE READ_GADGET_UNFORMATTED(ITER,NX,NY,NZ,T,ZETA,
     &             U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,ORIPA,
     &             N_DM,VAR,N_ST,UV,UM,HUBBLE_LITTLEH)
*********************************************************************
*      Reads data from generic input format
*********************************************************************

      IMPLICIT NONE

      INCLUDE 'input_files/asohf_parameters.dat'

      INTEGER NX,NY,NZ,ITER,NDXYZ
      REAL*4 T,AAA,BBB,CCC,MAP,ZETA,UV,UM,HUBBLE_LITTLEH
      INTEGER VAR !(=1: only DM; =2: DM+stars)

      INTEGER N_DM,N_ST,NBAS,NST0

      REAL*4 U2DM(PARTI_READ)
      REAL*4 U3DM(PARTI_READ)
      REAL*4 U4DM(PARTI_READ)
      REAL*4 MASAP(PARTI_READ)
      REAL*4 RXPA(PARTI_READ)
      REAL*4 RYPA(PARTI_READ)
      REAL*4 RZPA(PARTI_READ)
      INTEGER ORIPA(PARTI_READ)

      REAL CIO_MASS,CIO_SPEED,CIO_LENGTH,CIO_ALPHA,CIO_XC,CIO_YC,CIO_ZC
      COMMON /CONV_IO/ CIO_MASS,CIO_SPEED,CIO_LENGTH,CIO_ALPHA,CIO_XC,
     &                 CIO_YC,CIO_ZC
      REAL FACT_MASS,FACT_SPEED,FACT_LENGTH

*     Local variables
      INTEGER NTOT,IGAS0,IGAS1,IDM0,IDM1,IST0,IST1,IBH0,IBH1,I,J

*     IO variables
      CHARACTER*4 BLOCKLABEL
      INTEGER*4 BLOCKSIZE,NPP(6)
      REAL*8 MASS_ARR(6)
      REAL*8 TIME8,ZETA8
      REAL*8 CACA(5)
      REAL*8 BOXSIZE8,OMEGA_M8,OMEGA_LAMBDA8,HUBBLE_PARAM8

      REAL*4,ALLOCATABLE::SCR41(:)
      REAL*4,ALLOCATABLE::SCR42(:,:)
      INTEGER*4,ALLOCATABLE::SCRINT1(:)

      CHARACTER*3 ITER_STRING

***** FIX UNIT CONVERSIONS FOR GADGET *******
      !CIO_MASS=1.E10/HUBBLE_LITTLEH
      !CIO_SPEED=299792.458
      !CIO_LENGTH=1.0/HUBBLE_LITTLEH
      !CIO_ALPHA=0.5
*********************************************

*     READING DATA
      WRITE(ITER_STRING, '(I3.3)') ITER !For saving files to disk
      WRITE(*,*) ' > Reading iter',ITER_STRING

      OPEN(11, FILE='./simulation/snap_'//ITER_STRING,
     &     STATUS='UNKNOWN',ACTION='READ', FORM='UNFORMATTED')

*      Read the header ************************************************
       READ(11) BLOCKLABEL,BLOCKSIZE
       !WRITE(*,*) 'Found block ', BLOCKLABEL, ' with length', BLOCKSIZE
       READ(11) NPP,MASS_ARR,TIME8,ZETA8,CACA,BOXSIZE8,OMEGA_M8,
     &          OMEGA_LAMBDA8,HUBBLE_PARAM8

       !WRITE(*,*) NPP
       !WRITE(*,*) MASS_ARR
       !WRITE(*,*) 'Redshift:', ZETA8
       !WRITE(*,*) 'Box size (ckpc/h):', BOXSIZE8
       !WRITE(*,*) 'Om, Olambda, h:', OMEGA_M8, OMEGA_LAMBDA8,
       !&            HUBBLE_PARAM8

       NTOT=SUM(NPP)

       IGAS0=1
       IGAS1=NPP(1)
       IDM0=IGAS1+1
       IDM1=SUM(NPP(1:4))
       IST0=IDM1+1
       IST1=SUM(NPP(1:5))
       IBH0=IST1+1
       IBH1=SUM(NPP(1:6))
       !WRITE(*,*) 'Gas indices',IGAS0,IGAS1,IGAS1-IGAS0+1
       !WRITE(*,*) 'DM indices ',IDM0,IDM1,IDM1-IDM0+1
       !WRITE(*,*) 'ST indices ',IST0,IST1,IST1-IST0+1
       !WRITE(*,*) 'BH indices ',IBH0,IBH1,IBH1-IBH0+1

       N_DM=IDM1-IDM0+1
       IF (VAR.EQ.2) THEN
        N_ST=IST1-IST0+1
       ELSE
        N_ST=0
       END IF

*      Read the particle positions ************************************
       READ(11) BLOCKLABEL,BLOCKSIZE
       !WRITE(*,*) 'Found block ', BLOCKLABEL, ' with length', BLOCKSIZE
       ALLOCATE(SCR42(3,NTOT))
       READ(11) ((SCR42(J,I),J=1,3),I=1,NTOT) ! all particles
       !WRITE(*,*) '-X-', MINVAL(SCR42(1,:)), MAXVAL(SCR42(1,:))
       !WRITE(*,*) '-Y-', MINVAL(SCR42(2,:)), MAXVAL(SCR42(2,:))
       !WRITE(*,*) '-Z-', MINVAL(SCR42(3,:)), MAXVAL(SCR42(3,:))
       RXPA(1:N_DM)=SCR42(1,IDM0:IDM1)
       RYPA(1:N_DM)=SCR42(2,IDM0:IDM1)
       RZPA(1:N_DM)=SCR42(3,IDM0:IDM1)
       IF (VAR.EQ.2) THEN
        RXPA(N_DM+1:N_DM+N_ST)=SCR42(1,IST0:IST1)
        RYPA(N_DM+1:N_DM+N_ST)=SCR42(2,IST0:IST1)
        RZPA(N_DM+1:N_DM+N_ST)=SCR42(3,IST0:IST1)
       END IF

*      Read the particle velocities************************************
       READ(11) BLOCKLABEL,BLOCKSIZE
       !WRITE(*,*) 'Found block ', BLOCKLABEL, ' with length', BLOCKSIZE
       READ(11) ((SCR42(J,I),J=1,3),I=1,NTOT) ! all particles
       !WRITE(*,*) '-VX-', MINVAL(SCR42(1,:)), MAXVAL(SCR42(1,:))
       !WRITE(*,*) '-VY-', MINVAL(SCR42(2,:)), MAXVAL(SCR42(2,:))
       !WRITE(*,*) '-VZ-', MINVAL(SCR42(3,:)), MAXVAL(SCR42(3,:))
       U2DM(1:N_DM)=SCR42(1,IDM0:IDM1)
       U3DM(1:N_DM)=SCR42(2,IDM0:IDM1)
       U4DM(1:N_DM)=SCR42(3,IDM0:IDM1)
       IF (VAR.EQ.2) THEN
        U2DM(N_DM+1:N_DM+N_ST)=SCR42(1,IST0:IST1)
        U3DM(N_DM+1:N_DM+N_ST)=SCR42(2,IST0:IST1)
        U4DM(N_DM+1:N_DM+N_ST)=SCR42(3,IST0:IST1)
       END IF
       DEALLOCATE(SCR42)

*      Read the particle ids****************************************
       READ(11) BLOCKLABEL,BLOCKSIZE
       !WRITE(*,*) 'Found block ', BLOCKLABEL, ' with length', BLOCKSIZE
       ALLOCATE(SCRINT1(NTOT))
       READ(11) (SCRINT1(I),I=1,NTOT) ! all particles
       !WRITE(*,*) '-ID-', MINVAL(SCRINT1(:)), MAXVAL(SCRINT1(:))
       ORIPA(1:N_DM)=SCRINT1(IDM0:IDM1)
       IF (VAR.EQ.2) THEN
        ORIPA(N_DM+1:N_DM+N_ST)=SCRINT1(IST0:IST1)
       END IF
       DEALLOCATE(SCRINT1)

*      Read the particle masses****************************************
       READ(11) BLOCKLABEL,BLOCKSIZE
       !WRITE(*,*) 'Found block ', BLOCKLABEL, ' with length', BLOCKSIZE
       ALLOCATE(SCR41(NTOT))
       READ(11) (SCR41(I),I=1,NTOT) ! all particles
       !WRITE(*,*) '-M-', MINVAL(SCR41(:)), MAXVAL(SCR41(:))
       MASAP(1:N_DM)=SCR41(IDM0:IDM1)
       IF (VAR.EQ.2) THEN
        MASAP(N_DM+1:N_DM+N_ST)=SCR41(IST0:IST1)
       END IF
       DEALLOCATE(SCR41)

      CLOSE(11)

      IF (N_DM+N_ST.GT.PARTI_READ) THEN
       WRITE(*,*) 'WARNING: bad dimensioning of PARTI_READ',
     &            N_DM+N_ST,'>',PARTI_READ
       STOP
      END IF
      WRITE(*,*)
      WRITE(*,*) 'INPUT. DM x positions (min,max):',
     &            MINVAL(RXPA(1:N_DM)),MAXVAL(RXPA(1:N_DM))
      WRITE(*,*) 'INPUT. DM y positions (min,max):',
     &            MINVAL(RYPA(1:N_DM)),MAXVAL(RYPA(1:N_DM))
      WRITE(*,*) 'INPUT. DM z positions (min,max):',
     &            MINVAL(RZPA(1:N_DM)),MAXVAL(RZPA(1:N_DM))
      WRITE(*,*) 'INPUT. DM x velocities (min,max):',
     &            MINVAL(U2DM(1:N_DM)),MAXVAL(U2DM(1:N_DM))
      WRITE(*,*) 'INPUT. DM y velocities (min,max):',
     &            MINVAL(U3DM(1:N_DM)),MAXVAL(U3DM(1:N_DM))
      WRITE(*,*) 'INPUT. DM z velocities (min,max):',
     &            MINVAL(U4DM(1:N_DM)),MAXVAL(U4DM(1:N_DM))
      WRITE(*,*) 'INPUT. DM masses (min,max):',
     &            MINVAL(MASAP(1:N_DM)),MAXVAL(MASAP(1:N_DM))
      WRITE(*,*) 'INPUT. DM unique IDs (min,max):',
     &            MINVAL(ORIPA(1:N_DM)),MAXVAL(ORIPA(1:N_DM))
      WRITE(*,*)
      WRITE(*,*) 'TOTAL DM PARTICLES IN ITER=',N_DM


      IF (VAR.EQ.2) THEN
       WRITE(*,*)
       WRITE(*,*) 'INPUT. ST x positions (min,max):',
     &     MINVAL(RXPA(N_DM+1:N_DM+N_ST)),MAXVAL(RXPA(N_DM+1:N_DM+N_ST))
       WRITE(*,*) 'INPUT. ST y positions (min,max):',
     &     MINVAL(RYPA(N_DM+1:N_DM+N_ST)),MAXVAL(RYPA(N_DM+1:N_DM+N_ST))
       WRITE(*,*) 'INPUT. ST z positions (min,max):',
     &     MINVAL(RZPA(N_DM+1:N_DM+N_ST)),MAXVAL(RZPA(N_DM+1:N_DM+N_ST))
       WRITE(*,*) 'INPUT. ST x velocities (min,max):',
     &     MINVAL(U2DM(N_DM+1:N_DM+N_ST)),MAXVAL(U2DM(N_DM+1:N_DM+N_ST))
       WRITE(*,*) 'INPUT. ST y velocities (min,max):',
     &     MINVAL(U3DM(N_DM+1:N_DM+N_ST)),MAXVAL(U3DM(N_DM+1:N_DM+N_ST))
       WRITE(*,*) 'INPUT. ST z velocities (min,max):',
     &     MINVAL(U4DM(N_DM+1:N_DM+N_ST)),MAXVAL(U4DM(N_DM+1:N_DM+N_ST))
       WRITE(*,*) 'INPUT. ST masses (min,max):',
     &   MINVAL(MASAP(N_DM+1:N_DM+N_ST)),MAXVAL(MASAP(N_DM+1:N_DM+N_ST))
       WRITE(*,*) 'INPUT. ST unique IDs (min,max):',
     &   MINVAL(ORIPA(N_DM+1:N_DM+N_ST)),MAXVAL(ORIPA(N_DM+1:N_DM+N_ST))
       WRITE(*,*)
       WRITE(*,*) 'TOTAL STELLAR PARTICLES IN ITER=',N_ST
      END IF

      FACT_MASS=CIO_MASS/UM
      FACT_SPEED=(CIO_SPEED/UV)*(1+ZETA)**(CIO_ALPHA-1.0)
      FACT_LENGTH=CIO_LENGTH

!$OMP PARALLEL DO SHARED(N_DM,N_ST,RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,
!$OMP+                   FACT_LENGTH,FACT_SPEED,FACT_MASS,CIO_XC,CIO_YC,
!$OMP+                   CIO_ZC),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
      DO I=1,N_DM+N_ST
       RXPA(I)=(RXPA(I)-CIO_XC)*FACT_LENGTH
       RYPA(I)=(RYPA(I)-CIO_YC)*FACT_LENGTH
       RZPA(I)=(RZPA(I)-CIO_ZC)*FACT_LENGTH

       U2DM(I)=U2DM(I)*FACT_SPEED
       U3DM(I)=U3DM(I)*FACT_SPEED
       U4DM(I)=U4DM(I)*FACT_SPEED

       MASAP(I)=MASAP(I)*FACT_MASS
      END DO

      WRITE(*,*)
      WRITE(*,*) 'After unit conversion...'
      WRITE(*,*) 'x positions (min,max), in Mpc:',
     &     MINVAL(RXPA(1:N_DM+N_ST)),MAXVAL(RXPA(1:N_DM+N_ST))
      WRITE(*,*) 'y positions (min,max), in Mpc:',
     &     MINVAL(RYPA(1:N_DM+N_ST)),MAXVAL(RYPA(1:N_DM+N_ST))
      WRITE(*,*) 'z positions (min,max), in Mpc:',
     &     MINVAL(RZPA(1:N_DM+N_ST)),MAXVAL(RZPA(1:N_DM+N_ST))
      WRITE(*,*) 'x velocities (min,max), in c:',
     &     MINVAL(U2DM(1:N_DM+N_ST)),MAXVAL(U2DM(1:N_DM+N_ST))
      WRITE(*,*) 'y velocities (min,max), in c:',
     &     MINVAL(U3DM(1:N_DM+N_ST)),MAXVAL(U3DM(1:N_DM+N_ST))
      WRITE(*,*) 'z velocities (min,max), in c:',
     &     MINVAL(U4DM(1:N_DM+N_ST)),MAXVAL(U4DM(1:N_DM+N_ST))
      WRITE(*,*) 'Masses (min,max), in internal units:',
     &   MINVAL(MASAP(1:N_DM+N_ST)),MAXVAL(MASAP(1:N_DM+N_ST))
      WRITE(*,*) 'Unique IDs (min,max):',
     &   MINVAL(ORIPA(1:N_DM+N_ST)),MAXVAL(ORIPA(1:N_DM+N_ST))
      WRITE(*,*)

      ZETA=ZETA8

      RETURN
      END

*********************************************************************
       SUBROUTINE SORT_DM_PARTICLES(N_DM,NPART_ESP,N_ST,IR_KERN_STARS)
*********************************************************************
*      Reorders DM particles by species (assumes there are N_ESP
*       especies, each 8 times lighter than the previous one)
*********************************************************************
       USE PARTICLES
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER N_DM
       INTEGER NPART_ESP(0:N_ESP-1)
       INTEGER N_ST,IR_KERN_STARS

       INTEGER I,J,K,N,IESP,CONTA
       REAL MLOW,MHIGH,BAS,MAXMASS,MINMASS
       INTEGER,ALLOCATABLE::INDICES(:)
       REAL,ALLOCATABLE::SCR(:,:)
       INTEGER,ALLOCATABLE::SCRINT(:,:)

       WRITE(*,*) 'Sorting particles by mass'

       MAXMASS=MAXVAL(MASAP(1:N_DM))
       MINMASS=MINVAL(MASAP(1:N_DM))
       NPART_ESP=0

       CONTA=0
       ALLOCATE(INDICES(1:N_DM))
       DO IESP=0,N_ESP-1
        MHIGH=2*MAXMASS/8.0**IESP
        MLOW=0.5*MAXMASS/8.0**IESP
        IF (MHIGH.LT.MINMASS) EXIT

        DO I=1,N_DM
         BAS=MASAP(I)
         IF (BAS.LT.MHIGH) THEN
         IF (BAS.GT.MLOW) THEN
          CONTA=CONTA+1
          INDICES(CONTA)=I
         END IF
         END IF
        END DO
        IF (IESP.EQ.0) THEN
         NPART_ESP(IESP)=CONTA
        ELSE
         NPART_ESP(IESP)=CONTA-SUM(NPART_ESP(0:IESP-1))
        END IF
        IF (NPART_ESP(IESP).GT.0) THEN
         WRITE(*,*) 'Of species',IESP,', no. particles:',NPART_ESP(IESP)
        END IF
       END DO

       IF (CONTA.NE.N_DM.OR.SUM(NPART_ESP(0:N_ESP-1)).NE.N_DM) THEN
        WRITE(*,*) 'Wrong sorting, cannot continue',CONTA,N_DM
        STOP
       END IF

       ALLOCATE(SCR(7,N_DM),SCRINT(1,N_DM))

!$OMP PARALLEL DO SHARED(SCR,SCRINT,RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,
!$OMP+                   ORIPA,INDICES,N_DM),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
       DO I=1,N_DM
        SCR(1,I)=RXPA(INDICES(I))
        SCR(2,I)=RYPA(INDICES(I))
        SCR(3,I)=RZPA(INDICES(I))
        SCR(4,I)=U2DM(INDICES(I))
        SCR(5,I)=U3DM(INDICES(I))
        SCR(6,I)=U4DM(INDICES(I))
        SCR(7,I)=MASAP(INDICES(I))
        SCRINT(1,I)=ORIPA(INDICES(I))
       END DO

       DEALLOCATE(INDICES)

!$OMP PARALLEL DO SHARED(SCR,SCRINT,RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,
!$OMP+                   ORIPA,INDICES,N_DM),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
       DO I=1,N_DM
        RXPA(I)=SCR(1,I)
        RYPA(I)=SCR(2,I)
        RZPA(I)=SCR(3,I)
        U2DM(I)=SCR(4,I)
        U3DM(I)=SCR(5,I)
        U4DM(I)=SCR(6,I)
        MASAP(I)=SCR(7,I)
        ORIPA(I)=SCRINT(1,I)
       END DO

       DEALLOCATE(SCR,SCRINT)

       IF (N_ST.GT.0) THEN
        NPART_ESP(IR_KERN_STARS)=NPART_ESP(IR_KERN_STARS)+N_ST
        WRITE(*,*) 'Stars: Of species',IR_KERN_STARS,
     &             ', no. particles:',NPART_ESP(IR_KERN_STARS)
       END IF

C       WRITE(*,*) 'Checking...'
*      CHECK
C       BAS=MASAP(1)
C       DO I=2,N_DM
C        IF (MASAP(I).GT.1.0001*BAS) THEN
C         WRITE(*,*) 'Wrong, I=',I
C         STOP
C        END IF
C        BAS=MASAP(I)
C       END DO
C       stop

       RETURN
       END

*********************************************************************
       SUBROUTINE SORT_DM_PARTICLES_LOCALDENSITY(N_DM,NPART_ESP,N_ST,
     &                                           IR_KERN_STARS,RODO,RE0)
*********************************************************************
*      Assigns DM particles in species, according to local density,
*       and reorders them accordingly.
*********************************************************************
       USE PARTICLES
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER N_DM
       INTEGER NPART_ESP(0:N_ESP-1)
       INTEGER N_ST,IR_KERN_STARS
       REAL*4 RODO,RE0

       REAL*4 DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       INTEGER I,J,K,N,IESP,CONTA,NX,NY,NZ,IX,JY,KZ,PLEV
       REAL MLOW,MHIGH,BAS,MAXMASS,MINMASS
       REAL XMIN,YMIN,ZMIN
       INTEGER,ALLOCATABLE::INDICES(:)
       REAL,ALLOCATABLE::SCR(:,:)
       INTEGER,ALLOCATABLE::SCRINT(:,:)
       REAL,ALLOCATABLE::DENS(:,:,:)
       INTEGER,ALLOCATABLE::MOCKLEVEL(:)

       NX=NMAX
       NY=NMAY
       NZ=NMAZ

       ALLOCATE(DENS(NX,NY,NZ))
!$OMP PARALLEL DO SHARED(NX,NY,NZ,DENS),
!$OMP+            PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO KZ=1,NZ
       DO JY=1,NY
       DO IX=1,NX
        DENS(IX,JY,KZ)=0.0
       END DO
       END DO
       END DO

       XMIN=-NX*DX/2.0
       YMIN=-NY*DY/2.0
       ZMIN=-NZ*DZ/2.0
!$OMP PARALLEL DO SHARED(N_DM,RXPA,RYPA,RZPA,XMIN,YMIN,ZMIN,DX,DY,DZ,
!$OMP+                   MASAP,NX,NY,NZ),
!$OMP+            PRIVATE(IX,JY,KZ,I),
!$OMP+            REDUCTION(+:DENS), DEFAULT(NONE)
       DO I=1,N_DM
        IX=INT((RXPA(I)-XMIN)/DX)+1
        JY=INT((RYPA(I)-YMIN)/DY)+1
        KZ=INT((RZPA(I)-ZMIN)/DZ)+1
        IF (IX.LT.1) IX=1
        IF (IX.GT.NX) IX=NX
        IF (JY.LT.1) JY=1
        IF (JY.GT.NY) JY=NY
        IF (KZ.LT.1) KZ=1
        IF (KZ.GT.NZ) KZ=NZ
        DENS(IX,JY,KZ)=DENS(IX,JY,KZ)+MASAP(I)
       END DO

       BAS=DX*DY*DZ*RODO*RE0**3
!$OMP PARALLEL DO SHARED(NX,NY,NZ,DENS,BAS),
!$OMP+            PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO KZ=1,NZ
       DO JY=1,NY
       DO IX=1,NX
        DENS(IX,JY,KZ)=DENS(IX,JY,KZ)/BAS
       END DO
       END DO
       END DO

       ALLOCATE(MOCKLEVEL(N_DM))

!$OMP PARALLEL DO SHARED(N_DM,RXPA,RYPA,RZPA,XMIN,YMIN,ZMIN,DX,DY,DZ,
!$OMP+                   NX,NY,NZ,DENS,MOCKLEVEL),
!$OMP+            PRIVATE(IX,JY,KZ,I,BAS),
!$OMP+            DEFAULT(NONE)
       DO I=1,N_DM
        IX=INT((RXPA(I)-XMIN)/DX)+1
        JY=INT((RYPA(I)-YMIN)/DY)+1
        KZ=INT((RZPA(I)-ZMIN)/DZ)+1
        IF (IX.LT.1) IX=1
        IF (IX.GT.NX) IX=NX
        IF (JY.LT.1) JY=1
        IF (JY.GT.NY) JY=NY
        IF (KZ.LT.1) KZ=1
        IF (KZ.GT.NZ) KZ=NZ
        BAS=DENS(IX,JY,KZ)
        IF (BAS.GT.0.0) THEN
         MOCKLEVEL(I)=MAX(MIN(INT(LOG(BAS)/LOG(8.0)),N_ESP-1),0)
        ELSE
         MOCKLEVEL(I)=0
        END IF
       END DO

       DEALLOCATE(DENS)

       WRITE(*,*) 'Sorting particles by local density'

       NPART_ESP=0
       CONTA=0

       ALLOCATE(INDICES(1:N_DM))
       DO IESP=0,N_ESP-1
        DO I=1,N_DM
         PLEV=MOCKLEVEL(I)
         IF (PLEV.EQ.IESP) THEN
          CONTA=CONTA+1
          INDICES(CONTA)=I
         END IF
        END DO
        IF (IESP.EQ.0) THEN
         NPART_ESP(IESP)=CONTA
        ELSE
         NPART_ESP(IESP)=CONTA-SUM(NPART_ESP(0:IESP-1))
        END IF
        IF (NPART_ESP(IESP).GT.0) THEN
         WRITE(*,*) 'Of species',IESP,', no. particles:',NPART_ESP(IESP)
        END IF
       END DO

       DEALLOCATE(MOCKLEVEL)

       IF (CONTA.NE.N_DM.OR.SUM(NPART_ESP(0:N_ESP-1)).NE.N_DM) THEN
        WRITE(*,*) 'Wrong sorting, cannot continue',CONTA,N_DM
        STOP
       END IF

       ALLOCATE(SCR(7,N_DM),SCRINT(1,N_DM))

!$OMP PARALLEL DO SHARED(SCR,SCRINT,RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,
!$OMP+                   ORIPA,INDICES,N_DM),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
       DO I=1,N_DM
        SCR(1,I)=RXPA(INDICES(I))
        SCR(2,I)=RYPA(INDICES(I))
        SCR(3,I)=RZPA(INDICES(I))
        SCR(4,I)=U2DM(INDICES(I))
        SCR(5,I)=U3DM(INDICES(I))
        SCR(6,I)=U4DM(INDICES(I))
        SCR(7,I)=MASAP(INDICES(I))
        SCRINT(1,I)=ORIPA(INDICES(I))
       END DO

       DEALLOCATE(INDICES)

!$OMP PARALLEL DO SHARED(SCR,SCRINT,RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,
!$OMP+                   ORIPA,INDICES,N_DM),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
       DO I=1,N_DM
        RXPA(I)=SCR(1,I)
        RYPA(I)=SCR(2,I)
        RZPA(I)=SCR(3,I)
        U2DM(I)=SCR(4,I)
        U3DM(I)=SCR(5,I)
        U4DM(I)=SCR(6,I)
        MASAP(I)=SCR(7,I)
        ORIPA(I)=SCRINT(1,I)
       END DO

       DEALLOCATE(SCR,SCRINT)

       IF (N_ST.GT.0) THEN
        NPART_ESP(IR_KERN_STARS)=NPART_ESP(IR_KERN_STARS)+N_ST
        WRITE(*,*) 'Stars: Of species',IR_KERN_STARS,
     &             ', no. particles:',NPART_ESP(IR_KERN_STARS)
       END IF

C       WRITE(*,*) 'Checking...'
*      CHECK
C       BAS=MASAP(1)
C       DO I=2,N_DM
C        IF (MASAP(I).GT.1.0001*BAS) THEN
C         WRITE(*,*) 'Wrong, I=',I
C         STOP
C        END IF
C        BAS=MASAP(I)
C       END DO
C       stop

       RETURN
       END

*********************************************************************
       SUBROUTINE SORT_DM_PARTICLES_LOCALDENSITY_AND_MASS(N_DM,
     &                   NPART_ESP,N_ST,IR_KERN_STARS,RODO,RE0)
*********************************************************************
*      Assigns DM particles in species, according to local density,
*       and reorders them accordingly.
*********************************************************************
       USE PARTICLES
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER N_DM
       INTEGER NPART_ESP(0:N_ESP-1)
       INTEGER N_ST,IR_KERN_STARS
       REAL*4 RODO,RE0

       REAL*4 DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       INTEGER I,J,K,N,IESP,CONTA,NX,NY,NZ,IX,JY,KZ,PLEV,NN,II,JJ,KK
       INTEGER FAC_GRID,MAXLEV
       REAL MLOW,MHIGH,BAS,MAXMASS,MINMASS,FAC
       REAL XMIN,YMIN,ZMIN,DXPA,DYPA,DZPA
       INTEGER,ALLOCATABLE::INDICES(:)
       REAL,ALLOCATABLE::SCR(:,:)
       INTEGER,ALLOCATABLE::SCRINT(:,:)
       REAL,ALLOCATABLE::DENS(:,:,:)
       INTEGER,ALLOCATABLE::MOCKLEVEL(:)

       MAXMASS=MAXVAL(MASAP(1:N_DM))
       MINMASS=MINVAL(MASAP(1:N_DM))
       FAC_GRID=2
       MAXLEV=FAC_GRID

       NX=NMAX*FAC_GRID
       NY=NMAY*FAC_GRID
       NZ=NMAZ*FAC_GRID

       DXPA=DX/FAC_GRID
       DYPA=DY/FAC_GRID
       DZPA=DZ/FAC_GRID

       ALLOCATE(DENS(NX,NY,NZ))
!$OMP PARALLEL DO SHARED(NX,NY,NZ,DENS),
!$OMP+            PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO KZ=1,NZ
       DO JY=1,NY
       DO IX=1,NX
        DENS(IX,JY,KZ)=0.0
       END DO
       END DO
       END DO

       XMIN=-NX*DXPA/2.0
       YMIN=-NY*DYPA/2.0
       ZMIN=-NZ*DZPA/2.0
!$OMP PARALLEL DO SHARED(N_DM,RXPA,RYPA,RZPA,XMIN,YMIN,ZMIN,DXPA,DYPA,
!$OMP+                   DZPA,MASAP,NX,NY,NZ,MAXMASS,MAXLEV),
!$OMP+            PRIVATE(IX,JY,KZ,I,PLEV,II,JJ,KK,NN,FAC),
!$OMP+            REDUCTION(+:DENS), DEFAULT(NONE)
       DO I=1,N_DM
        IX=INT((RXPA(I)-XMIN)/DXPA)+1
        JY=INT((RYPA(I)-YMIN)/DYPA)+1
        KZ=INT((RZPA(I)-ZMIN)/DZPA)+1
        IF (IX.LT.1) IX=1
        IF (IX.GT.NX) IX=NX
        IF (JY.LT.1) JY=1
        IF (JY.GT.NY) JY=NY
        IF (KZ.LT.1) KZ=1
        IF (KZ.GT.NZ) KZ=NZ
        PLEV=INT(LOG(MAXMASS/MASAP(I))+0.01)
        NN=2**MAX(MAXLEV-PLEV,0)-1
        FAC=1./FLOAT((2*NN+1)**3)
        DO KK=-NN,NN
        DO JJ=-NN,NN
        DO II=-NN,NN
         DENS(IX,JY,KZ)=DENS(IX,JY,KZ)+MASAP(I)
        END DO
        END DO
        END DO


       END DO

       BAS=DX*DY*DZ*RODO*RE0**3
!$OMP PARALLEL DO SHARED(NX,NY,NZ,DENS,BAS),
!$OMP+            PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO KZ=1,NZ
       DO JY=1,NY
       DO IX=1,NX
        DENS(IX,JY,KZ)=DENS(IX,JY,KZ)/BAS
       END DO
       END DO
       END DO

       ALLOCATE(MOCKLEVEL(N_DM))

!$OMP PARALLEL DO SHARED(N_DM,RXPA,RYPA,RZPA,XMIN,YMIN,ZMIN,DXPA,DYPA,
!$OMP+                   DZPA,NX,NY,NZ,DENS,MOCKLEVEL,MAXMASS,MASAP),
!$OMP+            PRIVATE(IX,JY,KZ,I,BAS),
!$OMP+            DEFAULT(NONE)
       DO I=1,N_DM
        IX=INT((RXPA(I)-XMIN)/DXPA)+1
        JY=INT((RYPA(I)-YMIN)/DYPA)+1
        KZ=INT((RZPA(I)-ZMIN)/DZPA)+1
        IF (IX.LT.1) IX=1
        IF (IX.GT.NX) IX=NX
        IF (JY.LT.1) JY=1
        IF (JY.GT.NY) JY=NY
        IF (KZ.LT.1) KZ=1
        IF (KZ.GT.NZ) KZ=NZ
        BAS=DENS(IX,JY,KZ)
        IF (BAS.GT.0.0) THEN
         MOCKLEVEL(I)=MAX(MIN(
     &                   INT(LOG(BAS*MAXMASS/MASAP(I))/LOG(8.0)),
     &                N_ESP-1),0)
        ELSE
         MOCKLEVEL(I)=0
        END IF
       END DO

       DEALLOCATE(DENS)

       WRITE(*,*) 'Sorting particles by local density + particle mass'

       NPART_ESP=0
       CONTA=0

       ALLOCATE(INDICES(1:N_DM))
       DO IESP=0,N_ESP-1
        DO I=1,N_DM
         PLEV=MOCKLEVEL(I)
         IF (PLEV.EQ.IESP) THEN
          CONTA=CONTA+1
          INDICES(CONTA)=I
         END IF
        END DO
        IF (IESP.EQ.0) THEN
         NPART_ESP(IESP)=CONTA
        ELSE
         NPART_ESP(IESP)=CONTA-SUM(NPART_ESP(0:IESP-1))
        END IF
        IF (NPART_ESP(IESP).GT.0) THEN
         WRITE(*,*) 'Of species',IESP,', no. particles:',NPART_ESP(IESP)
        END IF
       END DO

       DEALLOCATE(MOCKLEVEL)

       IF (CONTA.NE.N_DM.OR.SUM(NPART_ESP(0:N_ESP-1)).NE.N_DM) THEN
        WRITE(*,*) 'Wrong sorting, cannot continue',CONTA,N_DM
        STOP
       END IF

       ALLOCATE(SCR(7,N_DM),SCRINT(1,N_DM))

!$OMP PARALLEL DO SHARED(SCR,SCRINT,RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,
!$OMP+                   ORIPA,INDICES,N_DM),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
       DO I=1,N_DM
        SCR(1,I)=RXPA(INDICES(I))
        SCR(2,I)=RYPA(INDICES(I))
        SCR(3,I)=RZPA(INDICES(I))
        SCR(4,I)=U2DM(INDICES(I))
        SCR(5,I)=U3DM(INDICES(I))
        SCR(6,I)=U4DM(INDICES(I))
        SCR(7,I)=MASAP(INDICES(I))
        SCRINT(1,I)=ORIPA(INDICES(I))
       END DO

       DEALLOCATE(INDICES)

!$OMP PARALLEL DO SHARED(SCR,SCRINT,RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,
!$OMP+                   ORIPA,INDICES,N_DM),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
       DO I=1,N_DM
        RXPA(I)=SCR(1,I)
        RYPA(I)=SCR(2,I)
        RZPA(I)=SCR(3,I)
        U2DM(I)=SCR(4,I)
        U3DM(I)=SCR(5,I)
        U4DM(I)=SCR(6,I)
        MASAP(I)=SCR(7,I)
        ORIPA(I)=SCRINT(1,I)
       END DO

       DEALLOCATE(SCR,SCRINT)

       IF (N_ST.GT.0) THEN
        NPART_ESP(IR_KERN_STARS)=NPART_ESP(IR_KERN_STARS)+N_ST
        WRITE(*,*) 'Stars: Of species',IR_KERN_STARS,
     &             ', no. particles:',NPART_ESP(IR_KERN_STARS)
       END IF

C       WRITE(*,*) 'Checking...'
*      CHECK
C       BAS=MASAP(1)
C       DO I=2,N_DM
C        IF (MASAP(I).GT.1.0001*BAS) THEN
C         WRITE(*,*) 'Wrong, I=',I
C         STOP
C        END IF
C        BAS=MASAP(I)
C       END DO
C       stop

       RETURN
       END
