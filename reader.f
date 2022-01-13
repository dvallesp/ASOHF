*********************************************************************
       SUBROUTINE READ_MASCLET(VAR,ITER,NX,NY,NZ,NDXYZ,T,ZETA,NL,NPATCH,
     &           PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &           PATCHRX,PATCHRY,PATCHRZ,MAP,
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

       CHARACTER*9 FILNOM1,FILNOM2,FILNOM4
       CHARACTER*10 FILNOM3
       CHARACTER*25 FIL1,FIL2,FIL4
       CHARACTER*26 FIL3

*      VARIABLES
       REAL*4 U1(NMAX,NMAY,NMAZ)
       REAL*4 U1G(NMAX,NMAY,NMAZ)
       REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL*4 U11G(NAMRX,NAMRY,NAMRZ,NPALEV)
       COMMON /VARIA/ U1,U11,U1G,U11G  !!,U12,U13,U14

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
       REAL*4 U2DM(PARTIRED)
       REAL*4 U3DM(PARTIRED)
       REAL*4 U4DM(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       REAL*4 RXPA(PARTIRED)
       REAL*4 RYPA(PARTIRED)
       REAL*4 RZPA(PARTIRED)

       INTEGER ORIPA(PARTIRED)

       REAL*4, ALLOCATABLE::SCR(:,:,:)

       REAL*4 UBAS(0:PARTIRED)
       INTEGER UBAS2(0:PARTIRED),CONTA,LOW1,LOW2


*      READING DATA
       CALL NOMFILE(ITER,FILNOM1,FILNOM2,FILNOM3,FILNOM4)
       WRITE(*,*) 'Reading iter',ITER,' ',FILNOM1,FILNOM2,FILNOM3,
     &            FILNOM4

       FIL1='simu_masclet/'//FILNOM1
       FIL2='simu_masclet/'//FILNOM2
       FIL3='simu_masclet/'//FILNOM3
       FIL4='simu_masclet/'//FILNOM4


       OPEN (33,FILE=FIL3,STATUS='UNKNOWN',ACTION='READ')
       OPEN (31,FILE=FIL1,
     &       STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')
       OPEN (32,FILE=FIL2,
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
        READ(31) (((U1G(I,J,K),I=1,N1),J=1,N2),K=1,N3)
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
        READ(31) (((SCR(IX,J,K),IX=1,N1),J=1,N2),K=1,N3)
           U11G(1:N1,1:N2,1:N3,I)=SCR(1:N1,1:N2,1:N3)
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
       SUBROUTINE READ_PARTICLES_MASCLET(ITER,NX,NY,NZ,T,ZETA,MAP,
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

       INTEGER I,J,K,IX,NL,IR,IRR,N1,N2,N3,N_DM,N_ST,NBAS,ARE_BH

       CHARACTER*9 FILNOM1,FILNOM2,FILNOM4
       CHARACTER*10 FILNOM3
       CHARACTER*25 FIL1,FIL2,FIL4
       CHARACTER*26 FIL3

       INTEGER,ALLOCATABLE::NPATCH(:)

       INTEGER,ALLOCATABLE::NPART(:),NPARTST(:),NPARTBH(:)
       REAL*4 U2DM(PARTIRED)
       REAL*4 U3DM(PARTIRED)
       REAL*4 U4DM(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       REAL*4 RXPA(PARTIRED)
       REAL*4 RYPA(PARTIRED)
       REAL*4 RZPA(PARTIRED)

       INTEGER ORIPA(PARTIRED)

       REAL*4 UBAS(0:PARTIRED)
       INTEGER UBAS2(0:PARTIRED),CONTA,LOW1,LOW2

       ARE_BH=0 ! Depends on MASCLET version (are there BHs??)

*      READING DATA
       CALL NOMFILE(ITER,FILNOM1,FILNOM2,FILNOM3,FILNOM4)
       WRITE(*,*) 'Reading iter',ITER,' ',FILNOM1,FILNOM2,FILNOM3,
     &                                    FILNOM4

       FIL1='simu_masclet/'//FILNOM1
       FIL2='simu_masclet/'//FILNOM2
       FIL3='simu_masclet/'//FILNOM3
       FIL4='simu_masclet/'//FILNOM4

*      GRID DATA
       OPEN (33,FILE=FIL3,STATUS='UNKNOWN',ACTION='READ')
       READ(33,*) IRR,T,NL,MAP
       READ(33,*) ZETA
       READ(33,*) IR,NDXYZ
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

       CLOSE(33)

**     DARK MATTER
       OPEN (32,FILE=FIL2,
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
***       END IF
       WRITE(*,*) 'TOTAL DM PARTICLES IN ITER=',CONTA
       IF (CONTA.NE.N_DM) THEN
        WRITE(*,*) 'WARNING: CONTA != N_DM',CONTA,N_DM
        STOP
       END IF

       IF (VAR.EQ.2) THEN
        N_ST=SUM(NPARTST(0:NL))+SUM(NPARTBH(0:NL))
        OPEN (34,FILE=FIL4,
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
       SUBROUTINE READ_PARTICLES_GENERAL(ITER,NX,NY,NZ,T,ZETA,NL,MAP,
     &            U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,ORIPA,LADO0,
     &            N_GAS,N_DM,N_PARTICLES)
*********************************************************************
*      Reads generic file with particles
*********************************************************************

       IMPLICIT NONE

       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NX,NY,NZ,ITER,NDXYZ
       REAL*4 T,AAA,BBB,CCC,MAP,ZETA,LADO0

       INTEGER I,J,K,IX,JY,KZ,NL,IR,IRR,VAR,N1,N2,N3
       INTEGER N_DM,N_GAS,N_PARTICLES
       INTEGER NBAS_DM,NBAS_GAS,NBAS_PARTICLES

       CHARACTER*9 FILNOM1,FILNOM2
       CHARACTER*10 FILNOM3
       CHARACTER*25 FIL1,FIL2
       CHARACTER*26 FIL3

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
       REAL*4 U2DM(PARTIRED)
       REAL*4 U3DM(PARTIRED)
       REAL*4 U4DM(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       REAL*4 RXPA(PARTIRED)
       REAL*4 RYPA(PARTIRED)
       REAL*4 RZPA(PARTIRED)

       INTEGER ORIPA(PARTIRED)

       REAL*4 UBAS(0:PARTIRED)
       INTEGER UBAS2(0:PARTIRED),CONTA,LOW1,LOW2

       OPEN (5,FILE='particle_list.dat',
     &               STATUS='UNKNOWN',ACTION='READ')

       READ(5,*) NBAS_PARTICLES,ZETA,T,NBAS_GAS,NBAS_DM

       IF (NBAS_PARTICLES.NE.N_PARTICLES.OR.
     &     NBAS_GAS.NE.N_GAS.OR.NBAS_DM.NE.N_DM) THEN
        WRITE(*,*) 'WARNING: number of particles badly specified'
        WRITE(*,*) NBAS_GAS,NBAS_DM,NBAS_PARTICLES
        WRITE(*,*) N_GAS,N_DM,N_PARTICLES
        STOP
       END IF

       WRITE(*,*)'N_DM,N_GAS=',N_DM,N_GAS

       DO I=1,N_DM
        READ(5,*) ORIPA(I),RXPA(I),RYPA(I),RZPA(I),
     &            U2DM(I),U3DM(I),U4DM(I),MASAP(I)
       END DO
       DO I=N_DM+1,N_DM+N_GAS
        READ(5,*) ORIPA(I),RXPA(I),RYPA(I),RZPA(I),
     &            U2DM(I),U3DM(I),U4DM(I),MASAP(I)
       END DO

       WRITE(*,*) 'ORIPA=',MAXVAL(ORIPA(1:N_PARTICLES)),
     &                      MINVAL(ORIPA(1:N_PARTICLES))

       CLOSE(5)

***********************
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

*      In [-L/2, L/2]
       RXPA(1:N_DM)=RXPA(1:N_DM)-LADO0*0.5
       RYPA(1:N_DM)=RYPA(1:N_DM)-LADO0*0.5
       RZPA(1:N_DM)=RZPA(1:N_DM)-LADO0*0.5

       WRITE(*,*)'CHECKING BOX SIDE', LADO0
       IX=0
       IX=COUNT(ABS(RXPA(1:N_DM)).GT.LADO0*0.5001)
       JY=0
       JY=COUNT(ABS(RYPA(1:N_DM)).GT.LADO0*0.5001)
       KZ=0
       KZ=COUNT(ABS(RZPA(1:N_DM)).GT.LADO0*0.5001)
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
       WRITE(*,*) MAXVAL(RXPA(N_DM+1:N_DM+N_GAS)),
     &            MINVAL(RXPA(N_DM+1:N_DM+N_GAS))
       WRITE(*,*) MAXVAL(RYPA(N_DM+1:N_DM+N_GAS)),
     &            MINVAL(RYPA(N_DM+1:N_DM+N_GAS))
       WRITE(*,*) MAXVAL(RZPA(N_DM+1:N_DM+N_GAS)),
     &            MINVAL(RZPA(N_DM+1:N_DM+N_GAS))

*      In [-L/2, L/2]
       RXPA(N_DM+1:N_DM+N_GAS)=RXPA(N_DM+1:N_DM+N_GAS)-LADO0*0.5
       RYPA(N_DM+1:N_DM+N_GAS)=RYPA(N_DM+1:N_DM+N_GAS)-LADO0*0.5
       RZPA(N_DM+1:N_DM+N_GAS)=RZPA(N_DM+1:N_DM+N_GAS)-LADO0*0.5

       WRITE(*,*)'CHECKING BOX SIDE', LADO0
       IX=0
       IX=COUNT(ABS(RXPA(N_DM+1:N_DM+N_GAS)).GT.LADO0*0.5001)
       JY=0
       JY=COUNT(ABS(RYPA(N_DM+1:N_DM+N_GAS)).GT.LADO0*0.5001)
       KZ=0
       KZ=COUNT(ABS(RZPA(N_DM+1:N_DM+N_GAS)).GT.LADO0*0.5001)
       IF(IX.GT.0.OR.JY.GT.0.OR.KZ.GT.0) THEN
        WRITE(*,*) IX,JY,KZ
        WRITE(*,*)'WARNING DM: LADO!!'
        STOP
       ENDIF

       WRITE(*,*) '///GAS///'
       WRITE(*,*) MAXVAL(RXPA(N_DM+1:N_DM+N_GAS)),
     &            MINVAL(RXPA(N_DM+1:N_DM+N_GAS))
       WRITE(*,*) MAXVAL(RYPA(N_DM+1:N_DM+N_GAS)),
     &            MINVAL(RYPA(N_DM+1:N_DM+N_GAS))
       WRITE(*,*) MAXVAL(RZPA(N_DM+1:N_DM+N_GAS)),
     &            MINVAL(RZPA(N_DM+1:N_DM+N_GAS))

       END IF !N_GAS

       WRITE(*,*) 'TOTAL PARTICLES IN ITER=',N_DM+N_GAS

       RETURN
       END


*********************************************************************
       SUBROUTINE SORT_DM_PARTICLES(U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,
     &                              RZPA,ORIPA,N_DM,NPART_ESP,N_ST,
     &                              IR_KERN_STARS)
*********************************************************************
*      Reorders DM particles by species (assumes there are N_ESP
*       especies, each 8 times lighter than the previous one)
*********************************************************************

       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       REAL*4 U2DM(PARTIRED),U3DM(PARTIRED),U4DM(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED)
       INTEGER ORIPA(PARTIRED)
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
