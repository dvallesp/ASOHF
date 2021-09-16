*********************************************************************
       SUBROUTINE LEER(VAR,ITER,NX,NY,NZ,NDXYZ,T,ZETA,NL,NPATCH,
     &           PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &           PATCHRX,PATCHRY,PATCHRZ,MAP,
     &           U2DM,U3DM, U4DM,MASAP,NPART,RXPA,RYPA,RZPA)
*********************************************************************
*      Reads MASCLET data: grids, gas density (clus files) and
*      DM particles information.
*      Must be checked depending on the version/flavour of MASCLET
*      the simulation has been run with
*********************************************************************

       IMPLICIT NONE

       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NX,NY,NZ,ITER,NDXYZ
       REAL*4 T,AAA,BBB,CCC,MAP,ZETA

       INTEGER I,J,K,IX,NL,IR,IRR,VAR,N1,N2,N3

       CHARACTER*9 FILNOM1,FILNOM2
       CHARACTER*10 FILNOM3
       CHARACTER*25 FIL1,FIL2
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

       INTEGER ORIPA1(PARTIRED),ORIPA2(PARTIRED)
       COMMON /PUNTEROS/ ORIPA1, ORIPA2

       REAL*4, ALLOCATABLE::SCR(:,:,:)

       REAL*4 UBAS(0:PARTI)
       INTEGER UBAS2(0:PARTI),CONTA,LOW1,LOW2


*      READING DATA
       CALL NOMFILE(ITER,FILNOM1,FILNOM2,FILNOM3)
       WRITE(*,*) 'Reading iter',ITER,' ',FILNOM1,FILNOM2,FILNOM3

       FIL1='simu_masclet/'//FILNOM1
       FIL2='simu_masclet/'//FILNOM2
       FIL3='simu_masclet/'//FILNOM3


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
       WRITE(*,*), 'NPATCH(IR), NPART(IR)',NPATCH(IR), NPART(IR)
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
C        READ(32) (ORIPA1(I),I=1,NDXYZ)   !OJO! las nuevas versioens de MASCLET no lo tienen
        READ(32) (ORIPA2(I),I=1,NDXYZ)    !particle ID
        CONTA=NDXYZ
        MASAP(1:NDXYZ)=MAP
        WRITE(*,*) 'ORIPA1=',MAXVAL(ORIPA1(1:NDXYZ)),
     &                       MINVAL(ORIPA1(1:NDXYZ))
        WRITE(*,*) 'ORIPA2=',MAXVAL(ORIPA2(1:NDXYZ)),
     &                       MINVAL(ORIPA2(1:NDXYZ))
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
c        READ(32) (UBAS2(IX),IX=1,NPART(IR))         !OJO! las nuevas versioens de MASCLET no lo tienen
c        IF (NPART(IR).GT.0)
c     &  ORIPA1(CONTA+1:CONTA+NPART(IR))=UBAS2(1:NPART(IR))
        READ(32) (UBAS2(IX),IX=1,NPART(IR))
        IF (NPART(IR).GT.0)
     &   ORIPA2(CONTA+1:CONTA+NPART(IR))=UBAS2(1:NPART(IR))

        IF (NPART(IR).GT.0) THEN
        WRITE(*,*) 'ORIPA1=',MAXVAL(ORIPA1(CONTA+1:CONTA+NPART(IR))),
     &                       MINVAL(ORIPA1(CONTA+1:CONTA+NPART(IR)))
        WRITE(*,*) 'ORIPA2=',MAXVAL(ORIPA2(CONTA+1:CONTA+NPART(IR))),
     &                       MINVAL(ORIPA2(CONTA+1:CONTA+NPART(IR)))
        END IF

        CONTA=CONTA+NPART(IR)
        WRITE(*,*) 'NPART(IR)=',IR,NPART(IR),CONTA


       END DO

       DEALLOCATE(SCR)

       CLOSE(32)
***       END IF

       CONTA=SUM(NPART(0:NL))
       WRITE(*,*) 'TOTAL PARTICLES IN ITER=',CONTA

       RETURN
       END
