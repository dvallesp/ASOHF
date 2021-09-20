************************************************************************
       SUBROUTINE MESHRENOEF(ITER,NX,NY,NZ,NL,COTA,NPATCH,
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



************************************
*      SEPARO LAS PARTI EN ESPECIES
************************************
       NDM_ESP(1)=NDMXYZ        !!!only DM particles
CX        NDM_ESP(1)=NPARTICULAS

       DO I=1, N_ESP
        WRITE(*,*) 'Particles of especie ', I,'= ',NDM_ESP(I)
       END DO

       IF (N_DM.NE.SUM(NDM_ESP(1:N_ESP))) THEN
        WRITE(*,*) 'WARNING: ESPECIES!'
        STOP
       END IF
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

!!!!! Esta rutina ha sido modificada (26/06/2019) para usar variables
!!!!! enteras en lugar de reales en la comprobación de qué parches
!!!!! solapan con otros.

!!!!! En este codigo CONTA2 es SOLAP. Y lo que era MARCA ya no se usa.
!!!!! Ademas, paso VECINO y NVECI en el argumento de la rutina

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

       INTEGER CR1,CR2,CR3,CR4,CR5,CR6
       INTEGER IR,I,J,IX,JY,KZ,II,JJ,KK
       INTEGER N1,N2,N3,L1,L2,L3, NL
       INTEGER NN1,NN2,NN3,LL1,LL2,LL3
       INTEGER KZ2,JY2,IX2,I2
       INTEGER NV,A2,B2,C2,K,LOW1, LOW2

       INTEGER VECINO(NPALEV,NPALEV),NVECI(NPALEV)
       INTEGER SOLAP(NAMRX,NAMRY,NAMRZ,NPALEV)

       REAL*4 A1,B1,C1,RIV1,RIV2,RIV3
       INTEGER CONTROL
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
! parallelize
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
