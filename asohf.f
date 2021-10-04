**********************************************************************
       PROGRAM ASOHF_v18
***********************************************************************
***********************************************************************
*      ASOHF IS AN ADAPTIVE SPHERICAL OVERDENSITY HALO FINDER.
*      IT FINDS HALOES AND SUBHALOES IN A COSMOLOGICAL SIMULATION,
*      BUILDS THEIR MERGING HISTORIES AND COMPUTES THEIR MAIN PHYSICAL
*      PROPERTIES.
***********************************************************************
*      For further details: Planelles & Quilis, 2010. A&A, 519:A94
***********************************************************************
*
*      Version: asohf_paralel_v8_NewVeins.f
*      (masclet1:/scratch/planelles/project_statistics/ASOHF/ASOHF_SA)
*
*---------------------GENERAL CONSIDERATIONS--------------------------*
*      In this version the merger tree is done every 2 iterations and
*      therefore, MAXITER=2 and MAXITER2=real num. of iters (MARK)
*
*      # INPUT FILES:
*           1)asohf_parameters.dat
*           2)asohf.dat
*
*      # OUTPUT FILES:
*           1)Families00000 ---> general information of all haloes
*           2)Merger_t00000 (optional, PLOT=2) ---> with %
*           3)Merger_r00000 (optional, PLOT=3) ---> main line
*           4)Grid_asohf00000 (only when FLAG_SA=0)
*
*           If REALCLUS(IFI,I)= -1  ------> HALO
*                             = 0   ------> RUBBISH (double or poor)
*                             = #>0 ------> SUBHALO within #
*
*           Proceeding way: 1)Halo refinement with DM particles
*                           2)Halo classification
*
*      # DATA TO READ:
*           1) MASCLET simulation in the "simu_masclet" directory
*           2) External file of particles: "particle_list.dat".
*              File in ascci containing:
*                #1st line        ---> NPARTICULAS,ZETA,T,N_GAS,N_DM
*                #1 line/particle ---> ID,X,Y,Z,VX,VY,VZ,MASS
*                 - Positions should go in units of Mpc within [0, L_box]
*                 - Velocities should go in units of c.
*                 - Masses should go in units of 9.1717e18 M_sun
*
*------------------------PREVIOUS SETTINGS----------------------------*
*      FLAG_SA     --> stand-alone halo finder(=0) or MASCLETs grid(=1)
*      FLAG_MASCLET--> MASCLET "as" stand-alone (=1)(it needs flag_sa=0)
*      FLAG_GAS    --> there is gas(=1) or only DM(=0) (see VAR!!)
*
*      If PLOT=2---> MERGER TREE WITH % (ONLY FOR HALOES AT IR=0)
*      If PLOT=3---> REDUCED MERGER TREE (ONLY FOR HALOES AT IR=0)
*
*      OJO!!! ONLY INTERNAL USE!!! hay que cambiar a mano:
*             1) El valor de las distintas COTAS
*             2) Rellenar las distintas especies de particulas
*             3) Tal y como esta, estamos pasando del gas al
*                construir la malla aunque si que lo leemos
**********************************************************************

       USE MTREE
       IMPLICIT NONE

       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER I,J,K,I2,J2
       INTEGER NX,NY,NZ,ITER,NDXYZ
       REAL*4 T,TEI

       REAL*4  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
       COMMON /GRID/   RADX,RADY,RADZ

       REAL*4  RX(0:NAMRX+1,NPALEV),RY(0:NAMRX+1,NPALEV),
     &         RZ(0:NAMRX+1,NPALEV)
       COMMON /GRIDAMR/ RX,RY,RZ

       REAL*4 PI,ACHE,T0,RE0,PI4ROD
       COMMON /DOS/ACHE,T0,RE0
       REAL*4 UNTERCIO,CGR,CGR2,ZI,RODO,ROI,REI,LADO,LADO0
       COMMON /CONS/PI4ROD,REI,CGR,PI
       REAL*4 OMEGA0

       REAL*4 RETE,HTE,ROTE
       COMMON /BACK/ RETE,HTE,ROTE

       REAL*4 DX,DY,DZ,H2
       COMMON /ESPACIADO/ DX,DY,DZ

       REAL*4 UV, UM

       INTEGER IX,JY,KZ,NL,IR,L1
       REAL*4 RX2,RY2,RZ2,A1,A2,B1,C1,A3,A4
       REAL*4 DXPA,DYPA,DZPA

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

       INTEGER NPART(0:NLEVELS)
       REAL*4 U2DM(PARTIRED),U3DM(PARTIRED),U4DM(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED)

       INTEGER CONTA2(NAMRX,NAMRY,NAMRZ,NPALEV) !this will be deleted
       INTEGER SOLAP(NAMRX,NAMRY,NAMRZ,NPALEV)

       INTEGER ORIPA1(PARTIRED),ORIPA2(PARTIRED)
       COMMON /PUNTEROS/ ORIPA1, ORIPA2

       INTEGER II,JJ,KK,KK1,KK2,KKK2,ICEN(3),IR2,III,KK3,KK4
       INTEGER IX1,IX2,N1,N2,N3,NTOT,VAR,IX3
       INTEGER NFILE,FIRST,EVERY,IFI,LAST,PLOT
       REAL*4 PRUEBAX,PRUEBAY,PRUEBAZ,DELTA_PRUEBA
       REAL*4 ZETA,AA1,AA2,LIM,BAS, BAS1, BAS2
       REAL*4 RRRR,R111,R222

       CHARACTER*13 FILE1, FILE4
       CHARACTER*14 FILE3, FILE7
       CHARACTER*30 FILERR3,FILERR7
       INTEGER*4 DATE(3), TIME(3), CONTAERR
       INTEGER CONTAITER,MARK(MAXITER2),NFILE2,IFI2

       INTEGER VECINO(NPALEV,NPALEV),NVECI(NPALEV),CEN(1),CEL
       INTEGER NV, NV2,FLAG,VID(NPALEV),NV3
       REAL*4 MAXIMO(NPALEV),MAXIMO2(NPALEV),CONTAP(NPALEV)
       INTEGER CEN2(1),CEL2

       REAL*4 UBAS1(NMAX,NMAY,NMAZ)
       REAL*4 RMIN,VOL,RSHELL,ESP, REF, RANT,VOL1,VOL2
       INTEGER SHELL, NSHELL,CONTA(NMAX,NMAY,NMAZ)
       INTEGER NSHELL_2
       REAL*4 DELTA, AA, VR, AADM, MASADM, MAP
       REAL*4 CONTRASTEC, CONTRASTEX, OMEGAZ
       REAL*4 VESC2, RS, CONCEN, VMAX2, F2, FC

       INTEGER NCLUS,ICMIN,ICMAX
       INTEGER SOLAPA(MAXNCLUS,NMAXSOLAP), NSOLAP(MAXNCLUS)
       INTEGER KONTA3,KONTA,KONTACAPA,KONTA1,KONTA2
       REAL*4 DIS, VKK, BASMAS, MASAKK
       REAL*4 CMX, CMY, CMZ, SOLMAS, MASAMIN, DELTA2
       REAL*4 VCMX, VCMY, VCMZ, VCM,VVV2, MASA2
       REAL*4 MASA(MAXNCLUS), RADIO(MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)

       REAL*8 CMX8,CMY8,CMZ8,VCMX8,VCMY8,VCMZ8

*      ---HALOS Y SUBHALOS---
       INTEGER NNCLUS(MAXITER)
       REAL*4 NR(MAXITER,NMAXNCLUS)
       REAL*4 CONCENTRA(MAXITER,NMAXNCLUS)
       REAL*4 ANGULARM(MAXITER,NMAXNCLUS)
       REAL*4 NMASA(MAXITER,NMAXNCLUS)
       REAL*4 NCLUSRX(MAXITER,NMAXNCLUS)
       REAL*4 NCLUSRY(MAXITER,NMAXNCLUS)
       REAL*4 NCLUSRZ(MAXITER,NMAXNCLUS)
       REAL*4 VMAXCLUS(MAXITER,NMAXNCLUS)
       REAL*4 VCM2(MAXITER,NMAXNCLUS)

       INTEGER IPLIP(MAXITER,NMAXNCLUS)
       INTEGER IPLIR(MAXITER,NMAXNCLUS)
       REAL*4 VX(MAXITER,NMAXNCLUS)
       REAL*4 VY(MAXITER,NMAXNCLUS)
       REAL*4 VZ(MAXITER,NMAXNCLUS)

       REAL*4 VCMAX(MAXITER,NMAXNCLUS)
       REAL*4 MCMAX(MAXITER,NMAXNCLUS)
       REAL*4 RCMAX(MAXITER,NMAXNCLUS)

       REAL*4 M200(MAXITER,NMAXNCLUS)
       REAL*4 R200(MAXITER,NMAXNCLUS)
       INTEGER FLAG_200

       INTEGER NCAPAS(NMAXNCLUS)
       INTEGER REALCLUS(MAXITER,MAXNCLUS)

*      ---PARTICULAS E ITERACIONES---
       INTEGER DMPITER(MAXITER)
       REAL*4 ZETAS(MAXITER),TIEMPO(MAXITER)
       INTEGER LIP(PARTIRED),LIR(PARTIRED),CONTADM(PARTIRED)
       INTEGER DMPCLUS(MAXITER,NMAXNCLUS)

*       ----- module mtree.f -----
*       INTEGER DMLIP(MAXITER,PARTIRED_PLOT)
*       INTEGER DMLIR(MAXITER,PARTIRED_PLOT)
*       INTEGER NHOST(MAXITER,0:ININL,PARTIRED)
*       INTEGER HHOST(MAXITER,NMAXDAD,0:ININL,PARTIRED)
*       REAL*4 NEW_MASAP(MAXITER,PARTIRED_PLOT)
*       REAL*4 RATIO(MAXITER,NMAXDAD,NMAXCLUS_PLOT)
*       INTEGER DAD(MAXITER,NMAXDAD,NMAXCLUS_PLOT)
*       INTEGER NDAD(MAXITER,NMAXNCLUS)

*      ---SUBSTRUCTURE---
       INTEGER LEVHAL(MAXNCLUS),NHALLEV(0:NLEVELS)
       INTEGER NLEVHAL(MAXITER,NMAXNCLUS)
       INTEGER SUBHALOS(NMAXNCLUS)

*      ---HALO SHAPE---
       REAL*4 INERTIA(3,3),EIGENVAL(MAXITER,3,NMAXNCLUS)
       REAL*4 BASEIGENVAL(3),AADMX(3)
       INTEGER DIMEN,NROT

*      ---STAND-ALONE HALO FINDER---
       INTEGER FLAG_SA,FLAG_GAS,FLAG_MASCLET,FLAG_WDM
       INTEGER N_GAS,N_DM,N_PARTICLES
       REAL*4 COTA(NCOTAS,0:NLEVELS)

*      ---UNBINDING---
       REAL*4 REF_MIN, REF_MAX
       INTEGER SALIDA,FAC
       REAL*4 DISTA(0:PARTIRED), DISTA2(0:PARTIRED)
       INTEGER QUIEN(PARTIRED), QUIEN2(PARTIRED)
       INTEGER INDICE(PARTIRED)

       REAL*4 DENSITOT(0:PARTIRED/10),RADIAL(0:PARTIRED/10)
       REAL*4 DENSA,DENSB,DENSC,NORMA

       REAL*4 KK_REAL,KK_REAL_2
       INTEGER KK_ENTERO,KK_ENTERO_2, NV_GOOD

       REAL*4, ALLOCATABLE::DDD(:)
       INTEGER, ALLOCATABLE::DDDP(:)
       INTEGER, ALLOCATABLE::DDDX(:)
       INTEGER, ALLOCATABLE::DDDY(:)
       INTEGER, ALLOCATABLE::DDDZ(:)
       INTEGER,ALLOCATABLE::INDICE2(:)
       REAL*4, ALLOCATABLE::DDD2(:)
       INTEGER, ALLOCATABLE::DDDP2(:)
       INTEGER, ALLOCATABLE::DDDX2(:)
       INTEGER, ALLOCATABLE::DDDY2(:)
       INTEGER, ALLOCATABLE::DDDZ2(:)

       REAL*4 BOUND
       INTEGER NX1,NX2,NY1,NY2,NZ1,NZ2, BORDES
       INTEGER LOW1,LOW2

       INTEGER PARCHLIM       ! =0, sin limite por nivel, =/ 0 limites
       INTEGER MPAPOLEV(100)  !maximo numero de parches por nivel, 100 niveles max.
       COMMON /LIMPARCH/ PARCHLIM,MPAPOLEV

*      ---PARALLEL---
       INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
       COMMON /PROCESADORES/ NUM

       INTEGER NMAXNCLUSBAS,PABAS,NPBAS,NLEVBAS,NUMPARTBAS,NBASPART_PLOT

       INTEGER, ALLOCATABLE:: IP_PARAL(:,:)
       INTEGER, ALLOCATABLE:: IR_PARAL(:,:)
       REAL*4, ALLOCATABLE:: MASAP_PARAL(:,:)
       INTEGER NN, IP, NH, IP2

       INTEGER TEMP0(NMAX,NMAY,NMAZ)
       INTEGER TEMP1(NAMRX,NAMRY,NAMRZ,NPALEV)

**************************************************************
*      OPENING FILES
**************************************************************
       OPEN(1,FILE='./input_files/asohf.dat',
     &              STATUS='UNKNOWN',ACTION='READ')


*      READING INITIAL DATA
****************************************************
*     NX,NY,NZ < or = NMAX,NMAY,NMAZ               *
****************************************************
       READ(1,*)
       READ(1,*)
       READ(1,*) FIRST,LAST,EVERY
       READ(1,*)
       READ(1,*) NX,NY,NZ
       READ(1,*)
       READ(1,*) N_GAS,N_DM
       READ(1,*)
       READ(1,*) ACHE,OMEGA0
       READ(1,*)
       READ(1,*) ZI,LADO0
       READ(1,*)
       READ(1,*) FLAG_PARALLEL,NUM
       READ(1,*)
       READ(1,*) NL
       READ(1,*)
       READ(1,*) PARCHLIM
       READ(1,*)
       READ(1,*) LIM
       READ(1,*)
       READ(1,*) VAR
       READ(1,*)
       READ(1,*) PLOT
       READ(1,*)
       READ(1,*) NSHELL
       READ(1,*)
       READ(1,*) BOUND, BORDES
       READ(1,*)
       READ(1,*) FLAG_SA,FLAG_MASCLET,FLAG_GAS
       READ(1,*)
       READ(1,*) FLAG_WDM

       CLOSE(1)

       N_PARTICLES=N_GAS+N_DM
       H2=ACHE
       MPAPOLEV(1:100)=LIM   !max. num. of patches per level, 100 levels max.

**************************************************************
*     ...PARALLEL RUNNING...
!$OMP PARALLEL SHARED(NUM)
!$OMP SINGLE
!$      NUM=OMP_GET_NUM_THREADS()
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL
**************************************************************

       CALL IDATE(DATE)
       CALL ITIME(TIME)
       WRITE(*,*) 'DATE=',DATE(1),'/',DATE(2),'/',DATE(3)
       WRITE(*,*) 'TIME=',TIME(1),':',TIME(2),':',TIME(3)

       WRITE(*,*) '************************************************'
       WRITE(*,*) '             GENERAL SETTINGS                   '
       WRITE(*,*) '************************************************'

       IF (FLAG_PARALLEL.EQ.1) THEN
         WRITE(*,*)'Running in PARALLEL in',NUM, 'processors'
       END IF
       IF (FLAG_PARALLEL.EQ.0)  WRITE(*,*)'Running in SERIAL...'

       IF(FLAG_SA.EQ.0.AND.FLAG_MASCLET.EQ.0)
     &       WRITE(*,*) 'ASOHF like stand-alone Halo Finder...'
       IF(FLAG_SA.EQ.0.AND.FLAG_MASCLET.EQ.1)
     &       WRITE(*,*) 'ASOHF reading MASCLETs PARTICLES...'
       IF(FLAG_SA.EQ.1)
     &       WRITE(*,*) 'ASOHF reading MASCLETs GRID...'

       IF(VAR.EQ.1) WRITE(*,*) 'Analysing only GAS'
       IF(VAR.EQ.2) WRITE(*,*) 'Analysing only DM'
       IF(VAR.EQ.0) WRITE(*,*) 'Analysing GAS+DM'

       WRITE(*,*) 'Min. number of particles per halo ', NUMPART

       WRITE(*,*) '************************************************'


***************************
*      GRID BUILDER
***************************
       LADO=LADO0-(LADO0/NX)
       CALL MALLA(NX,NY,NZ,LADO)

       WRITE(*,*) '************************************************'
       WRITE(*,*) '             GRID                               '
       WRITE(*,*) '************************************************'
       WRITE(*,*) 'SIDE LENGTH=',LADO
       WRITE(*,*) 'NX,DX,RADX(1),RADX(NX)=',NX,DX,RADX(1),RADX(NX)
       WRITE(*,*) 'NUMBER OF PATCHES PER LEVEL:'
       IF (PARCHLIM.EQ.0) WRITE(*,*) '  No limited patches per level!'
       IF (PARCHLIM.NE.0) THEN
             WRITE(*,*) '  Limit patches per level=', MPAPOLEV(1)
       END IF
       WRITE(*,*) '************************************************'



*********************************************************************
*      COSMOLOGICAL BACKGROUND
*********************************************************************
       PI=DACOS(-1.D0)
       UNTERCIO=1.D0/3.D0
       CGR=1.D0/(8.D0*PI)
       CGR2=2.D0*CGR
       ACHE=ACHE*3.66D-3
*      T0=364.298725        !T0=ACTUAL TIME
       T0=2.D0*UNTERCIO/ACHE
       RODO=OMEGA0*3.D0*ACHE**2
*      scale factor must be = 1Mpc  at z=0  in order to be consistent
*      with inipk.f and ini3d.f
*      in arbitrary units 1 ul=10.98 Mpc
       RE0=1.0/10.98
       ROI=RODO*(1.0+ZI)**3
       PI4ROD=4.D0*PI*ROI
       REI=RE0/(1.0+ZI)
       TEI=T0*(1.0+ZI)**(-1.5)   !TEI=INITIAL TIME

       UV=299792.458
       UM=9.1717E+18
       F2=LOG(3.0)-(2.0/3.0)

       HTE=ACHE !!!!!!!!!!!! DEBE PONER BIEN

********************************************************************
********************************************************************

       NFILE2=INT((LAST-FIRST)/EVERY) + 1
       NFILE=2    !SIEMPRE SON 2 ITERACIONES

       WRITE(*,*)'NFILE2=',NFILE2
       WRITE(*,*)'NFILE=',NFILE

!$OMP PARALLEL DO SHARED(NNCLUS,ZETAS,TIEMPO,NFILE,DMPITER),PRIVATE(I)
       DO I=1,NFILE
        NNCLUS(I)=0
        ZETAS(I)=0.0
        TIEMPO(I)=0.0
        DMPITER(I)=0
       END DO

       NMAXNCLUSBAS=NMAXNCLUS
!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,NMASA,NR,CONCENTRA,
!$OMP+                   ANGULARM,NCLUSRX,NCLUSRY,NCLUSRZ,
!$OMP+                   VCM2,VMAXCLUS,VX,VY,VZ),PRIVATE(I)
       DO I=1,NMAXNCLUSBAS
        NMASA(:,I)=0.0
        NR(:,I)=0.0
        CONCENTRA(:,I)=0.0
        ANGULARM(:,I)=0.0
        NCLUSRX(:,I)=0.0
        NCLUSRY(:,I)=0.0
        NCLUSRZ(:,I)=0.0
        VCM2(:,I)=0.0
        VMAXCLUS(:,I)=0.0

        VX(:,I)=0.0
        VY(:,I)=0.0
        VZ(:,I)=0.0
       END DO

!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,VCMAX,MCMAX,RCMAX,DMPCLUS,
!$OMP+                M200,R200,IPLIP,IPLIR,REALCLUS,NLEVHAL,
!$OMP+                EIGENVAL),PRIVATE(I)
       DO I=1,NMAXNCLUSBAS
        VCMAX(:,I)=0.0
        MCMAX(:,I)=0.0
        RCMAX(:,I)=0.0
        M200(:,I)=0.0
        R200(:,I)=0.0
        IPLIP(:,I)=0
        IPLIR(:,I)=0
        DMPCLUS(:,I)=0
        REALCLUS(:,I)=0    !de momento no hay halos
        NLEVHAL(:,I)=0
        EIGENVAL(:,:,I)=0.0
       END DO

       IF (PLOT.GT.1) THEN   !MERGER TREE

        ALLOCATE (DMLIP(MAXITER,PARTIRED_PLOT))
        ALLOCATE (DMLIR(MAXITER,PARTIRED_PLOT))
        ALLOCATE (NEW_MASAP(MAXITER,PARTIRED_PLOT))
        ALLOCATE (NHOST(MAXITER,0:ININL,PARTIRED))
        ALLOCATE (DAD(MAXITER,NMAXDAD,NMAXCLUS_PLOT))
        ALLOCATE (NDAD(MAXITER,NMAXNCLUS))
        ALLOCATE (HHOST(MAXITER,NMAXHOST,0:ININL,PARTIRED))

        NBASPART_PLOT=PARTIRED_PLOT
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,DMLIP,DMLIR,
!$OMP+             NEW_MASAP),PRIVATE(I)
        DO I=1,NBASPART_PLOT
         DMLIP(:,I)=0
         DMLIR(:,I)=0
         NEW_MASAP(:,I)=0.0
        END DO

        NBASPART_PLOT=PARTIRED
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,NHOST,HHOST),
!$OMP+            PRIVATE(I)
        DO I=1,NBASPART_PLOT
         NHOST(:,:,I)=0
         HHOST(:,:,:,I)=0
        END DO

        NBASPART_PLOT=NMAXCLUS_PLOT
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,DAD,NDAD),
!$OMP+            PRIVATE(I)
        DO I=1,NBASPART_PLOT
         NDAD(:,I)=0
         DAD(:,:,I)=0
        END DO

        IF (PLOT.EQ.2) THEN
         ALLOCATE(RATIO(MAXITER,NMAXDAD,NMAXCLUS_PLOT))
         NBASPART_PLOT=NMAXCLUS_PLOT
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,RATIO),
!$OMP+            PRIVATE(I)
         DO I=1,NBASPART_PLOT
          RATIO(:,:,I)=0.0
         END DO
        END IF  !PLOT=2
       END IF  !PLOT>1

       CONTAITER=0
       MARK(1:NFILE2)=0

       INERTIA=0.0

*///////// MAIN LOOP (ITERATIONS) /////////
*//////////////////////////////////////////
       DO IFI2=1, NFILE2                 !/
*//////////////////////////////////////////
*//////////////////////////////////////////

        CONTAITER=CONTAITER+1
        ITER=FIRST+EVERY*(IFI2-1)

        IFI=CONTAITER

        IF (MARK(IFI2).EQ.0) THEN

        MARK(IFI2)=1

        WRITE(*,*)'STARTING ITER', ITER, IFI2, IFI

        PATCHNX=0
        PATCHNY=0
        PATCHNZ=0
        PATCHX=0
        PATCHY=0
        PATCHZ=0
        PATCHRX=0.0
        PATCHRY=0.0
        PATCHRZ=0.0

        NPATCH=0
        PARE=0

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U1,U1G,UBAS1),PRIVATE(I,J,K)
        DO K=1,NZ
        DO J=1,NY
        DO I=1,NX
         U1(I,J,K)=-1.0        !valores minimos
         U1G(I,J,K)=-1.0
         UBAS1(I,J,K)=0.0
        END DO
        END DO
        END DO

**********
*     cambio especial para paralelizar!!
**********
        N1=NAMRX                !dimension max de todos los parches
        NPBAS=NPALEV            !numero total maximo de parches
        PABAS=PARTIRED
        NLEVBAS=NLEVELS         !numero maximo de niveles
        NMAXNCLUSBAS=MAXNCLUS   !num max de candidatos a halo
**********


!$OMP PARALLEL DO SHARED(NPBAS,N1,U11,U11G),
!$OMP+        PRIVATE(IX,JY,KZ,I)
        DO I=1,NPBAS
         DO KZ=1,N1
         DO JY=1,N1
         DO IX=1,N1
          U11(IX,JY,KZ,I)=-1.0
          U11G(IX,JY,KZ,I)=-1.0
         END DO
         END DO
         END DO
        END DO

        NHALLEV=0
        NPART=0

!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,MASA,RADIO,
!$OMP+                   CLUSRX,CLUSRY,CLUSRZ,LEVHAL,NSOLAP),
!$OMP+            PRIVATE(I)
        DO I=1,NMAXNCLUSBAS
         CLUSRX(I)=0.0
         CLUSRY(I)=0.0
         CLUSRZ(I)=0.0
         MASA(I)=0.0
         RADIO(I)=0.0
         LEVHAL(I)=0
         NSOLAP(I)=0
        END DO

        NCLUS=0
        CMX=0.0
        CMY=0.0
        CMZ=0.0
        MASA2=0.0
        ROTE=0.0
        RETE=0.0

        SUBHALOS=0
        SOLAPA=0

!$OMP PARALLEL DO SHARED(PABAS,U2DM,U3DM,U4DM,RXPA,RYPA,RZPA,
!$OMP+                   MASAP,ORIPA1,ORIPA2,LIP,LIR,CONTADM),
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

         LIP(I)=0         !si PARTI=PARTIRED!!!!!
         LIR(I)=0         !si PARTI=PARTIRED!!!!!
         CONTADM(I)=0     !si PARTI=PARTIRED!!!!!
        END DO

***************************************************
*     READING INPUT DATA
***************************************************

       IF (FLAG_SA.EQ.1) THEN
*       Reading MASCLET files directly
        CALL READ_MASCLET(VAR,ITER,NX,NY,NZ,NDXYZ,T,ZETA,NL,NPATCH,
     &            PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ,MAP,U2DM,U3DM,U4DM,MASAP,
     &            NPART,RXPA,RYPA,RZPA,N_DM)

        ! Background cosmology variables
        ROTE=RODO*(1.0+ZETA)**3
        RETE=RE0/(1.0+ZETA)
        ZETAS(IFI)=ZETA
        TIEMPO(IFI)=T

       ELSE
*       Reading external list of particles (either Masclet particles
*       or a general list of particles, depending on FLAG_MASCLET)
        IF (FLAG_MASCLET.EQ.1) THEN
         CALL READ_PARTICLES_MASCLET(ITER,NX,NY,NZ,T,ZETA,MAP,
     &                               U2DM,U3DM,U4DM,MASAP,RXPA,
     &                               RYPA,RZPA,N_DM)
        ELSE
         CALL READ_PARTICLES_GENERAL(ITER,NX,NY,NZ,T,ZETA,NL,MAP,
     &                               U2DM,U3DM,U4DM,MASAP,RXPA,
     &                               RYPA,RZPA,LADO0,N_GAS,N_DM,
     &                               N_PARTICLES)
        END IF

!      FIX THIS, REMOVE NPART (USELESS) FROM EVERYWHERE
       !NPART=0
       !NPART=N_DM
       ! for now, we will leave it like this for this to work temporarily
        NPART=0
        NPART(0)=N_PARTICLES

        ! Background cosmology variables
        ROTE=RODO*(1.0+ZETA)**3
        RETE=RE0/(1.0+ZETA)
        ZETAS(IFI)=ZETA
        TIEMPO(IFI)=T

        WRITE(*,*)'***********************'
        WRITE(*,*)'***** MESHRENOEF ******'
        WRITE(*,*)'***********************'

        IF (NL.GT.0) THEN
         WRITE(*,*)'==== Building the grid...', ITER, NL
         CALL CREATE_MESH(ITER,NX,NY,NZ,NL,NPATCH,PARE,PATCHNX,PATCHNY,
     &                    PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,
     &                    PATCHRZ,RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,
     &                    N_PARTICLES,N_DM,N_GAS,LADO0,T,ZETA)
         WRITE(*,*)'==== END building the grid...', ITER, NL
        END IF

        CALL INTERPOLATE_DENSITY(ITER,NX,NY,NZ,NL,NPATCH,PARE,
     &           PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &           PATCHRX,PATCHRY,PATCHRZ,RXPA,RYPA,RZPA,MASAP,
     &           N_PARTICLES,N_DM,N_GAS,LADO0,T,ZETA)

       WRITE(*,*)'***************************'
       WRITE(*,*)'***** END MESHRENOEF ******'
       WRITE(*,*)'***************************'

       END IF

****************************************************************
*      VIRIAL CONTRAST ! gets CONTRASTEC and OMEGAZ
*      (Bryan & Norman ApJ, 1998)
*      Delta_vir,c = 18*pi^2 + 82 x - 39 x^2; x=Omega_m(z)-1
****************************************************************
       CALL BRYAN_NORMAN_98(CONTRASTEC,OMEGAZ,OMEGA0,ZETA)

       WRITE(*,*) '************************************************'
       WRITE(*,*) '             "COSMOLOGICAL" PARAMETERS            '
       WRITE(*,*) '************************************************'
       WRITE(*,*) 'RETE=', RETE
       WRITE(*,*) 'ROTE=', ROTE
       WRITE(*,*) 'RODO,RE0,OMEGA0,OMEGAZ=', RODO,RE0,OMEGA0,OMEGAZ
       WRITE(*,*) 'Z=', ZETA
       WRITE(*,*) 'CONTRASTEC=',CONTRASTEC
       WRITE(*,*) '************************************************'

**************************************************************
*      Cleaning overlaps of patches
*      NOTE! we correct overlaps and not refinements because
*            we work within each level independentely
**************************************************************

       ! CONTA2: overlaps at level IR; (=1, keep), (=0, overlapped)
       DO IR=1,NL
        CALL VEINSGRID(IR,NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &                 PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &                 SOLAP,VECINO,NVECI)
       END DO

******* Compute CR0AMR *********************************************
c       CALL COMPUTE_CR0AMR(NL,NX,NY,NZ,NPATCH,PARE,PATCHNX,PATCHNY,
c     &                     PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
c     &                     PATCHRY,PATCHRZ,TEMP0,TEMP1,LADO0)
c
c       OPEN(99,FILE='output_files/density_asohf',STATUS='UNKNOWN',
c     &      FORM='UNFORMATTED')
c         write(99) (((u1(ix,jy,kz),ix=1,nx),jy=1,ny),kz=1,nz)
c         write(99) (((temp0(ix,jy,kz),ix=1,nx),jy=1,ny),kz=1,nz)
c         do i=1,sum(npatch(0:nl))
c          n1=patchnx(i)
c          n2=patchny(i)
c          n3=patchnz(i)
c          write(99) (((u11(ix,jy,kz,i),ix=1,n1),jy=1,n2),kz=1,n3)
c          write(99) (((temp1(ix,jy,kz,i),ix=1,n1),jy=1,n2),kz=1,n3)
c          write(99) (((conta2(ix,jy,kz,i),ix=1,n1),jy=1,n2),kz=1,n3)
c         end do
c       CLOSE(99)
*********************************************************************

       CALL CLEAN_OVERLAPS(NL,NPATCH,PATCHNX,PATCHNY,PATCHNZ,SOLAP,
     &                     U11)

**********************************************************************
******************************HALO FINDER*****************************
**********************************************************************

       WRITE(*,*) '************************************************'
       WRITE(*,*) '         HALO FINDING                           '
       WRITE(*,*) '************************************************'


**********************************************************
*      Looking for candidate haloes at the AMR levels
**********************************************************

       CALL HALOFIND_GRID(IFI,NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                    PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                    PATCHRY,PATCHRZ,PARE,NCLUS,MASA,RADIO,
     &                    CLUSRX,CLUSRY,CLUSRZ,REALCLUS,LEVHAL,
     &                    NSOLAP,SOLAPA,NHALLEV,BOUND,CONTRASTEC,RODO,
     &                    SOLAP,VECINO,NVECI)

       STOP

c!!! hacer paralelo
       CONTA2=1                !celdas de parches (1 valen, 0 no)
c!!!
       CONTAP(1:NPALEV)=1      !parches (1 valen, 0 no)
       ICEN=0
       CEL=0
       CEL2=0
       CEN=0
       CEN2=0
       VID=0

*      We proceed from top to bottom
       IF (NL.GE.1) THEN

       DO IR=NL, 1, -1    !!!!!!!!!!!!

       DXPA=0.0
       DYPA=0.0
       DZPA=0.0
       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))

       WRITE(*,*) 'RESOLUTION=',IR,DXPA

       REF=0.0
       REF=1.5*DXPA
       ESP=(LOG10(LADO*0.5)-LOG10(REF))/(NSHELL-1)

**************==estimacion==*****
       KONTA=0
!$OMP  PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,
!$OMP+                    U11,CONTRASTEC,LOW1,LOW2),
!$OMP+             PRIVATE(N1,N2,N3,I,IX,JY,KZ,UBAS1),
!$OMP+             REDUCTION(+:KONTA)
       DO I=LOW1,LOW2

        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)

        DO KZ=1,N3
        DO JY=1,N2
        DO IX=1,N1

        UBAS1(IX,JY,KZ)=0.0
        UBAS1(IX,JY,KZ)=U11(IX,JY,KZ,I)
        IF (UBAS1(IX,JY,KZ).GE.CONTRASTEC) KONTA=KONTA+1

        END DO
        END DO
        END DO

       END DO
       WRITE(*,*)'ESTIMATION:', IR,KONTA,CONTRASTEC
***************


*      Asi obtengo todos los vecinos de cada parche de IR
       VECINO=0
       NVECI=0
       CALL VEINSGRID(IR,NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &           PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &           CONTA2,VECINO,NVECI)

*      Empieza la busqueda de celdita en este nivel


!!!!!!OJO       CONTA2=1

       MAXIMO(1:NPATCH(IR))=0.0

!!!!hacer paralelo (I2)

       DO I2=1,NPATCH(IR)           !ID "local" de parche en IR
       I=SUM(NPATCH(0:IR-1)) + I2   !ID "absoluto" de parche en la simu.

        N1=0
        N2=0
        N3=0
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)

        MAXIMO(I2)=MAXVAL(U11(1:N1,1:N2,1:N3,I)*
     &                    CONTA2(1:N1,1:N2,1:N3,I))

        !SUSANA_QUESTION_1: YA NO ES NECESARIO *CONTA2, NO?
        ! creo que no pq U11 ya esta corregido mas arriba:
c!SUSANA        MAXIMO(I2)=MAXVAL(U11(1:N1,1:N2,1:N3,I))

        END DO
!!!!hacer paralelo

**////////**********
       KK_REAL=0.0
       KK_REAL=MAXVAL(MAXIMO(1:NPATCH(IR)))

C_VICENT       DO WHILE (KK_REAL.GE.CONTRASTEC)
       DO
       IF(KK_REAL.LT.CONTRASTEC) EXIT
C_VICENT

       CEN=MAXLOC(MAXIMO(1:NPATCH(IR)))    !numero "local" de parche en IR
       CEN(1)=SUM(NPATCH(0:IR-1))+CEN(1)   !numero "absoluto" de parche

        IF (CEN(1).GT.SUM(NPATCH(0:IR))) THEN
         WRITE(*,*) 'WARNING!!CEN1=',CEN(1),NPATCH(IR),SUM(NPATCH(0:IR))
         STOP
       END IF

       CEL=CEN(1)   !!numero "absoluto" de parche

**     TENGO QUE SABER TODA LA FAMILIA DE VECINOS DEL PARCHE ELEGIDO

       FLAG=0
       NV=0
       NV2=0
       NV3=0
       NV=1        !primer vecino: el mismo
       VID(NV)=CEL

       NV2=0
       NV3=NV

*       CONTAP(1:NPATCH(IR))=1
       CONTAP(LOW1:LOW2)=1    !(1 vale, 0 no)

       CONTAP(VID(NV))=0      !asi este no se puede autocontar

*****************8
C_VICENT       DO WHILE (FLAG.EQ.0)
       DO
       IF (FLAG.NE.0) EXIT
cC_VICENT
       DO II=NV2+1, NV3

       DO JJ=1, NVECI(VID(II))

       KK_ENTERO=0
       KK_ENTERO=CONTAP(VECINO(JJ,VID(II)))
       IF (KK_ENTERO.EQ.1) THEN

       NV=NV+1
       VID(NV)=VECINO(JJ,VID(II))
       CONTAP(VECINO(JJ,VID(II)))=0

       END IF

       END DO

       END DO

       IF (NV.EQ.NV3) THEN
        FLAG=1
       ELSE
        NV2=NV3
        NV3=NV
       END IF

       END DO  ! WHILE FLAG
******************8


*      VAMOS A LIQUIDAR TODA LA JERARQUIA DE VECINOS
*      voy a ordenar todas las celdas de los vecinos VECTORIZANDO!!
        KK_ENTERO=0
C_VICENT
!$OMP  PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,
!$OMP+                 U11,CONTRASTEC,NV,VID,CONTA2),
!$OMP+             PRIVATE(I,N1,N2,N3,KK_ENTERO_2),
!$OMP+             REDUCTION(+:KK_ENTERO)
C_VICENT
       DO I=1,NV
        N1=0
        N2=0
        N3=0
        N1=PATCHNX(VID(I))
        N2=PATCHNY(VID(I))
        N3=PATCHNZ(VID(I))

        KK_ENTERO_2=0
        KK_ENTERO_2=COUNT(U11(1:N1,1:N2,1:N3,VID(I))
     &              *CONTA2(1:N1,1:N2,1:N3,VID(I)).GE.CONTRASTEC)
c!SUSANA        KK_ENTERO_2=COUNT(U11(1:N1,1:N2,1:N3,VID(I)).GE.CONTRASTEC)
        KK_ENTERO=KK_ENTERO+KK_ENTERO_2

        !SUSANA_QUESTION_2: YA NO ES NECESARIO *CONTA2, NO?-->creo q no
        !SUSANA_QUESTION_2:no habria que poner *REAL(CONTA2)))
       END DO


       ALLOCATE(DDD(KK_ENTERO))
       ALLOCATE(DDDP(KK_ENTERO))
       ALLOCATE(DDDX(KK_ENTERO))
       ALLOCATE(DDDY(KK_ENTERO))
       ALLOCATE(DDDZ(KK_ENTERO))
       ALLOCATE(INDICE2(KK_ENTERO))
       ALLOCATE(DDD2(KK_ENTERO))
       ALLOCATE(DDDP2(KK_ENTERO))
       ALLOCATE(DDDX2(KK_ENTERO))
       ALLOCATE(DDDY2(KK_ENTERO))
       ALLOCATE(DDDZ2(KK_ENTERO))

C_VICENT
!$OMP  PARALLEL DO SHARED(KK_ENTERO,INDICE,DDD2,DDDP2,
!$OMP+       DDD,DDDP,DDDX,DDDY,DDDZ,DDDX2,DDDY2,DDDZ2),
!$OMP+       PRIVATE(I)
C_VICENT
        DO I=1, KK_ENTERO
         DDD(I)=0.0
         DDDP(I)=0
         DDDX(I)=0
         DDDY(I)=0
         DDDZ(I)=0
         INDICE(I)=0
         DDD2(I)=0.0
         DDDP2(I)=0
         DDDX2(I)=0
         DDDY2(I)=0
         DDDZ2(I)=0
       END DO

       II=0
       DO I=1,NV
        N1=PATCHNX(VID(I))
        N2=PATCHNY(VID(I))
        N3=PATCHNZ(VID(I))

        DO KZ=1, N3
        DO JY=1, N2
        DO IX=1, N1

        KK_REAL=0.0
        KK_REAL=U11(IX,JY,KZ,VID(I))*CONTA2(IX,JY,KZ,VID(I))
c!SUSANA        KK_REAL=U11(IX,JY,KZ,VID(I))

        IF(KK_REAL.GE.CONTRASTEC) THEN
         II=II+1
         DDD(II)=1.0/KK_REAL
         DDDP2(II)=VID(I)
         DDDX2(II)=IX
         DDDY2(II)=JY
         DDDZ2(II)=KZ
        END IF

        END DO
        END DO
        END DO
       END DO

       IF (II.NE.KK_ENTERO) THEN
         WRITE(*,*) 'WARNING,NV', II, KK_ENTERO
         STOP
       END IF

       CALL INDEXX(KK_ENTERO,DDD(1:KK_ENTERO),INDICE2(1:KK_ENTERO))
C_VICENT
!$OMP  PARALLEL DO SHARED(KK_ENTERO,INDICE2,DDD,DDD2),
!$OMP+       PRIVATE(I)
C_VICENT
       DO I=1,KK_ENTERO
        DDD2(I)=1.0/DDD(INDICE2(I))
       END DO
       !ahora ya estan ordenadas todas las celdas de todos los parches vecinos

!$OMP  PARALLEL DO SHARED(KK_ENTERO,INDICE2,DDD2,DDDP2,
!$OMP+       DDD,DDDP,DDDX,DDDY,DDDZ,DDDX2,DDDY2,DDDZ2),
!$OMP+       PRIVATE(I)
       DO I=1, KK_ENTERO
        DDD(I)=0.0
        DDDP(I)=0
        DDDX(I)=0
        DDDY(I)=0
        DDDZ(I)=0
        DDD(I)=DDD2(I)
        DDDP(I)=DDDP2(INDICE2(I))
        DDDX(I)=DDDX2(INDICE2(I))
        DDDY(I)=DDDY2(INDICE2(I))
        DDDZ(I)=DDDZ2(INDICE2(I))
       END DO

       NV_GOOD=0     !num. celdas parches vecinos ordenadas de > a <
       NV_GOOD=KK_ENTERO
       KK_ENTERO=0


       DO L1=1,NV_GOOD

*      de todos los vecinos de CEL este es el parche y celda elegidos!!

       CEL2=DDDP(L1)   !me da el parche correspondiente al vecino
       ICEN(1)=DDDX(L1)
       ICEN(2)=DDDY(L1)
       ICEN(3)=DDDZ(L1)

       KK_ENTERO=0
       KK_ENTERO=CONTA2(ICEN(1),ICEN(2),ICEN(3),CEL2)
       IF (KK_ENTERO.EQ.1) THEN

*      COORDENADAS DE ESTA CELDA DEL PARCHE
       RX2=0.0
       RY2=0.0
       RZ2=0.0
       RX2=PATCHRX(CEL2) - 0.5*DXPA + (ICEN(1)-1)*DXPA
       RY2=PATCHRY(CEL2) - 0.5*DYPA + (ICEN(2)-1)*DYPA
       RZ2=PATCHRZ(CEL2) - 0.5*DZPA + (ICEN(3)-1)*DZPA

       NCLUS=NCLUS+1
       LEVHAL(NCLUS)=IR
       REALCLUS(IFI,NCLUS)=-1

       NHALLEV(IR)=NHALLEV(IR)+1
       CLUSRX(NCLUS)=RX2
       CLUSRY(NCLUS)=RY2
       CLUSRZ(NCLUS)=RZ2

       PRUEBAX=0.0
       PRUEBAY=0.0
       PRUEBAZ=0.0
       DELTA_PRUEBA=0.0

*      EXTENDEMOS RADIALMENTE EL CLUSTER NCLUS
       RMIN=0.0
       DELTA=0.0
       SHELL=0
       BASMAS=0.0   !MASA TOTAL DEL CUMULO
       DELTA2=0.0   !DELTA TOTAL DEL CUMULO
       VOL2=0.0     !VOL TOTAL DEL CUMULO

       DELTA2=10.0*CONTRASTEC*ROTE

       DO WHILE(DELTA2.GT.(CONTRASTEC)*ROTE)

       IF(SHELL<NSHELL) THEN

          SHELL=SHELL+1
          DELTA=0.0
          RSHELL=0.0
          VOL=0.0

       ELSE
          WRITE(*,*) 'WARNING: DEMASIADAS CAPAS'
          WRITE(*,*) SHELL,DELTA,DELTA2,(CONTRASTEC)*ROTE,RSHELL
C          STOP
       END IF

       RRRR=0.0
       R111=0.0
       R222=0.0
       RRRR=LOG10(REF)+(SHELL-1)*ESP
       IF (SHELL.EQ.1) THEN
         R111=RRRR+0.5*ESP
         R222=0.0
         R111=10.0**R111
       ELSE
        R111=RRRR+0.5*ESP
        R222=RRRR-0.5*ESP
        R111=10.0**R111
        R222=10.0**R222
       END IF
       RSHELL=R111

       VOL=PI*(4.0/3.0)*((R111*RETE)**3-(R222*RETE)**3)

*      EXTENDEMOS buscando en todas las celdas de todos los vecinos...
C_VICENT ----PARALELIZAR EN II!!!!
!$OMP  PARALLEL DO SHARED(NV,PATCHRX,PATCHRY,PATCHRZ,VID,
!$OMP+       CONTA2,DXPA,DYPA,DZPA,CLUSRX,CLUSRY,CLUSRZ,NCLUS,
!$OMP+       R222,R111,U11),
!$OMP+       PRIVATE(II,N1,N2,N3,IX,JY,KZ,KK_ENTERO,
!$OMP+       RX2,RY2,RZ2,AA),
!$OMP+       REDUCTION(+:DELTA,PRUEBAX,PRUEBAY,PRUEBAZ,
!$OMP+                  DELTA_PRUEBA)
C_VICENT
       DO II=1,NV
        N1=0
        N2=0
        N3=0
        N1=PATCHNX(VID(II))
        N2=PATCHNY(VID(II))
        N3=PATCHNZ(VID(II))

       DO KZ=1, N3
       DO JY=1, N2
       DO IX=1, N1

       KK_ENTERO=0
       KK_ENTERO=CONTA2(IX,JY,KZ,VID(II))
       IF (KK_ENTERO.EQ.1) THEN

         RX2=0.0
         RY2=0.0
         RZ2=0.0
         RX2=PATCHRX(VID(II)) - 0.5*DXPA + (IX-1)*DXPA
         RY2=PATCHRY(VID(II)) - 0.5*DYPA + (JY-1)*DYPA
         RZ2=PATCHRZ(VID(II)) - 0.5*DZPA + (KZ-1)*DZPA

         AA=0.0
         AA=SQRT((CLUSRX(NCLUS)-RX2)**2+(CLUSRY(NCLUS)-RY2)**2
     &          +(CLUSRZ(NCLUS)-RZ2)**2)

         IF(AA.GE.R222.AND.AA.LE.R111) THEN

           DELTA=DELTA+U11(IX,JY,KZ,VID(II))*DXPA*DYPA*DZPA
           CONTA2(IX,JY,KZ,VID(II))=0

           PRUEBAX=PRUEBAX+RX2*U11(IX,JY,KZ,VID(II))
           PRUEBAY=PRUEBAY+RY2*U11(IX,JY,KZ,VID(II))
           PRUEBAZ=PRUEBAZ+RZ2*U11(IX,JY,KZ,VID(II))
           DELTA_PRUEBA=DELTA_PRUEBA+U11(IX,JY,KZ,VID(II))

         END IF

       END IF

       END DO
       END DO
       END DO

       END DO  ! NV

       BASMAS=BASMAS+DELTA*RODO*RE0**3
       VOL2=PI*(4.0/3.0)*(RSHELL*RETE)**3
       DELTA2=BASMAS/VOL2
       DELTA=(DELTA*RODO*RE0**3)/VOL
       RMIN=RSHELL

       END DO  !!EN DO WHILE DE CRECER EL HALO (DELTA2)

       MASA(NCLUS)=BASMAS
       RADIO(NCLUS)=RSHELL

       PRUEBAX=PRUEBAX/DELTA_PRUEBA
       PRUEBAY=PRUEBAY/DELTA_PRUEBA
       PRUEBAZ=PRUEBAZ/DELTA_PRUEBA

       CLUSRX(NCLUS)=PRUEBAX
       CLUSRY(NCLUS)=PRUEBAY
       CLUSRZ(NCLUS)=PRUEBAZ

       KK_REAL=0.0
       KK_REAL=MASA(NCLUS)
       IF (KK_REAL.LE.0.0) REALCLUS(IFI,NCLUS)=0

       END IF

       END DO  !NV_GOOD

       DEALLOCATE(DDD)
       DEALLOCATE(DDDP)
       DEALLOCATE(DDDX)
       DEALLOCATE(DDDY)
       DEALLOCATE(DDDZ)
       DEALLOCATE(INDICE2)
       DEALLOCATE(DDD2)
       DEALLOCATE(DDDP2)
       DEALLOCATE(DDDX2)
       DEALLOCATE(DDDY2)
       DEALLOCATE(DDDZ2)
*      Fin VECTORIZACION!!!

*       DO I=LOW1, LOW2
C_VICENT
!$OMP  PARALLEL DO SHARED(IR,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
!$OMP+              MAXIMO,U11),
!$OMP+       PRIVATE(I2,I,N1,N2,N3)
C_VICNET
        DO I2=1,NPATCH(IR)
         I=SUM(NPATCH(0:IR-1)) + I2

         N1=0
         N2=0
         N3=0
         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

*        MAXIMO(I)=MAXVAL(U11(1:N1,1:N2,1:N3,I)
*     &           *CONTA2(1:N1,1:N2,1:N3,I))
         MAXIMO(I2)=MAXVAL(U11(1:N1,1:N2,1:N3,I)
     &              *CONTA2(1:N1,1:N2,1:N3,I))
c!SUSANA         MAXIMO(I2)=MAXVAL(U11(1:N1,1:N2,1:N3,I))
       END DO

       KK_REAL=0.0
       KK_REAL=MAXVAL(MAXIMO(1:NPATCH(IR)))

       END DO !WHILE MAXVAL(MAXIMO(IR)).GT.CONTRASTEC
**////////**********


         WRITE(*,*) '==== First Estimate ==='
         WRITE(*,*)  IFI,IR,REF,ESP,NCLUS
         KONTA2=0
         KONTA2=COUNT(REALCLUS(IFI,1:NCLUS).EQ.-1)
         WRITE(*,*) 'POSSIBLE HALOS_0----->', KONTA2
         KONTA2=0
         KONTA2=COUNT(REALCLUS(IFI,1:NCLUS).EQ.0)
         WRITE(*,*) 'REMOVED HALOS_0----->', KONTA2


****************************************************
*      CORRECION DE SOLAPES:
*      VAMOS A VER QUE CUMULOS SOLAPAN EN IR
****************************************************

       CALL OVERLAPING(IFI,IR,NL,REF,ESP,BOUND,CONTA,CONTRASTEC,
     &                 NSHELL,RODO,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                 PATCHRX,PATCHRY,PATCHRZ,NX,NY,NZ,
     &                 NCLUS,MASA,RADIO,CLUSRX,CLUSRY,CLUSRZ,
     &                 REALCLUS,NSOLAP,SOLAPA,NHALLEV)

       WRITE(*,*) '==== After correcting overlaps ==='
       WRITE(*,*)  IFI,IR,REF,ESP,NCLUS
       KONTA2=0
       KONTA2=COUNT(REALCLUS(IFI,1:NCLUS).EQ.-1)
       WRITE(*,*) 'POSSIBLE HALOS_0----->', KONTA2
       KONTA2=0
       KONTA2=COUNT(REALCLUS(IFI,1:NCLUS).EQ.0)
       WRITE(*,*) 'REMOVED HALOS_0----->', KONTA2

****************************************************
*      FIN CORRECION DE SOLAPES EN IR
****************************************************

       WRITE(*,*)' FIN NIVEL ', IR
C       CALL ITIME(TIME)
C       WRITE(*,*) 'TIME=',TIME(1),':',TIME(2),':',TIME(3)

       END DO  !IR

       END IF !if HAY NIVELES






*************************************
*      CHECKING HALOES PER LEVEL
*************************************
       WRITE(*,*)'CHECKING HALOES PER LEVEL 1', NCLUS
       DO I=0, NL
        WRITE(*,*) I, NHALLEV(I)
       END DO
*************************************

*********************************************************
*      CHECKING....
*********************************************************
       WRITE(*,*) '---------------------------------'
       WRITE(*,*) 'CHECKING HIERARCHY_0:'
       WRITE(*,*) '---------------------------------'
       KONTA2=0
       KONTA2=COUNT(REALCLUS(IFI,1:NCLUS).EQ.-1)
       WRITE(*,*) 'POSSIBLE HALOS_0----->', KONTA2
       KONTA2=0
       KONTA2=COUNT(REALCLUS(IFI,1:NCLUS).EQ.0)
       WRITE(*,*) 'REMOVED HALOS_0----->', KONTA2
*******************************************************


*******************************************************
*      SORTING OUT ALL THE CLUSTERS
*******************************************************

       AA1=0.0
       DO I=1, NCLUS

       NNCLUS(IFI)=NNCLUS(IFI)+ABS(REALCLUS(IFI,I))

       KK_ENTERO=0
       KK_ENTERO=REALCLUS(IFI,I)
       IF(KK_ENTERO.EQ.-1) THEN

       KK_REAL=0.0
       KK_REAL=MASA(I)
       IF (KK_REAL.LE.0.0) WRITE(*,*) 'WAR:', I,MASA(I)*9.1717E18

          NCLUSRX(IFI,NNCLUS(IFI))=CLUSRX(I)
          NCLUSRY(IFI,NNCLUS(IFI))=CLUSRY(I)
          NCLUSRZ(IFI,NNCLUS(IFI))=CLUSRZ(I)
          NMASA(IFI,NNCLUS(IFI))=MASA(I)*9.1717E18
          NR(IFI,NNCLUS(IFI))=RADIO(I)
          AA1=AA1+NMASA(IFI,NNCLUS(IFI))
          NLEVHAL(IFI,NNCLUS(IFI))=LEVHAL(I)

       END IF

       END DO

       REALCLUS(IFI,1:NNCLUS(IFI))=-1

       WRITE(*,*)'MASSES MAX & MIN=', MAXVAL(NMASA(IFI,1:NNCLUS(IFI))),
     &            MINVAL(NMASA(IFI,1:NNCLUS(IFI)))

cv7       WRITE(*,*) 'ESTIMATION HALOES_0 (ONLY WITH GRID)'
cv7       WRITE(*,*)'=================================='
cv7       DO I=1, NNCLUS(IFI)
cv7CX       KK_ENTERO=0
cv7CX       KK_ENTERO=REALCLUS(IFI,I)
cv7CX       IF (KK_ENTERO.NE.0) THEN
cv7       WRITE(*,*) I, NMASA(IFI,I),NR(IFI,I),NCLUSRX(IFI,I),
cv7     &            NCLUSRY(IFI,I),NCLUSRZ(IFI,I),REALCLUS(IFI,I)
cv7CX       END IF
cv7       END DO
cv7       WRITE(*,*)'=================================='


       WRITE(*,*) 'NNCLUS=', NNCLUS(IFI)
       AA1=AA1/NNCLUS(IFI)
       WRITE(*,*)'MEAN MASS',AA1


*****
*******************************************************
**     HALOES AT THE EDGES OF THE BOX (JUST FOR CAUTION!)

       IF (BORDES.EQ.1) THEN

       KK_ENTERO=0
!!!!!paralelizar
       DO KK1=1, NNCLUS(IFI)

       IF((NCLUSRX(IFI,KK1)-NR(IFI,KK1)).GT.(LADO*0.5).OR.
     &    (NCLUSRX(IFI,KK1)-NR(IFI,KK1)).LT. (-LADO*0.5).OR.
     &    (NCLUSRX(IFI,KK1)+NR(IFI,KK1)).GT.(LADO*0.5).OR.
     &    (NCLUSRX(IFI,KK1)+NR(IFI,KK1)).LT. (-LADO*0.5).OR.
     &    (NCLUSRY(IFI,KK1)-NR(IFI,KK1)).GT.(LADO*0.5).OR.
     &    (NCLUSRY(IFI,KK1)-NR(IFI,KK1)).LT. (-LADO*0.5).OR.
     &    (NCLUSRY(IFI,KK1)+NR(IFI,KK1)).GT.(LADO*0.5).OR.
     &    (NCLUSRY(IFI,KK1)+NR(IFI,KK1)).LT. (-LADO*0.5).OR.
     &    (NCLUSRZ(IFI,KK1)-NR(IFI,KK1)).GT.(LADO*0.5).OR.
     &    (NCLUSRZ(IFI,KK1)-NR(IFI,KK1)).LT. (-LADO*0.5).OR.
     &    (NCLUSRZ(IFI,KK1)+NR(IFI,KK1)).GT.(LADO*0.5).OR.
     &    (NCLUSRZ(IFI,KK1)+NR(IFI,KK1)).LT. (-LADO*0.5)) THEN

       KK_ENTERO=KK_ENTERO+1
       END IF

       END DO

       WRITE(*,*) 'Haloes at the limits of the box?', KK_ENTERO

       END IF

************************************************************
**     COMPUTING CM AND PARTICLES FOR ALL THE HALOES
**     (We start here to work with partciles for the 1st time)
************************************************************

       WRITE(*,*) '============================================='
       WRITE(*,*) 'COMPUTING CM AND PARTICLES FOR ALL THE HALOES'
       WRITE(*,*) '============================================='


!$OMP  PARALLEL DO SHARED(NNCLUS,IFI,NLEVHAL,NCLUSRX,NCLUSRY,NCLUSRZ,
!$OMP+                   RXPA,RYPA,RZPA,MASAP,NL,NPART,NMASA,DMPCLUS,
!$OMP+                   VCM2,NR,U2DM,U3DM,U4DM),
!$OMP+   PRIVATE(I,MASADM,KONTA,BASMAS,MASAKK,VCMX,VCMY,VCMZ,VCM,
!$OMP+           REF_MIN,REF_MAX,LIP,LIR,BAS,IR,J,AADM,
!$OMP+           KK_REAL,KK1,KK2,CONTADM,CMX,CMY,CMZ,MASA2)
*****************************
       DO I=1, NNCLUS(IFI)
****************************

       MASADM=0.0
       KONTA=0
       BASMAS=0.0
       MASAKK=0.0

       VCMX=0.0
       VCMY=0.0
       VCMZ=0.0
       VCM=0.0
       MASA2=0.0

       REF_MIN=1.0E+10
       REF_MAX=-1.0

       LIP=0
       LIR=0

*****susana
C       IF (NLEVHAL(IFI,I).EQ.0) BAS=2.0
C       IF (NLEVHAL(IFI,I).EQ.1) BAS=1.25
C       IF (NLEVHAL(IFI,I).EQ.2) BAS=1.2
C       IF (NLEVHAL(IFI,I).GT.2) BAS=1.1
        BAS=1.1
*****

*      NIVEL BASE
       IR=0

       DO J=1, NPART(0)

       AADM=0.0
       AADM=SQRT((RXPA(J)-NCLUSRX(IFI,I))**2+
     &           (RYPA(J)-NCLUSRY(IFI,I))**2+
     &           (RZPA(J)-NCLUSRZ(IFI,I))**2)

       KK_REAL=0.0
       KK_REAL=NR(IFI,I)
       IF(AADM.LT.BAS*KK_REAL) THEN

        REF_MIN=MIN(REF_MIN,AADM)
        REF_MAX=MAX(REF_MAX,AADM)

        KONTA=KONTA+1
        MASADM=MASADM+MASAP(J)
        LIP(KONTA)=J
        LIR(KONTA)=IR

       END IF

       END DO


*      NIVELES AMR
       DO IR=1, NL

       KK1=0
       KK2=0
       KK1=SUM(NPART(0:IR-1))
       KK2=SUM(NPART(0:IR))

       DO J=KK1+1, KK2

       AADM=0.0
       AADM=SQRT((RXPA(J)-NCLUSRX(IFI,I))**2+
     &           (RYPA(J)-NCLUSRY(IFI,I))**2+
     &           (RZPA(J)-NCLUSRZ(IFI,I))**2)

       KK_REAL=0.0
       KK_REAL=NR(IFI,I)
       IF(AADM.LT.BAS*KK_REAL) THEN

        REF_MIN=MIN(REF_MIN,AADM)
        REF_MAX=MAX(REF_MAX,AADM)

        KONTA=KONTA+1
        MASADM=MASADM+MASAP(J)
        LIP(KONTA)=J
        LIR(KONTA)=IR

       END IF

       END DO
       END DO

       CONTADM(1:KONTA)=0

*      CHECKING BAD MASSES
       IF (MASADM.GT.0.0) THEN

         CALL CENTROMASAS_PART(KONTA,CONTADM,LIP,
     &          U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &          CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MASA2)

       END IF

*****////OJO: solo me quedo con las veloc. pq las coordenadas
***** las voy a refinar solo con la malla en PASO 2 si
***** fuera necesario

       VCM=SQRT(VCMX**2+VCMY**2+VCMZ**2)
       VCM2(IFI,I)=VCM
       MASAKK=MASADM*9.1717E18
       NMASA(IFI,I)=MASAKK
       DMPCLUS(IFI,I)=KONTA  !esto solo da un num inicial de part. a cada halo

*****************************
       END DO        !I
*****************************
       WRITE(*,*)'END COMPUTING CM AND PARTICLES FOR ALL THE HALOES'
       WRITE(*,*) '==================================================='

**     FIN VCM



************************************************
*     CLEANING THE SAMPLE OF HALOES
************************************************

*STEP 1) (1)************************************
************REMOVING POOR HALOES****************
************************************************

       KONTA2=0
       NUMPARTBAS=NUMPART
!$OMP  PARALLEL DO SHARED(NNCLUS,IFI,DMPCLUS,NUMPARTBAS,
!$OMP+                    REALCLUS),
!$OMP+             PRIVATE(I,KK_ENTERO)
!$OMP+             REDUCTION(+:KONTA2)
       DO I=1, NNCLUS(IFI)

        KK_ENTERO=0
        KK_ENTERO=DMPCLUS(IFI,I)
        IF (KK_ENTERO.LT.NUMPARTBAS) THEN
cx_test10        IF (KK_ENTERO.LT.1) THEN

        REALCLUS(IFI,I)=0
        KONTA2=KONTA2+1

       END IF

       END DO

       WRITE(*,*)'CHECKING POOR HALOS----->', KONTA2


*********************************************************
*      CHECKING....
*********************************************************
       WRITE(*,*) '---------------------------------'
       WRITE(*,*) 'CHECKING HIERARCHY:'
       WRITE(*,*) '---------------------------------'
       KONTA2=0
       KONTA2=COUNT(REALCLUS(IFI,1:NNCLUS(IFI)).EQ.-1)
       WRITE(*,*) 'REAL HALOS----->', KONTA2
       KONTA2=0
       KONTA2=COUNT(REALCLUS(IFI,1:NNCLUS(IFI)).EQ.0)
       WRITE(*,*) 'REMOVED HALOS----->', KONTA2
       KONTA2=0
       KONTA2=COUNT(REALCLUS(IFI,1:NNCLUS(IFI)).GT.0)
       WRITE(*,*) 'SUBSTRUCTURE----->', KONTA2
       WRITE(*,*) '---------------------------------'
*********************************************************


************************************************************
*      REFINING REAL HALOES WITH THE DM PARTICLES ONLY     *
************************************************************

       DIMEN=3   !DIMENSION DE LOS HALOS

       NCAPAS=0
       ESP=0.0
       REF=0.0

       WRITE(*,*)'Refining with DM particles...'
       WRITE(*,*)'=================================='

       DO I=0,NL
       WRITE(*,*)'Halos at level ', I,' =',
     &            COUNT(NLEVHAL(IFI,1:NNCLUS(IFI)).EQ.I),
     &            COUNT(REALCLUS(IFI,1:NNCLUS(IFI)).NE.0.AND.
     &                       NLEVHAL(IFI,1:NNCLUS(IFI)).EQ.I)
       END DO
       WRITE(*,*)'=================================='

CX
       KONTA2=0
       KONTA2=MAXVAL(DMPCLUS(IFI,1:NNCLUS(IFI)))
       WRITE(*,*) 'MAXVAL(DMPCLUS(IFI,1:NNCLUS(IFI)))=', KONTA2
       WRITE(*,*) 'MINVAL(DMPCLUS(IFI,1:NNCLUS(IFI)))=',
     &             MINVAL(DMPCLUS(IFI,1:NNCLUS(IFI)))
       !OJO: Debido a los valores de BAS anteriores y posteriores,
       !     a continuación pueden haber más particulas que antes!
       KONTA2=KONTA2+1000000 !!!!!PRUEBA!!!


       KONTA1=0
       KONTA1=NNCLUS(IFI)
       WRITE(*,*) 'NNCLUS(IFI)=', NNCLUS(IFI)

c_v6: .............
       NN=NUM*3      !NUMERO DE HALOS POR BUCLE (3/procesador cada vez)
       NH=NNCLUS(IFI)
       IP=0
       IP2=0
       WRITE(*,*)'NUMBER HALOS/PROCESSOR=', NN


        IF(PLOT.GT.1) THEN
         ALLOCATE(IP_PARAL(KONTA2,NN))   !CON KONTA1 ENTRE 1 Y NN
         ALLOCATE(IR_PARAL(KONTA2,NN))
         ALLOCATE(MASAP_PARAL(KONTA2,NN))
        END IF

       !+++++!
       DO
       !+++++!

       IP=IP+1
       IP2=IP+NN-1
       IF (IP2.GT.NH) IP2=NH

       IF (IP.GT.NH) EXIT
c_v6: .............

      IF(PLOT.GT.1) THEN
        NBASPART_PLOT=NN
!$OMP PARALLEL DO SHARED(NBASPART_PLOT,IP_PARAL,
!$OMP+                   IR_PARAL,MASAP_PARAL),
!$OMP+            PRIVATE(I)
        DO I=1,NBASPART_PLOT
           IP_PARAL(:,I)=0
           IR_PARAL(:,I)=0
           MASAP_PARAL(:,I)=0.0
        END DO
      END IF


       KONTA1=0
       KONTA2=0
       PABAS=PARTIRED_PLOT
       NUMPARTBAS=NUMPART
       !!!OJO: a continuacion cambio OMEGA0 por OMEGAZ

!$OMP  PARALLEL DO SHARED(NNCLUS,IFI,REALCLUS,
!$OMP+           NLEVHAL,NPART,RXPA,RYPA,RZPA,NCLUSRX,NCLUSRY,NCLUSRZ,
!$OMP+           NL,NR,MASAP,U2DM,U3DM,U4DM,VCM2,VX,VY,VZ,ACHE,
!$OMP+           PI,RETE,ROTE,VCMAX,MCMAX,RCMAX,M200,R200,CONTRASTEC,
!$OMP+           OMEGAZ,CGR,UM,UV,DMPCLUS,CONCENTRA,PLOT,
!$OMP+           ORIPA1,ORIPA2,ANGULARM,PABAS,IPLIP,IPLIR,DIMEN,
!$OMP+           EIGENVAL,NUMPARTBAS,IP_PARAL,IR_PARAL,MASAP_PARAL,
!$OMP+           IP,IP2,NH),
!$OMP+   PRIVATE(I,INERTIA,REF_MIN,REF_MAX,KK_ENTERO,MASADM,KONTA,
!$OMP+           BASMAS,MASAKK,DIS,VCM,VVV2,VR,LIP,LIR,CONCEN,RS,FC,
!$OMP+           VMAX2,KONTA2,BAS,IR,J,AADM,KK1,KK2,CONTADM,
!$OMP+           CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MASA2,DISTA,FAC,CONTAERR,
!$OMP+           NORMA,JJ,DENSITOT,RADIAL,SALIDA,BAS1,BAS2,FLAG_200,
!$OMP+           VOL,DELTA2,NCAPAS,RSHELL,KONTA3,NSHELL_2,KONTA1,
!$OMP+           DENSA,DENSB,DENSC,AADMX,VKK,AA,NROT,BASEIGENVAL)

*****************************
c_v6:       DO I=1, NNCLUS(IFI)
       DO I=IP, IP2
****************************

       KONTA1=0
       KONTA1=DMPCLUS(IFI,I)
       DMPCLUS(IFI,I)=0

       KK_ENTERO=0
       KK_ENTERO=REALCLUS(IFI,I)

       IF (KK_ENTERO.NE.0) THEN

       INERTIA=0.0
       REF_MIN=10.0e+10
       REF_MAX=-1.0

       MASADM=0.0
       KONTA=0
       BASMAS=0.0
       MASAKK=0.0

       DIS=1.0E+10    !DIS. DE LA PARTI MAS CERCANA AL CENTRO DEL HALO

       VCM=0.0
       VVV2=0.0    ! v propia de las particulas respecto al CM
       VR=0.0      ! v radial
       CMX=0.0     !just added:
       CMY=0.0
       CMZ=0.0
       VCMX=0.0
       VCMY=0.0
       VCMZ=0.0

       LIP=0
       LIR=0

       CONCEN=0.0       !concentracion perfil NFW
       RS=0.0           !radio de escala perfil NFW
       FC=0.0
       VMAX2=0.0
       KONTA2=0


**     COMPUTING VCM OF HALO I

*       BAS=1.0 + 2./(2.**(NLEVHAL(IFI,I)+1))
       IF (NLEVHAL(IFI,I).EQ.0) BAS=2.0
       IF (NLEVHAL(IFI,I).EQ.1) BAS=1.25
       IF (NLEVHAL(IFI,I).EQ.2) BAS=1.2
       IF (NLEVHAL(IFI,I).GT.2) BAS=1.1
*       BAS=1.2

*      NIVEL BASE
       IR=0

       DO J=1, NPART(0)
        AADM=SQRT((RXPA(J)-NCLUSRX(IFI,I))**2+
     &            (RYPA(J)-NCLUSRY(IFI,I))**2+
     &            (RZPA(J)-NCLUSRZ(IFI,I))**2)

        IF(AADM.LT.BAS*NR(IFI,I)) THEN
         REF_MIN=MIN(REF_MIN,AADM)
         REF_MAX=MAX(REF_MAX,AADM)

         KONTA=KONTA+1
         LIP(KONTA)=J
         LIR(KONTA)=IR

*        esta masa es solo una aproximacion para evitar divisiones por
*        cero. Se mantiene pero la correcta esta en CENTROMASAS_PART

         MASADM=MASADM+MASAP(J)
        END IF
       END DO


*      NIVELES AMR
       DO IR=1, NL

       KK1=0
       KK2=0
       KK1=SUM(NPART(0:IR-1))
       KK2=SUM(NPART(0:IR))

       DO J=KK1+1, KK2
        AADM=SQRT((RXPA(J)-NCLUSRX(IFI,I))**2+
     &            (RYPA(J)-NCLUSRY(IFI,I))**2+
     &            (RZPA(J)-NCLUSRZ(IFI,I))**2)


        IF(AADM<BAS*NR(IFI,I)) THEN
         REF_MIN=MIN(REF_MIN,AADM)
         REF_MAX=MAX(REF_MAX,AADM)

         KONTA=KONTA+1
         LIP(KONTA)=J
         LIR(KONTA)=IR

         MASADM=MASADM+MASAP(J)
        END IF
       END DO
       END DO

       CONTADM(1:KONTA)=0     !en principio todas estas ligadas

*      CHECKING BAD MASSES
cx       IF (MASADM.GT.0.0) THEN
       IF (MASADM.GT.0.0.AND.KONTA.GT.0) THEN

        CALL CENTROMASAS_PART(KONTA,CONTADM,LIP,
     &           U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &           CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MASA2)

       END IF

       VCM=SQRT(VCMX**2+VCMY**2+VCMZ**2)
       VCM2(IFI,I)=VCM
       MASAKK=MASADM*9.1717E18

*      ASIGNACION CM
       NCLUSRX(IFI,I)=CMX
       NCLUSRY(IFI,I)=CMY
       NCLUSRZ(IFI,I)=CMZ
       VX(IFI,I)=VCMX
       VY(IFI,I)=VCMY
       VZ(IFI,I)=VCMZ

**     FIN VCM

********************************************************************
*      UNBINDING:SCAPE VELOCITY
********************************************************************

       CONCEN=0.0
       RS=0.0
       FC=0.0
       VMAX2=0.0

       IF (MASADM.GT.0.0) THEN
        CONCEN=124.0*((MASADM*9.1717E18*(ACHE/3.66D-3)**(-1))**(-0.084))
        RS=NR(IFI,I)/CONCEN
       END IF

       CONCENTRA(IFI,I)=CONCEN

       FC=LOG(1.0+CONCEN)-(CONCEN/(1.0+CONCEN))
       IF (FC.GT.0.0) VMAX2=(CGR*MASADM*F2)/(RS*RETE*2.0*FC)

*       WRITE(*,*) 'PART.LIGADAS=', COUNT(CONTADM(1:KONTA).EQ.0)
       CONTAERR=0
       CONTAERR=COUNT(CONTADM(1:KONTA).EQ.0)

*      reordenar !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*      devuelve las particulas ordenadas de menor a mayor distancia
*      al centro (dista es la dist. al centro). Tambien reordena LIR y LIP

       DISTA=0.0
       CALL REORDENAR(KONTA,NCLUSRX(IFI,I),NCLUSRY(IFI,I),
     &   NCLUSRZ(IFI,I),RXPA,RYPA,RZPA,CONTADM,LIP,LIR,DISTA)

       FAC=0
       DO WHILE (CONTAERR.GT.0.OR.FAC.LT.4)

        FAC=FAC+1
        KONTA2=COUNT(CONTADM(1:KONTA).EQ.0)

        CALL UNBINDING4(FAC,IFI,I,REF_MIN,REF_MAX,DISTA,
     &           U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &           NR,NMASA,NCLUSRX,NCLUSRY,NCLUSRZ,
     &           LIP,LIR,KONTA,CONTADM,VX,VY,VZ)


        CALL REORDENAR(KONTA,NCLUSRX(IFI,I),NCLUSRY(IFI,I),
     &     NCLUSRZ(IFI,I),RXPA,RYPA,RZPA,CONTADM,LIP,LIR,DISTA)

        CONTAERR=abs(COUNT(CONTADM(1:KONTA).EQ.0)-KONTA2)


       END DO


       NR(IFI,I)=REF_MAX

*      AQUI: las particulas estan ordenadar de menor a mayor
*      distancia al centro!!!!!


**************************************************************
*      DENSITY PROFILE
**************************************************************

       KONTA2=COUNT(CONTADM(1:KONTA).EQ.0)

       REF_MIN=DISTA(1)
       REF_MAX=DISTA(KONTA2)

       NORMA=MAXVAL(MASAP)

       JJ=0
       DENSITOT=0.0
       RADIAL=0.0

*      perfil acumulado de masa

******************
       SALIDA=0
******************

*      SALIDA POR DENSIDAD MEDIA MENOR QUE THRESHOLD
       BAS=0.0
       BAS1=-1.0
       BAS2=0.0
*       VCMAX=-1.0

       FLAG_200=0

       DO J=1,KONTA2      !!!!! DEJO 80 por ciento de BINS DE SEGURIDAD

         VOL=PI*(4.0/3.0)*(DISTA(J)*RETE)**3

         BAS=BAS+(MASAP(LIP(J))/NORMA)

         DELTA2=NORMA*BAS/VOL/ROTE

         BAS2=BAS/DISTA(J)
         IF (BAS2.GT.BAS1) THEN
          VCMAX(IFI,I)=BAS2
          MCMAX(IFI,I)=BAS
          RCMAX(IFI,I)=DISTA(J)
          BAS1=BAS2
         END IF

CV2         IF (DELTA2.LE.200.0/OMEGA0.AND.FLAG_200.EQ.0) THEN
c_v8         IF (DELTA2.LE.200.0/OMEGA0) THEN
          IF (DELTA2.LE.200.0/OMEGAZ) THEN
          M200(IFI,I)=DELTA2*VOL*ROTE
          R200(IFI,I)=DISTA(J)
          FLAG_200=1
         END IF


       IF (DELTA2.LE.CONTRASTEC.AND.J.GT.INT(0.1*KONTA2)) THEN
          SALIDA=1
          EXIT
         END IF
       END DO            !!!!!!!!!!!!!!!

       IF (SALIDA.EQ.1) THEN
        NCAPAS(I)=J
        IF (J.EQ.KONTA2) THEN
         RSHELL=DISTA(J)
        ELSE
         RSHELL=DISTA(J+1)
        END IF

       END IF

       VCMAX(IFI,I)=VCMAX(IFI,I)*NORMA*CGR/RETE
       VCMAX(IFI,I)=SQRT(VCMAX(IFI,I))*UV
       MCMAX(IFI,I)=MCMAX(IFI,I)*NORMA*UM
       RCMAX(IFI,I)=RCMAX(IFI,I)   !*RETE
       M200(IFI,I)=M200(IFI,I)*UM
       R200(IFI,I)=R200(IFI,I)     !*RETE


       IF (SALIDA.NE.1.AND.KONTA2.NE.0) THEN   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
***    !!! DENSITOT TIENE QUE ESTAR DIMENSIONADO 0:KONTA2
       BAS=0.0
       FAC=INT(0.1*KONTA2)
       KONTA3=INT(REAL(KONTA2)/FAC)*FAC
       DO J=1,KONTA3
         BAS=BAS+(MASAP(LIP(J))/NORMA)
         IF (MOD(J,FAC).EQ.0) THEN
          JJ=JJ+1
          DENSITOT(JJ)=BAS
          RADIAL(JJ)=DISTA(J)
         END IF
       END DO
       DO J=KONTA3+1,KONTA2
         BAS=BAS+(MASAP(LIP(J))/NORMA)
       END DO
       JJ=JJ+1
       DENSITOT(JJ)=BAS
       RADIAL(JJ)=DISTA(KONTA2)

*      hasta aqui masa acumulada en bins de 10 particulas
       NSHELL_2=0
       NSHELL_2=JJ

COJO!!
       IF (NSHELL_2.GT.4) THEN
COJO!!
       DO J=INT(0.8*NSHELL_2),NSHELL_2       !!!!! DEJO 80 por cient de BINS DE SEGURIDA

         BAS=0.5*(RADIAL(J)+RADIAL(J-1))  !radio medio
         DENSA=(DENSITOT(J)-DENSITOT(J-1))/(RADIAL(J)-RADIAL(J-1))
         DENSA=DENSA/BAS/BAS

         BAS=0.5*(RADIAL(J+1)+RADIAL(J))  !radio medio
         DENSB=(DENSITOT(J+1)-DENSITOT(J))/(RADIAL(J+1)-RADIAL(J))
         DENSB=DENSB/BAS/BAS


         BAS=0.5*(RADIAL(J+2)+RADIAL(J+1))  !radio medio
         DENSC=(DENSITOT(J+2)-DENSITOT(J+1))/(RADIAL(J+2)-RADIAL(J+1))
         DENSC=DENSC/BAS/BAS

         IF (DENSC.GT.DENSB.AND.DENSB.GT.DENSA) THEN
           SALIDA=2
           EXIT
         END IF
       END DO            !!!!!!!!!!!!!!!

COJO!!
       END IF
COJO!!

*      SALIDA POR CAMBIO DE PENDIENTE DE LA DENSIDAD
       IF (SALIDA.EQ.2) THEN
        RSHELL=RADIAL(J+1)
        NCAPAS(I)=J
       END IF

       END IF       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

*      SIN CORTE
       IF (SALIDA.EQ.0.AND.KONTA2.NE.0) THEN
        J=J-1
        NCAPAS(I)=J
        RSHELL=REF_MAX
       END IF


*      YA CORTADO   !---------------------


***********************************************************
*      GUARDAMOS LAS PARTICULAS LIGADAS  DEL HALO I
***********************************************************

       KONTA2=0
       BASMAS=0.0
       DO J=1,KONTA
        IF (CONTADM(J).EQ.0) THEN

         AADMX(1)=RXPA(LIP(J))-NCLUSRX(IFI,I)
         AADMX(2)=RYPA(LIP(J))-NCLUSRY(IFI,I)
         AADMX(3)=RZPA(LIP(J))-NCLUSRZ(IFI,I)
         AADM=SQRT(AADMX(1)**2+AADMX(2)**2+AADMX(3)**2)

CV2
        VVV2=0.0
        VVV2=(U2DM(LIP(J))-VX(IFI,I))**2
     &      +(U3DM(LIP(J))-VY(IFI,I))**2
     &      +(U4DM(LIP(J))-VZ(IFI,I))**2
CV2


         IF (AADM.LE.RSHELL) THEN
          KONTA2=KONTA2+1
          BASMAS=BASMAS+(MASAP(LIP(J))/NORMA)

*********anyadido susana
          DMPCLUS(IFI,I)=DMPCLUS(IFI,I)+1
c          IF (PLOT.EQ.2) THEN
          IF (PLOT.GT.1) THEN
            IR_PARAL(DMPCLUS(IFI,I),I-IP+1)=ORIPA1(LIP(J))   !I-IP
            IP_PARAL(DMPCLUS(IFI,I),I-IP+1)=ORIPA2(LIP(J))
            MASAP_PARAL(DMPCLUS(IFI,I),I-IP+1)=MASAP(LIP(J))
          END IF

c          IF (PLOT.GT.1) THEN  ! We do ANY merger tree
c            DONDE(IFI,ORIPA1(LIP(J)),ORIPA2(LIP(J)))=I
c          ENDIF

          ANGULARM(IFI,I)=ANGULARM(IFI,I)+MASAP(LIP(J))*AADM*SQRT(VVV2)

**        INERTIA TENSOR
          DO JY=1, 3
          DO IX=1, 3
           INERTIA(IX,JY)=INERTIA(IX,JY)+AADMX(IX)*AADMX(JY)
          END DO
          END DO

**        VELOCITY OF THE FASTEST PARTICLE IN THE HALO
          VKK=0.0
          VKK=SQRT(U2DM(LIP(J))**2+U3DM(LIP(J))**2+
     &             U4DM(LIP(J))**2)

          IF (VKK.GT.VMAXCLUS(IFI,I)) THEN
           VMAXCLUS(IFI,I)=VKK
          END IF

**        CLOSEST PARTICLE TO THE CENTER OF THE HALO
          IF (AADM.LT.DIS) THEN
           DIS=AADM
           IPLIR(IFI,I)=ORIPA1(LIP(J))
           IPLIP(IFI,I)=ORIPA2(LIP(J))
          END IF

          IF(AADM.NE.0.0) THEN
           AA=0.0
           AA=((RXPA(LIP(J))-NCLUSRX(IFI,I))/AADM)*U2DM(LIP(J))+
     &        ((RYPA(LIP(J))-NCLUSRY(IFI,I))/AADM)*U3DM(LIP(J))+
     &        ((RZPA(LIP(J))-NCLUSRZ(IFI,I))/AADM)*U4DM(LIP(J))
           VR=VR+AA*MASAP(LIP(J))
          END IF


*********fin anyadido susana
         END IF    !AADM.LT.RSHELL
        END IF     !CONTADM
       END DO      !KONTA


*      CONTROL DE SEGURIDAD
       IF(DMPCLUS(IFI,I).NE.KONTA2) THEN
         WRITE(*,*) 'WARNING!', DMPCLUS(IFI,I),KONTA2
         STOP
       ENDIF


*******************************************************
*      SAVING MASSES, RADII, PROFILES AND SHAPES...
*******************************************************

       NMASA(IFI,I)=BASMAS*NORMA*9.1717E18    !!MASA2
       ANGULARM(IFI,I)=ANGULARM(IFI,I)/BASMAS
       NR(IFI,I)=RSHELL

       INERTIA(1:3,1:3)=INERTIA(1:3,1:3)/DMPCLUS(IFI,I)
       BASEIGENVAL(1:3)=0.0

       IF (DMPCLUS(IFI,I).GE.NUMPARTBAS) THEN
        CALL JACOBI(INERTIA,DIMEN,BASEIGENVAL,NROT)
        CALL SORT(BASEIGENVAL,DIMEN,DIMEN)
       END IF

       DO II=1, DIMEN
        EIGENVAL(IFI,II,I)=BASEIGENVAL(II)
        EIGENVAL(IFI,II,I)=SQRT(EIGENVAL(IFI,II,I))
       END DO

******************************************
************ FIN HALO I ******************
******************************************

*       IF (NMASA(IFI,I).EQ.0.0) KKKONTA=KKKONTA+1


       END IF   !IF REALCLUS.ne.0 (KK_ENTERO.NE.0)

*****************
       END DO   !NNCLUS I
****************


cx---------xc
       KONTA2=0
       KONTA2=SUM(DMPCLUS(IFI,1:NNCLUS(IFI)))
       PABAS=KONTA2
cx       WRITE(*,*) 'MAX. DMPITER(IFI)=', KONTA2
       IF(IP.EQ.1) WRITE(*,*) 'INITIAL MAX. DMPITER(IFI)=', KONTA2
C       IF(IP.GT.NH) WRITE(*,*) 'FINAL MAX. DMPITER(IFI)=', KONTA2
cx--------xc

*-----SOLO NECESARIO PARA EL MERGER TREE O PARA FLAG_WDM=1: ----*
       IF (PLOT.GT.1) THEN
*-----------------------------*

*****************************
c_v6:       DO I=1, NNCLUS(IFI)
       DO I=IP, IP2
****************************
       DO J=1,DMPCLUS(IFI,I)

          DMPITER(IFI)=DMPITER(IFI)+1

c          IF (FLAG_WDM.EQ.1.OR.PLOT.EQ.2) THEN
            DMLIR(IFI,DMPITER(IFI))=IR_PARAL(J,I-IP+1)
            DMLIP(IFI,DMPITER(IFI))=IP_PARAL(J,I-IP+1)
            NEW_MASAP(IFI,DMPITER(IFI))=MASAP_PARAL(J,I-IP+1)
            IX1=DMLIR(IFI,DMPITER(IFI))
            IX2=DMLIP(IFI,DMPITER(IFI))
            NHOST(IFI,IX1,IX2)=NHOST(IFI,IX1,IX2)+1
            IX3=NHOST(IFI,IX1,IX2)
            IF(IX3.GT.NMAXHOST) THEN
               WRITE(*,*)'WARNING!! NHOST>NMAXDAD'
               STOP
            END IF
            HHOST(IFI,IX3,IX1,IX2)=I
c          END IF
          !!COMMENT: tal y como esta, este if se necesita si PLOT>1!!
          !!OJO!!  PARTIRED_PLOT tiene que se >= PABAS
          IF (DMPITER(IFI).GT.PABAS) THEN
            WRITE(*,*)'WARNING! DMPITER(IFI)>',PABAS, PARTIRED_PLOT
            STOP
          END IF

       END DO  !J

****************************
       END DO      !I
****************************

*-----------------------------*
       ENDIF !PLOT
*-----------------------------*


c_v6: .............
        IP=IP2 !!! IP
        IF(IP.GE.NH) WRITE(*,*) 'FINAL MAX. DMPITER(IFI)=', KONTA2

       !+++++
       END DO
       !+++++

         IF (PLOT.GT.1) THEN
           DEALLOCATE(IP_PARAL)
           DEALLOCATE(IR_PARAL)
           DEALLOCATE(MASAP_PARAL)
        END IF
c_v6: .............


*************************************************
******** GENERAL CHECKING ***********************
*************************************************

       WRITE(*,*)'HALOES WITHOUT MASS=',
     &        COUNT(NMASA(IFI,1:NNCLUS(IFI)).LE.0.0)

       WRITE(*,*)'HALOES WITH MASS=0',
     &        COUNT(NMASA(IFI,1:NNCLUS(IFI)).EQ.0.0)
*       WRITE(*,*)'KKKONTA=',KKKONTA

       WRITE(*,*) 'Total number of particles within halos:',
     &            SUM(DMPCLUS(IFI,1:NNCLUS(IFI)))


       WRITE(*,*)'After refining with DM particles...'
       WRITE(*,*)'===================================='
       DO I=0,NL
       WRITE(*,*)'Haloes at level ', I,' =',
     &            COUNT(NLEVHAL(IFI,1:NNCLUS(IFI)).EQ.I.
     &            AND.REALCLUS(IFI,1:NNCLUS(IFI)).NE.0),
     &            COUNT(REALCLUS(IFI,1:NNCLUS(IFI)).EQ.-1)
       END DO
       WRITE(*,*)'===================================='

*************************************************
*************************************************

       WRITE(*,*) 'ESTIMATION HALOES_1 (AFTER UNBINDING)'
       WRITE(*,*)'=================================='
       DO I=1, NNCLUS(IFI)
       IF (REALCLUS(IFI,I).NE.0) THEN
       WRITE(*,*) I, NMASA(IFI,I),NR(IFI,I),NCLUSRX(IFI,I),
     &            NCLUSRY(IFI,I),NCLUSRZ(IFI,I)
       END IF
       END DO
       WRITE(*,*)'=================================='



*      REPETIMOS ESTO OTRA VEZ
*STEP 1) (2)************************************
************REMOVING POOR HALOES***************
************************************************
       KONTA2=0
       NUMPARTBAS=NUMPART
!$OMP  PARALLEL DO SHARED(NNCLUS,IFI,DMPCLUS,NUMPARTBAS,
!$OMP+             REALCLUS),PRIVATE(I,KK_ENTERO),
!$OMP+             REDUCTION(+:KONTA2)
       DO I=1, NNCLUS(IFI)

       KK_ENTERO=0
       KK_ENTERO=DMPCLUS(IFI,I)

       IF (KK_ENTERO.LT.NUMPARTBAS) THEN

       REALCLUS(IFI,I)=0
       KONTA2=KONTA2+1

       END IF

       END DO

       WRITE(*,*)'RE-CHECKING POOR HALOS----->', KONTA2


*PASO 2)****************************************
************RUBBISH (desbordamientos)***********
************************************************

       KONTA2=0

       DO I=1, NNCLUS(IFI)

       IF (REALCLUS(IFI,I).NE.0) THEN

       DO J=1, NNCLUS(IFI)


       IF (REALCLUS(IFI,J).NE.0.AND.NLEVHAL(IFI,J).GT.
     &      NLEVHAL(IFI,I)) THEN

       DIS=0.0
       DIS=SQRT((NCLUSRX(IFI,I)-NCLUSRX(IFI,J))**2+
     &          (NCLUSRY(IFI,I)-NCLUSRY(IFI,J))**2+
     &          (NCLUSRZ(IFI,I)-NCLUSRZ(IFI,J))**2)

       A1=0.0
       A1=MIN(NMASA(IFI,I),NMASA(IFI,J))/MAX(NMASA(IFI,I),NMASA(IFI,J))

       A2=0.0
       IF (MIN(ABS(VCM2(IFI,I)),ABS(VCM2(IFI,J))).NE.0.0) THEN
       A2=(ABS(VCM2(IFI,I)-VCM2(IFI,J)))/
     &                            MAX(ABS(VCM2(IFI,I)),ABS(VCM2(IFI,J)))
       END IF

       A3=0.0
       A3=MIN(NR(IFI,I),NR(IFI,J))

*       IF (DIS.LT.1.01*A3.AND.A1.GT.0.2.AND.A2.LT.3.0.OR.A1.GT.0.6) THEN
       IF (DIS.LT.1.01*A3.AND.A1.GT.0.2.AND.A2.LT.3.0) THEN
*       IF (DIS.LT.1.01*A3.AND.A1.GT.0.2.AND.A2.LT.5.0) THEN
*       IF (DIS.LT.1.01*A3.AND.A1.GT.0.2) THEN

       IF (NMASA(IFI,I).GT.NMASA(IFI,J)) THEN

       REALCLUS(IFI,J)=0
       KONTA2=KONTA2+1

       ELSE

       REALCLUS(IFI,I)=0
       KONTA2=KONTA2+1

       END IF

       END IF

       END IF   !realclus

       END DO
       END IF   !realclus
       END DO

       WRITE(*,*)'CHECKING RUBBISH----->', KONTA2


*PASO 3)********************************************
*******************SUBSTRUCTURE*********************
****************************************************

       SUBHALOS=0
       KONTA2=0
       DO I=1, NNCLUS(IFI)

       IF (REALCLUS(IFI,I).EQ.-1) THEN

       CONCEN=0.0
       RS=0.0
       FC=0.0
       VMAX2=0.0

       IF (NMASA(IFI,I).GT.0.0) THEN
       CONCEN=124.0*((NMASA(IFI,I)*(ACHE/3.66D-3)**(-1))**(-0.084))
       RS=NR(IFI,I)/CONCEN
       END IF

       CONCENTRA(IFI,I)=CONCEN

       FC=LOG(1.0+CONCEN)-(CONCEN/(1.0+CONCEN))

       IF (FC.GT.0.0) THEN
       VMAX2=(CGR*(NMASA(IFI,I)/9.1717E18)*F2)/(RS*RETE*2.0*FC)
       END IF


       DO J=1, NNCLUS(IFI)

       DIS=0.0
       DIS=SQRT((NCLUSRX(IFI,I)-NCLUSRX(IFI,J))**2+
     &          (NCLUSRY(IFI,I)-NCLUSRY(IFI,J))**2+
     &          (NCLUSRZ(IFI,I)-NCLUSRZ(IFI,J))**2)

       VVV2=0.0
       VVV2=(VCM2(IFI,I)-VCM2(IFI,J))**2

       VESC2=0.0
       IF (RS.NE.0.AND.DIS*F2/RS.NE.0) THEN
       VESC2=(4.0*VMAX2*LOG(1.0+DIS/RS))/(DIS*F2/RS)
       END IF

       A1=0.0
       A1=MIN(NMASA(IFI,I),NMASA(IFI,J))/MAX(NMASA(IFI,I),NMASA(IFI,J))

       IF (REALCLUS(IFI,J)==-1) THEN

       IF(NLEVHAL(IFI,J).GT.NLEVHAL(IFI,I).AND.
     &    DIS.LE.1.0*NR(IFI,I).AND.A1.LE.0.2.AND.VVV2.LE.VESC2) THEN


       REALCLUS(IFI,J)=I
       SUBHALOS(I)=SUBHALOS(I)+1
       IF (SUBHALOS(I).GT.NMAXSUB) THEN
       WRITE(*,*)'WARNING: DEMASIASOS SUBHALOS!!', IFI, I, SUBHALOS
       STOP
       END IF
       KONTA2=KONTA2+1

       END IF
       END IF

       END DO
       END IF
       END DO

       WRITE(*,*)'CHECKING SUBSTRUCTURE----->', KONTA2

       WRITE(*,*)'TOTAL NUMBER OF HALOS=',
     &            COUNT(REALCLUS(IFI,1:NNCLUS(IFI)).EQ.-1)


*********************************************************
*      CHECKING....2
*********************************************************
       WRITE(*,*) '---------------------------------'
       WRITE(*,*) 'RE-CHECKING HIERARCHY:'
       WRITE(*,*) '---------------------------------'
       KONTA2=0
       KONTA2=COUNT(REALCLUS(IFI,1:NNCLUS(IFI)).EQ.-1)
       WRITE(*,*) 'REAL HALOS----->', KONTA2
       KONTA2=0
       KONTA2=COUNT(REALCLUS(IFI,1:NNCLUS(IFI)).EQ.0)
       WRITE(*,*) 'REMOVED HALOS----->', KONTA2
       KONTA2=0
       KONTA2=COUNT(REALCLUS(IFI,1:NNCLUS(IFI)).GT.0)
       WRITE(*,*) 'SUBSTRUCTURE----->', KONTA2
       WRITE(*,*) '---------------------------------'

       WRITE(*,*)'At the end...'
       WRITE(*,*)'=================================='
       DO I=0,NL
       WRITE(*,*)'Haloes at level ', I,' =',
     &            COUNT(NLEVHAL(IFI,1:NNCLUS(IFI)).EQ.I.
     &            AND.REALCLUS(IFI,1:NNCLUS(IFI)).NE.0),
     &            COUNT(REALCLUS(IFI,1:NNCLUS(IFI)).EQ.-1)
       END DO
       WRITE(*,*)'=================================='


****************************************************
****************************************************
****************************************************


*************************************************
*************************************************
*===================Families===============
       KONTA2=0
c       KONTA2=COUNT(REALCLUS(IFI,1:NNCLUS(IFI)).EQ.-1)
       KONTA2=COUNT(REALCLUS(IFI,1:NNCLUS(IFI)).NE.0)
       CALL NOMFILE3(ITER,FILE3)
       FILERR3='./output_files/'//FILE3
       OPEN(3,FILE=FILERR3,STATUS='UNKNOWN')

       IF (FLAG_WDM.EQ.1) THEN
        CALL NOMFILE7(ITER,FILE7)
        FILERR7='./output_files/'//FILE7
        OPEN(4,FILE=FILERR7,STATUS='UNKNOWN',FORM='UNFORMATTED')
       END IF


       IF (FLAG_WDM.EQ.1) WRITE(4) KONTA2

       WRITE(3,*) '*********************NEW ITER*******************'
       WRITE(3,*) ITER, NNCLUS(IFI), KONTA2, ZETA
       WRITE(3,*) '************************************************'
       KONTA2=0

       DO I=1, NNCLUS(IFI)

       IF (REALCLUS(IFI,I).NE.0) THEN

         WRITE(3,135) I,NCLUSRX(IFI,I),NCLUSRY(IFI,I),NCLUSRZ(IFI,I),
     &         NMASA(IFI,I),NR(IFI,I),DMPCLUS(IFI,I),
     &         REALCLUS(IFI,I), NLEVHAL(IFI,I),SUBHALOS(I),
     &         EIGENVAL(IFI,1,I),EIGENVAL(IFI,2,I),EIGENVAL(IFI,3,I),
     &         VCM2(IFI,I)*UV,CONCENTRA(IFI,I),ANGULARM(IFI,I)*1.0E14,
     &         VCMAX(IFI,I),MCMAX(IFI,I),RCMAX(IFI,I),
     &         M200(IFI,I),R200(IFI,I),
     &         VX(IFI,I)*UV,VY(IFI,I)*UV,VZ(IFI,I)*UV


       IF (FLAG_WDM.EQ.1) THEN

       WRITE(4) I,REALCLUS(IFI,I),DMPCLUS(IFI,I)

       KK2=0
       KKK2=0
       KK2=SUM(DMPCLUS(IFI,1:I-1))+1
       KKK2=SUM(DMPCLUS(IFI,1:I))

       DO J2=KK2, KKK2

       IX1=0
       IX2=0
       IX1=DMLIP(IFI2,J2)

*      IR=0
       DO J=1, NPART(0)
        IX2=ORIPA2(J)
        IF (IX2.EQ.IX1) THEN
        WRITE(4) RXPA(J),RYPA(J),RZPA(J),MASAP(J),
     &           U2DM(J),U3DM(J),U4DM(J)
        END IF
       END DO

*      NIVELES AMR
       DO IR=1, NL

       KK1=0
       KK3=0
       KK1=SUM(NPART(0:IR-1))
       KK3=SUM(NPART(0:IR))

       DO J=KK1+1, KK3

       IX2=ORIPA2(J)
       IF (IX2.EQ.IX1) THEN
        WRITE(4) RXPA(J),RYPA(J),RZPA(J),MASAP(J),
     &           U2DM(J),U3DM(J),U4DM(J)
       END IF

       END DO
       END DO  !ir

       END DO  !J2
       END IF  !FLAG_WDM


       END IF  !realclus

       END DO


135    FORMAT (i8,'   ',3F15.6,'   ',e15.6,'   ',F15.6,'   ',i8,'    ',
     &         3i5,'   ',7F15.6,'   ',e15.6,'   ',F15.6,'   ',
     &         e15.6,'   ',F15.6,'   ',3F15.6)

       CLOSE(3)
       IF (FLAG_WDM.EQ.1) CLOSE(4)

*==========================================
*************************************************
*************************************************

       WRITE(*,*)'END ITER', ITER
       CALL IDATE(DATE)
       CALL ITIME(TIME)
       WRITE(*,*) 'DATE=',DATE(1),'/',DATE(2),'/',DATE(3)
       WRITE(*,*) 'TIME=',TIME(1),':',TIME(2),':',TIME(3)


*/////////////////////////////////////////////////
*      COMPLETE MERGER TREE (WITH % OF MASSES)
*////////////////////////////////////////////////

       IF (PLOT.EQ.2) THEN

       IF (CONTAITER.EQ.2) THEN

       IF (MARK(IFI2).EQ.1.AND.MARK(IFI2-1).EQ.1)THEN

       WRITE(*,*) 'MERGER TREE CHECKING:', CONTAITER
       DO I=1,NFILE
        WRITE(*,*) ZETAS(I),TIEMPO(I),NNCLUS(I),DMPITER(I)
       END DO

CCCCCCCCCCCCCCCCCCCC
       WRITE(*,*)'STARTING MERGER TREE', ITER-EVERY,'-',ITER
       CALL IDATE(DATE)
       CALL ITIME(TIME)
       WRITE(*,*) 'DATE=',DATE(1),'/',DATE(2),'/',DATE(3)
       WRITE(*,*) 'TIME=',TIME(1),':',TIME(2),':',TIME(3)
CCCCCCCCCCCCCCCCCCCCC

*/////////////CALLING MERGER TREE/////////////////////////*
         CALL MERGER(NFILE,NMASA,NNCLUS,DMPITER,
     &              DMPCLUS,REALCLUS)
*////////////////WRITING MERGER_TREE/////////////////////*

       CALL NOMFILE4(ITER-EVERY,FILE3)
       FILERR3='./output_files/'//FILE3
       OPEN(3,FILE=FILERR3,STATUS='UNKNOWN')

       DO I=NFILE, 2, -1
       WRITE(3,*)'================NEW ITER==================='
       KONTA2=0

c       KONTA2=COUNT(REALCLUS(I,1:NNCLUS(I)).EQ.-1.AND.
c     &        NLEVHAL(I,1:NNCLUS(I)).EQ.0)

C      NEW solo de halos reales
       KONTA2=COUNT(REALCLUS(I,1:NNCLUS(I)).EQ.-1)
       WRITE(3,*) I, ZETAS(I),TIEMPO(I),KONTA2
       KONTA2=0

c       KONTA2=COUNT(REALCLUS(I-1,1:NNCLUS(I-1)).EQ.-1.AND.
c     &        NLEVHAL(I-1,1:NNCLUS(I-1)).EQ.0)

C      NEW solo de halos reales
       KONTA2=COUNT(REALCLUS(I-1,1:NNCLUS(I-1)).EQ.-1)
       WRITE(3,*) I-1, ZETAS(I-1),TIEMPO(I-1),KONTA2
       KONTA2=0

       DO J=1, NNCLUS(I)
c       IF (REALCLUS(I,J).EQ.-1.AND.NLEVHAL(I,J).EQ.0) THEN

C      NEW solo de halos reales
       IF (REALCLUS(I,J).EQ.-1) THEN
       KONTA2=KONTA2+1
       WRITE(3,*) '---------------new clus----------------------'
       WRITE(3,*) KONTA2,J,NMASA(I,J),NDAD(I,J),NR(I,J)
       WRITE(3,*) NCLUSRX(I,J),NCLUSRY(I,J),NCLUSRZ(I,J),NLEVHAL(I,J)
       DO K=1, NDAD(I,J)
       WRITE(3,*) DAD(I,K,J),RATIO(I,K,J),NMASA(I-1,DAD(I,K,J)),
     &            NR(I-1,DAD(I,K,J))
       WRITE(3,*) NCLUSRX(I-1,DAD(I,K,J)),NCLUSRY(I-1,DAD(I,K,J)),
     &            NCLUSRZ(I-1,DAD(I,K,J)),NLEVHAL(I-1,DAD(I,K,J))
       END DO
       END IF
       END DO

       END DO

       CLOSE(3)

       CONTAITER=1

*      MACHACAMOS LAS VARIABLES DE ITER 1
C       WRITE(*,*)'1:',NNCLUS(1),ZETAS(1),TIEMPO(1),DMPITER(1)

cxcx       NMAXNCLUSBAS=NNCLUS(1)
       NMAXNCLUSBAS=NMAXNCLUS    !dimension original maxima
!$OMP  PARALLEL DO SHARED(NMAXNCLUSBAS,NMASA,NR,CONCENTRA,
!$OMP+                   ANGULARM,NCLUSRX,NCLUSRY,NCLUSRZ,
!$OMP+                   VCM2,VMAXCLUS,VX,VY,VZ,
!$OMP+                   VCMAX,MCMAX,RCMAX,DMPCLUS,
!$OMP+                   M200,R200,IPLIP,IPLIR,REALCLUS,
!$OMP+                   NLEVHAL,EIGENVAL),PRIVATE(I)
       DO I=1,NMAXNCLUSBAS
          NMASA(1,I)=0.0
          NR(1,I)=0.0
          CONCENTRA(1,I)=0.0
          ANGULARM(1,I)=0.0
          NCLUSRX(1,I)=0.0
          NCLUSRY(1,I)=0.0
          NCLUSRZ(1,I)=0.0
          VCM2(1,I)=0.0
          VMAXCLUS(1,I)=0.0
          VX(1,I)=0.0
          VY(1,I)=0.0
          VZ(1,I)=0.0
          VCMAX(1,I)=0.0
          MCMAX(1,I)=0.0
          RCMAX(1,I)=0.0
          M200(1,I)=0.0
          R200(1,I)=0.0
          IPLIP(1,I)=0
          IPLIR(1,I)=0
          DMPCLUS(1,I)=0
          REALCLUS(1,I)=0
          NLEVHAL(1,I)=0
          EIGENVAL(1,:,I)=0.0
        END DO

        IF (PLOT.GT.1) THEN
cx      NMAXNCLUSBAS=DMPITER(1)
        NMAXNCLUSBAS=PARTIRED_PLOT   !dimension original maxima
!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,DMLIP,DMLIR,
!$OMP+                   NEW_MASAP),
!$OMP+            PRIVATE(I)
       DO I=1,NMAXNCLUSBAS
         DMLIP(1,I)=0
         DMLIR(1,I)=0
         NEW_MASAP(1,I)=0.0
       END DO

       NBASPART_PLOT=PARTIRED
!$OMP PARALLEL DO SHARED(NBASPART_PLOT,NHOST,HHOST),
!$OMP+            PRIVATE(I)
      DO I=1,NBASPART_PLOT
         NHOST(1,:,I)=0
         HHOST(1,:,:,I)=0
      END DO


       NBASPART_PLOT=NMAXCLUS_PLOT
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,DAD,NDAD),
!$OMP+            PRIVATE(I)
        DO I=1,NBASPART_PLOT
           NDAD(1,I)=0
           DAD(1,:,I)=0
       END DO

      IF (PLOT.EQ.2) THEN
        NBASPART_PLOT=NMAXCLUS_PLOT
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,RATIO),
!$OMP+            PRIVATE(I)
       DO I=1,NBASPART_PLOT
          RATIO(1,:,I)=0.0
       END DO
      END IF


      END IF  !PLOT>1


       NNCLUS(1)=0
       ZETAS(1)=0.0
       TIEMPO(1)=0.0
       DMPITER(1)=0

*      GUARDAMOS LAS VARIABLES DE ITER 2 EN 1
C       WRITE(*,*)'1...:',NNCLUS(1),ZETAS(1),TIEMPO(1),DMPITER(1)

       NNCLUS(1)=NNCLUS(2)
       ZETAS(1)=ZETAS(2)
       TIEMPO(1)=TIEMPO(2)


       NMAXNCLUSBAS=NNCLUS(2)
!$OMP  PARALLEL DO SHARED(NMAXNCLUSBAS,NMASA,NR,CONCENTRA,
!$OMP+                   ANGULARM,NCLUSRX,NCLUSRY,NCLUSRZ,
!$OMP+                   VCM2,VMAXCLUS,VX,VY,VZ,
!$OMP+                   VCMAX,MCMAX,RCMAX,DMPCLUS,
!$OMP+                   M200,R200,IPLIP,IPLIR,REALCLUS,
!$OMP+                   NLEVHAL,EIGENVAL),PRIVATE(I)
       DO I=1,NMAXNCLUSBAS
          NMASA(1,I)=NMASA(2,I)
          NR(1,I)=NR(2,I)
          CONCENTRA(1,I)=CONCENTRA(2,I)
          ANGULARM(1,I)=ANGULARM(2,I)
          NCLUSRX(1,I)=NCLUSRX(2,I)
          NCLUSRY(1,I)=NCLUSRY(2,I)
          NCLUSRZ(1,I)=NCLUSRZ(2,I)
          VCM2(1,I)=VCM2(2,I)
          VMAXCLUS(1,I)=VMAXCLUS(2,I)
          VX(1,I)=VX(2,I)
          VY(1,I)=VY(2,I)
          VZ(1,I)=VZ(2,I)
          VCMAX(1,I)=VCMAX(2,I)
          MCMAX(1,I)=MCMAX(2,I)
          RCMAX(1,I)=RCMAX(2,I)
          M200(1,I)=M200(2,I)
          R200(1,I)=R200(2,I)
          IPLIP(1,I)=IPLIP(2,I)
          IPLIR(1,I)=IPLIR(2,I)
          DMPCLUS(1,I)=DMPCLUS(2,I)
          REALCLUS(1,I)=REALCLUS(2,I)
          NLEVHAL(1,I)=NLEVHAL(2,I)
          EIGENVAL(1,1:3,I)=EIGENVAL(2,1:3,I)
       END DO


       DMPITER(1)=DMPITER(2)

       IF (PLOT.GT.1) THEN

       NMAXNCLUSBAS=DMPITER(2)
!$OMP  PARALLEL DO SHARED(NMAXNCLUSBAS,DMLIP,DMLIR,
!$OMP+               NEW_MASAP),PRIVATE(I)
       DO I=1, NMAXNCLUSBAS
         DMLIP(1,I)=DMLIP(2,I)
         DMLIR(1,I)=DMLIR(2,I)
         NEW_MASAP(1,I)=NEW_MASAP(2,I)
        END DO

      NBASPART_PLOT=PARTIRED
!$OMP PARALLEL DO SHARED(NBASPART_PLOT,NHOST,HHOST),
!$OMP+            PRIVATE(I)
       DO I=1,NBASPART_PLOT
         NHOST(1,:,I)=NHOST(2,:,I)
         HHOST(1,:,:,I)=HHOST(2,:,:,I)
       END DO

       NBASPART_PLOT=NNCLUS(2)
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,DAD,NDAD),
!$OMP+            PRIVATE(I)
       DO I=1,NBASPART_PLOT
         NDAD(1,I)=NDAD(2,I)
c         DAD(1,1:NDAD(1,I),I)=DAD(2,1:NDAD(2,I),I)
         DAD(1,:,I)=DAD(2,:,I)
       END DO

       IF (PLOT.EQ.2) THEN

        NBASPART_PLOT=NNCLUS(2)
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,RATIO),
!$OMP+            PRIVATE(I)
         DO I=1,NBASPART_PLOT
c           RATIO(1,1:NDAD(1,I),I)=RATIO(2,1:NDAD(2,I),I)
           RATIO(1,:,I)=RATIO(2,:,I)
         END DO
       END IF


       END IF


C       WRITE(*,*)'AHORA 1 ES 2:',NNCLUS(1),ZETAS(1),TIEMPO(1),
C     &            DMPITER(1)

*      MACHACAMOS LAS VARIABLES DE ITER 2

cx       NMAXNCLUSBAS=NNCLUS(2)
        NMAXNCLUSBAS=NMAXNCLUS    !dimension original maxima
!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,NMASA,NR,CONCENTRA,
!$OMP+                   ANGULARM,NCLUSRX,NCLUSRY,NCLUSRZ,
!$OMP+                   VCM2,VMAXCLUS,VX,VY,VZ,
!$OMP+                   VCMAX,MCMAX,RCMAX,DMPCLUS,
!$OMP+                   M200,R200,IPLIP,IPLIR,REALCLUS,
!$OMP+                   NLEVHAL,EIGENVAL),PRIVATE(I)
       DO I=1, NMAXNCLUSBAS
         NMASA(2,I)=0.0
         NR(2,I)=0.0
         NCLUSRX(2,I)=0.0
         NCLUSRY(2,I)=0.0
         NCLUSRZ(2,I)=0.0
         VCM2(2,I)=0.0
         VMAXCLUS(2,I)=0.0
         VX(2,I)=0.0
         VY(2,I)=0.0
         VZ(2,I)=0.0
         VCMAX(2,I)=0.0
         MCMAX(2,I)=0.0
         RCMAX(2,I)=0.0
         CONCENTRA(2,I)=0.0
         ANGULARM(2,I)=0.0
         IPLIP(2,I)=0
         IPLIR(2,I)=0
         M200(2,I)=0.0
         R200(2,I)=0.0
         REALCLUS(2,I)=0
         NLEVHAL(2,I)=0
         DMPCLUS(2,I)=0
         EIGENVAL(2,:,I)=0.0
       END DO


       IF (PLOT.GT.1) THEN

cx       NMAXNCLUSBAS=DMPITER(2)
        NMAXNCLUSBAS=PARTIRED_PLOT   !dimension original maxima
!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,DMLIP,DMLIR,
!$OMP+              NEW_MASAP),PRIVATE(I)
       DO I=1, NMAXNCLUSBAS
         DMLIP(2,I)=0
         DMLIR(2,I)=0
         NEW_MASAP(2,I)=0.0
       END DO

       NBASPART_PLOT=PARTIRED
!$OMP PARALLEL DO SHARED(NBASPART_PLOT,NHOST,HHOST),
!$OMP+            PRIVATE(I)
      DO I=1,NBASPART_PLOT
         NHOST(2,:,I)=0
         HHOST(2,:,:,I)=0
      END DO

       NBASPART_PLOT=NMAXNCLUS
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,DAD,NDAD),
!$OMP+            PRIVATE(I)
       DO I=1,NBASPART_PLOT
          NDAD(2,I)=0
          DAD(2,:,I)=0
      END DO

      IF (PLOT.EQ.2) THEN

        NBASPART_PLOT=NMAXNCLUS
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,RATIO),
!$OMP+            PRIVATE(I)
         DO I=1,NBASPART_PLOT
          RATIO(2,:,I)=0.0
         END DO
      END IF


      END IF


       NNCLUS(2)=0
       ZETAS(2)=0.0
       TIEMPO(2)=0.0
       DMPITER(2)=0

C       WRITE(*,*)'DE MOMENTO 2 VACIO:',NNCLUS(2),ZETAS(2),TIEMPO(2),
C     &            DMPITER(2)

CCCCCCCCCCCCCCCCCCCCC
       WRITE(*,*)'END MERGER TREE', ITER-EVERY,'-',ITER
       CALL IDATE(DATE)
       CALL ITIME(TIME)
       WRITE(*,*) 'DATE=',DATE(1),'/',DATE(2),'/',DATE(3)
       WRITE(*,*) 'TIME=',TIME(1),':',TIME(2),':',TIME(3)
CCCCCCCCCCCCCCCCCCCCC

       END IF     !MARK
       END IF     !CONTAITER

*//////////////////////END WRITING////////////////////////

       END IF   !PLOT.EQ.2

*/////////////////////////////////////////////////
*      MERGER REDUCIDO (SOLO LINEA PRINCIPAL)
*////////////////////////////////////////////////

       IF (PLOT.EQ.3) THEN


       IF (CONTAITER.EQ.2) THEN

       IF (MARK(IFI2).EQ.1.AND.MARK(IFI2-1).EQ.1)THEN


C       WRITE(*,*)'AHORA 2 ES:',NNCLUS(2),ZETAS(2),TIEMPO(2),
C     &            DMPITER(2)

       WRITE(*,*) '---------------------------------------'
       WRITE(*,*) 'VAMOS A POR EL MERGER TREE:', CONTAITER
       DO I=1,NFILE
       WRITE(*,*) ZETAS(I),TIEMPO(I),NNCLUS(I),DMPITER(I)
       END DO

CCCCCCCCCCCCCCCCCCCCC
       WRITE(*,*)'STARTING MERGER TREE', ITER-EVERY,'-',ITER
       CALL IDATE(DATE)
       CALL ITIME(TIME)
       WRITE(*,*) 'DATE=',DATE(1),'/',DATE(2),'/',DATE(3)
       WRITE(*,*) 'TIME=',TIME(1),':',TIME(2),':',TIME(3)
CCCCCCCCCCCCCCCCCCCCC

*/////////////CALLING REDUCED MERGER TREE////////////////*

       CALL REDUCED_MERGER(NNCLUS,REALCLUS,IPLIR,IPLIP,
     &                     NCLUSRX,NCLUSRY,NCLUSRZ)

*////////////////WRITING REDUCED MERGER_TREE/////////////*

       CALL NOMFILE5(ITER-EVERY,FILE3)
       FILERR3='./output_files/'//FILE3
       OPEN(3,FILE=FILERR3,STATUS='UNKNOWN')

       DO I=NFILE, 2, -1
       WRITE(3,*)'================NEW ITER==================='
       KONTA2=0
C       KONTA2=COUNT(REALCLUS(I,1:NNCLUS(I)).EQ.-1.AND.
C     &        NLEVHAL(I,1:NNCLUS(I)).EQ.0)
C      NEW
       KONTA2=COUNT(REALCLUS(I,1:NNCLUS(I)).EQ.-1)
       WRITE(3,*) I, ZETAS(I),TIEMPO(I),KONTA2
       KONTA2=0
C       KONTA2=COUNT(REALCLUS(I-1,1:NNCLUS(I-1)).EQ.-1.AND.
C     &        NLEVHAL(I-1,1:NNCLUS(I-1)).EQ.0)
C      NEW
       KONTA2=COUNT(REALCLUS(I-1,1:NNCLUS(I-1)).EQ.-1)
       WRITE(3,*) I-1, ZETAS(I-1),TIEMPO(I-1),KONTA2
       KONTA2=0
       DO J=1, NNCLUS(I)
C       IF (REALCLUS(I,J).EQ.-1.AND.NLEVHAL(I,J).EQ.0) THEN
C      NEW
       IF (REALCLUS(I,J).EQ.-1) THEN
       KONTA2=KONTA2+1
       WRITE(3,*) '---------------new clus----------------------'
       WRITE(3,*) KONTA2,J,NMASA(I,J),NDAD(I,J),NR(I,J)
       WRITE(3,*) NCLUSRX(I,J),NCLUSRY(I,J),NCLUSRZ(I,J),NLEVHAL(I,J)
       DO K=1, NDAD(I,J)
       WRITE(3,*) DAD(I,K,J),NMASA(I-1,DAD(I,K,J)),NR(I-1,DAD(I,K,J))
       WRITE(3,*) NCLUSRX(I-1,DAD(I,K,J)),NCLUSRY(I-1,DAD(I,K,J)),
     &            NCLUSRZ(I-1,DAD(I,K,J)),NLEVHAL(I-1,DAD(I,K,J))
       END DO
       END IF
       END DO

       END DO

       CLOSE(3)

       CONTAITER=1


*      MACHACAMOS LAS VARIABLES DE ITER 1
C       WRITE(*,*)'1:',NNCLUS(1),ZETAS(1),TIEMPO(1),DMPITER(1)


cx      NMAXNCLUSBAS=NNCLUS(1)
        NMAXNCLUSBAS=NMAXNCLUS    !dimension original maxima
!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,NMASA,NR,CONCENTRA,
!$OMP+                   ANGULARM,NCLUSRX,NCLUSRY,NCLUSRZ,
!$OMP+                   VCM2,VMAXCLUS,VX,VY,VZ,
!$OMP+                   VCMAX,MCMAX,RCMAX,DMPCLUS,
!$OMP+                   M200,R200,IPLIP,IPLIR,REALCLUS,
!$OMP+                   NLEVHAL,EIGENVAL),PRIVATE(I)
       DO I=1,NMAXNCLUSBAS
          NMASA(1,I)=0.0
          NR(1,I)=0.0
          NCLUSRX(1,I)=0.0
          NCLUSRY(1,I)=0.0
          NCLUSRZ(1,I)=0.0
          VCM2(1,I)=0.0
          VMAXCLUS(1,I)=0.0
          VX(1,I)=0.0
          VY(1,I)=0.0
          VZ(1,I)=0.0
          VCMAX(1,I)=0.0
          MCMAX(1,I)=0.0
          RCMAX(1,I)=0.0
          CONCENTRA(1,I)=0.0
          ANGULARM(1,I)=0.0
          IPLIP(1,I)=0
          IPLIR(1,I)=0
          M200(1,I)=0.0
          R200(1,I)=0.0
          DMPCLUS(1,I)=0
          REALCLUS(1,I)=0
          NLEVHAL(1,I)=0
          EIGENVAL(1,:,I)=0.0
       END DO


       IF (PLOT.GT.1) THEN
cx       NMAXNCLUSBAS=DMPITER(1)
         NMAXNCLUSBAS=PARTIRED_PLOT   !dimension original maxima
!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,DMLIP,DMLIR,
!$OMP+              NEW_MASAP),PRIVATE(I)
       DO I=1, NMAXNCLUSBAS
         DMLIP(1,I)=0
         DMLIR(1,I)=0
         NEW_MASAP(1,I)=0.0
       END DO

      NBASPART_PLOT=PARTIRED
!$OMP PARALLEL DO SHARED(NBASPART_PLOT,NHOST,HHOST),
!$OMP+            PRIVATE(I)
       DO I=1,NBASPART_PLOT
         NHOST(1,:,I)=0
         HHOST(1,:,:,I)=0
       END DO


      NBASPART_PLOT=NMAXNCLUS
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,DAD,NDAD),
!$OMP+            PRIVATE(I)
       DO I=1,NBASPART_PLOT
         NDAD(1,I)=0
         DAD(1,:,I)=0
       END DO

        IF (PLOT.EQ.2) THEN
        NBASPART_PLOT=NMAXNCLUS
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,RATIO),
!$OMP+            PRIVATE(I)
         DO I=1,NBASPART_PLOT
            RATIO(1,:,I)=0.0
         END DO
        END IF

       END IF

       NNCLUS(1)=0
       ZETAS(1)=0.0
       TIEMPO(1)=0.0
       DMPITER(1)=0


*      GUARDAMOS LAS VARIABLES DE ITER 2 EN 1

C       WRITE(*,*)'1...:',NNCLUS(1),ZETAS(1),TIEMPO(1),DMPITER(1)

       NNCLUS(1)=NNCLUS(2)
       ZETAS(1)=ZETAS(2)
       TIEMPO(1)=TIEMPO(2)

      NMAXNCLUSBAS=NNCLUS(2)
!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,NMASA,NR,CONCENTRA,
!$OMP+                   ANGULARM,NCLUSRX,NCLUSRY,NCLUSRZ,
!$OMP+                   VCM2,VMAXCLUS,VX,VY,VZ,
!$OMP+                   VCMAX,MCMAX,RCMAX,DMPCLUS,
!$OMP+                   M200,R200,IPLIP,IPLIR,REALCLUS,
!$OMP+                   NLEVHAL,EIGENVAL),PRIVATE(I)
       DO I=1,NMAXNCLUSBAS
         NMASA(1,I)=NMASA(2,I)
         NR(1,I)=NR(2,I)
         NCLUSRX(1,I)=NCLUSRX(2,I)
         NCLUSRY(1,I)=NCLUSRY(2,I)
         NCLUSRZ(1,I)=NCLUSRZ(2,I)
         VCM2(1,I)=VCM2(2,I)
         VMAXCLUS(1,I)=VMAXCLUS(2,I)
         VX(1,I)=VX(2,I)
         VY(1,I)=VY(2,I)
         VZ(1,I)=VZ(2,I)
         VCMAX(1,I)=VCMAX(2,I)
         MCMAX(1,I)=MCMAX(2,I)
         RCMAX(1,I)=RCMAX(2,I)
         CONCENTRA(1,I)=CONCENTRA(2,I)
         ANGULARM(1,I)=ANGULARM(2,I)
         IPLIP(1,I)=IPLIP(2,I)
         IPLIR(1,I)=IPLIR(2,I)
         M200(1,I)=M200(2,I)
         R200(1,I)=R200(2,I)
         REALCLUS(1,I)=REALCLUS(2,I)
c         NDAD(1,I)=NDAD(2,I)
         NLEVHAL(1,I)=NLEVHAL(2,I)
         DMPCLUS(1,I)=DMPCLUS(2,I)
         EIGENVAL(1,1:3,I)=EIGENVAL(2,1:3,I)
        END DO

       DMPITER(1)=DMPITER(2)

       IF (PLOT.GT.1) THEN

       NMAXNCLUSBAS=DMPITER(2)
!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,DMLIP,DMLIR,
!$OMP+                 NEW_MASAP),PRIVATE(I)
       DO I=1,NMAXNCLUSBAS
         DMLIP(1,I)=DMLIP(2,I)
         DMLIR(1,I)=DMLIR(2,I)
         NEW_MASAP(1,I)=NEW_MASAP(2,I)
       END DO

      NBASPART_PLOT=PARTIRED
!$OMP PARALLEL DO SHARED(NBASPART_PLOT,NHOST,HHOST),
!$OMP+            PRIVATE(I)
      DO I=1,NBASPART_PLOT
        NHOST(1,:,I)=NHOST(2,:,I)
        HHOST(1,:,:,I)=HHOST(2,:,:,I)
      END DO

       NBASPART_PLOT=NNCLUS(2)
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,DAD,NDAD),
!$OMP+            PRIVATE(I)
        DO I=1,NBASPART_PLOT
           NDAD(1,I)=NDAD(2,I)
c           DAD(1,1:NDAD(1,I),I)=DAD(2,1:NDAD(2,I),I)
           DAD(1,:,I)=DAD(2,:,I)
        END DO

      IF(PLOT.EQ.2) THEN
        NBASPART_PLOT=NNCLUS(1)
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,RATIO),
!$OMP+            PRIVATE(I)
        DO I=1,NBASPART_PLOT
c           RATIO(1,1:NDAD(1,I),I)=RATIO(2,1:NDAD(2,I),I)
           RATIO(1,:,I)=RATIO(2,:,I)
        END DO

      END IF

      END IF


C       WRITE(*,*)'AHORA 1 ES 2:',NNCLUS(1),ZETAS(1),TIEMPO(1),
C     &            DMPITER(1)

*      MACHACAMOS LAS VARIABLES DE ITER 2

cx        NMAXNCLUSBAS=NNCLUS(2)
         NMAXNCLUSBAS=NMAXNCLUS    !dimension original maxima
!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,NMASA,NR,
!$OMP+     CONCENTRA,ANGULARM,NCLUSRX,NCLUSRY,NCLUSRZ,
!$OMP+             VCM2,VMAXCLUS,VX,VY,VZ,
!$OMP+             VCMAX,MCMAX,RCMAX,DMPCLUS,
!$OMP+             M200,R200,IPLIP,IPLIR,REALCLUS,
!$OMP+             NLEVHAL,EIGENVAL),PRIVATE(I)
        DO I=1,NMAXNCLUSBAS
         NMASA(2,I)=0.0
         NR(2,I)=0.0
         NCLUSRX(2,I)=0.0
         NCLUSRY(2,I)=0.0
         NCLUSRZ(2,I)=0.0
         VCM2(2,I)=0.0
         VMAXCLUS(2,I)=0.0
         VX(2,I)=0.0
         VY(2,I)=0.0
         VZ(2,I)=0.0
         VCMAX(2,I)=0.0
         MCMAX(2,I)=0.0
         RCMAX(2,I)=0.0
         CONCENTRA(2,I)=0.0
         ANGULARM(2,I)=0.0
         IPLIP(2,I)=0
         IPLIR(2,I)=0
         M200(2,I)=0.0
         R200(2,I)=0.0
         REALCLUS(2,I)=0
         NLEVHAL(2,I)=0
         DMPCLUS(2,I)=0
         EIGENVAL(2,:,I)=0.0
        END DO


        IF (PLOT.GT.1) THEN

cx        NMAXNCLUSBAS=DMPITER(2)
         NMAXNCLUSBAS=PARTIRED_PLOT   !dimension original maxima
!$OMP  PARALLEL DO SHARED(NMAXNCLUSBAS,DMLIP,DMLIR,
!$OMP+             NEW_MASAP),PRIVATE(I)
        DO I=1, NMAXNCLUSBAS
         DMLIP(2,I)=0
         DMLIR(2,I)=0
         NEW_MASAP(2,I)=0.0
        END DO


      NBASPART_PLOT=PARTIRED
!$OMP PARALLEL DO SHARED(NBASPART_PLOT,NHOST,HHOST),
!$OMP+            PRIVATE(I)
       DO I=1,NBASPART_PLOT
         NHOST(2,:,I)=0
         HHOST(2,:,:,I)=0
       END DO

       NBASPART_PLOT=NMAXNCLUS
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,DAD,NDAD),
!$OMP+            PRIVATE(I)
       DO I=1,NBASPART_PLOT
         NDAD(2,I)=0
         DAD(2,:,I)=0
       END DO

       IF (PLOT.EQ.2) THEN
        NBASPART_PLOT=NMAXNCLUS
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,RATIO),
!$OMP+            PRIVATE(I)
        DO I=1,NBASPART_PLOT
          RATIO(2,:,I)=0.0
        END DO


       END IF


       END IF

       NNCLUS(2)=0
       ZETAS(2)=0.0
       TIEMPO(2)=0.0
       DMPITER(2)=0

C       WRITE(*,*)'DE MOMENTO 2 VACIO:',NNCLUS(2),ZETAS(2),TIEMPO(2),
C     &            DMPITER(2)

CCCCCCCCCCCCCCCCCCCCC
       WRITE(*,*)'END MERGER TREE', ITER-EVERY,'-',ITER
       CALL IDATE(DATE)
       CALL ITIME(TIME)
       WRITE(*,*) 'DATE=',DATE(1),'/',DATE(2),'/',DATE(3)
       WRITE(*,*) 'TIME=',TIME(1),':',TIME(2),':',TIME(3)
CCCCCCCCCCCCCCCCCCCCC

       END IF   !MARK
       END IF   !CONTAITER

*//////////////////////END WRITING////////////////////////

       END IF   !PLOT

*/////////////////////////////////////////
       END IF    !MARK

*      Si no se hace merger_tree hay que inicializar todo lo que depende de ITER!!

       IF (CONTAITER.EQ.2) THEN

        CONTAITER=0

!$OMP PARALLEL DO SHARED(NNCLUS,ZETAS,TIEMPO,NFILE,DMPITER),PRIVATE(I)
       DO I=1,NFILE
        NNCLUS(I)=0
        ZETAS(I)=0.0
        TIEMPO(I)=0.0
        DMPITER(I)=0
       END DO

       NMAXNCLUSBAS=NMAXNCLUS
!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,NMASA,NR,CONCENTRA,
!$OMP+                   ANGULARM,NCLUSRX,NCLUSRY,NCLUSRZ,
!$OMP+                   VCM2,VMAXCLUS,VX,VY,VZ,IPLIP,IPLIR,
!$OMP+                   VCMAX,MCMAX,RCMAX,DMPCLUS,M200,R200,
!$OMP+                   REALCLUS,NLEVHAL,EIGENVAL),
!$OMP+              PRIVATE(I)
       DO I=1,NMAXNCLUSBAS
        NMASA(:,I)=0.0
        NR(:,I)=0.0
        CONCENTRA(:,I)=0.0
        ANGULARM(:,I)=0.0
        NCLUSRX(:,I)=0.0
        NCLUSRY(:,I)=0.0
        NCLUSRZ(:,I)=0.0
        VCM2(:,I)=0.0
        VMAXCLUS(:,I)=0.0
        VX(:,I)=0.0
        VY(:,I)=0.0
        VZ(:,I)=0.0
        VCMAX(:,I)=0.0
        MCMAX(:,I)=0.0
        RCMAX(:,I)=0.0
        M200(:,I)=0.0
        R200(:,I)=0.0
        IPLIP(:,I)=0
        IPLIR(:,I)=0
        DMPCLUS(:,I)=0
        REALCLUS(:,I)=0
        NLEVHAL(:,I)=0
        EIGENVAL(:,:,I)=0.0
        END DO


      IF (PLOT.GT.1) THEN

      NBASPART_PLOT=PARTIRED_PLOT
!$OMP PARALLEL DO SHARED(NBASPART_PLOT,DMLIP,DMLIR,
!$OMP+                   NEW_MASAP),
!$OMP+            PRIVATE(I)
       DO I=1,NBASPART_PLOT
         DMLIP(:,I)=0
         DMLIR(:,I)=0
         NEW_MASAP(:,I)=0.0
       END DO

      NBASPART_PLOT=PARTIRED
!$OMP PARALLEL DO SHARED(NBASPART_PLOT,NHOST,HHOST),
!$OMP+            PRIVATE(I)
       DO I=1,NBASPART_PLOT
          NHOST(:,:,I)=0
          HHOST(:,:,:,I)=0
       END DO


       NBASPART_PLOT=NMAXNCLUS
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,DAD,NDAD),
!$OMP+            PRIVATE(I)
       DO I=1,NBASPART_PLOT
          NDAD(:,I)=0
          DAD(:,:,I)=0
       END DO

       IF (PLOT.EQ.2) THEN

       NBASPART_PLOT=NMAXNCLUS
!$OMP  PARALLEL DO SHARED(NBASPART_PLOT,RATIO),
!$OMP+            PRIVATE(I)
       DO I=1,NBASPART_PLOT
          RATIO(:,:,I)=0.0
       END DO


       END IF

       END IF !PLOT>1


       END IF

       END DO    !FIN DE ITER
*/////////////////////////////////////////

        IF (PLOT.GT.1) THEN
         DEALLOCATE(DMLIP,DMLIR,NEW_MASAP)
         DEALLOCATE(HHOST,NHOST)
         DEALLOCATE(DAD,NDAD)
         IF (PLOT.EQ.2) DEALLOCATE(RATIO)
        END IF


       END

**//////////////////////////////////////////////////////////////

***********************************************************************
*      SUBROUTINES IN EXTERNAL FILES                                  *
***********************************************************************
*      Grid building
       INCLUDE 'grids.f'
*      I/O filenames
       INCLUDE 'nomfile.f'
*      Routines from 'Numerical Recipes in Fortran90', Press, Teukoslky et al.
       INCLUDE 'nr.f'
*      Halo finding procedures using the grid
       INCLUDE 'haloes_grids.f'
*      Halo finding procedures using particles
       INCLUDE 'haloes_particles.f'
*      Read MASCLET outputs (can be changed for other code outputs)
       INCLUDE 'reader.f'
*      Merger tree routines
       INCLUDE 'merger_tree.f'
