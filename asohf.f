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
*           If REALCLUS(I)= -1  ------> HALO
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
       INTEGER NPART_ESP(0:N_ESP-1)

       INTEGER CONTA2(NAMRX,NAMRY,NAMRZ,NPALEV) !this will be deleted
       INTEGER SOLAP(NAMRX,NAMRY,NAMRZ,NPALEV)

       INTEGER ORIPA2(PARTIRED)
       COMMON /PUNTEROS/ ORIPA2

       INTEGER II,JJ,KK,KK1,KK2,KKK2,ICEN(3),IR2,III,KK3,KK4
       INTEGER IX1,IX2,N1,N2,N3,NTOT,VAR,IX3
       INTEGER NFILE,FIRST,EVERY,LAST
       REAL*4 PRUEBAX,PRUEBAY,PRUEBAZ,DELTA_PRUEBA
       REAL*4 ZETA,AA2,LIM,BAS, BAS1, BAS2
       REAL*4 RRRR,R111,R222

       CHARACTER*13 FILE1, FILE4
       CHARACTER*14 FILE3, FILE7
       CHARACTER*30 FILERR3,FILERR7
       INTEGER*4 DATE(3), TIME(3), CONTAERR
       INTEGER MARK(MAXITER2),NFILE2,IFI2

       INTEGER VECINO(NPALEV,NPALEV),NVECI(NPALEV),CEN(1),CEL
       INTEGER NV, NV2,FLAG,VID(NPALEV),NV3
       REAL*4 MAXIMO(NPALEV),MAXIMO2(NPALEV),CONTAP(NPALEV)
       INTEGER CEN2(1),CEL2

       REAL*4 UBAS1(NMAX,NMAY,NMAZ)
       REAL*4 RMIN,VOL,RSHELL,ESP, REF, RANT,VOL1,VOL2
       INTEGER SHELL, NSHELL,CONTA(NMAX,NMAY,NMAZ)
       INTEGER NSHELL_2,BASINT
       REAL*4 DELTA, AA, VR, AADM, MASADM, MAP
       REAL*4 CONTRASTEC, CONTRASTEX, OMEGAZ
       REAL*4 VESC2, RS, CONCEN, VMAX2, F2, FC, RR

       INTEGER NCLUS,ICMIN,ICMAX
       INTEGER KONTA3,KONTA,KONTACAPA,KONTA1,KONTA2
       REAL*4 DIS, VKK, BASMAS, MASAKK
       REAL*4 CMX, CMY, CMZ, SOLMAS, MASAMIN, DELTA2
       REAL*4 VCMX, VCMY, VCMZ, VCM,VVV2, MASA2

       REAL*8 CMX8,CMY8,CMZ8,VCMX8,VCMY8,VCMZ8

*      ---HALOS AND SUBHALOS---
       REAL*4 MASA(MAXNCLUS), RADIO(MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
       INTEGER PATCHCLUS(MAXNCLUS)

       INTEGER REALCLUS(MAXNCLUS)
       INTEGER HALBORDERS(MAXNCLUS)

       REAL*4 CONCENTRA(NMAXNCLUS)
       REAL*4 ANGULARM(NMAXNCLUS)
       REAL*4 VMAXCLUS(NMAXNCLUS)
       REAL*4 VCM2(NMAXNCLUS)

       INTEGER IPLIP(NMAXNCLUS)
       REAL*4 VX(NMAXNCLUS)
       REAL*4 VY(NMAXNCLUS)
       REAL*4 VZ(NMAXNCLUS)

       REAL*4 VCMAX(NMAXNCLUS)
       REAL*4 MCMAX(NMAXNCLUS)
       REAL*4 RCMAX(NMAXNCLUS)

       REAL*4 M200C(NMAXNCLUS),R200C(NMAXNCLUS)
       REAL*4 M500C(NMAXNCLUS),R500C(NMAXNCLUS)
       REAL*4 M2500C(NMAXNCLUS),R2500C(NMAXNCLUS)
       REAL*4 M200M(NMAXNCLUS),R200M(NMAXNCLUS)
       REAL*4 M500M(NMAXNCLUS),R500M(NMAXNCLUS)
       REAL*4 M2500M(NMAXNCLUS),R2500M(NMAXNCLUS)

*      ---PARTICULAS E ITERACIONES---
       INTEGER DMPCLUS(NMAXNCLUS)

*      ---SUBSTRUCTURE---
       INTEGER LEVHAL(MAXNCLUS),NHALLEV(0:NLEVELS)
       INTEGER SUBHALOS(NMAXNCLUS)

*      ---HALO SHAPE---
       REAL*4 INERTIA(3,3),EIGENVAL(3,NMAXNCLUS)
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

       INTEGER PARCHLIM       ! =0: no limit patches per level, != 0: do limit
       INTEGER MPAPOLEV(NLEVELS)

       INTEGER REFINE_THR,MIN_PATCHSIZE,INTERP_DEGREE
       INTEGER BOR,BORAMR,BOR_OVLP
       REAL MINFRAC_REFINABLE,VOL_SOLAP_LOW

*      ---PARALLEL---
       INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
       COMMON /PROCESADORES/ NUM

       INTEGER NMAXNCLUSBAS,PABAS,NPBAS,NLEVBAS,NUMPARTBAS,NBASPART_PLOT

       INTEGER, ALLOCATABLE:: IP_PARAL(:,:)
       INTEGER, ALLOCATABLE:: IR_PARAL(:,:)
       REAL*4, ALLOCATABLE:: MASAP_PARAL(:,:)
       INTEGER NN, IP, NH, IP2

       INTEGER CR0AMR(NMAX,NMAY,NMAZ)
       INTEGER CR0AMR11(NAMRX,NAMRY,NAMRZ,NPALEV)

**************************************************************
*      OPENING FILES
**************************************************************
       OPEN(1,FILE='./input_files/asohf.dat',
     &              STATUS='UNKNOWN',ACTION='READ')


*      READING INITIAL DATA
****************************************************
*     NX,NY,NZ < or = NMAX,NMAY,NMAZ               *
****************************************************
       READ(1,*) !***********************************************************************
       READ(1,*) !*                         ASOHF PARAMETERS FILE                       *
       READ(1,*) !***********************************************************************
       READ(1,*) !*       General parameters block                                      *
       READ(1,*) !***********************************************************************
       READ(1,*) !Files: first, last, every -------------------------------------------->
       READ(1,*) FIRST,LAST,EVERY
       READ(1,*) !Cells per direction (NX,NY,NZ) --------------------------------------->
       READ(1,*) NX,NY,NZ
       READ(1,*) !GAS,DM particles (all levels) ---------------------------------------->
       READ(1,*) N_GAS,N_DM
       READ(1,*) !Hubble constant (h), omega matter ------------------------------------>
       READ(1,*) ACHE,OMEGA0
       READ(1,*) !Initial redshift, box size (Mpc) ------------------------------------->
       READ(1,*) ZI,LADO0
       READ(1,*) !Parallel(=1),serial(=0)/ Number of processors ------------------------>
       READ(1,*) FLAG_PARALLEL,NUM
       READ(1,*) !Reading flags: FLAG_SA,FLAG_MASCLET,FLAG_GAS ------------------------->
       READ(1,*) FLAG_SA,FLAG_MASCLET,FLAG_GAS
       READ(1,*) !***********************************************************************
       READ(1,*) !*       Mesh building parameters block                                *
       READ(1,*) !***********************************************************************
       READ(1,*) !Levels for the mesh (stand-alone) ------------------------------------>
       READ(1,*) NL
       IF (NL.GT.NLEVELS) THEN
        WRITE(*,*) 'Fatal ERROR: NLEVELS too small in parameters file',
     &             NL,NLEVELS
        STOP
       END IF
       READ(1,*) !PARCHLIM(=0 no limit patches/level,>0 limit) ------------------------->
       READ(1,*) PARCHLIM
       READ(1,*) !LIM=max num patches/level(needs PARCHLIM>0) -------------------------->
       IF (PARCHLIM.EQ.1) THEN
        READ(1,*) (MPAPOLEV(I),I=1,NL)
       ELSE
        READ(1,*)
       END IF
       READ(1,*) !Refinement threshold (num. part.), refinable fraction to extend ------>
       READ(1,*) REFINE_THR,MINFRAC_REFINABLE
       READ(1,*) !Minimum patch size (child cells) ------------------------------------->
       READ(1,*) MIN_PATCHSIZE
       READ(1,*) !Base grid refinement border, AMR grids refinement border ------------->
       READ(1,*) BOR,BORAMR
       READ(1,*) !Allow for addition overlap (to avoid losing signal) in the mesh ------>
       READ(1,*) BOR_OVLP
       READ(1,*) !Density interpolation kernel (1=linear, 2=quadratic) ----------------->
       READ(1,*) INTERP_DEGREE
       READ(1,*) !Variable for mesh halo finding: 1(gas), 2(dm), 0(gas+dm) ------------->
       READ(1,*) VAR
       READ(1,*) !***********************************************************************
       READ(1,*) !*       Halo finding parameters block                                 *
       READ(1,*) !***********************************************************************
       READ(1,*) !Max. reach around halos (Mpc), excluded cells in boundaries ---------->
       READ(1,*) BOUND, BORDES
       READ(1,*) !Minimum fraction of shared volume to merge (in grid search) ---------->
       READ(1,*) VOL_SOLAP_LOW
       READ(1,*) !FLAG_WDM (=1 write DM particles, =0 no) ------------------------------>
       READ(1,*) FLAG_WDM
       READ(1,*) !***********************************************************************
       READ(1,*) !*       Merger tree parameters block                                  *
       READ(1,*) !***********************************************************************
       READ(1,*) !Merger_tree: 1(no), 2(complete, with %), 3(main line) ---------------->
       READ(1,*) !PLOT


       CLOSE(1)

       N_PARTICLES=N_GAS+N_DM
       H2=ACHE

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

       NMAXNCLUSBAS=NMAXNCLUS
!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,CONCENTRA,
!$OMP+                   ANGULARM,CLUSRX,CLUSRY,CLUSRZ,
!$OMP+                   VCM2,VMAXCLUS,VX,VY,VZ),PRIVATE(I)
       DO I=1,NMAXNCLUSBAS
        CONCENTRA(I)=0.0
        ANGULARM(I)=0.0
        CLUSRX(I)=0.0
        CLUSRY(I)=0.0
        CLUSRZ(I)=0.0
        VCM2(I)=0.0
        VMAXCLUS(I)=0.0

        VX(I)=0.0
        VY(I)=0.0
        VZ(I)=0.0
       END DO

!$OMP PARALLEL DO SHARED(NMAXNCLUSBAS,VCMAX,MCMAX,RCMAX,DMPCLUS,M200C,
!$OMP+                   M500C,M2500C,M200M,M500M,M2500M,R200C,R500C,
!$OMP+                   R2500C,R200M,R500M,R2500M,IPLIP,REALCLUS,
!$OMP+                   LEVHAL,EIGENVAL),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
       DO I=1,NMAXNCLUSBAS
        VCMAX(I)=0.0
        MCMAX(I)=0.0
        RCMAX(I)=0.0
        M200C(I)=0.0
        M500C(I)=0.0
        M2500C(I)=0.0
        M200M(I)=0.0
        M500M(I)=0.0
        M2500M(I)=0.0
        R200C(I)=0.0
        R500C(I)=0.0
        R2500C(I)=0.0
        R200M(I)=0.0
        R500M(I)=0.0
        R2500M(I)=0.0
        IPLIP(I)=0
        LEVHAL(I)=0
        DMPCLUS(I)=0
        REALCLUS(I)=0    !de momento no hay halos
        EIGENVAL(:,I)=0.0
       END DO

       MARK(1:NFILE2)=0

       INERTIA=0.0

*///////// MAIN LOOP (ITERATIONS) /////////
*//////////////////////////////////////////
       DO IFI2=1, NFILE2                 !/
*//////////////////////////////////////////
*//////////////////////////////////////////

        ITER=FIRST+EVERY*(IFI2-1)

        WRITE(*,*)'STARTING ITER', ITER, IFI2

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
!$OMP+                   CLUSRX,CLUSRY,CLUSRZ,LEVHAL),
!$OMP+            PRIVATE(I)
        DO I=1,NMAXNCLUSBAS
         CLUSRX(I)=0.0
         CLUSRY(I)=0.0
         CLUSRZ(I)=0.0
         MASA(I)=0.0
         RADIO(I)=0.0
         LEVHAL(I)=0
        END DO

        NCLUS=0
        CMX=0.0
        CMY=0.0
        CMZ=0.0
        MASA2=0.0
        ROTE=0.0
        RETE=0.0

        SUBHALOS=0

!$OMP PARALLEL DO SHARED(PABAS,U2DM,U3DM,U4DM,RXPA,RYPA,RZPA,
!$OMP+                   MASAP,ORIPA2),
!$OMP+            PRIVATE(I)
        DO I=1,PABAS
         U2DM(I)=0.0
         U3DM(I)=0.0
         U4DM(I)=0.0
         RXPA(I)=0.0
         RYPA(I)=0.0
         RZPA(I)=0.0
         MASAP(I)=0.0
         ORIPA2(I)=0
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

        CALL SORT_DM_PARTICLES(U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &                         N_DM,NPART_ESP)

!      FIX THIS, REMOVE NPART (USELESS) FROM EVERYWHERE
       !NPART=0
       !NPART=N_DM
       ! for now, we will leave it like this for this to work temporarily
        NPART=0
        NPART(0)=N_PARTICLES

        ! Background cosmology variables
        ROTE=RODO*(1.0+ZETA)**3
        RETE=RE0/(1.0+ZETA)

        WRITE(*,*)'***********************'
        WRITE(*,*)'***** MESHRENOEF ******'
        WRITE(*,*)'***********************'

        IF (NL.GT.0) THEN
         WRITE(*,*)'==== Building the grid...', ITER, NL
         CALL CREATE_MESH(ITER,NX,NY,NZ,NL,NPATCH,PARE,PATCHNX,PATCHNY,
     &                    PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,
     &                    PATCHRZ,RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,
     &                    N_PARTICLES,N_DM,N_GAS,LADO0,T,ZETA,
     &                    REFINE_THR,MIN_PATCHSIZE,MINFRAC_REFINABLE,
     &                    BOR,BORAMR,BOR_OVLP)
         WRITE(*,*)'==== END building the grid...', ITER, NL
        END IF

c        WRITE(*,*) 'TSC density interpolation, levels min,max:',0,NL_TSC
c        CALL INTERPOLATE_DENSITY(ITER,NX,NY,NZ,NL_TSC,NPATCH,PARE,
c     &           PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
c     &           PATCHRX,PATCHRY,PATCHRZ,RXPA,RYPA,RZPA,MASAP,
c     &           N_PARTICLES,N_DM,N_GAS,LADO0,T,ZETA,NPART_ESP)
c
c        IF (NL_TSC.LT.NL) THEN
c         WRITE(*,*) 'Smooth density interpolation, levels min,max:',
c     &              NL_TSC+1,NL
c         CALL INTERPOLATE_DENSITY_KERNEL(ITER,NX,NY,NZ,NL_TSC,
c     &            NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,
c     &            PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,RXPA,RYPA,RZPA,
c     &            MASAP,N_PARTICLES,N_DM,N_GAS,LADO0,T,ZETA,NPART_ESP)
c        END IF

       CALL DENSITY(ITER,NX,NY,NZ,NL,NPATCH,PARE,PATCHNX,PATCHNY,
     &              PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,
     &              PATCHRZ,RXPA,RYPA,RZPA,MASAP,N_PARTICLES,N_DM,
     &              N_GAS,LADO0,T,ZETA,NPART_ESP,INTERP_DEGREE)

       WRITE(*,*)'***************************'
       WRITE(*,*)'***** END MESHRENOEF ******'
       WRITE(*,*)'***************************'

       END IF

c       DO IR=1,NL
c       DO I=SUM(NPATCH(0:IR-1))+1,SUM(NPATCH(0:IR))
c        N1=PATCHNX(I)
c        N2=PATCHNY(I)
c        N3=PATCHNZ(I)
c        WRITE(*,*) IR,I,SUM(U11(1:N1,1:N2,1:N3,I))/(N1*N2*N3),
c     &             MAXVAL(U11(1:N1,1:N2,1:N3,I))
c       END DO
c       END DO

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
       CALL COMPUTE_CR0AMR(NL,NX,NY,NZ,NPATCH,PARE,PATCHNX,PATCHNY,
     &                     PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                     PATCHRY,PATCHRZ,CR0AMR,CR0AMR11,LADO0)
c
       OPEN(99,FILE='output_files/density_asohf',STATUS='UNKNOWN',
     &      FORM='UNFORMATTED')
         write(99) (((u1(ix,jy,kz),ix=1,nx),jy=1,ny),kz=1,nz)
         write(99) (((cr0amr(ix,jy,kz),ix=1,nx),jy=1,ny),kz=1,nz)
         do i=1,sum(npatch(0:nl))
          n1=patchnx(i)
          n2=patchny(i)
          n3=patchnz(i)
          write(99) (((u11(ix,jy,kz,i),ix=1,n1),jy=1,n2),kz=1,n3)
          write(99) (((cr0amr11(ix,jy,kz,i),ix=1,n1),jy=1,n2),kz=1,n3)
          write(99) (((solap(ix,jy,kz,i),ix=1,n1),jy=1,n2),kz=1,n3)
         end do
       CLOSE(99)
*********************************************************************

c       CALL CLEAN_OVERLAPS(NL,NPATCH,PATCHNX,PATCHNY,PATCHNZ,SOLAP,
c     &                     U11)

**********************************************************************
******************************HALO FINDER*****************************
**********************************************************************

       WRITE(*,*) '************************************************'
       WRITE(*,*) '         HALO FINDING                           '
       WRITE(*,*) '************************************************'


**********************************************************
*      Looking for candidate haloes at the AMR levels
**********************************************************

       CALL HALOFIND_GRID(NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                    PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                    PATCHRY,PATCHRZ,PARE,NCLUS,MASA,RADIO,
     &                    CLUSRX,CLUSRY,CLUSRZ,REALCLUS,LEVHAL,
     &                    NHALLEV,BOUND,CONTRASTEC,RODO,
     &                    SOLAP,VECINO,NVECI,CR0AMR,CR0AMR11,PATCHCLUS,
     &                    VOL_SOLAP_LOW)

       open(55, file='./output_files/haloesgrids.res', status='unknown')
       do i=1,nclus
        write(55,*) clusrx(i),clusry(i),clusrz(i),radio(i),masa(i),
     &              levhal(i), realclus(i), patchclus(i)
       end do
       close(55)

*******************************************************
*      SORTING OUT ALL THE CLUSTERS
*******************************************************

       NHALLEV(:)=0
       J=0
       DO I=1,NCLUS
        IF (REALCLUS(I).EQ.0) CYCLE
        J=J+1
        CLUSRX(J)=CLUSRX(I)
        CLUSRY(J)=CLUSRY(I)
        CLUSRZ(J)=CLUSRZ(I)
        RADIO(J)=RADIO(I)
        MASA(J)=MASA(I)*UM ! From now on, in Solar Masses
        LEVHAL(J)=LEVHAL(I)
        PATCHCLUS(J)=PATCHCLUS(I)
        REALCLUS(J)=REALCLUS(I)
        NHALLEV(LEVHAL(J))=NHALLEV(LEVHAL(J))+1
       END DO
       NCLUS=J

       WRITE(*,*)'MASSES: MIN, MAX AND MEAN=', MINVAL(MASA(1:NCLUS)),
     &            MAXVAL(MASA(1:NCLUS)), SUM(MASA(1:NCLUS))/NCLUS
       WRITE(*,*) 'NCLUS=', NCLUS
       WRITE(*,*) 'MEAN MASS',SUM(MASA(1:NCLUS))/NCLUS

***************************************************************
**     HALOES AT THE EDGES OF THE BOX (JUST FOR CAUTION!)     *
***************************************************************

       IF (BORDES.EQ.1) THEN
        CALL HALOES_BORDER(NCLUS,CLUSRX,CLUSRY,CLUSRZ,RADIO,LADO0,
     &                     HALBORDERS)
        WRITE(*,*) 'Haloes close to the box borders:',
     &             SUM(HALBORDERS(1:NCLUS))
       END IF

************************************************************
**     Eliminating POOR haloes (less than a minimum number of particles)
**     (We start here to work with partciles for the 1st time)
************************************************************

       NUMPARTBAS=NUMPART
       CALL PRUNE_POOR_HALOES(NCLUS,CLUSRX,CLUSRY,CLUSRZ,RADIO,
     &                        REALCLUS,RXPA,RYPA,RZPA,N_DM,NUMPARTBAS,
     &                        DMPCLUS,1.0)

*********************************************************
*      CHECKING....
*********************************************************
       WRITE(*,*) '---------------------------------'
       WRITE(*,*) 'CHECKING HIERARCHY:'
       WRITE(*,*) '---------------------------------'
       KONTA2=0
       KONTA2=COUNT(REALCLUS(1:NCLUS).EQ.-1)
       WRITE(*,*) 'REAL HALOS ------------------>', KONTA2
       KONTA2=0
       KONTA2=COUNT(REALCLUS(1:NCLUS).EQ.0)
       WRITE(*,*) 'REMOVED HALOS --------------->', KONTA2
       KONTA2=0
       KONTA2=COUNT(REALCLUS(1:NCLUS).GT.0)
       WRITE(*,*) 'SUBSTRUCTURE candidates ----->', KONTA2
       WRITE(*,*) '---------------------------------'
*********************************************************

************************************************************
*      REFINING REAL HALOES WITH THE DM PARTICLES ONLY     *
************************************************************

       WRITE(*,*)'=================================='
       WRITE(*,*)'Refining with DM particles...'
       WRITE(*,*)'=================================='

       CALL HALOFIND_PARTICLES(NL,NCLUS,MASA,RADIO,CLUSRX,CLUSRY,
     &      CLUSRZ,REALCLUS,CONCENTRA,ANGULARM,VMAXCLUS,VCM2,IPLIP,
     &      VX,VCMAX,MCMAX,RCMAX,M200C,M500C,M2500C,M200M,M500M,M2500M,
     &      R200C,R500C,R2500C,R200M,R500M,R2500M,DMPCLUS,LEVHAL,
     &      EIGENVAL,N_DM,RXPA,RYPA,RZPA,MASAP,U2DM,U3DM,U4DM,ORIPA2,
     &      CONTRASTEC,OMEGAZ,UM,UV,F2)


*************************************************
******** GENERAL CHECKING ***********************
*************************************************

       WRITE(*,*)'HALOES WITHOUT MASS=',
     &        COUNT(MASA(1:NCLUS).LE.0.0)

       WRITE(*,*)'HALOES WITH MASS=0',
     &        COUNT(MASA(1:NCLUS).EQ.0.0)
*       WRITE(*,*)'KKKONTA=',KKKONTA

       WRITE(*,*) 'Total number of particles within halos:',
     &            SUM(DMPCLUS(1:NCLUS))


       WRITE(*,*)'After refining with DM particles...'
       WRITE(*,*)'===================================='
       DO I=0,NL
       WRITE(*,*)'Haloes at level ', I,' =',
     &            COUNT(LEVHAL(1:NCLUS).EQ.I
     &            .AND.REALCLUS(1:NCLUS).NE.0),
     &            COUNT(REALCLUS(1:NCLUS).EQ.-1)
       END DO
       WRITE(*,*)'===================================='

*************************************************
*************************************************

c       WRITE(*,*) 'ESTIMATION HALOES_1 (AFTER UNBINDING)'
c       WRITE(*,*)'=================================='
c       DO I=1, NCLUS
c       IF (REALCLUS(I).NE.0) THEN
c       WRITE(*,*) I, MASA(I),RADIO(I),CLUSRX(I),
c     &            CLUSRY(I),CLUSRZ(I)
c       END IF
c       END DO
c      WRITE(*,*)'=================================='



*      REPETIMOS ESTO OTRA VEZ
*STEP 1) (2)************************************
************REMOVING POOR HALOES***************
************************************************
       KONTA2=0
       NUMPARTBAS=NUMPART
!$OMP  PARALLEL DO SHARED(NCLUS,DMPCLUS,NUMPARTBAS,
!$OMP+             REALCLUS),PRIVATE(I,KK_ENTERO),
!$OMP+             REDUCTION(+:KONTA2)
       DO I=1, NCLUS

       KK_ENTERO=0
       KK_ENTERO=DMPCLUS(I)

       IF (KK_ENTERO.LT.NUMPARTBAS) THEN

       REALCLUS(I)=0
       KONTA2=KONTA2+1

       END IF

       END DO

       WRITE(*,*)'RE-CHECKING POOR HALOS----->', KONTA2


*PASO 2)****************************************
************RUBBISH (desbordamientos)***********
************************************************

       KONTA2=0

       DO I=1, NCLUS

       IF (REALCLUS(I).NE.0) THEN

       DO J=1, NCLUS


       IF (REALCLUS(J).NE.0.AND.LEVHAL(J).GT.
     &      LEVHAL(I)) THEN

       DIS=0.0
       DIS=SQRT((CLUSRX(I)-CLUSRX(J))**2+
     &          (CLUSRY(I)-CLUSRY(J))**2+
     &          (CLUSRZ(I)-CLUSRZ(J))**2)

       A1=0.0
       A1=MIN(MASA(I),MASA(J))/MAX(MASA(I),MASA(J))

       A2=0.0
       IF (MIN(ABS(VCM2(I)),ABS(VCM2(J))).NE.0.0) THEN
       A2=(ABS(VCM2(I)-VCM2(J)))/
     &                            MAX(ABS(VCM2(I)),ABS(VCM2(J)))
       END IF

       A3=0.0
       A3=MIN(RADIO(I),RADIO(J))

*       IF (DIS.LT.1.01*A3.AND.A1.GT.0.2.AND.A2.LT.3.0.OR.A1.GT.0.6) THEN
       IF (DIS.LT.1.01*A3.AND.A1.GT.0.2.AND.A2.LT.3.0) THEN
*       IF (DIS.LT.1.01*A3.AND.A1.GT.0.2.AND.A2.LT.5.0) THEN
*       IF (DIS.LT.1.01*A3.AND.A1.GT.0.2) THEN

       IF (MASA(I).GT.MASA(J)) THEN

       REALCLUS(J)=0
       KONTA2=KONTA2+1

       ELSE

       REALCLUS(I)=0
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
       DO I=1, NCLUS

       IF (REALCLUS(I).EQ.-1) THEN

       CONCEN=0.0
       RS=0.0
       FC=0.0
       VMAX2=0.0

       IF (MASA(I).GT.0.0) THEN
       CONCEN=124.0*((MASA(I)*(ACHE/3.66D-3)**(-1))**(-0.084))
       RS=RADIO(I)/CONCEN
       END IF

       CONCENTRA(I)=CONCEN

       FC=LOG(1.0+CONCEN)-(CONCEN/(1.0+CONCEN))

       IF (FC.GT.0.0) THEN
       VMAX2=(CGR*(MASA(I)/9.1717E18)*F2)/(RS*RETE*2.0*FC)
       END IF


       DO J=1, NCLUS

       DIS=0.0
       DIS=SQRT((CLUSRX(I)-CLUSRX(J))**2+
     &          (CLUSRY(I)-CLUSRY(J))**2+
     &          (CLUSRZ(I)-CLUSRZ(J))**2)

       VVV2=0.0
       VVV2=(VCM2(I)-VCM2(J))**2

       VESC2=0.0
       IF (RS.NE.0.AND.DIS*F2/RS.NE.0) THEN
       VESC2=(4.0*VMAX2*LOG(1.0+DIS/RS))/(DIS*F2/RS)
       END IF

       A1=0.0
       A1=MIN(MASA(I),MASA(J))/MAX(MASA(I),MASA(J))

       IF (REALCLUS(J)==-1) THEN

       IF(LEVHAL(J).GT.LEVHAL(I).AND.
     &    DIS.LE.1.0*RADIO(I).AND.A1.LE.0.2.AND.VVV2.LE.VESC2) THEN


       REALCLUS(J)=I
       SUBHALOS(I)=SUBHALOS(I)+1
       IF (SUBHALOS(I).GT.NMAXSUB) THEN
       WRITE(*,*)'WARNING: DEMASIASOS SUBHALOS!!', ITER, I, SUBHALOS
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
     &            COUNT(REALCLUS(1:NCLUS).EQ.-1)


*********************************************************
*      CHECKING....2
*********************************************************
       WRITE(*,*) '---------------------------------'
       WRITE(*,*) 'RE-CHECKING HIERARCHY:'
       WRITE(*,*) '---------------------------------'
       KONTA2=0
       KONTA2=COUNT(REALCLUS(1:NCLUS).EQ.-1)
       WRITE(*,*) 'REAL HALOS----->', KONTA2
       KONTA2=0
       KONTA2=COUNT(REALCLUS(1:NCLUS).EQ.0)
       WRITE(*,*) 'REMOVED HALOS----->', KONTA2
       KONTA2=0
       KONTA2=COUNT(REALCLUS(1:NCLUS).GT.0)
       WRITE(*,*) 'SUBSTRUCTURE----->', KONTA2
       WRITE(*,*) '---------------------------------'

       WRITE(*,*)'At the end...'
       WRITE(*,*)'=================================='
       DO I=0,NL
       WRITE(*,*)'Haloes at level ', I,' =',
     &            COUNT(LEVHAL(1:NCLUS).EQ.I.
     &            AND.REALCLUS(1:NCLUS).NE.0),
     &            COUNT(REALCLUS(1:NCLUS).EQ.-1)
       END DO
       WRITE(*,*)'=================================='


****************************************************
****************************************************
****************************************************


*************************************************
*************************************************
*===================Families===============
       KONTA2=0
c       KONTA2=COUNT(REALCLUS(1:NCLUS).EQ.-1)
       KONTA2=COUNT(REALCLUS(1:NCLUS).NE.0)
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
       WRITE(3,*) ITER, NCLUS, KONTA2, ZETA
       WRITE(3,*) '************************************************'
       KONTA2=0

       DO I=1, NCLUS

       IF (REALCLUS(I).NE.0) THEN

         WRITE(3,135) I,CLUSRX(I),CLUSRY(I),CLUSRZ(I),
     &         MASA(I),RADIO(I),DMPCLUS(I),
     &         REALCLUS(I), LEVHAL(I),SUBHALOS(I),
     &         EIGENVAL(1,I),EIGENVAL(2,I),EIGENVAL(3,I),
     &         VCM2(I)*UV,CONCENTRA(I),ANGULARM(I)*1.0E14,
     &         VCMAX(I),MCMAX(I),RCMAX(I),
     &         R200M(I),M200M(I),R200C(I),M200C(I),
     &         R500M(I),M500M(I),R500C(I),M500C(I),
     &         R2500M(I),M2500M(I),R2500C(I),M2500C(I),
     &         VX(I)*UV,VY(I)*UV,VZ(I)*UV


       IF (FLAG_WDM.EQ.1) THEN

       WRITE(4) I,REALCLUS(I),DMPCLUS(I)

       KK2=0
       KKK2=0
       KK2=SUM(DMPCLUS(1:I-1))+1
       KKK2=SUM(DMPCLUS(1:I))

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

*      Si no se hace merger_tree hay que inicializar todo lo que depende de ITER!!

       END DO    !FIN DE ITER
*/////////////////////////////////////////

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
       !INCLUDE 'merger_tree.f'
