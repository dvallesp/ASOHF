**********************************************************************
       PROGRAM ASOHF
***********************************************************************
*      ASOHF IS AN ADAPTIVE SPHERICAL OVERDENSITY HALO FINDER.
***********************************************************************
*      For further details, check:
*      - Planelles & Quilis, 2010. A&A, 519:A94
*      - Vallés-Pérez, Planelles & Quilis, 2022. A&A submitted
*      - Code documentation: https://asohf.github.io/
***********************************************************************

       USE PARTICLES
       IMPLICIT NONE

       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER I,J,K,I2,J2
       INTEGER NX,NY,NZ,ITER,NDXYZ
       REAL*4 T

       REAL*4  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
       COMMON /GRID/   RADX,RADY,RADZ

       REAL*4  RX(0:NAMRX+1,NPALEV),RY(0:NAMRX+1,NPALEV),
     &         RZ(0:NAMRX+1,NPALEV)
       COMMON /GRIDAMR/ RX,RY,RZ

       REAL*4 PI,ACHE,T0,RE0
       COMMON /DOS/ACHE,T0,RE0
       REAL*4 UNTERCIO,CGR,CGR2,RODO,LADO,LADO0
       COMMON /CONS/CGR,PI
       REAL*4 OMEGA0

       REAL*4 RETE,HTE,ROTE
       COMMON /BACK/ RETE,HTE,ROTE

       REAL*4 DX,DY,DZ,HUBBLE_LITTLEH
       COMMON /ESPACIADO/ DX,DY,DZ

       REAL*4 UV, UM

       INTEGER IX,JY,KZ,NL,IR,L1,NL_INPUT
       REAL*4 RX2,RY2,RZ2,A1,A2,B1,C1,A3,A4
       REAL*4 DXPA,DYPA,DZPA

*      VARIABLES
       REAL*4 U1(NMAX,NMAY,NMAZ)
       REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
       COMMON /VARIA/ U1,U11

c       REAL*4 POT(NMAX,NMAY,NMAZ)
c       REAL*4 POT1(NAMRX,NAMRY,NAMRZ,NPALEV)

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

       INTEGER NPART_ESP(0:N_ESP-1),NDMPART_X(0:NMAX)

       INTEGER SOLAP(NAMRX,NAMRY,NAMRZ,NPALEV)

       INTEGER FIRST,LAST,EVERY,NFILE,NFILE2,IFI2
       INTEGER N1,N2,N3,VAR,KONTA2,MAX_NUM_PART
       INTEGER NCLUS
       INTEGER*4 DATE(3), TIME(3)
       REAL ZETA,CONTRASTEC,F2,MAP,OMEGAZ

*      ---HALOS AND SUBHALOS---
       REAL*4 MASA(MAXNCLUS), RADIO(MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
       REAL*4 CLUSRXCM(MAXNCLUS),CLUSRYCM(MAXNCLUS),CLUSRZCM(MAXNCLUS)
       REAL*4 MSUB(MAXNCLUS),RSUB(MAXNCLUS)
       INTEGER PATCHCLUS(MAXNCLUS),REALCLUS(MAXNCLUS)
       INTEGER HALBORDERS(MAXNCLUS),LEVHAL(MAXNCLUS)
       REAL*4 CONCENTRA(NMAXNCLUS)
       REAL*4 ANGULARM(3,NMAXNCLUS)
       REAL*4 VMAXCLUS(NMAXNCLUS)
       INTEGER IPLIP(NMAXNCLUS)
       REAL*4 VX(NMAXNCLUS),VY(NMAXNCLUS),VZ(NMAXNCLUS)
       REAL*4 MEAN_VR(NMAXNCLUS)
       REAL*4 VCMAX(NMAXNCLUS),MCMAX(NMAXNCLUS),RCMAX(NMAXNCLUS)
       REAL*4 M200C(NMAXNCLUS),R200C(NMAXNCLUS)
       REAL*4 M500C(NMAXNCLUS),R500C(NMAXNCLUS)
       REAL*4 M2500C(NMAXNCLUS),R2500C(NMAXNCLUS)
       REAL*4 M200M(NMAXNCLUS),R200M(NMAXNCLUS)
       REAL*4 M500M(NMAXNCLUS),R500M(NMAXNCLUS)
       REAL*4 M2500M(NMAXNCLUS),R2500M(NMAXNCLUS)
       REAL*4 PROFILES(NBINS,2,NMAXNCLUS)
       REAL*4 VELOCITY_DISPERSION(NMAXNCLUS)
       REAL*4 RMAXSIGMA(NMAXNCLUS),MMAXSIGMA(NMAXNCLUS)
       REAL*4 KINETIC_E(NMAXNCLUS),POTENTIAL_E(NMAXNCLUS)
       REAL*4 FSUB(NMAXNCLUS) ! fraction of mass in substructures
       INTEGER NSUBS(NMAXNCLUS) ! number of substructures
       INTEGER INDCS_PARTICLES_PER_HALO(2,NMAXNCLUS)
       INTEGER DMPCLUS(MAXNCLUS)
       INTEGER NHALLEV(0:NLEVELS),SUBS_LEV(0:NLEVELS)
       INTEGER SUBHALOS(NMAXNCLUS)
       REAL*4 EIGENVAL(3,NMAXNCLUS)
       REAL*4 INERTIA_TENSOR(6,NMAXNCLUS)

*      ---STAND-ALONE HALO FINDER---
       INTEGER FLAG_SA,FLAG_GAS,FLAG_MASCLET,FLAG_WDM
       INTEGER N_DM,N_PARTICLES,N_ST,N_GAS,IR_KERN_STARS,MIN_NUM_PART
       INTEGER SPLIT_SPECIES,BORDES,PARCHLIM,MIN_NUM_PART_ST
       INTEGER MPAPOLEV(NLEVELS),MAX_PART_DSUM
       INTEGER REFINE_THR,MIN_PATCHSIZE,INTERP_DEGREE
       INTEGER BOR,BORAMR,BOR_OVLP
       INTEGER LOW1,LOW2
       REAL MINFRAC_REFINABLE,VOL_SOLAP_LOW,BOUND,FDM
       REAL XLDOM,XRDOM,YLDOM,YRDOM,ZLDOM,ZRDOM
       REAL CIO_MASS,CIO_SPEED,CIO_LENGTH,CIO_ALPHA,CIO_XC,CIO_YC,CIO_ZC
       COMMON /CONV_IO/ CIO_MASS,CIO_SPEED,CIO_LENGTH,CIO_ALPHA,CIO_XC,
     &                  CIO_YC,CIO_ZC
       REAL CIO_XC0,CIO_YC0,CIO_ZC0,LADO_BKP,LADO0_BKP

       INTEGER DO_DOMDECOMP
       REAL DDXL,DDXR,DDYL,DDYR,DDZL,DDZR
       COMMON /DOM_DECOMP/ DO_DOMDECOMP,DDXL,DDXR,DDYL,DDYR,DDZL,DDZR

*      ---PARALLEL---
       INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR,FLAG_PARALLEL
       COMMON /PROCESADORES/ NUM

       INTEGER FLAG_SUBS,FLAG_CENTRAL,DO_COMPUTE_ENERGIES,FLAG_STELLAR
       INTEGER FW1,FW2,FW3,FW4,FW5
       REAL STPAR_FACT_INC,STPAR_MAX_DIST,STPAR_MIN_OVERDENS,
     &      STPAR_MAX_R_PHYS

       INTEGER CR0AMR(NMAX,NMAY,NMAZ)
       INTEGER CR0AMR11(NAMRX,NAMRY,NAMRZ,NPALEV)

       CHARACTER*5 ITER_STRING

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
       READ(1,*) !Hubble constant (h), omega matter, fraction of DM to total mass ------>
       READ(1,*) ACHE,OMEGA0,FDM
       READ(1,*) !Max box sizelength (in length units specified below) ----------------->
       READ(1,*) LADO0
       READ(1,*) !Reading: IS_MASCLET (=0, no; =1, yes), GENERIC READER (see docs) ----->
       READ(1,*) FLAG_MASCLET,FLAG_SA
       READ(1,*) !Output flags: grid_asohf,density,haloes_grids,subs_grids,subs_part --->
       READ(1,*) FW1,FW2,FW3,FW4,FW5
       READ(1,*) !Input units: MASS (Msun; <0 for cMpc/h), LENGTH (cMpc; <0 for cMpc/h),
       READ(1,*) ! SPEED (km/s), ALPHA (v_input = a^alpha dx/dt; 1 is peculiar vel.) --->
       READ(1,*) CIO_MASS,CIO_LENGTH,CIO_SPEED,CIO_ALPHA
       READ(1,*) !Input domain (in input length units; x1,x2,y1,y2,z1,z2) -------------->
       READ(1,*) XLDOM,XRDOM,YLDOM,YRDOM,ZLDOM,ZRDOM
       READ(1,*) !***********************************************************************
       READ(1,*) !*       Domain decompose                                              *
       READ(1,*) !***********************************************************************
       READ(1,*) !Keep only particles inside a given domain (=0, no; =1, yes) ---------->
       READ(1,*) DO_DOMDECOMP
       READ(1,*) !Domain to keep particles (in input length units; x1,x2,y1,y2,z1,z2) -->
       READ(1,*) DDXL,DDXR,DDYL,DDYR,DDZL,DDZR
       READ(1,*) !***********************************************************************
       READ(1,*) !*       Mesh building parameters block                                *
       READ(1,*) !***********************************************************************
       READ(1,*) !Levels for the mesh (stand-alone) ------------------------------------>
       READ(1,*) NL_INPUT
       NL=NL_INPUT
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
       READ(1,*) !Variable for mesh halo finding: 1(dm), 2(dm+stars) ------------------->
       READ(1,*) VAR
       READ(1,*) !Kernel level for stars (if VAR=2) ------------------------------------>
       READ(1,*) IR_KERN_STARS
       READ(1,*) !Particle especies (0=there are different mass particles, 1=equal mass
       READ(1,*) !particles, use local density, 2=equal mass particles, do nothing) --->
       READ(1,*) SPLIT_SPECIES
       READ(1,*) !***********************************************************************
       READ(1,*) !*       Halo finding parameters block                                 *
       READ(1,*) !***********************************************************************
       READ(1,*) !Max. reach around halos (cMpc), excluded cells in boundaries --------->
       READ(1,*) BOUND, BORDES
       READ(1,*) !Minimum fraction of shared volume to merge (in grid search) ---------->
       READ(1,*) VOL_SOLAP_LOW
       READ(1,*) !FLAG_WDM (=1 write DM particles, =0 no) ------------------------------>
       READ(1,*) FLAG_WDM
       READ(1,*) !Search for substructure (=1 yes, =0 no) ------------------------------>
       READ(1,*) FLAG_SUBS
       !READ(1,*) !Search for cores (max sigma_v of bound particles) (=1 yes, =0 no) --->
       !READ(1,*) FLAG_CENTRAL
       READ(1,*) !Compute energies (=1 yes, =0 no), max num of part. for direct sum ---->
       READ(1,*) DO_COMPUTE_ENERGIES,MAX_PART_DSUM
       READ(1,*) !Minimum number of particles per halo --------------------------------->
       READ(1,*) MIN_NUM_PART
       READ(1,*) !***********************************************************************
       READ(1,*) !*       Stellar haloes (galaxies) block                               *
       READ(1,*) !***********************************************************************
       READ(1,*) !Look for stellar haloes (=1 yes, =0 no) ------------------------------>
       READ(1,*) FLAG_STELLAR
       READ(1,*) !Minimum number of stellar particles per stellar halo ----------------->
       READ(1,*) MIN_NUM_PART_ST
       READ(1,*) !Cut stellar halo if density increases more than factor from min ------>
       READ(1,*) STPAR_FACT_INC
       READ(1,*) !Cut stellar halo if radial distance of consecutive stars > (ckpc) ---->
       READ(1,*) STPAR_MAX_DIST
       READ(1,*) !Cut stellar halo if rho_* falls below this factor of rho_B ----------->
       READ(1,*) STPAR_MIN_OVERDENS
       READ(1,*) !Cut stellar halo at a maximum (>0, physical; <0, comoving) radius of
       READ(1,*) !(kpc) ---------------------------------------------------------------->
       READ(1,*) STPAR_MAX_R_PHYS

       CLOSE(1)

       N_PARTICLES=N_DM
       HUBBLE_LITTLEH=ACHE

       IF (CIO_MASS.LT.0.) CIO_MASS=-CIO_MASS/HUBBLE_LITTLEH
       IF (CIO_LENGTH.LT.0.) CIO_LENGTH=-CIO_LENGTH/HUBBLE_LITTLEH
       LADO0=LADO0*CIO_LENGTH
       ! center of the domain (in input length units)
       CIO_XC0=0.5*(XLDOM+XRDOM)
       CIO_YC0=0.5*(YLDOM+YRDOM)
       CIO_ZC0=0.5*(ZLDOM+ZRDOM)

       STPAR_MAX_DIST=STPAR_MAX_DIST/1000.0 ! to cMpc
       STPAR_MAX_R_PHYS=STPAR_MAX_R_PHYS/1000.0 ! to Mpc
**************************************************************
*     ...PARALLEL RUNNING...
!$OMP PARALLEL SHARED(NUM)
!$OMP SINGLE
      NUM=OMP_GET_NUM_THREADS()
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL
      FLAG_PARALLEL=0
      IF (NUM.GT.1) FLAG_PARALLEL=1
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

       IF (FLAG_MASCLET.EQ.1)
     &       WRITE(*,*) 'ASOHF reading MASCLET PARTICLES...'
       IF (FLAG_MASCLET.EQ.0.AND.FLAG_SA.EQ.0)
     &       WRITE(*,*) 'ASOHF reading GENERIC PARTICLE DATA...'
       IF (FLAG_MASCLET.EQ.0.AND.FLAG_SA.EQ.1)
     &       WRITE(*,*) 'ASOHF reading GADGET UNFORMATTED FILES...'

       IF(VAR.EQ.1) WRITE(*,*) 'Analysing only DM'
       IF(VAR.EQ.2) WRITE(*,*) 'Analysing DM+stars'

       WRITE(*,*) 'Min. number of particles per halo ',MIN_NUM_PART


***************************
*      GRID BUILDER
***************************
       LADO=LADO0-(LADO0/NX)
       CALL MALLA(NX,NY,NZ,LADO)

       WRITE(*,*)
       WRITE(*,*) '************************************************'
       WRITE(*,*) '                     GRID                       '
       WRITE(*,*) '************************************************'
       WRITE(*,*) 'SIDE LENGTH=',LADO
       WRITE(*,*) 'NX,DX,RADX(1),RADX(NX)=',NX,DX,RADX(1),RADX(NX)
       WRITE(*,*) 'NUMBER OF PATCHES PER LEVEL:'
       IF (PARCHLIM.EQ.0) WRITE(*,*) '  No limited patches per level!'
       IF (PARCHLIM.NE.0) THEN
             WRITE(*,*) '  Limit patches per level=', MPAPOLEV(1)
       END IF
       WRITE(*,*)

!      Backup these variables to avoid the domain decomposition messing
!       around with them.
       LADO_BKP=LADO
       LADO0_BKP=LADO0

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

       UV=299792.458
       UM=9.1717E+18
       F2=LOG(3.0)-(2.0/3.0)

       HTE=ACHE !!!!!!!!!!!! DEBE PONER BIEN

********************************************************************
********************************************************************

       NFILE2=INT((LAST-FIRST)/EVERY) + 1
       WRITE(*,*)'Number of iterations to analise=',NFILE2

*///////// MAIN LOOP (ITERATIONS) /////////
*//////////////////////////////////////////
       DO IFI2=1, NFILE2                 !/
*//////////////////////////////////////////
*//////////////////////////////////////////

        ITER=FIRST+EVERY*(IFI2-1)

        WRITE(*,*)
        WRITE(*,*) '************************************************'
        WRITE(*,*) '************************************************'
        WRITE(*,*) '* STARTING ITER', ITER, IFI2
        WRITE(*,*) '************************************************'
        WRITE(*,*) '************************************************'
        WRITE(*,*)
        WRITE(ITER_STRING, '(I5.5)') ITER !For saving files to disk

*       Recover these quantities from the backup, may have been change
*        if DOMDECOMP is performed.
        CIO_XC=CIO_XC0
        CIO_YC=CIO_YC0
        CIO_ZC=CIO_ZC0
        LADO=LADO_BKP
        LADO0=LADO0_BKP

        CALL INIT_OUTVARS(MASA,RADIO,CLUSRX,CLUSRY,CLUSRZ,CLUSRXCM,
     &           CLUSRYCM,CLUSRZCM,MSUB,RSUB,PATCHCLUS,REALCLUS,
     &           HALBORDERS,LEVHAL,CONCENTRA,ANGULARM,VMAXCLUS,IPLIP,VX,
     &           VY,VZ,MEAN_VR,VCMAX,MCMAX,RCMAX,M200C,M500C,M2500C,
     &           R200C,R500C,R2500C,M200M,M500M,M2500M,R200C,R500C,
     &           R2500C,PROFILES,VELOCITY_DISPERSION,RMAXSIGMA,
     &           MMAXSIGMA,KINETIC_E,POTENTIAL_E,FSUB,NSUBS,
     &           INDCS_PARTICLES_PER_HALO,DMPCLUS,NHALLEV,SUBS_LEV,
     &           SUBHALOS,EIGENVAL,INERTIA_TENSOR,NCLUS)

        CALL INIT_GRIDVARS(NX,NY,NZ,PATCHNX,PATCHNY,PATCHNZ,PATCHX,
     &           PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,PARE,NPATCH,U1,
     &           U11,ROTE,RETE)


***************************************************
*     READING INPUT DATA
***************************************************

*       Reading external list of particles (either Masclet particles
*       or a general list of particles, depending on FLAG_MASCLET)
        CALL READ_AND_ALLOC_PARTICLES(FLAG_MASCLET,FLAG_SA,ITER,NX,NY,
     &       NZ,T,ZETA,N_DM,VAR,N_ST,N_PARTICLES,UV,UM,HUBBLE_LITTLEH,
     &       LADO0,LADO)

        IF (SPLIT_SPECIES.EQ.0) THEN
         CALL SORT_DM_PARTICLES(N_DM,NPART_ESP,N_ST,IR_KERN_STARS)
        ELSE IF (SPLIT_SPECIES.EQ.1) THEN
         CALL SORT_DM_PARTICLES_LOCALDENSITY(N_DM,NPART_ESP,N_ST,
     &                                       IR_KERN_STARS,RODO,RE0)
        ELSE IF (SPLIT_SPECIES.EQ.2) THEN
         NPART_ESP(0)=N_DM
         NPART_ESP(1:N_ESP-1)=0
         IF (N_ST.GT.0) THEN
          IF (IR_KERN_STARS.GT.N_ESP-1) THEN
           WRITE(*,*) 'WARNING: IR_KERN_STARS > N_ESP-1',
     &                IR_KERN_STARS,N_ESP-1
           WRITE(*,*) 'Fix N_ESP to IR_KERN_STARS+1, at least'
           STOP
          END IF
          NPART_ESP(IR_KERN_STARS)=NPART_ESP(IR_KERN_STARS)+N_ST
          WRITE(*,*) 'Stars: Of species',IR_KERN_STARS,
     &               ', no. particles:',NPART_ESP(IR_KERN_STARS)
         END IF
        END IF !(SPLIT_SPECIES.EQ.0)

        ! Background cosmology variables
        ROTE=RODO*(1.0+ZETA)**3
        RETE=RE0/(1.0+ZETA)

        WRITE(*,*)
        WRITE(*,*)'***********************'
        WRITE(*,*)'***** MESHRENOEF ******'
        WRITE(*,*)'***********************'
        WRITE(*,*)

        NL=NL_INPUT ! Recover from backup!
        IF (NL.GT.0) THEN
         WRITE(*,*)'==== Building the grid...', ITER, NL
         CALL CREATE_MESH(ITER,NX,NY,NZ,NL,NPATCH,PARE,PATCHNX,PATCHNY,
     &                    PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,
     &                    PATCHRZ,N_PARTICLES,N_DM,N_GAS,LADO0,T,ZETA,
     &                    REFINE_THR,MIN_PATCHSIZE,MINFRAC_REFINABLE,
     &                    BOR,BORAMR,BOR_OVLP,NPART_ESP,FW1)
         WRITE(*,*)'==== END building the grid...', ITER, NL
        END IF

       CALL DENSITY(ITER,NX,NY,NZ,NL,NPATCH,PARE,PATCHNX,PATCHNY,
     &              PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,
     &              PATCHRZ,N_PARTICLES,N_DM,N_GAS,LADO0,T,ZETA,
     &              NPART_ESP,INTERP_DEGREE)

       WRITE(*,*)'***************************'
       WRITE(*,*)'***** END MESHRENOEF ******'
       WRITE(*,*)'***************************'

c       CALL POISSON(NL,NX,NY,NZ,DX,NPATCH,PARE,PATCHNX,PATCHNY,
c     &              PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,
c     &              PATCHRZ,RXPA,RYPA,RZPA,MASAP,N_PARTICLES,N_DM,
c     &              LADO0,POT,POT1)

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

c       WRITE(*,*) '***************************'
       WRITE(*,*) '* COSMOLOGICAL PARAMETERS *'
       WRITE(*,*) '***************************'
       WRITE(*,*) 'RETE=', RETE
       WRITE(*,*) 'ROTE=', ROTE
       WRITE(*,*) 'RODO,RE0,OMEGA0,OMEGAZ=', RODO,RE0,OMEGA0,OMEGAZ
       WRITE(*,*) 'Z=', ZETA
       WRITE(*,*) 'CONTRASTEC=',CONTRASTEC
c       WRITE(*,*) '***************************'

**************************************************************
*      Cleaning overlaps of patches
*      NOTE! we correct overlaps and not refinements because
*            we work within each level independentely
**************************************************************

       ! SOLAP: overlaps at level IR; (=1, keep), (=0, overlapped)
       DO IR=1,NL
        CALL VEINSGRID(IR,NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &                 PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &                 SOLAP)
       END DO

******* Compute CR0AMR *********************************************
       CALL COMPUTE_CR0AMR(NL,NX,NY,NZ,NPATCH,PARE,PATCHNX,PATCHNY,
     &                     PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                     PATCHRY,PATCHRZ,CR0AMR,CR0AMR11,LADO0)


       CALL RENORM_DENSITY(NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                     U1,U11,LADO0,RODO,RE0,FDM)

       IF (FW2.EQ.1) THEN
        OPEN(99,
     &       FILE='output_files/density_asohf'//ITER_STRING//'.res',
     &       STATUS='UNKNOWN',FORM='UNFORMATTED')
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
       END IF
*********************************************************************

**********************************************************************
******************************HALO FINDER*****************************
**********************************************************************

       WRITE(*,*)
       WRITE(*,*) '***************************'
       WRITE(*,*) '***    HALO FINDING     ***'
       WRITE(*,*) '***************************'

       CALL SORT_DM_PARTICLES_X(N_DM,NDMPART_X,NX,LADO0)

       open(55,file='./output_files/haloesgrids'//ITER_STRING//'.res',
     &        status='unknown')
       close(55)

**********************************************************
*      We proceed level by level, from coarse to fine
**********************************************************
       DO IR=0,NL

*       Looking for candidate haloes at each of the AMR levels

        CALL HALOFIND_GRID(IR,NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                     PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                     PATCHRY,PATCHRZ,PARE,NCLUS,MASA,RADIO,
     &                     CLUSRX,CLUSRY,CLUSRZ,REALCLUS,LEVHAL,
     &                     NHALLEV,BOUND,CONTRASTEC,RODO,
     &                     SOLAP,CR0AMR,CR0AMR11,PATCHCLUS,
     &                     VOL_SOLAP_LOW,CLUSRXCM,CLUSRYCM,CLUSRZCM)


*       Debug dumps

        IF (FW3.EQ.1) THEN
         open(55,
     &        file='./output_files/haloesgrids'//ITER_STRING//'.res',
     &        status='unknown',position='append')
         do i=1,nclus
          write(55,*) clusrx(i),clusry(i),clusrz(i),radio(i),masa(i),
     &                levhal(i), realclus(i), patchclus(i)
         end do
         close(55)
        END IF

        IF (NHALLEV(IR).EQ.0) CYCLE

*       Sorting out all the haloes

        CALL RE_SORT_HALOES(NCLUS,NHALLEV,REALCLUS,CLUSRX,CLUSRY,CLUSRZ,
     &                      RADIO,MASA,LEVHAL,PATCHCLUS,DMPCLUS,
     &                      HALBORDERS,CLUSRXCM,CLUSRYCM,CLUSRZCM,IR)

        LOW1=SUM(NHALLEV(0:IR-1))+1
        LOW2=SUM(NHALLEV(0:IR))
        WRITE(*,*)'... Masses: MIN, MAX and MEAN=',
     &             MINVAL(MASA(LOW1:LOW2))*UM,
     &             MAXVAL(MASA(LOW1:LOW2))*UM,
     &             SUM(MASA(LOW1:LOW2))/NHALLEV(IR)*UM
        !WRITE(*,*) 'NCLUS=', NCLUS

*       Haloes at the edge of the box (JUST FOR CAUTION!)

        IF (BORDES.EQ.1) THEN
         CALL HALOES_BORDER(NCLUS,CLUSRX,CLUSRY,CLUSRZ,RADIO,LADO0,
     &                      HALBORDERS,LOW1,LOW2)
         WRITE(*,*) '... Haloes close to the box borders:',
     &              SUM(HALBORDERS(LOW1:LOW2))
        END IF

*       Eliminating POOR haloes (less than a minimum number of particles)
*       (We start here to work with partciles for the 1st time)

        CALL PRUNE_POOR_HALOES(NCLUS,CLUSRX,CLUSRY,CLUSRZ,RADIO,
     &                         REALCLUS,N_DM,MIN_NUM_PART,DMPCLUS,
     &                         NDMPART_X,LADO0,1.0,1,LOW1,LOW2)
        CALL RE_SORT_HALOES(NCLUS,NHALLEV,REALCLUS,CLUSRX,CLUSRY,CLUSRZ,
     &                      RADIO,MASA,LEVHAL,PATCHCLUS,DMPCLUS,
     &                      HALBORDERS,CLUSRXCM,CLUSRYCM,CLUSRZCM,IR)

        LOW1=SUM(NHALLEV(0:IR-1))+1
        LOW2=SUM(NHALLEV(0:IR))

        !WRITE(*,*) '---------------------------------'
        !WRITE(*,*) 'CHECKING GRID FINDING AT LEVEL',IR
        KONTA2=COUNT(REALCLUS(LOW1:LOW2).EQ.-1)
        WRITE(*,*) 'Candidates at this level --->', KONTA2
        !WRITE(*,*) '---------------------------------'

*       REFINING REAL HALOES WITH THE DM PARTICLES ONLY

        WRITE(*,*)
        WRITE(*,*)'Refining with DM particles...'

        CALL HALOFIND_PARTICLES(NL,NCLUS,MASA,RADIO,CLUSRX,CLUSRY,
     &       CLUSRZ,REALCLUS,CONCENTRA,ANGULARM,VMAXCLUS,IPLIP,VX,VY,VZ,
     &       VCMAX,MCMAX,RCMAX,M200C,M500C,M2500C,M200M,M500M,M2500M,
     &       MSUB,R200C,R500C,R2500C,R200M,R500M,R2500M,RSUB,DMPCLUS,
     &       LEVHAL,EIGENVAL,N_DM,CONTRASTEC,OMEGAZ,UM,UV,LADO0,
     &       CLUSRXCM,CLUSRYCM,CLUSRZCM,MEAN_VR,INERTIA_TENSOR,NPATCH,
     &       PATCHCLUS,PROFILES,VELOCITY_DISPERSION,KINETIC_E,
     &       POTENTIAL_E,DO_COMPUTE_ENERGIES,INDCS_PARTICLES_PER_HALO,
     &       FLAG_WDM,ZETA,MIN_NUM_PART,NDMPART_X,VAR,MAX_PART_DSUM,
     &       LOW1,LOW2)

*       General CHECKING
        !WRITE(*,*)'HALOES WITHOUT MASS=',COUNT(MASA(LOW1:LOW2).LE.0.0)
        WRITE(*,*) '... Total number of particles within halos at IR:',
     &              SUM(DMPCLUS(LOW1:LOW2))

        CALL PRUNE_POOR_HALOES(NCLUS,CLUSRX,CLUSRY,CLUSRZ,RADIO,
     &                         REALCLUS,N_DM,MIN_NUM_PART,DMPCLUS,
     &                         NDMPART_X,LADO0,1.0,0,LOW1,LOW2)
        CALL CHECK_RUBISH(NCLUS,REALCLUS,CLUSRX,CLUSRY,CLUSRZ,VX,VY,VZ,
     &                    MASA,RADIO,LEVHAL,LOW1,LOW2)
        CALL ACCIDENTAL_SUBSTRUCTURE(NCLUS,REALCLUS,CLUSRX,CLUSRY,
     &                               CLUSRZ,VX,VY,VZ,MASA,RADIO,LEVHAL)


        WRITE(*,*) 'Finally, haloes found at this level --->',
     &             COUNT(REALCLUS(LOW1:LOW2).EQ.-1)

       END DO !IR=0,NL

       WRITE(*,*)
       WRITE(*,*) 'At the end...'
       WRITE(*,*) 'TOTAL NUMBER OF HALOS=',
     &            COUNT(REALCLUS(1:NCLUS).EQ.-1)
       WRITE(*,*)'=================================='
       DO I=0,NL
       WRITE(*,*)'Haloes at level ', I,' =',
     &            COUNT(LEVHAL(1:NCLUS).EQ.I.
     &            AND.REALCLUS(1:NCLUS).NE.0),
     &            COUNT(REALCLUS(1:NCLUS).EQ.-1)
       END DO
       WRITE(*,*)'=================================='
       SUBS_LEV(0)=NCLUS


****************************************************
****************************************************
****************************************************
       IF (FLAG_SUBS.EQ.1) THEN
        WRITE(*,*)
        WRITE(*,*) '***************************'
        WRITE(*,*) '** SUBSTRUCTURE SEARCH   **'
        WRITE(*,*) '***************************'

        IF (FW4.EQ.1) THEN
         OPEN(99, file='./output_files/substructuregrid'//
     &        ITER_STRING//'.res',status='unknown')
         CLOSE(99)
        END IF
        IF (FW5.EQ.1) THEN
         OPEN(100, file='./output_files/substructureparticles'//
     &        ITER_STRING//'.res',status='unknown')
         CLOSE(100)
        END IF

        DO IR=1,NL
         CALL SEARCH_SUBSTRUCTURE_GRID(IR,NL,NX,NY,NZ,NPATCH,PATCHNX,
     &                    PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                    PATCHRY,PATCHRZ,PARE,NCLUS,MASA,RADIO,CLUSRX,
     &                    CLUSRY,CLUSRZ,REALCLUS,LEVHAL,NHALLEV,BOUND,
     &                    CONTRASTEC,RODO,SOLAP,CR0AMR,CR0AMR11,
     &                    PATCHCLUS,VOL_SOLAP_LOW,CLUSRXCM,CLUSRYCM,
     &                    CLUSRZCM,RSUB,MSUB,SUBS_LEV,UM,PROFILES)

         IF (FW4.EQ.1) THEN
          open(99, file='./output_files/substructuregrid'//
     &        ITER_STRING//'.res',status='unknown',position='append')
          do i=subs_lev(0)+1,nclus
           write(99,*) clusrx(i),clusry(i),clusrz(i),msub(i),rsub(i),
     &               realclus(i)
          end do
          close(99)
         END IF

         CALL SUBSTRUCTURE_PARTICLES(IR,NL,NCLUS,MASA,RADIO,CLUSRX,
     &      CLUSRY,CLUSRZ,REALCLUS,CONCENTRA,ANGULARM,VMAXCLUS,IPLIP,VX,
     &      VY,VZ,VCMAX,MCMAX,RCMAX,M200C,M500C,M2500C,M200M,M500M,
     &      M2500M,MSUB,R200C,R500C,R2500C,R200M,R500M,R2500M,RSUB,
     &      DMPCLUS,LEVHAL,EIGENVAL,N_DM,CONTRASTEC,OMEGAZ,UM,UV,LADO0,
     &      CLUSRXCM,CLUSRYCM,CLUSRZCM,MEAN_VR,INERTIA_TENSOR,SUBS_LEV,
     &      PATCHCLUS,NPATCH,PROFILES,VELOCITY_DISPERSION,KINETIC_E,
     &      POTENTIAL_E,DO_COMPUTE_ENERGIES,INDCS_PARTICLES_PER_HALO,
     &      FLAG_WDM,ZETA,MIN_NUM_PART,MAX_NUM_PART,NDMPART_X,VAR,
     &      MAX_PART_DSUM)

         IF (FW5.EQ.1) THEN
          open(100, file='./output_files/substructureparticles'//
     &        ITER_STRING//'.res', status='unknown',position='append')
          do i=subs_lev(0)+1,nclus
           write(100,*) clusrx(i),clusry(i),clusrz(i),msub(i),rsub(i),
     &              realclus(i)
          end do
          close(100)
         END IF
        END DO

        WRITE(*,*) '===> TOTAL NUMBER OF SUBSTRCTURES:',
     &             COUNT(REALCLUS(1:NCLUS).GT.0),' <==='

        CALL FRACTION_MASS_SUBS(NCLUS,REALCLUS,MASA,MSUB,FSUB,NSUBS)
       END IF

****************************************************
****************************************************
****************************************************
*       IF (FLAG_CENTRAL.EQ.1) THEN
*       WRITE(*,*)
*       WRITE(*,*) '***************************'
*       WRITE(*,*) '**     CORE SEARCH       **'
*       WRITE(*,*) '***************************'
*
*       CALL CORE_SEARCH(NCLUS,MASA,RADIO,CLUSRX,CLUSRY,CLUSRZ,REALCLUS,
*     &                  MSUB,RSUB,SUBS_LEV,DMPCLUS,RMAXSIGMA,RXPA,RYPA,
*     &                  RZPA,MASAP,U2DM,U3DM,U4DM,N_DM,MMAXSIGMA,
*     &                  MAX_NUM_PART)
*       END IF

****************************************************
**********          STELLAR HALOES           *******
****************************************************

       IF (VAR.GT.1.AND.N_ST.GT.0.AND.FLAG_STELLAR.EQ.1) THEN
        WRITE(*,*)
        WRITE(*,*) '***************************'
        WRITE(*,*) '**    STELLAR HALOES     **'
        WRITE(*,*) '***************************'

        CALL STELLAR_HALOES(NCLUS,MASA,RADIO,MSUB,RSUB,REALCLUS,DMPCLUS,
     &                      CLUSRX,CLUSRY,CLUSRZ,N_DM,N_ST,NX,LADO0,
     &                      INDCS_PARTICLES_PER_HALO,UM,UV,
     &                      MIN_NUM_PART_ST,FLAG_WDM,ITER,ZETA,
     &                      STPAR_FACT_INC,STPAR_MAX_DIST,
     &                      STPAR_MIN_OVERDENS,STPAR_MAX_R_PHYS)
       END IF

*************************************************
*************************************************
*=================== OUTPUT FILES ===============
*************************************************
*************************************************

       CALL DECONVER_POSITIONS(NCLUS,REALCLUS,MAXNCLUS,MAXNCLUS,
     &                         CLUSRX,CLUSRY,CLUSRZ)
       CALL DECONVER_POSITIONS(NCLUS,REALCLUS,MAXNCLUS,MAXNCLUS,
     &                         CLUSRXCM,CLUSRYCM,CLUSRZCM)

       KONTA2=COUNT(REALCLUS(1:NCLUS).NE.0)

       OPEN(3,FILE='./output_files/families'//ITER_STRING,
     &      STATUS='UNKNOWN')
       IF (FLAG_WDM.EQ.1) THEN
        OPEN(4,FILE='./output_files/particles'//ITER_STRING,
     &       FORM='UNFORMATTED')
        WRITE(4) KONTA2
       END IF

       WRITE(3,*) '*********************NEW ITER*******************'
       WRITE(3,*) ITER, NCLUS, KONTA2, ZETA
       WRITE(3,*) '************************************************'

111    FORMAT(51A14)
112    FORMAT(2I14,3F14.6,E14.6,F14.6,E14.6,F14.6,2I14,3F14.6,3F14.6,
     &        6E14.6,3E14.6,F14.6,3F14.3,2F14.3,2E14.6,F14.3,E14.3,
     &        F14.6,F14.6,E14.6,F14.6,E14.6,F14.6,E14.6,F14.6,E14.6,
     &        F14.6,E14.6,F14.6,E14.6,F14.6,I14)


       WRITE(3,*) '=====================================================
     &==================================================================
     &==================================================================
     &==================================================================
     &==================================================================
     &================================================================='

       WRITE(3,111) 'Halo ID'     ,'Substr. of'  ,'Density peak',
     &            'coordinates' ,'(Mpc)'       ,'Virial mass' ,
     &            'Virial radi' ,'Substr. mass','Substr. radi',
     &            'Part. num.'  ,'Most bound'  ,'Center of'   ,
     &            'mass coords' ,'(Mpc)'       ,'Semiaxes'    ,
     &            '(Mpc)'       ,''            ,'Inertia'     ,
     &            'tensor'      ,'components'  ,'(cMpc^2)'    ,
     &            ''            ,''            ,'Spec. angul.',
     &            'momentum'    ,'(cMpc km/s)' ,'Veloc. disp.',
     &            'Bulk velocty','(km/s)'      ,''            ,
     &            'Max part vel','Mean V_rad'  ,'Kinetic E'   ,
     &            'Potential E' ,'Vcmax'       ,'Mass@Vcmax'  ,
     &            'r@Vcmax'     ,'R200m (Mpc)' ,'M200m (Msun)',
     &            'R200c (Mpc)' ,'M200c (Msun)','R500m (Mpc)' ,
     &            'M500m (Msun)','R500c (Mpc)' ,'M500c (Msun)',
     &            'R2500m (Mpc)','M2500m(Msun)','R2500c (Mpc)',
     &            'M2500c(Msun)','f_sub'       ,'N_subs'

       WRITE(3,111) ''            ,''            ,'x'           ,
     &            'y'           ,'z'           ,'(Msun)'      ,
     &            '(Mpc)'       ,'(Msun)'      ,'(Mpc)'       ,
     &            ''            ,'particle ID' ,'x'           ,
     &            'y'           ,'z'           ,'Major'       ,
     &            'Intermediate','Minor'       ,'Ixx'         ,
     &            'Ixy'         ,'Ixz'         ,'Iyy'         ,
     &            'Iyz'         ,'Izz'         ,'Lx'          ,
     &            'Ly'          ,'Lz'          ,'(km/s)'      ,
     &            'Vx'          ,'Vy'          ,'Vz'          ,
     &            '(km/s)'      ,'(km/s)'      ,'Msun(km/s)^2',
     &            'Msun(km/s)^2','(km/s)'      ,'(Msun)'      ,
     &            '(Mpc)'       ,''            ,''            ,
     &            ''            ,''            ,''            ,
     &            ''            ,''            ,''            ,
     &            ''            ,''            ,''            ,
     &            ''            ,''            ,''

       WRITE(3,*) '=====================================================
     &==================================================================
     &==================================================================
     &==================================================================
     &==================================================================
     &================================================================='

       KONTA2=0
       DO I=1,NCLUS

       IF (REALCLUS(I).NE.0) THEN
         WRITE(3,112) I,REALCLUS(I),CLUSRX(I),CLUSRY(I),CLUSRZ(I),
     &              MASA(I),RADIO(I),MSUB(I),RSUB(I),DMPCLUS(I),
     &              IPLIP(I),CLUSRXCM(I),CLUSRYCM(I),CLUSRZCM(I),
     &              (EIGENVAL(J,I),J=1,3),(INERTIA_TENSOR(J,I),J=1,6),
     &              (ANGULARM(J,I),J=1,3),VELOCITY_DISPERSION(I),
     &              VX(I)*UV,VY(I)*UV,VZ(I)*UV,VMAXCLUS(I),MEAN_VR(I),
     &              KINETIC_E(I),POTENTIAL_E(I),
     &              VCMAX(I),MCMAX(I),RCMAX(I),
     &              R200M(I),M200M(I),R200C(I),M200C(I),
     &              R500M(I),M500M(I),R500C(I),M500C(I),
     &              R2500M(I),M2500M(I),R2500C(I),M2500C(I),
     &              FSUB(I),NSUBS(I)
         IF (FLAG_WDM.EQ.1) THEN
          WRITE(4) I,(INDCS_PARTICLES_PER_HALO(J,I),J=1,2)
          KONTA2=MAX(KONTA2,INDCS_PARTICLES_PER_HALO(2,I))
         END IF
       END IF  !realclus

       END DO

       IF (FLAG_WDM.EQ.1) THEN
        WRITE(4) KONTA2
        WRITE(4) (PARTICLES_PER_HALO(J),J=1,KONTA2)
       END IF

       CLOSE(3)
       IF (FLAG_WDM.EQ.1) CLOSE(4)

*==========================================
*************************************************
*************************************************

       CALL DEALLOC

       WRITE(*,*)'END ITER', ITER
       CALL IDATE(DATE)
       CALL ITIME(TIME)
       WRITE(*,*) 'DATE=',DATE(1),'/',DATE(2),'/',DATE(3)
       WRITE(*,*) 'TIME=',TIME(1),':',TIME(2),':',TIME(3)

       END DO ! IFI2
**//////////////////////////////////////////////////////////////
       END
**//////////////////////////////////////////////////////////////

***********************************************************************
*      SUBROUTINES IN EXTERNAL FILES                                  *
***********************************************************************
*      Grid building
       INCLUDE 'grids.f'
*      Routines with basic numerical calculations stuff
       INCLUDE 'num.f'
*      Halo finding procedures using the grid
       INCLUDE 'haloes_grids.f'
*      Halo finding procedures using particles
       INCLUDE 'haloes_particles.f'
*      Substructure search
       INCLUDE 'substructure.f'
*      Read MASCLET outputs (can be changed for other code outputs)
       INCLUDE 'reader.f'
*      Calculate properties of stellar haloes
       INCLUDE 'stars.f'
*      Allocating and initialising stuff
       INCLUDE 'alloc.f'
*      Solve Poisson's equation for the gravitational potential
*       generated by DM
c       INCLUDE 'poisson.f'
