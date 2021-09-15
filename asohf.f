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
             
       INTEGER CONTA2(NAMRX,NAMRY,NAMRZ,NPALEV)
       
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
       INTEGER N_GAS,N_DM
       REAL*4 COTA(NCOTAS,0:NLEVELS)
       
*      ---UNBINDING---       
       REAL*4 REF_MIN, REF_MAX
       INTEGER SALIDA,FAC
       REAL*4 DISTA(0:PARTI), DISTA2(0:PARTI)
       INTEGER QUIEN(PARTI), QUIEN2(PARTI)
       INTEGER INDICE(PARTI)

       REAL*4 DENSITOT(0:PARTI/10),RADIAL(0:PARTI/10)
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
       READ(1,*) NDXYZ
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
       END IF        


       END IF  !PLOT>1


       CONTAITER=0
       MARK(1:NFILE2)=0
       
       INERTIA=0.0


*///////////////////////////////////////////      
       DO IFI2=1, NFILE2
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

       
cc!$OMP PARALLEL DO SHARED(NHALLEV,NPART,NLEVBAS),
cc!$OMP+        PRIVATE(IR)       
c       DO IR=0,NLEVBAS
        NHALLEV=0
        NPART=0
c       END DO


      
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
*     READING INITIAL DATA
***************************************************

       IF (FLAG_SA.EQ.1) THEN

*      Reading MASCLETs files directly
       CALL LEER(VAR,ITER,NX,NY,NZ,NDXYZ,T,ZETA,NL,NPATCH,
     &           PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &           PATCHRX,PATCHRY,PATCHRZ,MAP,
     &           U2DM,U3DM,U4DM,MASAP,NPART,RXPA,RYPA,RZPA)

       
       ROTE=RODO*(1.0+ZETA)**3   
       RETE=RE0/(1.0+ZETA)
       
       ZETAS(IFI)=ZETA
       TIEMPO(IFI)=T 
       
       ELSE
       
*      Reading external list of particles
       COTA=5.0           !OJO! hay que poner a mano el valor de las COTAS!
       WRITE(*,*)'COTA=', COTA
       
       WRITE(*,*)'***********************'
       WRITE(*,*)'***** MESHRENOEF ******'
       WRITE(*,*)'***********************'
       
       WRITE(*,*)'Building the grid...', ITER, IFI2, IFI
       CALL MESHRENOEF(ITER,NX,NY,NZ,NL,COTA,NPATCH,
     &             NPART,PATCHNX,PATCHNY,PATCHNZ,
     &             PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,
     &             PATCHRZ,PARE,U2DM,U3DM,U4DM,MASAP,MAP,
     &             RXPA,RYPA,RZPA,ZETA,T,LADO0,FLAG_MASCLET,PLOT)
       WRITE(*,*)'END building the grid...', ITER, IFI2, IFI

       ROTE=RODO*(1.0+ZETA)**3   
       RETE=RE0/(1.0+ZETA)
       
       ZETAS(IFI)=ZETA
       TIEMPO(IFI)=T 
       
!$OMP PARALLEL DO SHARED(NX,NY,NZ,U1,RETE,ROTE),PRIVATE(I,J,K)
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
       U1(I,J,K)=U1(I,J,K)/(RETE*RETE*RETE)
       U1(I,J,K)=(U1(I,J,K)/ROTE)-1.0
      END DO
      END DO
      END DO
       
       DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR)) 
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,
!$OMP+                   U11,RETE,ROTE,LOW1,LOW2),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I)
        DO I=LOW1,LOW2
         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)
         DO KZ=1,N3
         DO JY=1,N2
         DO IX=1,N1
          U11(IX,JY,KZ,I)=U11(IX,JY,KZ,I)/(RETE*RETE*RETE)
          U11(IX,JY,KZ,I)=(U11(IX,JY,KZ,I)/ROTE)-1.0                               
         END DO
         END DO
         END DO
        END DO
       END DO


       IF (FLAG_GAS.EQ.1) THEN 
       
!$OMP PARALLEL DO SHARED(NX,NY,NZ,U1G,RETE,ROTE),PRIVATE(I,J,K)
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
       U1G(I,J,K)=U1G(I,J,K)/(RETE*RETE*RETE)
       U1G(I,J,K)=(U1G(I,J,K)/ROTE)-1.0
      END DO
      END DO
      END DO

       DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR)) 
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,
!$OMP+                   U11G,RETE,ROTE,LOW1,LOW2),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I)
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        DO KZ=1,N3
        DO JY=1,N2
        DO IX=1,N1
         U11G(IX,JY,KZ,I)=U11G(IX,JY,KZ,I)/(RETE*RETE*RETE)
         U11G(IX,JY,KZ,I)=(U11G(IX,JY,KZ,I)/ROTE)-1.0                             
        END DO
        END DO
        END DO
       END DO
       END DO
       
       END IF 
       
       WRITE(*,*)'***************************'
       WRITE(*,*)'***** END MESHRENOEF ******'
       WRITE(*,*)'***************************'
       
       END IF 

****************************************************************
*      VIRIAL CONTRAST
*      (Bryan & Norman ApJ, 1998)
****************************************************************
       CONTRASTEX=0.0
       CONTRASTEC=0.0
       OMEGAZ=0.0
       CONTRASTEX=OMEGA0*(1.0+ZETA)**3
       CONTRASTEX=CONTRASTEX/(CONTRASTEX+1.0-OMEGA0)
       OMEGAZ=CONTRASTEX
       CONTRASTEX=CONTRASTEX-1.0
       CONTRASTEC=18.0*PI*PI+82.0*CONTRASTEX-39.0*CONTRASTEX*CONTRASTEX
       CONTRASTEC=CONTRASTEC/OMEGAZ       !!!este es el correcto
*****       CONTRASTEC=CONTRASTEC/OMEGA0    !!!!prueba!!!
*
*       CONTRASTEC=CONTRASTEC*1000.0
*
       WRITE(*,*) '************************************************'
       WRITE(*,*) '             "COSMOLOGICAL" PARAMETERS            '
       WRITE(*,*) '************************************************'
       WRITE(*,*) 'RETE=', RETE
       WRITE(*,*) 'ROTE=', ROTE
       WRITE(*,*) 'RODO,RE0,OMEGA0,OMEGAZ=', RODO,RE0,OMEGA0,OMEGAZ
       WRITE(*,*) 'Z=', ZETA
       WRITE(*,*) 'CONTRASTEC_OK=', CONTRASTEC*OMEGA0/OMEGAZ
       WRITE(*,*) 'CONTRASTEC=',CONTRASTEC,CONTRASTEC*OMEGA0,
     &                          CONTRASTEC*ROTE,CONTRASTEC/RODO
       WRITE(*,*)'CONTRASTEC_200=',200.0*ROTE*OMEGA0,200.0/OMEGAZ
       WRITE(*,*) '************************************************' 
       

****************************************************************
*      Joining all the densities in U1 and U11
****************************************************************
       WRITE(*,*) 'U1_000', MAXVAL(U1(1:NX,1:NY,1:NZ)), 
     &                      MINVAL(U1(1:NX,1:NY,1:NZ)) 
       WRITE(*,*) 'U1G_000',MAXVAL(U1G(1:NX,1:NY,1:NZ)), 
     &                      MINVAL(U1G(1:NX,1:NY,1:NZ)) 

*      --Only GAS--
       IF (VAR.EQ.1.AND.FLAG_GAS.EQ.1) THEN

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U1,U1G),PRIVATE(I,J,K)
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
       U1(I,J,K)=U1G(I,J,K)+1.0
      END DO
      END DO
      END DO

       DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR)) 
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,
!$OMP+                   U11,U11G,LOW1,LOW2),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I)
       DO I=LOW1,LOW1
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        DO KZ=1,N3
        DO JY=1,N2
        DO IX=1,N1
         U11(IX,JY,KZ,I)=U11G(IX,JY,KZ,I)+1.0
        END DO
        END DO
        END DO
       END DO
       END DO 


       END IF   !ONLY_GAS

*      --Only DM--       
       IF (VAR.EQ.2) THEN

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U1),PRIVATE(I,J,K)
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
       U1(I,J,K)=U1(I,J,K)+1.0
      END DO
      END DO
      END DO

      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR)) 
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,
!$OMP+                   U11,LOW1,LOW2),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I)
       DO I=LOW1, LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        DO KZ=1,N3
        DO JY=1,N2
        DO IX=1,N1
         U11(IX,JY,KZ,I)=U11(IX,JY,KZ,I)+1.0
        END DO
        END DO
        END DO
       END DO
       END DO 


       END IF    !ONLY_DM

*      --GAS+DM--       
       IF (VAR.EQ.0) THEN

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U1,U1G),PRIVATE(I,J,K)
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
       U1(I,J,K)=U1(I,J,K)+U1G(I,J,K)+1.0
      END DO
      END DO
      END DO


      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR)) 
!$OMP PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,
!$OMP+                   U11,U11G,LOW1,LOW2),
!$OMP+            PRIVATE(IX,JY,KZ,N1,N2,N3,I)
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        DO KZ=1,N3
        DO JY=1,N2
        DO IX=1,N1
         U11(IX,JY,KZ,I)=U11(IX,JY,KZ,I)+U11G(IX,JY,KZ,I)+1.0
        END DO
        END DO
        END DO
       END DO
       END DO 

       END IF  !GAS_&_DM
  
       
       WRITE(*,*) 'U1_0',MAXVAL(U1(1:NX,1:NY,1:NZ)), 
     &                   MINVAL(U1(1:NX,1:NY,1:NZ))
*************************************************************


**************************************************************
*      Cleaning overlaps of patches
*      NOTE! we correct overlaps and not refinements because
*            we work within each level independentely 
**************************************************************
       
       CONTA2=1   !AHORA HACE DE SOLAP (informacion solapada)
       DO IR=1,NL
       CALL VEINSGRID(IR,NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &           PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &           CONTA2,VECINO,NVECI)
       END DO
       
       
       DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))

!$OMP  PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,
!$OMP+                   U11,CONTA2,LOW1,LOW2),
!$OMP+            PRIVATE(N1,N2,N3,I,IX,JY,KZ)
       DO I=LOW1,LOW2
       N1=0
       N2=0
       N3=0
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)
       
        DO KZ=1,N3
        DO JY=1,N2
        DO IX=1,N1
        !Limpiamos U11 de solapes: 
       U11(IX,JY,KZ,I)=U11(IX,JY,KZ,I)*CONTA2(IX,JY,KZ,I)
       
       END DO
       END DO
       END DO
       
       END DO
       END DO
****************************************************************


**********************************************************************
******************************HALO FINDER*****************************
**********************************************************************

       WRITE(*,*) '************************************************'
       WRITE(*,*) '         HALO FINDING                           '
       WRITE(*,*) '************************************************'


**********************************************************
*      Looking for candidate haloes at the AMR levels
**********************************************************

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



**************************************************************
*      NIVEL BASE!!
**************************************************************

       IR=0
       REF=0.0
       REF=0.2*DX
       ESP=(LOG10(LADO*0.5)-LOG10(REF))/(NSHELL-1) 
       WRITE(*,*) 'RESOLUTION (IR,DX,NPART(0))=',IR,DX,NPART(0)
       
C_VICENT:PARALELIZAR       
       CONTA(1:NX,1:NY,1:NZ)=0        !celdas
       ICEN=0
       
       UBAS1(1:NX,1:NY,1:NZ)=0.0
       UBAS1(1:NX,1:NY,1:NZ)=U1(1:NX,1:NY,1:NZ)
       
       WRITE(*,*)'ESTIMATION_0:',IR,
     &            COUNT(UBAS1(1:NX,1:NY,1:NZ).GE.CONTRASTEC)
       
*      Ordenamos todas las celdas de los vecinos
       KK_ENTERO=0
       KK_ENTERO=COUNT(UBAS1(1:NX,1:NY,1:NZ).GE.CONTRASTEC)
       
       ALLOCATE(DDD(KK_ENTERO))
       ALLOCATE(DDDX(KK_ENTERO))
       ALLOCATE(DDDY(KK_ENTERO))
       ALLOCATE(DDDZ(KK_ENTERO)) 
       ALLOCATE(INDICE2(KK_ENTERO)) 
       ALLOCATE(DDD2(KK_ENTERO))
       ALLOCATE(DDDX2(KK_ENTERO))
       ALLOCATE(DDDY2(KK_ENTERO))
       ALLOCATE(DDDZ2(KK_ENTERO))
       
       DDD(1:KK_ENTERO)=0.0
       DDDX(1:KK_ENTERO)=0
       DDDY(1:KK_ENTERO)=0
       DDDZ(1:KK_ENTERO)=0
       INDICE2(1:KK_ENTERO)=0
       DDD2(1:KK_ENTERO)=0.0
       DDDX2(1:KK_ENTERO)=0
       DDDY2(1:KK_ENTERO)=0
       DDDZ2(1:KK_ENTERO)=0
       
       II=0
        
        DO KZ=1, NZ
        DO JY=1, NY
        DO IX=1, NX
       
        KK_REAL=0.0
        KK_REAL=UBAS1(IX,JY,KZ)
        
        IF(KK_REAL.GE.CONTRASTEC) THEN
         II=II+1
         DDD(II)=1.0/KK_REAL
         DDDX2(II)=IX
         DDDY2(II)=JY
         DDDZ2(II)=KZ
        END IF
        
        END DO
        END DO
        END DO
       
       IF (II.NE.KK_ENTERO) THEN
         WRITE(*,*) 'WARNING,NV', II, KK_ENTERO
         STOP
       END IF
       
       CALL INDEXX(KK_ENTERO,DDD(1:KK_ENTERO),INDICE2(1:KK_ENTERO)) 
       DO I=1,KK_ENTERO 
        DDD2(I)=1.0/DDD(INDICE2(I))
       END DO 
       !ahora ya estan ordenadas todas las celdas de IR=0

       DDD=0.0
       DDDX=0
       DDDY=0
       DDDZ=0
       
       DO I=1, KK_ENTERO 
        DDD(I)=DDD2(I)
        DDDX(I)=DDDX2(INDICE2(I))
        DDDY(I)=DDDY2(INDICE2(I))
        DDDZ(I)=DDDZ2(INDICE2(I))
       END DO

       NV_GOOD=0
       NV_GOOD=KK_ENTERO
       KK_ENTERO=0

       DO L1=1,NV_GOOD 

*      de todos los vecinos este es el parche elegido!! 

       ICEN(1)=DDDX(L1)
       ICEN(2)=DDDY(L1)
       ICEN(3)=DDDZ(L1)
       
       KK_ENTERO=0
       KK_ENTERO=CONTA(ICEN(1),ICEN(2),ICEN(3))
       
       IF(KK_ENTERO.EQ.0) THEN
                     
       NCLUS=NCLUS+1
       REALCLUS(IFI,NCLUS)=-1
       LEVHAL(NCLUS)=IR
       NHALLEV(IR)=NHALLEV(IR)+1
       
       IF(NCLUS.GT.MAXNCLUS) THEN
        WRITE(*,*) 'WARNING: NCLUS>MAXNCLUS!!!',NCLUS
        STOP       
       END IF
       
       CLUSRX(NCLUS)=RADX(ICEN(1))
       CLUSRY(NCLUS)=RADY(ICEN(2))
       CLUSRZ(NCLUS)=RADZ(ICEN(3))

*      posible "ALCANCE" del halo ---> mini box around it 
       NX1=0
       NX2=0
       NY1=0
       NY2=0
       NZ1=0
       NZ2=0
       
       PRUEBAX=0.0
       PRUEBAX=CLUSRX(NCLUS)+BOUND
       NX1=INT(((PRUEBAX-RADX(1))/DX)+0.5)+1
       IF (NX1.GT.NX) NX1=NX
       
       PRUEBAX=0.0
       PRUEBAX=CLUSRX(NCLUS)-BOUND
       NX2=INT(((PRUEBAX-RADX(1))/DX)+0.5)+1
       IF (NX2.LT.1) NX2=1
       
       PRUEBAX=0.0
       PRUEBAX=CLUSRY(NCLUS)+BOUND
       NY1=INT(((PRUEBAX-RADY(1))/DY)+0.5)+1
       IF (NY1.GT.NY) NY1=NY
       
       PRUEBAX=0.0
       PRUEBAX=CLUSRY(NCLUS)-BOUND
       NY2=INT(((PRUEBAX-RADY(1))/DY)+0.5)+1
       IF (NY2.LT.1) NY2=1
       
       PRUEBAX=0.0
       PRUEBAX=CLUSRZ(NCLUS)+BOUND
       NZ1=INT(((PRUEBAX-RADZ(1))/DZ)+0.5)+1
       IF (NZ1.GT.NZ) NZ1=NZ
       
       PRUEBAX=0.0
       PRUEBAX=CLUSRZ(NCLUS)-BOUND
       NZ2=INT(((PRUEBAX-RADZ(1))/DZ)+0.5)+1
       IF (NZ2.LT.1) NZ2=1
*
       
       PRUEBAX=0.0
       PRUEBAY=0.0
       PRUEBAZ=0.0
       DELTA_PRUEBA=0.0 

     
*      EXTENDEMOS RADIALMENTE EL CLUSTER
       
       RMIN=0.0       
       DELTA=0.0              
       SHELL=0
       BASMAS=0.0   !MASA TOTAL DEL CUMULO
       DELTA2=0.0   !DELTA TOTAL DEL CUMULO
       VOL2=0.0     !VOL TOTAL DEL CUMULO
       RANT=0.0              
       RSHELL=0.0

       DELTA2=10.0*CONTRASTEC*ROTE

       DO WHILE(DELTA2.GT.(CONTRASTEC)*ROTE)
                  
       IF(SHELL<NSHELL) THEN                 
          SHELL=SHELL+1
          DELTA=0.0
          RANT=RSHELL           
          RSHELL=0.0
          VOL=0.0          
       ELSE
          WRITE(*,*) 'WARNING: DEMASIADAS CAPAS'
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
       
           
*      Nivel base
       
***aki esto se puede acortar
 
       DO K=NZ2, NZ1
       DO J=NY2, NY1
       DO I=NX2, NX1
       
         AA=0.0
         AA=SQRT((RADX(I)-CLUSRX(NCLUS))**2+
     &           (RADY(J)-CLUSRY(NCLUS))**2+
     &           (RADZ(K)-CLUSRZ(NCLUS))**2)
         
         IF(AA.GE.R222.AND.AA.LE.R111) THEN
     
         DELTA=DELTA+U1(I,J,K)*DX*DY*DZ
         CONTA(I,J,K)=1 
         
         PRUEBAX=PRUEBAX+RADX(I)*U1(I,J,K)
         PRUEBAY=PRUEBAY+RADY(J)*U1(I,J,K)
         PRUEBAZ=PRUEBAZ+RADZ(K)*U1(I,J,K)
         DELTA_PRUEBA=DELTA_PRUEBA+U1(I,J,K)
         
         END IF
                         
       END DO
       END DO
       END DO 
       
       BASMAS=BASMAS+DELTA*RODO*RE0**3
       VOL2=PI*(4.0/3.0)*(RSHELL*RETE)**3
       DELTA2=BASMAS/VOL2              
       DELTA=(DELTA*RODO*RE0**3)/VOL  
       RMIN=RSHELL
                    
       END DO   ! while de crecer  (DELTA 2)
              
       RADIO(NCLUS)=RSHELL
       MASA(NCLUS)=BASMAS        
       
       PRUEBAX=PRUEBAX/DELTA_PRUEBA
       PRUEBAY=PRUEBAY/DELTA_PRUEBA
       PRUEBAZ=PRUEBAZ/DELTA_PRUEBA
         
       CLUSRX(NCLUS)=PRUEBAX
       CLUSRY(NCLUS)=PRUEBAY
       CLUSRZ(NCLUS)=PRUEBAZ
       
c8       ELSE
       
c8       UBAS1(ICEN(1),ICEN(2),ICEN(3))=0.0
       
       END IF       
       
       
       END DO     ! while de hallar halos
       

       DEALLOCATE(DDD)
       DEALLOCATE(DDDX)
       DEALLOCATE(DDDY)
       DEALLOCATE(DDDZ) 
       DEALLOCATE(INDICE2) 
       DEALLOCATE(DDD2)
       DEALLOCATE(DDDX2)
       DEALLOCATE(DDDY2)
       DEALLOCATE(DDDZ2)
C888888888888

       
****************************************************
*      CORRECION DE SOLAPES:
*      VAMOS A VER QUE CUMULOS SOLAPAN EN IR
**************************************************** 

       CALL OVERLAPING(IFI,IR,NL,REF,ESP,BOUND,CONTA,CONTRASTEC,
     &                 NSHELL,RODO,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                 PATCHRX,PATCHRY,PATCHRZ,NX,NY,NZ,
     &                 NCLUS,MASA,RADIO,CLUSRX,CLUSRY,CLUSRZ,
     &                 REALCLUS,NSOLAP,SOLAPA,NHALLEV)

**************************************************** 
*      FIN CORRECION DE SOLAPES EN IR     
****************************************************        

       WRITE(*,*)' FIN NIVEL BASE', 0
       
       

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


***************************************************************** 
       SUBROUTINE MERGER(NFILE,NMASA,NNCLUS,DMPITER, 
     &                   DMPCLUS,REALCLUS)
***************************************************************** 
      
       USE MTREE 
       IMPLICIT NONE
       
       INCLUDE 'input_files/asohf_parameters.dat'
       
       INTEGER IFI, NFILE, IFI2,MARCA,MARCA2
       INTEGER I,J,K,I2,J2,IX2,IX1,J1,IX3
       INTEGER KK, KK2, KKK, KKK2
       REAL*4 MASAP(PARTIRED)
       
       INTEGER NNCLUS(MAXITER)
       REAL*4 NMASA(MAXITER,NMAXNCLUS)       
       INTEGER REALCLUS(MAXITER,MAXNCLUS)
       
       INTEGER DMPITER(MAXITER)
       INTEGER DMPCLUS(MAXITER,NMAXNCLUS)
       INTEGER FLAG(NMAXNCLUS)
        
*       ----- module mtree.f -----       
*       INTEGER DMLIP(MAXITER,PARTIRED_PLOT)
*       INTEGER DMLIR(MAXITER,PARTIRED_PLOT)
*       INTEGER NHOST(MAXITER,0:ININL,PARTIRED)
*       INTEGER HHOST(MAXITER,NMAXDAD,0:ININL,PARTIRED)
*       REAL*4 NEW_MASAP(MAXITER,PARTIRED_PLOT)      
*       INTEGER NDAD(MAXITER,NMAXNCLUS)
*       REAL*4 RATIO(MAXITER,NMAXDAD,NMAXCLUS_PLOT)
*       INTEGER DAD(MAXITER,NMAXDAD,NMAXCLUS_PLOT)

              
       IX2=0
       IX1=0
       IX3=0
       IFI2=0
       
       IFI=1
       IFI2=IFI+1
       
       WRITE(*,*) 'ITERS=', IFI, IFI2
       WRITE(*,*) 'NNCLUS(ITER1)=',NNCLUS(IFI)
       WRITE(*,*) 'NNCLUS(ITER2)=',NNCLUS(IFI2)
       WRITE(*,*) 'DMPITER(IFI)=', DMPITER(IFI)
       WRITE(*,*) 'DMPITER(IFI2)=', DMPITER(IFI2)
       WRITE(*,*) '------------------------------------'
       
       DO I2=1, NNCLUS(IFI2)
       
       FLAG(1:NNCLUS(IFI))=0
       
       J1=0
       J1=REALCLUS(IFI2,I2)
       IF(J1.EQ.-1) THEN   !solo de halos reales en IFI2
       
*      particulas del cumulo I2               
       KK2=0
       KKK2=0
       KK2=SUM(DMPCLUS(IFI2,1:I2-1))+1
       KKK2=SUM(DMPCLUS(IFI2,1:I2))
                           
       WRITE(*,*)'I2,IFI2,KK2,KKK2=',I2,IFI2,KK2,KKK2

       DO J2=KK2, KKK2
       
       IX1=0
       IX2=0
       IX1=DMLIP(IFI2,J2)
       IX2=DMLIR(IFI2,J2)

       IX3=NHOST(IFI,IX2,IX1)    !cuantos host tiene esta parti en iter anterior

       DO K=1,IX3

       MARCA=0
       MARCA=HHOST(IFI,K,IX2,IX1)

       IF (MARCA.GT.0) THEN    !algun halo de la iter anterior

*      esto lo anyado para q solo lo haga entre (sub-)halos no eliminados!      
       MARCA2=0
       MARCA2=REALCLUS(IFI,MARCA)
       
       IF (MARCA2.NE.0) THEN

       IF (FLAG(MARCA).EQ.0) THEN
       
        NDAD(IFI2,I2)=NDAD(IFI2,I2)+1
        DAD(IFI2,NDAD(IFI2,I2),I2)=MARCA
        RATIO(IFI2,NDAD(IFI2,I2),I2)=RATIO(IFI2,NDAD(IFI2,I2),I2)
     &                               + NEW_MASAP(IFI2,J2)
        FLAG(MARCA)=NDAD(IFI2,I2)
       
       ELSE

       RATIO(IFI2,FLAG(MARCA),I2)=RATIO(IFI2,FLAG(MARCA),I2)
     &                           + NEW_MASAP(IFI2,J2)
       END IF   !FLAG(MARCA).EQ.0
       
       END IF   !halo o subhalo
       
       END IF   !(MARCA.GT.0) 

       END DO   !!K
              
       END DO   !!J2
       
       IF (NDAD(IFI2,I2).GT.NMAXDAD) THEN
        WRITE(*,*)'WARNING: NDAD!!!'
        STOP
       END IF 
       
       RATIO(IFI2,1:NDAD(IFI2,I2),I2)=RATIO(IFI2,1:NDAD(IFI2,I2),I2)
     &                               *9.1717E18*100.0/NMASA(IFI2,I2)
       
       END IF     !REALCLUS

       END DO
             
       RETURN                                                            
       END
**///////////////////////////////////////////       

****************************************************************** 
       SUBROUTINE REDUCED_MERGER(NNCLUS,REALCLUS,IPLIR,IPLIP,
     &                           NCLUSRX,NCLUSRY,NCLUSRZ)
***************************************************************** 
       
       USE MTREE
       IMPLICIT NONE
       
       INCLUDE 'input_files/asohf_parameters.dat'
       
       INTEGER IFI,IFI2,MARCA,MARCA2
       INTEGER I,J,K,I2,J2,IX2,IX1,J1,IX3,IX4
       REAL*4 DIS,DIS0,XX,YY,ZZ
       
       INTEGER NNCLUS(MAXITER)
       INTEGER REALCLUS(MAXITER,MAXNCLUS)
       
*       ----- module mtree.f -----  
*       INTEGER NHOST(MAXITER,0:ININL,PARTIRED)
*       INTEGER HHOST(MAXITER,NMAXDAD,0:ININL,PARTIRED)
*       INTEGER NDAD(MAXITER,NMAXNCLUS)
*       INTEGER DAD(MAXITER,NMAXDAD,NMAXCLUS_PLOT)
       
       INTEGER IPLIP(MAXITER,NMAXNCLUS)       
       INTEGER IPLIR(MAXITER,NMAXNCLUS)
       
       REAL*4 NCLUSRX(MAXITER,NMAXNCLUS)
       REAL*4 NCLUSRY(MAXITER,NMAXNCLUS)
       REAL*4 NCLUSRZ(MAXITER,NMAXNCLUS)
       
       
       IFI=1
       IFI2=IFI+1
       
       WRITE(*,*) 'ITERS=', IFI, IFI2
       WRITE(*,*) 'NNCLUS(ITER1)=',NNCLUS(IFI)
       WRITE(*,*) 'NNCLUS(ITER2)=',NNCLUS(IFI2)
       WRITE(*,*) '------------------------------------'
        
       DO I2=1, NNCLUS(IFI2)
       
       IF(REALCLUS(IFI2,I2).EQ.-1) THEN
       
       IX1=0
       IX2=0
       IX1=IPLIP(IFI2,I2)         !particula mas ligada de I2 en IFI2
       IX2=IPLIR(IFI2,I2)
       IX3=NHOST(IFI,IX2,IX1)     !num. hosts de esta parti. en IFI
       J1=0

       MARCA=0
     
       IF (IX3.GT.0) THEN

       IF (IX3.GT.1) THEN
         DIS0=100000000.0
         XX=0.0
         YY=0.0
         ZZ=0.0
         XX=NCLUSRX(IFI2,I2)
         YY=NCLUSRY(IFI2,I2)
         ZZ=NCLUSRZ(IFI2,I2)
         DO K=1, IX3
          IX4=HHOST(IFI,K,IX2,IX1)
          DIS=(NCLUSRX(IFI,IX4)-XX)**2+(NCLUSRY(IFI,IX4)-YY)**2+
     &        (NCLUSRZ(IFI,IX4)-ZZ)**2       
          DIS=SQRT(DIS)
          IF (DIS.LT.DIS0) THEN
            DIS0=DIS
            MARCA=IX4      
          END IF
         END DO
       ELSE
         MARCA=HHOST(IFI,IX3,IX2,IX1)
       END IF    

       END IF  !IX3.GT.0    

       IF (MARCA.GT.0) THEN
       
*      esto lo anyado para q solo lo haga entre (sub-)halos no eliminados!      
       MARCA2=0
       MARCA2=REALCLUS(IFI,MARCA)
       
       IF (MARCA2.NE.0) THEN

       NDAD(IFI2,I2)=1  !en este MT es siempre =1
       J1=1       
       DAD(IFI2,J1,I2)=MARCA

       END IF  !MARCA2

       END IF  !MARCA
       
       END IF     !REALCLUS
 
       END DO  !NNCLUS I2
             
       RETURN

       END
**///////////////////////////////////////////       
       

********************************************************************* 
       SUBROUTINE LEER(VAR,ITER,NX,NY,NZ,NDXYZ,T,ZETA,NL,NPATCH,
     &           PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &           PATCHRX,PATCHRY,PATCHRZ,MAP,
     &           U2DM,U3DM, U4DM,MASAP,NPART,RXPA,RYPA,RZPA)
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
                              
************************************************************************
        SUBROUTINE VEINSGRID(IR,NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &             PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,CONTA2,
     &             VECINO, NVECI)
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
       !INTEGER, ALLOCATABLE::VECINO(:,:)
       !INTEGER, ALLOCATABLE::NVECI(:)
       INTEGER NV,A2,B2,C2,K,LOW1, LOW2
       INTEGER VECINO(NPALEV,NPALEV),NVECI(NPALEV)

       INTEGER CONTA2(NAMRX,NAMRY,NAMRZ,NPALEV)
       INTEGER MARCA(NAMRX,NAMRY,NAMRZ,NPALEV)
      
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

       INTEGER IG1,IG2,JG1,JG2,KG1,KG2,IG3,JG3,KG3
       REAL*4 RXFIX,RYFIX,RZFIX
       INTEGER NPALEV2                  
       
       !INTEGER VECINO(NPALEV,NPALEV),NVECI(NPALEV)

       NPALEV2=MAX(100,INT(NPALEV/10))
       !ALLOCATE(VECINO(NPALEV2,NPATCH(IR)))
       !ALLOCATE(NVECI(NPATCH(IR)))
       

       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)

*      built auxiliar grid for comparison
       RXFIX=RADX(1) - DX*0.5 + 0.5*DXPA  
       RYFIX=RADY(1) - DY*0.5 + 0.5*DYPA 
       RZFIX=RADZ(1) - DZ*0.5 + 0.5*DZPA 

       !VECINO=0
       !NVECI=0
*!$OMP PARALLEL DO SHARED(IR,NPATCH,PARE,PATCHX,PATCHY,PATCHZ,
*!$OMP+    PATCHNX,PATCHNY,PATCHNZ,VECINO,NVECI),
*!$OMP+  PRIVATE(I,L1,L2,L3,N1,N2,N3,CR1,CR2,CR3,NV,J,LL1,
*!$OMP+         LL2,LL3,NN1,NN2,NN3,CR4,CR5,CR6,A1,A2,B1,B2,
*!$OMP+         C1,C2)       
       !DO I=1,NPATCH(IR)
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR)) 
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

         CONTA2(:,:,:,I)=1
          
         NV=0
         
         RX1=PATCHRX(I)-0.5*DXPA 
         RY1=PATCHRY(I)-0.5*DYPA 
         RZ1=PATCHRZ(I)-0.5*DZPA 
         RX2=PATCHRX(I)-0.5*DXPA+(N1-1)*DXPA 
         RY2=PATCHRY(I)-0.5*DYPA+(N2-1)*DYPA 
         RZ2=PATCHRZ(I)-0.5*DZPA+(N3-1)*DZPA 
         
!         DO J=1,NPATCH(IR)
!        DO J=SUM(NPATCH(0:IR-1))+1,SUM(NPATCH(0:IR))
         DO J=LOW1,LOW2
         IF (J.NE.I) THEN
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

          CONTROL=0

*         hay que mirar que alguno de los 8 vertices este dentro
*         vertice LL1,LL2,LL3
          RIV1=RXX1
          RIV2=RYY1
          RIV3=RZZ1

          IG1=INT(((RX1-RXFIX)/DXPA)+0.5) + 1
          JG1=INT(((RY1-RYFIX)/DYPA)+0.5) + 1
          KG1=INT(((RZ1-RZFIX)/DZPA)+0.5) + 1

          IG2=IG1 + N1 - 1
          JG2=JG1 + N2 - 1
          KG2=KG1 + N3 - 1

          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1

          IF (IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
           END IF

*         vertice LL1,LL2,CR6
          RIV1=RXX1
          RIV2=RYY1
          RIV3=RZZ2
          
          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1

          IF (CONTROL.EQ.0.AND.IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
          END IF

*         vertice LL1,CR5,LL3
          RIV1=RXX1
          RIV2=RYY2
          RIV3=RZZ1

          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1

          IF (CONTROL.EQ.0.AND.IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
          END IF          

*         vertice LL1,CR5,CR6
          RIV1=RXX1
          RIV2=RYY2
          RIV3=RZZ2

          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1
                                        
          IF (CONTROL.EQ.0.AND.IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
          END IF

*         vertice CR4,LL2,LL3
          RIV1=RXX2
          RIV2=RYY1
          RIV3=RZZ1

          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1

          IF (CONTROL.EQ.0.AND.IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
           END IF          

*         vertice CR4,LL2,CR6
          RIV1=RXX2
          RIV2=RYY1
          RIV3=RZZ2

          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1

          IF (CONTROL.EQ.0.AND.IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
           END IF          

*         vertice CR4,CR5,LL3
          RIV1=RXX2
          RIV2=RYY2
          RIV3=RZZ1

          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1

          IF (CONTROL.EQ.0.AND.IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
          END IF

*         vertice CR4,CR5,CR6
          RIV1=RXX2
          RIV2=RYY2
          RIV3=RZZ2

          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1

          IF (CONTROL.EQ.0.AND.IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
          END IF          


*        caso ----ll1 ---l1 -----cr1 -------cr4
          RIV1=RX1
          RIV2=RY1
          RIV3=RZ1


          IG1=INT(((RXX1-RXFIX)/DXPA)+0.5) + 1
          JG1=INT(((RYY1-RYFIX)/DYPA)+0.5) + 1
          KG1=INT(((RZZ1-RZFIX)/DZPA)+0.5) + 1
          
          IG2=IG1 + NN1 - 1
          JG2=JG1 + NN2 - 1
          KG2=KG1 + NN3 - 1
                    
          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1

          IF (CONTROL.EQ.0.AND.IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
          END IF

          RIV1=RX1
          RIV2=RY1
          RIV3=RZ2

          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1

          IF (CONTROL.EQ.0.AND.IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
          END IF

          RIV1=RX1
          RIV2=RY2
          RIV3=RZ1
          
          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1

          IF (CONTROL.EQ.0.AND.IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
          END IF

          RIV1=RX1
          RIV2=RY2
          RIV3=RZ2

          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1

          IF (CONTROL.EQ.0.AND.IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
          END IF

          RIV1=RX2
          RIV2=RY1
          RIV3=RZ1
          
          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1

          IF (CONTROL.EQ.0.AND.IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
          END IF
          RIV1=RX2
          RIV2=RY1
          RIV3=RZ2

          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1

          IF (CONTROL.EQ.0.AND.IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
          END IF


          RIV1=RX2
          RIV2=RY2
          RIV3=RZ1
          
          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1

          IF (CONTROL.EQ.0.AND.IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
          END IF

          RIV1=RX2
          RIV2=RY2
          RIV3=RZ2

          IG3=INT(((RIV1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RIV2-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RIV3-RZFIX)/DZPA)+0.5) + 1

          IF (CONTROL.EQ.0.AND.IG3.GE.IG1.AND.IG3.LE.IG2.AND.
     &        JG3.GE.JG1.AND.JG3.LE.JG2.AND.
     &        KG3.GE.KG1.AND.KG3.LE.KG2) THEN 
            NV=NV+1
            VECINO(NV,I2)=J
            CONTROL=1
          END IF


*          END IF
          END IF
         END DO
         NVECI(I2)=NV
**         write(*,*) ir,i,nv
       END DO  


*       DO I=1,NPATCH(IR)
*         WRITE(*,*) IR,I,NVECI(I)
*         DO J=1,NVECI(I)
*           WRITE(*,*) '        ','VECINOS DE',I,'  ',VECINO(J,I)
*         END DO
*       END DO
           
*      if there is dark matter, the location of neighbours cells can
*      be reused in VEINSDM

**       NGVECI(:)=NVECI(:)


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


         CORNX1=0
         CORNX2=0
         CORNXX1=0
         CORNXX2=0
         CORNY1=0
         CORNY2=0
         CORNYY1=0
         CORNYY2=0
         CORNZ1=0
         CORNZ2=0
         CORNZZ1=0
         CORNZZ2=0

  
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
         

*         MARCA(CORNXX1:CORNXX2,CORNYY1:CORNYY2,CORNZZ1:CORNZZ2,J)=1
*         SOLAP(CORNXX1:CORNXX2,CORNYY1:CORNYY2,CORNZZ1:CORNZZ2,
*     &         IR,J)=1
*         SOLAPP(CORNXX1:CORNXX2,CORNYY1:CORNYY2,CORNZZ1:CORNZZ2,
*      &         IR,J)=I



***        celdas madre del nivel inferior         
*           DO KK=CORNZZ1,CORNZZ2
*           DO JJ=CORNYY1,CORNYY2
*           DO II=CORNXX1,CORNXX2
*              IX=II-CORNXX1+CORNX1
*              JY=JJ-CORNYY1+CORNY1
*              KZ=KK-CORNZZ1+CORNZ1
**              IF (MARCA(IX,JY,KZ,I).EQ.0) THEN 
*              IF (SOLAP(IX,JY,KZ,I).EQ.1) THEN 
**               MARCA(II,JJ,KK,J)=1
****             the cell ii,jj,kk,ir,j overlaps ix,jy,kz,ir,i              
*               SOLAP(II,JJ,KK,J)=0
*            END IF 
*           END DO
*           END DO
*           END DO

***        celdas madre del nivel inferior         
           DO KK=CORNZZ1,CORNZZ2
           DO JJ=CORNYY1,CORNYY2
           DO II=CORNXX1,CORNXX2
              IX=II-CORNXX1+CORNX1
              JY=JJ-CORNYY1+CORNY1
              KZ=KK-CORNZZ1+CORNZ1
*              CONTA2(IX,JY,KZ,I)=CONTA2(IX,JY,KZ,I)+1 
              IF (CONTA2(IX,JY,KZ,I).EQ.1) CONTA2(II,JJ,KK,J)=0
           END DO
           END DO
           END DO

       END DO
       END DO

      
      RETURN
      END  


************************************************************************
      SUBROUTINE VEINSGRID_OLD(IR,NL,NPATCH,PARE,PATCHNX,PATCHNY,
     &           PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &           CONTA2,VECINO,NVECI)
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


      INTEGER CR1,CR2,CR3,CR4,CR5,CR6
      INTEGER IR,I,J,IX,JY,KZ,II,JJ,KK
      INTEGER N1,N2,N3,L1,L2,L3,NN1,NN2,NN3,LL1,LL2,LL3
      INTEGER KK2,JJ2,II2,KZ2,JY2,IX2,NL
      INTEGER VECINO(NPALEV,NPALEV),NVECI(NPALEV)
      INTEGER NV,A2,B2,C2,K

      INTEGER CONTA2(NAMRX,NAMRY,NAMRZ,NPALEV)  !es el SOLAP
      INTEGER MARCA(NAMRX,NAMRY,NAMRZ,NPALEV)
      !!!INTEGER*2 MARCA(NAMRX,NAMRY,NAMRZ,NPALEV)
      
      REAL*4 A1,B1,C1,RIV1,RIV2,RIV3
      INTEGER CONTROL
      INTEGER CORNX1,CORNXX1,CORNX2,CORNXX2
      INTEGER CORNY1,CORNYY1,CORNY2,CORNYY2
      INTEGER CORNZ1,CORNZZ1,CORNZ2,CORNZZ2
      REAL*4 RX1,RXX1,RX2,RXX2,RY1,RYY1,RY2,RYY2
      REAL*4 RZ1,RZZ1,RZ2,RZZ2
      INTEGER LOW1, LOW2, NPALEV2

      REAL*4 DXPA,DYPA,DZPA
      REAL*4 DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ
                 
      INTEGER NUM
      COMMON /PROCESADORES/ NUM
       

       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)
       
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       
       VECINO=0
       NVECI=0



c_SUSANA---> Quito PARE y NUM de SHARED porque no se usaba!!!
!$OMP   PARALLEL DO SHARED(PATCHX,PATCHY,PATCHZ,
!$OMP+             PATCHNX,PATCHNY,PATCHNZ,VECINO,NVECI,
!$OMP+             DXPA,DYPA,DZPA,PATCHRX,PATCHRY,PATCHRZ,LOW1,LOW2),
!$OMP+  PRIVATE(I,L1,L2,L3,N1,N2,N3,NV,J,LL1,
!$OMP+    LL2,LL3,NN1,NN2,NN3,A1,B1,C1,CONTROL,RX1,RY1,RZ1,RX2,RY2,RZ2,
!$OMP+    RXX1,RYY1,RZZ1,RXX2,RYY2,RZZ2,RIV1,RIV2,RIV3)
       DO I=LOW1,LOW2
         L1=PATCHX(I)
         L2=PATCHY(I)
         L3=PATCHZ(I)

         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         NV=0
         
         RX1=PATCHRX(I)-0.5*DXPA 
         RY1=PATCHRY(I)-0.5*DYPA 
         RZ1=PATCHRZ(I)-0.5*DZPA 
         RX2=PATCHRX(I)-0.5*DXPA+(N1-1)*DXPA 
         RY2=PATCHRY(I)-0.5*DYPA+(N2-1)*DYPA 
         RZ2=PATCHRZ(I)-0.5*DZPA+(N3-1)*DZPA 
         
         DO J=LOW1, LOW2
          IF (J.NE.I) THEN
*          IF (PARE(I).EQ.PARE(J)) THEN 
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
           
          CONTROL=0

*         hay que mirar que alguno de los 8 vertices este dentro
*         vertice LL1,LL2,LL3
          RIV1=RXX1
          RIV2=RYY1
          RIV3=RZZ1
          A1=(RIV1-RX1)*(RX2-RIV1)
          B1=(RIV2-RY1)*(RY2-RIV2) 
          C1=(RIV3-RZ1)*(RZ2-RIV3)
          IF (A1.GE.0.0.AND.B1.GE.0.0.AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF

*         vertice LL1,LL2,CR6
          RIV1=RXX1
          RIV2=RYY1
          RIV3=RZZ2
          A1=(RIV1-RX1)*(RX2-RIV1)
          B1=(RIV2-RY1)*(RY2-RIV2) 
          C1=(RIV3-RZ1)*(RZ2-RIV3)
          IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.
     &       AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF

*         vertice LL1,CR5,LL3
          RIV1=RXX1
          RIV2=RYY2
          RIV3=RZZ1
          A1=(RIV1-RX1)*(RX2-RIV1)
          B1=(RIV2-RY1)*(RY2-RIV2) 
          C1=(RIV3-RZ1)*(RZ2-RIV3)
          IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.
     &       AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF

*         vertice LL1,CR5,CR6
          RIV1=RXX1
          RIV2=RYY2
          RIV3=RZZ2
          A1=(RIV1-RX1)*(RX2-RIV1)
          B1=(RIV2-RY1)*(RY2-RIV2) 
          C1=(RIV3-RZ1)*(RZ2-RIV3)
          IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.
     &       AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF

*         vertice CR4,LL2,LL3
          RIV1=RXX2
          RIV2=RYY1
          RIV3=RZZ1
          A1=(RIV1-RX1)*(RX2-RIV1)
          B1=(RIV2-RY1)*(RY2-RIV2) 
          C1=(RIV3-RZ1)*(RZ2-RIV3)
          IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.
     &       AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF


*         vertice CR4,LL2,CR6
          RIV1=RXX2
          RIV2=RYY1
          RIV3=RZZ2
          A1=(RIV1-RX1)*(RX2-RIV1)
          B1=(RIV2-RY1)*(RY2-RIV2) 
          C1=(RIV3-RZ1)*(RZ2-RIV3)
          IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.
     &       AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF

*         vertice CR4,CR5,LL3
          RIV1=RXX2
          RIV2=RYY2
          RIV3=RZZ1
          A1=(RIV1-RX1)*(RX2-RIV1)
          B1=(RIV2-RY1)*(RY2-RIV2) 
          C1=(RIV3-RZ1)*(RZ2-RIV3)
          IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.
     &       AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF


*         vertice CR4,CR5,CR6
          RIV1=RXX2
          RIV2=RYY2
          RIV3=RZZ2
          A1=(RIV1-RX1)*(RX2-RIV1)
          B1=(RIV2-RY1)*(RY2-RIV2) 
          C1=(RIV3-RZ1)*(RZ2-RIV3)
          IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.
     &       AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF


*        caso ----ll1 ---l1 -----cr1 -------cr4
          RIV1=RX1
          RIV2=RY1
          RIV3=RZ1
          A1=(RIV1-RXX1)*(RXX2-RIV1)
          B1=(RIV2-RYY1)*(RYY2-RIV2) 
          C1=(RIV3-RZZ1)*(RZZ2-RIV3)
          IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.
     &        AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF

          RIV1=RX1
          RIV2=RY1
          RIV3=RZ2
          A1=(RIV1-RXX1)*(RXX2-RIV1)
          B1=(RIV2-RYY1)*(RYY2-RIV2) 
          C1=(RIV3-RZZ1)*(RZZ2-RIV3)
          IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.
     &        AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF

          RIV1=RX1
          RIV2=RY2
          RIV3=RZ1
          A1=(RIV1-RXX1)*(RXX2-RIV1)
          B1=(RIV2-RYY1)*(RYY2-RIV2) 
          C1=(RIV3-RZZ1)*(RZZ2-RIV3)
          IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.
     &        AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF

          RIV1=RX1
          RIV2=RY2
          RIV3=RZ2
          A1=(RIV1-RXX1)*(RXX2-RIV1)
          B1=(RIV2-RYY1)*(RYY2-RIV2) 
          C1=(RIV3-RZZ1)*(RZZ2-RIV3)
          IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.
     &        AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF

          RIV1=RX2
          RIV2=RY1
          RIV3=RZ1
          A1=(RIV1-RXX1)*(RXX2-RIV1)
          B1=(RIV2-RYY1)*(RYY2-RIV2) 
          C1=(RIV3-RZZ1)*(RZZ2-RIV3)
          IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.
     &        AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF

          RIV1=RX2
          RIV2=RY1
          RIV3=RZ2
          A1=(RIV1-RXX1)*(RXX2-RIV1)
          B1=(RIV2-RYY1)*(RYY2-RIV2) 
          C1=(RIV3-RZZ1)*(RZZ2-RIV3)
          IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.
     &        AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF

          RIV1=RX2
          RIV2=RY2
          RIV3=RZ1
          A1=(RIV1-RXX1)*(RXX2-RIV1)
          B1=(RIV2-RYY1)*(RYY2-RIV2) 
          C1=(RIV3-RZZ1)*(RZZ2-RIV3)
          IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.
     &        AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF

          RIV1=RX2
          RIV2=RY2
          RIV3=RZ2
          A1=(RIV1-RXX1)*(RXX2-RIV1)
          B1=(RIV2-RYY1)*(RYY2-RIV2) 
          C1=(RIV3-RZZ1)*(RZZ2-RIV3)
          IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.
     &        AND.C1.GE.0.0) THEN
            NV=NV+1
            VECINO(NV,I)=J
            CONTROL=1
          END IF

*          END IF
          END IF
         END DO
         NVECI(I)=NV
       END DO  


*      if there is dark matter, the location of neighbours cells can
*      be reused in VEINSDM

       MARCA=0
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
                           

         DO K=1,NVECI(I)
         J=VECINO(K,I)

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

         CORNX1=0
         CORNX2=0
         CORNXX1=0
         CORNXX2=0
         CORNY1=0
         CORNY2=0
         CORNYY1=0
         CORNYY2=0
         CORNZ1=0
         CORNZ2=0
         CORNZZ1=0
         CORNZZ2=0
  
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

**        celdas madre del nivel inferior         
           DO KK=CORNZZ1,CORNZZ2
           DO JJ=CORNYY1,CORNYY2
           DO II=CORNXX1,CORNXX2
              IX=II-CORNXX1+CORNX1
              JY=JJ-CORNYY1+CORNY1
              KZ=KK-CORNZZ1+CORNZ1
              IF (MARCA(IX,JY,KZ,I).EQ.0) THEN 
               MARCA(II,JJ,KK,J)=1
**             the cell ii,jj,kk,ir,j overlaps ix,jy,kz,ir,i              
               CONTA2(II,JJ,KK,J)=0
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
      
                                                              
***************************************************************
       SUBROUTINE NOMFILE(ITER,FILNOM1,FILNOM2,FILNOM3)
***************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*9 FILNOM1,FILNOM2
       CHARACTER*10 FILNOM3
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT
      
       CONTA=0
      
       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO
      
       FILNOM1='clus'//NOM
       FILNOM2='cldm'//NOM
       FILNOM3='grids'//NOM
      
       RETURN
       END 
       
       
       
*****************************************************************
       SUBROUTINE NOMFILE2(ITER,FILE1)
*****************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*14 FILE1
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT
      
       CONTA=0
      
       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO
      
       FILE1='subhalos'//NOM
      
       RETURN
       END
       

*****************************************************************
       SUBROUTINE NOMFILE3(ITER,FILE3)
*****************************************************************
       
       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*14 FILE3
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT
                                         
       CONTA=0
             
       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO
                                                        
       FILE3='families'//NOM 
              
       RETURN
       END
       
*****************************************************************
       SUBROUTINE NOMFILE4(ITER,FILE3)
*****************************************************************
       
       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*14 FILE3
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT
                                         
       CONTA=0
             
       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO
                                                        
       FILE3='merger_t'//NOM
             
       RETURN
       END 
       
*****************************************************************
       SUBROUTINE NOMFILE5(ITER,FILE3)
*****************************************************************
       
       IMPLICIT NONE
       
       INTEGER ITER
       CHARACTER*14 FILE3
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT
                                         
       CONTA=0
             
       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO
                                                        
       FILE3='merger_r'//NOM
             
       RETURN
       END
       
*****************************************************************
       SUBROUTINE NOMFILE6(ITER,FILE6)
*****************************************************************
       
       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*15 FILE6
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT
                                         
       CONTA=0
             
       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO
                                                        
       FILE6='grid_asohf'//NOM 
            
       RETURN
       END
       
*****************************************************************
       SUBROUTINE NOMFILE7(ITER,FILE7)
*****************************************************************
       
       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*14 FILE7
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT
                                         
       CONTA=0
             
       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO
                                                        
       FILE7='dm_parti'//NOM 
              
       RETURN
       END       
       
       
!*************************************************************
!* This subroutine computes all eigenvalues and eigenvectors *
!* of a real symmetric square matrix A(N,N). On output, ele- *
!* ments of A above the diagonal are destroyed. D(N) returns *
!* the eigenvalues of matrix A. V(N,N) contains, on output,  *
!* the eigenvectors of A by columns. THe normalization to    *
!* unity is made by main program before printing results.    *
!* NROT returns the number of Jacobi matrix rotations which  *
!* were required.                                            *
!* --------------------------------------------------------- *
!* Ref.:"NUMERICAL RECIPES, Cambridge University Press, 1986,*
!*       chap. 11, pages 346-348" [BIBLI 08].                *
!*************************************************************
       Subroutine JACOBI(A,N,D,NROT)
       
       implicit none
       integer N,NROT,ip,iq,ialloc,i,j
       real*4  A(1:N,1:N),D(1:N)  
       real*4, pointer :: B(:), Z(:)
       real*4  c,g,h,s,sm,t,tau,theta,tresh

       allocate(B(1:100))   !,stat=ialloc)
       allocate(Z(1:100))   !,stat=ialloc)

       do ip=1, N
         B(ip)=A(ip,ip)
         D(ip)=B(ip)
         Z(ip)=0.d0    
       end do
       NROT=0
       do i=1, 50
         sm=0.d0
         do ip=1, N-1           !sum off-diagonal elements
            do iq=ip+1, N
               sm=sm+ABS(A(ip,iq))
            end do
         end do
         if(sm==0.d0) return    !normal return
         if(i.lt.4) then
            tresh=0.2d0*sm**2
         else
            tresh=0.d0
         end if
         do ip=1, N-1
            do iq=ip+1, N
               g=100.d0*ABS(A(ip,iq))
!      after 4 sweeps,skip the rotation if the off-diag element is small
     
       if((i.gt.4).and.(ABS(D(ip))+g.eq.ABS(D(ip)))
     &     .and.(ABS(D(iq))+g.eq.ABS(D(iq)))) then
       
       A(ip,iq)=0.d0
       else if(ABS(A(ip,iq)).gt.tresh) then       
                  h=D(iq)-D(ip)
       if(ABS(h)+g.eq.ABS(h)) then
       t=A(ip,iq)/h
       else
       theta=0.5d0*h/A(ip,iq)  
            t=1.d0/(ABS(theta)+SQRT(1.d0+theta**2))
       if(theta.lt.0.d0) t=-t
       end if
       c=1.d0/SQRT(1.d0+t**2)
       s=t*c
          tau=s/(1.d0+c)
       h=t*A(ip,iq)
       Z(ip)=Z(ip)-h
       Z(iq)=Z(iq)+h
       D(ip)=D(ip)-h
       D(iq)=D(iq)+h
       A(ip,iq)=0.d0
       do j=1, ip-1
       g=A(j,ip)
       h=A(j,iq)
       A(j,ip)=g-s*(h+g*tau)
       A(j,iq)=h+s*(g-h*tau)
          end do
       do j=ip+1, iq-1
       g=A(ip,j)
       h=A(j,iq)
       A(ip,j)=g-s*(h+g*tau)
       A(j,iq)=h+s*(g-h*tau)
          end do      
       do j=iq+1, N
       g=A(ip,j)
       h=A(iq,j)
       A(ip,j)=g-s*(h+g*tau)
       A(iq,j)=h+s*(g-h*tau)
          end do  

          NROT=NROT+1
       end if                   !if ((i.gt.4)...
       end do                    !main iq loop
       end do                    !main ip loop
       do ip=1, N
       B(ip)=B(ip)+Z(ip)
       D(ip)=B(ip)
       Z(ip)=0.d0
       end do
       end do                    !main i loop
c       pause ' 50 iterations !'
       return
       END    

!     end of file ujacobi.f90

!*************************************************************
!* Esta subrutina ordena en orden descendente los autovalores*
!* obtenidos por JACOBI  (ver diagonalizando2.f)             *
!*************************************************************
       SUBROUTINE SORT(D,N,NP)
       
       IMPLICIT NONE
       
       integer N,NP,I,J,K
       real*4  D(NP)
       real*4 P
       
       DO I=1,N-1
       K=I
       P=D(I)
       DO J=I+1, N
               IF(D(J).GE.P) THEN
                  K=J
                  P=D(J)
               END IF
       END DO
       IF(K.NE.I) THEN
       D(K)=D(I)
       D(I)=P

       END IF
       
       END DO
       RETURN
       
       END  
       
       
***********************************************************
       SUBROUTINE UNBINDING4(FAC,IFI,I,REF_MIN,REF_MAX,DISTA,
     &           U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,     
     &           NR,NMASA,NCLUSRX,NCLUSRY,NCLUSRZ,
     &           LIP,LIR,KONTA,CONTADM,VX,VY,VZ)
***********************************************************       
       
       IMPLICIT NONE
                                                            
       INCLUDE 'input_files/asohf_parameters.dat'
       
       INTEGER I,J,K,IX,IMAX,JJ,FAC
       
       REAL*4 REI,CGR,PI,PI4ROD
       COMMON /CONS/PI4ROD,REI,CGR,PI
       
       REAL*4 RETE,HTE,ROTE                                              
       COMMON /BACK/ RETE,HTE,ROTE
       
       REAL*4 REF_MIN,REF_MAX
              
*      ---HALOS Y SUBHALOS---        
       REAL*4 NR(MAXITER,NMAXNCLUS)
       REAL*4 NMASA(MAXITER,NMAXNCLUS)       
       REAL*4 NCLUSRX(MAXITER,NMAXNCLUS)
       REAL*4 NCLUSRY(MAXITER,NMAXNCLUS)
       REAL*4 NCLUSRZ(MAXITER,NMAXNCLUS)
       REAL*4 VCM2(MAXITER,NMAXNCLUS)
       REAL*4 VX(MAXITER,NMAXNCLUS)
       REAL*4 VY(MAXITER,NMAXNCLUS)
       REAL*4 VZ(MAXITER,NMAXNCLUS)

*      ---PARTICULAS E ITERACIONES--- 
c       INTEGER LIP(PARTI), LIR(PARTI)
       INTEGER LIP(PARTIRED),LIR(PARTIRED),CONTADM(PARTIRED)
       INTEGER NPART(0:NLEVELS)
       REAL*4 U2DM(PARTIRED)  
       REAL*4 U3DM(PARTIRED)  
       REAL*4 U4DM(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       REAL*4 RXPA(PARTIRED) 
       REAL*4 RYPA(PARTIRED) 
       REAL*4 RZPA(PARTIRED)                           
       
       INTEGER KONTA,IFI,KONTA3,KONTA2
c       INTEGER CONTADM(PARTI)
       
       REAL*4 DISTA(0:PARTI)

       REAL*4 VVV2,VESC,AADMX(3),AADM,DR, AA, BB, CC 
       REAL*4 BAS
       REAL*4 CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MASA
       REAL*4 POTOK

*!!!!! ESPECIAL DOBLE PRECISON !!!!!!!!!!!!!!!!!!!!!
COJO       REAL*8 POT(KONTA)
       REAL*8 POT(0:KONTA)
       REAL*8 POT1
       REAL*8 BAS8
       REAL*8 MASA8,NORMA
       REAL*8 AA8
***********************************************

       POT=0.D0

*      particulas ya ordenada!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
*      SOLO LAS PARTICULAS DE 1 A KONTA2 CONTRIBUYEN 
       KONTA2=COUNT(CONTADM(1:KONTA).EQ.0)

*      masa maxima
       NORMA=MAXVAL(DBLE(MASAP))

       MASA8=0.D0
*      calculo del potencial en doble precision
COJO       MASA8=DBLE(MASAP(LIP(1)))/NORMA
       MASA8=DBLE(MASAP(1))/NORMA
       POT(1)=MASA8/DBLE(DISTA(1))
       
       DO J=2,KONTA2
         MASA8=MASA8+ DBLE(MASAP(LIP(J)))/NORMA
         BAS8=DISTA(J)-DISTA(J-1)
         POT(J)=POT(J-1)+ MASA8*BAS8/(DBLE(DISTA(J))**2)
       END DO

       POT1=POT(KONTA2) + MASA8/REF_MAX  ! INTEGRAL DE 0 A Rmax + CONSTANTE
       AA8=NORMA*DBLE(CGR/RETE)

       BB=2.0   !1.8
       IF (FAC.EQ.1) BB=8.0
       IF (FAC.EQ.2) BB=4.0
       IF (FAC.EQ.3) BB=2.0
CX       WRITE(*,*) ' ---> factor vesc:',bb
       
*      desligando particulas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO J=1,KONTA2
       
        POTOK=(POT(J) - POT1)*AA8
        
        VESC=SQRT(2.0*ABS(POTOK))
                           
        VVV2=(U2DM(LIP(J))-VX(IFI,I))**2
     &      +(U3DM(LIP(J))-VY(IFI,I))**2
     &      +(U4DM(LIP(J))-VZ(IFI,I))**2

        VVV2=SQRT(VVV2)
       
        IF (VVV2.GT.BB*VESC)  CONTADM(J)=1
       END DO     

CX       WRITE(*,*)'PARTICULAS NO LIGADAS=', 
CX     &            COUNT(CONTADM(1:KONTA).NE.0)
       

*      NEW CENTRO DE MASAS Y VELOCIDAD 
       CALL CENTROMASAS_PART(KONTA,CONTADM,LIP,
     &           U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,     
     &           CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MASA)

       KONTA3=COUNT(CONTADM(1:KONTA).EQ.0)

       IF (KONTA3.EQ.0) THEN 

       NCLUSRX(IFI,I)=0.0
       NCLUSRY(IFI,I)=0.0
       NCLUSRZ(IFI,I)=0.0
       
       VX(IFI,I)=0.0
       VY(IFI,I)=0.0
       VZ(IFI,I)=0.0

       NMASA(IFI,I)=0.0

CX       WRITE(*,*) 'PART.LIGADAS_3=', KONTA3
       
       NR(IFI,I)=0.0
       
       ELSE

       NCLUSRX(IFI,I)=CMX
       NCLUSRY(IFI,I)=CMY
       NCLUSRZ(IFI,I)=CMZ
       
       VX(IFI,I)=VCMX
       VY(IFI,I)=VCMY
       VZ(IFI,I)=VCMZ

       NMASA(IFI,I)=MASA*9.1717e+18


*      estimacion nuevo radio

       REF_MIN=10.0E+10
       REF_MAX=-1.0
       
       DO J=1,KONTA
       IF (CONTADM(J).EQ.0) THEN
        
        AADMX(1)=RXPA(LIP(J))-CMX
        AADMX(2)=RYPA(LIP(J))-CMY
        AADMX(3)=RZPA(LIP(J))-CMZ
       
        AADM=SQRT(AADMX(1)**2+AADMX(2)**2+AADMX(3)**2)
       
        REF_MIN=MIN(REF_MIN,AADM)
        REF_MAX=MAX(REF_MAX,AADM)
       END IF
       END DO
       
       
CX       write(*,*) 'new_r',NR(IFI,I),REF_MAX
       NR(IFI,I)=REF_MAX

       END IF

       RETURN
       END

***************************************************************************
      SUBROUTINE indexx(n,arr,indx)
***************************************************************************     


      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) then 
         write(*,*) 'NSTACK too small in indexx'
         stop
        endif 
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
**********************************************************************       

**********************************************************************
       SUBROUTINE REORDENAR(KONTA,NCLUSRX,NCLUSRY,NCLUSRZ,
     &                      RXPA,RYPA,RZPA,CONTADM,LIP,LIR,DISTA)
********************************************************************** 
       
       IMPLICIT NONE
            
       INCLUDE 'input_files/asohf_parameters.dat'
       
       INTEGER I,J,K,KONTA,KONTA2
       
       
*      ---HALOS Y SUBHALOS---        
       REAL*4 NCLUSRX
       REAL*4 NCLUSRY
       REAL*4 NCLUSRZ
     
*      ---PARTICULAS E ITERACIONES--- 
c       INTEGER LIP(PARTI), LIR(PARTI)
       INTEGER LIP(PARTIRED),LIR(PARTIRED),CONTADM(PARTIRED) 

       REAL*4 RXPA(PARTIRED) 
       REAL*4 RYPA(PARTIRED) 
       REAL*4 RZPA(PARTIRED)                           
       
c       INTEGER CONTADM(PARTI)
       
       REAL*4 DISTA(0:PARTI)
       
       INTEGER INDICE(KONTA)
       REAL*4 DISTA2(0:KONTA)
       INTEGER QUIEN(KONTA)
       INTEGER QUIEN2(KONTA)

       REAL*4 AADMX(3),AADM

*      reordenar !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

       QUIEN=0
       QUIEN2=0
       DISTA=1.E10
       INDICE=0
       DISTA2=0.0
       
       
       KONTA2=0
       DO J=1, KONTA
         IF (CONTADM(J).EQ.0) THEN 
         KONTA2=KONTA2+1
         
         AADMX(1)=RXPA(LIP(J))-NCLUSRX
         AADMX(2)=RYPA(LIP(J))-NCLUSRY
         AADMX(3)=RZPA(LIP(J))-NCLUSRZ
       
         AADM=SQRT(AADMX(1)**2+AADMX(2)**2+AADMX(3)**2)
       
         DISTA(KONTA2)=AADM
         QUIEN(KONTA2)=LIP(J)
         QUIEN2(KONTA2)=LIR(J)
         END IF
       END DO

       CALL INDEXX(KONTA2,DISTA(1:KONTA2),INDICE(1:KONTA2)) !las ordena todas
                   !no solo las seleccionadas
                   

       DO J=1,KONTA2
        DISTA2(J)=DISTA(INDICE(J))
       END DO 
       
       DISTA=1.E10
       LIP=0
       LIR=0
       
       CONTADM=1
       CONTADM(1:KONTA2)=0
       
       DO J=1,KONTA2
        DISTA(J)=DISTA2(J)
        LIP(J)=QUIEN(INDICE(J))
        LIR(J)=QUIEN2(INDICE(J))
       END DO

*      estan reordenados !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       RETURN
       END
       
***********************************************************
       SUBROUTINE CENTROMASAS_PART(N,CONTADM,LIP,
     &            U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,     
     &            CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MASA)
***********************************************************       
       
       IMPLICIT NONE
                                                            
       INCLUDE 'input_files/asohf_parameters.dat'
       
       INTEGER I,J,N
       
       REAL*4 U2DM(PARTIRED)
       REAL*4 U3DM(PARTIRED)
       REAL*4 U4DM(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       REAL*4 RXPA(PARTIRED) 
       REAL*4 RYPA(PARTIRED) 
       REAL*4 RZPA(PARTIRED)                           
       
c       INTEGER LIP(PARTI),CONTADM(PARTI)
       INTEGER LIP(PARTIRED),CONTADM(PARTIRED)

       REAL*4 CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MASA
       
*      ---- DOBLE PRECISION ------------------------       
       REAL*8 CMX8,CMY8,CMZ8,VCMX8,VCMY8,VCMZ8,MASA8
       REAL*8 NORMA,BAS
*      ---------------------------------------------       
              

       CMX8=0.D0
       CMY8=0.D0
       CMZ8=0.D0

       VCMX8=0.D0
       VCMY8=0.D0
       VCMZ8=0.D0
       
       MASA8=0.D0

       NORMA=DBLE(MAXVAL(MASAP))  ! NORMALIZACION MASA
       
       DO I=1,N
        IF (CONTADM(I).EQ.0) THEN 
         J=LIP(I)

         BAS=DBLE(MASAP(J))/NORMA

         CMX8=CMX8 + DBLE(RXPA(J))*BAS
         CMY8=CMY8 + DBLE(RYPA(J))*BAS
         CMZ8=CMZ8 + DBLE(RZPA(J))*BAS
       
         VCMX8=DBLE(U2DM(J))*BAS + VCMX8
         VCMY8=DBLE(U3DM(J))*BAS + VCMY8
         VCMZ8=DBLE(U4DM(J))*BAS + VCMZ8

         MASA8=MASA8+BAS

        END IF
       END DO

       IF (MASA8.GT.0.D0) THEN 
       
       CMX8=(CMX8/MASA8)
       CMY8=(CMY8/MASA8)
       CMZ8=(CMZ8/MASA8)

       VCMX8=(VCMX8/MASA8)
       VCMY8=(VCMY8/MASA8)
       VCMZ8=(VCMZ8/MASA8)

       MASA8=MASA8*NORMA

       ELSE
       
       CMX8=0.D0
       CMY8=0.D0
       CMZ8=0.D0

       VCMX8=0.D0
       VCMY8=0.D0
       VCMZ8=0.D0

       MASA8=0.D0
       
       END IF
       
*      TRANSFORMACION   A REAL*4

       CMX=CMX8
       CMY=CMY8
       CMZ=CMZ8

       VCMX=VCMX8
       VCMY=VCMY8
       VCMZ=VCMZ8

       MASA=MASA8

       RETURN
       END
********************************************************************        
        

********************************************************************
       SUBROUTINE OVERLAPING(IFI,IR,NL,REF,ESP,BOUND,CONTA,CONTRASTEC,
     &                       NSHELL,RODO,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                       PATCHRX,PATCHRY,PATCHRZ,NX,NY,NZ,
     &                       NCLUS,MASA,RADIO,CLUSRX,CLUSRY,CLUSRZ,
     &                       REALCLUS,NSOLAP,SOLAPA,NHALLEV)
********************************************************************       
                      
       
       IMPLICIT NONE
                                                            
       INCLUDE 'input_files/asohf_parameters.dat'
       
       REAL*4 PI,ACHE,T0,RE0,PI4ROD
       COMMON /DOS/ACHE,T0,RE0                

       REAL*4 UNTERCIO,CGR,CGR2,ZI,RODO,ROI,REI,LADO,LADO0
       COMMON /CONS/PI4ROD,REI,CGR,PI
     
       REAL*4 RETE,HTE,ROTE                                              
       COMMON /BACK/ RETE,HTE,ROTE

       REAL*4 DX,DY,DZ,H2
       COMMON /ESPACIADO/ DX,DY,DZ
       
       REAL*4  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)                           
       COMMON /GRID/   RADX,RADY,RADZ        
       
*      VARIABLES
       REAL*4 U1(NMAX,NMAY,NMAZ)
       REAL*4 U1G(NMAX,NMAY,NMAZ)
       REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL*4 U11G(NAMRX,NAMRY,NAMRZ,NPALEV)  
       COMMON /VARIA/ U1,U11,U1G,U11G   !!,U12,U13,U14  !,U2,U3,U4 
      
       INTEGER NPATCH(0:NLEVELS), PARE(NPALEV)                          
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
       REAL*4 PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV) 
       
       INTEGER I,J,K,IR,NL,II,JJ,KK,IR2,IFI
       INTEGER NX,NY,NZ,KK_ENTERO,KK_ENTERO_2
       REAL*4 KK_REAL,KK_REAL_2, DIS, AA, AA2
       
       INTEGER NCLUS,SHELL,NSHELL,NHALLEV(0:NLEVELS)
       INTEGER SOLAPA(MAXNCLUS,NMAXSOLAP),NSOLAP(MAXNCLUS)
       INTEGER REALCLUS(MAXITER,MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
       
       INTEGER ICMIN,ICMAX,KONTA,KONTACAPA
       INTEGER KONTA1,KONTA2,KONTA3      
       REAL*4 RADIO(MAXNCLUS),MASA(MAXNCLUS)
       REAL*4 CMX,CMY,CMZ,SOLMAS,MASAMIN,BASMAS,MASAKK
       REAL*4 DELTA,VOL,RANT,RSHELL,DELTA2,CONTRASTEC
       REAL*4 RX2,RY2,RZ2,DXPA,DYPA,DZPA
       REAL*4 REF,ESP,VOL2,RMIN,RRRR,R111,R222
       INTEGER KK1,KK2,KK3,KK4,N1,N2,N3
       INTEGER IX,JY,KZ,CONTA(NMAX,NMAY,NMAZ) 
       REAL*4 PRUEBAX1,PRUEBAY1,PRUEBAZ1,BOUND
       REAL*4 PRUEBAX2,PRUEBAY2,PRUEBAZ2
       REAL*4 RX1,RY1,RZ1
       INTEGER NX1,NX2,NY1,NY2,NZ1,NZ2,LOW1,LOW2
       INTEGER CONTROL
       REAL*4 RIV1,RIV2,RIV3,A1,B1,C1, BOUND2


       WRITE(*,*) '== HALOES OVERLAPPING IN IR =', IR
       WRITE(*,*) 'NCLUS=', NCLUS, NHALLEV(IR)

       DXPA=0.0
       DYPA=0.0
       DZPA=0.0
       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)
       
       WRITE(*,*) 'RESOLUCION=',IR,DXPA
       
       KK2=0
       KK3=0
       IF (IR.EQ.NL) THEN
        KK2=0
        KK3=NHALLEV(IR)
       ELSE
        KK2=SUM(NHALLEV(IR+1:NL))
        KK3=NCLUS
       END IF 
       
*      La primera vez barremos todos los halos con todos! 
       DO I=KK2+1, KK3
       DO J=I+1, KK3
       
       DIS=0.0
       DIS=SQRT((CLUSRX(I)-CLUSRX(J))**2+
     &          (CLUSRY(I)-CLUSRY(J))**2+
     &          (CLUSRZ(I)-CLUSRZ(J))**2)
       
       KK_REAL=0.0
       KK_REAL=RADIO(I)+RADIO(J) 
       IF(KK_REAL.GT.DIS) THEN

          NSOLAP(I)=NSOLAP(I)+1
          NSOLAP(J)=NSOLAP(J)+1
          KK_ENTERO=0
          KK_ENTERO=NSOLAP(I)

          IF (KK_ENTERO.GT.NMAXSOLAP) THEN
          WRITE(*,*) 'WARNING!: MUCHOS SOLAPES', IFI, I, IR
          STOP
          END IF
          
          SOLAPA(J,NSOLAP(J))=I
          SOLAPA(I,NSOLAP(I))=J
       END IF
              
       END DO
       END DO
       

*      CALCULAMOS SOLAPAMIENTOS Y REDEFINIMOS
       KK2=0
       KK3=0
       IF (IR.EQ.NL) THEN
        KK2=0
        KK3=NHALLEV(IR)
       ELSE
        KK2=SUM(NHALLEV(IR+1:NL))
        KK3=NCLUS
       END IF
       
       IR2=IR
       
       DO I=KK2+1, KK3

       KONTA1=0
       KONTA2=0
       KONTA2=NSOLAP(I)
       
333    DO KK=KONTA1+1, KONTA2
       J=SOLAPA(I,KK)

       SOLMAS=0.0
       MASAMIN=0.0
       CMX=0.0
       CMY=0.0
       CMZ=0.0 
       
       IF (J.GT.0) THEN
       
       KK_ENTERO=0
       KK_ENTERO_2=0
       KK_ENTERO=REALCLUS(IFI,I)
       KK_ENTERO_2=REALCLUS(IFI,J)
       
       IF(KK_ENTERO.EQ.-1.AND.KK_ENTERO_2.EQ.-1) THEN

CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX
c       IF(IR.NE.0) THEN
CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX
     

       ! EXTENSION MAX Y MIN DE LOS HALOS:
       ! Lo calculo solo una vez tanto para IR=0 como IR>0
       PRUEBAX1=0.0
       PRUEBAX2=0.0
       PRUEBAY1=0.0
       PRUEBAY2=0.0
       PRUEBAZ1=0.0
       PRUEBAZ2=0.0

       BOUND2=0.0
       BOUND2=5.5*(RADIO(I)+RADIO(J))
       PRUEBAX1=MIN(CLUSRX(I),CLUSRX(J)) - BOUND2
       PRUEBAX2=MAX(CLUSRX(I),CLUSRX(J)) + BOUND2
       PRUEBAY1=MIN(CLUSRY(I),CLUSRY(J)) - BOUND2
       PRUEBAY2=MAX(CLUSRY(I),CLUSRY(J)) + BOUND2
       PRUEBAZ1=MIN(CLUSRZ(I),CLUSRZ(J)) - BOUND2
       PRUEBAZ2=MAX(CLUSRZ(I),CLUSRZ(J)) + BOUND2      
 
CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX
       IF(IR.NE.0) THEN
CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX


       LOW1=SUM(NPATCH(0:IR2-1))+1
       LOW2=SUM(NPATCH(0:IR2))

       DO II=LOW1,LOW2
       
       N1=PATCHNX(II) 
       N2=PATCHNY(II) 
       N3=PATCHNZ(II) 
      
*xc    --> Es necesario barrer el parche II?

       !PARCHE DE BORDE A BORDE
       RX1=0.0   !extension minima
       RY1=0.0
       RZ1=0.0
       RX2=0.0   !extension minima                                             
       RY2=0.0                                             
       RZ2=0.0            
       RX1= PATCHRX(II) - 0.5*DXPA 
       RY1= PATCHRY(II) - 0.5*DYPA 
       RZ1= PATCHRZ(II) - 0.5*DZPA 
       RX2= RX1 + (N1-1)*DXPA
       RY2= RY1 + (N2-1)*DYPA
       RZ2= RZ1 + (N3-1)*DZPA       

       CONTROL=0
*      Hay que mirar que alguno de los 8 vertices este dentro:

       !! Vertice LL1,LL2,LL3
       RIV1=RX1
       RIV2=RY1
       RIV3=RZ1
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)
       IF (A1.GE.0.0.AND.B1.GE.0.0.AND.C1.GE.0.0) CONTROL=1

       !! Vertice LL1,LL2,CR6
       RIV1=RX1
       RIV2=RY1
       RIV3=RZ2
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)
       IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.AND.
     &                      C1.GE.0.0) CONTROL=1       

       !! Vertice LL1,CR5,LL3
       RIV1=RX1
       RIV2=RY2
       RIV3=RZ1
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)
       IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.AND.
     &                      C1.GE.0.0) CONTROL=1       

       !! Vertice LL1,CR5,CR6
       RIV1=RX1
       RIV2=RY2
       RIV3=RZ2
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)
       IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.AND.
     &                      C1.GE.0.0) CONTROL=1

       !! Vertice CR4,LL2,LL3
       RIV1=RX2
       RIV2=RY1
       RIV3=RZ1
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)
       IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.AND.
     &                     C1.GE.0.0) CONTROL=1       

       !! Vertice CR4,LL2,CR6
       RIV1=RX2
       RIV2=RY1
       RIV3=RZ2
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)
       IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.AND.
     &                      C1.GE.0.0) CONTROL=1

       !!V ertice CR4,CR5,LL3
       RIV1=RX2
       RIV2=RY2
       RIV3=RZ1
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)
       IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.AND.
     &                      C1.GE.0.0) CONTROL=1

       !! Vertice CR4,CR5,CR6
       RIV1=RX2
       RIV2=RY2
       RIV3=RZ2
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)
       IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.AND.
     &                      C1.GE.0.0) CONTROL=1

       IF (CONTROL.NE.0) THEN

*xc
       DO KZ=1, N3  !!!miramos si la celda (ix,jy,kz) esta compartida
       DO JY=1, N2 
       DO IX=1, N1
              
         RX2=0.0
         RY2=0.0
         RZ2=0.0
         AA=0.0
         AA2=0.0
         
         RX2=PATCHRX(II) - 0.5*DXPA + (IX-1)*DXPA
         RY2=PATCHRY(II) - 0.5*DYPA + (JY-1)*DYPA
         RZ2=PATCHRZ(II) - 0.5*DZPA + (KZ-1)*DZPA
  
         AA=SQRT((CLUSRX(I)-RX2)**2+(CLUSRY(I)-RY2)**2
     &          +(CLUSRZ(I)-RZ2)**2)
     
         AA2=SQRT((CLUSRX(J)-RX2)**2+(CLUSRY(J)-RY2)**2
     &           +(CLUSRZ(J)-RZ2)**2)
       
        KK_REAL=0.0
        KK_REAL_2=0.0
        KK_REAL=RADIO(I)
        KK_REAL_2=RADIO(J)
       
        IF (AA.LT.KK_REAL.AND.AA2.LT.KK_REAL_2) THEN
          SOLMAS=SOLMAS+U11(IX,JY,KZ,II)*DXPA*DYPA*DZPA
        END IF
       
       END DO
       END DO
       END DO 

       END IF  !!CONTROL

       END DO  !!PARCHE II

CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX
       ELSE
CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX

*      posible "extension" maxima del halo 
       NX1=0
       NX2=0
       NY1=0
       NY2=0
       NZ1=0
       NZ2=0
       
       NX1=INT(((PRUEBAX1-RADX(1))/DX)+0.5)+1
       IF (NX1.LT.1) NX1=1  

       NX2=INT(((PRUEBAX2-RADX(1))/DX)+0.5)+1
       IF (NX2.GT.NX) NX2=NX

       NY1=INT(((PRUEBAY1-RADY(1))/DY)+0.5)+1
       IF (NY1.LT.1) NY1=1

       NY2=INT(((PRUEBAY2-RADY(1))/DY)+0.5)+1
       IF (NY2.GT.NY) NY2=NY

       NZ1=INT(((PRUEBAZ1-RADZ(1))/DZ)+0.5)+1
       IF (NZ1.LT.1) NZ1=1 

       NZ2=INT(((PRUEBAZ2-RADZ(1))/DZ)+0.5)+1
       IF (NZ2.GT.NZ) NZ2=NZ

         DO KZ=NZ1, NZ2
         DO JY=NY1, NY2
         DO IX=NX1, NX2 

           AA=0.0
           AA=SQRT((RADX(IX)-CLUSRX(I))**2+
     &             (RADY(JY)-CLUSRY(I))**2+
     &             (RADZ(KZ)-CLUSRZ(I))**2)
     
           AA2=0.0
           AA2=SQRT((RADX(IX)-CLUSRX(J))**2+
     &              (RADY(JY)-CLUSRY(J))**2+
     &              (RADZ(KZ)-CLUSRZ(J))**2)
             
           KK_REAL=0.0
           KK_REAL_2=0.0
           KK_REAL=RADIO(I)
           KK_REAL_2=RADIO(J) 
           IF (AA.LT.KK_REAL.AND.AA2.LT.KK_REAL_2) THEN
                SOLMAS=SOLMAS+U1(IX,JY,KZ)*DX*DY*DZ
           END IF
       
          END DO
          END DO
          END DO

       
CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX
       END IF     !IR.NE.0
CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX   

       END IF  !!IF(KK_ENTERO.EQ.-1.AND.KK_ENTERO_2.EQ.-1) 

       
       SOLMAS=SOLMAS*RODO*(RE0**3)

       KK_REAL=0.0
       KK_REAL_2=0.0
       KK_REAL=MASA(I)
       KK_REAL_2=MASA(J)
       
       IF(KK_REAL.LT.KK_REAL_2) THEN
        MASAMIN=MASA(I)
        ICMIN=I
        ICMAX=J
       ELSE
        MASAMIN=MASA(J)
        ICMIN=J
        ICMAX=I
       END IF
                     

*      Caso 1: redefinicion de halos solapados
*      EN ESTE CASO REDEFININMOS UN UNICO HALO
       
       KK_ENTERO=0
       KK_ENTERO_2=0
       KK_ENTERO=REALCLUS(IFI,ICMAX)
       KK_ENTERO_2=REALCLUS(IFI,ICMIN)
       
       IF(SOLMAS.GT.0.5*MASAMIN.AND.SOLMAS.LT.0.8*MASAMIN.AND.
     &    KK_ENTERO.EQ.-1.AND.KK_ENTERO_2.EQ.-1)THEN
     
       MASA(ICMIN)=MASA(ICMIN)-SOLMAS
       
       KK_REAL=0.0
       KK_REAL=MASA(ICMIN)
       IF(KK_REAL.LT.0.0) MASA(ICMIN)=0.0


*      REDEFINIMOS LOS CUMULOS SOLAPADOS
       CMX=0.0
       CMY=0.0
       CMZ=0.0
       
       CMX=CLUSRX(I)*MASA(I)+CLUSRX(J)*MASA(J)
       CMY=CLUSRY(I)*MASA(I)+CLUSRY(J)*MASA(J)
       CMZ=CLUSRZ(I)*MASA(I)+CLUSRZ(J)*MASA(J)
       
       CMX=CMX/(MASA(I)+MASA(J))
       CMY=CMY/(MASA(I)+MASA(J))
       CMZ=CMZ/(MASA(I)+MASA(J))
       
*      Se elimina el cumulo de menor masa: ICMIN
       REALCLUS(IFI,ICMIN)=0

*      Se redefinen las coordenadas del de mayor masa: ICMAX
       CLUSRX(ICMAX)=CMX
       CLUSRY(ICMAX)=CMY
       CLUSRZ(ICMAX)=CMZ
       
*      Se extienden capas entorno al nuevo cumulo ICMAX
       RMIN=0.0       
       DELTA=0.0
       AA=0.0
       RANT=0.0
              
       SHELL=0
       BASMAS=0.0
       DELTA2=0.0
       VOL2=0.0
       RSHELL=0.0
                     
       DELTA2=10.0*CONTRASTEC*ROTE         

       DO WHILE(DELTA2.GT.(CONTRASTEC)*ROTE)
             
       IF(SHELL<NSHELL) THEN
                 
          SHELL=SHELL+1
          DELTA=0.0
          RANT=RSHELL           
          RSHELL=0.0
          VOL=0.0
       ELSE
          WRITE(*,*) 'WARNING:DEMASIADAS CAPAS'
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


*      EXTENDEMOS...
       
CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX
CX       IF(IR.NE.0) THEN
CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX
             
       !Region aprox. ocupada por el halo
       PRUEBAX1=0.0
       PRUEBAX2=0.0
       PRUEBAY1=0.0
       PRUEBAY2=0.0
       PRUEBAZ1=0.0
       PRUEBAZ2=0.0

       BOUND2=0.0
       BOUND2=5.5*(RADIO(ICMAX)+RADIO(ICMIN))
       PRUEBAX1=CLUSRX(ICMAX) - BOUND2
       PRUEBAX2=CLUSRX(ICMAX) + BOUND2
       PRUEBAY1=CLUSRY(ICMAX) - BOUND2
       PRUEBAY2=CLUSRY(ICMAX) + BOUND2
       PRUEBAZ1=CLUSRZ(ICMAX) - BOUND2
       PRUEBAZ2=CLUSRZ(ICMAX) + BOUND2

CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX
       IF(IR.NE.0) THEN
CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX

       LOW1=SUM(NPATCH(0:IR2-1))+1
       LOW2=SUM(NPATCH(0:IR2))
       
       DO II=LOW1,LOW2
       
       N1=PATCHNX(II)
       N2=PATCHNY(II)
       N3=PATCHNZ(II)

*xc    --> Es necesario barrer el parche II?
       RX1=0.0   !extension minima
       RY1=0.0
       RZ1=0.0
       RX2=0.0   !extension minima                                             
       RY2=0.0
       RZ2=0.0
       RX1=PATCHRX(II) - 0.5*DXPA
       RY1=PATCHRY(II) - 0.5*DYPA
       RZ1=PATCHRZ(II) - 0.5*DZPA
       RX2= RX1 + (N1-1)*DXPA
       RY2= RY1 + (N2-1)*DYPA
       RZ2= RZ1 + (N3-1)*DZPA

       CONTROL=0
*      Hay que mirar que alguno de los 8 vertices este dentro:

       !! Vertice LL1,LL2,LL3
       RIV1=RX1
       RIV2=RY1
       RIV3=RZ1
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)
       IF (A1.GE.0.0.AND.B1.GE.0.0.AND.C1.GE.0.0) CONTROL=1

       !! Vertice LL1,LL2,CR6
       RIV1=RX1
       RIV2=RY1
       RIV3=RZ2
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)
       IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.AND.
     &                      C1.GE.0.0) CONTROL=1

       !! Vertice LL1,CR5,LL3
       RIV1=RX1
       RIV2=RY2
       RIV3=RZ1
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)
       IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.AND.
     &                      C1.GE.0.0) CONTROL=1

       !! Vertice LL1,CR5,CR6
       RIV1=RX1
       RIV2=RY2
       RIV3=RZ2
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)
       IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.AND.
     &                      C1.GE.0.0) CONTROL=1

       !! Vertice CR4,LL2,LL3
       RIV1=RX2
       RIV2=RY1
       RIV3=RZ1
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)
       IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.AND.
     &                     C1.GE.0.0) CONTROL=1

       !! Vertice CR4,LL2,CR6
       RIV1=RX2
       RIV2=RY1
       RIV3=RZ2
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)
       IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.AND.
     &                      C1.GE.0.0) CONTROL=1


       !!Vertice CR4,CR5,LL3
       RIV1=RX2
       RIV2=RY2
       RIV3=RZ1
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)       
       IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.AND.
     &                      C1.GE.0.0) CONTROL=1

       !! Vertice CR4,CR5,CR6
       RIV1=RX2
       RIV2=RY2
       RIV3=RZ2
       A1=(RIV1-PRUEBAX1)*(PRUEBAX2-RIV1)
       B1=(RIV2-PRUEBAY1)*(PRUEBAY2-RIV2)
       C1=(RIV3-PRUEBAZ1)*(PRUEBAZ2-RIV3)
       IF (CONTROL.EQ.0.AND.A1.GE.0.0.AND.B1.GE.0.0.AND.
     &                      C1.GE.0.0) CONTROL=1       

       IF (CONTROL.NE.0) THEN

*xc

       DO KZ=1, N3 
       DO JY=1, N2 
       DO IX=1, N1 
       
         RX2=0.0
         RY2=0.0
         RZ2=0.0
         
         RX2=PATCHRX(II) - 0.5*DXPA + (IX-1)*DXPA
         RY2=PATCHRY(II) - 0.5*DYPA + (JY-1)*DYPA
         RZ2=PATCHRZ(II) - 0.5*DZPA + (KZ-1)*DZPA
  
         AA=0.0
         AA=SQRT((CLUSRX(ICMAX)-RX2)**2+(CLUSRY(ICMAX)-RY2)**2
     &          +(CLUSRZ(ICMAX)-RZ2)**2)
         
          IF(AA.GE.R222.AND.AA.LE.R111) THEN     
           DELTA=DELTA+U11(IX,JY,KZ,II)*DXPA*DYPA*DZPA
         END IF
       
       END DO
       END DO
       END DO   


       END IF   !! EXTENSION PARCHE Y HALOS

       END DO   !! PARCHE II

CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX       
       ELSE
CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX
*      Nivel base

*      posible "ALCANCE" del halo 
       NX1=0
       NX2=0
       NY1=0
       NY2=0
       NZ1=0
       NZ2=0
       
       NX1=INT(((PRUEBAX1-RADX(1))/DX)+0.5)+1
       IF (NX1.LT.1) NX1=1
       
       NX2=INT(((PRUEBAX2-RADX(1))/DX)+0.5)+1
       IF (NX2.GT.NX) NX2=NX
       
       NY1=INT(((PRUEBAY1-RADY(1))/DY)+0.5)+1
       IF (NY1.LT.1) NY1=1
       
       NY2=INT(((PRUEBAY2-RADY(1))/DY)+0.5)+1
       IF (NY2.GT.NY) NY2=NY
       
       NZ1=INT(((PRUEBAZ1-RADZ(1))/DZ)+0.5)+1
       IF (NZ1.LT.1) NZ1=1
       
       NZ2=INT(((PRUEBAZ2-RADZ(1))/DZ)+0.5)+1
       IF (NZ2.GT.NZ) NZ2=NZ
       
      
       DO KZ=NZ1, NZ2
       DO JY=NY1, NY2
       DO IX=NX1, NX2

       AA=0.0
       AA=SQRT((RADX(IX)-CLUSRX(ICMAX))**2+
     &         (RADY(JY)-CLUSRY(ICMAX))**2+
     &         (RADZ(KZ)-CLUSRZ(ICMAX))**2)
  
       IF(AA.GE.R222.AND.AA.LE.R111) THEN 
         DELTA=DELTA+U1(IX,JY,KZ)*DX*DY*DZ
         CONTA(IX,JY,KZ)=1                      
       END IF
       
       END DO
       END DO
       END DO 


CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX       
       END IF
CXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCXCX

       BASMAS=BASMAS+DELTA*RODO*RE0**3
       VOL2=PI*(4.0/3.0)*(RSHELL*RETE)**3
       DELTA2=BASMAS/VOL2              
       DELTA=(DELTA*RODO*RE0**3)/VOL  
       RMIN=RSHELL
                    
       END DO      !END DO WHILE DE CRECER EL HALO (DELTA2)
       
       RADIO(ICMAX)=RSHELL
       MASA(ICMAX)=BASMAS
       
  
*      HAY QUE REDEFINIR LOS SOLAPES DEL NUEVO ICMAX 
       KK1=0
       KK4=0
       IF (IR.EQ.NL) THEN
        KK1=0
        KK4=NHALLEV(IR)
       ELSE
        KK1=SUM(NHALLEV(IR+1:NL))
        KK4=NCLUS
       END IF 
       
       DO K=KK1+1, KK4     !todos los halos del nivel
       
       KONTA=1  !se supone que no solapan  
       DO JJ=1, NSOLAP(ICMAX)
        KK_ENTERO=0
        KK_ENTERO=SOLAPA(ICMAX,JJ)
        IF (KK_ENTERO.EQ.K) KONTA=0  !ya solapaban
        IF (KONTA.EQ.0) EXIT
       END DO
       IF (K.EQ.ICMAX) KONTA=0  !ICMAX solapa con el mismo
       
       IF (KONTA.EQ.1) THEN     !K es un posible nuevo solape de ICMAX
       
       DIS=0.0
       DIS=SQRT((CLUSRX(ICMAX)-CLUSRX(K))**2+
     &          (CLUSRY(ICMAX)-CLUSRY(K))**2+
     &          (CLUSRZ(ICMAX)-CLUSRZ(K))**2)
       
       KK_REAL=0.0
       KK_REAL=RADIO(ICMAX)+RADIO(K)
        
       IF(KK_REAL.GT.DIS) THEN
       
          NSOLAP(ICMAX)=NSOLAP(ICMAX)+1
          NSOLAP(K)=NSOLAP(K)+1

          KK_ENTERO=0
          KK_ENTERO_2=0
          KK_ENTERO=NSOLAP(ICMAX)
          KK_ENTERO_2=NSOLAP(K)
          
          IF (KK_ENTERO.GT.NMAXSOLAP) THEN
           WRITE(*,*) 'WARNING! MUCHOS SOLAPES', IFI,I,IR
           STOP
          END IF 
          IF (KK_ENTERO_2.GT.NMAXSOLAP) THEN
           WRITE(*,*) 'WARNING! MUCHOS SOLAPES', IFI,I,IR
           STOP
          END IF
          
          SOLAPA(ICMAX,NSOLAP(ICMAX))=K
          SOLAPA(K,NSOLAP(K))=ICMAX

       END IF   !!KK_REAL
              
       END IF   !!KONTA
       
       END DO

****   Fin nuevos solapes de ICMAX

       END IF    !caso 1: redefinicion de halos solapados
       
*      Caso 2: se elimina un halo menor solapado por otro mayor      
*      elimino el cumulo de menor masa: ICMIN       
       
       KK_ENTERO=0
       KK_ENTERO_2=0
       KK_ENTERO=REALCLUS(IFI,ICMAX)
       KK_ENTERO_2=REALCLUS(IFI,ICMIN)
       IF(SOLMAS.GE.0.8*MASAMIN.AND.KK_ENTERO.EQ.-1.AND.
     &    KK_ENTERO_2.EQ.-1)THEN
       
          REALCLUS(IFI,ICMIN)=0 
       END IF    !caso 2
       
       END IF   !J.GT.0
       
       END DO    !NSOLAP(I)=KONTA          
       
       KK_ENTERO=0
       KK_ENTERO=NSOLAP(I)
       IF (KONTA2.LT.KK_ENTERO) THEN
        KONTA1=KONTA2
        KONTA2=NSOLAP(I)
       GOTO 333
       END IF
       
       END DO    !HALO I
       

*      FIN CORRECION DE SOLAPES
        

       RETURN
       END 
       
************************************************************************
       SUBROUTINE MESHRENOEF(ITER,NX,NY,NZ,NL,COTA,NPATCH,
     &                   NPART,PATCHNX,PATCHNY,PATCHNZ,
     &                   PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,
     &                   PATCHRZ,PARE,U2DM,U3DM,U4DM,MASAP,MAP,
     &                   RXPA,RYPA,RZPA,ZETA,T,LADO0,FLAG_MASCLET,PLOT)
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
