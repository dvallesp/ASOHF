**********************************************************************
       SUBROUTINE RE_SORT_HALOES(NCLUS,NHALLEV,REALCLUS,CLUSRX,CLUSRY,
     &                           CLUSRZ,RADIO,MASA,LEVHAL,PATCHCLUS,
     &                           DMPCLUS)
**********************************************************************
*      Resorts the clusters, getting rid of REALCLUS=0 ones
**********************************************************************
        IMPLICIT NONE
        INCLUDE 'input_files/asohf_parameters.dat'

        INTEGER NCLUS
        INTEGER NHALLEV(0:NLEVELS)
        INTEGER REALCLUS(MAXNCLUS)
        REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
        REAL*4 MASA(MAXNCLUS), RADIO(MAXNCLUS)
        INTEGER LEVHAL(MAXNCLUS),PATCHCLUS(MAXNCLUS)
        INTEGER DMPCLUS(NMAXNCLUS)

        INTEGER I,J
        INTEGER,ALLOCATABLE::RESORT(:)

        NHALLEV(:)=0
        J=0
        ALLOCATE(RESORT(NCLUS))
        DO I=1,NCLUS
         IF (REALCLUS(I).EQ.0) CYCLE
         J=J+1
         RESORT(I)=J
         CLUSRX(J)=CLUSRX(I)
         CLUSRY(J)=CLUSRY(I)
         CLUSRZ(J)=CLUSRZ(I)
         RADIO(J)=RADIO(I)
         MASA(J)=MASA(I)
         LEVHAL(J)=LEVHAL(I)
         PATCHCLUS(J)=PATCHCLUS(I)
         DMPCLUS(J)=DMPCLUS(I)
         IF (REALCLUS(I).LE.0) THEN
          REALCLUS(J)=REALCLUS(I)
         ELSE
          REALCLUS(J)=RESORT(REALCLUS(I))
         END IF
         NHALLEV(LEVHAL(J))=NHALLEV(LEVHAL(J))+1
        END DO

        NCLUS=J
        DEALLOCATE(RESORT)

        RETURN
        END


**********************************************************************
        SUBROUTINE HALOES_BORDER(NCLUS,CLUSRX,CLUSRY,CLUSRZ,RADIO,
     &                           LADO0,HALBORDERS)
**********************************************************************
*       Identifies the haloes at the border of the box for special
*        treatment if necessary
**********************************************************************

        IMPLICIT NONE

        INCLUDE 'input_files/asohf_parameters.dat'

        INTEGER NCLUS
        REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
        REAL*4 RADIO(MAXNCLUS)
        REAL LADO0
        INTEGER HALBORDERS(MAXNCLUS)

        INTEGER I
        REAL XL,XR,X1,X2,Y1,Y2,Z1,Z2,XC,YC,ZC,RC

        XL=-LADO0/2
        XR=LADO0/2

!$OMP PARALLEL DO SHARED(NCLUS,HALBORDERS),PRIVATE(I),DEFAULT(NONE)
        DO I=1,NCLUS
         HALBORDERS(I)=0
        END DO

!$OMP PARALLEL DO SHARED(NCLUS,RADIO,CLUSRX,CLUSRY,CLUSRZ,XL,XR,
!$OMP+                   HALBORDERS),
!$OMP+            PRIVATE(I,RC,XC,YC,ZC,X1,X2,Y1,Y2,Z1,Z2),
!$OMP+            DEFAULT(NONE)
        DO I=1,NCLUS
         RC=RADIO(I)
         XC=CLUSRX(I)
         YC=CLUSRY(I)
         ZC=CLUSRZ(I)
         X1=XC-1.25*RC
         X2=XC+1.25*RC
         Y1=YC-1.25*RC
         Y2=YC+1.25*RC
         Z1=ZC-1.25*RC
         Z2=ZC+1.25*RC
         IF (X1.LT.XL.OR.X2.GT.XR.OR.
     &       Y1.LT.XL.OR.Y2.GT.XR.OR.
     &       Z1.LT.XL.OR.Z2.GT.XR) THEN
          HALBORDERS(I)=1
         END IF
        END DO

        RETURN
        END


**********************************************************************
        SUBROUTINE COUNT_PARTICLES_HALO(RXPA,RYPA,RZPA,N_DM,
     &                                  CMX,CMY,CMZ,R,NUM_PART_HAL)
**********************************************************************
*       Counts number of particles inside a halo
**********************************************************************
         IMPLICIT NONE
         INCLUDE 'input_files/asohf_parameters.dat'

         REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED)
         INTEGER N_DM
         REAL CMX,CMY,CMZ,R
         INTEGER NUM_PART_HAL

         INTEGER I
         REAL BASR,R2

         NUM_PART_HAL=0
         R2=R**2

         DO I=1,N_DM
          BASR=(RXPA(I)-CMX)**2 + (RYPA(I)-CMY)**2 + (RZPA(I)-CMZ)**2
          IF (BASR.LT.R2) NUM_PART_HAL=NUM_PART_HAL+1
         END DO

         RETURN
         END

**********************************************************************
       SUBROUTINE PRUNE_POOR_HALOES(NCLUS,CLUSRX,CLUSRY,CLUSRZ,RADIO,
     &                              REALCLUS,RXPA,RYPA,RZPA,N_DM,
     &                              MIN_NUM_PART,DMPCLUS,FACRAD,
     &                              DO_NEED_TO_COUNT)
**********************************************************************
*       Removes haloes with too few particles
**********************************************************************

        IMPLICIT NONE
        INCLUDE 'input_files/asohf_parameters.dat'

        INTEGER NCLUS
        REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
        REAL*4 RADIO(MAXNCLUS)
        INTEGER REALCLUS(MAXNCLUS)
        REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED)
        INTEGER N_DM,MIN_NUM_PART
        INTEGER DMPCLUS(NMAXNCLUS)
        REAL FACRAD
        INTEGER DO_NEED_TO_COUNT

        INTEGER I,BASINT,KONTA2
        REAL CMX,CMY,CMZ,RR

        KONTA2=0


        IF (DO_NEED_TO_COUNT.EQ.1) THEN
!$OMP PARALLEL DO SHARED(NCLUS,CLUSRX,CLUSRY,CLUSRZ,RADIO,N_DM,
!$OMP+                   MIN_NUM_PART,REALCLUS,RXPA,RYPA,RZPA,FACRAD,
!$OMP+                   DMPCLUS),
!$OMP+            PRIVATE(I,CMX,CMY,CMZ,RR,BASINT),
!$OMP+            REDUCTION(+:KONTA2)
!$OMP+            DEFAULT(NONE)
****************************
         DO I=1,NCLUS
****************************
          CMX=CLUSRX(I)
          CMY=CLUSRY(I)
          CMZ=CLUSRZ(I)
          RR=FACRAD*RADIO(I)

          CALL COUNT_PARTICLES_HALO(RXPA,RYPA,RZPA,N_DM,CMX,CMY,CMZ,RR,
     &                              BASINT)

          IF (BASINT.LT.MIN_NUM_PART) THEN
           REALCLUS(I)=0
           KONTA2=KONTA2+1
           !WRITE(*,*) I,BASINT,RR
          END IF

          DMPCLUS(I)=BASINT

*****************************
         END DO        !I
*****************************
         WRITE(*,*)'CHECKING POOR HALOS----->', KONTA2
        ELSE
!$OMP  PARALLEL DO SHARED(NCLUS,DMPCLUS,MIN_NUM_PART,
!$OMP+             REALCLUS),PRIVATE(I,BASINT),
!$OMP+             REDUCTION(+:KONTA2)
         DO I=1, NCLUS
          BASINT=DMPCLUS(I)
          IF (BASINT.LT.MIN_NUM_PART) THEN
           REALCLUS(I)=0
           KONTA2=KONTA2+1
          END IF
         END DO
         WRITE(*,*)'RE-CHECKING POOR HALOS----->', KONTA2
        END IF

        RETURN
        END

**********************************************************************
       SUBROUTINE CHECK_RUBISH(NCLUS,REALCLUS,CLUSRX,CLUSRY,CLUSRZ,VX,
     &                         VY,VZ,MASA,RADIO,LEVHAL)
**********************************************************************
*      Looks for overlapping haloes
**********************************************************************
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NCLUS
       INTEGER REALCLUS(MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
       REAL*4 VX(NMAXNCLUS),VY(NMAXNCLUS),VZ(NMAXNCLUS)
       REAL*4 RADIO(MAXNCLUS),MASA(MAXNCLUS)
       INTEGER LEVHAL(MAXNCLUS)

       INTEGER KONTA,I,J
       REAL A1,A2,A3,DIS,VVV1,VVV2

       KONTA=0

       DO I=1, NCLUS
        IF (REALCLUS(I).NE.0) THEN
         DO J=1, NCLUS
          IF (REALCLUS(J).NE.0.AND.LEVHAL(J).GT.
     &        LEVHAL(I)) THEN

           DIS=SQRT((CLUSRX(I)-CLUSRX(J))**2+
     &              (CLUSRY(I)-CLUSRY(J))**2+
     &              (CLUSRZ(I)-CLUSRZ(J))**2)

           A1=MIN(MASA(I),MASA(J))/MAX(MASA(I),MASA(J))

           VVV1=SQRT(VX(I)**2+VY(I)**2+VZ(I)**2)
           VVV2=SQRT(VX(J)**2+VY(J)**2+VZ(J)**2)
           IF (MIN(VVV1,VVV2).NE.0.0) THEN
            A2=(ABS(VVV1-VVV2))/MAX(VVV1,VVV2)
           END IF

           A3=MIN(RADIO(I),RADIO(J))

           IF (DIS.LT.1.01*A3.AND.A1.GT.0.2.AND.A2.LT.3.0) THEN
            IF (MASA(I).GT.MASA(J)) THEN
             REALCLUS(J)=0
             KONTA=KONTA+1
            ELSE
             REALCLUS(I)=0
             KONTA=KONTA+1
            END IF
           END IF

          END IF   !realclus
         END DO
        END IF   !realclus
       END DO

       WRITE(*,*)'CHECKING RUBBISH----->', KONTA


       RETURN
       END

**********************************************************************
       SUBROUTINE ACCIDENTAL_SUBSTRUCTURE(NCLUS,REALCLUS,CLUSRX,CLUSRY,
     &                                    CLUSRZ,VX,VY,VZ,MASA,RADIO,
     &                                    LEVHAL)
**********************************************************************
*      Looks for overlapping haloes
**********************************************************************
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NCLUS
       INTEGER REALCLUS(MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
       REAL*4 VX(NMAXNCLUS),VY(NMAXNCLUS),VZ(NMAXNCLUS)
       REAL*4 RADIO(MAXNCLUS),MASA(MAXNCLUS)
       INTEGER LEVHAL(MAXNCLUS)

       INTEGER KONTA,I,J
       REAL A1,A2,A3,DIS,VVV1,VVV2

       KONTA=0

       DO I=1, NCLUS
        IF (REALCLUS(I).EQ.-1) THEN
         DO J=1, NCLUS
          IF (REALCLUS(J)==-1) THEN
           DIS=SQRT((CLUSRX(I)-CLUSRX(J))**2+
     &              (CLUSRY(I)-CLUSRY(J))**2+
     &              (CLUSRZ(I)-CLUSRZ(J))**2)
           IF(LEVHAL(J).GT.LEVHAL(I).AND.
     &        DIS+RADIO(J).LE.1.0*RADIO(I)) THEN
            REALCLUS(J)=0
            KONTA=KONTA+1

           END IF
          END IF

         END DO
        END IF
       END DO

       WRITE(*,*)'CHECKING ACCIDENTAL SUBSTRUCTURE----->', KONTA

       RETURN
       END

**********************************************************************
        SUBROUTINE RECENTER_DENSITY_PEAK_PARTICLES(CX,CY,CZ,R,RXPA,RYPA,
     &                                             RZPA,MASAP,N_DM,
     &                                             DXPAMIN,MAX_NUM_PART)
**********************************************************************
*       Recenters density peak using particles
**********************************************************************
        IMPLICIT NONE
        INCLUDE 'input_files/asohf_parameters.dat'

        REAL CX,CY,CZ,R
        REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED),
     &         MASAP(PARTIRED),DXPAMIN
        INTEGER N_DM,MAX_NUM_PART

        INTEGER KONTA,FLAG_LARGER,I,NN,IX,JY,KZ,IP,MAX_NUM_PART_LOCAL
        INTEGER INMAX(3),KONTA2,FLAG_ITER,NUMPARTMIN,WELL_ALLOCATED
        REAL RADIO,BAS,XL,YL,ZL,DDXX,BASX,BASY,BASZ
        REAL,ALLOCATABLE::DENS(:,:,:)
        INTEGER,ALLOCATABLE::LIP(:)

        NUMPARTMIN=32 !4**3/2

        MAX_NUM_PART_LOCAL=MAX_NUM_PART
        WELL_ALLOCATED=0

        DO WHILE (WELL_ALLOCATED.EQ.0)
         IF (ALLOCATED(LIP)) DEALLOCATE(LIP)
         ALLOCATE(LIP(MAX_NUM_PART_LOCAL))
         WELL_ALLOCATED=1
         write(*,*) cx,cy,cz,'allocated with',max_num_part_local

         FLAG_LARGER=1
         RADIO=MAX(0.05*R,DXPAMIN)
         DO WHILE (FLAG_LARGER.EQ.1)
          KONTA=0
          DO I=1,N_DM
           IF (CX-RADIO.LT.RXPA(I).AND.RXPA(I).LT.CX+RADIO.AND.
     &        CY-RADIO.LT.RYPA(I).AND.RYPA(I).LT.CY+RADIO.AND.
     &        CZ-RADIO.LT.RZPA(I).AND.RZPA(I).LT.CZ+RADIO) THEN
            KONTA=KONTA+1
            IF (KONTA.GT.MAX_NUM_PART_LOCAL) THEN
             WELL_ALLOCATED=0
             EXIT
            END IF
            LIP(KONTA)=I
           END IF
          END DO !I=1,N_DM

          IF (WELL_ALLOCATED.EQ.0) EXIT

          IF (KONTA.GT.NUMPARTMIN) THEN
           FLAG_LARGER=0
          ELSE
           RADIO=RADIO*1.2
          END IF
         END DO ! WHILE (FLAG_LARGER.EQ.1)

         IF (WELL_ALLOCATED.EQ.0) THEN
          MAX_NUM_PART_LOCAL=2*MAX_NUM_PART_LOCAL
         END IF

        END DO ! WHILE (WELL_ALLOCATED.EQ.0)
        write(*,*) cx,cy,cz,'good allocation',max_num_part_local,konta

c        WRITE(*,*) RADIO,KONTA,CX,CY,CZ
        NN=4

        ALLOCATE(DENS(NN,NN,NN))
        DDXX=2.0*RADIO/FLOAT(NN)
        XL=CX-RADIO
        YL=CY-RADIO
        ZL=CZ-RADIO
        FLAG_ITER=1
        DO WHILE (FLAG_ITER.EQ.1)
         DO KZ=1,NN
         DO JY=1,NN
         DO IX=1,NN
          DENS(IX,JY,KZ)=0.0
         END DO
         END DO
         END DO

         DO I=1,KONTA
          IP=LIP(I)
          IX=INT((RXPA(IP)-XL)/DDXX)+1
          JY=INT((RYPA(IP)-YL)/DDXX)+1
          KZ=INT((RZPA(IP)-ZL)/DDXX)+1
          !IF (JY.EQ.0) WRITE(*,*) (RYPA(IP)-YL)/DDXX
          DENS(IX,JY,KZ)=DENS(IX,JY,KZ)+MASAP(IP)
         END DO

         INMAX=MAXLOC(DENS)
         IX=INMAX(1)
         JY=INMAX(2)
         KZ=INMAX(3)
         CX=XL+(IX-0.5)*DDXX
         CY=YL+(JY-0.5)*DDXX
         CZ=ZL+(KZ-0.5)*DDXX
         RADIO=RADIO/2.0
         XL=CX-RADIO
         YL=CY-RADIO
         ZL=CZ-RADIO
         DDXX=DDXX/2.0

         KONTA2=0
         DO I=1,KONTA
          IP=LIP(I)
          IF (CX-RADIO.LT.RXPA(IP).AND.RXPA(IP).LT.CX+RADIO.AND.
     &        CY-RADIO.LT.RYPA(IP).AND.RYPA(IP).LT.CY+RADIO.AND.
     &        CZ-RADIO.LT.RZPA(IP).AND.RZPA(IP).LT.CZ+RADIO) THEN
           KONTA2=KONTA2+1
           LIP(KONTA2)=IP
          END IF
         END DO

         KONTA=KONTA2

         IF (KONTA2.LT.NUMPARTMIN) FLAG_ITER=0

c         WRITE(*,*) RADIO,KONTA2,CX,CY,CZ,DENS(IX,JY,KZ)/DDXX**3,
c     &              IX,JY,KZ,FLAG_ITER
        END DO

        DEALLOCATE(DENS)

        BAS=0.0
        BASX=0.0
        BASY=0.0
        BASZ=0.0
        DO I=1,KONTA
         IP=LIP(I)
         BAS=BAS+MASAP(IP)
         BASX=BASX+RXPA(IP)*MASAP(IP)
         BASY=BASY+RYPA(IP)*MASAP(IP)
         BASZ=BASZ+RZPA(IP)*MASAP(IP)
        END DO

        CX=BASX/BAS
        CY=BASY/BAS
        CZ=BASZ/BAS

        RETURN
        END


**********************************************************************
       SUBROUTINE HALOFIND_PARTICLES(NL,NCLUS,MASA,RADIO,CLUSRX,CLUSRY,
     &      CLUSRZ,REALCLUS,CONCENTRA,ANGULARM,VMAXCLUS,IPLIP,VX,VY,VZ,
     &      VCMAX,MCMAX,RCMAX,M200C,M500C,M2500C,M200M,M500M,M2500M,
     &      MSUB,R200C,R500C,R2500C,R200M,R500M,R2500M,RSUB,DMPCLUS,
     &      LEVHAL,EIGENVAL,N_DM,RXPA,RYPA,RZPA,MASAP,U2DM,U3DM,U4DM,
     &      ORIPA,CONTRASTEC,OMEGAZ,UM,UV,LADO0,CLUSRXCM,CLUSRYCM,
     &      CLUSRZCM,MEAN_VR,INERTIA_TENSOR,NPATCH,PATCHCLUS,PROFILES,
     &      VELOCITY_DISPERSION,KINETIC_E,POTENTIAL_E,
     &      DO_COMPUTE_ENERGIES,PARTICLES_PER_HALO,
     &      INDCS_PARTICLES_PER_HALO,FLAG_WDM,ZETA,MIN_NUM_PART)
**********************************************************************
*      Refines halo identification with DM particles
**********************************************************************
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NL,NCLUS
       REAL*4 MASA(MAXNCLUS),RADIO(MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
       INTEGER REALCLUS(MAXNCLUS)
       REAL*4 CONCENTRA(NMAXNCLUS),ANGULARM(3,NMAXNCLUS)
       REAL*4 VMAXCLUS(NMAXNCLUS)
       INTEGER IPLIP(NMAXNCLUS)
       REAL*4 VX(NMAXNCLUS),VY(NMAXNCLUS),VZ(NMAXNCLUS)
       REAL*4 VCMAX(NMAXNCLUS),MCMAX(NMAXNCLUS),RCMAX(NMAXNCLUS)
       REAL*4 M200C(NMAXNCLUS),R200C(NMAXNCLUS)
       REAL*4 M500C(NMAXNCLUS),R500C(NMAXNCLUS)
       REAL*4 M2500C(NMAXNCLUS),R2500C(NMAXNCLUS)
       REAL*4 M200M(NMAXNCLUS),R200M(NMAXNCLUS)
       REAL*4 M500M(NMAXNCLUS),R500M(NMAXNCLUS)
       REAL*4 M2500M(NMAXNCLUS),R2500M(NMAXNCLUS)
       REAL*4 MSUB(MAXNCLUS),RSUB(MAXNCLUS)
       INTEGER DMPCLUS(NMAXNCLUS),LEVHAL(MAXNCLUS)
       REAL*4 EIGENVAL(3,NMAXNCLUS)
       INTEGER N_DM
       REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED)
       REAL*4 U2DM(PARTIRED),U3DM(PARTIRED),U4DM(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       INTEGER ORIPA(PARTIRED)
       REAL*4 CONTRASTEC,OMEGAZ,UM,UV,LADO0
       REAL*4 CLUSRXCM(MAXNCLUS),CLUSRYCM(MAXNCLUS),CLUSRZCM(MAXNCLUS)
       REAL*4 MEAN_VR(NMAXNCLUS),INERTIA_TENSOR(6,NMAXNCLUS)
       INTEGER NPATCH(0:NLEVELS),PATCHCLUS(MAXNCLUS)
       REAL*4 PROFILES(NBINS,2,NMAXNCLUS),VELOCITY_DISPERSION(NMAXNCLUS)
       REAL*4 KINETIC_E(NMAXNCLUS),POTENTIAL_E(NMAXNCLUS)
       INTEGER DO_COMPUTE_ENERGIES
       INTEGER PARTICLES_PER_HALO(PARTIRED)
       INTEGER INDCS_PARTICLES_PER_HALO(2,NMAXNCLUS),FLAG_WDM
       REAL*4 ZETA
       INTEGER MIN_NUM_PART

       REAL*4 PI,ACHE,T0,RE0,PI4ROD
       COMMON /DOS/ACHE,T0,RE0
       REAL*4 UNTERCIO,CGR,CGR2,ZI,RODO,ROI,REI
       COMMON /CONS/PI4ROD,REI,CGR,PI
       REAL*4 OMEGA0

       REAL*4 RETE,HTE,ROTE
       COMMON /BACK/ RETE,HTE,ROTE

       REAL*4 DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

*      Local variables
       INTEGER DIMEN,KONTA1,KONTA2,I,KK_ENTERO,NROT,II
       INTEGER KONTA,IR,J,CONTAERR,JJ,SALIDA,KONTA3,NSHELL_2,KONTA2PREV
       INTEGER IX,JY,KK1,KK2,FAC,ITER_SHRINK,COUNT_1,COUNT_2,JJCORE
       INTEGER FLAG200C,FLAG500C,FLAG2500C,FLAG200M,FLAG500M,FLAG2500M
       INTEGER FLAGVIR,EACH_PROF,IPATCH,IRR,MOST_BOUND_IDX,LOWP1,LOWP2
       INTEGER NCAPAS(NMAXNCLUS),MAX_NUM_PART,WELL_ALLOCATED,IDX_VIR
       REAL ESP,REF,REF_MIN,REF_MAX,DIS,VCM,MINOVERDENS
       REAL VVV2,CONCEN,RS,BAS,AADM,CMX,CMY,CMZ,VCMX,VCMY,CX,CY,CZ
       REAL VCMZ,MASA2,NORMA,BAS1,BAS2,VOL,DELTA2,RSHELL,RCLUS,BASVEC(3)
       REAL DENSA,DENSB,DENSC,VKK,AA,XP,YP,ZP,MP,t1,t2
       REAL INERTIA(3,3),BASEIGENVAL(3),RADII_ITER(7),BASVECCM(3)
       REAL BASVCM(3),EPOT,GCONS
       REAL DENSITOT(0:1000),RADIAL(0:1000),DENSR(0:1000)
       REAL LOGDERIV(0:1000)

       INTEGER,ALLOCATABLE::LIP(:),CONTADM(:)
       REAL,ALLOCATABLE::DISTA(:)

*      DOUBLE PRECISION VARIABLES
       REAL*8 MASADM,BASMAS,VR,BASX,BASY,BASZ,BAS8,BASVX,BASVY,BASVZ
       REAL*8 SIGMA_HALO,EKIN,INERTIA8(3,3)

       ! For writing DM particles
       INTEGER,ALLOCATABLE::PARTICLES_PROC(:,:),PROC_NPARTICLES(:)
       INTEGER,ALLOCATABLE::HALOES_PROC(:,:)
       INTEGER NUM_PROC,ID_PROC,IPART_PROC,OMP_GET_THREAD_NUM
       COMMON /PROCESADORES/ NUM_PROC

       PI=DACOS(-1.D0)
       GCONS=4.301E-9 ! in Msun^-1*Mpc*km^2*s^-2

       DIMEN=3   !DIMENSION DE LOS HALOS
       NCAPAS=0
       ESP=0.0
       REF=0.0
       KONTA1=0
       KONTA2=0

       MINOVERDENS=MIN(200.0,CONTRASTEC)

c       DO I=0,NL
c       WRITE(*,*)'Halos at level ', I,' =',
c     &            COUNT(LEVHAL(1:NCLUS).EQ.I),
c     &            COUNT(REALCLUS(1:NCLUS).NE.0.AND.
c     &                       LEVHAL(1:NCLUS).EQ.I)
c       END DO
c       WRITE(*,*)'=================================='

       MAX_NUM_PART=MAXVAL(DMPCLUS(1:NCLUS))
       WRITE(*,*) 'Max num. of part. in a halo=',
     &            MAX_NUM_PART
       KONTA1=100000000
       DO I=1,NCLUS
        IF (REALCLUS(I).NE.0) KONTA1=MIN(KONTA1,DMPCLUS(I))
       END DO
       WRITE(*,*) 'Min num. of part. in a halo=',KONTA1
       WRITE(*,*) 'NCLUS=', NCLUS

       NORMA=MAXVAL(MASAP)

       IF (FLAG_WDM.EQ.1) THEN
        ALLOCATE(PARTICLES_PROC(PARTIRED,NUM_PROC),
     &           HALOES_PROC(3,NCLUS),
     &           PROC_NPARTICLES(NUM_PROC))

        PROC_NPARTICLES(1:NUM_PROC)=0
       END IF

       MAX_NUM_PART=MIN(INT(MAX(1.05*MAX_NUM_PART,
     &                   MAX_NUM_PART*(CONTRASTEC/MINOVERDENS)**(1.0))),
     &                  PARTIRED)

!$OMP  PARALLEL DO SHARED(NCLUS,REALCLUS,PATCHCLUS,NPATCH,
!$OMP+           LEVHAL,RXPA,RYPA,RZPA,CLUSRX,CLUSRY,CLUSRZ,NL,MASAP,
!$OMP+           U2DM,U3DM,U4DM,VX,VY,VZ,ACHE,PI,RETE,ROTE,VCMAX,
!$OMP+           MCMAX,RCMAX,CONTRASTEC,OMEGAZ,CGR,UM,UV,DMPCLUS,
!$OMP+           CONCENTRA,ORIPA,ANGULARM,IPLIP,DIMEN,EIGENVAL,
!$OMP+           MIN_NUM_PART,RADIO,MASA,VMAXCLUS,N_DM,NORMA,R200M,
!$OMP+           R500M,R2500M,R200C,R500C,R2500C,M200M,M500M,M2500M,
!$OMP+           M200C,M500C,M2500C,RSUB,MSUB,MINOVERDENS,CLUSRXCM,
!$OMP+           CLUSRYCM,CLUSRZCM,DX,MEAN_VR,INERTIA_TENSOR,PROFILES,
!$OMP+           VELOCITY_DISPERSION,GCONS,KINETIC_E,POTENTIAL_E,
!$OMP+           DO_COMPUTE_ENERGIES,PARTICLES_PROC,HALOES_PROC,
!$OMP+           PROC_NPARTICLES,FLAG_WDM,ZETA),
!$OMP+   PRIVATE(I,INERTIA,REF_MIN,REF_MAX,KK_ENTERO,MASADM,KONTA,
!$OMP+           BASMAS,DIS,VCM,VVV2,VR,LIP,CONCEN,RS,KONTA2,BAS,IR,J,
!$OMP+           AADM,KK1,KK2,CONTADM,CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MASA2,
!$OMP+           DISTA,FAC,CONTAERR,JJ,DENSITOT,RADIAL,SALIDA,BAS1,BAS2,
!$OMP+           VOL,DELTA2,NCAPAS,RSHELL,KONTA3,NSHELL_2,KONTA1,
!$OMP+           DENSA,DENSB,DENSC,BASVEC,BASVECCM,VKK,AA,NROT,
!$OMP+           BASEIGENVAL,BASX,BASY,BASZ,XP,YP,ZP,MP,RCLUS,COUNT_1,
!$OMP+           COUNT_2,KONTA2PREV,FLAG200C,FLAG200M,FLAG500C,FLAG500M,
!$OMP+           FLAG2500C,FLAG2500M,FLAGVIR,EACH_PROF,DENSR,LOGDERIV,
!$OMP+           CX,CY,CZ,JJCORE,RADII_ITER,BASVCM,IPATCH,IRR,BASVX,
!$OMP+           BASVY,BASVZ,SIGMA_HALO,EKIN,EPOT,MOST_BOUND_IDX,
!$OMP+           ID_PROC,IPART_PROC,BAS8,INERTIA8,MAX_NUM_PART,
!$OMP+           WELL_ALLOCATED,IDX_VIR),
!$OMP+   SCHEDULE(DYNAMIC), DEFAULT(NONE)
*****************************
       DO I=1,NCLUS
****************************
       KK_ENTERO=REALCLUS(I)
       IF (KK_ENTERO.NE.0) THEN

        MAX_NUM_PART=INT(1.5*DMPCLUS(I))
        IF (MAX_NUM_PART.LT.4000) MAX_NUM_PART=4000
        IF (MAX_NUM_PART.GT.PARTIRED) MAX_NUM_PART=PARTIRED

        !write(*,*) max_num_part

*********************************************************************
*       RECENTERING
*********************************************************************
        ! Find level used for recentering
        IPATCH=PATCHCLUS(I)
        DO IRR=1,NL
         IF (SUM(NPATCH(0:IRR-1))+1.LE.IPATCH.AND.
     &       IPATCH.LE.SUM(NPATCH(0:IRR))) EXIT
        END DO

        CX=CLUSRX(I)
        CY=CLUSRY(I)
        CZ=CLUSRZ(I)
        BAS=RADIO(I)

        CALL RECENTER_DENSITY_PEAK_PARTICLES(CX,CY,CZ,BAS,RXPA,RYPA,
     &                                       RZPA,MASAP,N_DM,
     &                                       DX/2.0**IRR,MAX_NUM_PART)

        BAS=(CLUSRX(I)-CX)**2+(CLUSRY(I)-CY)**2+(CLUSRZ(I)-CZ)**2
        BAS=SQRT(BAS)
c        WRITE(*,*) 'Recentering shift', i, bas, bas/radio(i)
        CLUSRX(I)=CX
        CLUSRY(I)=CY
        CLUSRZ(I)=CZ

*********************************************************************
*       ALLOCATING PARTICLES
*********************************************************************
        WELL_ALLOCATED=0
        DO WHILE (WELL_ALLOCATED.EQ.0)
         IF (ALLOCATED(LIP)) DEALLOCATE(CONTADM,LIP,DISTA)
         ALLOCATE(LIP(MAX_NUM_PART),CONTADM(MAX_NUM_PART))
         ALLOCATE (DISTA(0:MAX_NUM_PART))
         WELL_ALLOCATED=1

         INERTIA8=0.D0
         REF_MIN=10.0e+10
         REF_MAX=-1.0
         MASADM=0.D0
         KONTA=0
         BASMAS=0.D0
         DIS=1.0E+10    !Distance (to the center) of the most central particle
         VCM=0.0
         VVV2=0.0    ! v propia de las particulas respecto al CM
         VR=0.D0      ! v radial
         CMX=0.0
         CMY=0.0
         CMZ=0.0
         VCMX=0.0
         VCMY=0.0
         VCMZ=0.0
         LIP=0
         CONCEN=0.0       !concentration NFW profile
         RS=0.0           !scale radius NFW profile
         KONTA2=0
         BASX=0.D0
         BASY=0.D0
         BASZ=0.D0

         CX=CLUSRX(I)
         CY=CLUSRY(I)
         CZ=CLUSRZ(I)
         RCLUS=RADIO(I)
         DELTA2=100.0*MINOVERDENS ! just to ensure it enters the loop
         DO WHILE (DELTA2.GT.0.99*MINOVERDENS)
          KONTA=0
          MASADM=0.D0

          BAS=-1.0
          DO J=1,N_DM
           XP=RXPA(J)
           YP=RYPA(J)
           ZP=RZPA(J)
           MP=MASAP(J)
           AADM=SQRT((XP-CX)**2+(YP-CY)**2+(ZP-CZ)**2)
           IF(AADM.LT.RCLUS) THEN
            KONTA=KONTA+1
            IF (KONTA.GT.MAX_NUM_PART) THEN
             WELL_ALLOCATED=0
             EXIT
            END IF
            LIP(KONTA)=J
            MASADM=MASADM+MP
            BAS=MAX(BAS,AADM)
           END IF
          END DO

          IF (WELL_ALLOCATED.EQ.0) EXIT

          DELTA2=MASADM/(ROTE*RETE**3*(4*PI/3)*BAS**3)

          IF (DELTA2.GT.0.9*MINOVERDENS) RCLUS=1.25*RCLUS
         END DO ! WHILE (DELTA2.GT.0.99*MINOVERDENS)

         IF (WELL_ALLOCATED.EQ.0) THEN
          !WRITE(*,*) 'WARNING: konta>max_num_part',KONTA,MAX_NUM_PART,I
          MAX_NUM_PART=MAX_NUM_PART*2
         END IF
        END DO ! WHILE (WELL_ALLOCATED.EQ.0)


        IF (MASADM.LE.0.D0.OR.KONTA.EQ.0) THEN
         REALCLUS(I)=0
         CYCLE
        END IF

        CONTADM(1:KONTA)=0     !en principio todas estas ligadas

        CONTAERR=KONTA
        DISTA=0.0
        KONTA2=0
        CALL REORDENAR(KONTA,CX,CY,CZ,RXPA,RYPA,RZPA,CONTADM,LIP,
     &                 DISTA,KONTA2,1,MAX_NUM_PART,IDX_VIR)
        REF_MAX=DISTA(KONTA2)
        REF_MIN=DISTA(1)

        RADIO(I)=RCLUS
        MASA(I)=MASADM*UM
        DMPCLUS(I)=KONTA

        CALL FIND_IDX_VIR(DISTA,MASAP,LIP,MAX_NUM_PART,ROTE,RETE,PI,
     &                    CONTRASTEC,IDX_VIR)
c        WRITE(*,*) '*',sqrt((cmx-cx)**2+(cmy-cy)**2+(cmz-cz)**2)
c        write(*,*) '**',vcmx,vcmy,vcmz,vcm

*********************************************************************
*       END RECENTERING AND COMPUTING VCM OF HALO I (SHRINKING SPHERE)
*********************************************************************

********************************************************************
*      UNBINDING:SCAPE VELOCITY
********************************************************************

        CALL CENTROMASAS_PART(IDX_VIR,CONTADM,LIP,U2DM,U3DM,U4DM,MASAP,
     &                        RXPA,RYPA,RZPA,CMX,CMY,CMZ,VCMX,VCMY,VCMZ,
     &                        MASA2,MAX_NUM_PART)
        VCM=SQRT(VCMX**2+VCMY**2+VCMZ**2)
        CLUSRXCM(I)=CMX
        CLUSRYCM(I)=CMY
        CLUSRZCM(I)=CMZ
        VX(I)=VCMX
        VY(I)=VCMY
        VZ(I)=VCMZ

        FAC=0
        DO WHILE (CONTAERR.GT.0.OR.FAC.LT.3)
         FAC=FAC+1
         KONTA2PREV=KONTA2
         CALL UNBINDING8(FAC,I,REF_MIN,REF_MAX,DISTA,U2DM,U3DM,U4DM,
     &                   MASAP,RXPA,RYPA,RZPA,RADIO,MASA,CLUSRXCM,
     &                   CLUSRYCM,CLUSRZCM,LIP,KONTA,CONTADM,VX,VY,VZ,
     &                   REALCLUS,KONTA2,MAX_NUM_PART,IDX_VIR)
         CALL REORDENAR(KONTA,CX,CY,CZ,RXPA,RYPA,RZPA,CONTADM,LIP,
     &                  DISTA,KONTA2,0,MAX_NUM_PART,IDX_VIR)
         REF_MAX=DISTA(KONTA2)
         REF_MIN=DISTA(1)
         CONTAERR=KONTA2PREV-KONTA2
        END DO

        count_1=konta-konta2
        count_2=konta2 !backup
c        write(*,*) 'Unbinding V_ESC',i,'. ',konta,'-->',konta2,
c     &             '. Pruned:',count_1,'. Iters:', FAC

********************************************************************
*      UNBINDING: PHASE SPACE
********************************************************************

        FAC=0
        CONTAERR=KONTA2
        DO WHILE (CONTAERR.GT.0.OR.FAC.LT.4)
         FAC=FAC+1
         KONTA2PREV=KONTA2
         CALL UNBINDING_SIGMA(FAC,I,REF_MIN,REF_MAX,U2DM,U3DM,U4DM,RXPA,
     &                        RYPA,RZPA,MASAP,RADIO,MASA,CLUSRXCM,
     &                        CLUSRYCM,CLUSRZCM,LIP,KONTA,CONTADM,VX,
     &                        VY,VZ,REALCLUS,KONTA2,MAX_NUM_PART,
     &                        IDX_VIR)
         CALL REORDENAR(KONTA,CX,CY,CZ,RXPA,RYPA,RZPA,CONTADM,LIP,
     &                  DISTA,KONTA2,0,MAX_NUM_PART,IDX_VIR)
         REF_MAX=DISTA(KONTA2)
         REF_MIN=DISTA(1)
         CONTAERR=KONTA2PREV-KONTA2
         !write(*,*) 'sigma unbinding: iter,unbound',fac,contaerr
        END DO

        count_2=count_2-konta2
c        write(*,*) 'Unbinding SIGMA',i,'. ',konta,'-->',konta2,
c     &             '. Pruned:',count_2,'. Iters:', FAC
c        write(*,*) '--'

********************************************************************
*      DISCARD POOR HALOES
********************************************************************
        IF (KONTA2.LT.MIN_NUM_PART) THEN
         REALCLUS(I)=0
        ELSE
**************************************************************
*      DENSITY PROFILE
**************************************************************
         REF_MIN=DISTA(1)
         REF_MAX=DISTA(KONTA2) ! DISTANCE TO THE FURTHERST BOUND PARTICLE
         JJ=0
         DENSITOT=0.0
         RADIAL=0.0
****************** stopping condition
         SALIDA=0
******************
         BAS8=0.D0 ! to store cummulative mass profile
         BAS1=-1.0 ! to store vcmax
         BAS2=0.0 ! to store vc(r)
         BASVX=0.D0 ! to compute CM velocity
         BASVY=0.D0
         BASVZ=0.D0
*        VCMAX=-1.0
         ! Initialise flags (whether each overdensity has been found)
         FLAG200C=0
         FLAG500C=0
         FLAG2500C=0
         FLAG200M=0
         FLAG500M=0
         FLAG2500M=0
         FLAGVIR=0

         FAC=MAX(100,INT(0.05*KONTA2))
         KONTA3=INT(REAL(KONTA2)/FAC)*FAC
         NSHELL_2=0
         JJCORE=INT(MAX(10.0,KONTA2/100.0))
         DO J=1,JJCORE
          JJ=LIP(J)
          BAS8=BAS8+(MASAP(JJ)/NORMA)
          BASVX=BASVX+(MASAP(JJ)/NORMA)*U2DM(JJ)
          BASVY=BASVY+(MASAP(JJ)/NORMA)*U3DM(JJ)
          BASVZ=BASVZ+(MASAP(JJ)/NORMA)*U4DM(JJ)
         END DO
         DO J=JJCORE+1,KONTA2      !!!!! DEJO 80 por ciento de BINS DE SEGURIDAD
          JJ=LIP(J)
          VOL=PI*(4.0/3.0)*(DISTA(J)*RETE)**3
          BAS8=BAS8+(MASAP(JJ)/NORMA)
          BASVX=BASVX+(MASAP(JJ)/NORMA)*U2DM(JJ)
          BASVY=BASVY+(MASAP(JJ)/NORMA)*U3DM(JJ)
          BASVZ=BASVZ+(MASAP(JJ)/NORMA)*U4DM(JJ)

          DELTA2=NORMA*BAS8/VOL/ROTE ! overdensity

          BAS2=BAS8/DISTA(J) ! vc(r)
          IF (BAS2.GT.BAS1) THEN
           VCMAX(I)=BAS2
           MCMAX(I)=BAS8
           RCMAX(I)=DISTA(J)
           BAS1=BAS2
          END IF

          IF (FLAG2500C.EQ.0) THEN
           IF (DELTA2.LE.2500.0/OMEGAZ) THEN
            M2500C(I)=DELTA2*VOL*ROTE*UM
            R2500C(I)=DISTA(J)
            FLAG2500C=1
            IF (J.EQ.JJCORE+1) THEN
             M2500C(I)=-1.0
             R2500C(I)=-1.0
            END IF
           END IF
          END IF

          IF (FLAG500C.EQ.0) THEN
           IF (DELTA2.LE.500.0/OMEGAZ) THEN
            M500C(I)=DELTA2*VOL*ROTE*UM
            R500C(I)=DISTA(J)
            FLAG500C=1
            IF (J.EQ.JJCORE+1) THEN
             M500C(I)=-1.0
             R500C(I)=-1.0
            END IF
           END IF
          END IF

          IF (FLAG200C.EQ.0) THEN
           IF (DELTA2.LE.200.0/OMEGAZ) THEN
            M200C(I)=DELTA2*VOL*ROTE*UM
            R200C(I)=DISTA(J)
            FLAG200C=1
            IF (J.EQ.JJCORE+1) THEN
             M200C(I)=-1.0
             R200C(I)=-1.0
            END IF
           END IF
          END IF

          IF (FLAG2500M.EQ.0) THEN
           IF (DELTA2.LE.2500.0) THEN
            M2500M(I)=DELTA2*VOL*ROTE*UM
            R2500M(I)=DISTA(J)
            FLAG2500M=1
            IF (J.EQ.JJCORE+1) THEN
             M2500M(I)=-1.0
             R2500M(I)=-1.0
            END IF
           END IF
          END IF

          IF (FLAG500M.EQ.0) THEN
           IF (DELTA2.LE.500.0) THEN
            M500M(I)=DELTA2*VOL*ROTE*UM
            R500M(I)=DISTA(J)
            FLAG500M=1
            IF (J.EQ.JJCORE+1) THEN
             M500M(I)=-1.0
             R500M(I)=-1.0
            END IF
           END IF
          END IF

          IF (FLAG200M.EQ.0) THEN
           IF (DELTA2.LE.200.0) THEN
            M200M(I)=DELTA2*VOL*ROTE*UM
            R200M(I)=DISTA(J)
            FLAG200M=1
            IF (J.EQ.JJCORE+1) THEN
             M200M(I)=-1.0
             R200M(I)=-1.0
            END IF
           END IF
          END IF

          ! STOPPING CONDITIONS
          ! 1. Below virial overdensity
          IF (FLAGVIR.EQ.0) THEN
           IF (DELTA2.LE.CONTRASTEC.AND.J.GT.INT(0.1*KONTA2)) THEN
            SALIDA=1 ! this means we've found the virial radius
                     ! but we don't exit straightaway because we still
                     ! want to find 200m
            FLAGVIR=1
            !MASA(I)=DELTA2*VOL*ROTE*UM
            !RADIO(I)=DISTA(J)
            BAS8=BAS8-(MASAP(JJ)/NORMA)
            BASVX=BASVX-(MASAP(JJ)/NORMA)*U2DM(JJ)
            BASVY=BASVY-(MASAP(JJ)/NORMA)*U3DM(JJ)
            BASVZ=BASVZ-(MASAP(JJ)/NORMA)*U4DM(JJ)

            MASA(I)=BAS8*NORMA*UM
            VOL=(NORMA*BAS8)/(CONTRASTEC*ROTE)
            RADIO(I)=(((3*VOL)/(4*PI))**(1.0/3.0))/RETE
            !WRITE(*,*) DISTA(J-1),RADIO(I),DISTA(J)

            VX(I)=BASVX/BAS8
            VY(I)=BASVY/BAS8
            VZ(I)=BASVZ/BAS8

            NCAPAS(I)=J-1
            !IF (J.EQ.KONTA2) THEN
            ! RSHELL=DISTA(J)
            !ELSE
            ! RSHELL=DISTA(J+1)
            !END IF
            RSHELL=RADIO(I)

           END IF
          END IF

          ! profile
          IF (MOD(J,FAC).EQ.0) THEN
           NSHELL_2=NSHELL_2+1
           DENSITOT(NSHELL_2)=NORMA*BAS8*UM
           RADIAL(NSHELL_2)=DISTA(J)
          END IF

          IF (SALIDA.EQ.1.AND.FLAG200M.EQ.1) EXIT
         END DO ! J=1,KONTA2

         !WRITE(*,*) RADIAL(1:NSHELL_2)
         !WRITE(*,*) DENSITOT(1:NSHELL_2)

         IF (MOD(J,FAC).NE.0) THEN
          NSHELL_2=NSHELL_2+1
          DENSITOT(NSHELL_2)=MASA(I)
          RADIAL(NSHELL_2)=RSHELL
         END IF

c         WRITE(*,*) 'HALO I,KONTA2,NSHELL_2,KK_ENTERO=',
c     &               I,KONTA2,NSHELL_2,KK_ENTERO
c         WRITE(*,*) CLUSRX(I),CLUSRY(I),CLUSRZ(I)
c         WRITE(*,*) R2500C(I),R500C(I),R200C(I),R2500M(I),R500M(I),
c     &              R200M(I),RADIO(I)
c         WRITE(*,*) M2500C(I),M500C(I),M200C(I),M2500M(I),M500M(I),
c     &              M200M(I),MASA(I)
c         WRITE(*,*) '---'

         BAS=VCMAX(I)*NORMA*CGR/RETE
         VCMAX(I)=SQRT(BAS)*UV
         MCMAX(I)=MCMAX(I)*NORMA*UM
         RCMAX(I)=RCMAX(I)   !*RETE

         IF (KK_ENTERO.EQ.-1.AND.SALIDA.NE.1) THEN
          SALIDA=1
          FLAGVIR=1
          !MASA(I)=DELTA2*VOL*ROTE*UM
          !RADIO(I)=DISTA(J)
          BAS8=BAS8!-(MASAP(JJ)/NORMA)
          BASVX=BASVX!-(MASAP(JJ)/NORMA)*U2DM(JJ)
          BASVY=BASVY!-(MASAP(JJ)/NORMA)*U3DM(JJ)
          BASVZ=BASVZ!-(MASAP(JJ)/NORMA)*U4DM(JJ)

          MASA(I)=BAS8*NORMA*UM
          VOL=(NORMA*BAS8)/(CONTRASTEC*ROTE)
          RADIO(I)=(((3*VOL)/(4*PI))**(1.0/3.0))/RETE
          !WRITE(*,*) DISTA(J-1),RADIO(I),DISTA(J)

          VX(I)=BASVX/BAS8
          VY(I)=BASVY/BAS8
          VZ(I)=BASVZ/BAS8

          NCAPAS(I)=J-1
          !IF (J.EQ.KONTA2) THEN
          ! RSHELL=DISTA(J)
          !ELSE
          ! RSHELL=DISTA(J+1)
          !END IF
          RSHELL=RADIO(I)

          IF (FLAG200M.EQ.0) THEN
           M200M(I)=MASA(I)
           R200M(I)=RADIO(I)*(CONTRASTEC/200.0)**(1.0/3.0)
          END IF

          WRITE(*,*) 'POSSIBLE PROBLEM WITH HALO',I,DELTA2,CMX,CMY,CMZ,
     &                NSHELL_2,RADIAL(NSHELL_2)

          IF (MOD(J,FAC).NE.0) THEN
           NSHELL_2=NSHELL_2+1
           DENSITOT(NSHELL_2)=MASA(I)
           RADIAL(NSHELL_2)=RSHELL
          END IF
         END IF

***********************************************************
*      GUARDAMOS LAS PARTICULAS LIGADAS  DEL HALO I
***********************************************************
         KONTA=NCAPAS(I)
         IF (KONTA.LT.MIN_NUM_PART) THEN
          REALCLUS(I)=0
          CYCLE
         END IF
         KONTA2=0
         BASMAS=0.D0
         DMPCLUS(I)=0

         VCMX=VX(I)
         VCMY=VY(I)
         VCMZ=VZ(I)

         BASX=0.D0
         BASY=0.D0
         BASZ=0.D0
         INERTIA8=0.D0
         SIGMA_HALO=0.D0
         EKIN=0.D0
         EPOT=0.0

         DIS=1000000.0
         VMAXCLUS(I)=-1.0
         EACH_PROF=NCAPAS(I)/NBINS
         DO J=1,KONTA
          IF (CONTADM(J).EQ.0) THEN
           JJ=LIP(J)
           BASVEC(1)=RXPA(JJ)-CX
           BASVEC(2)=RYPA(JJ)-CY
           BASVEC(3)=RZPA(JJ)-CZ
           AADM=SQRT(BASVEC(1)**2+BASVEC(2)**2+BASVEC(3)**2)
           IF (AADM.LE.RSHELL) THEN
            KONTA2=KONTA2+1
            BASMAS=BASMAS+(MASAP(JJ)/NORMA)

            IF (MOD(J,EACH_PROF).EQ.0) THEN
             IF (J/EACH_PROF.LE.NBINS) THEN
              PROFILES(J/EACH_PROF,1,I)=AADM
              PROFILES(J/EACH_PROF,2,I)=BASMAS*NORMA*UM
             END IF
            END IF

            BASVECCM(1)=RXPA(JJ)-CMX
            BASVECCM(2)=RYPA(JJ)-CMY
            BASVECCM(3)=RZPA(JJ)-CMZ

            BASVCM(1)=U2DM(JJ)-VCMX
            BASVCM(2)=U3DM(JJ)-VCMY
            BASVCM(3)=U4DM(JJ)-VCMZ

            VVV2=BASVCM(1)**2+BASVCM(2)**2+BASVCM(3)**2
            SIGMA_HALO=SIGMA_HALO+VVV2
            EKIN=EKIN+MASAP(JJ)*VVV2

**          ANGULAR MOMENTUM
            BASX=BASX+MASAP(JJ)*(BASVECCM(2)*BASVCM(3)
     &                          -BASVECCM(3)*BASVCM(2))
            BASY=BASY+MASAP(JJ)*(BASVECCM(3)*BASVCM(1)
     &                          -BASVECCM(1)*BASVCM(3))
            BASZ=BASZ+MASAP(JJ)*(BASVECCM(1)*BASVCM(2)
     &                          -BASVECCM(2)*BASVCM(1))

**          INERTIA TENSOR
            DO JY=1,3
            DO IX=1,3
              INERTIA8(IX,JY)=INERTIA8(IX,JY)
     &                         +MASAP(JJ)*BASVECCM(IX)*BASVECCM(JY)
            END DO
            END DO

**        VELOCITY OF THE FASTEST PARTICLE IN THE HALO
            IF (VVV2.GT.VMAXCLUS(I)) THEN
             VMAXCLUS(I)=VVV2
            END IF

**        CLOSEST PARTICLE TO THE CENTER OF THE HALO
C            IF (AADM.LT.DIS) THEN
C             DIS=AADM
C             IPLIP(I)=ORIPA(JJ)
C            END IF

            IF(AADM.NE.0.0) THEN
             AA=(BASVEC(1)/AADM)*BASVCM(1)+
     &          (BASVEC(2)/AADM)*BASVCM(2)+
     &          (BASVEC(3)/AADM)*BASVCM(3)
             VR=VR+AA*MASAP(JJ)
            END IF
           ELSE      !AADM.LT.RSHELL
            CONTADM(J)=1
           END IF    !AADM.LT.RSHELL
          END IF     !CONTADM
         END DO      !KONTA

*******************************************************
*      SAVING MASSES, RADII, PROFILES AND SHAPES...
*******************************************************
         DMPCLUS(I)=KONTA2
         IF (KONTA2.LT.MIN_NUM_PART) THEN
          REALCLUS(I)=0
          CYCLE
         END IF

         IF (FLAG_WDM.EQ.1) THEN
          ID_PROC=OMP_GET_THREAD_NUM()+1
          IPART_PROC=PROC_NPARTICLES(ID_PROC)
          HALOES_PROC(1,I)=ID_PROC
          HALOES_PROC(2,I)=IPART_PROC+1
          DO J=1,KONTA
           IF (CONTADM(J).EQ.0) THEN
            JJ=LIP(J)
            IPART_PROC=IPART_PROC+1
            PARTICLES_PROC(IPART_PROC,ID_PROC)=ORIPA(JJ)
           END IF
          END DO
          PROC_NPARTICLES(ID_PROC)=IPART_PROC
          HALOES_PROC(3,I)=IPART_PROC
         END IF

         IF (DO_COMPUTE_ENERGIES.EQ.1) THEN
          CALL COMPUTE_EPOT(KONTA,KONTA2,LIP,RXPA,RYPA,RZPA,MASAP,
     &                      CONTADM,EPOT,MOST_BOUND_IDX,MAX_NUM_PART)
          EPOT=EPOT*UM**2*GCONS ! Gravitational Energy in Msun * km^2 * s^-2
          EKIN=0.5*EKIN*UM*UV**2 ! Kinetic Energy in Msun * km^2 * s^-2
          KINETIC_E(I)=EKIN
          POTENTIAL_E(I)=EPOT*(1+ZETA)
          IPLIP(I)=ORIPA(MOST_BOUND_IDX)
c          WRITE(*,*) 'halo',i,'most bound particle',
c     & RXPA(MOST_BOUND_IDX),RYPA(MOST_BOUND_IDX),RZPA(MOST_BOUND_IDX),
c     & CLUSRX(I),CLUSRY(I),CLUSRZ(I)
         ELSE
          POTENTIAL_E(I)=0.0
          IPLIP(I)=ORIPA(LIP(1)) ! we take the centralmost particle as the most bound...
         END IF

         !MASA(I)=BASMAS*NORMA*UM
         !RADIO(I)=RSHELL
         ANGULARM(1,I)=BASX*UV / (BASMAS*NORMA)
         ANGULARM(2,I)=BASY*UV / (BASMAS*NORMA)
         ANGULARM(3,I)=BASZ*UV / (BASMAS*NORMA)
         MEAN_VR(I)=VR/(BASMAS*NORMA)*UV

         ! to simple precision
         INERTIA(1:3,1:3)=INERTIA8(1:3,1:3)/(BASMAS*NORMA)
         INERTIA(1,2)=INERTIA(2,1)
         INERTIA(1,3)=INERTIA(3,1)
         INERTIA(2,3)=INERTIA(3,2)
         INERTIA_TENSOR(1,I)=INERTIA(1,1)
         INERTIA_TENSOR(2,I)=INERTIA(1,2)
         INERTIA_TENSOR(3,I)=INERTIA(1,3)
         INERTIA_TENSOR(4,I)=INERTIA(2,2)
         INERTIA_TENSOR(5,I)=INERTIA(2,3)
         INERTIA_TENSOR(6,I)=INERTIA(3,3)

         VELOCITY_DISPERSION(I)=SQRT(SIGMA_HALO/FLOAT(KONTA2))*UV

         VMAXCLUS(I)=SQRT(VMAXCLUS(I))*UV

         BASEIGENVAL(1:3)=0.0
         IF (DMPCLUS(I).GE.MIN_NUM_PART) THEN
          CALL JACOBI(INERTIA,DIMEN,BASEIGENVAL,NROT)
          CALL SORT(BASEIGENVAL,DIMEN,DIMEN)
         END IF

         DO II=1,DIMEN
          EIGENVAL(II,I)=SQRT(5.0*BASEIGENVAL(II))
         END DO

C         WRITE(*,*)
C         WRITE(*,*)I,CLUSRX(I),CLUSRY(I),CLUSRZ(I),
C     &         MASA(I),RADIO(I),DMPCLUS(I),
C     &         REALCLUS(I), LEVHAL(I),(INERTIA_TENSOR(IX,I),IX=1,6),
C     &         EIGENVAL(1,I),EIGENVAL(2,I),EIGENVAL(3,I),MEAN_VR(I)*UV,
C     &         CONCENTRA(I),(ANGULARM(J,I),J=1,3),
C     &         VCMAX(I),MCMAX(I),RCMAX(I),
C     &         R200M(I),M200M(I),R200C(I),M200C(I),
C     &         R500M(I),M500M(I),R500C(I),M500C(I),
C     &         R2500M(I),M2500M(I),R2500C(I),M2500C(I),
C     &         VX(I)*UV,VY(I)*UV,VZ(I)*UV

******************************************
************ FIN HALO I ******************
******************************************
        END IF ! (ELSE of KONTA2.LT.MIN_NUM_PART) (poor haloes after unbinding)
       END IF ! (realclus(i).ne.0)

*****************
       END DO   !I=IP,IP2
****************

       IF (FLAG_WDM.EQ.1) THEN
        J=0

        DO I=1,NCLUS
         IF (REALCLUS(I).EQ.0) THEN
          INDCS_PARTICLES_PER_HALO(1,I)=-1
          INDCS_PARTICLES_PER_HALO(2,I)=-1
          CYCLE
         END IF
         INDCS_PARTICLES_PER_HALO(1,I)=J+1

         ID_PROC=HALOES_PROC(1,I)
         LOWP1=HALOES_PROC(2,I)
         LOWP2=HALOES_PROC(3,I)

         DO IPART_PROC=LOWP1,LOWP2
          J=J+1
          PARTICLES_PER_HALO(J)=PARTICLES_PROC(IPART_PROC,ID_PROC)
         END DO

         INDCS_PARTICLES_PER_HALO(2,I)=J
        END DO

        DEALLOCATE(HALOES_PROC, PARTICLES_PROC, PROC_NPARTICLES)
       END IF

       RETURN
       END

**********************************************************************
       SUBROUTINE FIND_IDX_VIR(DISTA,MASAP,LIP,MAX_NUM_PART,ROTE,RETE,
     &                         PI,CONTRASTEC,IDX_VIR)
**********************************************************************
*      Estimates the virial radius (before unbounding)
**********************************************************************
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       REAL DISTA(0:MAX_NUM_PART)
       REAL MASAP(PARTIRED)
       INTEGER LIP(MAX_NUM_PART)
       REAL ROTE,RETE,PI,CONTRASTEC
       INTEGER MAX_NUM_PART,IDX_VIR

       INTEGER J,JJCORE
       REAL DELTA
       REAL*8 MENC

       MENC=0.0
       JJCORE=MIN(IDX_VIR/10, 50)
       DO J=1,JJCORE
         MENC=MENC+DBLE(MASAP(LIP(J)))
       END DO

       DO J=JJCORE+1,IDX_VIR
         MENC=MENC+DBLE(MASAP(LIP(J)))
         DELTA=MENC/(ROTE*RETE**3*(4*PI/3)*DISTA(J)**3)
         IF (DELTA.LT.CONTRASTEC) EXIT
       END DO

       IF (J.EQ.IDX_VIR+1) THEN
        WRITE(*,*) 'NOT PRE-FOUND VIRIAL RADIUS!'
        STOP
       END IF

       IDX_VIR=J

       RETURN
       END

**********************************************************************
       SUBROUTINE COMPUTE_EPOT(KONTA,KONTA2,LIP,RXPA,RYPA,RZPA,MASAP,
     &                         CONTADM,EPOT,MOST_BOUND_IDX,MAX_NUM_PART)
**********************************************************************
*      Computes the gravitational potential energy of a halo, taking
*       into account only bound particles, by direct summation (small
*       halos) or does an estimation by sampling
**********************************************************************

       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER KONTA,KONTA2
       INTEGER LIP(MAX_NUM_PART)
       REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       INTEGER CONTADM(MAX_NUM_PART)
       REAL EPOT
       INTEGER MOST_BOUND_IDX
       INTEGER MAX_NUM_PART

       INTEGER J,JJ,K,KK,NSAMPLE!,NPAIRS
       INTEGER FLAG,JJJ,KKK
       REAL X1,X2,Y1,Y2,Z1,Z2,M1,M2,U,LARGEST_ENERGY

*      DOUBLE PRECISION / LONG LOCAL VARIABLES
       REAL*8 BAS,EPOT8
       REAL*8,ALLOCATABLE::EPOT_PART(:)
       INTEGER(8) NPAIRS,NPAIRS_TOTAL,NPAIRS_EXPECT,NBAS64 !this has to be 8-byte because of possible OVERFLOWS
       !REAL NPAIRS_TOTAL,NPAIRS_EXPECT

       EPOT8=0.D0


       IF (KONTA2.LT.1000) THEN ! direct summation
        ALLOCATE(EPOT_PART(KONTA))

        DO J=1,KONTA
         EPOT_PART(J)=0.D0
        END DO

        DO J=1,KONTA
         JJ=LIP(J)
         IF (CONTADM(J).EQ.1) CYCLE
         X1=RXPA(JJ)
         Y1=RYPA(JJ)
         Z1=RZPA(JJ)
         M1=MASAP(JJ)
         DO K=J+1,KONTA
          KK=LIP(K)
          IF (CONTADM(K).EQ.1) CYCLE
          X2=RXPA(KK)
          Y2=RYPA(KK)
          Z2=RZPA(KK)
          M2=MASAP(KK)
          BAS=-M1*M2/SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
          EPOT_PART(J)=EPOT_PART(J)+BAS
          EPOT_PART(K)=EPOT_PART(K)+BAS
          EPOT8=EPOT8+BAS
         END DO
        END DO

        DO J=1,KONTA
         IF (CONTADM(J).EQ.1) THEN
          EPOT_PART(J)=100000.D0
         ELSE
          EPOT_PART(J)=EPOT_PART(J)/MASAP(LIP(J)) ! gravitational energy per unit mass
         END IF
        END DO
        MOST_BOUND_IDX=LIP(MINLOC(EPOT_PART,1))

        DEALLOCATE(EPOT_PART)
       ELSE ! sampling
        NSAMPLE=MAX(1000,INT(0.01*KONTA2))
        NBAS64=NSAMPLE
        NPAIRS_EXPECT=NBAS64*(NBAS64-1)/2
        NPAIRS=0

        CALL RANDOM_SEED() ! set the seed for random numbers

        LARGEST_ENERGY=1000000.0
        DO JJJ=1,NSAMPLE
         ! Generate a random, valid particle
         FLAG=0
         DO WHILE (FLAG.EQ.0)
          CALL RANDOM_NUMBER(U)
          J=1+FLOOR(U*KONTA)
          IF (CONTADM(J).EQ.0) FLAG=1
         END DO
         JJ=LIP(J)
         X1=RXPA(JJ)
         Y1=RYPA(JJ)
         Z1=RZPA(JJ)
         M1=MASAP(JJ)

         BAS=0.D0

         DO KKK=JJJ+1,NSAMPLE
          FLAG=0
          DO WHILE (FLAG.EQ.0)
           CALL RANDOM_NUMBER(U)
           K=1+FLOOR(U*KONTA)
           IF (CONTADM(K).EQ.0) THEN
            FLAG=1
            IF (K.EQ.J) FLAG=0
           END IF
          END DO
          KK=LIP(K)
          X2=RXPA(KK)
          Y2=RYPA(KK)
          Z2=RZPA(KK)
          M2=MASAP(KK)

          NPAIRS=NPAIRS+1
          BAS=BAS-M1*M2/SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
         END DO

         EPOT8=EPOT8+BAS

         BAS=BAS/M1 ! gravitational energy per unit mass
         IF (BAS.LT.LARGEST_ENERGY) THEN
          LARGEST_ENERGY=BAS
          MOST_BOUND_IDX=JJ
         END IF
        END DO
        NBAS64=KONTA
        NPAIRS_TOTAL=NBAS64*(NBAS64-1)/2
        ! here NPAIRS_TOTAL is real, to avoid overflows with 4-byte integers
        !TEMP=EPOT
        EPOT8=EPOT8*FLOAT(NPAIRS_TOTAL)/FLOAT(NPAIRS)
       END IF !(KONTA2.LT.1000)

       EPOT=EPOT8

       RETURN
       END


**********************************************************************
       SUBROUTINE REORDENAR(KONTA,CX,CY,CZ,RXPA,RYPA,RZPA,CONTADM,LIP,
     &                      DISTA,KONTA2,DO_SORT,MAX_NUM_PART,IDX_VIR)
**********************************************************************
*      Sorts the particles with increasing distance to the center of
*      the halo. Only particles with CONTADM=0 are sorted (the others
*      are already pruned particles and therefore they are ignored)
*      Do_sort=1: particles are sorted by distance
*      Do_sort=0: particles are assumed to be sorted; just prune unbound
**********************************************************************

       IMPLICIT NONE

       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER I,J,K,KONTA,KONTA2,JJ,DO_SORT,MAX_NUM_PART,IDX_VIR

*      ---HALOS Y SUBHALOS---
       REAL*4 CX,CY,CZ
       INTEGER RELOC_VIR

*      ---PARTICULAS E ITERACIONES---
       INTEGER LIP(MAX_NUM_PART),CONTADM(MAX_NUM_PART)

       REAL*4 RXPA(PARTIRED)
       REAL*4 RYPA(PARTIRED)
       REAL*4 RZPA(PARTIRED)

c       INTEGER CONTADM(PARTI)

       REAL*4 DISTA(0:MAX_NUM_PART)

       INTEGER INDICE(KONTA)
       REAL*4 DISTA2(0:KONTA)
       INTEGER QUIEN(KONTA)

       REAL*4 AADMX,AADMY,AADMZ,AADM

       IF (DO_SORT.EQ.1) THEN ! we have to sort the particles
        QUIEN=0
        DISTA=1.E10
        INDICE=0
        DISTA2=0.0

        KONTA2=0
        DO J=1,KONTA
         IF (CONTADM(J).EQ.0) THEN
          KONTA2=KONTA2+1
          JJ=LIP(J)
          AADMX=RXPA(JJ)-CX
          AADMY=RYPA(JJ)-CY
          AADMZ=RZPA(JJ)-CZ
          AADM=SQRT(AADMX**2+AADMY**2+AADMZ**2)

          DISTA(KONTA2)=AADM
          QUIEN(KONTA2)=JJ
         END IF
        END DO

        CALL INDEXX(KONTA2,DISTA(1:KONTA2),INDICE(1:KONTA2))

        DO J=1,KONTA2
         DISTA2(J)=DISTA(J)
        END DO

        DISTA=1.E10
        LIP=0

        DO J=1,KONTA2
         JJ=INDICE(J)
         DISTA(J)=DISTA2(JJ)
         LIP(J)=QUIEN(JJ)
        END DO

        IDX_VIR=KONTA2 ! Yet to be determined
       ELSE ! they are already sorted (since the center does not change)
        KONTA2=0
        RELOC_VIR=0
        DO J=1,KONTA
         IF (CONTADM(J).EQ.0) THEN
          KONTA2=KONTA2+1
          DISTA(KONTA2)=DISTA(J)
          LIP(KONTA2)=LIP(J)

          IF (RELOC_VIR.EQ.0) THEN
           IF (J.GE.IDX_VIR) THEN
            RELOC_VIR=1
            IDX_VIR=KONTA2
           END IF
          END IF

         END IF
        END DO

        IF (RELOC_VIR.EQ.0) IDX_VIR=KONTA2
       END IF

       CONTADM=1
       CONTADM(1:KONTA2)=0

       RETURN
       END

***********************************************************
       SUBROUTINE CENTROMASAS_PART(N,CONTADM,LIP,
     &            U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &            CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MASA,MAX_NUM_PART)
***********************************************************
*      Computes the center of mass of the particles
*      within the halo
***********************************************************

       IMPLICIT NONE

       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER I,J,N,MAX_NUM_PART

       REAL*4 U2DM(PARTIRED)
       REAL*4 U3DM(PARTIRED)
       REAL*4 U4DM(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       REAL*4 RXPA(PARTIRED)
       REAL*4 RYPA(PARTIRED)
       REAL*4 RZPA(PARTIRED)

       INTEGER LIP(MAX_NUM_PART),CONTADM(MAX_NUM_PART)

       REAL*4 CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MASA

*      ---- DOUBLE PRECISION -----------------------
       REAL*8 CMX8,CMY8,CMZ8,VCMX8,VCMY8,VCMZ8,MASA8
       REAL*8 BAS!,NORMA
*      ---------------------------------------------


       CMX8=0.D0
       CMY8=0.D0
       CMZ8=0.D0

       VCMX8=0.D0
       VCMY8=0.D0
       VCMZ8=0.D0

       MASA8=0.D0

       !NORMA=DBLE(MAXVAL(MASAP))  ! NORMALIZACION MASA

       DO I=1,N
        IF (CONTADM(I).EQ.0) THEN
         J=LIP(I)

         BAS=DBLE(MASAP(J))!/NORMA

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

       MASA8=MASA8!*NORMA

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

***********************************************************
       SUBROUTINE UNBINDING8(FAC,I,REF_MIN,REF_MAX,DISTA,
     &           U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &           RADIO,MASA,CLUSRXCM,CLUSRYCM,CLUSRZCM,
     &           LIP,KONTA,CONTADM,VX,VY,VZ,REALCLUS,KONTA2,
     &           MAX_NUM_PART,IDX_VIR)
***********************************************************
*      Finds and discards the unbound particles (those
*      with speed larger than the scape velocity).
*      Potential is computed in double precision.
***********************************************************

       IMPLICIT NONE

       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER I,J,K,IX,IMAX,JJ,FAC,MAX_NUM_PART,IDX_VIR

       REAL*4 REI,CGR,PI,PI4ROD
       COMMON /CONS/PI4ROD,REI,CGR,PI

       REAL*4 RETE,HTE,ROTE
       COMMON /BACK/ RETE,HTE,ROTE

       REAL*4 REF_MIN,REF_MAX

*      ---HALOS Y SUBHALOS---
       REAL*4 RADIO(MAXNCLUS)
       REAL*4 MASA(MAXNCLUS)
       REAL*4 CLUSRXCM(MAXNCLUS)
       REAL*4 CLUSRYCM(MAXNCLUS)
       REAL*4 CLUSRZCM(MAXNCLUS)
       REAL*4 VX(NMAXNCLUS)
       REAL*4 VY(NMAXNCLUS)
       REAL*4 VZ(NMAXNCLUS)
       INTEGER REALCLUS(MAXNCLUS)

*      ---PARTICULAS E ITERACIONES---
       INTEGER LIP(MAX_NUM_PART),CONTADM(MAX_NUM_PART)
       INTEGER NPART(0:NLEVELS)
       REAL*4 U2DM(PARTIRED)
       REAL*4 U3DM(PARTIRED)
       REAL*4 U4DM(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       REAL*4 RXPA(PARTIRED)
       REAL*4 RYPA(PARTIRED)
       REAL*4 RZPA(PARTIRED)

       INTEGER KONTA,KONTA3,KONTA2
c       INTEGER CONTADM(PARTI)

       REAL*4 DISTA(0:MAX_NUM_PART)

       REAL*4 VVV2,VESC2,AADMX(3),AADM,DR, AA, BB, CC
       REAL*4 BAS
       REAL*4 CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MMM
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

       IF (KONTA2.GT.0) THEN
*      Max mass
        NORMA=DBLE(MAXVAL(MASAP))
        MASA8=DBLE(MASAP(1))/NORMA

        !POT(1)=MASA8/DBLE(DISTA(1))
        !WRITE(*,*) 'IN UNBINDING, KONTA2=',KONTA2
        DO J=1,KONTA2
         IF (DISTA(J).GT.0.01*REF_MAX) EXIT
        END DO
        JJ=J
        MASA8=0.D0
        DO J=1,JJ
         MASA8=MASA8+DBLE(MASAP(LIP(J)))/NORMA
        END DO
        DO J=1,JJ
         POT(J)=MASA8/DISTA(JJ)
        END DO

        !WRITE(*,*) 'KONTA2,JJ=',konta2,jj

        DO J=JJ+1,KONTA2
          MASA8=MASA8+DBLE(MASAP(LIP(J)))/NORMA
          IF (DISTA(J).NE.DISTA(J-1)) THEN
           BAS8=DISTA(J)-DISTA(J-1)
          ELSE
           BAS8=0.D0
          END IF
          POT(J)=POT(J-1)+MASA8*BAS8/(DBLE(DISTA(J))**2)
        END DO

        POT1=POT(KONTA2) + MASA8/REF_MAX
        !POT1 is the constant to be subtracted to the computed potential
        !so that the potential origin is located at infinity

        AA8=NORMA*DBLE(CGR/RETE)

        BB=2.0
        IF (FAC.EQ.1) BB=8.0
        IF (FAC.EQ.2) BB=4.0

        BB=BB**2 !(we compare the squared velocities)

*      Find particles able to escape the potential well
        DO J=1,KONTA2

         POTOK=(POT(J)-POT1)*AA8
         VESC2=2.0*ABS(POTOK)

         VVV2=(U2DM(LIP(J))-VX(I))**2
     &       +(U3DM(LIP(J))-VY(I))**2
     &       +(U4DM(LIP(J))-VZ(I))**2

         IF (VVV2.GT.BB*VESC2)  CONTADM(J)=1
        END DO

*      NEW CENTER OF MASS AND ITS VELOCITY
        CALL CENTROMASAS_PART(IDX_VIR,CONTADM,LIP,
     &           U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &           CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MMM,MAX_NUM_PART)

       END IF

       KONTA3=COUNT(CONTADM(1:KONTA2).EQ.0)

       IF (KONTA3.EQ.0) THEN

        CLUSRXCM(I)=0.0
        CLUSRYCM(I)=0.0
        CLUSRZCM(I)=0.0

        VX(I)=0.0
        VY(I)=0.0
        VZ(I)=0.0

        MASA(I)=0.0
        RADIO(I)=0.0

        REALCLUS(I)=0

       ELSE

        CLUSRXCM(I)=CMX
        CLUSRYCM(I)=CMY
        CLUSRZCM(I)=CMZ

        VX(I)=VCMX
        VY(I)=VCMY
        VZ(I)=VCMZ

        MASA(I)=MMM*9.1717e+18

       END IF

       RETURN
       END


***********************************************************
       SUBROUTINE UNBINDING_SIGMA(FAC,I,REF_MIN,REF_MAX,U2DM,U3DM,U4DM,
     &                            RXPA,RYPA,RZPA,MASAP,RADIO,MASA,
     &                            CLUSRXCM,CLUSRYCM,CLUSRZCM,LIP,KONTA,
     &                            CONTADM,VX,VY,VZ,REALCLUS,KONTA2,
     &                            MAX_NUM_PART,IDX_VIR)
***********************************************************
*      Finds and discards the unbound particles (those
*      with speed larger than the scape velocity).
*      Potential is computed in double precision.
***********************************************************

       IMPLICIT NONE

       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER FAC
       INTEGER I
       REAL REF_MIN,REF_MAX
       REAL*4 U2DM(PARTIRED),U3DM(PARTIRED),U4DM(PARTIRED)
       REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       REAL*4 RADIO(MAXNCLUS),MASA(MAXNCLUS)
       REAL*4 CLUSRXCM(MAXNCLUS),CLUSRYCM(MAXNCLUS),
     &        CLUSRZCM(MAXNCLUS)
       INTEGER LIP(MAX_NUM_PART),CONTADM(MAX_NUM_PART)
       INTEGER KONTA
       REAL*4 VX(NMAXNCLUS),VY(NMAXNCLUS),VZ(NMAXNCLUS)
       INTEGER REALCLUS(MAXNCLUS),KONTA2,MAX_NUM_PART,IDX_VIR

       REAL CMX,CMY,CMZ,VXCM,VYCM,VZCM,BAS,BB,AADM,AADMX,AADMY
       REAL AADMZ,MMM
       INTEGER J,JJ,KONTA3
       REAL,ALLOCATABLE::DESV2(:)

*      DOUBLE PRECISION VARIABLES
       REAL*8 SIGMA2

       BB=MAX(6.0-1.0*(FAC-1), 3.0)
       BB=BB**2 ! This is because we compare velocities squared

       CMX=CLUSRXCM(I)
       CMY=CLUSRYCM(I)
       CMZ=CLUSRZCM(I)
       VXCM=VX(I)
       VYCM=VY(I)
       VZCM=VZ(I)

       ALLOCATE(DESV2(1:KONTA2))

       IF (KONTA2.GT.0) THEN
        SIGMA2=0.D0
        DO J=1,KONTA2
         JJ=LIP(J)
         BAS=(U2DM(JJ)-VXCM)**2+(U3DM(JJ)-VYCM)**2+(U4DM(JJ)-VZCM)**2
         DESV2(J)=BAS
         SIGMA2=SIGMA2+BAS
        END DO

        IF (KONTA2.GT.1) SIGMA2=SIGMA2/(KONTA2-1)

*       Find particles with too large relative velocity
        DO J=1,KONTA2
         IF (DESV2(J).GT.BB*SIGMA2) CONTADM(J)=1
        END DO

*       NEW CENTER OF MASS AND ITS VELOCITY
        CALL CENTROMASAS_PART(IDX_VIR,CONTADM,LIP,
     &           U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &           CMX,CMY,CMZ,VXCM,VYCM,VZCM,MMM,MAX_NUM_PART)

       END IF

       KONTA3=COUNT(CONTADM(1:KONTA2).EQ.0)

       IF (KONTA3.EQ.0) THEN

        CLUSRXCM(I)=0.0
        CLUSRYCM(I)=0.0
        CLUSRZCM(I)=0.0

        VX(I)=0.0
        VY(I)=0.0
        VZ(I)=0.0

        MASA(I)=0.0
        RADIO(I)=0.0

        REALCLUS(I)=0

       ELSE

        CLUSRXCM(I)=CMX
        CLUSRYCM(I)=CMY
        CLUSRZCM(I)=CMZ

        VX(I)=VXCM
        VY(I)=VYCM
        VZ(I)=VZCM

        MASA(I)=MMM*9.1717e+18

       END IF

       DEALLOCATE(DESV2)


       RETURN
       END
