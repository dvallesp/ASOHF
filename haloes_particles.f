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
     &                              NUMPARTBAS,DMPCLUS,FACRAD)
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
        INTEGER N_DM,NUMPARTBAS
        INTEGER DMPCLUS(NMAXNCLUS)
        REAL FACRAD

        INTEGER I,BASINT,KONTA2
        REAL CMX,CMY,CMZ,RR

        KONTA2=0

!$OMP PARALLEL DO SHARED(NCLUS,CLUSRX,CLUSRY,CLUSRZ,RADIO,N_DM,
!$OMP+                   NUMPARTBAS,REALCLUS,RXPA,RYPA,RZPA,FACRAD,
!$OMP+                   DMPCLUS),
!$OMP+            PRIVATE(I,CMX,CMY,CMZ,RR,BASINT),
!$OMP+            REDUCTION(+:KONTA2)
!$OMP+            DEFAULT(NONE)
*****************************
        DO I=1,NCLUS
****************************
         CMX=CLUSRX(I)
         CMY=CLUSRY(I)
         CMZ=CLUSRZ(I)
         RR=FACRAD*RADIO(I)

         CALL COUNT_PARTICLES_HALO(RXPA,RYPA,RZPA,N_DM,CMX,CMY,CMZ,RR,
     &                             BASINT)

         IF (BASINT.LT.NUMPARTBAS) THEN
          REALCLUS(I)=0
          KONTA2=KONTA2+1
          !WRITE(*,*) I,BASINT,RR
         END IF

         DMPCLUS(I)=BASINT

*****************************
        END DO        !I
*****************************

        WRITE(*,*)'CHECKING POOR HALOS----->', KONTA2

        RETURN
        END


**********************************************************************
       SUBROUTINE HALOFIND_PARTICLES(NL,NCLUS,MASA,RADIO,CLUSRX,CLUSRY,
     &      CLUSRZ,REALCLUS,CONCENTRA,ANGULARM,VMAXCLUS,VCM2,IPLIP,
     &      VX,VCMAX,MCMAX,RCMAX,M200,R200,DMPCLUS,LEVHAL,EIGENVAL,
     &      N_DM,RXPA,RYPA,RZPA,MASAP,U2DM,U3DM,U4DM,ORIPA2,CONTRASTEC,
     &      OMEGAZ,UM,UV,F2)
**********************************************************************
*      Refines halo identification with DM particles
**********************************************************************
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NL,NCLUS
       REAL*4 MASA(MAXNCLUS),RADIO(MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
       INTEGER REALCLUS(MAXNCLUS)
       REAL*4 CONCENTRA(NMAXNCLUS),ANGULARM(NMAXNCLUS)
       REAL*4 VMAXCLUS(NMAXNCLUS),VCM2(NMAXNCLUS)
       INTEGER IPLIP(NMAXNCLUS)
       REAL*4 VX(NMAXNCLUS),VY(NMAXNCLUS),VZ(NMAXNCLUS)
       REAL*4 VCMAX(NMAXNCLUS),MCMAX(NMAXNCLUS),RCMAX(NMAXNCLUS)
       REAL*4 M200(NMAXNCLUS),R200(NMAXNCLUS)
       INTEGER DMPCLUS(NMAXNCLUS),LEVHAL(MAXNCLUS)
       REAL*4 EIGENVAL(3,NMAXNCLUS)
       INTEGER N_DM
       REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED)
       REAL*4 U2DM(PARTIRED),U3DM(PARTIRED),U4DM(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       INTEGER ORIPA2(PARTIRED)
       REAL*4 CONTRASTEC,OMEGAZ,UM,UV,F2

       REAL*4 PI,ACHE,T0,RE0,PI4ROD
       COMMON /DOS/ACHE,T0,RE0
       REAL*4 UNTERCIO,CGR,CGR2,ZI,RODO,ROI,REI,LADO,LADO0
       COMMON /CONS/PI4ROD,REI,CGR,PI
       REAL*4 OMEGA0

       REAL*4 RETE,HTE,ROTE
       COMMON /BACK/ RETE,HTE,ROTE

*      Local variables
       INTEGER DIMEN,KONTA1,KONTA2,I,PABAS,NUMPARTBAS,KK_ENTERO,NROT,II
       INTEGER KONTA,IR,J,CONTAERR,JJ,SALIDA,FLAG_200,KONTA3,NSHELL_2
       INTEGER IX,JY,KK1,KK2,FAC,ITER_SHRINK,IS_SUB,COUNT_1,COUNT_2
       INTEGER KONTA2PREV
       INTEGER NCAPAS(NMAXNCLUS)
       INTEGER LIP(PARTIRED),CONTADM(PARTIRED)
       REAL ESP,REF,REF_MIN,REF_MAX,MASADM,BASMAS,DIS,VCM
       REAL VVV2,VR,CONCEN,RS,BAS,AADM,CMX,CMY,CMZ,VCMX,VCMY
       REAL VCMZ,MASA2,NORMA,BAS1,BAS2,VOL,DELTA2,RSHELL,RCLUS
       REAL DENSA,DENSB,DENSC,VKK,AA,BASX,BASY,BASZ,XP,YP,ZP,MP
       REAL INERTIA(3,3),BASEIGENVAL(3),AADMX(3),RADII_ITER(7)
       REAL DENSITOT(0:PARTIRED/10),RADIAL(0:PARTIRED/10)
       REAL DISTA(0:PARTIRED)

       PI=DACOS(-1.D0)

       DIMEN=3   !DIMENSION DE LOS HALOS
       NCAPAS=0
       ESP=0.0
       REF=0.0
       KONTA1=0
       KONTA2=0

       DO I=0,NL
       WRITE(*,*)'Halos at level ', I,' =',
     &            COUNT(LEVHAL(1:NCLUS).EQ.I),
     &            COUNT(REALCLUS(1:NCLUS).NE.0.AND.
     &                       LEVHAL(1:NCLUS).EQ.I)
       END DO
       WRITE(*,*)'=================================='

       WRITE(*,*) 'Max num. of part. in a halo=',
     &            MAXVAL(DMPCLUS(1:NCLUS))
       KONTA1=1000000
       DO I=1,NCLUS
        IF (REALCLUS(I).NE.0) KONTA1=MIN(KONTA1,DMPCLUS(I))
       END DO
       WRITE(*,*) 'Min num. of part. in a halo=',
     &            KONTA1
       WRITE(*,*) 'NCLUS=', NCLUS

       PABAS=PARTIRED_PLOT
       NUMPARTBAS=NUMPART
       NORMA=MAXVAL(MASAP)

!$OMP  PARALLEL DO SHARED(NCLUS,REALCLUS,
!$OMP+           LEVHAL,RXPA,RYPA,RZPA,CLUSRX,CLUSRY,CLUSRZ,
!$OMP+           NL,MASAP,U2DM,U3DM,U4DM,VCM2,VX,VY,VZ,ACHE,
!$OMP+           PI,RETE,ROTE,VCMAX,MCMAX,RCMAX,M200,R200,CONTRASTEC,
!$OMP+           OMEGAZ,CGR,UM,UV,DMPCLUS,CONCENTRA,
!$OMP+           ORIPA2,ANGULARM,PABAS,IPLIP,DIMEN,
!$OMP+           EIGENVAL,NUMPARTBAS,
!$OMP+           RADIO,MASA,F2,VMAXCLUS,N_DM,NORMA),
!$OMP+   PRIVATE(I,INERTIA,REF_MIN,REF_MAX,KK_ENTERO,MASADM,KONTA,
!$OMP+           BASMAS,DIS,VCM,VVV2,VR,LIP,CONCEN,RS,
!$OMP+           KONTA2,BAS,IR,J,AADM,KK1,KK2,CONTADM,
!$OMP+           CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MASA2,DISTA,FAC,CONTAERR,
!$OMP+           JJ,DENSITOT,RADIAL,SALIDA,BAS1,BAS2,FLAG_200,
!$OMP+           VOL,DELTA2,NCAPAS,RSHELL,KONTA3,NSHELL_2,KONTA1,
!$OMP+           DENSA,DENSB,DENSC,AADMX,VKK,AA,NROT,BASEIGENVAL,BASX,
!$OMP+           BASY,BASZ,XP,YP,ZP,MP,RCLUS,RADII_ITER,IS_SUB,COUNT_1,
!$OMP+           COUNT_2,KONTA2PREV),
!$OMP+   SCHEDULE(DYNAMIC,2), DEFAULT(NONE), IF(.FALSE.)
*****************************
       DO I=1,NCLUS
****************************
       KK_ENTERO=REALCLUS(I)
       IF (KK_ENTERO.NE.0) THEN

        INERTIA=0.0
        REF_MIN=10.0e+10
        REF_MAX=-1.0
        MASADM=0.0
        KONTA=0
        BASMAS=0.0
        DIS=1.0E+10    !Distance (to the center) of the most central particle
        VCM=0.0
        VVV2=0.0    ! v propia de las particulas respecto al CM
        VR=0.0      ! v radial
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
        BASX=0.0
        BASY=0.0
        BASZ=0.0

*********************************************************************
*       RECENTERING AND COMPUTING VCM OF HALO I (SHRINKING SPHERE)
*********************************************************************

        RADII_ITER=(/2.0,1.8,1.6,1.4,1.2,1.1,1.0/)
        CMX=CLUSRX(I)
        CMY=CLUSRY(I)
        CMZ=CLUSRZ(I)

        IS_SUB=REALCLUS(I)
        IF (IS_SUB.EQ.-1) THEN
         IS_SUB=0
        ELSE IF (IS_SUB.GT.0) THEN
         IS_SUB=1
        END IF

*       First iteration
        KONTA=0
        MASADM=0.0
        BASX=0.0
        BASY=0.0
        BASZ=0.0
        BAS=RADII_ITER(1)
        RCLUS=BAS*RADIO(I)
        DO J=1,N_DM
         XP=RXPA(J)
         YP=RYPA(J)
         ZP=RZPA(J)
         MP=MASAP(J)
         AADM=SQRT((XP-CMX)**2+(YP-CMY)**2+(ZP-CMZ)**2)
         IF(AADM.LT.RCLUS) THEN
          KONTA=KONTA+1
          LIP(KONTA)=J
          MASADM=MASADM+MP

          BASX=BASX+MP*XP
          BASY=BASY+MP*YP
          BASZ=BASZ+MP*ZP
         END IF
        END DO

        CMX=BASX/MASADM
        CMY=BASY/MASADM
        CMZ=BASZ/MASADM

        DELTA2=MASADM/(ROTE*RETE**3*(4*PI/3)*RCLUS**3)
        !WRITE(*,*) I,is_sub,1,RCLUS,CMX,CMY,CMZ,DELTA2,CONTRASTEC
*       Subsequent iterations (shrinking the sphere)
        IF (DELTA2.LT.CONTRASTEC.OR.IS_SUB.EQ.1) THEN
         loop_shrink: DO ITER_SHRINK=2,7
          BAS=RADII_ITER(ITER_SHRINK)
          RCLUS=BAS*RADIO(I)
          KONTA2=0
          MASADM=0.0
          BASX=0.0
          BASY=0.0
          BASZ=0.0
          DO JJ=1,KONTA
           J=LIP(JJ)
           XP=RXPA(J)
           YP=RYPA(J)
           ZP=RZPA(J)
           MP=MASAP(J)
           AADM=SQRT((XP-CMX)**2+(YP-CMY)**2+(ZP-CMZ)**2)
           IF(AADM.LT.RCLUS) THEN
            REF_MIN=MIN(REF_MIN,AADM)
            REF_MAX=MAX(REF_MAX,AADM)
            KONTA2=KONTA2+1
            LIP(KONTA2)=J
            MASADM=MASADM+MP

            BASX=BASX+MP*XP
            BASY=BASY+MP*YP
            BASZ=BASZ+MP*ZP
           END IF
          END DO
          KONTA=KONTA2
          CMX=BASX/MASADM
          CMY=BASY/MASADM
          CMZ=BASZ/MASADM

          DELTA2=MASADM/(ROTE*RETE**3*(4*PI/3)*RCLUS**3)
    !      WRITE(*,*) I,is_sub,ITER_SHRINK,RCLUS,CMX,CMY,CMZ,DELTA2,
    ! & CONTRASTEC
          IF (DELTA2.GT.CONTRASTEC.AND.IS_SUB.EQ.0)
     &     EXIT loop_shrink

         END DO loop_shrink
        END IF !(DELTA2.LT.CONTRASTEC.OR.REALCLUS(I).GT.0)

        IF (MASADM.LE.0.0.OR.KONTA.EQ.0) THEN
         REALCLUS(I)=0
         CYCLE
        END IF

        CONTADM(1:KONTA)=0     !en principio todas estas ligadas
        CALL CENTROMASAS_PART(KONTA,CONTADM,LIP,
     &           U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &           CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MASA2)

        VCM=SQRT(VCMX**2+VCMY**2+VCMZ**2)

        CLUSRX(I)=CMX
        CLUSRY(I)=CMY
        CLUSRZ(I)=CMZ
        VX(I)=VCMX
        VY(I)=VCMY
        VZ(I)=VCMZ
        VCM2(I)=VCM
        RADIO(I)=RCLUS
        MASA(I)=MASADM*UM
        !write(*,*) '**',vcmx,vcmy,vcmz,vcm

*       Last, find the particles a more generous radius
        BAS=1.5
        RCLUS=BAS*RADIO(I)
        KONTA=0
        REF_MIN=10.0e+10
        REF_MAX=-1.0
        DO J=1,N_DM
         XP=RXPA(J)
         YP=RYPA(J)
         ZP=RZPA(J)
         AADM=SQRT((XP-CMX)**2+(YP-CMY)**2+(ZP-CMZ)**2)
         IF(AADM.LT.RCLUS) THEN
          REF_MIN=MIN(REF_MIN,AADM)
          REF_MAX=MAX(REF_MAX,AADM)
          KONTA=KONTA+1
          LIP(KONTA)=J
         END IF
        END DO
        DMPCLUS(I)=KONTA
        CONTADM(1:KONTA)=0

*********************************************************************
*       END RECENTERING AND COMPUTING VCM OF HALO I (SHRINKING SPHERE)
*********************************************************************

********************************************************************
*      UNBINDING:SCAPE VELOCITY
********************************************************************

        CONTAERR=KONTA
        DISTA=0.0
        KONTA2=0
        CALL REORDENAR(KONTA,CMX,CMY,CMZ,RXPA,RYPA,RZPA,CONTADM,LIP,
     &                 DISTA,KONTA2)

        FAC=0
        DO WHILE (CONTAERR.GT.0.OR.FAC.LT.3)
         FAC=FAC+1
         KONTA2PREV=KONTA2
         CALL UNBINDING8(FAC,I,REF_MIN,REF_MAX,DISTA,U2DM,U3DM,U4DM,
     &                   MASAP,RXPA,RYPA,RZPA,RADIO,MASA,CLUSRX,CLUSRY,
     &                   CLUSRZ,LIP,KONTA,CONTADM,VX,VY,VZ,REALCLUS)
         CALL REORDENAR(KONTA,CMX,CMY,CMZ,RXPA,RYPA,RZPA,CONTADM,LIP,
     &                  DISTA,KONTA2)
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
     &                        RYPA,RZPA,MASAP,RADIO,MASA,CLUSRX,CLUSRY,
     &                        CLUSRZ,LIP,KONTA,CONTADM,VX,VY,VZ,
     &                        REALCLUS)
         CALL REORDENAR(KONTA,CMX,CMY,CMZ,RXPA,RYPA,RZPA,CONTADM,LIP,
     &                  DISTA,KONTA2)
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
        IF (KONTA2.LT.NUMPARTBAS) THEN
         REALCLUS(I)=0
        ELSE
**************************************************************
*      DENSITY PROFILE
**************************************************************
         REF_MIN=DISTA(1)
         REF_MAX=DISTA(KONTA2)
         JJ=0
         DENSITOT=0.0
         RADIAL=0.0
****************** stopping condition
         SALIDA=0
******************
         BAS=0.0 ! to store cummulative mass profile
         BAS1=-1.0 ! to store vcmax
         BAS2=0.0 ! to store vc(r)
*        VCMAX=-1.0
         FLAG_200=0 ! whether r200c has been found

         DO J=1,KONTA2      !!!!! DEJO 80 por ciento de BINS DE SEGURIDAD
          VOL=PI*(4.0/3.0)*(DISTA(J)*RETE)**3
          BAS=BAS+(MASAP(LIP(J))/NORMA)

          DELTA2=NORMA*BAS/VOL/ROTE
          BAS2=BAS/DISTA(J) ! vc(r)

          IF (BAS2.GT.BAS1) THEN
           VCMAX(I)=BAS2
           MCMAX(I)=BAS
           RCMAX(I)=DISTA(J)
           BAS1=BAS2
          END IF

          IF (FLAG_200.EQ.0) THEN
           IF (DELTA2.LE.200.0/OMEGAZ) THEN
            M200(I)=DELTA2*VOL*ROTE
            R200(I)=DISTA(J)
            FLAG_200=1
           END IF
          END IF

          IF (DELTA2.LE.CONTRASTEC.AND.J.GT.INT(0.1*KONTA2)) THEN
           SALIDA=1
           EXIT
          END IF
         END DO ! J=1,KONTA2

         VCMAX(I)=VCMAX(I)*NORMA*CGR/RETE
         VCMAX(I)=SQRT(VCMAX(I))*UV
         MCMAX(I)=MCMAX(I)*NORMA*UM
         RCMAX(I)=RCMAX(I)   !*RETE
         M200(I)=M200(I)*UM
         R200(I)=R200(I)     !*RETE

         IF (SALIDA.EQ.1) THEN
          NCAPAS(I)=J
          IF (J.EQ.KONTA2) THEN
           RSHELL=DISTA(J)
          ELSE
           RSHELL=DISTA(J+1)
          END IF
         END IF

         IF (SALIDA.NE.1.AND.KONTA2.NE.0) THEN   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

*         hasta aqui masa acumulada en bins de 10 particulas
          NSHELL_2=JJ
          write(*,*) 'nshell2=',nshell_2

          IF (NSHELL_2.GT.4) THEN
           DO J=INT(0.8*NSHELL_2),NSHELL_2       !!!!! DEJO 80 por cient de BINS DE SEGURIDA
            BAS=0.5*(RADIAL(J)+RADIAL(J-1))  !radio medio
            DENSA=(DENSITOT(J)-DENSITOT(J-1))/(RADIAL(J)-RADIAL(J-1))
            DENSA=DENSA/BAS/BAS

            BAS=0.5*(RADIAL(J+1)+RADIAL(J))  !radio medio
            DENSB=(DENSITOT(J+1)-DENSITOT(J))/(RADIAL(J+1)-RADIAL(J))
            DENSB=DENSB/BAS/BAS

            BAS=0.5*(RADIAL(J+2)+RADIAL(J+1))  !radio medio
            DENSC=(DENSITOT(J+2)-DENSITOT(J+1))/
     &                                       (RADIAL(J+2)-RADIAL(J+1))
            DENSC=DENSC/BAS/BAS

            IF (DENSC.GT.DENSB.AND.DENSB.GT.DENSA) THEN
             SALIDA=2
             EXIT
            END IF
           END DO            !!!!!!!!!!!!!!!
          END IF ! NSHELL_2

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
         DMPCLUS(I)=0
         DO J=1,KONTA
          IF (CONTADM(J).EQ.0) THEN

           AADMX(1)=RXPA(LIP(J))-CLUSRX(I)
           AADMX(2)=RYPA(LIP(J))-CLUSRY(I)
           AADMX(3)=RZPA(LIP(J))-CLUSRZ(I)
           AADM=SQRT(AADMX(1)**2+AADMX(2)**2+AADMX(3)**2)

CV2
           VVV2=0.0
           VVV2=(U2DM(LIP(J))-VX(I))**2
     &         +(U3DM(LIP(J))-VY(I))**2
     &         +(U4DM(LIP(J))-VZ(I))**2
CV2

           IF (AADM.LE.RSHELL) THEN
            KONTA2=KONTA2+1
            BASMAS=BASMAS+(MASAP(LIP(J))/NORMA)

*********anyadido susana
            DMPCLUS(I)=DMPCLUS(I)+1

            ANGULARM(I)=ANGULARM(I)+MASAP(LIP(J))*AADM*SQRT(VVV2)

**          INERTIA TENSOR
            DO JY=1, 3
            DO IX=1, 3
              INERTIA(IX,JY)=INERTIA(IX,JY)+AADMX(IX)*AADMX(JY)
            END DO
            END DO

**        VELOCITY OF THE FASTEST PARTICLE IN THE HALO
            VKK=0.0
            VKK=SQRT(U2DM(LIP(J))**2+U3DM(LIP(J))**2+
     &               U4DM(LIP(J))**2)

            IF (VKK.GT.VMAXCLUS(I)) THEN
             VMAXCLUS(I)=VKK
            END IF

**        CLOSEST PARTICLE TO THE CENTER OF THE HALO
            IF (AADM.LT.DIS) THEN
             DIS=AADM
             IPLIP(I)=ORIPA2(LIP(J))
            END IF

            IF(AADM.NE.0.0) THEN
             AA=0.0
             AA=((RXPA(LIP(J))-CLUSRX(I))/AADM)*U2DM(LIP(J))+
     &          ((RYPA(LIP(J))-CLUSRY(I))/AADM)*U3DM(LIP(J))+
     &          ((RZPA(LIP(J))-CLUSRZ(I))/AADM)*U4DM(LIP(J))
             VR=VR+AA*MASAP(LIP(J))
            END IF
           END IF    !AADM.LT.RSHELL
          END IF     !CONTADM
         END DO      !KONTA


*       CONTROL DE SEGURIDAD
         IF(DMPCLUS(I).NE.KONTA2) THEN
          WRITE(*,*) 'WARNING!', DMPCLUS(I),KONTA2
          STOP
         ENDIF


*******************************************************
*      SAVING MASSES, RADII, PROFILES AND SHAPES...
*******************************************************

         MASA(I)=BASMAS*NORMA*9.1717E18    !!MASA2
         ANGULARM(I)=ANGULARM(I)/BASMAS
         RADIO(I)=RSHELL

         INERTIA(1:3,1:3)=INERTIA(1:3,1:3)/DMPCLUS(I)
         BASEIGENVAL(1:3)=0.0

         IF (DMPCLUS(I).GE.NUMPARTBAS) THEN
          CALL JACOBI(INERTIA,DIMEN,BASEIGENVAL,NROT)
          CALL SORT(BASEIGENVAL,DIMEN,DIMEN)
         END IF

         DO II=1, DIMEN
          EIGENVAL(II,I)=BASEIGENVAL(II)
          EIGENVAL(II,I)=SQRT(EIGENVAL(II,I))
         END DO

******************************************
************ FIN HALO I ******************
******************************************
        END IF ! (ELSE of KONTA2.LT.NUMPARTBAS) (poor haloes after unbinding)
       END IF ! (realclus(i).ne.0)

*****************
       END DO   !I=IP,IP2
****************


       RETURN
       END




**********************************************************************
       SUBROUTINE REORDENAR(KONTA,CMX,CMY,CMZ,
     &                      RXPA,RYPA,RZPA,CONTADM,LIP,DISTA,KONTA2)
**********************************************************************
*      Sorts the particles with increasing distance to the center of
*      the halo. Only particles with CONTADM=0 are sorted (the others
*      are already pruned particles and therefore they are ignored)
**********************************************************************

       IMPLICIT NONE

       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER I,J,K,KONTA,KONTA2,JJ

*      ---HALOS Y SUBHALOS---
       REAL*4 CMX,CMY,CMZ

*      ---PARTICULAS E ITERACIONES---
       INTEGER LIP(PARTIRED),CONTADM(PARTIRED)

       REAL*4 RXPA(PARTIRED)
       REAL*4 RYPA(PARTIRED)
       REAL*4 RZPA(PARTIRED)

c       INTEGER CONTADM(PARTI)

       REAL*4 DISTA(0:PARTIRED)

       INTEGER INDICE(KONTA)
       REAL*4 DISTA2(0:KONTA)
       INTEGER QUIEN(KONTA)

       REAL*4 AADMX,AADMY,AADMZ,AADM

       QUIEN=0
       DISTA=1.E10
       INDICE=0
       DISTA2=0.0

       KONTA2=0
       DO J=1,KONTA
        IF (CONTADM(J).EQ.0) THEN
         KONTA2=KONTA2+1
         JJ=LIP(J)
         AADMX=RXPA(JJ)-CMX
         AADMY=RYPA(JJ)-CMY
         AADMZ=RZPA(JJ)-CMZ
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

       CONTADM=1
       CONTADM(1:KONTA2)=0

       RETURN
       END

***********************************************************
       SUBROUTINE CENTROMASAS_PART(N,CONTADM,LIP,
     &            U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &            CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MASA)
***********************************************************
*      Computes the center of mass of the particles
*      within the halo
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

       INTEGER LIP(PARTIRED),CONTADM(PARTIRED)

       REAL*4 CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MASA

*      ---- DOUBLE PRECISION -----------------------
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

***********************************************************
       SUBROUTINE UNBINDING8(FAC,I,REF_MIN,REF_MAX,DISTA,
     &           U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &           RADIO,MASA,CLUSRX,CLUSRY,CLUSRZ,
     &           LIP,KONTA,CONTADM,VX,VY,VZ,REALCLUS)
***********************************************************
*      Finds and discards the unbound particles (those
*      with speed larger than the scape velocity).
*      Potential is computed in double precision.
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
       REAL*4 RADIO(NMAXNCLUS)
       REAL*4 MASA(NMAXNCLUS)
       REAL*4 CLUSRX(NMAXNCLUS)
       REAL*4 CLUSRY(NMAXNCLUS)
       REAL*4 CLUSRZ(NMAXNCLUS)
       REAL*4 VCM2(NMAXNCLUS) !!?
       REAL*4 VX(NMAXNCLUS)
       REAL*4 VY(NMAXNCLUS)
       REAL*4 VZ(NMAXNCLUS)
       INTEGER REALCLUS(MAXNCLUS)

*      ---PARTICULAS E ITERACIONES---
       INTEGER LIP(PARTIRED),CONTADM(PARTIRED)
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

       REAL*4 DISTA(0:PARTIRED)

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

       KONTA2=COUNT(CONTADM(1:KONTA).EQ.0)

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
          BAS8=DISTA(J)-DISTA(J-1)
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
        CALL CENTROMASAS_PART(KONTA,CONTADM,LIP,
     &           U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &           CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MMM)

       END IF

       KONTA3=COUNT(CONTADM(1:KONTA).EQ.0)

       IF (KONTA3.EQ.0) THEN

        CLUSRX(I)=0.0
        CLUSRY(I)=0.0
        CLUSRZ(I)=0.0

        VX(I)=0.0
        VY(I)=0.0
        VZ(I)=0.0

        MASA(I)=0.0
        RADIO(I)=0.0

        REALCLUS(I)=0

       ELSE

        CLUSRX(I)=CMX
        CLUSRY(I)=CMY
        CLUSRZ(I)=CMZ

        VX(I)=VCMX
        VY(I)=VCMY
        VZ(I)=VCMZ

        MASA(I)=MMM*9.1717e+18


*       estimacion nuevo radio

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


CX       write(*,*) 'new_r',RADIO(I),REF_MAX
        RADIO(I)=REF_MAX

       END IF

       RETURN
       END


***********************************************************
       SUBROUTINE UNBINDING_SIGMA(FAC,I,REF_MIN,REF_MAX,U2DM,U3DM,U4DM,
     &                            RXPA,RYPA,RZPA,MASAP,RADIO,MASA,
     &                            CLUSRX,CLUSRY,CLUSRZ,LIP,KONTA,
     &                            CONTADM,VX,VY,VZ,REALCLUS)
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
       REAL*4 RADIO(NMAXNCLUS),MASA(NMAXNCLUS)
       REAL*4 CLUSRX(NMAXNCLUS),CLUSRY(NMAXNCLUS),CLUSRZ(NMAXNCLUS)
       INTEGER LIP(PARTIRED),CONTADM(PARTIRED)
       INTEGER KONTA
       REAL*4 VX(NMAXNCLUS),VY(NMAXNCLUS),VZ(NMAXNCLUS)
       INTEGER REALCLUS(MAXNCLUS)

       REAL CMX,CMY,CMZ,VXCM,VYCM,VZCM,BAS,SIGMA2,BB,AADM,AADMX,AADMY
       REAL AADMZ,MMM
       INTEGER J,JJ,KONTA2,KONTA3
       REAL,ALLOCATABLE::DESV2(:)

       BB=MAX(6.0-1.0*(FAC-1), 3.0)
       BB=BB**2 ! This is because we compare velocities squared

       CMX=CLUSRX(I)
       CMY=CLUSRY(I)
       CMZ=CLUSRZ(I)
       VXCM=VX(I)
       VYCM=VY(I)
       VZCM=VZ(I)

       KONTA2=COUNT(CONTADM(1:KONTA).EQ.0)

       ALLOCATE(DESV2(1:KONTA2))

       IF (KONTA2.GT.0) THEN
        SIGMA2=0.0
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
        CALL CENTROMASAS_PART(KONTA,CONTADM,LIP,
     &           U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &           CMX,CMY,CMZ,VXCM,VYCM,VZCM,MMM)

       END IF

       KONTA3=COUNT(CONTADM(1:KONTA2).EQ.0)

       IF (KONTA3.EQ.0) THEN

        CLUSRX(I)=0.0
        CLUSRY(I)=0.0
        CLUSRZ(I)=0.0

        VX(I)=0.0
        VY(I)=0.0
        VZ(I)=0.0

        MASA(I)=0.0
        RADIO(I)=0.0

        REALCLUS(I)=0

       ELSE

        CLUSRX(I)=CMX
        CLUSRY(I)=CMY
        CLUSRZ(I)=CMZ

        VX(I)=VXCM
        VY(I)=VYCM
        VZ(I)=VZCM

        MASA(I)=MMM*9.1717e+18

*       estimacion nuevo radio

        REF_MIN=10.0E+10
        REF_MAX=-1.0

        DO J=1,KONTA2
         IF (CONTADM(J).EQ.0) THEN
           AADMX=RXPA(LIP(J))-CMX
           AADMY=RYPA(LIP(J))-CMY
           AADMZ=RZPA(LIP(J))-CMZ

           AADM=SQRT(AADMX**2+AADMY**2+AADMZ**2)

           REF_MIN=MIN(REF_MIN,AADM)
           REF_MAX=MAX(REF_MAX,AADM)
         END IF
        END DO

        RADIO(I)=REF_MAX

       END IF

       DEALLOCATE(DESV2)


       RETURN
       END
