********************************************************************
      SUBROUTINE STELLAR_HALOES(NCLUS,MASA,RADIO,MSUB,RSUB,REALCLUS,
     &                          DMPCLUS,CLUSRX,CLUSRY,CLUSRZ,RXPA,RYPA,
     &                          RZPA,MASAP,U2DM,U3DM,U4DM,ORIPA,N_DM,
     &                          N_ST,NX,LADO0,PARTICLES_PER_HALO,
     &                          INDCS_PARTICLES_PER_HALO,UM,
     &                          MIN_NUM_PART_ST)
********************************************************************

      IMPLICIT NONE
      INCLUDE 'input_files/asohf_parameters.dat'

      INTEGER NCLUS
      REAL*4 MASA(MAXNCLUS),RADIO(MAXNCLUS)
      REAL*4 MSUB(MAXNCLUS),RSUB(MAXNCLUS)
      INTEGER REALCLUS(MAXNCLUS), DMPCLUS(MAXNCLUS)
      REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
      REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED)
      REAL*4 MASAP(PARTIRED)
      REAL*4 U2DM(PARTIRED),U3DM(PARTIRED),U4DM(PARTIRED)
      INTEGER ORIPA(PARTIRED)
      INTEGER N_DM,N_ST,NX
      REAL*4 LADO0
      INTEGER PARTICLES_PER_HALO(PARTIRED)
      INTEGER INDCS_PARTICLES_PER_HALO(2,NMAXNCLUS)
      REAL*4 UM
      INTEGER MIN_NUM_PART_ST

      INTEGER NSTPART_X(0:NMAX),I,LOWP1,LOWP2,J,JJ,NST_HALO,NDM_HALO
      INTEGER MAX_NUM_PART_LOCAL,WELL_ALLOCATED,MINORIPA,MAXORIPA
      INTEGER NPART_HALO,BASINT,KONTA,KONTA2,FAC,CONTAERR
      INTEGER COUNT_1,COUNT_2,KONTA2PREV
      REAL XLDOM,CX,CY,CZ,RCLUS,RCLUS2,XP,YP,ZP,CMX,CMY,CMZ
      REAL VCMX,VCMY,VCMZ,BASMAS,REF_MIN,REF_MAX

      INTEGER,ALLOCATABLE::ORIPADM_LOT(:)
      INTEGER,ALLOCATABLE::LIP(:),CONTADM(:),LIPST(:)
      REAL,ALLOCATABLE::DISTA(:),DISTAST(:)

      WRITE(*,*) 'DM, stellar particles:', N_DM, N_ST

      XLDOM=-LADO0/2.0

**********************************************************************
*     Sort stellar particles
**********************************************************************
      CALL SORT_STELLAR_PARTICLES_X(U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &                              ORIPA,N_DM,N_ST,NSTPART_X,NX,LADO0)

**********************************************************************
*     Build ORIPA_DM Look-up table to get particles from the oripas
**********************************************************************
      MINORIPA=MINVAL(ORIPA(1:N_DM))
      MAXORIPA=MAXVAL(ORIPA(1:N_DM))
      ALLOCATE(ORIPADM_LOT(MINORIPA:MAXORIPA))
!$OMP PARALLEL DO SHARED(MINORIPA,MAXORIPA,ORIPADM_LOT),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
      DO I=MINORIPA,MAXORIPA
       ORIPADM_LOT(I)=-1
      END DO

!$OMP PARALLEL DO SHARED(N_DM,ORIPADM_LOT,ORIPA),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
      DO I=1,N_DM
       ORIPADM_LOT(ORIPA(I))=I
      END DO

      WRITE(*,*) 'ORIPA LOT DONE',MINORIPA,MAXORIPA

**********************************************************************
*     Main loop through DM haloes
**********************************************************************

      DO I=1,NCLUS
       IF (REALCLUS(I).EQ.0) CYCLE

       !***********************************************
       !!! IDENTIFY STELLAR PARTICLES INSIDE THE HALO
       !***********************************************
       CX=CLUSRX(I)
       CY=CLUSRY(I)
       CZ=CLUSRZ(I)
       IF (REALCLUS(I).EQ.-1) THEN
        RCLUS=RADIO(I)
       ELSE
        RCLUS=RSUB(I)
       END IF
       RCLUS2=RCLUS**2

       CALL FIND_PARTICLE_INDICES(CX,RCLUS,XLDOM,NSTPART_X,LOWP1,LOWP2)

       MAX_NUM_PART_LOCAL=0
       DO J=LOWP1,LOWP2
        IF ((RXPA(J)-CX)**2+(RYPA(J)-CY)**2+(RZPA(J)-CZ)**2.LT.RCLUS2)
     &    MAX_NUM_PART_LOCAL=MAX_NUM_PART_LOCAL+1
       END DO  !J=LOWP1,LOWP2

C       WRITE(*,*) 'HALO, NUM STARS:',I,MAX_NUM_PART_LOCAL,LOWP1,LOWP2
       IF (MAX_NUM_PART_LOCAL.EQ.0) CYCLE

       ALLOCATE(LIPST(MAX_NUM_PART_LOCAL))

       JJ=0
       DO J=LOWP1,LOWP2
        IF ((RXPA(J)-CX)**2+(RYPA(J)-CY)**2+(RZPA(J)-CZ)**2
     &      .LT.RCLUS2) THEN
         JJ=JJ+1
         LIPST(JJ)=J
        END IF
       END DO  !J=LOWP1,LOWP2

       IF (JJ.NE.MAX_NUM_PART_LOCAL) THEN
        WRITE(*,*) 'Wrong allocation of stars',JJ,MAX_NUM_PART_LOCAL
        STOP
       END IF

       IF (JJ.LT.MIN_NUM_PART_ST) THEN
        DEALLOCATE(LIPST)
        CYCLE
       END IF

       !***********************************************
       !!! RESCUE DM PARTICLES
       !***********************************************
       NDM_HALO=DMPCLUS(I)
       NST_HALO=MAX_NUM_PART_LOCAL
       NPART_HALO=NDM_HALO+NST_HALO

       ALLOCATE(LIP(NPART_HALO),CONTADM(NPART_HALO),DISTA(0:NPART_HALO))

       LOWP1=INDCS_PARTICLES_PER_HALO(1,I)
       LOWP2=INDCS_PARTICLES_PER_HALO(2,I)
       JJ=0
       DO J=LOWP1,LOWP2
        JJ=JJ+1
        LIP(JJ)=ORIPADM_LOT(PARTICLES_PER_HALO(J))
       END DO
       DO J=1,NST_HALO
        JJ=JJ+1
        LIP(JJ)=LIPST(J)
       END DO
       IF (JJ.NE.NPART_HALO) THEN
        WRITE(*,*) 'Problem with halo',I,JJ,NPART_HALO
        STOP
       END IF
c       WRITE(*,*) I,MINVAL(DISTA(1:NPART_HALO)),
c     &            MAXVAL(DISTA(1:NPART_HALO)),RCLUS

       !***********************************************
       !!! SORT DM AND STELLAR PARTICLES ALTOGETHER
       !***********************************************
       ! Stellar particles will have LIP>N_DM
       KONTA=NPART_HALO
       BASINT=KONTA
       CONTADM(1:KONTA)=0
       CALL REORDENAR(KONTA,CX,CY,CZ,RXPA,RYPA,RZPA,CONTADM,LIP,
     &                DISTA,KONTA2,1,NPART_HALO,BASINT)

C       JJ=-1
C       DO J=2,NPART_HALO
C        IF (DISTA(J).LT.DISTA(J-1)) THEN
C         WRITE(*,*) 'HALO',I,'UNSORTED PARTICLES!',J,DISTA(J),DISTA(J-1)
C         EXIT
C        END IF
C       END DO

C       IF (COUNT(LIP.GT.N_DM).NE.NST_HALO) THEN
C        WRITE(*,*) 'HAVE LOST STARS,',COUNT(LIP.GT.N_DM),NST_HALO
C        STOP
C       END IF

       !***********************************************
       !!! UNBINDING: SCAPE VELOCITY
       !***********************************************

       REF_MIN=DISTA(1)
       REF_MAX=DISTA(NPART_HALO)

       CALL CENTROMASAS_PART(KONTA,CONTADM,LIP,U2DM,U3DM,U4DM,MASAP,
     &                       RXPA,RYPA,RZPA,CMX,CMY,CMZ,VCMX,VCMY,VCMZ,
     &                       BASMAS,NPART_HALO)
       WRITE(*,*) VCMX,VCMY,VCMZ

       FAC=0
       DO WHILE (CONTAERR.GT.0.OR.FAC.LT.3)
        FAC=FAC+1
        KONTA2PREV=KONTA2
        CALL UNBINDING8_STARS(FAC,REF_MIN,REF_MAX,DISTA,U2DM,U3DM,
     &                  U4DM,MASAP,RXPA,RYPA,RZPA,LIP,KONTA,
     &                  CONTADM,KONTA2,NPART_HALO,UM,VCMX,VCMY,VCMZ,
     &                  N_DM)
        BASINT=KONTA
        CALL REORDENAR(KONTA,CX,CY,CZ,RXPA,RYPA,RZPA,CONTADM,LIP,
     &                 DISTA,KONTA2,0,NPART_HALO,BASINT)
        REF_MAX=DISTA(KONTA2)
        REF_MIN=DISTA(1)
        CONTAERR=KONTA2PREV-KONTA2
       END DO

       count_1=konta-konta2
       count_2=konta2 !backup
C       write(*,*) 'Unbinding V_ESC',i,'. ',konta-ndm_halo,'-->',
C     &             konta2-ndm_halo,'. Pruned:',count_1,'. Iters:', FAC

       !***********************************************
       !!! GET RID OF DM PARTICLES (WE NO LONGER WANT THEM)
       !***********************************************
       NST_HALO=COUNT(LIP.GT.N_DM.AND.CONTADM.EQ.0)
       IF (NST_HALO.LT.MIN_NUM_PART_ST) THEN
         DEALLOCATE(LIPST,LIP,CONTADM,DISTA)
         CYCLE
       END IF
       DEALLOCATE(LIPST)
       ALLOCATE(DISTAST(0:NST_HALO),LIPST(NST_HALO))

       JJ=0
       DO J=1,KONTA2
        IF (LIP(J).GT.N_DM.AND.CONTADM(J).EQ.0) THEN
         JJ=JJ+1
         LIPST(JJ)=LIP(J)
         DISTAST(JJ)=DISTA(J)
        END IF
       END DO
       KONTA2=JJ
C       write(*,*) 'now we have stellar particles:',KONTA2,NST_HALO

       DEALLOCATE(LIP,DISTA,CONTADM)
       ALLOCATE(CONTADM(NST_HALO))
       CONTADM=1
       CONTADM(1:KONTA2)=0

       CALL CENTROMASAS_PART(KONTA2,CONTADM,LIPST,
     &          U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &          CMX,CMY,CMZ,VCMX,VCMY,VCMZ,BASMAS,NST_HALO)

       !***********************************************
       !!! UNBINDING: PHASE SPACE
       !***********************************************
       FAC=0
       CONTAERR=KONTA2
       KONTA=KONTA2
       DO WHILE (CONTAERR.GT.0.OR.FAC.LT.4)
        FAC=FAC+1
        KONTA2PREV=KONTA2
        CALL UNBINDING_SIGMA_STARS(FAC,REF_MIN,REF_MAX,U2DM,U3DM,U4DM,
     &               RXPA,RYPA,RZPA,MASAP,LIPST,CONTADM,KONTA2,
     &               NST_HALO,UM,VCMX,VCMY,VCMZ,N_DM)
        BASINT=KONTA
        CALL REORDENAR(KONTA,CX,CY,CZ,RXPA,RYPA,RZPA,CONTADM,LIPST,
     &                 DISTAST,KONTA2,0,NST_HALO,BASINT)
        REF_MAX=DISTAST(KONTA2)
        REF_MIN=DISTAST(1)
        CONTAERR=KONTA2PREV-KONTA2
        !write(*,*) 'sigma unbinding: iter,unbound',fac,contaerr
       END DO

       count_2=konta-konta2
C       write(*,*) 'Unbinding SIGMA',i,'. ',konta,'-->',konta2,
C     &            '. Pruned:',count_2,'. Iters:', FAC
C       write(*,*) '--'

       IF (KONTA2.LT.MIN_NUM_PART_ST) THEN
        DEALLOCATE(CONTADM,LIPST,DISTAST)
        CYCLE
       END IF

       WRITE(*,*) 'ACCEPTED STELLAR HALO',I,KONTA2




       DEALLOCATE(LIPST,CONTADM,DISTAST)
      END DO

      RETURN
      END

***********************************************************
       SUBROUTINE UNBINDING_SIGMA_STARS(FAC,REF_MIN,REF_MAX,U2DM,U3DM,
     &              U4DM,RXPA,RYPA,RZPA,MASAP,LIP,CONTADM,
     &              KONTA2,MAX_NUM_PART,UM,VCMX,VCMY,VCMZ,N_DM)
***********************************************************
*      Finds and discards the unbound particles (those
*      with speed larger than the scape velocity).
*      VERSION FOR STARS (dark matter particles have already been treated)
***********************************************************

       IMPLICIT NONE

       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER FAC
       REAL REF_MIN,REF_MAX
       REAL*4 U2DM(PARTIRED),U3DM(PARTIRED),U4DM(PARTIRED)
       REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       INTEGER LIP(MAX_NUM_PART),CONTADM(MAX_NUM_PART)
       INTEGER KONTA2,MAX_NUM_PART
       REAL*4 UM,VCMX,VCMY,VCMZ
       INTEGER N_DM

       REAL CMX,CMY,CMZ,BAS,BB,AADM,AADMX,AADMY
       REAL AADMZ,MMM
       INTEGER J,JJ,KONTA3
       REAL,ALLOCATABLE::DESV2(:)

*      DOUBLE PRECISION VARIABLES
       REAL*8 SIGMA2

       BB=MAX(6.0-1.0*(FAC-1), 3.0)
       BB=BB**2 ! This is because we compare velocities squared

       IF (KONTA2.GT.0) THEN

        ALLOCATE(DESV2(1:KONTA2))

        SIGMA2=0.D0
        DO J=1,KONTA2
         JJ=LIP(J)
         BAS=(U2DM(JJ)-VCMX)**2+(U3DM(JJ)-VCMY)**2+(U4DM(JJ)-VCMZ)**2
         DESV2(J)=BAS
         SIGMA2=SIGMA2+BAS
        END DO

        IF (KONTA2.GT.1) SIGMA2=SIGMA2/(KONTA2-1)

*       Find particles with too large relative velocity
        DO J=1,KONTA2
         IF (DESV2(J).GT.BB*SIGMA2) CONTADM(J)=1
        END DO

*       NEW CENTER OF MASS AND ITS VELOCITY
        CALL CENTROMASAS_PART(KONTA2,CONTADM,LIP,
     &           U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &           CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MMM,MAX_NUM_PART)

        DEALLOCATE(DESV2)

       END IF


       RETURN
       END

***********************************************************
       SUBROUTINE UNBINDING8_STARS(FAC,REF_MIN,REF_MAX,DISTA,
     &           U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &           LIP,KONTA,CONTADM,KONTA2,MAX_NUM_PART,UM,
     &           VCMX,VCMY,VCMZ,N_DM)
***********************************************************
*      Finds and discards the unbound particles (those
*      with speed larger than the scape velocity).
*      Potential is computed in double precision.
*      VERSION FOR STARS (dark matter particles have already been treated)
***********************************************************

       IMPLICIT NONE

       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER FAC
       REAL*4 REF_MIN,REF_MAX
       REAL*4 DISTA(0:MAX_NUM_PART)
       REAL*4 U2DM(PARTIRED),U3DM(PARTIRED),U4DM(PARTIRED)
       REAL*4 MASAP(PARTIRED)
       REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED)
       INTEGER LIP(MAX_NUM_PART),CONTADM(MAX_NUM_PART)
       INTEGER KONTA,KONTA2,MAX_NUM_PART
       REAL*4 UM
       REAL*4 VCMX,VCMY,VCMZ
       INTEGER N_DM

       INTEGER J,K,IX,IMAX,JJ

       REAL*4 REI,CGR,PI,PI4ROD
       COMMON /CONS/PI4ROD,REI,CGR,PI

       REAL*4 RETE,HTE,ROTE
       COMMON /BACK/ RETE,HTE,ROTE

       INTEGER KONTA3
       REAL*4 VVV2,VESC2,AADMX(3),AADM,DR, AA, BB, CC
       REAL*4 BAS
       REAL*4 CMX,CMY,CMZ,MMM
       REAL*4 POTOK

*!!!!! ESPECIAL DOBLE PRECISON !!!!!!!!!!!!!!!!!!!!!
COJO       REAL*8 POT(KONTA)
       REAL*8 POT(0:KONTA)
       REAL*8 POT1
       REAL*8 BAS8
       REAL*8 MASA8
       REAL*8 AA8
***********************************************

       POT=0.D0

       IF (KONTA2.GT.0) THEN
*      Max mass
       MASA8=DBLE(MASAP(1))

       !POT(1)=MASA8/DBLE(DISTA(1))
       !WRITE(*,*) 'IN UNBINDING, KONTA2=',KONTA2
       DO J=1,KONTA2
        IF (DISTA(J).GT.0.01*REF_MAX) EXIT
       END DO
       JJ=J
       MASA8=0.D0
       DO J=1,JJ
        MASA8=MASA8+DBLE(MASAP(LIP(J)))
       END DO
       DO J=1,JJ
        POT(J)=MASA8/DISTA(JJ)
       END DO

       DO J=JJ+1,KONTA2
         MASA8=MASA8+DBLE(MASAP(LIP(J)))
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

       AA8=DBLE(CGR/RETE)

       BB=2.0
       IF (FAC.EQ.1) BB=8.0
       IF (FAC.EQ.2) BB=4.0

       BB=BB**2 !(we compare the squared velocities)

*      Find particles able to escape the potential well
       DO J=1,KONTA2
        IF (LIP(J).LE.N_DM) CYCLE

        POTOK=(POT(J)-POT1)*AA8
        VESC2=2.0*ABS(POTOK)

        VVV2=(U2DM(LIP(J))-VCMX)**2
     &      +(U3DM(LIP(J))-VCMY)**2
     &      +(U4DM(LIP(J))-VCMZ)**2

C        WRITE(*,*) SQRT(VVV2),SQRT(VESC2),DISTA(J)/REF_MAX,
C     &              DISTA(J),REF_MAX
        IF (VVV2.GT.BB*VESC2)  CONTADM(J)=1
       END DO

*      NEW CENTER OF MASS AND ITS VELOCITY
       CALL CENTROMASAS_PART(KONTA2,CONTADM,LIP,
     &          U2DM,U3DM,U4DM,MASAP,RXPA,RYPA,RZPA,
     &          CMX,CMY,CMZ,VCMX,VCMY,VCMZ,MMM,MAX_NUM_PART)

       END IF

       RETURN
       END


*********************************************************************
       SUBROUTINE SORT_STELLAR_PARTICLES_X(U2DM,U3DM,U4DM,MASAP,RXPA,
     &                  RYPA,RZPA,ORIPA,N_DM,N_ST,NSTPART_X,NX,LADO0)
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
       INTEGER N_DM,N_ST,NSTPART_X(0:NMAX),NX
       REAL*4 LADO0

       INTEGER I,CONTA,IX,IXLAST
       REAL XL,DX,RADXR(0:NMAX)
       INTEGER,ALLOCATABLE::INDICES(:)
       REAL,ALLOCATABLE::SCR(:,:)
       INTEGER,ALLOCATABLE::SCRINT(:,:)

       WRITE(*,*) 'Sorting stellar particles by X coordinate',
     &            ' (for faster search)'

       DO I=1,NX
        NSTPART_X(I)=0
       END DO

       ALLOCATE(INDICES(1:N_ST))
       CALL INDEXX(N_ST,RXPA(N_DM+1:N_DM+N_ST),INDICES)

       ALLOCATE(SCR(1:7,1:N_ST), SCRINT(1,1:N_ST))

!$OMP PARALLEL DO SHARED(SCR,SCRINT,RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,
!$OMP+                   ORIPA,INDICES,N_DM,N_ST),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
       DO I=1,N_ST
        SCR(1,I)=RXPA(N_DM+INDICES(I))
        SCR(2,I)=RYPA(N_DM+INDICES(I))
        SCR(3,I)=RZPA(N_DM+INDICES(I))
        SCR(4,I)=U2DM(N_DM+INDICES(I))
        SCR(5,I)=U3DM(N_DM+INDICES(I))
        SCR(6,I)=U4DM(N_DM+INDICES(I))
        SCR(7,I)=MASAP(N_DM+INDICES(I))
        SCRINT(1,I)=ORIPA(N_DM+INDICES(I))
       END DO

       DEALLOCATE(INDICES)

!$OMP PARALLEL DO SHARED(SCR,SCRINT,RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,
!$OMP+                   ORIPA,INDICES,N_DM,N_ST),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
       DO I=1,N_ST
        RXPA(N_DM+I)=SCR(1,I)
        RYPA(N_DM+I)=SCR(2,I)
        RZPA(N_DM+I)=SCR(3,I)
        U2DM(N_DM+I)=SCR(4,I)
        U3DM(N_DM+I)=SCR(5,I)
        U4DM(N_DM+I)=SCR(6,I)
        MASAP(N_DM+I)=SCR(7,I)
        ORIPA(N_DM+I)=SCRINT(1,I)
       END DO

       DEALLOCATE(SCR,SCRINT)

       XL=-LADO0/2
       DX=LADO0/NX
       DO I=0,NX
        !RADXL(I)=XL+(I-1)*DX ! X left interface of cell I
        RADXR(I)=XL+I*DX ! X right interface of cell I
       END DO

       NSTPART_X(0)=N_DM
       IX=1
       DO I=N_DM+1,N_DM+N_ST
        IF (RXPA(I).GT.RADXR(IX)) THEN
         NSTPART_X(IX)=I-1
         IX=IX+1
         DO WHILE (RXPA(I).GT.RADXR(IX))
          NSTPART_X(IX)=I-1
          IX=IX+1
         END DO
        END IF
       END DO

       IXLAST=IX
       DO IX=IXLAST,NX
        NSTPART_X(IX)=I-1
       END DO

C       WRITE(*,*) 'CHECKING'
C       DO IX=1,NX
C        WRITE(*,*) IX,':',NSTPART_X(IX-1)+1,NSTPART_X(IX),':',
C     &             radxr(ix-1),radxr(ix),NSTPART_X(IX)-NSTPART_X(IX-1)
C        if (ix.eq.1) then
C         WRITE(*,*) ' --> just before left',' first particle'
C        else
C         WRITE(*,*) ' --> just before left',RXPA(NSTPART_X(IX-1))
C        end if
C
C        WRITE(*,*) ' --> just after left', RXPA(NSTPART_X(IX-1)+1)
C        WRITE(*,*) ' --> just before right',RXPA(NSTPART_X(IX))
C        if (ix.eq.nx) then
C         WRITE(*,*) ' --> just after right',' last particle'
C        else
C         WRITE(*,*) ' --> just after right',RXPA(NSTPART_X(IX)+1)
C        end if
C
C       END DO
C
C       WRITE(*,*)

       RETURN
       END
