*****************************************************************
       SUBROUTINE MERGER(NFILE,NMASA,NNCLUS,DMPITER,
     &                   DMPCLUS,REALCLUS)
*****************************************************************
*      Computes the complete merger tree (used if PLOT=2),
*      i.e., finds all haloes at the previous output which share
*      particles with a given halo at the following output.
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
*      Computes the reduced merger tree (used if PLOT=3),
*      i.e., finds the main 'parent' halo of each halo, at its
*      previous iteration
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
