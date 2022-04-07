**********************************************************************
       SUBROUTINE DIAGONALISE(INPUT_MATRIX,EIGENVALUES_OUT)
**********************************************************************
*      Diagonalise using the Jacobi rotations method
**********************************************************************

       IMPLICIT NONE
       INTEGER N
       REAL*4  INPUT_MATRIX(1:3,1:3),EIGENVALUES_OUT(1:3)
       REAL*8 MATRIX(1:3,1:3),EIGEN_LOC(1:3)
       INTEGER III,JJJ,I,J
       REAL*8 AUX_ARR_1(100),AUX_ARR_2(100)
       REAL*8 COSINUS,VAR1,VAR2,SINUS,SUMA
       REAL*8 TANGENT,TANHALF,ROT_ANG,UPPER_LIM,SUM_ELEMENTS

       N=3

       ! WE'LL WORK WITH A, INPUT_MATRIX IS UNMODIFIED
       DO III=1,N
       DO JJJ=1,N
        MATRIX(III,JJJ)=INPUT_MATRIX(III,JJJ)
       END DO
       END DO

       DO III=1,N
         AUX_ARR_1(III)=MATRIX(III,III) ! diagonal
         EIGEN_LOC(III)=AUX_ARR_1(III) ! this will contain EIGENVALUES_OUT
         AUX_ARR_2(III)=0.d0
       END DO

!      SUM OF ELEMENTS, TO CONTROL ERRORS
       SUM_ELEMENTS=0.0
       DO III=1,N
       DO JJJ=III,N
        SUM_ELEMENTS=SUM_ELEMENTS+ABS(MATRIX(III,JJJ))
       END DO
       END DO

!      JACOBI LOOP, 100 ITERATIONS MAX
       DO I=1,100
        SUMA=0.D0
        DO III=1,N-1
        DO JJJ=III+1,N
         SUMA=SUMA+ABS(MATRIX(III,JJJ))
        END DO
        END DO
        IF (SUMA.LT.1.e-4*SUM_ELEMENTS) EXIT ! WE HAVE CONVERGED

        IF(I.LT.4) THEN
           UPPER_LIM=0.2D0*SUMA**2
        ELSE
           UPPER_LIM=0.D0
        END IF
        DO III=1,N-1
         DO JJJ=III+1, N
          VAR1=100.D0*ABS(MATRIX(III,JJJ))

          IF ((I.GT.4).AND.
     &        (ABS(EIGEN_LOC(III))+VAR1.EQ.ABS(EIGEN_LOC(III))).AND.
     &        (ABS(EIGEN_LOC(JJJ))+VAR1.EQ.ABS(EIGEN_LOC(JJJ)))) THEN

           MATRIX(III,JJJ)=0.D0
          ELSE IF (ABS(MATRIX(III,JJJ)).GT.UPPER_LIM) THEN
           VAR2=EIGEN_LOC(JJJ)-EIGEN_LOC(III)
           IF (ABS(VAR2)+VAR1.EQ.ABS(VAR2)) THEN
            TANGENT=MATRIX(III,JJJ)/VAR2
           ELSE
            ROT_ANG=0.5D0*VAR2/MATRIX(III,JJJ)
            TANGENT=1.D0/(ABS(ROT_ANG)+SQRT(1.D0+ROT_ANG**2))
            IF (ROT_ANG.lt.0.D0) TANGENT=-TANGENT
           END IF
           COSINUS=1.D0/SQRT(1.D0+TANGENT**2)
           SINUS=TANGENT*COSINUS
           TANHALF=SINUS/(1.D0+COSINUS)
           VAR2=TANGENT*MATRIX(III,JJJ)
           AUX_ARR_2(III)=AUX_ARR_2(III)-VAR2
           AUX_ARR_2(JJJ)=AUX_ARR_2(JJJ)+VAR2
           EIGEN_LOC(III)=EIGEN_LOC(III)-VAR2
           EIGEN_LOC(JJJ)=EIGEN_LOC(JJJ)+VAR2
           MATRIX(III,JJJ)=0.D0
           DO J=1,III-1
            VAR1=MATRIX(J,III)
            VAR2=MATRIX(J,JJJ)
            MATRIX(J,III)=VAR1-SINUS*(VAR2+VAR1*TANHALF)
            MATRIX(J,JJJ)=VAR2+SINUS*(VAR1-VAR2*TANHALF)
           END DO
           DO J=III+1,JJJ-1
            VAR1=MATRIX(III,J)
            VAR2=MATRIX(J,JJJ)
            MATRIX(III,J)=VAR1-SINUS*(VAR2+VAR1*TANHALF)
            MATRIX(J,JJJ)=VAR2+SINUS*(VAR1-VAR2*TANHALF)
           END DO
           DO J=JJJ+1,N
            VAR1=MATRIX(III,J)
            VAR2=MATRIX(JJJ,J)
            MATRIX(III,J)=VAR1-SINUS*(VAR2+VAR1*TANHALF)
            MATRIX(JJJ,J)=VAR2+SINUS*(VAR1-VAR2*TANHALF)
          END DO

          END IF  !(I.GT.4)
         END DO !JJJ=III+1,N
        END DO !III=1,N-1
        DO III=1,N
         AUX_ARR_1(III)=AUX_ARR_1(III)+AUX_ARR_2(III)
         EIGEN_LOC(III)=AUX_ARR_1(III)
         AUX_ARR_2(III)=0.D0
        END DO
       END DO !I=1,100

!       MOVE THE EIGENVALUES_OUT TO THE OUTPUT VARIABLE
       DO III=1,N
        EIGENVALUES_OUT(III)=EIGEN_LOC(III)
       END DO

        RETURN
       END

**********************************************************************
       SUBROUTINE SORT_EIGEN(EIGENVALUES,N)
**********************************************************************
*      Sorts the eigenvalues in increasing order
**********************************************************************
       IMPLICIT NONE

       INTEGER N
       REAL*4  EIGENVALUES(N)

       INTEGER II,JJ,KK
       REAL*4 VALOR

       DO II=1,N-1
        KK=II
        VALOR=EIGENVALUES(II)
        DO JJ=II+1,N
         IF(EIGENVALUES(JJ).GE.VALOR) THEN
          KK=JJ
          VALOR=EIGENVALUES(JJ)
         END IF
        END DO
        IF(KK.NE.II) THEN
         EIGENVALUES(KK)=EIGENVALUES(II)
         EIGENVALUES(II)=VALOR
        END IF
       END DO

       RETURN
       END

**********************************************************************
       SUBROUTINE ARGSORT(N,ARRAY,INDX)
**********************************************************************
*      Return the indices which make ARRAY sorted in increasing order.
*      Uses MERGE SORT
*      Rewritten from the idea of the algorithm (LGPL-3.0 License) in
*       https://github.com/Astrokiwi/simple_fortran_argsort/blob/master/LICENSE
**********************************************************************
       INTEGER N
       REAL*4 ARRAY(N)
       INTEGER INDX(N)

       INTEGER TEMPINT(N)

       INTEGER DELTA
       INTEGER I,J,III,K,KKK

       DO I=1,N
        INDX(I)=I
       END DO

       IF (N.EQ.1) RETURN

       DELTA=1
       DO WHILE (DELTA.LT.N)
        DO III=1,N-DELTA,DELTA*2
         I1=III
         I2=III+DELTA
         KKK=MIN(DELTA*2,N-III+1)
         K=1

         DO WHILE (I1.LT.III+DELTA.AND.I2.LT.III+KKK)
          IF (ARRAY(INDX(I1)).LT.ARRAY(INDX(I2))) THEN
           TEMPINT(K)=INDX(I1)
           I1=I1+1
           K=K+1
          ELSE
           TEMPINT(K)=INDX(I2)
           I2=I2+1
           K=K+1
          END IF
         END DO

         IF (I1.LT.III+DELTA) THEN
          TEMPINT(K:KKK)=INDX(I1:III+DELTA-1)
         ELSE
          TEMPINT(K:KKK)=INDX(I2:III+KKK-1)
         ENDIf
          INDX(III:III+KKK-1)=TEMPINT(1:KKK)
         END DO
         DELTA=DELTA*2
       END DO

       RETURN
      END SUBROUTINE
