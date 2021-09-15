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
