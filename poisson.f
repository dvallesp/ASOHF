************************************************************************
      SUBROUTINE POISSON(NL,NX,NY,NZ,DX,NPATCH,PARE,PATCHNX,PATCHNY,
     &                   PATCHNZ,PATCHX,PATCHY,PATCHZ,U1,U11,POT,POT1)
************************************************************************

      IMPLICIT NONE

      INCLUDE 'input_files/asohf_parameters.dat'

      INTEGER NL,NX,NY,NZ
      REAL DX
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PARE(NPALEV)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      REAL*4 U1(NMAX,NMAY,NMAZ)
      REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
      REAL*4 POT(NMAX,NMAY,NMAZ)
      REAL*4 POT1(NAMRX,NAMRY,NAMRZ,NPALEV)

*     LOCAL VARIABLES
      REAL KKK(NMAX,NMAY,NMAZ)
      INTEGER I,MAXIT,IX,JY,KZ,N1,N2,N3
      REAL PRECIS

      INTEGER CR3AMR1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1X(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Y(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Z(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)

*     level 0
      CALL MOMENTO(DX,NX,NY,NZ,KKK)
      CALL POFFT3D(NX,NY,NZ,KKK,U1,POT)

      WRITE(*,*) 'Potential throgh i=64,j=64,k=1 to 128'
      WRITE(*,*) (POT(64,64,I),I=1,128)

*     levels > 0
      CALL COMPUTE_CR3AMR(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                    PATCHX,PATCHY,PATCHZ,PARE,CR3AMR1,CR3AMR1X,
     &                    CR3AMR1Y,CR3AMR1Z)

c      WRITE(*,*) PATCHX(100),PATCHY(100),PATCHZ(100)
c      WRITE(*,*) CR3AMR1(1,1,1,100),CR3AMR1X(1,1,1,100),
c     &           CR3AMR1Y(1,1,1,100),CR3AMR1Z(1,1,1,100)
c
c      WRITE(*,*) PARE(401),PATCHX(401),PATCHY(401),PATCHZ(401)
c      WRITE(*,*) CR3AMR1(1,1,1,401),CR3AMR1X(1,1,1,401),
c     &           CR3AMR1Y(1,1,1,401),CR3AMR1Z(1,1,1,401)


      ! THEN WE SHOULD CALL TO POISSON, BUT CAUTION WITH COMMON VARIABLES
      ! (MAKE THEM PARAMETERS) AND WITH DIMENSIONS OF POTS, DENSITIES.
      ! WE WILL HAVE TO MAKE A LOCAL DENSITY AS WELL TO EXTEND THE SOURCE,
      ! SINCE WE DON'T WANT TO CHANGE THE DIMENSIONS OF THE OUTSIDE VARIABLE
      MAXIT=500
      PRECIS=1E-6
      CALL POTAMR(NL,NX,NY,NZ,DX,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &            PATCHX,PATCHY,PATCHZ,U11,POT1,U1,POT,PRECIS,MAXIT,
     &            CR3AMR1,CR3AMR1X,CR3AMR1Y,CR3AMR1Z)


      DO I=1,SUM(NPATCH(0:NL))
       WRITE(*,*) I,
     & MINVAL(POT1(1:PATCHNX(I),1:PATCHNY(I),1:PATCHNZ(I),I)),
     & MAXVAL(POT1(1:PATCHNX(I),1:PATCHNY(I),1:PATCHNZ(I),I))
      END DO

      OPEN(99,FILE='output_files/potential_asohf',STATUS='UNKNOWN',
     &      FORM='UNFORMATTED')
        write(99) (((pot(ix,jy,kz),ix=1,nx),jy=1,ny),kz=1,nz)
        do i=1,sum(npatch(0:nl))
         n1=patchnx(i)
         n2=patchny(i)
         n3=patchnz(i)
         write(99) (((pot1(ix,jy,kz,i),ix=1,n1),jy=1,n2),kz=1,n3)
        end do
      CLOSE(99)


      RETURN
      END

************************************************************************
      SUBROUTINE MOMENTO(DX,NX,NY,NZ,KKK)
************************************************************************
*     Computes the momentum indices (i.e., the Green function) to
*     solve Poisson's equation in Fourier space, and stores them in
*     the array KKK
************************************************************************

      IMPLICIT NONE

      INCLUDE 'input_files/asohf_parameters.dat'

      INTEGER NX,NY,NZ,I,J,K
      real PI, DELTA,CONSTA,KX2,KY2,KZ2,DEL2
      INTEGER KX,KY,KZ,NX2,NY2,NZ2

      real KKK(NMAX,NMAY,NMAZ)
      real DX        !size of coarse grid cell

      DELTA = DX

      PI=ACOS(-1.0)

**    KERNEL SIN
**    CONSTA=(2.0*PI/(DELTA*NX))*(DELTA/2.0)

      CONSTA=PI/NX
      DEL2=(DELTA*0.5)**2

      NX2=NX/2+2
      NY2=NY/2+2
      NZ2=NZ/2+2

**  IN THIS SECTION MOMENTUM SPACE INDECES ARE COMPUTED FOR POISSON EQ.
      DO I=1,NX
        IF (I.GE.NX2) THEN
          KX=I-NX-1
        ELSE
          KX=I-1
        END IF
        KX2=CONSTA*KX
        KX2=(SIN(KX2))**2
      DO J=1,NY
        IF (J.GE.NY2) THEN
          KY=J-NY-1
        ELSE
          KY=J-1
        END IF
        KY2=CONSTA*KY
        KY2=(SIN(KY2))**2
      DO K=1,NZ
        IF (K.GE.NZ2) THEN
          KZ=K-NZ-1
        ELSE
          KZ=K-1
        END IF
        KZ2=CONSTA*KZ
        KZ2=(SIN(KZ2))**2

        KKK(I,J,K)=KX2+KY2+KZ2
      END DO
      END DO
      END DO

**     CENTRAL CONDITION
**     KKK(1,1,1)=0.D0
      KKK(1,1,1)=1.0E30

      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
        KKK(I,J,K)=-DEL2/KKK(I,J,K)
      END DO
      END DO
      END DO


      RETURN
      END

************************************************************************
      SUBROUTINE POFFT3D(NX,NY,NZ,KKK,U1,POT)
************************************************************************
*     Solve Poisson's equation in Fourier space, assuming periodic
*     boundary conditions.
************************************************************************

      IMPLICIT NONE

      INCLUDE 'input_files/asohf_parameters.dat'

      INTEGER NX,NY,NZ,I,J,K

      REAL  U1(NMAX,NMAY,NMAZ)    ! source in poisson equation
      REAL  POT(NMAX,NMAY,NMAZ)   ! field to solve

*     FFT variables
      REAL KKK(NMAX,NMAY,NMAZ)
      REAL DATA1(2*NMAX*NMAY*NMAZ)

      INTEGER I1,I2,IJK,NYZ,NXYZ,NNN(3)
      INTEGER NZ2,NZ3,NX2,NX3,NY2,NY3


      NYZ=NY*NZ
      NXYZ=NX*NY*NZ

**    FFT METHOD
      DO I=1,NX
      DO J=1,NY
      DO K=1,NZ
        IJK=(I-1)*NYZ + (J-1)*NZ + K
        I1= 2*IJK -1
        I2= 2*IJK
        DATA1(I1)=U1(I,J,K)    ! real part
        DATA1(I2)=0.0          ! imaginary part
      END DO
      END DO
      END DO

      NNN(1)=NX
      NNN(2)=NY
      NNN(3)=NZ

      CALL FOURN(DATA1,NNN,3,1)

**    POISSON EQUATION IN MOMENTUM SPACE
      DO I=1,NX
      DO J=1,NY
      DO K=1,NZ
         IJK=(I-1)*NYZ + (J-1)*NZ + K
         I1= 2*IJK -1
         I2= 2*IJK
         DATA1(I1)=DATA1(I1)*KKK(I,J,K)
         DATA1(I2)=DATA1(I2)*KKK(I,J,K)
      END DO
      END DO
      END DO

      CALL FOURN(DATA1,NNN,3,-1)

      DO I=1,NX
      DO J=1,NY
      DO K=1,NZ
         IJK=(I-1)*NYZ + (J-1)*NZ + K
         I1= 2*IJK -1
         POT(I,J,K)=DATA1(I1)/NXYZ
      END DO
      END DO
      END DO

      RETURN
      END

***********************************************************************
      SUBROUTINE FOURN(DATA,NN,NDIM,ISIGN)
***********************************************************************
*     Computes the FFT (in NDIM dimensions). These function has been
*     adapted from the following reference:
*---------------------------------------------------------------------*
*     Ref.: Numerical Recipes in FORTRAN 77: Volume 1                 *
*     W.H. Press, B.P. Flannery, S.A. Teukolsky, W.T. Vetterling      *
*     1992, Cambridge University Press                                *
***********************************************************************

      IMPLICIT real(A-H,O-Z)

      INCLUDE 'input_files/asohf_parameters.dat'

      DOUBLE PRECISION WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION NN(NDIM),DATA(*)
      NTOT=1
      DO 11 IDIM=1,NDIM
        NTOT=NTOT*NN(IDIM)
11    CONTINUE
      NPREV=1
      DO 18 IDIM=1,NDIM
        N=NN(IDIM)
        NREM=NTOT/(N*NPREV)
        IP1=2*NPREV
        IP2=IP1*N
        IP3=IP2*NREM
        I2REV=1
        DO 14 I2=1,IP2,IP1
          IF(I2.LT.I2REV)THEN
            DO 13 I1=I2,I2+IP1-2,2
              DO 12 I3=I1,IP3,IP2
                I3REV=I2REV+I3-I2
                TEMPR=DATA(I3)
                TEMPI=DATA(I3+1)
                DATA(I3)=DATA(I3REV)
                DATA(I3+1)=DATA(I3REV+1)
                DATA(I3REV)=TEMPR
                DATA(I3REV+1)=TEMPI
12            CONTINUE
13          CONTINUE
          ENDIF
          IBIT=IP2/2
1         IF ((IBIT.GE.IP1).AND.(I2REV.GT.IBIT)) THEN
            I2REV=I2REV-IBIT
            IBIT=IBIT/2
          GO TO 1
          ENDIF
          I2REV=I2REV+IBIT
14      CONTINUE
        IFP1=IP1
2       IF(IFP1.LT.IP2)THEN
          IFP2=2*IFP1
          THETA=ISIGN*6.28318530717959D0/(IFP2/IP1)
          WPR=-2.D0*DSIN(0.5D0*THETA)**2
          WPI=DSIN(THETA)
          WR=1.D0
          WI=0.D0
          DO 17 I3=1,IFP1,IP1
            DO 16 I1=I3,I3+IP1-2,2
              DO 15 I2=I1,IP3,IFP2
                K1=I2
                K2=K1+IFP1
                TEMPR=SNGL(WR)*DATA(K2)-SNGL(WI)*DATA(K2+1)
                TEMPI=SNGL(WR)*DATA(K2+1)+SNGL(WI)*DATA(K2)
                DATA(K2)=DATA(K1)-TEMPR
                DATA(K2+1)=DATA(K1+1)-TEMPI
                DATA(K1)=DATA(K1)+TEMPR
                DATA(K1+1)=DATA(K1+1)+TEMPI
15            CONTINUE
16          CONTINUE
            WTEMP=WR
            WR=WR*WPR-WI*WPI+WR
            WI=WI*WPR+WTEMP*WPI+WI
17        CONTINUE
          IFP1=IFP2
        GO TO 2
        ENDIF
        NPREV=N*NPREV
18    CONTINUE
      RETURN
      END

************************************************************************
      SUBROUTINE COMPUTE_CR3AMR(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &                          PATCHNZ,PATCHX,PATCHY,PATCHZ,PARE,
     &                          CR3AMR1,CR3AMR1X,CR3AMR1Y,CR3AMR1Z)
************************************************************************
*     Build the AMR grid: cells center and interface positions for each
*     refinement patch. Also computes the CR3AMR variables, which
*     contain the "parent" cell of a given one which is well-inside its
*     parent patch.
************************************************************************

      IMPLICIT NONE

      INCLUDE 'input_files/asohf_parameters.dat'

      INTEGER NX,NY,NZ,NL1,NL

      INTEGER NPATCH(0:NLEVELS)
      INTEGER PARE(NPALEV)
      INTEGER PATCHNX(NPALEV)
      INTEGER PATCHNY(NPALEV)
      INTEGER PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV)
      INTEGER PATCHY(NPALEV)
      INTEGER PATCHZ(NPALEV)

      INTEGER CR3AMR1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1X(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Y(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Z(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)

      INTEGER I,J,K,IX,JY,KZ,I1,J1,K1,IR,IPALE,BOR
      INTEGER NEF,NCELL,PAX1,PAX2,PAY1,PAY2,PAZ1,PAZ2
      INTEGER IPATCH,II,JJ,KK,N1,N2,N3
      INTEGER INMAX(3)
      INTEGER NBAS,IPA2,NPX,NPY,NPZ
      real DXPA,DYPA,DZPA

      real XX,YY,ZZ,XX1,YY1,ZZ1,XX2,YY2,ZZ2,BAS1,BAS2,BAS3

      INTEGER CONTROL,IP,P1,P2,P3,NP1,NP2,NP3
      INTEGER L1,L2,L3,LL1,LL2,LL3,CR1,CR2,CR3,CR4,CR5,CR6
      INTEGER LN1,LN2,LN3,LNN1,LNN2,LNN3
      INTEGER LOW1,LOW2,LOW3,LOW4
      INTEGER KR1,KR2,KR3,ABUELO,MARCA

      DO IR=NL,2,-1
      LOW1=SUM(NPATCH(0:IR-1))+1
      LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(IR,NL,LOW1,LOW2,PATCHX,PATCHY,PATCHZ,PATCHNX,
!$OMP+                   PATCHNY,PATCHNZ,CR3AMR1,NX,NY,NZ,CR3AMR1X,
!$OMP+                   CR3AMR1Y,CR3AMR1Z,PARE),
!$OMP+            PRIVATE(I,L1,L2,L3,IX,JY,KZ,KR1,KR2,KR3,MARCA,
!$OMP+                    CR1,CR2,CR3,ABUELO),
!$OMP+            DEFAULT(NONE)
      DO I=LOW1,LOW2
       L1=PATCHX(I)
       L2=PATCHY(I)
       L3=PATCHZ(I)
       DO KZ=-2,PATCHNZ(I)+3
       DO JY=-2,PATCHNY(I)+3
       DO IX=-2,PATCHNX(I)+3

       MARCA=0

       CR1=L1-1+INT((IX+1)/2)
       CR2=L2-1+INT((JY+1)/2)
       CR3=L3-1+INT((KZ+1)/2)

       IF (IX.LT.-1) CR1=INT((IX+1)/2)+L1-2
       IF (JY.LT.-1) CR2=INT((JY+1)/2)+L2-2
       IF (KZ.LT.-1) CR3=INT((KZ+1)/2)+L3-2

       KR1=CR1
       KR2=CR2
       KR3=CR3

       ABUELO=PARE(I)

       DO WHILE (MARCA.EQ.0)  !...................
        IF (CR1.LT.2.OR.CR1.GT.PATCHNX(ABUELO)-1.OR.   !----------------
     &      CR2.LT.2.OR.CR2.GT.PATCHNY(ABUELO)-1.OR.
     &      CR3.LT.2.OR.CR3.GT.PATCHNZ(ABUELO)-1) THEN
*        the progenitor cell of the cell ix,jy,kz
*        is in the boundary of its patch. we need to recursively look
*        for a cell not in the boundary

         CR1=PATCHX(ABUELO)-1+INT((KR1+1)/2)
         CR2=PATCHY(ABUELO)-1+INT((KR2+1)/2)
         CR3=PATCHZ(ABUELO)-1+INT((KR3+1)/2)

         IF (KR1.LT.-1) CR1=PATCHX(ABUELO)-2+INT((KR1+1)/2)
         IF (KR2.LT.-1) CR2=PATCHY(ABUELO)-2+INT((KR2+1)/2)
         IF (KR3.LT.-1) CR3=PATCHZ(ABUELO)-2+INT((KR3+1)/2)

         KR1=CR1
         KR2=CR2
         KR3=CR3

         ABUELO=PARE(ABUELO)

         IF (ABUELO.EQ.0) THEN
          MARCA=1
          CR3AMR1(IX,JY,KZ,I)=0
          CR3AMR1X(IX,JY,KZ,I)=KR1
          CR3AMR1Y(IX,JY,KZ,I)=KR2
          CR3AMR1Z(IX,JY,KZ,I)=KR3
         END IF
        ELSE                                           !----------------
         CR3AMR1(IX,JY,KZ,I)=ABUELO
         CR3AMR1X(IX,JY,KZ,I)=KR1
         CR3AMR1Y(IX,JY,KZ,I)=KR2
         CR3AMR1Z(IX,JY,KZ,I)=KR3
*        this cell is well inside its parent patch!!
         MARCA=1                                       !----------------
        END IF
       END DO                 !...................

       END DO
       END DO
       END DO
      END DO
      END DO


      IR=1
!$OMP PARALLEL DO SHARED(IR,NL,PATCHX,PATCHY,PATCHZ,PATCHNX,PATCHNY,
!$OMP+                   PATCHNZ,CR3AMR1,NX,NY,NZ,CR3AMR1X,CR3AMR1Y,
!$OMP+                   CR3AMR1Z,PARE),
!$OMP+            PRIVATE(I,L1,L2,L3,IX,JY,KZ,KR1,KR2,KR3,CR1,CR2,CR3),
!$OMP+            DEFAULT(NONE)
      DO I=1,NPATCH(IR)

       L1=PATCHX(I)
       L2=PATCHY(I)
       L3=PATCHZ(I)

       DO KZ=-2,PATCHNZ(I)+3
       DO JY=-2,PATCHNY(I)+3
       DO IX=-2,PATCHNX(I)+3

         CR1=L1-1+INT((IX+1)/2)
         CR2=L2-1+INT((JY+1)/2)
         CR3=L3-1+INT((KZ+1)/2)

         IF (IX.LT.-1) CR1=INT((IX+1)/2)+L1-2
         IF (JY.LT.-1) CR2=INT((JY+1)/2)+L2-2
         IF (KZ.LT.-1) CR3=INT((KZ+1)/2)+L3-2

         KR1=CR1
         KR2=CR2
         KR3=CR3

         CR3AMR1(IX,JY,KZ,I)=0
         CR3AMR1X(IX,JY,KZ,I)=KR1
         CR3AMR1Y(IX,JY,KZ,I)=KR2
         CR3AMR1Z(IX,JY,KZ,I)=KR3

       END DO
       END DO
       END DO
      END DO

      RETURN
      END

************************************************************************
      SUBROUTINE POTAMR(NL,NX,NY,NZ,DX,NPATCH,PARE,PATCHNX,PATCHNY,
     &                  PATCHNZ,PATCHX,PATCHY,PATCHZ,AU11,APOT1,
     &                  U1,POT,PRECIS,MAXIT,CR3AMR1,CR3AMR1X,CR3AMR1Y,
     &                  CR3AMR1Z)
************************************************************************
*     Solves Poisson equation for the refinement patches, taking into
*     account the boundary conditions imposed by the coarser cells.
*     It uses a SOR method with Chebyshev acceleration procedure to
*     set the overrelaxation parameter. It uses 3 ficticious cell at
*     each boundary. Boundary conditions are enforced in the ficticious
*     boundary (cells -2 and N+3).
************************************************************************

      IMPLICIT NONE

      INCLUDE 'input_files/asohf_parameters.dat'

      INTEGER N1,N2,N3,NP1,NP2,NP3,L1,L2,L3
      INTEGER I,J,K,IX,JY,KZ,I1,J2,K3,II,IR,BOR
      INTEGER NX,NY,NZ,NL

      REAL  AU11(NAMRX,NAMRY,NAMRZ,NPALEV)
      REAL  APOT1(NAMRX,NAMRY,NAMRZ,NPALEV)

      REAL  U1(NMAX,NMAY,NMAZ)
      REAL  POT(NMAX,NMAZ,NMAZ)

      INTEGER NPATCH(0:NLEVELS)
      INTEGER PARE(NPALEV)
      INTEGER PATCHNX(NPALEV)
      INTEGER PATCHNY(NPALEV)
      INTEGER PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV)
      INTEGER PATCHY(NPALEV)
      INTEGER PATCHZ(NPALEV)

*     Here pot1 and u1 are local !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL POT1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3)
      REAL U11(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3)
*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real BAS,ERROR,ERRMAX,BASS,PI
      real SSS,DXPA,DX,WWW,ERR,ERRTOT
      real RADIUS,SNOR,RESNOR,PRECIS
      INTEGER CR1,CR2,CR3,MARCA,LOW1,LOW2,MAXIT
      real AAA,BBB,CCC

      REAL*4  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
      COMMON /GRID/ RADX,RADY,RADZ

      REAL*4  ARX(0:NAMRX+1,NPALEV),ARY(0:NAMRX+1,NPALEV),
     &        ARZ(0:NAMRX+1,NPALEV)
      COMMON /GRIDAMR/ ARX,ARY,ARZ

*     Here RX,RY,RZ are local !!!!!!!!!!!!!!!!!!!!!!!!!!!1
      REAL*4  RX(-2:NAMRX+3,NPALEV),RY(-2:NAMRX+3,NPALEV),
     &        RZ(-2:NAMRX+3,NPALEV)
*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

      real UBAS(3,3,3),RXBAS(3),RYBAS(3),RZBAS(3),FUIN
      INTEGER KARE,KR1,KR2,KR3

      INTEGER CR3AMR1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1X(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Y(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Z(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
*     ---------------------------------------------------------------------------

      PI=ACOS(-1.0)
      BOR=2

      DO IR=1,NL
       DXPA=DX/(2.0**IR)
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)

        DO IX=0,N1+1
         RX(IX,I)=ARX(IX,I)
        END DO
        BAS=RX(0,I)
        RX(-1,I)=BAS-DXPA
        RX(-2,I)=BAS-2.0*DXPA
        BAS=RX(N1+1,I)
        RX(N1+2,I)=BAS+DXPA
        RX(N1+3,I)=BAS+2.0*DXPA

        DO JY=0,N2+1
         RY(JY,I)=ARY(JY,I)
        END DO
        BAS=RY(0,I)
        RY(-1,I)=BAS-DXPA
        RY(-2,I)=BAS-2.0*DXPA
        BAS=RY(N2+1,I)
        RY(N2+2,I)=BAS+DXPA
        RY(N2+3,I)=BAS+2.0*DXPA

        DO KZ=0,N3+1
         RZ(KZ,I)=ARZ(KZ,I)
        END DO
        BAS=RZ(0,I)
        RZ(-1,I)=BAS-DXPA
        RZ(-2,I)=BAS-2.0*DXPA
        BAS=RZ(N3+1,I)
        RZ(N3+2,I)=BAS+DXPA
        RZ(N3+3,I)=BAS+2.0*DXPA
       END DO
      END DO

      IR=1
      DXPA=DX/(2.0**IR)

!$OMP PARALLEL DO SHARED(IR,DX,DXPA,POT,APOT1,NX,NY,NZ,NPATCH,
!$OMP+         PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
!$OMP+         AU11,PI,PRECIS,CR3AMR1X,CR3AMR1Y,CR3AMR1Z,U1,MAXIT,BOR),
!$OMP+      PRIVATE(I,N1,N2,N3,L1,L2,L3,NP1,NP2,NP3,JY,KZ,IX,
!$OMP+         I1,J2,K3,II,SSS,BAS,BASS,SNOR,RESNOR,ERR,ERRTOT,
!$OMP+         WWW,RADIUS,ERROR,ERRMAX,CR1,CR2,CR3,POT1,MARCA,
!$OMP+         UBAS,FUIN,KR1,KR2,KR3,U11),
!$OMP+      DEFAULT(NONE)
      DO I=1,NPATCH(IR)

       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       L1=PATCHX(I)
       L2=PATCHY(I)
       L3=PATCHZ(I)

       NP1=NX
       NP2=NY
       NP3=NZ

       ! initialize the potential (real and ficticious cells) by interpolation
       DO KZ=-2,N3+3
       DO JY=-2,N2+3
       DO IX=-2,N1+3
        KR1=CR3AMR1X(IX,JY,KZ,I)
        KR2=CR3AMR1Y(IX,JY,KZ,I)
        KR3=CR3AMR1Z(IX,JY,KZ,I)

        UBAS(1:3,1:3,1:3)=POT(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
        CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
        POT1(IX,JY,KZ)=FUIN
       END DO
       END DO
       END DO

       ! extend the source (ficticious cells) by interpolation
       DO KZ=1-BOR,N3+BOR
       DO JY=1-BOR,N2+BOR
       DO IX=1-BOR,N1+BOR
        IF(IX.LT.1.OR.IX.GT.N1.OR.
     &     JY.LT.1.OR.jy.GT.N2.OR.
     &     KZ.LT.1.OR.kz.GT.N3) THEN
         KR1=CR3AMR1X(IX,JY,KZ,I)
         KR2=CR3AMR1Y(IX,JY,KZ,I)
         KR3=CR3AMR1Z(IX,JY,KZ,I)

         UBAS(1:3,1:3,1:3)=U1(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
         CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
         U11(IX,JY,KZ)=FUIN
        ELSE
         U11(IX,JY,KZ)=AU11(IX,JY,KZ,I)
        END IF

       END DO
       END DO
       END DO

       SNOR=0.0
       DO KZ=1-BOR,N3+BOR
       DO JY=1-BOR,N2+BOR
       DO IX=1-BOR,N1+BOR
         SSS=DXPA*DXPA*U11(IX,JY,KZ)
         SNOR=SNOR+ABS(SSS)
       END DO
       END DO
       END DO


       RADIUS=COS(PI/((N1+N2+N3+12)/3.0))
       WWW=1.0
       RESNOR=SNOR
       II=0
       MARCA=0
       DO WHILE (MARCA.EQ.0.OR.II.LT.2) ! SOR iteration
        II=II+1
        RESNOR=0.0

        ERRTOT=-1.0
        DO KZ=1-BOR,N3+BOR
        DO JY=1-BOR,N2+BOR
        DO IX=1-BOR,N1+BOR
         IF (MOD((IX+JY+KZ),2).EQ.MOD((II+1),2)) THEN
*         POISSON EQUATION
          SSS=DXPA*DXPA*U11(IX,JY,KZ)

          ERR=POT1(IX,JY,KZ)
          BAS=POT1(IX+1,JY,KZ)+POT1(IX-1,JY,KZ)+POT1(IX,JY+1,KZ)
     &       +POT1(IX,JY-1,KZ)+POT1(IX,JY,KZ+1)+POT1(IX,JY,KZ-1)
     &       -6.0*POT1(IX,JY,KZ) - SSS

          RESNOR=RESNOR+ABS(BAS)

          POT1(IX,JY,KZ)=POT1(IX,JY,KZ)+WWW*BAS/6.0
          IF (ERR.NE.0.0) ERR=POT1(IX,JY,KZ)/ERR - 1.0
          ERRTOT=MAX(ERRTOT,ABS(ERR))
         END IF
        END DO
        END DO
        END DO
        IF (RESNOR.LT.(PRECIS*SNOR)) MARCA=1
*
        IF (ERRTOT.LE.0.1*PRECIS.AND.II.GT.2) MARCA=1

        WWW=1.0/(1.0-0.25*WWW*RADIUS**2)
        IF (II.EQ.1) WWW=1.0/(1.0-0.5*RADIUS**2)
*
        IF (II.GT.MAXIT) MARCA=1
       END DO

       APOT1(1:N1,1:N2,1:N3,I)=POT1(1:N1,1:N2,1:N3)
      END DO


*     OTHER LEVELS
      DO IR=2,NL

      DXPA=DX/(2**IR)

      LOW1=SUM(NPATCH(0:IR-1))+1
      LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(IR,DX,DXPA,POT,APOT1,NX,NY,NZ,NPATCH,
!$OMP+        PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
!$OMP+        AU11,PI,PRECIS,LOW1,LOW2,
!$OMP+        CR3AMR1,CR3AMR1X,CR3AMR1Y,CR3AMR1Z,PARE,
!$OMP+        RX,RY,RZ,RADX,RADY,RADZ,U1,MAXIT,BOR),
!$OMP+     PRIVATE(I,N1,N2,N3,L1,L2,L3,NP1,NP2,NP3,JY,KZ,IX,MARCA,
!$OMP+        I1,J2,K3,II,SSS,BAS,BASS,RESNOR,SNOR,ERR,ERRTOT,
!$OMP+        WWW,RADIUS,ERROR,ERRMAX,CR1,CR2,CR3,POT1,
!$OMP+        KARE,KR1,KR2,KR3,UBAS,RXBAS,RYBAS,RZBAS,FUIN,
!$OMP+        AAA,BBB,CCC,U11),
!$OMP+     DEFAULT(NONE)
      DO I=LOW1,LOW2
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       L1=PATCHX(I)
       L2=PATCHY(I)
       L3=PATCHZ(I)

       NP1=PATCHNX(PARE(I))
       NP2=PATCHNY(PARE(I))
       NP3=PATCHNZ(PARE(I))


       DO KZ=-2,N3+3
       DO JY=-2,N2+3
       DO IX=-2,N1+3

        KARE=CR3AMR1(IX,JY,KZ,I)
        KR1=CR3AMR1X(IX,JY,KZ,I)
        KR2=CR3AMR1Y(IX,JY,KZ,I)
        KR3=CR3AMR1Z(IX,JY,KZ,I)

        AAA=RX(IX,I)
        BBB=RY(JY,I)
        CCC=RZ(KZ,I)

        IF (KARE.GT.0) THEN
          UBAS(1:3,1:3,1:3)=
     &        APOT1(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)

          RXBAS(1:3)=RX(KR1-1:KR1+1,KARE)
          RYBAS(1:3)=RY(KR2-1:KR2+1,KARE)
          RZBAS(1:3)=RZ(KR3-1:KR3+1,KARE)

          CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                   RXBAS,RYBAS,RZBAS,UBAS,FUIN)
          POT1(IX,JY,KZ)=FUIN
        ELSE
          UBAS(1:3,1:3,1:3)=
     &        POT(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)

          RXBAS(1:3)=RADX(KR1-1:KR1+1)
          RYBAS(1:3)=RADY(KR2-1:KR2+1)
          RZBAS(1:3)=RADZ(KR3-1:KR3+1)

          CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                   RXBAS,RYBAS,RZBAS,UBAS,FUIN)
          POT1(IX,JY,KZ)=FUIN
        ENDIF

       END DO
       END DO
       END DO

       DO KZ=1-BOR,N3+BOR
       DO JY=1-BOR,N2+BOR
       DO IX=1-BOR,N1+BOR

        IF(IX.LT.1.OR.IX.GT.N1.OR.
     &     JY.LT.1.OR.jy.GT.N2.OR.
     &     KZ.LT.1.OR.kz.GT.N3) THEN

          KARE=CR3AMR1(IX,JY,KZ,I)
          KR1=CR3AMR1X(IX,JY,KZ,I)
          KR2=CR3AMR1Y(IX,JY,KZ,I)
          KR3=CR3AMR1Z(IX,JY,KZ,I)

          AAA=RX(IX,I)
          BBB=RY(JY,I)
          CCC=RZ(KZ,I)

          IF (KARE.GT.0) THEN
            UBAS(1:3,1:3,1:3)=
     &        AU11(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)

            RXBAS(1:3)=RX(KR1-1:KR1+1,KARE)
            RYBAS(1:3)=RY(KR2-1:KR2+1,KARE)
            RZBAS(1:3)=RZ(KR3-1:KR3+1,KARE)

            CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                   RXBAS,RYBAS,RZBAS,UBAS,FUIN)
            U11(IX,JY,KZ)=FUIN
          ELSE
            UBAS(1:3,1:3,1:3)=
     &        U1(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)

            RXBAS(1:3)=RADX(KR1-1:KR1+1)
            RYBAS(1:3)=RADY(KR2-1:KR2+1)
            RZBAS(1:3)=RADZ(KR3-1:KR3+1)

            CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                   RXBAS,RYBAS,RZBAS,UBAS,FUIN)
            U11(IX,JY,KZ)=FUIN
          ENDIF
        ELSE
         U11(IX,JY,KZ)=AU11(IX,JY,KZ,I)
        END IF

       END DO
       END DO
       END DO


       SNOR=0.0
       DO KZ=1-BOR,N3+BOR
       DO JY=1-BOR,N2+BOR
       DO IX=1-BOR,N1+BOR
         SSS=DXPA*DXPA*U11(IX,JY,KZ)
         SNOR=SNOR+ABS(SSS)
       END DO
       END DO
       END DO


       RADIUS=COS(PI/((N1+N2+N3+12)/3.0))
*
       WWW=1.0
       II=0
       RESNOR=SNOR
       MARCA=0
       DO WHILE (MARCA.EQ.0.OR.II.LT.2)
        II=II+1
        RESNOR=0.0

        ERRTOT=-1.0
        DO KZ=1-BOR,N3+BOR
        DO JY=1-BOR,N2+BOR
        DO IX=1-BOR,N1+BOR
         IF (MOD((IX+JY+KZ),2).EQ.MOD((II+1),2)) THEN
*         POISSON EQUATION
          SSS=DXPA*DXPA*U11(IX,JY,KZ)

          ERR=POT1(IX,JY,KZ)
          BAS=POT1(IX+1,JY,KZ)+POT1(IX-1,JY,KZ)+POT1(IX,JY+1,KZ)
     &       +POT1(IX,JY-1,KZ)+POT1(IX,JY,KZ+1)+POT1(IX,JY,KZ-1)
     &       -6.0*POT1(IX,JY,KZ) - SSS

          RESNOR=RESNOR+ABS(BAS)
          POT1(IX,JY,KZ)=POT1(IX,JY,KZ)+WWW*BAS/6.0
          IF (ERR.NE.0.0) ERR=POT1(IX,JY,KZ)/ERR - 1.0
          ERRTOT=MAX(ERRTOT,ABS(ERR))
         END IF
        END DO
        END DO
        END DO
        IF (RESNOR.LT.(PRECIS*SNOR)) MARCA=1
*
        IF (ERRTOT.LE.0.1*PRECIS.AND.II.GT.2) MARCA=1
        WWW=1.0/(1.0-0.25*WWW*RADIUS**2)
        IF (II.EQ.1) WWW=1.0/(1.0-0.5*RADIUS**2)
*
        IF (II.GT.MAXIT) MARCA=1
       END DO

       APOT1(1:N1,1:N2,1:N3,I)=POT1(1:N1,1:N2,1:N3)

      END DO
      END DO

      RETURN
      END


***********************************************************************
      SUBROUTINE LININT52D_NEW(IX,JY,KZ,U,FUIN)
***********************************************************************
*     Linear interpolation of a l=1 cell from the base grid
************************************************************************

      IMPLICIT NONE

      INCLUDE 'input_files/asohf_parameters.dat'

      INTEGER II,JJ,KK,IX,JY,KZ
      real XX,YY,ZZ,SIGNO
      real FUIN, U(3,3,3)

      real DXX,DYY,DZZ,DXMAS,DYMAS,DZMAS,DXMIN,DYMIN,DZMIN
      real LIM,LIMA
      real DXCEN,DYCEN,DZCEN

      LIM=8.0

      II=-1
      JJ=-1
      KK=-1
      IF (MOD(IX,2).EQ.0) II=1
      IF (MOD(JY,2).EQ.0) JJ=1
      IF (MOD(KZ,2).EQ.0) KK=1

      DXX=0.0
      DYY=0.0
      DZZ=0.0

*     DU/DX
      DXMAS=U(3,2,2)-U(2,2,2)
      DXMIN=U(2,2,2)-U(1,2,2)
      DXCEN=0.5*(DXMAS+DXMIN)
      LIMA=ABS(U(3,2,2)-U(1,2,2))/
     &     MAX(1.E-30,ABS(MIN(U(3,2,2),U(1,2,2))))

      IF ((DXMIN*DXMAS).GT.0.0) THEN
       DXX=MIN(ABS(DXCEN),ABS(DXMIN),ABS(DXMAS))
       SIGNO=1.0
       IF (DXCEN.LT.0.0) SIGNO=-1.0
       DXX=DXX*SIGNO
      ELSE
       DXX=0.0
      END IF
      IF (LIMA.GT.LIM) DXX=0.0

*     DU/DY
      DYMAS=U(2,3,2)-U(2,2,2)
      DYMIN=U(2,2,2)-U(2,1,2)
      DYCEN=0.5*(DYMAS+DYMIN)
      LIMA=ABS(U(2,3,2)-U(2,1,2))/
     &     MAX(1.E-30,ABS(MIN(U(2,3,2),U(2,1,2))))

      IF ((DYMIN*DYMAS).GT.0.0) THEN
       DYY=MIN(ABS(DYCEN),ABS(DYMIN),ABS(DYMAS))
       SIGNO=1.0
       IF (DYCEN.LT.0.0) SIGNO=-1.0
       DYY=DYY*SIGNO
      ELSE
       DYY=0.0
      END IF
      IF (LIMA.GT.LIM) DYY=0.0

*     DU/DZ
      DZMAS=U(2,2,3)-U(2,2,2)
      DZMIN=U(2,2,2)-U(2,2,1)
      DZCEN=0.5*(DZMAS+DZMIN)
      LIMA=ABS(U(2,2,3)-U(2,2,1))/
     &     MAX(1.E-30,ABS(MIN(U(2,2,3),U(2,2,1))))

      IF ((DZMIN*DZMAS).GT.0.0) THEN
       DZZ=MIN(ABS(DZCEN),ABS(DZMIN),ABS(DZMAS))
       SIGNO=1.0
       IF (DZCEN.LT.0.0) SIGNO=-1.0
       DZZ=DZZ*SIGNO
      ELSE
       DZZ=0.0
      END IF
      IF (LIMA.GT.LIM) DZZ=0.0

      XX=0.25
      YY=0.25
      ZZ=0.25
      IF (II.LT.0) XX=-0.25
      IF (JJ.LT.0) YY=-0.25
      IF (KK.LT.0) ZZ=-0.25

      FUIN=U(2,2,2) + XX*DXX + YY*DYY + ZZ*DZZ

      RETURN
      END

***********************************************************************
      SUBROUTINE LININT52D_NEW_REAL(XX,YY,ZZ,RX,RY,RZ,U,FUIN)
***********************************************************************
*     Linear interpolation from refined levels
************************************************************************

      IMPLICIT NONE

      INCLUDE 'input_files/asohf_parameters.dat'

      INTEGER II,JJ,KK,IX,JY,KZ
      real XX,YY,ZZ,XXX,YYY,ZZZ,SIGNO
      real FUIN
      real U(3,3,3),RX(3),RY(3),RZ(3)

      real DXX,DYY,DZZ,DXMAS,DYMAS,DZMAS,DXMIN,DYMIN,DZMIN
      real LIM,LIMA
      real DXCEN,DYCEN,DZCEN

      LIM=8.0

      DXX=0.0
      DYY=0.0
      DZZ=0.0

*     DU/DX
      DXMAS=(U(3,2,2)-U(2,2,2))/(RX(3)-RX(2))
      DXMIN=(U(2,2,2)-U(1,2,2))/(RX(2)-RX(1))
      DXCEN=0.5*(DXMAS+DXMIN)

      LIMA=ABS(U(3,2,2)-U(1,2,2))/
     &     MAX(1.E-30,ABS(MIN(U(3,2,2),U(1,2,2))))

      IF ((DXMIN*DXMAS).GT.0.0) THEN
       DXX = ABS(DXCEN)
       SIGNO=1.0
       IF (DXCEN.LT.0.0) SIGNO=-1.0
       DXX=DXX*SIGNO
      ELSE
       DXX=0.0
      END IF
      IF (LIMA.GT.LIM) DXX=0.0

*     DU/DY
      DYMAS=(U(2,3,2)-U(2,2,2))/(RY(3)-RY(2))
      DYMIN=(U(2,2,2)-U(2,1,2))/(RY(2)-RY(1))
      DYCEN=0.5*(DYMAS+DYMIN)

      LIMA=ABS(U(2,3,2)-U(2,1,2))/
     &     MAX(1.E-30,ABS(MIN(U(2,3,2),U(2,1,2))))

      IF ((DYMIN*DYMAS).GT.0.0) THEN
       DYY = ABS(DYCEN)
       SIGNO=1.0
       IF (DYCEN.LT.0.0) SIGNO=-1.0
       DYY=DYY*SIGNO
      ELSE
       DYY=0.0
      END IF
      IF (LIMA.GT.LIM) DYY=0.0

*     DU/DZ
      DZMAS=(U(2,2,3)-U(2,2,2))/(RZ(3)-RZ(2))
      DZMIN=(U(2,2,2)-U(2,2,1))/(RZ(2)-RZ(1))
      DZCEN=0.5*(DZMAS+DZMIN)

      LIMA=ABS(U(2,2,3)-U(2,2,1))/
     &     MAX(1.E-30,ABS(MIN(U(2,2,3),U(2,2,1))))

      IF ((DZMIN*DZMAS).GT.0.0) THEN
       DZZ = ABS(DZCEN)
       SIGNO=1.0
       IF (DZCEN.LT.0.0) SIGNO=-1.0
       DZZ=DZZ*SIGNO
      ELSE
       DZZ=0.0
      END IF
      IF (LIMA.GT.LIM) DZZ=0.0

      XXX=XX-RX(2)
      YYY=YY-RY(2)
      ZZZ=ZZ-RZ(2)

      FUIN=U(2,2,2) + XXX*DXX + YYY*DYY + ZZZ*DZZ

      RETURN
      END
