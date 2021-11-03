********************************************************************
       SUBROUTINE BRYAN_NORMAN_98(CONTRASTEC,OMEGAZ,OMEGA0,ZETA)
****************************************************************
*      VIRIAL CONTRAST (in terms of the BACKGROUND MATTER density)
*      (Bryan & Norman ApJ, 1998)
*      Delta_vir,c = 18*pi^2 + 82 x - 39 x^2; x=Omega_m(z)-1
****************************************************************
        IMPLICIT NONE

        REAL CONTRASTEC,OMEGAZ ! intent:out
        REAL OMEGA0,ZETA ! intent in
        REAL CONTRASTEX,PI,OMEGALAMBDA0,BAS

        PI=DACOS(-1.D0)
        OMEGALAMBDA0=1.0-OMEGA0
        BAS=OMEGA0*(1+ZETA)**3
        OMEGAZ=BAS/(BAS+OMEGALAMBDA0)
        CONTRASTEX=OMEGAZ - 1.0
        CONTRASTEC= 18.0*PI**2 + 82.0*CONTRASTEX - 39.0*CONTRASTEX**2
        CONTRASTEC=CONTRASTEC/OMEGAZ

        RETURN

       END

********************************************************************
       SUBROUTINE oldOVERLAPING(IFI,IR,NL,REF,ESP,BOUND,CONTA,
     &                       CONTRASTEC,NSHELL,RODO,NPATCH,PATCHNX,
     &                       PATCHNY,PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,NX,
     &                       NY,NZ,NCLUS,MASA,RADIO,CLUSRX,CLUSRY,
     &                       CLUSRZ,REALCLUS,NSOLAP,SOLAPA,NHALLEV)
********************************************************************
*      Accounts for the overlaps between haloes found within the
*      grid (prior to refining with particles)
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

********************************************************************
       SUBROUTINE OVERLAPPING(IFI,IR,NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                        PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                        PATCHRY,PATCHRZ,PARE,NCLUS,MASA,RADIO,
     &                        CLUSRX,CLUSRY,CLUSRZ,REALCLUS,LEVHAL,
     &                        NHALLEV,BOUND,CONTRASTEC,RODO,SOLAP,
     &                        VECINO,NVECI,CR0AMR,CR0AMR11,
     &                        VOL_SOLAP_LOW)
********************************************************************
*      Detect and correct overlaps on the cluster catalogue
********************************************************************

       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

*      I/O DATA
       INTEGER IFI,IR,NL,NX,NY,NZ
       INTEGER NPATCH(0:NLEVELS)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
       REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
       INTEGER PARE(NPALEV)
       INTEGER NCLUS
       REAL MASA(MAXNCLUS),RADIO(MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
       INTEGER REALCLUS(MAXITER,MAXNCLUS),LEVHAL(MAXNCLUS)
       INTEGER NHALLEV(0:NLEVELS)
       REAL BOUND,CONTRASTEC,RODO
       INTEGER SOLAP(NAMRX,NAMRY,NAMRZ,NPALEV)
       INTEGER VECINO(NPALEV,NPALEV),NVECI(NPALEV)
       INTEGER CR0AMR(NMAX,NMAY,NMAZ)
       INTEGER CR0AMR11(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL VOL_SOLAP_LOW

*      GLOBAL VARIABLES
       REAL*4 DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       REAL*4  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
       COMMON /GRID/ RADX,RADY,RADZ

       REAL*4  RX(0:NAMRX+1,NPALEV),RY(0:NAMRX+1,NPALEV),
     &         RZ(0:NAMRX+1,NPALEV)
       COMMON /GRIDAMR/ RX,RY,RZ

       REAL*4 RETE,HTE,ROTE
       COMMON /BACK/ RETE,HTE,ROTE

       REAL*4 U1(NMAX,NMAY,NMAZ)
       REAL*4 U1G(NMAX,NMAY,NMAZ)
       REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL*4 U11G(NAMRX,NAMRY,NAMRZ,NPALEV)
       COMMON /VARIA/ U1,U11,U1G,U11G

       REAL*4 ACHE,T0,RE0
       COMMON /DOS/ ACHE,T0,RE0

*      LOCAL VARIABLES
       INTEGER,ALLOCATABLE::SOLAPA(:,:),NSOLAP(:)
       INTEGER IRR,I,J,K,IX,JY,KZ,IPATCH,LOWH1,LOWH2,LOWH3,BASINT
       INTEGER IMAXCLUS,IMINCLUS,CONTA,FLAG_ITER,NUM_ITERS,BASINT2
       REAL BAS,BASX,BASY,BASZ,R1,R2,DIST,X1,Y1,Z1,X2,Y2,Z2,VOL1,VOL2
       REAL VINT,XI,PI,SOLAP_LOWER_THR,M1,M2

       SOLAP_LOWER_THR=VOL_SOLAP_LOW

       PI=DACOS(-1.D0)

       WRITE(*,*) '== HALOES OVERLAPPING IN IR =', IR
       WRITE(*,*) 'NCLUS,NHALLEV(IR)=', NCLUS, NHALLEV(IR)

       LOWH1=SUM(NHALLEV(0:IR-1))+1
       LOWH2=SUM(NHALLEV(0:IR))

       ALLOCATE(NSOLAP(LOWH1:LOWH2), SOLAPA(LOWH1:LOWH2,NMAXSOLAP))

!$OMP PARALLEL DO SHARED(NSOLAP,SOLAPA), PRIVATE(I,J), DEFAULT(NONE)
       DO I=LOWH1,LOWH2
        NSOLAP(I)=0
        DO J=1,NMAXSOLAP
         SOLAPA(I,J)=0
        END DO
       END DO

       NUM_ITERS=0
       FLAG_ITER=1

       DO WHILE (FLAG_ITER.EQ.1)
        NUM_ITERS=NUM_ITERS+1
        CONTA=0

*       1. Look for clusters completely included in other ones, and delete
*        the smallest ones.
        DO I=LOWH1,LOWH2
         BASINT=REALCLUS(IFI,I)
         IF (BASINT.NE.0) THEN
          R1=RADIO(I)
          X1=CLUSRX(I)
          Y1=CLUSRY(I)
          Z1=CLUSRZ(I)
          DO J=I+1,LOWH2
           BASINT2=REALCLUS(IFI,J)
           IF (BASINT2.NE.0) THEN
            R2=RADIO(J)
            X2=CLUSRX(J)
            Y2=CLUSRY(J)
            Z2=CLUSRZ(J)
            DIST=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
            IF (ABS(R1-R2).GE.DIST) THEN
             ! A sphere is completely included in the other one.
             ! We shall remove the smallest one.
             IMAXCLUS=I
             IMINCLUS=J
             IF (R2.GT.R1) THEN
              IMAXCLUS=J
              IMINCLUS=I
             END IF !(R2.GT.R1)

             REALCLUS(IFI,IMINCLUS)=0 !The smallest halo is removed
             CONTA=CONTA+1
            END IF !(ABS(R1-R2).GE.DIST)
           END IF !(BASINT2.NE.0)
          END DO !J=I+1,LOWH2
         END IF !(BASINT.NE.0) THEN
        END DO !I=LOWH1,LOWH2

*       2. Look for clusters which overlap a large amount of their
*        volume, and merge them.
        DO I=LOWH1,LOWH2
         BASINT=REALCLUS(IFI,I)
         IF (BASINT.NE.0) THEN
          R1=RADIO(I)
          X1=CLUSRX(I)
          Y1=CLUSRY(I)
          Z1=CLUSRZ(I)
          M1=MASA(I)
          DO J=I+1,LOWH2
           BASINT2=REALCLUS(IFI,J)
           IF (BASINT2.NE.0) THEN
            R2=RADIO(J)
            X2=CLUSRX(J)
            Y2=CLUSRY(J)
            Z2=CLUSRZ(J)
            M2=MASA(J)

            ! Compute the volume overlap
            DIST=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
            BAS=DIST/(R1+R2)
            IF (BAS.LT.1.0) THEN
             XI=BAS
             VINT=(PI/12.0) * (R1+R2)**3 * ((1-XI)**2/XI) *
     &               (XI**2 + 2*XI - 3*(R1-R2)**2/(R1+R2)**2)
             VINT=VINT/(4*PI/3*MIN(R1,R2)**3) ! fraction of the smallest cluster volume overlapped

             IF (VINT.GT.SOLAP_LOWER_THR) THEN
              IMAXCLUS=I
              IMINCLUS=J
              IF (R2.GT.R1) THEN
               IMAXCLUS=J
               IMINCLUS=I
              END IF !(R2.GT.R1)

              REALCLUS(IFI,IMINCLUS)=0 !The smallest halo is removed
              CONTA=CONTA+1

              RADIO(IMAXCLUS)=(R1**3+R2**3)**(1.0/3.0)
              CLUSRX(IMAXCLUS)=(M1*X1+M2*X2)/(M1+M2)
              CLUSRY(IMAXCLUS)=(M1*Y1+M2*Y2)/(M1+M2)
              CLUSRZ(IMAXCLUS)=(M1*Z1+M2*Z2)/(M1+M2)

C              write(*,*) '--------------'
C              write(*,*) i,j,dist,bas
C              write(*,*) x1,y1,z1,r1,m1
C              write(*,*) x2,y2,z2,r2,m2
C              write(*,*) vint
C              write(*,*) '-->',CLUSRX(IMAXCLUS),CLUSRY(IMAXCLUS),
C     &                         CLUSRZ(IMAXCLUS),RADIO(IMAXCLUS)

             END IF !(VINT.GT.SOLAP_LOWER_THR)
            END IF !(BAS.LT.1.0)
           END IF !(BASINT2.NE.0)
          END DO !J=I+1,LOWH2
         END IF !(BASINT.NE.0) THEN
        END DO !I=LOWH1,LOWH2

        IF (CONTA.EQ.0) FLAG_ITER=0
        !WRITE(*,*) 'OVERLAPPING',IR,NUM_ITERS,CONTA
       END DO !WHILE (FLAG_ITER.EQ.1)

       BASINT=COUNT(REALCLUS(IFI,LOWH1:LOWH2).EQ.0)
       WRITE(*,*) 'REMOVED HALOS_0----->', BASINT
       BASINT=COUNT(REALCLUS(IFI,LOWH1:LOWH2).NE.0)
       WRITE(*,*) 'POSSIBLE HALOS_0----->', BASINT
       BASINT=COUNT(REALCLUS(IFI,LOWH1:LOWH2).EQ.-1)
       BASINT2=COUNT(REALCLUS(IFI,LOWH1:LOWH2).GT.0)
       WRITE(*,*) '--> Of which free,substructure:',BASINT,BASINT2
       WRITE(*,*)

       RETURN
       END

********************************************************************
       SUBROUTINE HALOFIND_GRID(IFI,NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                          PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                          PATCHRY,PATCHRZ,PARE,NCLUS,MASA,RADIO,
     &                          CLUSRX,CLUSRY,CLUSRZ,REALCLUS,LEVHAL,
     &                          NHALLEV,BOUND,CONTRASTEC,RODO,SOLAP,
     &                          VECINO,NVECI,CR0AMR,CR0AMR11,PATCHCLUS,
     &                          VOL_SOLAP_LOW)
********************************************************************
*      Pipeline for tentative halo finding over the grid
********************************************************************

       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

*      I/O DATA
       INTEGER IFI,NL,NX,NY,NZ
       INTEGER NPATCH(0:NLEVELS)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
       REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
       INTEGER PARE(NPALEV)
       INTEGER NCLUS
       REAL MASA(MAXNCLUS),RADIO(MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
       INTEGER REALCLUS(MAXITER,MAXNCLUS),LEVHAL(MAXNCLUS)
       INTEGER PATCHCLUS(MAXNCLUS)
       INTEGER NHALLEV(0:NLEVELS)
       REAL BOUND,CONTRASTEC,RODO
       INTEGER SOLAP(NAMRX,NAMRY,NAMRZ,NPALEV)
       INTEGER VECINO(NPALEV,NPALEV),NVECI(NPALEV)
       INTEGER CR0AMR(NMAX,NMAY,NMAZ)
       INTEGER CR0AMR11(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL VOL_SOLAP_LOW

*      GLOBAL VARIABLES
       REAL*4 DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       REAL*4  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
       COMMON /GRID/ RADX,RADY,RADZ

       REAL*4  RX(0:NAMRX+1,NPALEV),RY(0:NAMRX+1,NPALEV),
     &         RZ(0:NAMRX+1,NPALEV)
       COMMON /GRIDAMR/ RX,RY,RZ

       REAL*4 RETE,HTE,ROTE
       COMMON /BACK/ RETE,HTE,ROTE

       REAL*4 U1(NMAX,NMAY,NMAZ)
       REAL*4 U1G(NMAX,NMAY,NMAZ)
       REAL*4 U11(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL*4 U11G(NAMRX,NAMRY,NAMRZ,NPALEV)
       COMMON /VARIA/ U1,U11,U1G,U11G

       REAL*4 ACHE,T0,RE0
       COMMON /DOS/ ACHE,T0,RE0

*      LOCAL VARIABLES
       INTEGER CONTA(NMAX,NMAY,NMAZ)
       INTEGER CONTA1(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL UBAS1(NMAX,NMAY,NMAZ)

       REAL MAXIMO(NPALEV)
       INTEGER VID(NLEVELS,NPALEV),NVID(NLEVELS)
       INTEGER RELEVANT_PATCHES(NPALEV),NRELEVANT_PATCHES(NLEVELS)

       INTEGER IR,IX,JY,KZ,I,J,K,II,JJ,KK,IPATCH,ICEN(3),NV_GOOD
       INTEGER L1,L2,L3,NX1,NX2,NY1,NY2,NZ1,NZ2,KK_ENTERO,ITER_GROW
       INTEGER N1,N2,N3,KONTA,LOW1,LOW2,I2,ICEN1(1),ICEN4(4),CEL,I1,J1
       INTEGER K1,BORAMR,IRR,BASINT,FLAG_ITER
       REAL PRUEBAX,PRUEBAY,PRUEBAZ,RMIN,BASMASS_SHELL,BASMASS,DELTA
       REAL ESP,ESP_LOG,BAS,KK_REAL,RSHELL,R_INT,R_EXT,RANT
       REAL BASDELTA,AA,PI,VOLCELL,BASX,BASY,BASZ,BASVOL
       REAL X1,X2,Y1,Y2,Z1,Z2,DXPA,DYPA,DZPA,BASXX,BASYY,BASZZ
       REAL XCEN,YCEN,ZCEN,BOUNDIR,X3,Y3,Z3,X4,Y4,Z4,MINDERIV
       REAL VECDENS(1000),VECRAD(1000),DERIVATIVE(1000),BASVOL_SHELL

       REAL*4, ALLOCATABLE::DDD(:)
       INTEGER, ALLOCATABLE::DDDX(:),DDDY(:),DDDZ(:),DDDP(:)

**************************************************************
*      NIVEL BASE!!
**************************************************************

       PI=DACOS(-1.D0)

       IR=0
       ESP=0.2*DX
       ESP_LOG=1.05
       WRITE(*,*) '<-------- BASE GRID -------->'
       WRITE(*,*) '--> IR,DX:',IR,DX

!$OMP PARALLEL DO SHARED(NX,NY,NZ,CONTA,UBAS1,U1),
!$OMP+            PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO KZ=1,NZ
       DO JY=1,NY
       DO IX=1,NX
        CONTA(IX,JY,KZ)=0
        UBAS1(IX,JY,KZ)=U1(IX,JY,KZ)
       END DO
       END DO
       END DO

       KK_ENTERO=0
!$OMP PARALLEL DO SHARED(NX,NY,NZ,UBAS1,CONTRASTEC),
!$OMP+            PRIVATE(IX,JY,KZ),
!$OMP+            REDUCTION(+:KK_ENTERO),
!$OMP+            DEFAULT(NONE)
       DO KZ=1,NZ
       DO JY=1,NY
       DO IX=1,NX
        IF (UBAS1(IX,JY,KZ).GE.CONTRASTEC) KK_ENTERO=KK_ENTERO+1
       END DO
       END DO
       END DO
       WRITE(*,*) 'Clean cells above virial contrast',KK_ENTERO
*      Ordenamos todas las celdas de los vecinos
       ALLOCATE(DDD(KK_ENTERO))
       ALLOCATE(DDDX(KK_ENTERO))
       ALLOCATE(DDDY(KK_ENTERO))
       ALLOCATE(DDDZ(KK_ENTERO))

!$OMP PARALLEL DO SHARED(DDD,DDDX,DDDY,DDDZ,KK_ENTERO),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
       DO I=1,KK_ENTERO
        DDD(I)=0.0
        DDDX(I)=0
        DDDY(I)=0
        DDDZ(I)=0
       END DO

       II=0
       DO KZ=1,NZ
       DO JY=1,NY
       DO IX=1,NX
        BAS=UBAS1(IX,JY,KZ)
        IF (BAS.GE.CONTRASTEC) THEN
         II=II+1
         DDD(II)=BAS
         DDDX(II)=IX
         DDDY(II)=JY
         DDDZ(II)=KZ
        END IF
       END DO
       END DO
       END DO

       IF (II.NE.KK_ENTERO) THEN
         WRITE(*,*) 'WARNING: bad allocation of DDD', II, KK_ENTERO
         STOP
       END IF

       CALL SORT_CELLS(KK_ENTERO,DDD,DDDX,DDDY,DDDZ)

       NV_GOOD=KK_ENTERO
       KK_ENTERO=0

*      Go through all the center candidates
       DO L1=1,NV_GOOD
        ICEN(1)=DDDX(L1)
        ICEN(2)=DDDY(L1)
        ICEN(3)=DDDZ(L1)

        KK_ENTERO=CONTA(ICEN(1),ICEN(2),ICEN(3))
        IF(KK_ENTERO.EQ.0) THEN ! this means this peak is not inside a halo yet
c         WRITE(*,*) U1(ICEN(1),ICEN(2),ICEN(3))
         NCLUS=NCLUS+1
         REALCLUS(IFI,NCLUS)=-1
         LEVHAL(NCLUS)=IR
         NHALLEV(IR)=NHALLEV(IR)+1
         PATCHCLUS(NCLUS)=0

         IF(NCLUS.GT.MAXNCLUS) THEN
          WRITE(*,*) 'WARNING: NCLUS>MAXNCLUS!!!',NCLUS,MAXNCLUS
          STOP
         END IF

         CLUSRX(NCLUS)=RADX(ICEN(1))
         CLUSRY(NCLUS)=RADY(ICEN(2))
         CLUSRZ(NCLUS)=RADZ(ICEN(3))

*        tentative reach of the halo ---> build a mini box around it
*        (by excess, set as a parameter in asohf.dat)
         BASX=CLUSRX(NCLUS)-BOUND
         NX1=INT(((BASX-RADX(1))/DX)+0.5)+1
         IF (NX1.LT.1) NX1=1

         BASX=CLUSRX(NCLUS)+BOUND
         NX2=INT(((BASX-RADX(1))/DX)+0.5)+1
         IF (NX2.GT.NX) NX2=NX

         BASY=CLUSRY(NCLUS)-BOUND
         NY1=INT(((BASY-RADY(1))/DY)+0.5)+1
         IF (NY1.LT.1) NY1=1

         BASY=CLUSRY(NCLUS)+BOUND
         NY2=INT(((BASY-RADY(1))/DY)+0.5)+1
         IF (NY2.GT.NY) NY2=NY

         BASZ=CLUSRZ(NCLUS)-BOUND
         NZ1=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
         IF (NZ1.LT.1) NZ1=1

         BASZ=CLUSRZ(NCLUS)+BOUND
         NZ2=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
         IF (NZ2.GT.NZ) NZ2=NZ

*        Now, we extend the cluster radially from its center
         BASX=0.0
         BASY=0.0
         BASZ=0.0
         BASDELTA=0.0

         BASMASS_SHELL=0.0
         BASMASS=0.0   !TOTAL MASS OF THE CLUSTER
         DELTA=0.0    !TOTAL CONTRAST OF THE CLUSTER
         BASVOL=0.0     !TOTAL VOLUME OF THE CLUSTER (sphere)

         R_INT=0.0
         R_EXT=0.5*DX

*        increase the radius until density falls below the virial value
         DELTA=10.0*CONTRASTEC*ROTE ! this is to ensure we enter the loop
         ITER_GROW=0
         DO WHILE(DELTA.GT.CONTRASTEC*ROTE)
          ITER_GROW=ITER_GROW+1

          IF (ITER_GROW.GT.1) THEN
           IF (R_EXT.LE.BOUND) THEN
            R_INT=R_EXT
            R_EXT=MAX(R_EXT+ESP, R_EXT*ESP_LOG)
           ELSE
            WRITE(*,*) 'WARNING: growing not converged', r_int, r_ext,
     &                                                   bound
            STOP
           END IF
          END IF

          VOLCELL=DX*DY*DZ
          BASMASS_SHELL=0.0
          II=0
          DO K=NZ1,NZ2
          DO J=NY1,NY2
          DO I=NX1,NX2
           AA=SQRT((RADX(I)-CLUSRX(NCLUS))**2 +
     &             (RADY(J)-CLUSRY(NCLUS))**2 +
     &             (RADZ(K)-CLUSRZ(NCLUS))**2)

           IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
            CONTA(I,J,K)=1 ! do not try to find an additional halo here (at this level)

            BASVOL=BASVOL+VOLCELL

            BAS=U1(I,J,K)*VOLCELL !U1 is not density contrast, but 1+delta = rho/rho_B!!!
            BASMASS_SHELL=BASMASS_SHELL+BAS

            BASX=BASX+RADX(I)*BAS
            BASY=BASY+RADY(J)*BAS
            BASZ=BASZ+RADZ(K)*BAS
            BASDELTA=BASDELTA+BAS

            II=II+1
           END IF

          END DO
          END DO
          END DO

          BASMASS=BASMASS+BASMASS_SHELL*RODO*RE0**3
          DELTA=BASMASS/(BASVOL*RETE**3)
c          WRITE(*,*) DELTA/ROTE, II
         END DO   ! do while (DELTA)

         RADIO(NCLUS)=R_EXT
         MASA(NCLUS)=BASMASS

         CLUSRX(NCLUS)=BASX/BASDELTA
         CLUSRY(NCLUS)=BASY/BASDELTA
         CLUSRZ(NCLUS)=BASZ/BASDELTA

C         WRITE(*,*) CLUSRX(NCLUS),CLUSRY(NCLUS),CLUSRZ(NCLUS),
C     &             RADIO(NCLUS),MASA(NCLUS)*9.1717E18,IR

        END IF ! KK_ENTERO.EQ.0
       END DO ! I=1,NV_GOOD

       DEALLOCATE(DDD,DDDX,DDDY,DDDZ)

       WRITE(*,*) 'At level',IR,', num. haloes:', NHALLEV(IR)

****************************************************
*      CORRECION DE SOLAPES:
*      VAMOS A VER QUE CUMULOS SOLAPAN EN IR
****************************************************

c       CALL OVERLAPPING(IFI,IR,NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
c     &                        PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
c     &                        PATCHRY,PATCHRZ,PARE,NCLUS,MASA,RADIO,
c     &                        CLUSRX,CLUSRY,CLUSRZ,REALCLUS,LEVHAL,
c     &                        NHALLEV,BOUND,CONTRASTEC,RODO,SOLAP,
c     &                        VECINO,NVECI,CR0AMR,CR0AMR11)

****************************************************
*      FIN CORRECION DE SOLAPES EN IR
****************************************************
       WRITE(*,*) 'End of base level', 0

       IF (NL.GT.0) THEN
        WRITE(*,*)
        WRITE(*,*) '<-------- Now proceeding with the',NL,
     &             'AMR levels -------->'
*      Find the l=1 cells covered by l=0 haloes (they will
*      potentially host substructures)
!$OMP PARALLEL DO SHARED(NPATCH,PATCHNX,PATCHNY,PATCHNZ,CONTA1),
!$OMP+            PRIVATE(I,N1,N2,N3,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
         DO I=1,NPATCH(1)
          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)
          DO IX=1,N1
          DO JY=1,N2
          DO KZ=1,N3
           CONTA1(IX,JY,KZ,I)=1
          END DO
          END DO
          END DO
         END DO

!$OMP PARALLEL DO SHARED(NPATCH,PATCHNX,PATCHNY,PATCHNZ,PATCHRX,PATCHRY,
!$OMP+                   PATCHRZ,DX,DY,DZ,CLUSRX,CLUSRY,CLUSRZ,RADIO,
!$OMP+                   RX,RY,RZ,CONTA1,NCLUS),
!$OMP+            PRIVATE(I,N1,N2,N3,X1,X2,Y1,Y2,Z1,Z2,BASX,BASY,BASZ,
!$OMP+                    BAS,IX,JY,KZ,AA,X3,Y3,Z3,X4,Y4,Z4),
!$OMP+            DEFAULT(NONE)
         DO I=1,NPATCH(1)
          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)
          X1=PATCHRX(I)-DX/2.0
          Y1=PATCHRY(I)-DY/2.0
          Z1=PATCHRZ(I)-DZ/2.0
          X2=X1+N1*DX/2.0
          Y2=Y1+N2*DY/2.0
          Z2=Z1+N3*DZ/2.0
          DO II=1,NCLUS
           BASX=CLUSRX(II)
           BASY=CLUSRY(II)
           BASZ=CLUSRZ(II)
           BAS=RADIO(II)
           X3=BASX-BAS
           Y3=BASY-BAS
           Z3=BASZ-BAS
           X4=BASX+BAS
           Y4=BASY+BAS
           Z4=BASZ+BAS
           IF (X1.LE.X4.AND.X3.LE.X2.AND.
     &         Y1.LE.Y4.AND.Y3.LE.Y2.AND.
     &         Z1.LE.Z4.AND.Z3.LE.Z2) THEN
c           WRITE(*,*) BASX,BASY,BASZ,BAS
            DO IX=1,N1
            DO JY=1,N2
            DO KZ=1,N3
             AA=SQRT((RX(IX,I)-BASX)**2 +
     &               (RY(JY,I)-BASY)**2 +
     &               (RZ(KZ,I)-BASZ)**2)
             IF (AA.LE.BAS) CONTA1(IX,JY,KZ,I)=-1
            END DO
            END DO
            END DO
           END IF
          END DO
         END DO

*       And mark the centers of haloes (to avoid identifying haloes at same levels as substructure)
!$OMP PARALLEL DO SHARED(NPATCH,DX,DY,DZ,PATCHNX,PATCHNY,PATCHNZ,
!$OMP+                   PATCHRX,PATCHRY,PATCHRZ,NCLUS,CLUSRX,CLUSRY,
!$OMP+                   CLUSRZ,CONTA1),
!$OMP+            PRIVATE(I,N1,N2,N3,X1,X2,Y1,Y2,Z1,Z2,II,BASX,BASY,
!$OMP+                    BASZ,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
        DO I=1,NPATCH(1)
         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)
         X1=PATCHRX(I)-DX/2.0
         Y1=PATCHRY(I)-DY/2.0
         Z1=PATCHRZ(I)-DZ/2.0
         X2=X1+N1*DX/2.0
         Y2=Y1+N2*DY/2.0
         Z2=Z1+N3*DZ/2.0
         DO II=1,NCLUS
          BASX=CLUSRX(II)
          BASY=CLUSRY(II)
          BASZ=CLUSRZ(II)
          BASX=(BASX-X1)*(X2-BASX)
          BASY=(BASY-Y1)*(Y2-BASY)
          BASZ=(BASZ-Z1)*(Z2-BASZ)
          IF (BASX.GT.0.AND.BASY.GT.0.AND.BASZ.GT.0) THEN
           BASX=CLUSRX(II)-X1
           BASY=CLUSRY(II)-Y1
           BASZ=CLUSRZ(II)-Z1
           IX=INT(BASX/(DX/2.0))+1
           JY=INT(BASY/(DY/2.0))+1
           KZ=INT(BASZ/(DZ/2.0))+1
           DO I1=IX-1,IX+1
           DO J1=JY-1,JY+1
           DO K1=KZ-1,KZ+1
            IF (I1.GE.1.AND.I1.LE.N1.AND.
     &          J1.GE.1.AND.J1.LE.N2.AND.
     &          K1.GE.1.AND.K1.LE.N3) THEN
             CONTA1(I1,J1,K1,I)=-2
            END IF
           END DO
           END DO
           END DO
          END IF
         END DO
        END DO

        ! clean conta1 from overlaps
        CALL CLEAN_OVERLAPS_INT(NL,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                          SOLAP,CONTA1)

*       CONTA1=0 --> the cell is overlapping another IR cell, and slave (not to be counted)
*       CONTA1=1 --> outside any halo
*       CONTA1=2 --> inside a halo at this same level IR
*       CONTA1=-1 --> outside haloes at this same level, but inside an IR-1 halo
*       CONTA1=-2 --> outside haloes at this same level, but center of a IR-1 halo

        DO IR=1,NL
        ! then proceed with the halo finding as written in the asohf.f now, but with the new modifications
        ! especially, cycling not only around the vecinos, but also around the pares....
         DXPA=DX/(2.0**IR)
         DYPA=DY/(2.0**IR)
         DZPA=DZ/(2.0**IR)

         WRITE(*,*) '----> NEW LEVEL: IR,DX=',IR,DXPA,'<----'

         LOW1=SUM(NPATCH(0:IR-1))+1
         LOW2=SUM(NPATCH(0:IR))

         ESP=0.2*DXPA
         ESP_LOG=1.05
         BORAMR=1
         BOUNDIR=BOUND/1.75**IR

*        estimation to allocate
         KK_ENTERO=0
!$OMP  PARALLEL DO SHARED(PATCHNX,PATCHNY,PATCHNZ,U11,CONTRASTEC,LOW1,
!$OMP+                    LOW2,CONTA1,BORAMR),
!$OMP+             PRIVATE(N1,N2,N3,I,IX,JY,KZ),
!$OMP+             REDUCTION(+:KK_ENTERO),
!$OMP+             DEFAULT(NONE)
         DO I=LOW1,LOW2
          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)
          DO KZ=1+BORAMR,N3-BORAMR
          DO JY=1+BORAMR,N2-BORAMR
          DO IX=1+BORAMR,N1-BORAMR
           IF (U11(IX,JY,KZ,I).GE.CONTRASTEC.AND.
     &         CONTA1(IX,JY,KZ,I).NE.0) KK_ENTERO=KK_ENTERO+1
          END DO
          END DO
          END DO
         END DO
         WRITE(*,*) 'Clean cells above virial contrast',KK_ENTERO

         ALLOCATE(DDD(KK_ENTERO))
         ALLOCATE(DDDX(KK_ENTERO))
         ALLOCATE(DDDY(KK_ENTERO))
         ALLOCATE(DDDZ(KK_ENTERO))
         ALLOCATE(DDDP(KK_ENTERO))

!$OMP PARALLEL DO SHARED(DDD,DDDX,DDDY,DDDZ,DDDP,KK_ENTERO),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
         DO I=1,KK_ENTERO
          DDD(I)=0.0
          DDDX(I)=0
          DDDY(I)=0
          DDDZ(I)=0
          DDDP(I)=0
         END DO

         II=0
         DO I=LOW1,LOW2
          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)
          DO KZ=1+BORAMR,N3-BORAMR
          DO JY=1+BORAMR,N2-BORAMR
          DO IX=1+BORAMR,N1-BORAMR
           IF (U11(IX,JY,KZ,I).GE.CONTRASTEC.AND.
     &         CONTA1(IX,JY,KZ,I).NE.0) THEN
            II=II+1
            DDD(II)=U11(IX,JY,KZ,I)
            DDDX(II)=IX
            DDDY(II)=JY
            DDDZ(II)=KZ
            DDDP(II)=I
           END IF
          END DO
          END DO
          END DO
         END DO

         !WRITE(*,*) 'CHECK KK_ENTERO,II=',KK_ENTERO,II

         CALL SORT_CELLS_AMR(KK_ENTERO,DDD,DDDX,DDDY,DDDZ,DDDP)

c         DO IX=1,KK_ENTERO
c          WRITE(*,*) DDD(IX),DDDP(IX),DDDX(IX),DDDY(IX),DDDZ(IX),
c     &               CONTA1(DDDX(IX),DDDY(IX),DDDZ(IX),DDDP(IX))
c         END DO

         NV_GOOD=KK_ENTERO
         KK_ENTERO=0

         ! cycle through the cells, and see how to grow the halo, depending
         ! on whether it is virialised halo or substructure
         ! problem: we might need to check whether candidates to substructure
         ! are bonafide density maxima, since excluding the center cell
         ! does not ensure the code will not pick a dense cell just next to it.
*        Go through all the center candidates
         DO L1=1,NV_GOOD
          ICEN4(1)=DDDX(L1)
          ICEN4(2)=DDDY(L1)
          ICEN4(3)=DDDZ(L1)
          ICEN4(4)=DDDP(L1)
          IX=ICEN4(1)
          JY=ICEN4(2)
          KZ=ICEN4(3)
          IPATCH=ICEN4(4)
          XCEN=RX(IX,IPATCH)
          YCEN=RY(JY,IPATCH)
          ZCEN=RZ(KZ,IPATCH)

          KK_ENTERO=CONTA1(IX,JY,KZ,IPATCH)
          IF (KK_ENTERO.EQ.1) THEN ! this means this peak is not inside a halo yet
           NCLUS=NCLUS+1
           REALCLUS(IFI,NCLUS)=-1
           LEVHAL(NCLUS)=IR
           NHALLEV(IR)=NHALLEV(IR)+1
           PATCHCLUS(NCLUS)=IPATCH

           IF(NCLUS.GT.MAXNCLUS) THEN
            WRITE(*,*) 'WARNING: NCLUS>MAXNCLUS!!!',NCLUS,MAXNCLUS
            STOP
           END IF

           CLUSRX(NCLUS)=XCEN
           CLUSRY(NCLUS)=YCEN
           CLUSRZ(NCLUS)=ZCEN

*          tentative reach of the base grid
           BASX=CLUSRX(NCLUS)-BOUNDIR
           NX1=INT(((BASX-RADX(1))/DX)+0.5)+1
           IF (NX1.LT.1) NX1=1

           BASX=CLUSRX(NCLUS)+BOUNDIR
           NX2=INT(((BASX-RADX(1))/DX)+0.5)+1
           IF (NX2.GT.NX) NX2=NX

           BASY=CLUSRY(NCLUS)-BOUNDIR
           NY1=INT(((BASY-RADY(1))/DY)+0.5)+1
           IF (NY1.LT.1) NY1=1

           BASY=CLUSRY(NCLUS)+BOUNDIR
           NY2=INT(((BASY-RADY(1))/DY)+0.5)+1
           IF (NY2.GT.NY) NY2=NY

           BASZ=CLUSRZ(NCLUS)-BOUNDIR
           NZ1=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
           IF (NZ1.LT.1) NZ1=1

           BASZ=CLUSRZ(NCLUS)+BOUNDIR
           NZ2=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
           IF (NZ2.GT.NZ) NZ2=NZ

*          patches that could contain this halo
           CALL PATCHES_SPHERE(NPATCH,PATCHRX,PATCHRY,PATCHRZ,PATCHNX,
     &                         PATCHNY,PATCHNZ,XCEN,YCEN,ZCEN,BOUNDIR,
     &                         IR,NL,RELEVANT_PATCHES,NRELEVANT_PATCHES)

c           WRITE(*,*) CLUSRX(NCLUS),NX1,NX2,RADX(NX1),RADX(NX2)
c           DO I=1,NRELEVANT_PATCHES(1)
c            WRITE(*,*) I,PATCHRX(RELEVANT_PATCHES(I))-DX/2.0,
c     &                 PATCHRX(RELEVANT_PATCHES(I))+
c     &                 (PATCHNX(RELEVANT_PATCHES(I))-1)*DX/2.0
c           END DO

*          Now, we extend the cluster radially from its center
           BASX=0.0
           BASY=0.0
           BASZ=0.0
           BASDELTA=0.0

           BASMASS_SHELL=0.0
           BASMASS=0.0   !TOTAL MASS OF THE CLUSTER
           DELTA=0.0    !TOTAL CONTRAST OF THE CLUSTER
           BASVOL=0.0     !TOTAL VOLUME OF THE CLUSTER (sphere)

           R_INT=0.0
           R_EXT=0.5*DXPA

*          increase the radius until density falls below the virial value
           DELTA=10.0*CONTRASTEC*ROTE ! this is to ensure we enter the loop
           ITER_GROW=0
           II=0
           JJ=0

           DO WHILE(DELTA.GT.CONTRASTEC*ROTE)
            ITER_GROW=ITER_GROW+1

            IF (ITER_GROW.GT.1) THEN
             IF (R_EXT.LE.BOUNDIR) THEN
              R_INT=R_EXT
              R_EXT=MAX(R_EXT+ESP, R_EXT*ESP_LOG)
             ELSE
              WRITE(*,*) 'WARNING: growing not converged', r_int, r_ext,
     &                    boundir,iter_grow,delta/rote,kk_entero
              STOP
             END IF
            END IF

            BASMASS_SHELL=0.0

            VOLCELL=DX*DY*DZ
            DO K=NZ1,NZ2
            DO J=NY1,NY2
            DO I=NX1,NX2
             IF (CR0AMR(I,J,K).EQ.1) THEN
              AA=SQRT((RADX(I)-CLUSRX(NCLUS))**2 +
     &                (RADY(J)-CLUSRY(NCLUS))**2 +
     &                (RADZ(K)-CLUSRZ(NCLUS))**2)

              IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
               !CONTA(I,J,K)=1 ! do not try to find an additional halo here (at this level)
               II=II+1
               BASVOL=BASVOL+VOLCELL

               BAS=U1(I,J,K)*VOLCELL !U1 is not density contrast, but 1+delta = rho/rho_B!!!
               BASMASS_SHELL=BASMASS_SHELL+BAS

               BASX=BASX+RADX(I)*BAS
               BASY=BASY+RADY(J)*BAS
               BASZ=BASZ+RADZ(K)*BAS
               BASDELTA=BASDELTA+BAS
              END IF
             END IF
            END DO
            END DO
            END DO

            ! AND NOW LOOP UP TO IR-1, AND THEN ANOTHER LOOP ON IR, OVER RELEVANT_PATCHES
            DO IRR=1,IR-1
             DXPA=DX/2.0**IRR
             DYPA=DY/2.0**IRR
             DZPA=DZ/2.0**IRR
             LOW1=SUM(NRELEVANT_PATCHES(1:IRR-1))+1
             LOW2=SUM(NRELEVANT_PATCHES(1:IRR))
             VOLCELL=DXPA*DYPA*DZPA
             DO BASINT=LOW1,LOW2
              IPATCH=RELEVANT_PATCHES(BASINT)
              N1=PATCHNX(IPATCH)
              N2=PATCHNY(IPATCH)
              N3=PATCHNZ(IPATCH)
              DO KZ=1,N3
              DO JY=1,N2
              DO IX=1,N1
               IF (CR0AMR11(IX,JY,KZ,IPATCH).EQ.1) THEN
                IF (CONTA1(IX,JY,KZ,IPATCH).NE.0) THEN
                 ! DO STUFF
                 AA=SQRT((RX(IX,IPATCH)-CLUSRX(NCLUS))**2 +
     &                   (RY(JY,IPATCH)-CLUSRY(NCLUS))**2 +
     &                   (RZ(KZ,IPATCH)-CLUSRZ(NCLUS))**2)

                 IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
                  II=II+1
                  BASVOL=BASVOL+VOLCELL

                  BAS=U11(IX,JY,KZ,IPATCH)*VOLCELL !U1 is not density contrast, but 1+delta = rho/rho_B!!!
                  BASMASS_SHELL=BASMASS_SHELL+BAS

                  BASX=BASX+RX(IX,IPATCH)*BAS
                  BASY=BASY+RY(JY,IPATCH)*BAS
                  BASZ=BASZ+RZ(KZ,IPATCH)*BAS
                  BASDELTA=BASDELTA+BAS
                 END IF
                END IF
               END IF
              END DO
              END DO
              END DO
             END DO
            END DO

            IRR=IR
            DXPA=DX/(2.0**IRR)
            DYPA=DY/(2.0**IRR)
            DZPA=DZ/(2.0**IRR)
            LOW1=SUM(NRELEVANT_PATCHES(1:IRR-1))+1
            LOW2=SUM(NRELEVANT_PATCHES(1:IRR))
            VOLCELL=DXPA*DYPA*DZPA
            DO BASINT=LOW1,LOW2
             IPATCH=RELEVANT_PATCHES(BASINT)
             !if (ir.gt.1) write(*,*) 'ipatch=',ipatch,low1,low2,basint
             N1=PATCHNX(IPATCH)
             N2=PATCHNY(IPATCH)
             N3=PATCHNZ(IPATCH)
             DO KZ=1,N3
             DO JY=1,N2
             DO IX=1,N1
              IF (CONTA1(IX,JY,KZ,IPATCH).NE.0) THEN

               AA=SQRT((RX(IX,IPATCH)-CLUSRX(NCLUS))**2 +
     &                 (RY(JY,IPATCH)-CLUSRY(NCLUS))**2 +
     &                 (RZ(KZ,IPATCH)-CLUSRZ(NCLUS))**2)

               IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
                JJ=JJ+1
                CONTA1(IX,JY,KZ,IPATCH)=2 ! do not try to find an additional halo here (at this level)
                BASVOL=BASVOL+VOLCELL

                BAS=U11(IX,JY,KZ,IPATCH)*VOLCELL !U1 is not density contrast, but 1+delta = rho/rho_B!!!
                BASMASS_SHELL=BASMASS_SHELL+BAS

                BASX=BASX+RX(IX,IPATCH)*BAS
                BASY=BASY+RY(JY,IPATCH)*BAS
                BASZ=BASZ+RZ(KZ,IPATCH)*BAS
                BASDELTA=BASDELTA+BAS
               END IF
              END IF
             END DO
             END DO
             END DO
            END DO

            BASMASS=BASMASS+BASMASS_SHELL*RODO*RE0**3
            DELTA=BASMASS/(BASVOL*RETE**3)

c            write(*,*) iter_grow,r_ext,
c     &                 clusrx(nclus),clusry(nclus),clusrz(nclus),
c     &                 delta/rote
           END DO   ! do while (DELTA)

           RADIO(NCLUS)=R_EXT
           MASA(NCLUS)=BASMASS

           CLUSRX(NCLUS)=BASX/BASDELTA
           CLUSRY(NCLUS)=BASY/BASDELTA
           CLUSRZ(NCLUS)=BASZ/BASDELTA

c           WRITE(*,*) IR,CLUSRX(NCLUS),CLUSRY(NCLUS),CLUSRZ(NCLUS),
c     &             RADIO(NCLUS),MASA(NCLUS)*9.1717E18,II,JJ,DELTA/ROTE,
c     &             NCLUS
          ELSE IF(KK_ENTERO.EQ.-1) THEN ! this mean this peak must be a substructure
           BASX =U11(IX+1,JY,KZ,IPATCH)-U11(IX,JY,KZ,IPATCH)
           BASY =U11(IX,JY+1,KZ,IPATCH)-U11(IX,JY,KZ,IPATCH)
           BASZ =U11(IX,JY,KZ+1,IPATCH)-U11(IX,JY,KZ,IPATCH)
           BASXX=U11(IX,JY,KZ,IPATCH)  -U11(IX-1,JY,KZ,IPATCH)
           BASYY=U11(IX,JY,KZ,IPATCH)  -U11(IX,JY-1,KZ,IPATCH)
           BASZZ=U11(IX,JY,KZ,IPATCH)  -U11(IX,JY,KZ-1,IPATCH)
           IF (BASX.LT.0) THEN
           IF (BASY.LT.0) THEN
           IF (BASZ.LT.0) THEN
           IF (BASXX.GT.0) THEN
           IF (BASYY.GT.0) THEN
           IF (BASZZ.GT.0) THEN ! then it's a local maximum
C            WRITE(*,*) XCEN,YCEN,ZCEN,'substructure',ir
C            WRITE(*,*) 'CHECK:',
C     &                 MINVAL(U11(IX-1:IX+1,JY-1:JY+1,KZ-1:KZ+1,IPATCH))

            NCLUS=NCLUS+1
            REALCLUS(IFI,NCLUS)=-1
            LEVHAL(NCLUS)=IR
            NHALLEV(IR)=NHALLEV(IR)+1
            PATCHCLUS(NCLUS)=IPATCH

            IF(NCLUS.GT.MAXNCLUS) THEN
             WRITE(*,*) 'WARNING: NCLUS>MAXNCLUS!!!',NCLUS,MAXNCLUS
             STOP
            END IF

            CLUSRX(NCLUS)=XCEN
            CLUSRY(NCLUS)=YCEN
            CLUSRZ(NCLUS)=ZCEN

*           tentative reach of the base grid
            BASX=CLUSRX(NCLUS)-BOUNDIR
            NX1=INT(((BASX-RADX(1))/DX)+0.5)+1
            IF (NX1.LT.1) NX1=1

            BASX=CLUSRX(NCLUS)+BOUNDIR
            NX2=INT(((BASX-RADX(1))/DX)+0.5)+1
            IF (NX2.GT.NX) NX2=NX

            BASY=CLUSRY(NCLUS)-BOUNDIR
            NY1=INT(((BASY-RADY(1))/DY)+0.5)+1
            IF (NY1.LT.1) NY1=1

            BASY=CLUSRY(NCLUS)+BOUNDIR
            NY2=INT(((BASY-RADY(1))/DY)+0.5)+1
            IF (NY2.GT.NY) NY2=NY

            BASZ=CLUSRZ(NCLUS)-BOUNDIR
            NZ1=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
            IF (NZ1.LT.1) NZ1=1

            BASZ=CLUSRZ(NCLUS)+BOUNDIR
            NZ2=INT(((BASZ-RADZ(1))/DZ)+0.5)+1
            IF (NZ2.GT.NZ) NZ2=NZ

*           patches that could contain this halo
            CALL PATCHES_SPHERE(NPATCH,PATCHRX,PATCHRY,PATCHRZ,PATCHNX,
     &                          PATCHNY,PATCHNZ,XCEN,YCEN,ZCEN,BOUNDIR,
     &                          IR,NL,RELEVANT_PATCHES,
     &                          NRELEVANT_PATCHES)

*           Now, we extend the cluster radially from its center
            BASX=0.0
            BASY=0.0
            BASZ=0.0
            BASDELTA=0.0

            BASMASS_SHELL=0.0
            BASMASS=0.0   !TOTAL MASS OF THE CLUSTER
            DELTA=0.0    !TOTAL CONTRAST OF THE CLUSTER
            BASVOL=0.0     !TOTAL VOLUME OF THE CLUSTER (sphere)

            R_INT=0.0
            R_EXT=0.5*DXPA

*           increase the radius until density falls below the virial value
            DELTA=10.0*CONTRASTEC*ROTE ! this is to ensure we enter the loop
            ITER_GROW=0
            FLAG_ITER=1
            JJ=0
            MINDERIV=1000.0

            DO WHILE(DELTA.GT.CONTRASTEC*ROTE.AND.FLAG_ITER.EQ.1)
             ITER_GROW=ITER_GROW+1
             II=0

             IF (ITER_GROW.GT.1) THEN
              IF (R_EXT.LE.BOUNDIR) THEN
               R_INT=R_EXT
               R_EXT=MAX(R_EXT+ESP, R_EXT*ESP_LOG)
              ELSE
               WRITE(*,*) 'WARNING: growing not converged', r_int,
     &                    r_ext,boundir,iter_grow,delta/rote,kk_entero
               STOP
              END IF
             END IF

             BASMASS_SHELL=0.0
             BASVOL_SHELL=0.0

             VOLCELL=DX*DY*DZ
             DO K=NZ1,NZ2
             DO J=NY1,NY2
             DO I=NX1,NX2
              IF (CR0AMR(I,J,K).EQ.1) THEN
               AA=SQRT((RADX(I)-CLUSRX(NCLUS))**2 +
     &                 (RADY(J)-CLUSRY(NCLUS))**2 +
     &                 (RADZ(K)-CLUSRZ(NCLUS))**2)

               IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
                !CONTA(I,J,K)=1 ! do not try to find an additional halo here (at this level)
                II=II+1
                BASVOL=BASVOL+VOLCELL

                BAS=U1(I,J,K)*VOLCELL !U1 is not density contrast, but 1+delta = rho/rho_B!!!
                BASMASS_SHELL=BASMASS_SHELL+BAS
                BASVOL_SHELL=BASVOL_SHELL+VOLCELL

                BASX=BASX+RADX(I)*BAS
                BASY=BASY+RADY(J)*BAS
                BASZ=BASZ+RADZ(K)*BAS
                BASDELTA=BASDELTA+BAS
               END IF
              END IF
             END DO
             END DO
             END DO

             ! AND NOW LOOP UP TO IR-1, AND THEN ANOTHER LOOP ON IR, OVER RELEVANT_PATCHES
             DO IRR=1,IR-1
              DXPA=DX/2.0**IRR
              DYPA=DY/2.0**IRR
              DZPA=DZ/2.0**IRR
              LOW1=SUM(NRELEVANT_PATCHES(1:IRR-1))+1
              LOW2=SUM(NRELEVANT_PATCHES(1:IRR))
              VOLCELL=DXPA*DYPA*DZPA
              DO BASINT=LOW1,LOW2
               IPATCH=RELEVANT_PATCHES(BASINT)
               N1=PATCHNX(IPATCH)
               N2=PATCHNY(IPATCH)
               N3=PATCHNZ(IPATCH)
               DO KZ=1,N3
               DO JY=1,N2
               DO IX=1,N1
                IF (CR0AMR11(IX,JY,KZ,IPATCH).EQ.1) THEN
                 IF (CONTA1(IX,JY,KZ,IPATCH).NE.0) THEN
                  ! DO STUFF
                  AA=SQRT((RX(IX,IPATCH)-CLUSRX(NCLUS))**2 +
     &                    (RY(JY,IPATCH)-CLUSRY(NCLUS))**2 +
     &                    (RZ(KZ,IPATCH)-CLUSRZ(NCLUS))**2)

                  IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
                   II=II+1
                   BASVOL=BASVOL+VOLCELL

                   BAS=U11(IX,JY,KZ,IPATCH)*VOLCELL !U1 is not density contrast, but 1+delta = rho/rho_B!!!
                   BASMASS_SHELL=BASMASS_SHELL+BAS
                   BASVOL_SHELL=BASVOL_SHELL+VOLCELL

                   BASX=BASX+RX(IX,IPATCH)*BAS
                   BASY=BASY+RY(JY,IPATCH)*BAS
                   BASZ=BASZ+RZ(KZ,IPATCH)*BAS
                   BASDELTA=BASDELTA+BAS
                  END IF
                 END IF
                END IF
               END DO
               END DO
               END DO
              END DO
             END DO

             IRR=IR
             DXPA=DX/(2.0**IRR)
             DYPA=DY/(2.0**IRR)
             DZPA=DZ/(2.0**IRR)
             LOW1=SUM(NRELEVANT_PATCHES(1:IRR-1))+1
             LOW2=SUM(NRELEVANT_PATCHES(1:IRR))
             VOLCELL=DXPA*DYPA*DZPA
             DO BASINT=LOW1,LOW2
              IPATCH=RELEVANT_PATCHES(BASINT)
              !if (ir.gt.1) write(*,*) 'ipatch=',ipatch,low1,low2,basint
              N1=PATCHNX(IPATCH)
              N2=PATCHNY(IPATCH)
              N3=PATCHNZ(IPATCH)
              DO KZ=1,N3
              DO JY=1,N2
              DO IX=1,N1
               IF (CONTA1(IX,JY,KZ,IPATCH).NE.0) THEN

                AA=SQRT((RX(IX,IPATCH)-CLUSRX(NCLUS))**2 +
     &                  (RY(JY,IPATCH)-CLUSRY(NCLUS))**2 +
     &                  (RZ(KZ,IPATCH)-CLUSRZ(NCLUS))**2)

                IF (AA.GE.R_INT.AND.AA.LT.R_EXT) THEN
                 II=II+1
                 CONTA1(IX,JY,KZ,IPATCH)=2 ! do not try to find an additional halo here (at this level)
                 BASVOL=BASVOL+VOLCELL

                 BAS=U11(IX,JY,KZ,IPATCH)*VOLCELL !U1 is not density contrast, but 1+delta = rho/rho_B!!!
                 BASMASS_SHELL=BASMASS_SHELL+BAS
                 BASVOL_SHELL=BASVOL_SHELL+VOLCELL

                 BASX=BASX+RX(IX,IPATCH)*BAS
                 BASY=BASY+RY(JY,IPATCH)*BAS
                 BASZ=BASZ+RZ(KZ,IPATCH)*BAS
                 BASDELTA=BASDELTA+BAS
                END IF
               END IF
              END DO
              END DO
              END DO
             END DO

             BASMASS=BASMASS+BASMASS_SHELL*RODO*RE0**3
             DELTA=BASMASS/(BASVOL*RETE**3)

             IF (II.GT.0) THEN
              JJ=JJ+1
              VECRAD(JJ)=0.5*(R_INT+R_EXT)
              VECDENS(JJ)=BASMASS_SHELL*RODO*RE0**3 /
     &                     (BASVOL_SHELL*RETE**3)
              IF (JJ.EQ.1) THEN
               DERIVATIVE(JJ)=1000.0
              ELSE
               BAS=LOG(VECDENS(JJ))-LOG(VECDENS(JJ-1))
               DERIVATIVE(JJ)=BAS/(LOG(VECRAD(JJ))-LOG(VECRAD(JJ-1)))
               IF (DERIVATIVE(JJ).GT.MINDERIV) FLAG_ITER=0
              END IF
              BAS=DERIVATIVE(JJ)
              MINDERIV=MIN(MINDERIV,BAS+0.1*ABS(BAS))
C              WRITE(*,*) JJ,DERIVATIVE(JJ),II,DELTA/ROTE
             END IF

c            write(*,*) iter_grow,r_ext,
c     &                 clusrx(nclus),clusry(nclus),clusrz(nclus),
c     &                 delta/rote
            END DO   ! do while (DELTA.AND.FLAG_ITER)

            RADIO(NCLUS)=R_EXT
            MASA(NCLUS)=BASMASS

            CLUSRX(NCLUS)=BASX/BASDELTA
            CLUSRY(NCLUS)=BASY/BASDELTA
            CLUSRZ(NCLUS)=BASZ/BASDELTA

*           Find its parent (from higher to lower level)
            busca_pare: DO I1=IR-1,0,-1
             LOW1=SUM(NHALLEV(0:I1-1))+1
             LOW2=SUM(NHALLEV(0:I1))
             X1=CLUSRX(NCLUS)
             Y1=CLUSRY(NCLUS)
             Z1=CLUSRZ(NCLUS)
             DO II=LOW1,LOW2
              X2=CLUSRX(II)
              Y2=CLUSRY(II)
              Z2=CLUSRZ(II)
              BAS=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
              IF (BAS/(1.05*RADIO(II)).LT.1) THEN
               REALCLUS(IFI,NCLUS)=II
               EXIT busca_pare
              END IF
             END DO
            END DO busca_pare

c            WRITE(*,*) IR,CLUSRX(NCLUS),CLUSRY(NCLUS),CLUSRZ(NCLUS),
c     &             RADIO(NCLUS),MASA(NCLUS)*9.1717E18,DELTA/ROTE,
c     &             NCLUS,REALCLUS(IFI,NCLUS)

           END IF
           END IF
           END IF
           END IF
           END IF
           END IF
          END IF ! KK_ENTERO.EQ.1 (FREE HALO) OR .EQ.-1 (SUBSTRUCTURE)
         END DO ! L1=1,NV_GOOD

         DEALLOCATE(DDD,DDDX,DDDY,DDDZ,DDDP)

         WRITE(*,*) 'At level', IR,', no. haloes:', NHALLEV(IR)
         LOW1=SUM(NHALLEV(0:IR-1))+1
         LOW2=SUM(NHALLEV(0:IR))
         WRITE(*,*) '--> Of which, potential substructure:',
     &               COUNT(REALCLUS(IFI,LOW1:LOW2).GT.0)
         WRITE(*,*)

c         CALL OVERLAPPING(IFI,IR,NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
c     &                        PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
c     &                        PATCHRY,PATCHRZ,PARE,NCLUS,MASA,RADIO,
c     &                        CLUSRX,CLUSRY,CLUSRZ,REALCLUS,LEVHAL,
c     &                        NHALLEV,BOUND,CONTRASTEC,RODO,SOLAP,
c     &                        VECINO,NVECI,CR0AMR,CR0AMR11)

         IF (IR.LT.NL) THEN
          LOW1=SUM(NPATCH(0:IR))+1
          LOW2=SUM(NPATCH(0:IR+1))

!$OMP PARALLEL DO SHARED(NPATCH,PATCHNX,PATCHNY,PATCHNZ,CONTA1,LOW1,
!$OMP+                   LOW2),
!$OMP+            PRIVATE(I,N1,N2,N3,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
          DO I=LOW1,LOW2
           N1=PATCHNX(I)
           N2=PATCHNY(I)
           N3=PATCHNZ(I)
           DO IX=1,N1
           DO JY=1,N2
           DO KZ=1,N3
            CONTA1(IX,JY,KZ,I)=1
           END DO
           END DO
           END DO
          END DO

          DXPA=DX/(2.0**(IR+1))
          DYPA=DY/(2.0**(IR+1))
          DZPA=DZ/(2.0**(IR+1))

!$OMP PARALLEL DO SHARED(NPATCH,PATCHNX,PATCHNY,PATCHNZ,PATCHRX,PATCHRY,
!$OMP+                   PATCHRZ,DX,DY,DZ,CLUSRX,CLUSRY,CLUSRZ,RADIO,
!$OMP+                   RX,RY,RZ,CONTA1,NCLUS,IR,LOW1,LOW2,DXPA,DYPA,
!$OMP+                   DZPA),
!$OMP+            PRIVATE(I,N1,N2,N3,X1,X2,Y1,Y2,Z1,Z2,BASX,BASY,BASZ,
!$OMP+                    BAS,IX,JY,KZ,AA,X3,Y3,Z3,X4,Y4,Z4),
!$OMP+            DEFAULT(NONE)
          DO I=LOW1,LOW2
           N1=PATCHNX(I)
           N2=PATCHNY(I)
           N3=PATCHNZ(I)
           X1=PATCHRX(I)-DXPA
           Y1=PATCHRY(I)-DYPA
           Z1=PATCHRZ(I)-DZPA
           X2=X1+N1*DXPA
           Y2=Y1+N2*DYPA
           Z2=Z1+N3*DZPA
           DO II=1,NCLUS
            BASX=CLUSRX(II)
            BASY=CLUSRY(II)
            BASZ=CLUSRZ(II)
            BAS=RADIO(II)
            X3=BASX-BAS
            Y3=BASY-BAS
            Z3=BASZ-BAS
            X4=BASX+BAS
            Y4=BASY+BAS
            Z4=BASZ+BAS
            IF (X1.LE.X4.AND.X3.LE.X2.AND.
     &          Y1.LE.Y4.AND.Y3.LE.Y2.AND.
     &          Z1.LE.Z4.AND.Z3.LE.Z2) THEN
c           WRITE(*,*) BASX,BASY,BASZ,BAS
             DO IX=1,N1
             DO JY=1,N2
             DO KZ=1,N3
              AA=SQRT((RX(IX,I)-BASX)**2 +
     &                (RY(JY,I)-BASY)**2 +
     &                (RZ(KZ,I)-BASZ)**2)
              IF (AA.LE.BAS) CONTA1(IX,JY,KZ,I)=-1
             END DO
             END DO
             END DO
            END IF
           END DO
          END DO

*       And mark the centers of haloes (to avoid identifying haloes at same levels as substructure)
!$OMP PARALLEL DO SHARED(NPATCH,DX,DY,DZ,PATCHNX,PATCHNY,PATCHNZ,
!$OMP+                   PATCHRX,PATCHRY,PATCHRZ,NCLUS,CLUSRX,CLUSRY,
!$OMP+                   CLUSRZ,CONTA1,LOW1,LOW2,DXPA,DYPA,DZPA),
!$OMP+            PRIVATE(I,N1,N2,N3,X1,X2,Y1,Y2,Z1,Z2,II,BASX,BASY,
!$OMP+                    BASZ,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
          DO I=LOW1,LOW2
           N1=PATCHNX(I)
           N2=PATCHNY(I)
           N3=PATCHNZ(I)
           X1=PATCHRX(I)-DXPA
           Y1=PATCHRY(I)-DYPA
           Z1=PATCHRZ(I)-DZPA
           X2=X1+N1*DXPA
           Y2=Y1+N2*DYPA
           Z2=Z1+N3*DZPA
           DO II=1,NCLUS
            BASX=CLUSRX(II)
            BASY=CLUSRY(II)
            BASZ=CLUSRZ(II)
            BASX=(BASX-X1)*(X2-BASX)
            BASY=(BASY-Y1)*(Y2-BASY)
            BASZ=(BASZ-Z1)*(Z2-BASZ)
            IF (BASX.GT.0.AND.BASY.GT.0.AND.BASZ.GT.0) THEN
             BASX=CLUSRX(II)-X1
             BASY=CLUSRY(II)-Y1
             BASZ=CLUSRZ(II)-Z1
             IX=INT(BASX/DXPA)+1
             JY=INT(BASY/DYPA)+1
             KZ=INT(BASZ/DZPA)+1
             DO I1=IX-1,IX+1
             DO J1=JY-1,JY+1
             DO K1=KZ-1,KZ+1
              IF (I1.GE.1.AND.I1.LE.N1.AND.
     &            J1.GE.1.AND.J1.LE.N2.AND.
     &            K1.GE.1.AND.K1.LE.N3) THEN
               CONTA1(I1,J1,K1,I)=-2
              END IF
             END DO
             END DO
             END DO
            END IF
           END DO
          END DO

          ! clean conta1 from overlaps
          CALL CLEAN_OVERLAPS_INT(NL,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                            SOLAP,CONTA1)

         END IF
        END DO !(IR=1,NL)

       END IF !(NL.GT.0)

       DO IR=0,NL
        IF (NHALLEV(IR).GT.0) THEN
         CALL OVERLAPPING(IFI,IR,NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                        PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                        PATCHRY,PATCHRZ,PARE,NCLUS,MASA,RADIO,
     &                        CLUSRX,CLUSRY,CLUSRZ,REALCLUS,LEVHAL,
     &                        NHALLEV,BOUND,CONTRASTEC,RODO,SOLAP,
     &                        VECINO,NVECI,CR0AMR,CR0AMR11,
     &                        VOL_SOLAP_LOW)
        END IF
       END DO

       BASINT=COUNT(REALCLUS(IFI,1:NCLUS).EQ.-1)
       WRITE(*,*) '====> TOTAL NUMBER OF TENTATIVE FREE HALOES',BASINT
       BASINT=COUNT(REALCLUS(IFI,1:NCLUS).GT.0)
       WRITE(*,*) '====> TOTAL NUMBER OF TENTATIVE SUBSTRUCTURES',BASINT

       RETURN
       END

********************************************************************
       SUBROUTINE PATCHES_SPHERE(NPATCH,PATCHRX,PATCHRY,PATCHRZ,
     &                           PATCHNX,PATCHNY,PATCHNZ,
     &                           XCEN,YCEN,ZCEN,REACH,MAXLEVEL,NL,
     &                           RELEVANT_PATCHES,NRELEVANT_PATCHES)
********************************************************************
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER NPATCH(0:NLEVELS)
       REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       REAL XCEN,YCEN,ZCEN,REACH
       INTEGER MAXLEVEL,NL
       INTEGER RELEVANT_PATCHES(NPALEV),NRELEVANT_PATCHES(NLEVELS)

       REAL*4 DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       REAL X1,X2,X3,X4,Y1,Y2,Y3,Y4,Z1,Z2,Z3,Z4,DXPA,DYPA,DZPA
       INTEGER I,J,K,IX,JY,KZ,II,JJ,KK,LOW1,LOW2,IR,N1,N2,N3,LOW3,CONTA

       DO IR=1,NL
        NRELEVANT_PATCHES(IR)=0
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2
         RELEVANT_PATCHES(I)=0
        END DO
       END DO

       X1=XCEN-REACH
       X2=XCEN+REACH
       Y1=YCEN-REACH
       Y2=YCEN+REACH
       Z1=ZCEN-REACH
       Z2=ZCEN+REACH

       LOW3=0
       DO IR=1,MAXLEVEL
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))

        DXPA=DX/2.0**IR
        DYPA=DY/2.0**IR
        DZPA=DZ/2.0**IR

        CONTA=0
        DO I=LOW1,LOW2
         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         X3=PATCHRX(I)-DXPA
         Y3=PATCHRY(I)-DYPA
         Z3=PATCHRZ(I)-DZPA

         X4=X3+N1*DXPA
         Y4=Y3+N2*DYPA
         Z4=Z3+N3*DZPA

         IF (X1.LE.X4.AND.X3.LE.X2) THEN
         IF (Y1.LE.Y4.AND.Y3.LE.Y2) THEN
         IF (Z1.LE.Z4.AND.Z3.LE.Z2) THEN
          CONTA=CONTA+1
          RELEVANT_PATCHES(LOW3+CONTA)=I
         END IF
         END IF
         END IF
        END DO
        NRELEVANT_PATCHES(IR)=CONTA
        LOW3=LOW3+NRELEVANT_PATCHES(IR)
       END DO


       RETURN
       END

********************************************************************
       SUBROUTINE FIND_NEIGHBOURING_PATCHES(IR,NL,CEL,PARE,VECINO,
     &                                      NVECI,NPATCH,VID,NVID)
********************************************************************
*      Finds recursively all the patches touching a given one,
*      and does the same with their parents up to level 1
*      UNUSED NOW / NOT WORKING
********************************************************************
       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

       INTEGER IR,NL,CEL
       INTEGER VECINO(NPALEV,NPALEV),NVECI(NPALEV),NPATCH(0:NLEVELS)
       INTEGER PARE(NPALEV),VID(NLEVELS,NPALEV),NVID(NLEVELS)

       INTEGER IRR,FLAG,NV,NV2,NV3,II,JJ,LOW1,LOW2,KK_ENTERO,IPATCH
       INTEGER CONTAP(NPALEV)

       DO II=1,NL
        NVID(II)=0
        DO JJ=1,NPALEV
         VID(II,JJ)=0
        END DO
       END DO

       IPATCH=CEL

       DO IRR=IR,1,-1

        FLAG=0

        NV=1        !first neighbour: the patch itself
        VID(IRR,NV)=IPATCH

        NV2=0
        NV3=NV

        LOW1=SUM(NPATCH(0:IRR-1))+1
        LOW2=SUM(NPATCH(0:IRR))
        write(*,*) irr,low1,low2
        CONTAP(LOW1:LOW2)=1    !(1 vale, 0 no)
        CONTAP(VID(IRR,NV))=0      !asi este no se puede autocontar
        write(*,*) irr,ipatch,nveci(vid(irr,1))
        DO WHILE (FLAG.EQ.0)
         DO II=NV2+1, NV3
          DO JJ=1, NVECI(VID(IRR,II))
           KK_ENTERO=0
           KK_ENTERO=CONTAP(VECINO(JJ,VID(IRR,II)))
           if (irr.eq.1) write(*,*) ii,jj,VECINO(JJ,VID(IRR,II)),
     &                              VID(IRR,II),kk_entero
           IF (KK_ENTERO.EQ.1) THEN
            NV=NV+1
            VID(IRR,NV)=VECINO(JJ,VID(IRR,II))
            CONTAP(VECINO(JJ,VID(IRR,II)))=0
           END IF
          END DO
         END DO

         IF (NV.EQ.NV3) THEN
          FLAG=1
         ELSE
          NV2=NV3
          NV3=NV
         END IF

        END DO  ! WHILE (FLAG)

        NVID(IRR)=NV

        IF (IRR.GT.1) IPATCH=PARE(IPATCH)

       END DO !IRR=IR,1,-1

       RETURN
       END

********************************************************************
       SUBROUTINE SORT_CELLS(KK_ENTERO,DDD,DDDX,DDDY,DDDZ)
********************************************************************
*      Sorts cells decreasingly in density (DDD)
********************************************************************

       IMPLICIT NONE

       INTEGER KK_ENTERO
       REAL DDD(KK_ENTERO)
       INTEGER DDDX(KK_ENTERO),DDDY(KK_ENTERO),DDDZ(KK_ENTERO)

       REAL DDD2(KK_ENTERO)
       INTEGER DDDX2(KK_ENTERO),DDDY2(KK_ENTERO),DDDZ2(KK_ENTERO)
       INTEGER INDICE2(KK_ENTERO)
       INTEGER I

!$OMP PARALLEL DO SHARED(KK_ENTERO,DDD,DDDX,DDDY,DDDZ,DDD2,DDDX2,
!$OMP+                   DDDY2,DDDZ2),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
       DO I=1,KK_ENTERO
        DDD2(I)=1.0/DDD(I)
        DDDX2(I)=DDDX(I)
        DDDY2(I)=DDDY(I)
        DDDZ2(I)=DDDZ(I)
       END DO

       CALL INDEXX(KK_ENTERO,DDD2,INDICE2)

!$OMP PARALLEL DO SHARED(KK_ENTERO,DDD,DDDX,DDDY,DDDZ,DDD2,DDDX2,
!$OMP+                   DDDY2,DDDZ2,INDICE2),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
       DO I=1,KK_ENTERO
        DDD(I)=1.0/DDD2(INDICE2(I))
        DDDX(I)=DDDX2(INDICE2(I))
        DDDY(I)=DDDY2(INDICE2(I))
        DDDZ(I)=DDDZ2(INDICE2(I))
       END DO

       RETURN
       END

********************************************************************
       SUBROUTINE SORT_CELLS_AMR(KK_ENTERO,DDD,DDDX,DDDY,DDDZ,DDDP)
********************************************************************
*      Sorts cells decreasingly in density (DDD) (including patch number)
********************************************************************

       IMPLICIT NONE

       INTEGER KK_ENTERO
       REAL DDD(KK_ENTERO)
       INTEGER DDDX(KK_ENTERO),DDDY(KK_ENTERO),DDDZ(KK_ENTERO)
       INTEGER DDDP(KK_ENTERO)

       REAL DDD2(KK_ENTERO)
       INTEGER DDDX2(KK_ENTERO),DDDY2(KK_ENTERO),DDDZ2(KK_ENTERO)
       INTEGER INDICE2(KK_ENTERO),DDDP2(KK_ENTERO)
       INTEGER I

!$OMP PARALLEL DO SHARED(KK_ENTERO,DDD,DDDX,DDDY,DDDZ,DDD2,DDDX2,
!$OMP+                   DDDY2,DDDZ2,DDDP,DDDP2),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
       DO I=1,KK_ENTERO
        DDD2(I)=1.0/DDD(I)
        DDDX2(I)=DDDX(I)
        DDDY2(I)=DDDY(I)
        DDDZ2(I)=DDDZ(I)
        DDDP2(I)=DDDP(I)
       END DO

       CALL INDEXX(KK_ENTERO,DDD2,INDICE2)

!$OMP PARALLEL DO SHARED(KK_ENTERO,DDD,DDDX,DDDY,DDDZ,DDD2,DDDX2,
!$OMP+                   DDDY2,DDDZ2,INDICE2,DDDP,DDDP2),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
       DO I=1,KK_ENTERO
        DDD(I)=1.0/DDD2(INDICE2(I))
        DDDX(I)=DDDX2(INDICE2(I))
        DDDY(I)=DDDY2(INDICE2(I))
        DDDZ(I)=DDDZ2(INDICE2(I))
        DDDP(I)=DDDP2(INDICE2(I))
       END DO

       RETURN
       END
