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
       SUBROUTINE OVERLAPING(IFI,IR,NL,REF,ESP,BOUND,CONTA,CONTRASTEC,
     &                       NSHELL,RODO,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                       PATCHRX,PATCHRY,PATCHRZ,NX,NY,NZ,
     &                       NCLUS,MASA,RADIO,CLUSRX,CLUSRY,CLUSRZ,
     &                       REALCLUS,NSOLAP,SOLAPA,NHALLEV)
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
       SUBROUTINE HALOFIND_GRID(IFI,NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                          PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,NCLUS,
     &                          MASA,RADIO,CLUSRX,CLUSRY,CLUSRZ,
     &                          REALCLUS,LEVHAL,NSOLAP,SOLAPA,NHALLEV,
     &                          BOUND,CONTRASTEC,RODO)
********************************************************************
*      Pipeline for tentative halo finding over the grid
********************************************************************

       IMPLICIT NONE
       INCLUDE 'input_files/asohf_parameters.dat'

*      I/O DATA
       INTEGER IFI,NL,NX,NY,NZ
       INTEGER NPATCH(0:NLEVELS)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
       INTEGER NCLUS
       REAL MASA(MAXNCLUS),RADIO(MAXNCLUS)
       REAL*4 CLUSRX(MAXNCLUS),CLUSRY(MAXNCLUS),CLUSRZ(MAXNCLUS)
       INTEGER SOLAPA(MAXNCLUS,NMAXSOLAP),NSOLAP(MAXNCLUS)
       INTEGER REALCLUS(MAXITER,MAXNCLUS),LEVHAL(MAXNCLUS)
       INTEGER NHALLEV(0:NLEVELS)
       REAL BOUND,CONTRASTEC,RODO

*      GLOBAL VARIABLES
       REAL*4 DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       REAL*4  RADX(0:NMAX+1),RADY(0:NMAY+1),RADZ(0:NMAZ+1)
       COMMON /GRID/ RADX,RADY,RADZ

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
       INTEGER CONTA2(NAMRX,NAMRY,NAMRZ,NPALEV)
       REAL UBAS1(NMAX,NMAY,NMAZ)

       INTEGER IR,NSHELL,IX,JY,KZ,I,J,K,II,IPATCH,ICEN(3),NV_GOOD
       INTEGER L1,NX1,NX2,NY1,NY2,NZ1,NZ2,KK_ENTERO,ITER_GROW
       REAL PRUEBAX,PRUEBAY,PRUEBAZ,RMIN,BASMASS_SHELL,BASMASS,DELTA
       REAL REF,ESP,ESP_LOG,BAS,KK_REAL,RSHELL,R_INT,R_EXT,RANT
       REAL BASDELTA,AA,PI,VOLCELL,BASX,BASY,BASZ,BASVOL
       REAL*4, ALLOCATABLE::DDD(:)
       INTEGER, ALLOCATABLE::DDDX(:),DDDY(:),DDDZ(:)

**************************************************************
*      NIVEL BASE!!
**************************************************************

       PI=DACOS(-1.D0)

       IR=0
       REF=0.2*DX !THIS WILL BE GOTTEN RID OF
       ESP=0.2*DX
       ESP_LOG=1.05
       NSHELL=50 !THIS WILL BE GOTTEN RID OF
       WRITE(*,*) 'RESOLUTION (IR,DX)=',IR,DX

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
       WRITE(*,*)'ESTIMATION_0:',IR,KK_ENTERO
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
         WRITE(*,*) U1(ICEN(1),ICEN(2),ICEN(3))
         NCLUS=NCLUS+1
         REALCLUS(IFI,NCLUS)=-1
         LEVHAL(NCLUS)=IR
         NHALLEV(IR)=NHALLEV(IR)+1

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
          WRITE(*,*) DELTA/ROTE, II
         END DO   ! do while (DELTA)

         RADIO(NCLUS)=R_EXT
         MASA(NCLUS)=BASMASS

         CLUSRX(NCLUS)=BASX/BASDELTA
         CLUSRY(NCLUS)=BASY/BASDELTA
         CLUSRZ(NCLUS)=BASZ/BASDELTA

         WRITE(*,*) CLUSRX(NCLUS),CLUSRY(NCLUS),CLUSRZ(NCLUS),
     &             RADIO(NCLUS),MASA(NCLUS)*9.1717E18

        END IF ! KK_ENTERO.EQ.0
       END DO ! I=1,NV_GOOD

       DEALLOCATE(DDD,DDDX,DDDY,DDDZ)

****************************************************
*      CORRECION DE SOLAPES:
*      VAMOS A VER QUE CUMULOS SOLAPAN EN IR
****************************************************

       write(*,*) IFI,IR,NL,REF,ESP,BOUND,CONTRASTEC,
     &                 NSHELL,RODO,NPATCH,NX,NY,NZ,
     &                 NCLUS,
     &                 NHALLEV
       CALL OVERLAPING(IFI,IR,NL,REF,ESP,BOUND,CONTA,CONTRASTEC,
     &                 NSHELL,RODO,NPATCH,PATCHNX,PATCHNY,PATCHNZ,
     &                 PATCHRX,PATCHRY,PATCHRZ,NX,NY,NZ,
     &                 NCLUS,MASA,RADIO,CLUSRX,CLUSRY,CLUSRZ,
     &                 REALCLUS,NSOLAP,SOLAPA,NHALLEV)

****************************************************
*      FIN CORRECION DE SOLAPES EN IR
****************************************************

       WRITE(*,*)' FIN NIVEL BASE', 0



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
