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
       real*4  c,g,h,s,sm,t,tau,theta,tresh,sum_elements

       allocate(B(1:100))   !,stat=ialloc)
       allocate(Z(1:100))   !,stat=ialloc)

       do ip=1, N
         B(ip)=A(ip,ip)
         D(ip)=B(ip)
         Z(ip)=0.d0
       end do
       NROT=0

       sum_elements=0.0
       do ip=1,N
       do iq=ip,N
        sum_elements=sum_elements+abs(A(ip,iq))
       end do
       end do

       do i=1, 50
         sm=0.d0
         do ip=1, N-1           !sum off-diagonal elements
            do iq=ip+1, N
               sm=sm+ABS(A(ip,iq))
            end do
         end do
         if(sm.lt.1.e-4*sum_elements) return    !normal return
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

**********************************************************************
      SUBROUTINE indexx(n,arr,indx)
**********************************************************************
*     Sorts a list of n elements (arr), without modifying it,
*     increasingly. The sorted indices are returned in the variable
*     indx.
**********************************************************************


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
