      subroutine genBORN4(q2,shat,r,p,wt,*)
c----generate phase space weight and vectors p(i,4) for i=1,2,3,4,5,6
C    q2 and shat are input
C    
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'phasemin.f'
      integer nu

      double precision r(mxdim),sqrts,wt,wt34,
     . p(mxpart,4),p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      double precision pswt,xjac,xx(2),tau,tau0,y,q2,yq,shat

      common/energy/sqrts
      common/xx0/xx

      wt=0d0

      tau=shat/sqrts**2
      tau0=q2/sqrts**2
      y=0.5d0*dlog(tau0)*(1d0-2d0*r(10))
      xjac=dabs(dlog(tau0))

      xx(1)=dsqrt(tau0)*dexp(+y)
      xx(2)=dsqrt(tau0)*dexp(-y)

c---if x's out of normal range alternative return
      if   ((xx(1) .gt. 1d0) 
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) return 99

      p1(4)=-xx(1)*sqrts*half
      p1(1)=zip
      p1(2)=zip
      p1(3)=-xx(1)*sqrts*half

      p2(4)=-xx(2)*sqrts*half
      p2(1)=zip
      p2(2)=zip
      p2(3)=+xx(2)*sqrts*half


      call phase4(r,p1,p2,p3,p4,p5,p6,pswt,*99) 

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      enddo 
      wt=xjac*pswt
      return
 99   continue
      wt=0d0

      return
      end
