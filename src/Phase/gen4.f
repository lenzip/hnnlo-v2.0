      subroutine gen4(r,p,wt4,*)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'debug.f'
      include 'process.f'
      include 'phasemin.f'
      integer nu
      double precision r(mxdim)
      double precision wt4,p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      double precision p(mxpart,4)
      double precision pswt,xjac,p1ext(4),p2ext(4)
      double precision xx(2),tau,x1mx2,surd
      double precision lntaum
      common/pext/p1ext,p2ext
      common/x1x2/xx


      wt4=0d0

      lntaum=dlog(taumin)
      tau=dexp(lntaum*(one-r(9)))
      xjac=-lntaum*tau

c      tau=(one-taumin)*r(9)**2+taumin
c      xjac=2*r(9)*(one-taumin)

      x1mx2=two*r(10)-one
      surd=dsqrt(x1mx2**2+four*tau) 
           
      xx(1)=half*(+x1mx2+surd)
      xx(2)=half*(-x1mx2+surd)

      xjac=xjac*two/surd

      if   ((xx(1) .gt. 1d0) 
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) return 1 

      do nu=1,4
      p1(nu)=xx(1)*p1ext(nu)
      p2(nu)=xx(2)*p2ext(nu)
      enddo

  
      call phase4(r,p1,p2,p3,p4,p5,p6,pswt,*999) 
      

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      p(7,nu)=0d0

      enddo 

      wt4=xjac*pswt
      
      if (debug) write(6,*) 'wt4 in gen4',wt4
      return

 999  return 1
      end

