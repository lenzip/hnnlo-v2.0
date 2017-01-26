c--- Generates 2->4 phase space with (3456) a Breit-Wigner around
c--- the Higgs mass
      subroutine gen4h(r,p,wt4,*)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'masses.f'
      include 'phasemin.f'
      include 'limits.f'
      integer nu
      double precision r(mxdim)
      double precision wt4,p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      double precision p(mxpart,4),sqrts,rtshat
      double precision pswt,xjac
      double precision xx(2),s3456,wt3456,ymax,yave
      common/energy/sqrts
      common/x1x2/xx

      wt4=0d0

      call breitw(r(9),wsqmin,wsqmax,hmass,hwidth,s3456,wt3456)
            
      rtshat=dsqrt(s3456)
      ymax=dlog(sqrts/rtshat)
      yave=ymax*(two*r(10)-1d0)
      xjac=two*ymax*wt3456
           
      xx(1)=rtshat/sqrts*exp(+yave)
      xx(2)=rtshat/sqrts*exp(-yave)

      if   ((xx(1) .gt. 1d0) 
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) then
      write(6,*) 'problems with xx(1),xx(2) in gen4h',xx(1),xx(2)  
      return 1 
      endif

      p1(4)=-0.5d0*xx(1)*sqrts
      p1(1)=0d0
      p1(2)=0d0
      p1(3)=-0.5d0*xx(1)*sqrts
      
      p2(4)=-0.5d0*xx(2)*sqrts
      p2(1)=0d0
      p2(2)=0d0
      p2(3)=+0.5d0*xx(2)*sqrts

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

      wt4=xjac*pswt/sqrts**2
      
      return

 999  return 1
      end
