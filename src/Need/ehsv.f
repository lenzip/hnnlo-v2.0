      double complex function ehsva4(s,t,u,mfsq1,mfsq2)
C     ehsv:EqnA.8
      implicit none
      double precision s,t,u,mfsq1,mfsq2
      double complex ehsvb4
      ehsva4=ehsvb4(s,t,u,mfsq1,mfsq2)
     .      +ehsvb4(u,s,t,mfsq1,mfsq2)
     .      +ehsvb4(t,u,s,mfsq1,mfsq2)
      return 
      end

      double complex function ehsva2(s,t,u,mfsq1,mfsq2)
C     ehsv:EqnA.9
      implicit none
      double precision s,t,u,mfsq1,mfsq2
      double complex ehsvb2
      ehsva2=ehsvb2(s,t,u,mfsq1,mfsq2)
     .      +ehsvb2(s,u,t,mfsq1,mfsq2)
      return 
      end
CH--  change of ehsvb4 (adding dependence on the fermion mass)
      double complex function ehsvb4f(s,t,u,mfsq)
      implicit none
C     ehsv:EqnA.10
      include 'masses.f'
      double precision hmass2,s,t,u,mfsq
      double complex w2f,w3f
      hmass2=hmass**2
      ehsvb4f=mfsq/hmass2*(-2d0/3d0
     . +(mfsq/hmass2-0.25d0)*(w2f(t,mfsq)-w2f(hmass2,mfsq)
     . +w3f(s,t,u,hmass2,mfsq)))
      return 
      end
CH--  
       double complex FUNCTION ehsvb4(s,t,u,mfsq1,mfsq2)
       implicit none
       double complex ehsvb4f
       double precision s,t,u,mfsq1,mfsq2
       integer approxim
       common/flag/approxim
       ehsvb4=(approxim-1)*ehsvb4f(s,t,u,mfsq1)+ehsvb4f(s,t,u,mfsq2)
ch       ehsvb4=ehsvb4f(s,t,u,mfsq2)
       return
       end
CH--  change of ehsvb2 (adding dependence on the fermion mass)
      double complex function ehsvb2f(s,t,u,mfsq)
C     ehsv:EqnA.11
      implicit none
      include 'masses.f'
      double precision hmass2,s,t,u,mfsq
      double complex w1f,w2f,w3f
      hmass2=hmass**2
      ehsvb2f=mfsq/hmass2**2*(s*(u-s)/(s+u)
     . +2d0*u*t*(u+2d0*s)/(s+u)**2*(w1f(t,mfsq)-w1f(hmass2,mfsq))
     . +(mfsq-0.25d0*s)
     . *(0.5d0*w2f(s,mfsq)+0.5d0*w2f(hmass2,mfsq)-w2f(t,mfsq)
     . +w3f(s,t,u,hmass2,mfsq))
     . +s**2*(2d0*mfsq/(s+u)**2-0.5d0/(s+u))*
     . (w2f(t,mfsq)-w2f(hmass2,mfsq))
     . +0.5d0*u*t/s*(w2f(hmass2,mfsq)-2d0*w2f(t,mfsq))
     . +0.125d0*(s-12d0*mfsq-4d0*u*t/s)*w3f(t,s,u,hmass2,mfsq))
      return 
      end
CH--
       double complex FUNCTION ehsvb2(s,t,u,mfsq1,mfsq2)
       implicit none
       double complex ehsvb2f
       double precision s,t,u,mfsq1,mfsq2
       integer approxim
       common/flag/approxim
       ehsvb2=(approxim-1)*ehsvb2f(s,t,u,mfsq1)+ehsvb2f(s,t,u,mfsq2)
ch       ehsvb2=ehsvb2f(s,t,u,mfsq2)
       return
       end
CH--  change of ehsva5 (adding dependence on the fermion mass)
      double complex function ehsva5f(s,t,u,mfsq)
C     ehsv:EqnA.14
      implicit none
      include 'masses.f'
      double precision hmass2,s,t,u,mfsq
      double complex w1f,w2f
      hmass2=hmass**2
      ehsva5f=mfsq/hmass2*(4d0+4d0*s/(u+t)*(w1f(s,mfsq)
     . -w1f(hmass2,mfsq))
     . +(1d0-4d0*mfsq/(u+t))*(w2f(s,mfsq)-w2f(hmass2,mfsq)))
      return 
      end
CH--
       double complex FUNCTION ehsva5(s,t,u,mfsq1,mfsq2)
       implicit none
       double complex ehsva5f
       double precision s,t,u,mfsq1,mfsq2
       integer approxim
       common/flag/approxim
       ehsva5=(approxim-1)*ehsva5f(s,t,u,mfsq1)+ehsva5f(s,t,u,mfsq2)
ch       ehsva5=ehsva5f(s,t,u,mfsq2)
       return
       end

CH--  change of w1 (adding dependence on the fermion mass)
      double complex function w1f(s,mfsq)
C     ehsv:EqnA.19
      implicit none
      include 'constants.f'
      include 'masses.f'
      double precision rat,s,temp,acosh,asinh,mfsq
      rat=4d0*mfsq/s
      temp=dsqrt(dabs(1d0/rat))
      if (rat .lt. 0d0) then
          w1f=2d0*dsqrt(1d0-rat)*asinh(temp)
      elseif (rat .gt. 1d0) then
          w1f=2d0*dsqrt(rat-1d0)*asin(temp)
      else 
          temp=2d0*acosh(temp)
          w1f=dsqrt(1d0-rat)*dcmplx(temp,-pi)
       endif
      return
      end
CH--  change of w2 (adding dependence on the fermion mass)
      double complex function w2f(s,mfsq)
C     ehsv:EqnA.20
      implicit none
      include 'constants.f'
      include 'masses.f'
      double precision rat,s,tempr,tempi,acosh,asinh,mfsq
      rat=s/(4d0*mfsq)
      tempr=dsqrt(dabs(rat))
      if (rat .lt. 0d0) then
          tempr=asinh(tempr)
          w2f=4d0*tempr**2
      elseif (rat .gt. 1d0) then
          tempr=acosh(tempr)
          tempi=-4d0*tempr*pi
          tempr=+4d0*tempr**2-pi**2
          w2f=dcmplx(tempr,tempi)
      else 
          tempr=asin(tempr)
          w2f=-4d0*tempr**2
      endif
      return
      end
CH--  change of w3 (adding dependence on the fermion mass)
      double complex function w3f(s,t,u,varg,mfsq)
C     ehsv:EqnA.17
      implicit none
      double complex i3f,w2f
      double precision s,t,u,mfsq,varg
      w3f=i3f(s,t,u,varg,mfsq)-i3f(s,t,u,s,mfsq)
     .                        -i3f(s,t,u,u,mfsq)
      return
      end

CH--  change of i3 (adding dependence on the fermion mass)
      double complex function i3f(s,t,u,varg,mfsq)
C     ehsv:EqnA.21
      implicit none
      include 'constants.f'
      include 'masses.f'
      double precision s,t,u,varg,rat,al,be,ga,r,theta,phi,mfsq,DDILOG
      double precision arg1,arg2,arg3,arg4,targ
      double complex li2,zth,zph
      rat=4d0*mfsq/varg
      if (rat .lt. 0d0) then
           be=0.5d0*(1d0+dsqrt(1d0+4d0*t*mfsq/(u*s)))
           if (be .eq. 1d0)then
             be=1.00000001d0
           endif
           ga=0.5d0*(1d0+dsqrt(1d0-rat))
           arg1=ga/(ga+be-1d0)
           arg2=(ga-1d0)/(ga+be-1d0)
           arg3=(be-ga)/be
           arg4=(be-ga)/(be-1d0)
           i3f=2d0/(2d0*be-1d0)
     .     *(-DDILOG(arg1)+DDILOG(arg2)+DDILOG(arg3)-DDILOG(arg4)
     .     +0.5d0*(dlog(be)**2-dlog(be-1d0)**2)
     .     +dlog(ga)*dlog((ga+be-1d0)/be)
     .     +dlog(ga-1d0)*dlog((be-1d0)/(ga+be-1d0)))
      elseif (rat .gt. 1d0) then
           be=0.5d0*(1d0+dsqrt(1d0+4d0*t*mfsq/(u*s)))
           al=dsqrt(rat-1d0)
           r=dsqrt((al**2+1d0)/(al**2+(2d0*be-1d0)**2))
CC
           targ=r*(al**2+2d0*be-1d0)/(1d0+al**2)
           if(targ.le.1d0) then
            phi=acos(targ)
           else
            phi=0d0
           endif
C           phi=acos(r*(al**2+2d0*be-1d0)/(1d0+al**2))
           theta=acos(r*(al**2-2d0*be+1d0)/(1d0+al**2))
           zth=r*dcmplx(cos(theta),sin(theta))
           zph=r*dcmplx(cos(phi),sin(phi))
           i3f=2d0/(2d0*be-1d0)
     .     *(2d0*dble(li2(zth))-2d0*dble(li2(zph))
     .     +(phi-theta)*(phi+theta-pi))
      else
           be=0.5d0*(1d0+dsqrt(1d0+4d0*t*mfsq/(u*s)))
           ga=0.5d0*(1d0+dsqrt(1d0-rat))
           arg1=ga/(ga+be-1d0)
           arg2=(ga-1d0)/(ga+be-1d0)
           arg3=ga/(ga-be)
           arg4=(ga-1d0)/(ga-be)
           i3f=2d0/(2d0*be-1d0)
     .     *(-DDILOG(arg1)+DDILOG(arg2)+DDILOG(arg3)-DDILOG(arg4)
     .     +dlog(ga/(1d0-ga))*dlog((ga+be-1d0)/(be-ga))
     .     -im*pi*dlog((ga+be-1d0)/(be-ga)))
      endif

      return
      end

      double precision function acosh(y)
      implicit none
      double precision y      
      acosh=dlog(y+dsqrt(y**2-1d0))
      return
      end

      double precision function asinh(y)
      implicit none
      double precision y      
      asinh=dlog(y+dsqrt(y**2+1d0))
      return
      end

