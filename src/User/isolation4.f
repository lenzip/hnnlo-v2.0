C     Lepton isolation for H->ZZ->4l

      logical function isolation4(p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),Ris,pt,eps,Eti
      double precision pt7,pt8,R
      integer i

      isolation4=.true.

c     Leptons are isolated if total transverse energy in a cone
c     of radius Ris is smaller than eps*pt

      Ris=0.2d0
      eps=0.05d0

      pt7=pt(7,p)
      pt8=pt(8,p)
      
C     Loop over leptons

      do i=3,6
      
      Eti=0d0

C     Transverse energy around lepton i

      if(pt8.eq.0d0) then
       if(r(p,i,7).lt.Ris)Eti=pt7
      else
       if((r(p,i,7).lt.Ris).and.(r(p,i,8).lt.Ris)) then
        Eti=pt7+pt8
       elseif ((r(p,i,7).lt.Ris).and.(r(p,i,8).gt.Ris)) then
        Eti=pt7
       elseif ((r(p,i,7).gt.Ris).and.(r(p,i,8).lt.Ris)) then
        Eti=pt8
       endif
      endif


      if(Eti.gt.eps*pt(i,p))then
       isolation4=.false.
       return
      endif

      enddo


      return
      end
