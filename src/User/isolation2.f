C     Photon isolation

      logical function isolation2(p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),Ris,R35,R36,R45,R46,pt
      double precision pt5,pt6,ptm,R,Etmax,Et3,Et4

      isolation2=.true.

c     Photons are isolated if total transverse energy in a cone
c     of radius Ris is smaller than Etmax

      Ris=0.3d0
      Etmax=6d0

      pt5=pt(5,p)
      pt6=pt(6,p)
      ptm=max(pt5,pt6)
 

      if(pt5+pt6.lt.Etmax) return
      

C     One parton

      if(pt6.eq.0d0) then
       if(pt5.lt.Etmax.or.(r(p,3,5).gt.Ris.and.r(p,4,5).gt.Ris))then
        return
       else
        isolation2=.false.
        return
       endif
      endif


CC    Two partons

C     Transverse energy around photon 3

      if((r(p,3,5).lt.Ris).and.(r(p,3,6).lt.Ris)) then
       Et3=pt5+pt6
      elseif ((r(p,3,5).lt.Ris).and.(r(p,3,6).gt.Ris)) then
       Et3=pt5
      elseif ((r(p,3,5).gt.Ris).and.(r(p,3,6).lt.Ris)) then
       Et3=pt6
      else
       Et3=0d0
      endif

C     Transverse energy around photon 4

      if((r(p,4,5).lt.Ris).and.(r(p,4,6).lt.Ris)) then
       Et4=pt5+pt6
      elseif ((r(p,4,5).lt.Ris).and.(r(p,4,6).gt.Ris)) then
       Et4=pt5
      elseif ((r(p,4,5).gt.Ris).and.(r(p,4,6).lt.Ris)) then
       Et4=pt6
      else
       Et4=0d0
      endif

      if(Et3.lt.Etmax.and.Et4.lt.Etmax) then
       return
      else
       isolation2=.false.
       return
      endif

      write(6,*)'Error in isolation subroutine'
      stop

      return
      end
