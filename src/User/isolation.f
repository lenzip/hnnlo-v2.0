C     Check if leptons (photons) are isolated

C     higgsdec=1: handled with isolation2
C     higgsdec=2: isolation not implemented
C     higgsdec=31,32: handled with isolation4

      logical function isolation(p)
      implicit none

      include 'constants.f'
      double precision p(mxpart,4)
      integer higgsdec,ndec
      common/higgsdec/higgsdec,ndec

      logical isolation2,isolation4

      if(higgsdec.eq.1) then
       isolation=isolation2(p)
      elseif(higgsdec.eq.2) then
       isolation=.true.
      elseif(higgsdec/10.eq.3) then
       isolation=isolation4(p)
      endif

      return
      end
