      logical function cuts(pjet,njets)
      implicit none
      include 'constants.f'
      include 'masses.f'
      integer i,j,k,njets
      double precision pjet(mxpart,4),etvec(4)
      double precision pt,etarap
CC
      double precision pt3,pt4,eta3,eta4,pt34,y34,yraptwo
      double precision pt5,eta5,eta6,ptmiss,ptjetmax,m45,pt45,mth
      double precision cosphi45,deltaphi,cosphiLLMet
      double precision ppt(1:4),pto(1:4),mz(1:4),dmz(1:4),tmp(1:2)

      integer l(4),m(4),i1,i2,i3,i4

      logical isol
      common/isol/isol

      integer higgsdec,ndec
      common/higgsdec/higgsdec,ndec


CC
      cuts=.false.


      if(higgsdec.eq.1) then

CC    H->2gamma


CC    Decide here if photons should be isolated (default is not)

      
      isol=.false.


CC    Insert here cuts

      pt3=dsqrt(pjet(3,1)**2+pjet(3,2)**2)
      pt4=dsqrt(pjet(4,1)**2+pjet(4,2)**2)      
      eta3=etarap(3,pjet)
      eta4=etarap(4,pjet)


      pt34=dsqrt((pjet(3,1)+pjet(4,1))**2+(pjet(3,2)+pjet(4,2))**2)

CCC
ch      if(pt34.lt.100d0) cuts=.true.
CCC

      y34=yraptwo(3,4,pjet)

CC     Acceptance cuts on the photons

c       if((min(pt3,pt4).lt.35d0).or.(max(pt3,pt4).lt.40d0).
c     #   or.(dabs(eta3).gt.2.5d0).or.(dabs(eta4).gt.2.5d0)) then 
c        cuts=.true.
c       endif



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      elseif(higgsdec.eq.2) then

CC    H->2W->lnulnu

CC    Selection cuts from hep-ph/0402218

CC    Leptons pt's and eta

      pt4=dsqrt(pjet(4,1)**2+pjet(4,2)**2)
      pt5=dsqrt(pjet(5,1)**2+pjet(5,2)**2)      
      eta4=etarap(4,pjet)
      eta5=etarap(5,pjet)
      pt45=dsqrt((pjet(4,1)+pjet(5,1))**2+(pjet(4,2)+pjet(5,2))**2)
C   leptons in acceptance
      if((dabs(eta4).gt.2.5).or.(dabs(eta5).gt.2.5)) cuts=.true.      
c pt > 20
      if(min(pt4,pt5).lt.20d0) cuts=.true.
c ptll > 30      
      if (pt45.lt.30) cuts=.true.

c      if((max(pt4,pt5).lt.35d0).or.(max(pt4,pt5).gt.50d0)) cuts=.true.


CC    Ptmiss > 20

       ptmiss=dsqrt((pjet(3,1)+pjet(6,1))**2+(pjet(3,2)+pjet(6,2))**2)

      if(ptmiss.lt.20d0) cuts=.true.


CC    mll > 50

      m45=dsqrt((pjet(4,4)+pjet(5,4))**2
     &    -(pjet(4,1)+pjet(5,1))**2
     &    -(pjet(4,2)+pjet(5,2))**2
     &    -(pjet(4,3)+pjet(5,3))**2)


      if(m45.lt.50d0) cuts=.true. 

CC     Jet Veto

       ptjetmax=max(pt(7,pjet),pt(8,pjet))

c      if(ptjetmax.gt.30d0) cuts=.true. 


C     azimuthal angle of (4,5)

      cosphi45=(pjet(4,1)* pjet(5,1)+ 
     &    pjet(4,2)* pjet(5,2))/(pt(4,pjet)*pt(5,pjet))

C     deltaphi

      deltaphi=dabs(dacos(cosphi45))


C Convert to degrees

      deltaphi=deltaphi*180/3.141592653589793d0

c      if(deltaphi.gt.45d0) cuts=.true.  


c mth > 60
      cosphiLLMet=((pjet(3,1)+pjet(6,1))*(pjet(4,1)+pjet(5,1)) + 
     &             (pjet(3,2)+pjet(6,2))*(pjet(4,2)+pjet(5,2)) ) /
     &             (ptmiss*pt45)
      mth=dsqrt(2*ptmiss*pt45*cosphiLLMet)

      if(mth.lt.60d0) cuts=.true.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CC    H->2Z->4l

      elseif(higgsdec/10.eq.3) then

C     Leptons should be isolated

c      isol=.true.

C     Standard rapidity cut

      eta3=dabs(etarap(3,pjet))
      eta4=dabs(etarap(4,pjet))
      eta5=dabs(etarap(5,pjet))
      eta6=dabs(etarap(6,pjet))

c      if(max(eta3,eta4,eta5,eta6).gt.2.5d0) cuts=.true.
      

CC    Order lepton pt's

      do i=1,4
      ppt(i)=dsqrt(pjet(2+i,1)**2+pjet(2+i,2)**2)
      enddo

      pto(1)=max(ppt(1),ppt(2),ppt(3),ppt(4))
      pto(4)=min(ppt(1),ppt(2),ppt(3),ppt(4))

      j=0

      do i=1,4
       if((ppt(i).ne.pto(1)).and.(ppt(i).ne.pto(4))) then
       j=j+1
       tmp(j)=ppt(i)
       endif
      enddo

      pto(2)=max(tmp(1),tmp(2))
      pto(3)=min(tmp(1),tmp(2))


c      if((pto(1).lt.30d0).or.(pto(2).lt.25d0).or.(pto(3).lt.15d0)
c     & .or.(pto(4).lt.7d0)) cuts=.true.


C Cuts on ll invariant masses

C Only in the case H->2Z->4e
       
      if(higgsdec.eq.32) then

      l(1)=3
      m(1)=4

      l(2)=5
      m(2)=6

      l(3)=3
      m(3)=6

      l(4)=4
      m(4)=5

      do i=1,4

      mz(i)=dsqrt((pjet(l(i),4)+pjet(m(i),4))**2
     &     -(pjet(l(i),1)+pjet(m(i),1))**2
     &     -(pjet(l(i),2)+pjet(m(i),2))**2
     &     -(pjet(l(i),3)+pjet(m(i),3))**2)


      dmz(i)=dabs(mz(i)-zmass)

      enddo



C     Find the closest and next-to closest to mz and put cuts


      do i=1,4       
       if(dmz(i).eq.min(dmz(1),dmz(2),dmz(3),dmz(4))) then
        i1=i
       elseif(dmz(i).eq.max(dmz(1),dmz(2),dmz(3),dmz(4))) then
        i4=i
       else 
        i2=i
       endif
      enddo

     
      i3=10-i1-i2-i4

      if(dmz(i2).gt.dmz(i3))then
       i2=i3 
       i3=10-i1-i2-i4
      endif


c      if(mz(i1).gt.101d0.or.mz(i1).lt.81d0) cuts=.true.
      
c      if(mz(i2).gt.110d0.or.mz(i2).lt.40d0) cuts=.true.



      endif




      endif


      return     
      end
 
 
 
 
