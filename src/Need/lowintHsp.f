C Full H2st included

C Scale dependence added up to NNLO

C Effect of spin correlations included

      double precision function lowintHst(r,wgt)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'limits.f'
      include 'npart.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'noglue.f'
      include 'process.f'
      include 'efficiency.f'
      include 'phasemin.f'

C
      include 'qcdcouple.f'
      include 'rescoeff.f'

c --- DSW.
      integer pflav,pbarflav,approxim
c --- To use VEGAS random number sequence :
      double precision ran2
      integer ih1,ih2,j,k,nvec,sgnj,sgnk,flgq
      double precision r(mxdim),W,sqrts,xmsq,val,
     . fx10(-nf:nf),fx20(-nf:nf),p(mxpart,4),pjet(mxpart,4),
     . pswt,rscalestart,fscalestart
      double precision wgt,msqc(-nf:nf,-nf:nf),msqb,m3,m4,m5,xmsqjk
      double precision xx(2),flux,vol,vol_mass,vol3_mass,BrnRat
      logical bin,first,includedipole
CC
      logical cuts
      double precision x1p,x2p,fx1p(-nf:nf),fx2p(-nf:nf)
      double precision asopi,z1,z2,alfa,beta,cut,diff
      double precision tdelta,tH1st,tH1stF,xx10,xx20,tH2st,th1stH
      double precision tgaga,tcga,tgamma2
      double precision diff10,diff20,diffc10,diffc20,diffg10,diffg20
      double precision diff1f,diff2f,diffg1f,diffg2f,diffc1f,diffc2f
      double precision Pggreg,D0int,D1int,Cgq,Pgq,LF,LR,H2st,H2stgq
      double precision H2ggREG,dot,q2,Ggq,Ggg,spgq1,spgq2,tH2sp
      double precision Pggggreg,Pgggq,Pgqqg,Pgqqq
      double precision CgqPqq,CgqPqg,P2gg,P2gq
CH
      double precision mtop,mbot,sigmatb,sigmat
CH
C
      double precision beta1,H2ggdelta,H2ggD0
      common/Hstcoeff/beta1,H2ggdelta,H2ggD0

      external Pggreg,D0int,Cgq,Pgq,H2stgq,H2ggREG,P2gg,P2gq,
     &         Pggggreg,Pgqqg,Pgqqq,Pgggq,CgqPqq,CgqPqg

      integer order
      common/nnlo/order
CC
      integer jets,higgsdec,ndec
      common/parts_int/jets
      common/higgsdec/higgsdec,ndec
CC
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/x1x2/xx
      common/BrnRat/BrnRat
      common/flag/approxim
CH
      common/Born/sigmatb,sigmat
      common/tbmass/mtop,mbot
CH
      data p/48*0d0/
      data first/.true./
      save first,rscalestart,fscalestart
      
      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale
      endif

      ntotshot=ntotshot+1
      lowintHst=0d0

C     The number of jets is zero for this peice

      jets=0      

C
      W=sqrts**2


      if (higgsdec.eq.1) then   
       npart=2
       call gen2(r,p,pswt,*999)
      else
       npart=4
       call gen4h(r,p,pswt,*999)
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      nvec=npart+2
      call dotem(nvec,p,s)

      call masscuts(s,*999)
      
C     Here there is no need to cut sij

C      call smalls(s,npart,*999)                                                 

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out

      if(cuts(p,0) .eqv. .true.) then
        goto 999
      endif
      
      
      xx(1)=-2d0*p(1,4)/sqrts
      xx(2)=-2d0*p(2,4)/sqrts


c--- Calculate the required gg->H matrix element      
      if(higgsdec.eq.1 .and. approxim.eq.0) then
        call gg_h(p,msqc)
        q2=2*dot(p,3,4)
       elseif(higgsdec.eq.1 .and. approxim.ne.0)then
        call gg_hexact(p,msqc)
        q2=2*dot(p,3,4)
      elseif(higgsdec.eq.2 .and. approxim.eq.0) then
       call qqb_hww(p,msqc)
        q2=2*(dot(p,3,4)+dot(p,3,5)+dot(p,3,6)
     &   +dot(p,4,5)+dot(p,4,6)+dot(p,5,6))
      elseif(higgsdec.eq.2 .and. approxim.ne.0) then
       call qqb_hwwexact(p,msqc)
       q2=2*(dot(p,3,4)+dot(p,3,5)+dot(p,3,6)
     &   +dot(p,4,5)+dot(p,4,6)+dot(p,5,6))
      elseif(approxim.eq.0)then
       call qqb_hzz(p,msqc)
       q2=2*(dot(p,3,4)+dot(p,3,5)+dot(p,3,6)
     &   +dot(p,4,5)+dot(p,4,6)+dot(p,5,6))
      elseif(approxim.ne.0)then
       call qqb_hzzexact(p,msqc)
       q2=2*(dot(p,3,4)+dot(p,3,5)+dot(p,3,6)
     &   +dot(p,4,5)+dot(p,4,6)+dot(p,5,6))
      endif

C    Only gg channel is non vanishing

      msqb=msqc(0,0)
      
            
      flux=fbGeV2/(2d0*xx(1)*xx(2)*W)
      
c--- initialize a PDF set here, if calculating errors
  777 continue    
      xmsq=0d0


      asopi=ason2pi*2



C     Compute Q2

      LF=dlog(q2/facscale**2)
      LR=dlog(q2/scale**2)


C Scaled momentum fractions

      cut=1d-7
   
C ndim is here 6 for H->2gamma and 12 for H->WW

      beta=cut+(1-cut)*r(ndim-1)
      alfa=cut+(1-cut)*r(ndim)


      xx10=xx(1)
      xx20=xx(2)

      z1=xx10**beta
      z2=xx20**alfa


c--- calculate PDF's  

      call fdist(ih1,xx10,facscale,fx10)
      call fdist(ih2,xx20,facscale,fx20)

      call fdist(ih1,xx10**(1-beta),facscale,fx1p)
      call fdist(ih2,xx20**(1-alfa),facscale,fx2p)

       if(noglue) then
        fx10(0)=0d0
        fx20(0)=0d0
        fx1p(0)=0d0
        fx2p(0)=0d0
       endif

       flgq=1
       if(gqonly)flgq=0


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC Start calculation

      tdelta=0d0
      tH1st=0d0
      tH1stH=0d0
      tH1stF=0d0
      tH2st=0d0

      diff10=0d0
      diff20=0d0
      diff1f=0d0
      diff2f=0d0

      diffc10=0d0
      diffc1f=0d0
      diffc20=0d0
      diffc2f=0d0

      diffg10=0d0
      diffg1f=0d0
      diffg20=0d0
      diffg2f=0d0

      spgq1=0d0
      spgq2=0d0

      tgamma2=0d0
      tgaga=0d0
      tcga=0d0

      tH2sp=0d0


CC    First diagonal part, then loop over j=1,nf

      if(gqonly)goto 71

C     Simplest term without convolutions
  

      tdelta=tdelta+fx10(0)*fx20(0)*msqb  

      if(order.eq.0)goto 72

C     H1st delta term

      tH1st=tH1st+2*C1ggdeltaex*fx10(0)*fx20(0)*msqb
      tH1stH=tH1stH+2*C1ggdelta*fx10(0)*fx20(0)*msqb

c     gammagg: non delta terms, first leg    


      diff=-dlog(xx10)*((fx1p(0)-fx10(0)*xx10**beta)*3/(1-z1)
     &    +fx1p(0)*Pggreg(z1))

      tH1stF=tH1stF+diff*fx20(0)*msqb
      tH1stF=tH1stF-3*D0int(xx10)*fx10(0)*fx20(0)*msqb

c     gammagg: non delta terms, second leg    

      diff=-dlog(xx20)*((fx2p(0)-fx20(0)*xx20**alfa)*3/(1-z2)
     &    +fx2p(0)*Pggreg(z2))

      tH1stF=tH1stF+diff*fx10(0)*msqb
      tH1stF=tH1stF-3*D0int(xx20)*fx10(0)*fx20(0)*msqb

c     gammagg: delta term, both legs


      tH1stF=tH1stF+2*beta0*fx10(0)*fx20(0)*msqb

c     H2st, gg channel: D0(z), first leg

      diff=-dlog(xx10)*((fx1p(0)-fx10(0)*xx10**beta)*H2ggD0/(1-z1)
     &    +fx1p(0)*H2ggreg(z1))

      tH2st=tH2st+0.5d0*diff*fx20(0)*msqb
      tH2st=tH2st-0.5d0*H2ggD0*D0int(xx10)*fx10(0)*fx20(0)*msqb

c     H2st, gg channel: D0(z), second leg
      
      diff=-dlog(xx20)*((fx2p(0)-fx20(0)*xx20**alfa)*H2ggD0/(1-z2)
     &    +fx2p(0)*H2ggreg(z2))

      tH2st=tH2st+0.5d0*diff*fx10(0)*msqb
      tH2st=tH2st-0.5d0*H2ggD0*D0int(xx20)*fx10(0)*fx20(0)*msqb


CCCCCC  Terms required for NNLO Scale dependence



CC    (gamma+gamma)*(gamma+gamma) term

C     First part: one gamma for each leg: gluon channel

      diff10=-dlog(xx10)
     &  *(fx1p(0)-fx10(0)*xx10**beta)/(1-z1)
     &  - D0int(xx10)*fx10(0)

      diff10=3*diff10+beta0*fx10(0)
     &     -dlog(xx10)*fx1p(0)*Pggreg(z1)

      diff20=-dlog(xx20)
     &  *(fx2p(0)-fx20(0)*xx20**alfa)/(1-z2)
     &  - D0int(xx20)*fx20(0)

      diff20=3*diff20+beta0*fx20(0)
     &     -dlog(xx20)*fx2p(0)*Pggreg(z2)


C    Second: gamma*gamma: gluon channel

c    First leg

      diff=-dlog(xx10)*((fx1p(0)-fx10(0)*xx10**beta)
     &    *(D0gggg/(1-z1)+D1gggg*dlog(1-z1)/(1-z1))
     &      +fx1p(0)*(Pggggreg(z1)+Pgqqg(z1)))
     &    +(Deltagggg-D0gggg*D0int(xx10)-D1gggg*D1int(xx10))*fx10(0)


      tgaga=tgaga+diff*flgq*fx20(0)*msqb


c    Second leg

      diff=-dlog(xx20)*((fx2p(0)-fx20(0)*xx20**alfa)
     &    *(D0gggg/(1-z2)+D1gggg*dlog(1-z2)/(1-z2))
     &      +fx2p(0)*(Pggggreg(z2)+Pgqqg(z2)))
     &    +(Deltagggg-D0gggg*D0int(xx20)-D1gggg*D1int(xx20))*fx20(0)

      tgaga=tgaga+diff*flgq*fx10(0)*msqb



CC    Start  (C+C)*(gamma+gamma) term: diagonal part


c    gamma first leg, C second

      diffg10=-dlog(xx10)
     &  *((fx1p(0)-fx10(0)*xx10**beta)*3d0/(1-z1)+Pggreg(z1)*fx1p(0))
     &    +fx10(0)*(beta0-3*D0int(xx10))

      diffc20=C1ggdelta*fx20(0)

c    gamma second leg, C first

      diffg20=-dlog(xx20)
     &  *((fx2p(0)-fx20(0)*xx20**alfa)*3d0/(1-z2)+Pggreg(z2)*fx2p(0))
     &    +fx20(0)*(beta0-3*D0int(xx20))


      diffc10=C1ggdelta*fx10(0)

c    C*gamma: first leg (ignore delta term in Cgg: taken into account in H1stf)

      tcga=tcga+CgqPqg(z1)*(-dlog(xx10))*flgq*fx1p(0)*fx20(0)*msqb 

c    C*gamma: second leg (ignore delta term in Cgg: taken into account in H1stf)

      tcga=tcga+CgqPqg(z2)*(-dlog(xx20))*flgq*fx2p(0)*fx10(0)*msqb 

c    End of (C+C)*(gamma+gamma)


CC    gamma2: diagonal part

c     First leg

      diff=-dlog(xx10)
     &  *(fx1p(0)-fx10(0)*xx10**beta)/(1-z1)
     &  - D0int(xx10)*fx10(0)  

      tgamma2=tgamma2+(1.5d0*Kappa*diff-dlog(xx10)*P2gg(z1)*fx1p(0))
     &                *flgq*fx20(0)*msqb


c     Second leg

      diff=-dlog(xx20)
     &  *(fx2p(0)-fx20(0)*xx20**alfa)/(1-z2)
     &  - D0int(xx20)*fx20(0)  

      tgamma2=tgamma2+(1.5d0*Kappa*diff-dlog(xx20)*P2gg(z2)*fx2p(0))
     &                *flgq*fx10(0)*msqb




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 71   continue

      do j=1,nf


C     H1st: Cgq, first leg


      tH1st=tH1st+(fx1p(j)+fx1p(-j))*Cgq(z1)
     & *(-dlog(xx10))*fx20(0)*msqb
      tH1stH=tH1stH+(fx1p(j)+fx1p(-j))*Cgq(z1)
     & *(-dlog(xx10))*fx20(0)*msqb


C     H1st: Cgq, second leg

      
      tH1st=tH1st+(fx2p(j)+fx2p(-j))*Cgq(z2)
     & *(-dlog(xx20))*fx10(0)*msqb
      tH1stH=tH1stH+(fx2p(j)+fx2p(-j))*Cgq(z2)
     & *(-dlog(xx20))*fx10(0)*msqb

    

C     H1st: muf dependence: Pgq, first leg


      tH1stF=tH1stF+(-dlog(xx10))
     & *(fx1p(j)+fx1p(-j))*Pgq(z1)*fx20(0)*msqb


C     H1st: muf dependence: Pgq, second leg


      tH1stF=tH1stF+(-dlog(xx20))
     & *(fx2p(j)+fx2p(-j))*Pgq(z2)*fx10(0)*msqb
  

C NEW: Add H2 qg contribution

c First leg


      tH2st=tH2st+(fx1p(j)+fx1p(-j))*H2stgq(z1)
     & *(-dlog(xx10))*fx20(0)*msqb

c Second leg


      tH2st=tH2st+(fx2p(j)+fx2p(-j))*H2stgq(z2)
     & *(-dlog(xx20))*fx10(0)*msqb


CCCCC Terms needed for NNLO scale dependence


CC    Now (gamma+gamma)*(gamma+gamma) term

C     First part: one gamma for each leg


      diff1f=diff1f-dlog(xx10)*Pgq(z1)*(fx1p(j)+fx1p(-j))

      diff2f=diff2f-dlog(xx20)*Pgq(z2)*(fx2p(j)+fx2p(-j))


C     Second part: gamma*gamma terms

c     Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
c              + Pijjk(z) + Deltaijjk delta(1-z)

      tgaga=tgaga-dlog(xx10)*(Pgqqq(z1)+Pgggq(z1))*(fx1p(j)+fx1p(-j))
     &            *fx20(0)*msqb

      tgaga=tgaga-dlog(xx20)*(Pgqqq(z2)+Pgggq(z2))*(fx2p(j)+fx2p(-j))
     &            *fx10(0)*msqb


C    End of (gamma+gamma)*(gamma+gamma) term


C    Start (C+C)*(gamma+gamma)


c    Gamma first leg

      diffg1f=diffg1f-dlog(xx10)*(fx1p(j)+fx1p(-j))*Pgq(z1)

c    C second leg

      diffc2f=diffc2f-dlog(xx20)*(fx2p(j)+fx2p(-j))*Cgq(z2)

c    Gamma second leg

      diffg2f=diffg2f-dlog(xx20)*(fx2p(j)+fx2p(-j))*Pgq(z2)

c    C first leg

      diffc1f=diffc1f-dlog(xx10)*(fx1p(j)+fx1p(-j))*Cgq(z1)

c    C*gamma: first leg (ignore delta term in Cgg: taken into account in H1stf)


      tcga=tcga+CgqPqq(z1)*(-dlog(xx10))*(fx1p(j)+fx1p(-j))*fx20(0)*msqb 


c    C*gamma: second leg (ignore delta term in Cgg: taken into account in H1stf)


      tcga=tcga+CgqPqq(z2)*(-dlog(xx20))*(fx2p(j)+fx2p(-j))*fx10(0)*msqb 


CC    gamma2: qg channel


c    First leg

      tgamma2=tgamma2
     &       -dlog(xx10)*P2gq(z1)*(fx1p(j)+fx1p(-j))*fx20(0)*msqb 

c    Second leg

      tgamma2=tgamma2
     &       -dlog(xx20)*P2gq(z2)*(fx2p(j)+fx2p(-j))*fx10(0)*msqb 



CC   Effect of spin correlations

c     First leg

      spgq1=spgq1-dlog(xx10)*(fx1p(j)+fx1p(-j))*Ggq(z1)

c    Second leg

      spgq2=spgq2-dlog(xx20)*(fx2p(j)+fx2p(-j))*Ggq(z2)

      enddo


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c    Check it !      

      tgaga=tgaga+2*msqb
     # *(flgq*diff10*diff20+flgq*diff1f*diff2f
     #   +diff10*diff2f+diff1f*diff20) 

c

c    gamma first leg, C second leg

      tcga=tcga+msqb*
     # (flgq*diffg10*diffc20+flgq*diffg1f*diffc2f
     #          +diffg10*diffc2f+diffg1f*diffc20)

c    gamma second leg, C first leg

      tcga=tcga+msqb*
     # (flgq*diffg20*diffc10+flgq*diffg2f*diffc1f
     #          +diffg20*diffc1f+diffg2f*diffc10)



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c    qq contribution to H2 (C1C1)


      tH2st=tH2st+msqb*diffc1f*diffc2f*flgq


c    spin correlations

c    gg channel

      
      tH2sp=msqb*dlog(xx10)*dlog(xx20)*
     . fx1p(0)*fx2p(0)*Ggg(z1)*Ggg(z2)*flgq

c    gq+qg channel

      tH2sp=tH2sp-msqb*dlog(xx10)*fx1p(0)*Ggg(z1)*spgq2
      tH2sp=tH2sp-msqb*dlog(xx20)*fx2p(0)*Ggg(z2)*spgq1

c    qq channel


      tH2sp=tH2sp+msqb*spgq1*spgq2*flgq


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 72   xmsq=tdelta


      if(order.eq.1)then
       xmsq=xmsq+asopi*(tH1st+LF*tH1stF-2*beta0*LR*tdelta)
      endif     


      if(order.eq.2 .and. approxim.eq.0) then
       xmsq=xmsq+asopi*(tH1st+LF*tH1stF-2*beta0*LR*tdelta)
     &          +asopi**2*(tdelta*H2ggdelta+tH2st+tH2sp)
      
CC     add scale dependence at NNLO
  
       xmsq=xmsq+asopi**2*(0.5d0*beta0*LF**2*tH1stF
     &               +tgamma2*LF
     &               -2*beta0*LR*(tH1st+LF*tH1stF-2*beta0*LR*tdelta)
     &               -beta0*LF*LR*tH1stF
     &               +LF*tcga-beta0*LR*tH1st+0.5d0*LF**2*tgaga       
     &               -2*(0.5d0*(beta0*LR)**2+beta1*LR)*tdelta)


c     Include missing delta term from C*gamma (no factor 2 here !)

      xmsq=xmsq+asopi**2*(LF*C1ggdelta*tH1stF)

C     Include missing term from contact term in 2 loop AP

      xmsq=xmsq+asopi**2*(2*Delta2gg*tdelta)*LF
      
      elseif(order.eq.2 .and. approxim .ne.0) then
       xmsq=xmsq+asopi*(tH1st+LF*tH1stF-2*beta0*LR*tdelta)
     &          +asopi**2*(tdelta*H2ggdelta+tH2st+tH2sp)
     &          *sigmat/sigmatb
      
CC     add scale dependence at NNLO
  
       xmsq=xmsq+asopi**2*(0.5d0*beta0*LF**2*tH1stF
     &               +tgamma2*LF
     &               -2*beta0*LR*(tH1stH+LF*tH1stF-2*beta0*LR*tdelta)
     &               -beta0*LF*LR*tH1stF
     &               +LF*tcga-beta0*LR*tH1stH+0.5d0*LF**2*tgaga       
     &               -2*(0.5d0*(beta0*LR)**2+beta1*LR)*tdelta)
     &               *sigmat/sigmatb


c     Include missing delta term from C*gamma (no factor 2 here !)

      xmsq=xmsq+asopi**2*(LF*C1ggdelta*tH1stF)*sigmat/sigmatb

C     Include missing term from contact term in 2 loop AP

      xmsq=xmsq+asopi**2*(2*Delta2gg*tdelta)*LF*sigmat/sigmatb


      endif



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      lowintHst=flux*pswt*xmsq/BrnRat


      call getptildejet(0,pjet)
      
      call dotem(nvec,pjet,s)

    

      val=lowintHst*wgt

c--- update the maximum weight so far, if necessary
c---  but not if we are already unweighting ...
c      if ((.not.unweight) .and. (dabs(val) .gt. wtmax)) then
c        wtmax=dabs(val)
c      endif

      if (bin) then
        val=val/dfloat(itmx)
CC      call plotter(pjet,val,0)
        call plotter(p,val,0)
      endif

      return

 999  continue
      ntotzero=ntotzero+1
      
      return
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




      function H2stgq(z)
      implicit none
      real *8 H2stgq,H2stgq1,H2stgq2,H2stgq3,H2stgqfull
      real *8 Pi,Z3,myli2,myli3,z
      integer nf

      external myli2,myli3

      Pi=3.14159265358979d0
      Z3=1.20205690316d0

      nf=5


c    Full result: Check it !

            H2stgqfull=
     #      -(12072-224*nf-396*Pi**2-12444*z+224*nf*z+432*Pi**2*z+
     #      2064*z**2-52*nf*z**2-432*Pi**2*z**2-1824*z**3+
     #      72*Pi**2*z**3+1584*dlog(1-z)-120*nf*dlog(1-z)-
     #      162*Pi**2*dlog(1-z)-1584*z*dlog(1-z)+120*nf*z*dlog(1-z)-
     #      324*Pi**2*z*dlog(1-z)+414*z**2*dlog(1-z) -
     #      24*nf*z**2*dlog(1-z)-162*Pi**2*z**2*dlog(1-z) +
     #      378*dlog(1-z)**2-36*nf*dlog(1-z)**2-378*z*dlog(1-z)**2 +
     #      36*nf*z*dlog(1-z)**2+99*z**2*dlog(1-z)**2 -
     #      18*nf*z**2*dlog(1-z)**2+60*dlog(1-z)**3 -
     #      60*z*dlog(1-z)**3+30*z**2*dlog(1-z)**3+1296*dlog(z) +
     #      6318*z*dlog(z)-288*z**2*dlog(z)+1584*z**3*dlog(z) +
     #      2376*dlog(1-z)*dlog(z)-2592*z*dlog(1-z)*dlog(z) +
     #      972*z**2*dlog(1-z)*dlog(z)-432*z**3*dlog(1-z)*dlog(z) +
     #      324*dlog(1-z)**2*dlog(z)+1620*z*dlog(1-z)**2*dlog(z) +
     #      486*z**2*dlog(1-z)**2*dlog(z)-900*z*dlog(z)**2 -
     #      189*z**2*dlog(z)**2-216*z**3*dlog(z)**2 -
     #      324*dlog(1-z)*dlog(z)**2+324*z*dlog(1-z)*dlog(z)**2 -
     #      162*z**2*dlog(1-z)*dlog(z)**2+84*z*dlog(z)**3 +
     #      66*z**2*dlog(z)**3+270*Pi**2*dlog(1+z)+
     #      432*Pi**2*z*dlog(1+z)+216*Pi**2*z**2*dlog(1+z)-
     #      324*z**2*dlog(z)*dlog(1+z)-1296*dlog(1-z)*dlog(z)*dlog(1+z)-
     #      2592*z*dlog(1-z)*dlog(z)*dlog(1+z) -
     #      1296*z**2*dlog(1-z)*dlog(z)*dlog(1+z) +
     #      324*dlog(z)**2*dlog(1+z)+324*z*dlog(z)**2*dlog(1+z) +
     #      162*z**2*dlog(z)**2*dlog(1+z)+648*dlog(z)*dlog(1+z)**2+
     #      1296*z*dlog(z)*dlog(1+z)**2+648*z**2*dlog(z)*dlog(1+z)**2 -
     #      216*dlog(1+z)**3-216*z*dlog(1+z)**3 -
     #      108*z**2*dlog(1+z)**3 -
     #      324*(z**2+2*(1+z)**2*dlog(1-z)+(2+2*z+z**2)*dlog(z) -
     #      2*dlog(1+z)-4*z*dlog(1+z)-2*z**2*dlog(1+z))*myli2(-z)-108*
     #       (-22+24*z-6*z**2+4*z**3-6*(1+z)**2*dlog(1-z) +
     #         3*(6-2*z+3*z**2)*dlog(z)+6*dlog(1+z)+12*z*dlog(1+z)+
     #         6*z**2*dlog(1+z))*myli2(z) +
     #      648*dlog(1-z)*myli2((1-z)/(1+z)) +
     #      1296*z*dlog(1-z)*myli2((1-z)/(1+z)) +
     #      648*z**2*dlog(1-z)*myli2((1-z)/(1+z)) -
     #      648*dlog(1+z)*myli2((1-z)/(1+z)) -
     #      1296*z*dlog(1+z)*myli2((1-z)/(1+z)) -
     #      648*z**2*dlog(1+z)*myli2((1-z)/(1+z)) -
     #      648*dlog(1-z)*myli2((-1+z)/(1+z)) -
     #      1296*z*dlog(1-z)*myli2((-1+z)/(1+z)) -
     #      648*z**2*dlog(1-z)*myli2((-1+z)/(1+z)) +
     #      648*dlog(1+z)*myli2((-1+z)/(1+z)) +
     #      1296*z*dlog(1+z)*myli2((-1+z)/(1+z)) +
     #      648*z**2*dlog(1+z)*myli2((-1+z)/(1+z)) +
     #      1944*myli3(-z)+1944*z*myli3(-z) +
     #      972*z**2*myli3(-z)+3240*myli3(z)-648*z*myli3(z) +
     #      1620*z**2*myli3(z)+1296*myli3(1/(1+z)) +
     #      1296*z*myli3(1/(1+z))+648*z**2*myli3(1/(1+z)) -
     #      5184*Z3+3240*z*Z3-2592*z**2*Z3)/(324d0*z)


C     Subtract Ggg*Ggq term
    
      H2stgq=H2stgqfull-(-4)*(2*(1-z)+(1+z)*dlog(z))/4d0

      return
      end


C     Regular part of H2st in gg channel

      function H2ggREG(z)
      implicit none
      real *8 H2ggREG
      real *8 Pi,Z3,myli2,myli3,z
      integer nf

      external myli2,myli3

      Pi=3.14159265358979d0
      Z3=1.20205690316d0

      nf=5

c      H2ggREG=(16164-226*nf+396*Pi**2-17790*z-52*nf*z+1512*Pi**2*z- 
c     -  5064*z**2 + 390*nf*z**2 + 684*Pi**2*z**2 + 5784*z**3 + 
c     -  314*nf*z**3 - 972*Pi**2*z**3 - 11100*z**4 - 164*nf*z**4 - 
c     -  1080*Pi**2*z**4 + 12006*z**5 - 262*nf*z**5 - 540*Pi**2*z**5- 
c     -  9648*dlog(1-z) + 864*Pi**2*dlog(1-z) + 7992*z*dlog(1-z) - 
c     -  1728*Pi**2*z*dlog(1-z) + 7434*z**2*dlog(1-z) + 
c     -  18*nf*z**2*dlog(1-z) - 4176*z**3*dlog(1-z) + 
c     -  864*Pi**2*z**3*dlog(1-z) + 2214*z**4*dlog(1-z) - 
c     -  18*nf*z**4*dlog(1-z) - 864*Pi**2*z**4*dlog(1-z) - 
c     -  3816*z**5*dlog(1-z) + 864*Pi**2*z**5*dlog(1-z) - 
c     -  432*dlog(1-z)**3-2592*z*dlog(1-z)**3-864*z**2*dlog(1-z)**3+ 
c     -  2160*z**3*dlog(1-z)**3 + 1296*z**4*dlog(1-z)**3 + 
c     -  432*z**5*dlog(1-z)**3 + 3060*dlog(z) - 216*Pi**2*dlog(z)+ 
c     -  693*z*dlog(z) - 222*nf*z*dlog(z) + 432*Pi**2*z*dlog(z) - 
c     -  4203*z**2*dlog(z)-204*nf*z**2*dlog(z)+2727*z**3*dlog(z) + 
c     -  222*nf*z**3*dlog(z)-216*Pi**2*z**3*dlog(z)+1143*z**4*dlog(z)+ 
c     -  204*nf*z**4*dlog(z)+216*Pi**2*z**4*dlog(z)-3420*z**5*dlog(z)- 
c     -  216*Pi**2*z**5*dlog(z) + 2376*dlog(1-z)*dlog(z)- 
c     -  2592*z*dlog(1-z)*dlog(z) + 216*z**2*dlog(1-z)*dlog(z)+ 
c     -  216*z**3*dlog(1-z)*dlog(z) - 2592*z**4*dlog(1-z)*dlog(z)+ 
c     -  2376*z**5*dlog(1-z)*dlog(z) - 1620*dlog(1-z)**2*dlog(z) - 
c     -  7452*z*dlog(1-z)**2*dlog(z) - 2916*z**2*dlog(1-z)**2*dlog(z)+ 
c     -  6156*z**3*dlog(1-z)**2*dlog(z)+4212*z**4*dlog(1-z)**2*dlog(z)+ 
c     -  972*z**5*dlog(1-z)**2*dlog(z)
c     -  -1188*dlog(z)**2-2295*z*dlog(z)**2- 
c     -  54*nf*z*dlog(z)**2 - 783*z**2*dlog(z)**2-30*nf*z**2*dlog(z)**2+ 
c     -  891*z**3*dlog(z)**2+54*nf*z**3*dlog(z)**2+1971*z**4*dlog(z)**2+ 
c     -  30*nf*z**4*dlog(z)**2+1404*z**5*dlog(z)**2 + 
c     -  1620*dlog(1-z)*dlog(z)**2 + 6804*z*dlog(1-z)*dlog(z)**2 + 
c     -  2268*z**2*dlog(1-z)*dlog(z)**2-6156*z**3*dlog(1-z)*dlog(z)**2- 
c     -  4212*z**4*dlog(1-z)*dlog(z)**2 - 972*z**5*dlog(1-z)*dlog(z)**2- 
c     -  432*dlog(z)**3 - 756*z*dlog(z)**3 - 8*nf*z*dlog(z)**3 - 
c     -  216*z**2*dlog(z)**3 - 8*nf*z**2*dlog(z)**3+756*z**3*dlog(z)**3+ 
c     -  8*nf*z**3*dlog(z)**3+648*z**4*dlog(z)**3+8*nf*z**4*dlog(z)**3+ 
c     -  108*z**5*dlog(z)**3 - 108*Pi**2*dlog(1+z) - 
c     -  108*Pi**2*z*dlog(1+z) - 108*Pi**2*z**2*dlog(1+z) + 
c     -  108*Pi**2*z**3*dlog(1+z) + 108*Pi**2*z**4*dlog(1+z) + 
c     -  108*Pi**2*z**5*dlog(1+z) + 324*dlog(z)**2*dlog(1+z) + 
c     -  324*z*dlog(z)**2*dlog(1+z) + 324*z**2*dlog(z)**2*dlog(1+z) - 
c     -  324*z**3*dlog(z)**2*dlog(1+z) - 324*z**4*dlog(z)**2*dlog(1+z) - 
c     -  324*z**5*dlog(z)**2*dlog(1+z) - 648*dlog(z)*dlog(1+z)**2 - 
c     -  648*z*dlog(z)*dlog(1+z)**2 - 648*z**2*dlog(z)*dlog(1+z)**2 + 
c     -  648*z**3*dlog(z)*dlog(1+z)**2 + 648*z**4*dlog(z)*dlog(1+z)**2 + 
c     -  648*z**5*dlog(z)*dlog(1+z)**2 + 216*dlog(1+z)**3 + 
c     -  216*z*dlog(1+z)**3 + 216*z**2*dlog(1+z)**3 - 
c     -  216*z**3*dlog(1+z)**3 - 216*z**4*dlog(1+z)**3 - 
c     -  216*z**5*dlog(1+z)**3 + 
c     -  648*(-1 + z)*(1 + z + z**2)**2*dlog(z)*myli2(-z) + 
c     -  216*(-11 - 42*z - 19*z**2 + 27*z**3 + 30*z**4 + 15*z**5 + 
c     -  12*(-1 - 6*z - 2*z**2 + 5*z**3 + 3*z**4 + z**5)*dlog(1-z) - 
c     -  3*(-1 - 9*z - z**2 + 9*z**3 + 5*z**4 + z**5)*dlog(z))*myli2(z)
c     -  + 648*myli3(-z) + 648*z*myli3(-z) + 
c     -  648*z**2*myli3(-z) - 648*z**3*myli3(-z) - 
c     -  648*z**4*myli3(-z) - 648*z**5*myli3(-z) + 
c     -  1944*myli3(z) + 12312*z*myli3(z) + 5832*z**2*myli3(z) - 
c     -  8424*z**3*myli3(z) - 4536*z**4*myli3(z) - 
c     -  3240*z**5*myli3(z) + 2592*myli3(z/(-1 + z)) + 
c     -  15552*z*myli3(z/(-1 + z)) + 5184*z**2*myli3(z/(-1 + z)) - 
c     -  12960*z**3*myli3(z/(-1 + z)) - 
c     -  7776*z**4*myli3(z/(-1 + z)) - 2592*z**5*myli3(z/(-1 + z))- 
c     -  1296*myli3(z/(1 + z)) - 1296*z*myli3(z/(1 + z)) - 
c     -  1296*z**2*myli3(z/(1 + z)) + 1296*z**3*myli3(z/(1 + z)) + 
c     -  1296*z**4*myli3(z/(1 + z)) + 1296*z**5*myli3(z/(1 + z)) - 
c     -  7776*Z3 + 4212*z*Z3 - 4212*z**2*Z3 - 
c     -  648*z**3*Z3 + 10368*z**4*Z3 - 5832*z**5*Z3)/
c     -  (72d0*z*(-1 + z**2))



        H2ggREG=-(-10776+226*nf+396*Pi**2+12408*z+52*nf*z-432*Pi**2*z+ 
     -   1656*z**2 - 390*nf*z**2 + 36*Pi**2*z**2 - 2388*z**3 - 
     -   314*nf*z**3 + 36*Pi**2*z**3 + 9120*z**4 + 164*nf*z**4 - 
     -   432*Pi**2*z**4 - 10020*z**5 + 262*nf*z**5 + 396*Pi**2*z**5 + 
     -   54*z**2*dlog(1-z) - 18*nf*z**2*dlog(1-z) - 54*z**4*dlog(1-z)+ 
     -   18*nf*z**4*dlog(1-z)-648*dlog(z)-6957*z*dlog(z)+ 
     -   222*nf*z*dlog(z)-693*z**2*dlog(z)+204*nf*z**2*dlog(z)+ 
     -   2133*z**3*dlog(z)-222*nf*z**3*dlog(z)+1341*z**4*dlog(z) - 
     -   204*nf*z**4*dlog(z)+4824*z**5*dlog(z)-2376*dlog(1-z)*dlog(z)+ 
     -   2592*z*dlog(1-z)*dlog(z)-216*z**2*dlog(1-z)*dlog(z)- 
     -   216*z**3*dlog(1-z)*dlog(z)+2592*z**4*dlog(1-z)*dlog(z)- 
     -   2376*z**5*dlog(1-z)*dlog(z) + 324*dlog(1-z)**2*dlog(z)- 
     -   324*z*dlog(1-z)**2*dlog(z)+324*z**2*dlog(1-z)**2*dlog(z) + 
     -   324*z**3*dlog(1-z)**2*dlog(z)-324*z**4*dlog(1-z)**2*dlog(z)+ 
     -   324*z**5*dlog(1-z)**2*dlog(z) + 675*z*dlog(z)**2+ 
     -   54*nf*z*dlog(z)**2-297*z**2*dlog(z)**2+30*nf*z**2*dlog(z)**2+ 
     -   513*z**3*dlog(z)**2-54*nf*z**3*dlog(z)**2+297*z**4*dlog(z)**2- 
     -   30*nf*z**4*dlog(z)**2 - 1188*z**5*dlog(z)**2 + 
     -   324*dlog(1-z)*dlog(z)**2 - 324*z*dlog(1-z)*dlog(z)**2 + 
     -   324*z**2*dlog(1-z)*dlog(z)**2 + 324*z**3*dlog(1-z)*dlog(z)**2- 
     -   324*z**4*dlog(1-z)*dlog(z)**2 + 324*z**5*dlog(1-z)*dlog(z)**2- 
     -   108*z*dlog(z)**3 + 8*nf*z*dlog(z)**3-216*z**2*dlog(z)**3 + 
     -   8*nf*z**2*dlog(z)**3+108*z**3*dlog(z)**3-8*nf*z**3*dlog(z)**3+ 
     -   216*z**4*dlog(z)**3-8*nf*z**4*dlog(z)**3-108*z**5*dlog(z)**3+ 
     -   108*Pi**2*dlog(1+z) + 108*Pi**2*z*dlog(1+z) + 
     -   108*Pi**2*z**2*dlog(1+z) - 108*Pi**2*z**3*dlog(1+z) - 
     -   108*Pi**2*z**4*dlog(1+z) - 108*Pi**2*z**5*dlog(1+z) - 
     -   324*dlog(z)**2*dlog(1+z) - 324*z*dlog(z)**2*dlog(1+z) - 
     -   324*z**2*dlog(z)**2*dlog(1+z) + 324*z**3*dlog(z)**2*dlog(1+z)+ 
     -   324*z**4*dlog(z)**2*dlog(1+z) + 324*z**5*dlog(z)**2*dlog(1+z)+ 
     -   648*dlog(z)*dlog(1+z)**2 + 648*z*dlog(z)*dlog(1+z)**2 + 
     -   648*z**2*dlog(z)*dlog(1+z)**2 - 648*z**3*dlog(z)*dlog(1+z)**2- 
     -   648*z**4*dlog(z)*dlog(1+z)**2 - 648*z**5*dlog(z)*dlog(1+z)**2- 
     -   216*dlog(1+z)**3 - 216*z*dlog(1+z)**3 - 216*z**2*dlog(1+z)**3+ 
     -   216*z**3*dlog(1+z)**3 + 216*z**4*dlog(1+z)**3 + 
     -   216*z**5*dlog(1+z)**3 - 
     -   648*(-1 + z)*(1 + z + z**2)**2*dlog(z)*myli2(-z) + 
     -   216*(-((-1 + z)**2*(11 + 10*z + 10*z**2 + 11*z**3)) + 
     -   3*(3 - z + 3*z**2 + z**3 - 3*z**4 + z**5)*dlog(z))*myli2(z)
     -    - 648*myli3(-z) - 648*z*myli3(-z) - 
     -   648*z**2*myli3(-z) + 648*z**3*myli3(-z) + 
     -   648*z**4*myli3(-z) + 648*z**5*myli3(-z) - 
     -   3240*myli3(z) + 648*z*myli3(z) - 3240*z**2*myli3(z) - 
     -   648*z**3*myli3(z) + 3240*z**4*myli3(z) - 
     -   648*z**5*myli3(z) + 1296*myli3(z/(1 + z)) + 
     -   1296*z*myli3(z/(1 + z)) + 1296*z**2*myli3(z/(1 + z)) - 
     -   1296*z**3*myli3(z/(1 + z)) - 1296*z**4*myli3(z/(1 + z)) - 
     -   1296*z**5*myli3(z/(1 + z)) + 3888*Z3 - 6804*z*Z3 + 
     -   1620*z**2*Z3 + 4536*z**3*Z3 - 3888*z**4*Z3 + 
     -   4536*z**5*Z3)/(72d0*z*(-1 + z**2))

C     Subtract Ggg*Ggg

      H2ggREG=H2ggREG-(-9)*(2*(1-z)+(1+z)*dlog(z))/4d0

      return
      end

CC Spin correlations

      function Ggq(z)
      real *8 z,Ggq,CF

      CF=4d0/3
      
      Ggq=CF*(1-z)/z
    

      return
      end

      function Ggg(z)
      real *8 z,Ggg,CA

      CA=3d0
      
      Ggg=CA*(1-z)/z
    

      return
      end
