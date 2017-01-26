      double precision function virtint(r,wgt)
      implicit none
      include 'constants.f'
      include 'noglue.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'npart.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'agq.f'
      include 'PR_new.f'
      include 'PR_cs_new.f'
      include 'PR_twojet.f'
      include 'msq_cs.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'clustering.f'
      include 'efficiency.f'
      include 'lc.f'
      include 'process.f'
      include 'maxwt.f'
      include 'limits.f'
      include 'b0.f'
      double precision mqq(0:2,fn:nf,fn:nf)
      double precision msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      double precision msqx_cs(0:2,-nf:nf,-nf:nf)
      double precision AP(-1:1,-1:1,3)
CC
CC    Variables to be passed to the counterterm
CC      
      double precision qt2,qq2,shat,dot
      common/count/qt2,qq2,shat
CC

      integer ih1,ih2,j,k,cs,nvec,is,ia,ib,ic
      double precision p(mxpart,4),pjet(mxpart,4),r(mxdim),W,sqrts,xmsq,
     . val,fx1(-nf:nf),fx2(-nf:nf),fx1z(-nf:nf),fx2z(-nf:nf)
      double precision pswt,xjac,rscalestart,fscalestart,
     . wgt,msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),msqvdk(-nf:nf,-nf:nf),
     . msq_qq,msq_aa,msq_aq,msq_qa,msq_qg,msq_gq,epcorr,
     . msqr(-nf:nf,-nf:nf)
      double precision xx(2),z,x1onz,x2onz,flux,omz,
     . BrnRat,xmsq_old,tmp,resfac
      integer nshot,rvcolourchoice,sgnj,sgnk
      logical bin,first,includedipole
      character*4 mypart
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/x1x2/xx
      common/BrnRat/BrnRat
      common/rvcolourchoice/rvcolourchoice
      common/mypart/mypart
      integer higgsdec,ndec
      common/higgsdec/higgsdec,ndec
      common/rescale/resfac

      data p/48*0d0/
      data nshot/1/
      data first/.true./
      save first,rscalestart,fscalestart
      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale
      endif

      ntotshot=ntotshot+1
      virtint=0d0

      W=sqrts**2


      

      if(higgsdec.eq.1) then   
       npart=3   
       call gen3(r,p,pswt,*999)
       qq2=2*dot(p,3,4)
       qt2=p(5,1)**2+p(5,2)**2  
      else
       npart=5
       call gen5(r,p,pswt,*999)
       qq2=2*(dot(p,3,4)+dot(p,3,5)+dot(p,3,6)
     &     +dot(p,4,5)+dot(p,4,6)+dot(p,5,6))
       qt2=p(7,1)**2+p(7,2)**2
      endif

      shat=2*dot(p,1,2)

      nvec=npart+2
      call dotem(nvec,p,s)

c---impose mass cuts on final state
      call masscuts(s,*999)
c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)
         
c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      if (includedipole(0,p) .eqv. .false.) then
        goto 999
      endif
      
     
      z=r(ndim)**2
      if (nshot .eq. 1) z=0.95d0
      xjac=two*dsqrt(z)

      omz=1d0-z

      flux=fbGeV2/(2d0*xx(1)*xx(2)*W)

c--- to test poles, we need colourchoice=0, but save real value
      if (nshot .eq. 1) then
        rvcolourchoice=colourchoice
        colourchoice=0
      endif
      
   12 continue
c--- point to restart from when checking epsilon poles

c--- correction to epinv from AP subtraction when mu_FAC != mu_REN,
c--- corresponding to subtracting -1/epinv*Pab*log(musq_REN/musq_FAC)
      epcorr=epinv+2d0*dlog(scale/facscale)
          

      AP(q,q,1)=+ason2pi*Cf*1.5d0*epcorr
      AP(q,q,2)=+ason2pi*Cf*(-1d0-z)*epcorr
      AP(q,q,3)=+ason2pi*Cf*2d0/omz*epcorr
      AP(a,a,1)=+ason2pi*Cf*1.5d0*epcorr
      AP(a,a,2)=+ason2pi*Cf*(-1d0-z)*epcorr
      AP(a,a,3)=+ason2pi*Cf*2d0/omz*epcorr

      AP(q,g,1)=0d0
      AP(q,g,2)=ason2pi*Tr*(z**2+omz**2)*epcorr
      AP(q,g,3)=0d0
      AP(a,g,1)=0d0
      AP(a,g,2)=ason2pi*Tr*(z**2+omz**2)*epcorr
      AP(a,g,3)=0d0

      AP(g,q,1)=0d0
      AP(g,q,2)=ason2pi*Cf*(1d0+omz**2)/z*epcorr
      AP(g,q,3)=0d0
      AP(g,a,1)=0d0
      AP(g,a,2)=ason2pi*Cf*(1d0+omz**2)/z*epcorr
      AP(g,a,3)=0d0

      AP(g,g,1)=+ason2pi*b0*epcorr
      AP(g,g,2)=+ason2pi*xn*2d0*(1d0/z+z*omz-2d0)*epcorr
      AP(g,g,3)=+ason2pi*xn*2d0/omz*epcorr

      
      do ia=-1,+1
      do ib=-1,+1
      do ic=-1,+1
      do is=1,3
        Q1(ia,ib,ic,is)=0d0
        Q2(ia,ib,ic,is)=0d0
      do cs=0,2
        R1(ia,ib,ic,cs,is)=0d0
        R2(ia,ib,ic,cs,is)=0d0
      do j=1,8
        S1(ia,ib,ic,j,cs,is)=0d0
        S2(ia,ib,ic,j,cs,is)=0d0
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
     
c--- Calculate the required matrix elements      

      if(higgsdec.eq.1) then
        call gg_hg(p,msq)
        call gg_hg_rescaled(p,msqr)
        call gg_hg_v(p,msqv)
        call gg_hg_z(p,z)
      elseif(higgsdec.eq.2) then
        call qqb_hww_g(p,msq)
        call qqb_hww_g_rescaled(p,msqr)
        call gg_hwwg_v(p,msqv)
        call gg_hwwg_z(p,z)
      else
        call qqb_hzz_g(p,msq)
        call qqb_hzz_g_rescaled(p,msqr)
        call gg_hzzg_v(p,msqv)
        call gg_hzzg_z(p,z)
      endif      

      msqr=msqr*resfac

            
  777 continue    
      xmsq=0d0

      call fdist(ih1,xx(1),facscale,fx1)
      call fdist(ih2,xx(2),facscale,fx2)

      do j=-nf,nf
      fx1z(j)=0d0
      fx2z(j)=0d0
      enddo
            
      if (z .gt. xx(1)) then
         x1onz=xx(1)/z
         call fdist(ih1,x1onz,facscale,fx1z)
      endif
      if (z .gt. xx(2)) then
         x2onz=xx(2)/z
         call fdist(ih2,x2onz,facscale,fx2z)
      endif         
      
      do j=-nf,nf
      do k=-nf,nf

      
      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) goto 20
      endif

      if (gqonly) then
      if (((j.eq.0).and.(k.eq.0)) .or. ((j.ne.0).and.(k.ne.0))) goto 20
      endif      
      
      if (noglue) then 
      if ((j.eq.0) .or. (k.eq.0)) goto 20
      endif


      tmp=xmsq

c--- The variables R1 and R2 provide the Regular and Plus pieces associated
c--- with radiation from leg 1 (R1(a,b,c,cs,is)) and leg 2 (R2(a,b,c,cs,is))
c--- In each case the parton labelling is using the normal QM notation of 
c--- putting everything backward
c---       emitted line after emission =    a
c---       emitter before emission     =    b
c---       spectator                   =    c
c--- There is no label for he or she who is emitted.
c--- Note that in general each piece will be composed of many different
c--- dipole contributions


c--- SUM BY TOTAL MATRIX ELEMENTS: everything else
C--QQ
      if     ((j .gt. 0) .and. (k.gt.0)) then
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*one+msqr(j,k)*(AP(q,q,1)-AP(q,q,3)+Q1(q,q,q,1)
     &             -Q1(q,q,q,3)
     &             +AP(q,q,1)-AP(q,q,3)+Q2(q,q,q,1)-Q2(q,q,q,3)))
     &             *fx1(j)*fx2(k)
     & +(msqr(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,q,2)+Q1(q,q,q,3))
     & + msqr(g,k)*(AP(g,q,2)+Q1(g,q,q,2)))*fx1z(j)/z*fx2(k)
     & +(msqr(j,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,q,2)+Q2(q,q,q,3))
     & + msqr(j,g)*(AP(g,q,2)+Q2(g,q,q,2)))*fx1(j)*fx2z(k)/z
C--QbarQbar
      elseif ((j .lt. 0) .and. (k.lt.0)) then
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*one+msqr(j,k)*(AP(a,a,1)-AP(a,a,3)+Q1(a,a,a,1)
     &              -Q1(a,a,a,3)
     &              +AP(a,a,1)-AP(a,a,3)+Q2(a,a,a,1)-Q2(a,a,a,3)))
     &              *fx1(j)*fx2(k)
     & +(msqr(j,k)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,a,2)+Q1(a,a,a,3))
     & + msqr(g,k)*(AP(g,a,2)+Q1(g,a,a,2)))*fx1z(j)/z*fx2(k)
     & +(msqr(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,a,2)+Q2(a,a,a,3))
     & + msqr(j,g)*(AP(g,a,2)+Q2(g,a,a,2)))*fx1(j)*fx2z(k)/z
C--QQbar
      elseif ((j .gt. 0) .and. (k.lt.0)) then
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*one+msqr(j,k)*(AP(q,q,1)-AP(q,q,3)+Q1(q,q,a,1)
     &              -Q1(q,q,a,3)
     &              +AP(a,a,1)-AP(a,a,3)+Q2(a,a,q,1)-Q2(a,a,q,3)))
     &              *fx1(j)*fx2(k)
     & +(msqr(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,a,3)+Q1(q,q,a,2))
     & + msqr(g,k)*(AP(g,q,2)+Q1(g,q,a,2)))*fx1z(j)/z*fx2(k)
     & +(msqr(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,q,3)+Q2(a,a,q,2))
     & + msqr(j,g)*(AP(g,a,2)+Q2(g,a,q,2)))*fx1(j)*fx2z(k)/z

      elseif ((j .lt. 0) .and. (k.gt.0)) then
C--QbarQ
      xmsq=xmsq+(msqv(j,k)
     & +msq(j,k)*one+msqr(j,k)*(AP(a,a,1)-AP(a,a,3)+Q1(a,a,q,1)
     &             -Q1(a,a,q,3)
     &             +AP(q,q,1)-AP(q,q,3)+Q2(q,q,a,1)-Q2(q,q,a,3)))
     &             *fx1(j)*fx2(k)
     & +(msqr(j,k)*(AP(a,a,3)+AP(a,a,2)+Q1(a,a,q,3)+Q1(a,a,q,2))
     & + msqr(g,k)*(AP(g,a,2)+Q1(g,a,q,2)))*fx1z(j)/z*fx2(k)
     & +(msqr(j,k)*(AP(q,q,3)+AP(q,q,2)+Q2(q,q,a,3)+Q2(q,q,a,2))
     & + msqr(j,g)*(AP(g,q,2)+Q2(g,q,a,2)))*fx1(j)*fx2z(k)/z

      elseif ((j .eq. g) .and. (k.eq.g)) then
C--gg
    
       msq_qg=msqr(+5,g)+msqr(+4,g)+msqr(+3,g)+msqr(+2,g)+msqr(+1,g)
     &       +msqr(-5,g)+msqr(-4,g)+msqr(-3,g)+msqr(-2,g)+msqr(-1,g)
       msq_gq=msqr(g,+5)+msqr(g,+4)+msqr(g,+3)+msqr(g,+2)+msqr(g,+1)
     &       +msqr(g,-5)+msqr(g,-4)+msqr(g,-3)+msqr(g,-2)+msqr(g,-1)
       xmsq=xmsq+(msqv(g,g)
     &  +msq(g,g)*one+msqr(g,g)*(AP(g,g,1)-AP(g,g,3)+Q1(g,g,g,1)
     &                -Q1(g,g,g,3)
     &                +AP(g,g,1)-AP(g,g,3)+Q2(g,g,g,1)-Q2(g,g,g,3)))
     &                *fx1(g)*fx2(g)
     &  +(msqr(g,g)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,g,2)+Q1(g,g,g,3))
     &  +   msq_qg*(AP(q,g,2)+Q1(q,g,g,2)))*fx1z(g)/z*fx2(g)
     &  +(msqr(g,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,g,2)+Q2(g,g,g,3))
     &  +   msq_gq*(AP(q,g,2)+Q2(q,g,g,2)))*fx1(g)*fx2z(g)/z

      elseif (j .eq. g) then
C--gQ
       if    (k .gt. 0) then
       msq_aq=msqr(-1,k)+msqr(-2,k)+msqr(-3,k)+msqr(-4,k)+msqr(-5,k)
       msq_qq=msqr(+1,k)+msqr(+2,k)+msqr(+3,k)+msqr(+4,k)+msqr(+5,k)
       xmsq=xmsq+(msqv(g,k)
     & +msq(g,k)*one+msqr(g,k)*(AP(g,g,1)-AP(g,g,3)+Q1(g,g,q,1)
     &               -Q1(g,g,q,3)
     &               +AP(q,q,1)-AP(q,q,3)+Q2(q,q,g,1)-Q2(q,q,g,3)))
     &               *fx1(g)*fx2(k)
     & +(msqr(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,q,2)+Q1(g,g,q,3))
     & +   msq_aq*(AP(a,g,2)+Q1(a,g,q,2))
     & +   msq_qq*(AP(q,g,2)+Q1(q,g,q,2)))*fx1z(g)/z*fx2(k)
     & +(msqr(g,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,g,2)+Q2(q,q,g,3))
     & + msqr(g,g)*(AP(g,q,2)+Q2(g,q,g,2)))*fx1(g)*fx2z(k)/z
C--gQbar

       elseif (k.lt.0) then
       msq_qa=msqr(+1,k)+msqr(+2,k)+msqr(+3,k)+msqr(+4,k)+msqr(+5,k)
       msq_aa=msqr(-1,k)+msqr(-2,k)+msqr(-3,k)+msqr(-4,k)+msqr(-5,k)
       xmsq=xmsq+(msqv(g,k)
     & +msq(g,k)*one+msqr(g,k)*(AP(g,g,1)-AP(g,g,3)+Q1(g,g,a,1)
     &               -Q1(g,g,a,3)
     &               +AP(a,a,1)-AP(a,a,3)+Q2(a,a,g,1)-Q2(a,a,g,3)))
     &               *fx1(g)*fx2(k)
     & +(msqr(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,a,2)+Q1(g,g,a,3))
     & +   msq_qa*(AP(q,g,2)+Q1(q,g,a,2))
     & +   msq_aa*(AP(a,g,2)+Q1(a,g,a,2)))*fx1z(g)/z*fx2(k)
     & +(msqr(g,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,g,2)+Q2(a,a,g,3))
     & + msqr(g,g)*(AP(g,a,2)+Q2(g,a,g,2)))*fx1(g)*fx2z(k)/z
       endif
C--Qg
      elseif (k .eq. g) then
       if     (j.gt.0) then
       msq_qa=msqr(j,-1)+msqr(j,-2)+msqr(j,-3)+msqr(j,-4)+msqr(j,-5)
       msq_qq=msqr(j,+1)+msqr(j,+2)+msqr(j,+3)+msqr(j,+4)+msqr(j,+5)
       xmsq=xmsq+(msqv(j,g)
     & +msq(j,g)*one+msqr(j,g)*(AP(q,q,1)-AP(q,q,3)+Q1(q,q,g,1)
     &               -Q1(q,q,g,3)
     &               +AP(g,g,1)-AP(g,g,3)+Q2(g,g,q,1)-Q2(g,g,q,3)))
     &               *fx1(j)*fx2(g)
     & +(msqr(j,g)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,g,2)+Q1(q,q,g,3))
     & + msqr(g,g)*(AP(g,q,2)+Q1(g,q,g,2)))*fx1z(j)/z*fx2(g)
     & +(msqr(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,q,2)+Q2(g,g,q,3))
     & +   msq_qa*(AP(a,g,2)+Q2(a,g,q,2))
     & +   msq_qq*(AP(q,g,2)+Q2(q,g,q,2)))*fx1(j)*fx2z(g)/z
C--Qbarg
       elseif (j.lt.0) then
       msq_aq=msqr(j,+1)+msqr(j,+2)+msqr(j,+3)+msqr(j,+4)+msqr(j,+5)
       msq_aa=msqr(j,-1)+msqr(j,-2)+msqr(j,-3)+msqr(j,-4)+msqr(j,-5)
       xmsq=xmsq+(msqv(j,g)
     & +msq(j,g)*one+msqr(j,g)*(AP(a,a,1)-AP(a,a,3)+Q1(a,a,g,1)
     &               -Q1(a,a,g,3)
     &               +AP(g,g,1)-AP(g,g,3)+Q2(g,g,a,1)-Q2(g,g,a,3)))
     &               *fx1(j)*fx2(g)
     & +(msqr(j,g)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,g,2)+Q1(a,a,g,3))
     & + msqr(g,g)*(AP(g,a,2)+Q1(g,a,g,2)))*fx1z(j)/z*fx2(g)
     & +(msqr(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,a,3)+Q2(g,g,a,2))
     & + msq_aq*(AP(q,g,2)+Q2(q,g,a,2))
     & + msq_aa*(AP(a,g,2)+Q2(a,g,a,2)))*fx1(j)*fx2z(g)/z
       endif
      
      endif

      if     (j .gt. 0) then
        sgnj=+1
      elseif (j .lt. 0) then
        sgnj=-1
      else
        sgnj=0
      endif
      if     (k .gt. 0) then
        sgnk=+1
      elseif (k .lt. 0) then
        sgnk=-1
      else
        sgnk=0
      endif

      
 20   continue

      enddo
      enddo
      

      virtint=flux*xjac*pswt*xmsq/BrnRat

c--- code to check that epsilon poles cancel      
      if (nshot .eq. 1) then
        if (xmsq .eq. 0d0) goto 999
        xmsq_old=xmsq
        nshot=nshot+1
        epinv=0d0
        epinv2=0d0
c        epinv=1d0
c        epinv2=1d0
        goto 12
      elseif (nshot .eq. 2) then
        nshot=nshot+1

CC        if (abs(xmsq_old/xmsq-1d0) .gt. 1d-6) then
          if (abs(xmsq_old/xmsq-1d0) .gt. 1d-4) then
          write(6,*) 'epsilon fails to cancel'
          write(6,*) 'xmsq (epinv=large) = ',xmsq_old
          write(6,*) 'xmsq (epinv=zero ) = ',xmsq
c          stop
        else
c          write(6,*) 'Poles cancelled!'
          colourchoice=rvcolourchoice
        endif
      endif

      call getptildejet(0,pjet)
      
      call dotem(nvec,pjet,s)


      val=virtint*wgt 
c--- update the maximum weight so far, if necessary
      if (val .gt. wtmax) then
        wtmax=val
      endif

      if (bin) then
        val=val/dfloat(itmx) 
        call plotter(pjet,val,0)
      endif

      return

 999  continue
      ntotzero=ntotzero+1
      
      return
      end


