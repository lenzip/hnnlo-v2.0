      subroutine gg_hg_v(p,msq)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scheme.f'
C     (Taken from Ravindran, Smith, van Neerven hep-ph/0201114)
C     Modified by overall factors
      integer iglue,j,k,approxim
      double precision p(4,mxpart),msq(fn:nf,fn:nf)
      double precision ss,tt,uu,s34,
     . virtgg,virtqa,virtaq,virtqg,virtgq,hdecay,Asq,fac
      double precision mtop,mbot,sigmatb,sigmat
      common/tbmass/mtop,mbot
      common/flag/approxim
      common/Born/sigmatb,sigmat
      parameter(iglue=5)

      scheme='tH-V'

      call dotem(iglue,p,s)
      ss=s(1,2)
      tt=s(1,iglue)
      uu=s(2,iglue)

      Asq=(as/(3d0*pi))**2/vevsq

C   Deal with Higgs decay to b-bbar
      s34=s(3,4)+2d0*mb**2
      hdecay=xn*gwsq*mbot**2/(4d0*wmass**2)*2d0*(s34-4d0*mb**2)
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
      if (approxim.ne.0) then
        fac=ason2pi*Asq*gsq*hdecay*sigmat
      elseif (approxim.eq.0) then
        fac=ason2pi*Asq*gsq*hdecay
      endif
      call hjetfill(ss,tt,uu,virtgg,virtqa,virtaq,virtqg,virtgq)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      if ((j.eq.0).and.(k.eq.0)) msq(j,k)=avegg*fac*virtgg
      if ((j.gt.0).and.(k.eq.-j)) msq(j,k)=aveqq*fac*virtqa
      if ((j.lt.0).and.(k.eq.-j)) msq(j,k)=aveqq*fac*virtaq
      if ((j.eq.0).and.(k.ne.0)) msq(j,k)=aveqg*fac*virtgq
      if ((j.ne.0).and.(k.eq.0)) msq(j,k)=aveqg*fac*virtqg
      enddo
      enddo

      return
      end
