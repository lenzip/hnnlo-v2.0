      subroutine gg_hwwgg(p,msq)
      implicit none

CC    NEW from gg_hgg: replace H->bbar into H->WW

c---Matrix element squared averaged over initial colors and spins
c
c     g(-p1)+g(-p2) --> H(nu(p3)+e^+(p4)+e^-(p5)+nubar(p6))+g(p7)+g(p8)

      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j,k,iglue1,iglue2
      double precision p(mxpart,4),Asq,fac,sh
      double precision Hgggg,Hqagg,Haqgg,Hgqgq,Hgaga,Hqgqg,Hagag,Hggqa
      double precision 
     . Hqrqr,Hqqqq,
     . Habab,Haaaa,
     . Hqarb,Hqaqa,Hqbqb,
     . Haqbr,Haqaq,Hbqbq
      double precision msq(-nf:nf,-nf:nf),hdecay,s34
      double precision mtop,mbot
      common/tbmass/mtop,mbot

      parameter(iglue1=7,iglue2=8)


C---fill spinor products upto maximum number
      call spinoru(iglue2,p,za,zb)  

C   Deal with Higgs decay to b-bbar

c      s34=s(3,4)+2d0*mb**2
c      hdecay=xn*gwsq*mbot**2/(4d0*wmass**2)*2d0*(s34-4d0*mb**2) 
c      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
c      Asq=(as/(3d0*pi))**2/vevsq
c      fac=gsq**2*Asq*hdecay

C     Higgs virtuality

      sh=s(3,4)+s(3,5)+s(3,6)+s(4,5)+s(4,6)+s(5,6)

      hdecay=gwsq**3*wmass**2*s(3,5)*s(4,6)
      hdecay=hdecay/((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
      hdecay=hdecay/((s(5,6)-wmass**2)**2+(wmass*wwidth)**2)
      hdecay=hdecay/((sh-hmass**2)**2+(hmass*hwidth)**2)


      Asq=(as/(3d0*pi))**2/vevsq
      fac=gsq**2*Asq*hdecay

C--four gluon terms
      call h4g(1,2,iglue1,iglue2,Hgggg)

C--two quark two gluon terms
      call hqqgg(1,2,iglue1,iglue2,Hqagg)
c      call hqqgg(2,1,iglue1,iglue2,Haqgg)
C====symmetric in first two arguments
      Haqgg=Hqagg

      call hqqgg(1,iglue1,2,iglue2,Hqgqg)
c      call hqqgg(iglue1,1,2,iglue2,Hagag)
C====symmetric in first two arguments
      Hagag=Hqgqg

      call hqqgg(2,iglue1,1,iglue2,Hgqgq)
c      call hqqgg(iglue1,2,1,iglue2,Hgaga)
C====symmetric in first two arguments
      Hgaga=Hgqgq

      call hqqgg(iglue2,iglue1,1,2,Hggqa)

C---four quark terms
      call H4qn(1,2,iglue1,iglue2,Hqrqr)
      call H4qi(1,2,iglue1,iglue2,Hqqqq)
C---four anti-quark terms
c      call H4qn(iglue1,iglue2,1,2,Habab)
c      call H4qi(iglue1,iglue2,1,2,Haaaa)
      Habab=Hqrqr
      Haaaa=Hqqqq

C-qqb
      call H4qn(1,iglue2,2,iglue1,Hqarb)
      call H4qi(1,iglue2,2,iglue1,Hqaqa)

      call H4qn(1,iglue2,iglue1,2,Hqbqb)


C-qbq
      Haqbr=Hqarb
      call H4qi(2,iglue2,1,iglue1,Haqaq)
      call H4qn(2,iglue2,iglue1,1,Hbqbq)

      do j=fn,nf
      do k=fn,nf
      msq(j,k)=0d0
      if ((j.gt.0).and.(k.gt.0)) then 
      if (j.eq.k) then
      msq(j,k)=0.5d0*aveqq*fac*Hqqqq
      else
      msq(j,k)=aveqq*fac*Hqrqr
      endif
      endif
      if ((j.lt.0).and.(k.lt.0)) then 
      if (j.eq.k) then
      msq(j,k)=0.5d0*aveqq*fac*Haaaa
      else
      msq(j,k)=aveqq*fac*Habab
      endif
      endif

      if ((j.gt.0).and.(k.lt.0)) then
      if (j.eq.-k) then
      msq(j,k)=aveqq*fac*(0.5d0*Hqagg+Hqaqa+(nf-1)*Hqarb)
      else
      msq(j,k)=aveqq*fac*Hqbqb
      endif
      endif

      if ((j.lt.0).and.(k.gt.0)) then
      if (j.eq.-k) then
      msq(j,k)=aveqq*fac*(0.5d0*Haqgg+Haqaq+dfloat(nf-1)*Haqbr)
      else
      msq(j,k)=aveqq*fac*Hbqbq
      endif
      endif

      if ((j.gt.0).and.(k.eq.0)) msq(j,0)=aveqg*fac*Hqgqg
      if ((j.lt.0).and.(k.eq.0)) msq(j,0)=aveqg*fac*Hagag

      if ((j.eq.0).and.(k.gt.0)) msq(0,k)=aveqg*fac*Hgqgq
      if ((j.eq.0).and.(k.lt.0)) msq(0,k)=aveqg*fac*Hgaga

      if ((j.eq.0).and.(k.eq.0)) then
      msq(0,0)=avegg*fac*(0.5d0*Hgggg+dfloat(nf)*Hggqa)
      endif
      
      enddo
      enddo

      return
      end

 
