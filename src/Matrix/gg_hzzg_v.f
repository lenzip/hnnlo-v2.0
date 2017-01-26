CC    NEW: gg->Hg virtual contribution with H->ZZ decay


      subroutine gg_hzzg_v(p,msq)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scheme.f'
      include 'zcouple.f'

C     (Taken from Ravindran, Smith, van Neerven hep-ph/0201114)
C     Modified by overall factors
      integer iglue,j,k
      double precision p(4,mxpart),msq(fn:nf,fn:nf)
      double precision ss,tt,uu,s34,
     . virtgg,virtqa,virtaq,virtqg,virtgq,hdecay,Asq,fac,shiggs
      double precision dec,decay,interf,num,den,resfac
      common/Rescale/resfac

      parameter(iglue=7)  ! gluon label

      logical int
      common/int/int

      scheme='tH-V'

      call dotem(iglue,p,s)
      ss=s(1,2)
      tt=s(1,iglue)
      uu=s(2,iglue)

      Asq=(as/(3d0*pi))**2/vevsq



      shiggs=s(3,4)+s(3,5)+s(3,6)+s(4,5)+s(4,6)+s(5,6)
      interf=0d0

    

      decay=(((l1*l2)**2+(r1*r2)**2)*s(3,5)*s(4,6)
     .  +((r1*l2)**2+(r2*l1)**2)*s(3,6)*s(4,5))


      decay=decay/((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)
      decay=decay/((s(5,6)-zmass**2)**2+(zmass*zwidth)**2)

C     Here only H->ZZ->(34)(56): diagram with (1<->3) accounted for
C     by adding a factor 2

      if(int.eqv..false.)goto 39

      decay=2*decay


C Interference contribution (check if factor 2 is there !)

      interf=2*((l1*l2)**2+(r1*r2)**2)*s(3,5)*s(4,6)

      num=((s(3,4)-zmass**2)*(s(5,6)-zmass**2)*(s(4,5)-zmass**2)*
     .   (s(3,6)-zmass**2)+(zmass*zwidth)**4+
     .   (zmass*zwidth)**2*(2*zmass**4-zmass**2*
     .   (s(3,4)+s(5,6)+s(4,5)+s(3,6))+s(3,4)*s(3,6)+s(3,4)*s(4,5)+
     .   s(3,6)*s(5,6)+s(4,5)*s(5,6)-s(3,6)*s(4,5)-s(3,4)*s(5,6)))

      den=((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)*
     .   ((s(5,6)-zmass**2)**2+(zmass*zwidth)**2)*
     .   ((s(4,5)-zmass**2)**2+(zmass*zwidth)**2)*
     .   ((s(3,6)-zmass**2)**2+(zmass*zwidth)**2)


      interf=interf*num/den
            
 39   continue

      dec=gwsq**3*zmass**2*4d0*xw**2/(one-xw)*
     .      (decay+interf)/((shiggs-hmass**2)**2+(hmass*hwidth)**2)


C     In case of identical particles add 1/4 symmetry factor

      if(int) dec=dec/4

      fac=ason2pi*Asq*gsq*dec*resfac

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
