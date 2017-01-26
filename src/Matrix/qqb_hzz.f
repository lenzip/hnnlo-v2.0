      subroutine qqb_hzz(p,msq)
      implicit none
C----Lowest order matrix element for H production
C----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
C     g(-p1)+g(-p2)->H->ZZ->(mu-(p3)+mu+(p4)+e-(p5)+e+(p6))

C    Includes interference contribution in case mu->e


      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      integer j,k
      double precision abisq,msq(-nf:nf,-nf:nf),p(mxpart,4),s,shiggs
      double precision dec,decay,gg,Asq,interf,num,den
      logical int
      common/int/int

      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

c---set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

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

CH--  multiplication with the abisq
      Asq=(as/(3d0*pi))**2/vevsq
ch      Asq=(as/(3d0*pi))**2/vevsq*abisq(mt**2/(hmass**2),mbsq/(hmass**2))
      gg=0.5d0*V*Asq*shiggs**2

c---calculate propagators
      msq(0,0)=avegg*gg*dec


      return
      end
