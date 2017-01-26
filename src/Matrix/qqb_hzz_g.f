      subroutine qqb_hzz_g(p,msq)
      implicit none
c----NLO matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H -->  Z(e^-(p3)+e^+(p4)) + Z(mu^-(p5)+mu^+(p6))
c    +g(p7)
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision dec,interf,num,den,shiggs
      double precision sh,ss,tt,uu,decay,s(mxpart,mxpart)
      double precision aw,qqb,qg,gq,gg,ehsvm3,ehsvm4

      logical int
      common/int/int

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      aw=gwsq/(4d0*pi)
      call dotem(7,p,s)


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

C     Interference contribution

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

c--   calculate propagators
      ss=s(1,2)
      tt=s(1,7)
      uu=s(2,7)
      sh=s(1,2)+s(1,7)+s(2,7)
      if(ss.eq.0d0.or.tt.eq.0d0.or.uu.eq.0d0)then
       gg=0d0
       qqb=0d0
       gq=0d0
       qg=0d0
      else

ch      gg=aw*as**3*4d0*V/9d0*xn*(sh**4+ss**4+tt**4+uu**4)
ch     . /(ss*tt*uu*wmass**2)
ch      qqb=aw*as**3*2d0*V/9d0*(tt**2+uu**2)/(ss*wmass**2)
ch      gq=-aw*as**3*2d0*V/9d0*(ss**2+tt**2)/(uu*wmass**2)
ch      qg=-aw*as**3*2d0*V/9d0*(ss**2+uu**2)/(tt*wmass**2)

      gg=ehsvm3(ss,tt,uu)
      qqb=ehsvm4(ss,tt,uu)
      gq=-ehsvm4(uu,tt,ss)
      qg=-ehsvm4(tt,ss,uu)


      gg=avegg*gg*dec
      gq=aveqg*gq*dec
      qg=aveqg*qg*dec
      qqb=aveqq*qqb*dec
      endif


c--set msq=0 to initialize
ch      do j=-nf,nf
ch      do k=-nf,nf
ch      if ((k .eq. -j) .and. (j .ne. 0)) then
ch      msq(j,k)=qqb
ch      elseif ((j .eq. 0) .and. (k .ne. 0)) then
ch      msq(j,k)=gq
ch      elseif ((j .ne. 0) .and. (k .eq. 0)) then
ch      msq(j,k)=qg
ch      elseif ((k .eq. 0) .and. (j .eq. 0)) then
ch      msq(j,k)=gg
ch      endif
ch      enddo
ch      enddo
ch      return
ch      end


CH      mbsq=origmbsq

      do j=-nf,nf    
      do k=-nf,nf
      msq(j,k)=0d0

      if ((j.eq. 0) .or. (k.eq.0)) then
           if ((j.eq. 0) .and. (k.eq.0)) then
                msq(j,k)=gg
           elseif ((j.eq.0).and.(k.ne.0)) then
                msq(j,k)=gq
           elseif ((j.ne.0).and.(k.eq.0)) then
                msq(j,k)=qg
           endif
      elseif ((j.eq.-k).and. (j.ne.0)) then
           msq(j,k)=qqb
      endif

      enddo
      enddo

      return
      end
