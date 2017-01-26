      subroutine gg_hzzg_gvec(p,n,in,msq)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
C  in is the label of the momentum contracted with n
      integer j,k,in,iglue
      double precision msq(-nf:nf,-nf:nf),s,shiggs,interf,num,den
      double precision n(4),p(mxpart,4),dot,decay,dec,fac,
     . qqghn,ggghn,p1p2(-1:1,-1:1)
      parameter(iglue=7)

      logical int
      common/int/int

      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

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


C Interference contribution

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

  
      fac=dec

      do j=-1,+1
      do k=-1,+1
      p1p2(j,k)=0d0
      enddo
      enddo

      if (in .eq. 1) then
      p1p2(0,-1)=-aveqg*fac*qqghn(2,iglue,1,p,n)
      p1p2(0,+1)=-aveqg*fac*qqghn(2,iglue,1,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(iglue,2,1,p,n)
      elseif (in .eq. 2) then
      p1p2(+1,0)=-aveqg*fac*qqghn(1,iglue,2,p,n)
      p1p2(-1,0)=-aveqg*fac*qqghn(iglue,1,2,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(1,iglue,2,p,n)
      elseif (in .eq. 7) then     
      p1p2(1,-1)=+aveqq*fac*qqghn(1,2,iglue,p,n)
      p1p2(-1,1)=+aveqq*fac*qqghn(2,1,iglue,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(1,2,iglue,p,n)
      endif

      do j=-nf,nf
      do k=-nf,nf
      if     ((j .gt. 0) .and. (k .eq. -j)) then
          msq(j,k)=p1p2(1,-1)
      elseif ((j .lt. 0) .and. (k .eq. -j)) then
          msq(j,k)=p1p2(-1,1)
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
          msq(j,k)=p1p2(0,0)
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &    p1p2(+1,0)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &    p1p2(-1,0)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=
     &    p1p2(0,+1)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=
     &    p1p2(0,-1)
      endif
      enddo
      enddo
 
      return
      end

