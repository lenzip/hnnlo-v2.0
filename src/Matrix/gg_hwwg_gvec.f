      subroutine gg_hwwg_gvec(p,n,in,msq)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
C  in is the label of the momentum contracted with n
      integer j,k,in
      double precision msq(-nf:nf,-nf:nf)
      double precision mtop,mbot
      common/tbmass/mtop,mbot
      double precision n(4),p(mxpart,4),dot,hdecay,sh,fac,
     . qqghn,ggghn,p1p2(-1:1,-1:1)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

C   Deal with Higgs decay to b-bbar
c      s34=2d0*Dot(p,3,4)+2d0*mb**2
c      hdecay=xn*gwsq*mbot**2/(4d0*wmass**2)*2d0*(s34-4d0*mb**2) 
c      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)


C     Higgs virtuality

c      sh=s(3,4)+s(3,5)+s(3,6)+s(4,5)+s(4,6)+s(5,6)

      sh=2*(Dot(p,3,4)+Dot(p,3,5)+Dot(p,3,6)
     &     +Dot(p,4,5)+Dot(p,4,6)+Dot(p,5,6))

      hdecay=gwsq**3*wmass**2*(4*Dot(p,3,5)*Dot(p,4,6))
      hdecay=hdecay/((2*Dot(p,3,4)-wmass**2)**2+(wmass*wwidth)**2)
      hdecay=hdecay/((2*Dot(p,5,6)-wmass**2)**2+(wmass*wwidth)**2)
      hdecay=hdecay/((sh-hmass**2)**2+(hmass*hwidth)**2)
      


      fac=hdecay

      do j=-1,+1
      do k=-1,+1
      p1p2(j,k)=0d0
      enddo
      enddo

      if (in .eq. 1) then
      p1p2(0,-1)=-aveqg*fac*qqghn(2,7,1,p,n)
      p1p2(0,+1)=-aveqg*fac*qqghn(2,7,1,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(7,2,1,p,n)
      elseif (in .eq. 2) then
      p1p2(+1,0)=-aveqg*fac*qqghn(1,7,2,p,n)
      p1p2(-1,0)=-aveqg*fac*qqghn(7,1,2,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(1,7,2,p,n)
      elseif (in .eq. 7) then     
      p1p2(1,-1)=+aveqq*fac*qqghn(1,2,7,p,n)
      p1p2(-1,1)=+aveqq*fac*qqghn(2,1,7,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(1,2,7,p,n)
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


