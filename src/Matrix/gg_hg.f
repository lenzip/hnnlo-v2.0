       subroutine gg_hg(p,msq)
       implicit none
c---Matrix element squared averaged over initial colors and spins
c     f(-p1) + f(-p2) --> H + f(p5)
c                         |
c                         --> b(p3)+bbar(p4)
c                            
c--all momenta incoming
c
c--- Matrix elements are taken from:
c--- R.~K.~Ellis, I.~Hinchliffe, M.~Soldate and J.~J.~van der Bij,
c--- %``Higgs Decay To Tau+ Tau-: A Possible Signature Of Intermediate
c--- % Mass Higgs Bosons At The SSC,''
c--- Nucl.\ Phys.\ B {\bf 297}, 221 (1988).
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),hbb
      double precision ehsvm3,ehsvm4,origmbsq,gg,qg,gq,qq
      double precision mtop,mbot
      common/tbmass/mtop,mbot
  
      call dotem(5,p,s)

c--- matrix element for H -> bbbar

      if(s(1,2).eq.0d0.or.s(1,5).eq.0d0.or.s(2,5).eq.0d0)then
        gg=0d0
        qq=0d0
        qg=0d0
        gq=0d0
      else
      hbb=xn*gwsq*0.5d0*mbot**2/wmass**2*(s(3,4)-2d0*mb**2)/
     .     ((s(3,4)-hmass**2+2d0*mb**2)**2+(hmass*hwidth)**2)

CH      origmbsq=mbsq
CH      mbsq=mt**2

      gg=+avegg*ehsvm3(s(1,2),s(1,5),s(2,5))*hbb
      qq=+aveqq*ehsvm4(s(1,2),s(1,5),s(2,5))*hbb
      qg=-aveqg*ehsvm4(s(1,5),s(1,2),s(2,5))*hbb
      gq=-aveqg*ehsvm4(s(2,5),s(1,5),s(1,2))*hbb
      endif
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
           msq(j,k)=qq
      endif

      enddo
      enddo

      return
      end
      
      
      double precision function ehsvm3(s,t,u)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
c---Matrix element squared Eqn 2.2 of EHSV
      double complex ehsva2,ehsva4
      double precision s,t,u
      double precision mtop,mbot,fac
      common/tbmass/mtop,mbot
      integer approxim
      common/flag/approxim
c--- approx TRUE uses the heavy fermion approximation to Msq


      if (approxim.eq.0) then
      ehsvm3=(as/(3d0*pi))**2/vevsq*gsq*xn*V*(
     .        hmass**8+s**4+t**4+u**4)/s/t/u
      else
      ehsvm3=
     .  abs(ehsva2(s,t,u,mbot**2,mtop**2))**2
     . +abs(ehsva2(u,s,t,mbot**2,mtop**2))**2
     . +abs(ehsva2(t,u,s,mbot**2,mtop**2))**2
     . +abs(ehsva4(s,t,u,mbot**2,mtop**2))**2 
      ehsvm3=(as/pi)**2/vevsq*gsq*xn*V*hmass**8
     .           /(s*t*u)*ehsvm3

      endif
      
      return
      end
      
      
      double precision function ehsvm4(s,t,u)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
c---Matrix element squared Eqn 2.6 of EHSV
      double complex ehsva5
      double precision s,t,u
      double precision mtop,mbot,fac
      common/tbmass/mtop,mbot
      integer approxim
      common/flag/approxim
c--- approx TRUE uses the heavy fermion approximation to Msq


      if (approxim.eq.0) then
      ehsvm4=(as/(3d0*pi))**2/vevsq*gsq*V/2d0*(u**2+t**2)/s
      else
      ehsvm4=abs(ehsva5(s,t,u,mbot**2,mtop**2))**2
      ehsvm4=(as/(2d0*pi))**2/vevsq*gsq*V/2d0*(u**2+t**2)
     .  /s
     . *hmass**4/(u+t)**2*ehsvm4
      

      endif
      
      return
      end
