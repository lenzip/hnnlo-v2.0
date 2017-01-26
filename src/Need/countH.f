CC    Counterterm to be subtracted from real+virt to get a finite
CC    cross section at qt->0

CC    Version allowing to isolate qg+gq contribution

CC    NEW: full scale dependence included (to be checked !)

CC    NEW: Q2,qt2 and shat obtained trough a common block:
CC         no need to generate real event every time

      double precision function countint(vector,wgt)
      implicit none
      include 'constants.f'
      include 'realonly.f'
      include 'virtonly.f'
      include 'noglue.f'
      include 'vegas_common.f'
      include 'ptilde.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'zerowidth.f'
      include 'efficiency.f'
      include 'masses.f'
      include 'limits.f'
C
      include 'jetlabel.f'
      include 'qcdcouple.f'
      include 'phasemin.f'
      include 'rescoeff.f'
C
      integer ih1,ih2,j,nd,nvec,order
      double precision vector(mxdim),val,xint
      double precision sqrts
      double precision p(mxpart,4),pjet(mxpart,4)
      double precision pswt,rscalestart,fscalestart
      double precision s(mxpart,mxpart),wgt
      double precision msqc(-nf:nf,-nf:nf),xmsq(0:maxd)
      double precision BrnRat,xreal,xreal2
CC
      logical cuts
      double precision ptrans(mxpart,4)
      double precision q2,qt2,shat,Itilde
      double precision fx10(-nf:nf),fx20(-nf:nf)
      double precision fx1p(-nf:nf),fx2p(-nf:nf)
      double precision alfa,beta,diff,Pggreg,Pgq,Cgq
      double precision xjacq2,xjacqt2,xjcut,xth,x3,almin,almax
      double precision xmio,qtcut,fluxborn
      double precision shad,zmax,tauh,Vol,pswt0
      double precision xx0(2),xx10,xx20
      double precision sig1,sig2,LR,LF
      double precision sig11,sig12
      double precision sig21,sig22,sig23,sig24
      double precision tdelta,tH1st,tH1stF,tgaga,tcga,tgamma2
      double precision LL1,LL2,LL3,LL4
      double precision z1,z2,diff10,diff20,diff1f,diff2f
      double precision diffg10,diffg20,diffc10,diffc20,cut
      double precision diffg1f,diffg2f,diffc1f,diffc2f
      double precision D0int,D1int
      double precision Pggggreg,Pgggq,Pgqqg,Pgqqq
      double precision CgqPqq,CgqPqg
      double precision msqb,abisqnnlo,mtop,mbot,abisq
      double precision P2gq,P2gg,sigmatb,sigmat
      external Itilde,Pggreg,Pgq,Cgq,D0int,D1int,
     &         Pggggreg,Pgggq,Pgqqg,Pgqqq,CgqPqq,CgqPqg

      common/xmio/xmio
      common/qtcut/qtcut
      common/xx0/xx0
      common/nnlo/order
      common/tbmass/mtop,mbot
      common/Born/sigmatb,sigmat
CC
CC    Variables passed from virtint or lowint
CC
      common/count/qt2,q2,shat

      integer higgsdec,ndec,approxim
      common/higgsdec/higgsdec,ndec
      integer n2,n3,flgq
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3
      common/xreal/xreal,xreal2
      common/flag/approxim
      logical bin,first
      logical incldip(0:maxd)
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/BrnRat/BrnRat
      common/incldip/incldip

      data first/.true./
      save first,rscalestart,fscalestart
      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale
      endif
      ntotshot=ntotshot+1
      pswt=0d0
      countint=0d0      


      do nd=0,1
      xmsq(nd)=0d0
      enddo
      

CC Check if q2 is the proper interval

      if(q2.lt.wsqmin.or.q2.gt.wsqmax) goto 999   
      
CC Cut on qt2

      if(dsqrt(qt2).lt.qtcut) goto 999


      shad=sqrts**2

      xmio=dsqrt(qt2/q2)

      npart=ndec+1
      nvec=npart+2
      

CC   LR,LF
    
      LR=dlog(q2/scale**2)
      LF=dlog(q2/facscale**2) 


      xth=vector(2)
      x3=vector(1)
      if(higgsdec.ne.1)then
       xth=vector(12)
       x3=vector(13)
      endif
     
CC   jacobian for qt2

    
      xjacqt2=(shat-q2)**2/shat*dabs(1-2*xth)


CC   jacobian for Q2

      if (zerowidth) then
      xjacq2=pi*hmass*hwidth
      else
      almin=datan(-hmass/hwidth)
      almax=datan((shat-hmass**2)/hmass/hwidth)
      xjacq2=hmass*hwidth
     & *(almax-almin)/(dcos((almax-almin)*x3+almin))**2
      endif



CC   Volume

      zmax=1+2*xmio**2-2*xmio*dsqrt(1+xmio**2)
      tauh=q2/shad
      Vol=-(dlog(zmax)-dlog(tauh))/dlog(taumin)     
   

  

CC   LL1,LL2,LL3,LL4: large log (squared) corresponding to eq. (136) 
CC   In this way normalization is fixed to dsigma/dqt2


      LL1=Itilde(1)/q2**2
      LL2=Itilde(2)/q2**2
      LL3=Itilde(3)/q2**2
      LL4=Itilde(4)/q2**2

  

CC Generate BORN momenta for counterterm
      
      if(higgsdec.eq.1) then
       call genBORN2(q2,shat,vector,ptrans,pswt0,*999)
      else
       call genBORN4(q2,shat,vector,ptrans,pswt0,*999)   
      endif


      call storeptilde(1,ptrans)

    

CC Here we have to check if the counterevent passes the cuts

       jets=0
       incldip(1)=cuts(ptrans,0)
       if (incldip(1)) goto 999


CC Compute Born matrix element


      if(higgsdec.eq.1 .and. approxim.eq.0) then
        call gg_h(ptrans,msqc)
       elseif(higgsdec.eq.1 .and. approxim.ne.0)then
        call gg_hexact(ptrans,msqc)
      elseif(higgsdec.eq.2 .and. approxim.eq.0) then
       call qqb_hww(ptrans,msqc)
      elseif(higgsdec.eq.2 .and. approxim.ne.0)then
       call qqb_hwwexact(ptrans,msqc)
      elseif(approxim.eq.0)then
       call qqb_hzz(ptrans,msqc)
      elseif(approxim.ne.0)then
       call qqb_hzzexact(ptrans,msqc)
      endif
  
      msqb=msqc(0,0)



      do nd=0,ndmax
      xmsq(nd)=0d0
      enddo

C Scaled momentum fractions

      cut=1d-8
   
      if(higgsdec.eq.1) then
       beta=cut+(1-cut)*vector(8)
       alfa=cut+(1-cut)*vector(9)  
      else
       beta=cut+(1-cut)*vector(14)
       alfa=cut+(1-cut)*vector(15)
      endif

      xx10=xx0(1)
      xx20=xx0(2)

      z1=xx10**beta
      z2=xx20**alfa


c    To be included later

      xjcut=(1-cut)**2

            
c--- calculate PDF's  

      call fdist(ih1,xx10,facscale,fx10)
      call fdist(ih2,xx20,facscale,fx20)

      call fdist(ih1,xx10**(1-beta),facscale,fx1p)
      call fdist(ih2,xx20**(1-alfa),facscale,fx2p)




CC Flux for Born cross section


   
       fluxborn=fbGeV2/(2*q2)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC Start construction of the counterterm

        flgq=1
        if(gqonly)flgq=0

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

        tdelta=0d0
        tH1st=0d0
        tH1stF=0d0
        tgaga=0d0
        tcga=0d0
        tgamma2=0d0


        sig11=0d0
        sig12=0d0
        sig21=0d0      
        sig22=0d0
        sig23=0d0
        sig24=0d0


CC    First diagonal part, then loop over j=1,nf

      if(gqonly)goto 71

C     Simplest term without convolutions
  
      tdelta=tdelta+fx10(0)*fx20(0)*msqb     



C     H1st delta term

      tH1st=tH1st+2*C1ggdelta*fx10(0)*fx20(0)*msqb
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

 71   if(order.eq.1) goto 74

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


 74   continue


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           
      do j=1,nf


C     H1st: Cgq, first leg


      tH1st=tH1st+(fx1p(j)+fx1p(-j))*Cgq(z1)
     & *(-dlog(xx10))*fx20(0)*msqb


C     H1st: Cgq, second leg

      
      tH1st=tH1st+(fx2p(j)+fx2p(-j))*Cgq(z2)
     & *(-dlog(xx20))*fx10(0)*msqb
    

C     H1st: muf dependence: Pgq, first leg


      tH1stF=tH1stF+(-dlog(xx10))
     & *(fx1p(j)+fx1p(-j))*Pgq(z1)*fx20(0)*msqb


C     H1st: muf dependence: Pgq, second leg


      tH1stF=tH1stF+(-dlog(xx20))
     & *(fx2p(j)+fx2p(-j))*Pgq(z2)*fx10(0)*msqb
  

CC    End of H1st

      if(order.eq.1) goto 75

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


 75   continue


      enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

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
      

CC   First order

      sig12=-0.5d0*A1g*tdelta
      sig11=-B1g*tdelta-tH1stF


CC   Second order

       sig24=(A1g)**2/8*tdelta
       
       sig23=-beta0*A1g/3*tdelta-0.5d0*A1g*sig11

       sig22=0.5d0*(beta0*A1g*LR-A2g)*tdelta
     &     -0.5d0*A1g*(tH1st+LF*tH1stF)
     &     -0.5d0*(B1g-beta0)*sig11
     &     +0.5d0*B1g*tH1stF
     &     +0.5d0*tgaga
C
C    Add mur dependence from H1st
C
     &     +beta0*A1g*LR*tdelta


       sig21=-beta0*LR*sig11-B1g*(tH1st+LF*tH1stF)
     &     -LF*tgaga-B2g*tdelta+beta0*tH1st-tcga-tgamma2



c     Include missing delta term from C*gamma (no factor 2 here !)

       sig21=sig21-C1ggdelta*tH1stF


C     Include missing term from contact term in 2 loop AP

       sig21=sig21-2*Delta2gg*tdelta

C     NEW: Include mur dependence from H1st

       sig21=sig21+2*beta0*LR*(B1g*tdelta+tH1stF)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CC Include as/pi factors and sum O(as) and O(as^2) contributions

           sig1=sig12*LL2+sig11*LL1
           sig2=sig24*LL4+sig23*LL3+sig22*LL2+sig21*LL1
         

           sig1=sig1*ason2pi*2
           sig2=sig2*(ason2pi*2)**2

           if(order.eq.1)then
            xmsq(1)=-sig1
            else 
              xmsq(1)=-(sig1+sig2)
           endif

           if (order.eq.2 .and. approxim.ne.0) then
            sig2=sig2*sigmat
            sig2=sig2/sigmatb
            xmsq(1)=-(sig1+sig2)
           endif


CC Include iacobians

      xmsq(1)=xmsq(1)*xjacqt2*xjacq2*q2/shad/Vol


      countint=0d0
      xint=0d0


C Multiply by BORN phase space weight

        xmsq(1)=xmsq(1)*fluxborn*pswt0/BrnRat/2d0
      

c---Add to total

        xint=xmsq(1)
        val=xmsq(1)*wgt
                

c---if we're binning, add to histo too
        if (bin) then
          call getptildejet(1,pjet)
          call dotem(nvec,pjet,s)
          val=val/dfloat(itmx)
           call plotter(ptrans,val,1)        
C           call plotter(p,val,0)
        endif

     
      countint=xint
     
      xreal=xreal+xint*wgt/dfloat(itmx)
      xreal2=xreal2+(xint*wgt)**2/dfloat(itmx)
      


      return

 999  countint=0d0
      ntotzero=ntotzero+1
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC gg splitting function (not normalized !)

      function Pgg(x)
      implicit none
      real *8 Pgg,x

      Pgg=1d0/(1-x)+1d0/x-2+x*(1-x)

      return
      end


CC gg splitting function: regular part (with asopi normalization)


      function Pggreg(z)
      implicit none
      real *8 Pggreg,z
      Pggreg=3*((1-2*z)/z+z*(1-z))
      return
      end

CC gq splitting function (with asopi normalization)

      function Pgq(z)
      implicit none
      real *8 Pgq,z
      Pgq=2d0/3*(1+(1-z)**2)/z
      return
      end




CC Cgq coefficient (with asopi normalization)

      function Cgq(z)
      implicit none
      real *8 Cgq,z
      Cgq=2d0/3*z
      return
      end



CC Integral of 1/(1-x) from 0 to z

      function D0int(z)
      implicit none
      real *8 D0int,z
      D0int=-dlog(1-z)
      return
      end

CC Integral of log(1-x)/(1-x) from 0 to z

      function D1int(z)
      implicit none
      real *8 D1int,z
      D1int=-0.5d0*dlog(1-z)**2
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                P*P convolutions
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      function Pggggreg(x) ! checked !
      implicit none
      real *8 Pggggreg,x
      real *8 Pggreg,beta0
      integer nf
      external Pggreg

      nf=5
      beta0=(33-2*nf)/12d0

      Pggggreg=-9*dlog(x)/(1-x)+2*beta0*Pggreg(x)
     &   +9*(3*(1-x)+11d0/3/x*(x**3-1d0)+2d0/3*dlog(1-x)*Pggreg(x)
     &   +dlog(x)*(x**2-3*x-1d0/x))

      return
      end


      function Pgggq(x) ! checked !
      implicit none
      real *8 Pgggq,x
      real *8 Pgq,beta0
      integer nf
      external Pgq

      nf=5
      beta0=(33-2*nf)/12d0

      Pgggq=2*((1+(1-x)**2)/x*dlog(1-x)-2*(1+x+1d0/x)*dlog(x)
     &    +4-31d0/6/x+x/2+2d0/3*x**2)+beta0*Pgq(x)

      return
      end

      function Pgqqq(x) ! checked !
      implicit none
      real *8 Pgqqq,x
      
      Pgqqq=4d0/9*((2-x)*dlog(x)+dlog(1-x)*(2*x+4/x-4)+2-x/2)
      
      return
      end

      function Pgqqg(x) ! checked !
      implicit none
      real *8 Pgqqg,x
      integer nf

      nf=5 
      Pgqqg=1d0/6*(1+4d0/3/x-x-4*x**2/3+2*(1+x)*dlog(x))

      Pgqqg=2*nf*Pgqqg

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                C*P convolutions
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      function CgqPqq(x) ! checked !
      implicit none
      real *8 CgqPqq,x
      
      CgqPqq=2d0/9*(2+x+4*x*dlog(1-x)-2*x*dlog(x))

      return
      end



      function CgqPqg(x) ! checked !
      implicit none
      real *8 CgqPqg,x
      integer nf

      nf=5       
      CgqPqg=1d0/6*(1+x-2*x**2+2*x*dlog(x))

      CgqPqg=2*nf*CgqPqg

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C           Two loop AP:  pqq of ESW is my 3/2 Pqq
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      function P2gg(x)
      implicit none
      real *8 P2gg,S2,pgg,x,pi
      integer nf
      external pgg,S2

      nf=5

      pi=3.14159265358979d0

      P2gg=2d0/3*nf*(-16+8*x+20d0/3*x**2+4d0/3/x-(6+10*x)*dlog(x)
     &    -(2+2*x)*dlog(x)**2)+1.5d0*nf*(2-2*x+26d0/9*(x**2-1d0/x)
     &    -4d0/3*(1+x)*dlog(x)-20d0/9*(1/x-2+x*(1-x)))
     &    +9*(27d0/2*(1-x)+67d0/9*(x**2-1d0/x)
     &        -(25d0/3-11d0/3*x+44d0/3*x**2)*dlog(x)
     &        +4*(1+x)*dlog(x)**2+2*pgg(-x)*S2(x)
     &    +(67d0/9-pi**2/3)*(1/x-2+x*(1-x))
     &     +(-4*dlog(x)*dlog(1-x)+dlog(x)**2)*pgg(x))
     

      P2gg=P2gg/4d0

      return
      end


      function P2gq(x)
      implicit none
      real *8 P2gq,Pgq,S2,x,logx,logomx,pi
      external Pgq,S2
      integer nf
     
      pi=3.14159265358979d0

      nf=5

      logx=dlog(x)
      logomx=dlog(1-x)

      P2gq=16d0/9*(-2.5d0-3.5d0*x+(2+3.5d0*x)*logx
     &     -(1-0.5d0*x)*logx**2-2*x*logomx
     &     -(3*logomx+logomx**2)*1.5d0*Pgq(x))
     &     +4d0*(28d0/9+65d0/18*x+44d0/9*x**2-(12+5*x+8d0/3*x**2)*logx
     &     +(4+x)*logx**2+2*x*logomx+S2(x)*1.5d0*Pgq(-x)
     &     +(0.5d0-2*logx*logomx+0.5d0*logx**2+11d0/3*logomx+logomx**2
     &     -pi**2/6)*1.5d0*Pgq(x))
     &     +2d0/3*nf*(-4d0*x/3-(20d0/9+4d0/3*logomx)*1.5d0*Pgq(x))


      P2gq=P2gq/4

      return
      end



C    S2: Eq. (4.114) ESW

      function S2(x)
      implicit none
      real *8 x,pi,S2,myli2
      external myli2      
      pi=3.14159265358979d0

      S2=-2*myli2(-x)+0.5d0*dlog(x)**2-2*dlog(x)*dlog(1+x)-pi**2/6
      return
      end
