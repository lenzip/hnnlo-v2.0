      subroutine coupling
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'nlooprun.f'
      include 'process.f'
      include 'ewinput.f'
      include 'b0.f'
      character*4 part,mypart
      common/part/part
      integer i,order,approxim
      double precision aemmz,alphas,amz,cmass,bmass
      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      double precision xsq,topwidth
      double precision mtop,mbot
      character*3 inlabel(10)
      common/cabib/Vud,Vus,Vub,
     &             Vcd,Vcs,Vcb
      common/qmass/cmass,bmass
      common/em/aemmz
      common/couple/amz
      common/mypart/mypart
      common/nnlo/order
      common/tbmass/mtop,mbot
      common/flag/approxim
ch      common/flag/appr


c--- This is the MCFM default, corresponding to an effective
c--- field theory approach valid for scales below the top-mass
C--- (see Georgi, Nucl. Phys. B 363 (1991) 301).
c--- There are 4 inputs here instead of the usual 3 ...

         Gf = Gf_inp
         aemmz  = aemmz_inp
         wmass  = wmass_inp
         zmass  = zmass_inp
        
c--- ... and as result, both xw and mtop are derived
c--- (using wmass=zmass*dsqrt(rho)*cos(theta_w)
c---    and rho=1+3d0*aemmz/16d0/pi/xw*(mt/wmass)**2 )
         xw  = fourpi*aemmz/(8d0*wmass**2*Gf/rt2)

C      Even when approxim=0 mtop is taken from input

C      if (approxim.eq.0) then
C         mtop  = dsqrt(16d0*pisq/3d0/rt2/Gf*(
C     .          wmass**2/zmass**2/(1d0-xw)-1d0))
C      endif
     

c--- Now set up the other derived parameters
      gwsq=fourpi*aemmz/xw
      esq=gwsq*xw
      gw=dsqrt(gwsq)
      call couplz(xw)

c--- Calculate the appropriate Higgs vacuum expectation value.
c--- This vevsq is defined so that gwsq/(4*wmass**2)=Gf*rt2=1/vevsq
c--- (ie differs from definition in ESW)
      vevsq=1d0/rt2/Gf

c--- Set-up twidth, using LO formula except when including radiation in decay
      xsq=(wmass/mtop)**2

c--- set up the beta-function
      b0=(xn*11d0-2d0*nf)/6d0

c--- initialize the pdf set
      nlooprun=0
      call pdfset      

      cmass=dsqrt(mcsq)
      bmass=mbot
      musq=scale**2
 
c--- set the number of loops to use in the running of alpha_s
c--- if it hasn't been set by pdfset already
      if (nlooprun.eq.0) then
       nlooprun=order+1
      endif

c--- initialize alpha_s
      as=alphas(abs(scale),amz,nlooprun)

      ason2pi=as/twopi
      ason4pi=as/fourpi
      gsq=fourpi*as

***************************************
      write(6,99)gsq/fourpi,amz,nlooprun
      

 99   format(' CCCCCCCCCCC  Strong coupling, alpha_s   CCCCCCCCCCCC'/,
     .      ' C                                                  C'/,
     .      ' C   alpha_s (mur)=',f8.5,'                         C'/,                 
     .      ' C   alpha_s (mz) =',f8.5,'                         C'/,     
     .      ' C   (using ',i1,
     .      '-loop running for alpha_s)             C'/, 
     .      ' C                                                  C') 
    

      
      return
      end


      block data wsalam1
      implicit none
      include 'constants.f'
      include 'ewcharge.f'
      data Q(-5)/+0.333333333333333d0/
      data Q(-4)/-0.666666666666667d0/
      data Q(-3)/+0.333333333333333d0/
      data Q(-2)/-0.666666666666667d0/
      data Q(-1)/+0.333333333333333d0/
      data Q(0)/+0d0/
      data Q(+1)/-0.333333333333333d0/
      data Q(+2)/+0.666666666666667d0/
      data Q(+3)/-0.333333333333333d0/
      data Q(+4)/+0.666666666666667d0/
      data Q(+5)/-0.333333333333333d0/
      data tau/1d0,-1d0,1d0,-1d0,1d0,0d0,-1d0,1d0,-1d0,1d0,-1d0/
      end 



