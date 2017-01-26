      subroutine bookplot(n,tag,titlex,var,wt,xmin,xmax,dx,llplot) 
      implicit none
      include 'nplot.f'
      integer n,higgsdec
      character*(*) titlex
      character*3 llplot
      character*4 tag,mypart
      double precision var,wt,xmin,xmax,dx
      common/mypart/mypart

      if (tag.eq.'book') then
          call mbook(n,titlex,dx,xmin,xmax)
          call mbook(20+n,titlex,dx,xmin,xmax)
          call mbook(40+n,titlex,dx,xmin,xmax)
          call mbook(60+n,titlex,dx,xmin,xmax)
          call mbook(80+n,titlex,dx,xmin,xmax)
          call mbook(100+n,titlex,dx,xmin,xmax)
          call mbook(120+n,titlex,dx,xmin,xmax)
      elseif (tag .eq. 'plot') then
          call mfill(n,var,wt)
        linlog(n)=llplot
        titlearray(n)=titlex
      endif

      return
      end

    
      subroutine plotter(p,wt,switch)
      implicit none
      include 'clustering.f'
      include 'constants.f'
      include 'cutoff.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'mxdim.f'
      include 'process.f'
      include 'removebr.f'
      include 'masses.f'

      integer n,switch,nplotmax
      character tag*4
  
      double precision dot
      double precision wt
      double precision p(mxpart,4)
      double precision pt,yraptwo,pttwo
      integer eventpart

      logical first,jetmerge
      common/nplotmax/nplotmax
      common/jetmerge/jetmerge
      integer order,higgsdec,ndec
      common/higgsdec/higgsdec,ndec
      common/nnlo/order
      data first/.true./
      save first

      double precision m45
      double precision ptj1,ptj2,pjm
      double precision pto(1:4)  
      double precision damping
      double precision xjets,xjets_incl
      double precision eta3,pt3,eta4,pt4,eta5,pt5,eta6,pt6
      double precision eta7,pt7,eta8,pt8
      double precision deltaphi,cosphi45,ptmiss,pt45 
      double precision mjj,detajj,cosphiLLMet,mth
      double precision etarap

      wt=wt  !  Convert to Powheg-Box units.

      if (first) then
         tag='book'
c--- ensure we initialize all possible histograms
         eventpart=npart+3
         
         ptj1=0d0
         ptj2=0d0
         pjm=0d0
         m45=0d0
         xjets=0d0
         xjets_incl=0d0
         
         jetmerge=.true.
CC      Here set jets to the maximum number of jets 
CC      to book the necessary histograms: 0 at LO, 1 at NLO and 2 at NNLO
         jets=order
CC
         goto 99
      else
         tag='plot'
      endif


      eta3=etarap(3,p)
      pt3=pt(3,p)
      eta4=etarap(4,p)
      pt4=pt(4,p)

      eta5=etarap(5,p)
      pt5=pt(5,p)
      eta6=etarap(6,p)
      pt6=pt(6,p)

      pt7=pt(7,p)
      pt8=pt(8,p)

      pto(1)=max(pt4,pt5)
      pto(2)=min(pt4,pt5) 



      ptj1=pt(7,p)
      ptj2=pt(8,p)
      if(ptj1.lt.ptj2) then
        pjm=ptj2
        ptj2=ptj1
        ptj1=pjm
      endif   
     
      if(ptj1.gt.0d0.and.ptj2.gt.0d0) then
       xjets=2d0
      elseif(max(ptj1,ptj2).gt.0d0.and.min(ptj1,ptj2).eq.0d0) then
       xjets=1d0
      else
       xjets=0d0
      endif

      mjj=-1d0
      detajj=-1d0
      if (xjets.ge.2d0) then
        mjj=dsqrt((p(7,4)+p(8,4))**2
     &    -(p(7,1)+p(8,1))**2
     &    -(p(7,2)+p(8,2))**2
     &    -(p(7,3)+p(8,3))**2)
        detajj=dabs(etarap(7,p)-etarap(8,p))
      endif

      m45=dsqrt((p(4,4)+p(5,4))**2
     &    -(p(4,1)+p(5,1))**2
     &    -(p(4,2)+p(5,2))**2
     &    -(p(4,3)+p(5,3))**2)

      if (xjets.ge.2d0) then
       if ((mjj.gt.500d0).and.(detajj.gt.3.5d0)) then
          xjets=3d0
       endif
      endif

      cosphi45=(p(4,1)*p(5,1)+p(4,2)*p(5,2))/(pt4*pt5)
      deltaphi=dabs(dacos(cosphi45))
     &      *180/3.141592653589793d0

      ptmiss=dsqrt((p(3,1)+p(6,1))**2+(p(3,2)+p(6,2))**2)
      pt45=dsqrt((p(4,1)+p(5,1))**2+(p(4,2)+p(5,2))**2)
      cosphiLLMet=((p(3,1)+p(6,1))*(p(4,1)+p(5,1)) +
     &             (p(3,2)+p(6,2))*(p(4,2)+p(5,2)) ) /
     &             (ptmiss*pt45) 

      mth=dsqrt(2*ptmiss*pt45*cosphiLLMet)
 

 99   continue
      n=1                  
      call bookplot(n,tag,'eta3',eta3,wt,-4d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt3',pt3,wt,0d0,150d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'eta4',eta4,wt,-4d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt4',pt4,wt,0d0,150d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'eta5',eta5,wt,-5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt5',pt5,wt,0d0,500d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'eta6',eta6,wt,-5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt6',pt6,wt,0d0,500d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'pto1',pto(1),wt,0d0,200d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'pto2',pto(2),wt,0d0,200d0,5d0,'lin')
      n=n+1 
      call bookplot(n,tag,'mjj',mjj,wt,0d0,5000d0, 100d0,'lin')
      n=n+1
      call bookplot(n,tag,'detajj',detajj,wt,0d0,10d0, 1d0,'lin')
      n=n+1
      call bookplot(n,tag,'ptj1',ptj1,wt,0d0,300d0,3d0,'lin')
      n=n+1
      call bookplot(n,tag,'ptj2',ptj2,wt,0d0,300d0,3d0,'lin')
      n=n+1
      call bookplot(n,tag,'category',xjets,wt,-0.5d0,3.5d0,1d0,'lin')
      n=n+1
      call bookplot(n,tag,'mll',m45,wt,0d0,500d0,10d0,'lin')
      n=n+1
      call bookplot(n,tag,'mth',mth,wt,0d0,500d0,10d0,'lin')
      n=n+1
      call bookplot(n,tag,'ptll',pt45,wt,0d0,500d0,10d0,'lin')
      n=n+1
      call bookplot(n,tag,'ptmiss',ptmiss,wt,0d0,500d0,10d0,'lin')
      n=n+1
      if (xjets.ge.0d0) then
        xjets_incl=0d0
        call bookplot(n,tag,'incljets',xjets_incl,wt,-0.5d0,3.5d0,1d0,
     & 'lin') 
      endif
      if (xjets.ge.1d0) then
        xjets_incl=1d0
        call bookplot(n,tag,'incljets',xjets_incl,wt,-0.5d0,3.5d0,1d0,
     & 'lin')
      endif
      if (xjets.ge.2d0) then
        xjets_incl=2d0
        call bookplot(n,tag,'incljets',xjets_incl,wt,-0.5d0,3.5d0,1d0,
     & 'lin')
      endif
      if (xjets.ge.3d0) then  
        xjets_incl=3d0
        call bookplot(n,tag,'incljets',xjets_incl,wt,-0.5d0,3.5d0,1d0,
     & 'lin')
      endif
 


      if (n.gt.20) then
        write(6,*) 'WARNING - TOO MANY HISTOGRAMS!'
        write(6,*) n,' > 20, which is the built-in maximum'
        stop
      endif

c--- set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
      

      return 
      end
      

      
