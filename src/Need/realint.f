      double precision function realint(vector,wgt)
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
      include 'efficiency.f'
      include 'maxwt.f'
      include 'process.f'
      integer ih1,ih2,j,k,nd,nmax,nmin,nvec
      double precision vector(mxdim),W,val,xint
      double precision sqrts,fx1(-nf:nf),fx2(-nf:nf)
      double precision p(mxpart,4),pjet(mxpart,4),p1ext(4),p2ext(4)
      double precision pswt,rscalestart,fscalestart
      double precision s(mxpart,mxpart),wgt,msq(-nf:nf,-nf:nf)
      double precision msqc(maxd,-nf:nf,-nf:nf),xmsq(0:maxd),xmsqjk
      double precision flux,BrnRat,xreal,xreal2,resfac
      double precision xx1,xx2,q(mxpart,4)
      integer n2,n3,sgnj,sgnk
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3
      common/xreal/xreal,xreal2
      logical bin,first,failed
      logical incldip(0:maxd),includedipole,includereal
      external gg_Hg,gg_Hgg,gg_Hg_gs
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/Pext/p1ext,p2ext
      common/nmax/nmax
      common/BrnRat/BrnRat
      common/nmin/nmin
      common/incldip/incldip
      integer higgsdec,ndec
      common/higgsdec/higgsdec,ndec
      common/rescale/resfac


      data p/48*0d0/
      data first/.true./
      save first,rscalestart,fscalestart
      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale
      endif
      ntotshot=ntotshot+1
      pswt=0d0
      realint=0d0      

      W=sqrts**2
      
      if (first) then
         write(6,*)
         write(6,*) 'nmin=',nmin,',nmax=',nmax
         write(6,*)
         first=.false.
      endif
      
      if(higgsdec.eq.1) then
       npart=4
       call gen4(vector,p,pswt,*999)
      else
       npart=6
       call gen6h(vector,p,pswt,*999)
      endif
      
      nvec=npart+2
      call dotem(nvec,p,s)
      
c---impose cuts on final state
      call masscuts(s,*999)


c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)
      

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead set to zero
      includereal=includedipole(0,p)
      incldip(0)=includereal 

      if (includereal .eqv. .false.) then
        do j=-nf,nf
        do k=-nf,nf
          msq(j,k)=0d0
        enddo
        enddo
      endif
      
      
c---- generate collinear points that satisy the jet cuts (for checking)
c      call singgen(p,s,*998)
            
c----calculate the x's for the incoming partons from generated momenta

      xx1=two*(p(1,4)*p2ext(4)-p(1,3)*p2ext(3))/W
      xx2=two*(p(2,4)*p1ext(4)-p(2,3)*p1ext(3))/W

      
      if ((xx1 .gt. 1d0) .or. (xx2 .gt. 1d0)) then
         realint=0d0
         return
      endif

c--- Calculate the required matrix elements    

      if(higgsdec.eq.1) then
       if (includereal) call gg_hgg(p,msq)
       call gg_hg_gs(p,msqc)
      elseif(higgsdec.eq.2) then
       if (includereal) call gg_hwwgg(p,msq)
       call gg_hwwg_gs(p,msqc)
      else
       if (includereal) call gg_hzzgg(p,msq)
       call gg_hzzg_gs(p,msqc)
      endif

      msq=msq*resfac
      msqc=msqc*resfac


      do nd=0,ndmax
      xmsq(nd)=0d0
      enddo
      
            
      flux=fbGeV2/(two*xx1*xx2*W)


  777 continue    
      do nd=0,ndmax
      xmsq(nd)=0d0
      enddo

            
c--- calculate PDF's  
      call fdist(ih1,xx1,facscale,fx1)
      call fdist(ih2,xx2,facscale,fx2)
      
      do j=-nf,nf
      do k=-nf,nf

      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) goto 20
      endif

      if (gqonly) then
      if (((j.eq.0).and.(k.eq.0)) .or. ((j.ne.0).and.(k.ne.0))) goto 20
      endif      
      
      if (noglue) then 
      if ((j.eq.0) .or. (k.eq.0)) goto 20
      endif

      if (realonly) then 
        xmsq(0)=xmsq(0)+fx1(j)*fx2(k)*msq(j,k)
        do nd=1,ndmax
        xmsq(nd)=0d0
        enddo
      elseif (virtonly) then
         xmsq(0)=0d0
         do nd=1,ndmax
           xmsq(nd)=xmsq(nd)+fx1(j)*fx2(k)*(-msqc(nd,j,k))
         enddo
      else

         if     (j .gt. 0) then
           sgnj=+1
         elseif (j .lt. 0) then
           sgnj=-1
         else
           sgnj=0
         endif
         if     (k .gt. 0) then
           sgnk=+1
         elseif (k .lt. 0) then
           sgnk=-1
         else
           sgnk=0
         endif

         xmsqjk=fx1(j)*fx2(k)*msq(j,k)
         xmsq(0)=xmsq(0)+xmsqjk
       
         do nd=1,ndmax
           xmsqjk=fx1(j)*fx2(k)*(-msqc(nd,j,k))
           xmsq(nd)=xmsq(nd)+xmsqjk
         enddo
         
      endif
 20   continue

      enddo
      enddo


      realint=0d0
      xint=0d0

c---trial with weight of real alone
c---first set up all dipole contributions
c---this is the value of integral including subtractions
      do nd=0,ndmax
        xmsq(nd)=xmsq(nd)*flux*pswt/BrnRat
        failed=.false.
        
c--- if this dipole has no contribution, go to end of loop
c        if (xmsq(nd) .eq. 0d0) goto 997         
         
        if (nd .eq. 0) then
c---if there's no real contribution, record the event as failing to pass cuts
          if (xmsq(nd) .eq. 0d0) then
             failed=.true.
             goto 996
          endif
        else
c--- if this dipole has no contribution, go to end of loop
          if (xmsq(nd) .eq. 0d0) goto 997         
c---check whether each counter-event passes the cuts
          do j=1,mxpart
          do k=1,4
          q(j,k)=ptilde(nd,j,k)
          enddo
          enddo
          incldip(nd)=includedipole(nd,q)
          if (incldip(nd) .eqv. .false.) failed=.true.
        endif

 996    if (failed) then
          if (nd .eq. 0) then
            ncutzero=ncutzero+1
            ntotzero=ntotzero+1
          endif
          call dotem(nvec,p,s)
          xmsq(nd)=0d0
          goto 997         
        endif
c---if it does, add to total
        xint=xint+xmsq(nd)

        val=xmsq(nd)*wgt
                
c--- update the maximum weight so far, if necessary
        if (dabs(val) .gt. wtmax) then
          wtmax=dabs(val)
        endif

c---if we're binning, add to histo too
        if (bin) then
          call getptildejet(nd,pjet)
          call dotem(nvec,pjet,s)
          val=val/dfloat(itmx)
          if (nd .eq. 0) then
            call plotter(pjet,val,0)
          else
            call plotter(pjet,val,1)
          endif
        endif
c---otherwise, skip contribution
 997    continue
      enddo

      call dotem(nvec,p,s)

c 998  continue


      realint=xint

      xreal=xreal+xint*wgt/dfloat(itmx)
      xreal2=xreal2+(xint*wgt)**2/dfloat(itmx)
      

      return

 999  realint=0d0
      ntotzero=ntotzero+1
 
      return
      end
















