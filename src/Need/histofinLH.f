      block data linlog_data
      implicit none
      include 'nplot.f'
      data linlog/150*'lin'/
      end

      subroutine histofin(xsec,xsec_err,itno,itmx)
c--- This outputs the final histograms for itno=0
c--- For itno>0, this is an intermediate result only
      implicit none
      include 'nplot.f'
      include 'histo.f'
      integer j,nlength,itno,itmx,nplotmax,L
      logical fin,snd
      common/fin/fin
      character*30 runstring
      character*72 runname,outfiledat,outfiletop,outfileerr
      character*3 oldbook
      double precision xsec,xsec_err
      double precision EHIST(4,40,100)   
      integer IHISTOMATCH(100),ICOUNTHISTO 
      common/runstring/runstring      
      common/runname/runname
      common/nlength/nlength
      common/nplotmax/nplotmax
      COMMON/EHISTO/EHIST,IHISTOMATCH,ICOUNTHISTO
      character*4 part,mypart
      common/mypart/mypart
      common/part/part
      integer order
      common/nnlo/order

C     Part is what it is computing
C     mypart is what is set at the beginning

C     Fin is true at the very last iteration
C     Snd is true when computing real after virt

      snd=.false.

      if((part.eq.'real').and.(mypart.eq.'tota')) snd=.true.


c
c     Accumula i valori e i valori al quadrato per l'analisi statistica,
c     e svuota l'istogramma di accumulo.
c

      do j=1,nplotmax
       call mopera(j,'A',j+20,j+40,1d0,1d0)
      enddo

    
      if(itno.lt.itmx) then

c      do j=1,nplotmax

c      call flush(6)

c--- ensure that MFINAL doesn't turn off booking for intermediate results
c      oldbook=book(j)
c      call mfinal(j)
c      if (itno .gt. 0) then
c      book(j)=oldbook
c      endif
c      enddo


C     Generate intermediate topdrawer file (no errors !)

      do j=1,nplotmax
      call flush(6)
      
      call strcat(runstring,'.top',outfiletop)
      open(unit=99,file=outfiletop,status='unknown')

CC    Rescale by itmx/itno and put in j+120

      call mopera(j+20,'F',2,j+120,dfloat(itmx)/dfloat(itno),1d0)

CC    If needed, combine with final result for virtual

      if(snd) then
       call mopera(j+120,'+',j+80,j+120,1d0,1d0)
      endif   

      call mfinal(j+120)
      call mtop(j+120,150,'x','y',linlog(j))
      
      enddo
      close (unit=99)


CC If we are at the last step

      else

      call strcat(runstring,'.top',outfiletop)
      
      open(unit=99,file=outfiletop,status='unknown')
      
c--- write out run info to top of files

      call writeinfo(99,xsec,xsec_err,itno)      


C     Complete statistical analysis


      do j=1,nplotmax

c accumula j+1 riscalato in j e pone in j+2 la stima dell'errore
c

      call mopera(j+20,'E',j+40,j,1d0,1d0)


c
c accumula l'errore in quadratura
c
      call mopera(j+40,'Q',j+60,j+60,1d0,1d0)

      enddo



C If this is the end rescale by itmx and put results,error in j,j+60

C If not rescale by itmx and put results,error in j+80,j+100
C Then clear j+20,j+40,j+60,j+80 and number of entries


      do j=1,nplotmax

      if(fin) then
       call mopera(j,'F',2,j,dfloat(itmx),1d0)
       call mopera(j+60,'F',2,j+60,dfloat(itmx),1d0)
      else
       call mopera(j,'F',2,j+80,dfloat(itmx),1d0)
       call mopera(j+60,'F',2,j+100,dfloat(itmx),1d0)
       call mopera(j,'F',2,j,0d0,1d0)          ! clear
       call mopera(j+20,'F',2,j+20,0d0,1d0)    ! clear
       call mopera(j+40,'F',2,j+20,0d0,1d0)    ! clear
       call mopera(j+60,'F',2,j+60,0d0,1d0)    ! clear
       ihis(j,1)=0
       ihis(j+20,1)=0
       ihis(j+40,1)=0
       ihis(j+60,1)=0
      endif     
 
      enddo


C If needed combine the histograms (used when computing real+virt)

      do j=1,nplotmax

      if(fin.and.(mypart.eq.'tota').and.(order.ne.0)) then
       call mopera(j,'+',j+80,j,1d0,1d0) ! sum real+virtual
C Combine errors in quadrature
       call mopera(j+60,'S',j,j+60,1d0,1d0) ! square it
       call mopera(j+100,'S',j,j+100,1d0,1d0)! square it
       call mopera(j+60,'+',j+100,j+60,1d0,1d0) ! sum squares
       call mopera(j+60,'R',j,j+60,1d0,1d0) ! takes square root 
      endif

      enddo

C     Generate final topdrawer file

      do j=1,nplotmax
      call flush(6)  

      if(fin) then    
       call mfinal(j)
       call mtop(j,60+j,'x','y',linlog(j))
      else
       call mfinal(j+80)
       call mtop(j+80,100+j,'x','y',linlog(j))      
      endif

      enddo

      endif

      close (unit=99)


      return
      end

