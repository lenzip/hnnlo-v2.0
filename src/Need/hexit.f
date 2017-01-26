      subroutine hexit(xinteg,xinteg_err)
      implicit none
      double precision xinteg,xinteg_err
      integer itmx1,ncall1,itmx2,ncall2
      common/iterat/itmx1,ncall1,itmx2,ncall2

      
c--- Print-out the value of the integral and its error
      write(6,*) 
      write(6,53)xinteg,xinteg_err
      write(6,*)
     
   53 format('Cross section is',f13.3,' +/-',f10.3,' fb')
 

      call histofin(xinteg,xinteg_err,itmx2,itmx2)
      

      return
      
      end
