C     This program computes the Higgs boson production cross section
C     in pp or ppbar collisions at LO, NLO and NNLO accuracy
C
C     Version 1.0 (July 2007)
C
C     Please refer to: 
C     S. Catani, M. Grazzini, Phys. Rev. Lett. 98 (2007) 222002

C     This is the main program

      program hnnlo

      implicit none
      include 'constants.f'
      include 'gridinfo.f'
      integer itmx1,ncall1,itmx2,ncall2
      double precision integ,integ_err
      double precision p(mxpart,4),wt
      common/iterat/itmx1,ncall1,itmx2,ncall2
    

CC    Initialization

      call hinit


CC    Warm up

      if(readin.eqv. .false.)
     & call integrate(0,itmx1,ncall1,.false.,integ,integ_err)

      
CC    Main run    

      call integrate(1,itmx2,ncall2,.true.,integ,integ_err)
    


CC    Final processing and print-out

      call hexit(integ,integ_err)
      
      stop
      end
       
