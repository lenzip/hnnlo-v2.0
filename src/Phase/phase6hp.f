CC    NEW: generate p12->(p3+p4+p5+p6)+p7+p8
CC    where p3,p4,p5,p6 come from the decay of the Higgs boson

      subroutine phase6hp(r,p1,p2,p3,p4,p5,p6,p7,p8,wt,*)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'zerowidth.f'
      include 'process.f'

      logical oldzerowidth
      integer n2,n3
      double precision r(mxdim)
      double precision p1(4),p2(4),p5(4),p6(4),p3(4),p4(4),p7(4),p8(4)
      double precision p12(4),p3456(4),p34(4),p56(4),p78(4),smin
      double precision wt,wt0,wt12,wt3456,wt78,wt34,wt56


      integer j
      parameter(wt0=1d0/twopi**4)
      wt=0d0
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo

CC    p12->p3456+p78      
   
      call phi1_2m(hmass,r(2),r(3),r(4),1d-10,p12,p3456,p78,wt12,*99)
      
      
CC    p3456->p34+p56
    
      call phi1_2(r(5),r(6),r(7),r(8),p3456,p34,p56,wt3456,*99)

CC    p34->p3+p4

      call phi3m0(r(11),r(12),p34,p3,p4,wt34,*99)

CC    p56->p5+p6

      call phi3m0(r(13),r(14),p56,p5,p6,wt56,*99)

CC    p78->p7+p8

      call phi3m0(r(15),r(16),p78,p7,p8,wt78,*99)

      wt=wt0*wt12*wt3456*wt34*wt56*wt78
      
      return
 99   wt=0d0
      return 1
      end

