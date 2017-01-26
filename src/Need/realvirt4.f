CC    Attempt to put real and virtual together for H->3456 at NNLO

      double precision function realvirt4(vector,wgt)
      implicit none
      include 'vegas_common.f'
      integer j,order
      double precision wgt,realint,virtint,lowint,countint
      double precision vector(mxdim),vv(mxdim),tmpr,tmpv,tmpc
      external realint,virtint,lowint,countint
     
      common/nnlo/order

      do j=1,mxdim
      vv(j)=vector(j)
      enddo

CC    Mapping real <-> virtual for H->3456

C     x1,x2

      vv(9)=vector(9)
      vv(10)=vector(10)

C     pH->p34+p56

      vv(1)=vector(5)
      vv(2)=vector(6)
      vv(3)=vector(7)
      vv(4)=vector(8)

C     p34->p3+p4

      vv(7)=vector(11)
      vv(8)=vector(12)

C     p56->p5+p6

      vv(5)=vector(13)
      vv(6)=vector(14)


C     p12->p3456+p7

      vv(13)=vector(1)
      vv(12)=vector(3)
      vv(11)=vector(4)

C     Additional dimension for virtual

      vv(14)=vector(2)

C     Dummies

      vv(15)=vector(15)
      vv(16)=vector(16)



      tmpr=0d0      
      tmpv=0d0


      if(order.eq.2) then
      tmpr=realint(vector,wgt)
      tmpv=virtint(vv,wgt)
      else
      tmpv=lowint(vv,wgt)


      endif 
      
CC    Counterterm

      tmpc=countint(vv,wgt) 

 
      realvirt4=tmpr+tmpv+tmpc

      return
      end
