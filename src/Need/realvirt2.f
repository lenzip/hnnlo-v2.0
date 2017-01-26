CC    Attempt to put real and virtual together: H->2gamma

      double precision function realvirt2(vector,wgt)
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

CC    Mapping real <-> virtual

      vv(1)=vector(2)
      vv(2)=vector(3)
      vv(3)=vector(4)
      vv(4)=vector(7)
      vv(5)=vector(8)
      vv(6)=vector(9)
      vv(7)=vector(10)
CC    Additional dimension for virtual
      vv(8)=vector(1)
CC    Dummies
      vv(9)=vector(5)
      vv(10)=vector(6)

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
 
      realvirt2=tmpr+tmpv+tmpc


      return
      end
