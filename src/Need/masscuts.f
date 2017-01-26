      subroutine masscuts(s,*)
      implicit none
      include 'constants.f'
      include 'npart.f'
      include 'limits.f'
      logical first
      double precision s(mxpart,mxpart),scut
      integer nqcdjets,nqcdstart
      common/nqcdjets/nqcdjets,nqcdstart
      integer higgsdec,ndec
      common/higgsdec/higgsdec,ndec
      data first/.true./
      save first

      if(ndec.eq.2) then
       scut=s(3,4)
      else
       scut=s(3,4)+s(3,5)+s(3,6)+s(4,5)+s(4,6)+s(5,6)
      endif


      if (scut.lt.wsqmin.or.scut.gt.wsqmax) return 1

            

     
   98 format(' *      ',f8.2,'  <   ',a12,'  < ',f8.2,'      *')
   99 format(' *          ',f8.2,'  <   ',a3,'  < ',f8.2,'           *')
     
      return
      end

