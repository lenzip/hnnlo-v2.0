      logical function includedipole(nd,ptrans)
      implicit none
      include 'constants.f'
      include 'clustering.f'
      include 'npart.f'
      include 'ptilde.f'
      include 'jetlabel.f'
      double precision ptrans(mxpart,4),pjet(mxpart,4),rcut,pt34
      double precision qt2,qtcut
      integer i,j,nd,nqcdjets,nqcdstart,notag,isub
      logical cuts,failedcuts,makecuts,isolation,isol

      common/nqcdjets/nqcdjets,nqcdstart
      common/rcut/rcut
      common/makecuts/makecuts
      common/notag/notag
      common/qtcut/qtcut
      common/isol/isol
      
      integer higgsdec,ndec
      common/higgsdec/higgsdec,ndec

      includedipole=.true.

      if (nd .gt. 0) then
        isub=1
      else
        isub=0
      endif

CC    Compute qt2

      if(higgsdec.eq.1) then
       qt2=(ptrans(3,1)+ptrans(4,1))**2+(ptrans(3,2)+ptrans(4,2))**2      
      else
       qt2=(ptrans(3,1)+ptrans(4,1)+ptrans(5,1)+ptrans(6,1))**2
     &    +(ptrans(3,2)+ptrans(4,2)+ptrans(5,2)+ptrans(6,2))**2
      endif


      if(dsqrt(qt2).lt.qtcut) then
       includedipole=.false.
      else

       call genclust2(ptrans,rcut,pjet,isub)
       do j=1,4
        do i=1,npart+2
        ptildejet(nd,i,j)=pjet(i,j)
        enddo
       enddo


CC    Insert here isolation cut

      if(isol) then
       if(isolation(ptrans).eqv..false.) includedipole=.false.
      endif 
 
c--- check the lepton cuts

        if (makecuts) then
          failedcuts=cuts(pjet,jets)
          if (failedcuts) includedipole=.false.
        endif
 
      endif
      
      return
      end
            
