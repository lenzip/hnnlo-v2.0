      subroutine higgsw(br)
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'
      double precision wff,mfsq,br
      double precision mtop,mbot
      common/tbmass/mtop,mbot
c-----approximate form for the width of the standard model higgs
c-----valid for low masses
      wff(mfsq)=dsqrt(2d0)/8d0/pi*gf*hmass*mfsq
     & *(1d0-4d0*mfsq/hmass**2)**1.5d0

      hwidth=3d0*(wff(mbot**2)+wff(mcsq))+wff(mtausq)
      write(6,*) 'hmass,hwidth',hmass,hwidth
      write(6,*) 'mtausq,mcsq,mbot**2',mtausq,mcsq,mbot**2
      write(6,*) 
      br=3d0*wff(mbot**2)/hwidth
      return
      end
