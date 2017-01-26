c.....Function BK(n,z)
c.....BK(n,z) is the n-derivative of BesselK[nu,z]
c.....with respect to nu in nu=1

c.....Itilde defined as in the paper

      function Itilde(m)
      implicit none
      double precision zBK,argum,Itilde
      double precision Eulergamma,b0,z2,z3,logx
      double precision xmio 
      integer m
      common/xmio/xmio

      external zBK

      Eulergamma=0.577215664902d0
      z2=1.64493406685d0
      z3=1.20205690316d0
      b0=2*dexp(-Eulergamma)

C

      argum=b0*xmio
      logx=dlog(xmio)
      
      if (m.eq.1) then
         Itilde=-zBK(0,argum)/xmio**2
      elseif (m.eq.2) then
         Itilde=2d0/xmio**2*(zBK(0,argum)*logx-zBK(1,argum))
      elseif (m.eq.3) then
         Itilde=-3/xmio**2*(zBK(0,argum)*(logx**2-z2)
     &          -2*zBK(1,argum)*logx+zBK(2,argum))
      elseif (m.eq.4) then
         Itilde=4/xmio**2*(zBK(0,argum)*(logx**3-3*z2*logx+2*z3)
     &         -3*zBK(1,argum)*(logx**2-z2)+3*zBK(2,argum)*logx
     &         -zBK(3,argum))
      endif


      return
      end

C     n-derivative of the function BesselK[nu,z] 
C     with respect to nu for nu=1
C     NOTE: IT IS MULTIPLIED by z

      function zBK(n,z)
      implicit none
      double precision zbk,fb,errest,z,zz,max,adpint
      double precision zbesselk0,zbesselk1,zbesselk2,zbesselk3
      integer n,nn,ifail
      common/nuorder/nn
      common/zz/zz
      external fb,zbesselk0,zbesselk1,zbesselk2,zbesselk3
      nn=n
      zz=z
      max=10d0

C     Use approximated form only for z<1.5d0

      if(z.gt.1.5d0) then
      zbk=z*adpint(fb,0d0,max,1d-10,1d-5,errest,ifail)
      goto 99
      endif

      if(n.eq.0) then
      zbk=zbesselk0(z)
      elseif(n.eq.1) then
      zbk=zbesselk1(z)
      elseif(n.eq.2) then
      zbk=zbesselk2(z)      
      elseif(n.eq.3) then
      zbk=zbesselk3(z)
      endif

 99   return
      end

C     Approximated forms of BesselK[nu,z] for nu=1 and
C     derivatives with respect to nu of BesselK[nu,z] at nu=1
C     All functions multiplied by z

      function zbesselk0(z)
      implicit none
      double precision zbesselk0,z,zm,loz,egamma
      egamma=0.577215664902d0
      zm=z/2d0
      loz=dlog(zm)

      zbesselk0=1d0+z*zm*(loz-0.5d0*(1-2*egamma))
     &         +zm**4*(loz-0.5d0*(2.5d0-2*egamma))
     &         +zm**6/6d0*(loz-0.5d0*(10d0/3-2*egamma))
     &         +zm**8/72d0*(loz-0.5d0*(47d0/12-2*egamma))
     &         +zm**10/1440d0*(loz-0.5d0*(131d0/30-2*egamma))
      
      return
      end


      function zbesselk1(z)
      implicit none
      double precision zbesselk1,z,zm,loz,egamma
      egamma=0.577215664902d0
      zm=z/2d0
      loz=dlog(zm)
      zbesselk1=-(loz+egamma)-zm**2*(loz-1+egamma)
     &         -0.25d0*zm**4*(loz-1.5d0+egamma)
     &         -zm**6/36d0*(loz-11d0/6+egamma)
     &         -zm**8/576d0*(loz-25d0/12+egamma)
      

      return
      end


      function zbesselk2(z)
      implicit none
      double precision zbesselk2,z,a(0:13),loz,zm
      data a(0) / 1.15443132980306572d0/
      data a(1) / 1.97811199065594511d0/
      data a(2) / 0.154431329803065721d0/
      data a(3) / 4.801792651508824500d0/
      data a(4) / 0.806235643470665767d0/
      data a(5) /-0.672784335098467139d0/
      data a(6) / 3.285072828402112960d0/
      data a(7) /-1.945338757678943440d0/
      data a(8) /-0.181575166960855634d0/
      data a(9) / 0.694195147571435559d0/
      data a(10)/-0.607655744858515573d0/
      data a(11)/-0.019182189839330562d0/
      data a(12)/ 0.068894530444636532d0/
      data a(13)/-0.070514317816328185d0/


      zm=z/2
      loz=dlog(zm)
      
      zbesselk2=loz**2+a(0)*loz+a(1)
     &       +zm**2*(2*loz**3/3d0+a(2)*loz**2+a(3)*loz+a(4))
     &       +zm**4*(loz**3/3d0+a(5)*loz**2+a(6)*loz+a(7))
     &       +zm**6*(loz**3/18d0+a(8)*loz**2+a(9)*loz+a(10))
     &       +zm**8*(loz**3/216d0+a(11)*loz**2+a(12)*loz+a(13))
      return
      end


      function zbesselk3(z)
      implicit none
      double precision zbesselk3,z,b(0:14),loz,zm

      data b(0) / 1.731646994704598580d0/
      data b(1) / 5.934335971967835330d0/
      data b(2) / 5.444874456485317730d0/
      data b(3) /-1.268353005295401420d0/ 
      data b(4) / 8.471041982558638170d0/
      data b(5) /-3.026167526073320430d0/
      data b(6) /-0.692088251323850355d0/ 
      data b(7) / 2.809848746963509900d0/
      data b(8) /-2.161466255000085060d0/
      data b(9) /-0.104676472369316706d0/
      data b(10)/ 0.381989731242156681d0/
      data b(11)/-0.367492827636283900d0/
      data b(12)/-0.007844362856415627d0/
      data b(13)/ 0.027796539630842606d0/
      data b(14)/-0.029917436634978395d0/


      zm=z/2
      loz=dlog(zm)
      
      zbesselk3=loz**3+b(0)*loz**2+b(1)*loz+b(2)
     &       +zm**2*(loz**3+b(3)*loz**2+b(4)*loz+b(5))
     &       +zm**4*(loz**3/4d0+b(6)*loz**2+b(7)*loz+b(8))
     &       +zm**6*(loz**3/36d0+b(9)*loz**2+b(10)*loz+b(11))
     &       +zm**8*(loz**3/576d0+b(12)*loz**2+b(13)*loz+b(14))
      zbesselk3=-zbesselk3

      return
      end

      
      function fb(t)
      implicit none
      integer nn,nu
      double precision fb,t,zz
      common/nuorder/nn
      common/zz/zz
      nu=1
      if(nn.eq.0) then
         fb=dexp(-zz*dcosh(t))*dcosh(nu*t)
      elseif(nn.eq.1) then
         fb=dexp(-zz*dcosh(t))*t*dsinh(nu*t)
      elseif(nn.eq.2) then
         fb=dexp(-zz*dcosh(t))*t*t*dcosh(nu*t)
      elseif(nn.eq.3) then
         fb=dexp(-zz*dcosh(t))*t*t*t*dsinh(nu*t)
      endif      
      return
      end
      
