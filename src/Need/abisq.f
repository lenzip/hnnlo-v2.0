       FUNCTION ABI(X)
       IMPLICIT NONE
       double precision X,FR,FI,root,dasin,PI,etap,etam,rlog
       double precision AR,AI
       COMPLEX*16 ABI
       PI=4.D0*DATAN(1.D0)
       IF(X.GE.0.25D0) THEN
                 FR=-2.D0*(DASIN(0.5D0/DSQRT(X)))**2
                 FI=0.D0
       ELSEIF(X.LT.0.25D0) THEN
                 ROOT=DSQRT(0.25D0-X)
                 ETAP=0.5D0+ROOT
                 ETAM=0.5D0-ROOT
                 RLOG=DLOG(ETAP/ETAM)
                 FR=0.5D0*(RLOG**2-PI**2)
                 FI=PI*RLOG
                 ENDIF
       AR=2.D0*X+X*(4.D0*X-1.D0)*FR
       AI=       X*(4.D0*X-1.D0)*FI
       ABI=3D0*(AR+(0,1)*AI)
       RETURN
       END


       double precision FUNCTION ABISQ(x,y)
       implicit none
       complex *16 abi,t1,t2
       double precision x,y
       integer approxim
       common/flag/approxim
       external ABI
       t1=abi(x)+(approxim-1)*abi(y)
ch       t1=abi(x)
       t2=t1*dconjg(t1)
       abisq=dreal(t2)
       return
       end
