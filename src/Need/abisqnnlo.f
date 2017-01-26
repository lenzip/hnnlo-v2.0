       FUNCTION ABINNLO(X)
       IMPLICIT NONE
       double precision X,FR,FI,root,dasin,PI,etap,etam,rlog
       double precision AR,AI
       COMPLEX*16 ABINNLO
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
       ABINNLO=3D0*(AR+(0,1)*AI)
       RETURN
       END


       double precision FUNCTION ABISQNNLO(x)
       implicit none
       complex *16 abinnlo,t1,t2
       double precision x,y
       integer approxim
       external ABINNLO
       t1=abinnlo(x)
       t2=t1*dconjg(t1)
       abisqnnlo=dreal(t2)
       return
       end
