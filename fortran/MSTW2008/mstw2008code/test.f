      PROGRAM example

      IMPLICIT NONE
      DOUBLE PRECISION x,q,xf,GetOnePDF
      CHARACTER prefix*50

      prefix = "../Grids/mstw2008lo"

      x=1.d-2
      q=12.d0

      xf = GetOnePDF(prefix,0,x,q,0) ! central value of gluon momentum density


      WRITE(6,*) 'x=',x,'   q=',q,'GeV   xf=',xf,'   f=',xf/x

      END
