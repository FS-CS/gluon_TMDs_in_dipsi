c**********
c This file contains the different TMD convolutions as functions of bt or integrated over it, as well as definitions of the gluon PDF using the MSTW2008 data. The convolutions are nested integrals computed successively: the fully unintegrated convolution is integrated using the 'dcadredo' package and the result is passed to the next integration function using a global-scope variable.
c**********

c*****************   gMSTW   *****************

      double precision function gMSTWbar(xpdf,bt)
      implicit none
      double precision xpdf,bt
      gMSTWbar = GetOnePDF(prefix,0,xpdf,mubar(bt),0)/xpdf
      end

      double precision function gMSTWfroz(xpdf,bt)
      implicit none
      double precision xpdf,bt
      gMSTWfroz = GetOnePDF(prefix,0,xpdf,mubfrozen(bt),0)/xpdf
      end

      double precision function gMSTWstar(xpdf,bt)
      implicit none
      double precision xpdf,bt
      gMSTWstar = GetOnePDF(prefix,0,xpdf,mub(bt),0)/xpdf
      end

c*****************   Cff   *****************

      double precision function Cff()
      implicit none
      double precision btlow,bthigh
      real*8 err1
      integer ierr1
c     bt-integration 
      btlow = 0.00001
      bthigh = 30d0
      Cff = dcadredo(Cff_bt,btlow,bthigh,0.d0,precision/10.,err1,ierr1)
      end

c Remark: Sastar(0)=Infinity => error BUT exp(-Infinity)=0, so very small bt doesn't contribute to the integral value, hence the very small but non-zero btlow
c Safroz(0)=0 => exp(0)=1 but starting at bt=0 or bt sligthly greater than 0 doesn't change results (tested)
 
      double precision function Cff_bt(bt)
      implicit none
      double precision bt
      Cff_bt = 1d0/(2d0*pi)*bt*bessel_j0(bt*qt)*exp(-Sa(bt)-
     &Snp(bt))*gMSTW(x1,bt)*gMSTW(x2,bt)
      end

c*****************   Chh   *****************

      double precision function Chh()
      implicit none
      double precision btlow,bthigh
      real*8 err3
      integer ierr3
c     bt-integration 
      btlow = 0.00001
      bthigh = 30d0
      Chh = dcadredo3(Chh_bt,btlow,bthigh,0.d0,precision/10.,err3,ierr3)
      end

      double precision function Chh_bt(bt)
      implicit none
      double precision y1min,y1max,bt
      real*8 err2
      integer ierr2
      vbt=bt
c     y1-integration 
      y1min = x1
      y1max = 1d0
      Chh_bt = dcadredo2(Chh_bt_y1,y1min,y1max,0.d0,precision/10.,
     &err2,ierr2)
      end

      double precision function Chh_bt_y1(y1)
      implicit none
      double precision y2min,y2max,y1
      real*8 err1
      integer ierr1
      vy1=y1
c     y2-integration 
      y2min = x2
      y2max = 1d0
      Chh_bt_y1 = dcadredo(Chh_bt_y1_y2,y2min,y2max,0.d0,precision/10.,
     &err1,ierr1)
      end

      double precision function Chh_bt_y1_y2(y2)
      implicit none
      double precision y2
      Chh_bt_y1_y2 = 1d0/(2d0*pi)*vbt*bessel_j0(vbt*qt)*exp(-Sa(vbt)-
     &Snp(vbt))*(Nc*alphas(vbt)/pi)**2*(1d0/x1-1d0/vy1)*
     &gMSTW(vy1,vbt)*(1d0/x2-1d0/y2)*gMSTW(y2,vbt)
      end

c*****************   Cw3afh   *****************

      double precision function Cw3afh()
      implicit none
      double precision btlow,bthigh
      real*8 err2
      integer ierr2
c     bt-integration 
      btlow = 0.00001
      bthigh = 30d0
      Cw3afh = dcadredo2(Cw3afh_bt,btlow,bthigh,0.d0,precision/10.,
     &err2,ierr2)
      end

      double precision function Cw3afh_bt(bt)
      implicit none
      double precision y1min,y1max,bt
      real*8 err1
      integer ierr1
      vbt=bt
c     y1-integration 
      y1min = x1
      y1max = 1d0
      Cw3afh_bt = dcadredo(Cw3afh_bt_y1,y1min,y1max,0.d0,precision/10.,
     &err1,ierr1)
      end

      double precision function Cw3afh_bt_y1(y1)
      implicit none
      double precision y1
      Cw3afh_bt_y1=1d0/(2d0*pi)*vbt*bessel_jn(2,vbt*qt)*exp(-Sa(vbt)-
     &Snp(vbt))*(Nc*alphas(vbt)/pi)*(1d0/x1-1d0/y1)*gMSTW(y1,vbt)
     &*gMSTW(x2,vbt)
      end

c*****************   Cw3bhf   *****************

      double precision function Cw3bhf()
      implicit none
      double precision btlow,bthigh
      real*8 err2
      integer ierr2
c     bt-integration 
      btlow = 0.00001
      bthigh = 30d0
      Cw3bhf = dcadredo2(Cw3bhf_bt,btlow,bthigh,0.d0,precision/10.,
     &err2,ierr2)
      end

      double precision function Cw3bhf_bt(bt)
      implicit none
      double precision y2min,y2max,bt
      real*8 err1
      integer ierr1
      vbt=bt
c     y2-integration 
      y2min = x2
      y2max = 1d0
      Cw3bhf_bt = dcadredo(Cw3bhf_bt_y2,y2min,y2max,0.d0,precision/10.,
     &err1,ierr1)
      end

      double precision function Cw3bhf_bt_y2(y2)
      implicit none
      double precision y2
      Cw3bhf_bt_y2=1d0/(2d0*pi)*vbt*bessel_jn(2,vbt*qt)*exp(-Sa(vbt)-
     &Snp(vbt))*(Nc*alphas(vbt)/pi)*(1d0/x2-1d0/y2)*gMSTW(y2,vbt)
     &*gMSTW(x1,vbt)
      end

c*****************   Cw4hh   *****************

      double precision function Cw4hh()
      implicit none
      double precision btlow,bthigh
      real*8 err3
      integer ierr3
c     bt-integration 
      btlow = 0.00001
      bthigh = 30d0
      Cw4hh = dcadredo3(Cw4hh_bt,btlow,bthigh,0.d0,precision/10.,err3,
     &ierr3)
      end

      double precision function Cw4hh_bt(bt)
      implicit none
      double precision y1min,y1max,bt
      real*8 err2
      integer ierr2
      vbt=bt
c     y1-integration 
      y1min = x1
      y1max = 1d0
      Cw4hh_bt = dcadredo2(Cw4hh_bt_y1,y1min,y1max,0.d0,precision/10.,
     &err2,ierr2)
      end

      double precision function Cw4hh_bt_y1(y1)
      implicit none
      double precision y2min,y2max,y1
      real*8 err1
      integer ierr1
      vy1=y1
c     y2-integration 
      y2min = x2
      y2max = 1d0
      Cw4hh_bt_y1 = dcadredo(Cw4hh_bt_y1_y2,y2min,y2max,0.d0,
     &precision/10.,err1,ierr1)
      end

      double precision function Cw4hh_bt_y1_y2(y2)
      implicit none
      double precision y2
      Cw4hh_bt_y1_y2 = 1d0/(2d0*pi)*vbt*bessel_jn(4,vbt*qt)*exp(-
     &Sa(vbt)-Snp(vbt))*(Nc*alphas(vbt)/pi)**2*(1d0/x1-1d0/vy1)*
     &gMSTW(vy1,vbt)*(1d0/x2-1d0/y2)*gMSTW(y2,vbt)
      end
