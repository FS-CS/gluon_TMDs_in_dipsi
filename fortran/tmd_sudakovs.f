c**********
c This file contains various definitions of the perturbative and nonperturbative Sudakov factors
c**********

      double precision function SaO2(bt)
      implicit none
      double precision mu,bt
      mu = mubar(bt)
c   This version takes into account the O(alpha_s^2) corrections to Sa and uses the simple scale mubar=b0/btbar, with btbar that has both upper and lower cuts
      SaO2 = (12*Nc*(Log(mu**2) - Log(Q**2) + ((-67 + 3*Pi**2 + 20*
     &nf*Tf)*Log(mu**2/Q**2)*Log(Q**2))/(3.*(-33 + 2*nf)*(2*Log(Lambda)
     & - Log(mu**2))*(-2*Log(Lambda) + Log(Q**2))) - Log(Lambda)*Log(1/
     &(-2*Log(Lambda) + Log(mu**2))) + Log(Lambda)*Log(-2*Log(Lambda) +
     & Log(mu**2)) + ((-67 + 3*Pi**2 + 20*nf*Tf)*(-((4*Log(Lambda)**2 +
     & Log(mu**2)*Log(Q**2))*Log((-2*Log(Lambda) + Log(mu**2))/(-2*
     &Log(Lambda) + Log(Q**2)))) + 2*Log(Lambda)*(Log(mu**2)*(-1 + Log(
     &-2*Log(Lambda) + Log(mu**2)) - Log(-2*Log(Lambda) + Log(Q**2))) +
     & Log(Q**2)*(1 + Log(-2*Log(Lambda) + Log(mu**2)) - Log(-2*Log(
     &Lambda) + Log(Q**2))))))/(3.*(-33 + 2*nf)*(2*Log(Lambda) - Log(
     &mu**2))*(-2*Log(Lambda) + Log(Q**2))) - 2*Log(Lambda)*Log(-2*Log(
     &Lambda) + Log(Q**2)) + ((-11 + (2*nf)/Nc)*(-Log(-2*Log(Lambda) + 
     &Log(mu**2)) + Log(-2*Log(Lambda) + Log(Q**2))))/6. + Log(Q**2)*(-
     &Log(-2*Log(Lambda) + Log(mu**2)) + Log(-2*Log(Lambda) + Log(Q**2)
     &))))/(33 - 2*nf)
      end

      double precision function Snpbc1(bt)    ! becomes 0 at bt=1.3GeV-1 for Q=8 GeV
      implicit none
      double precision bt
      Snpbc1 = 2*log(Q)*(bc(bt))**2
      end

      double precision function Snpbc2(bt)    ! becomes 0 at bt=2GeV-1 for Q=8 GeV
      implicit none
      double precision bt
      Snpbc2 = 0.64*log(Q)*(bc(bt))**2
      end

      double precision function Snpbc4(bt)    ! becomes 0 at bt=4GeV-1 for Q=8 GeV
      implicit none
      double precision bt
      Snpbc4 = 0.16*log(Q)*(bc(bt))**2
      end

      double precision function Snpbc8(bt)    ! becomes 0 at bt=8GeV-1 for Q=8 GeV
      implicit none
      double precision bt
      Snpbc8 = 0.04*log(Q)*(bc(bt))**2
      end

c**********

      double precision function Safroz(bt)
      implicit none
      double precision mu,muint,bt
      mu = mubfrozen(bt)
      muint = mub(bt)
c     at bt=0, mub(0)=Infinity => log(log(muint**2/Lambda**2))=Infinity BUT log(Q**2/mubfrozen(0)**2)=log(1)=0 => log(Q**2/mu**2)*log(log(muint**2/Lambda**2)=0 because the log(Q**2/mu**2) goes to 0 faster than the log(log(muint**2/Lambda**2)) goes to Infinity
      Safroz=-1d0/(2d0*nf-33d0)*(6d0*Nc*log(mu**2)+(2d0*nf-11d0*Nc-12d0*
     &Nc*log(Lambda))*log(log(Q**2)-2d0*log(Lambda))-6d0*Nc*log(Lambda)*
     &log(1d0/(log(mu**2)-2d0*log(Lambda)))+6d0*Nc*log(Q**2)*(log(log(
     &Q**2)-2d0*log(Lambda))-log(log(mu**2)-2d0*log(Lambda))-1d0)+(11d0
     &*Nc-2d0*nf+6d0*Nc*log(Lambda))*log(log(mu**2)-2d0*log(Lambda))-
     &6d0*Nc*log(Q**2/mu**2)*log(log(muint**2/Lambda**2)/log(mu**2/
     &Lambda**2)))
      end

      double precision function Sastar(bt)
      implicit none
      double precision mu,muint,bt
      mu = mub(bt)
      muint = mub(bt)
      Sastar=-1d0/(2d0*nf-33d0)*(6d0*Nc*log(mu**2)+(2d0*nf-11d0*Nc-12d0*
     &Nc*log(Lambda))*log(log(Q**2)-2d0*log(Lambda))-6d0*Nc*log(Lambda)*
     &log(1d0/(log(mu**2)-2d0*log(Lambda)))+6d0*Nc*log(Q**2)*(log(log(
     &Q**2)-2d0*log(Lambda))-log(log(mu**2)-2d0*log(Lambda))-1d0)+(11d0
     &*Nc-2d0*nf+6d0*Nc*log(Lambda))*log(log(mu**2)-2d0*log(Lambda))-
     &6d0*Nc*log(Q**2/mu**2)*log(log(muint**2/Lambda**2)/log(mu**2/
     &Lambda**2)))
      end

      double precision function Sacss(bt)
      implicit none
      double precision mu,muint,bt
      mu = mub(bt)
      Sacss=36d0/(2d0*nf-33d0)*(log(Q**2/mu**2)+log(Q**2/Lambda**2)*log(
     &1d0-log(Q**2/mu**2)/log(Q**2/Lambda**2))+(11d0-2d0*nf/Nc)/6d0*log(
     &log(Q**2/Lambda**2)/log(mu**2/Lambda**2)))
      end

      double precision function SnpBLNY(bt)
      implicit none
      double precision bt
c     Ca/Cf=3/(4/3)=9/4 : : the fit of Snp was made for quarks, so we scale it by a factor Ca/Cf to apply it to gluons
      SnpBLNY = 9d0/4d0*(0.68d0*log(Q/(2d0*1.6d0))+0.21d0*(1d0-0.6d0*   ! BLNY (to use with btmax = 0.5 GeV^-1)
     &log(100*xxpdf)))*bt**2
      end


      double precision function SnpKN(bt)
      implicit none
      double precision bt
c     Ca/Cf=3/(4/3)=9/4 : : the fit of Snp was made for quarks, so we scale it by a factor Ca/Cf to apply it to gluons
      SnpKN = 9d0/4d0*(0.184d0*log(Q/(2d0*1.6d0))+0.21d0*(1d0-0.6d0*   ! Konychev-Nadolsky (to use with btmax = 1.5 GeV^-1)
     &log(100*xxpdf)))*bt**2
      end

      double precision function SnpAR(bt)
      implicit none
      double precision bt
c     Ca/Cf=3/(4/3)=9/4 : : the fit of Snp was made for quarks, so we scale it by a factor Ca/Cf to apply it to gluons
c     Aybat-Rogers (to use with btmax = 1.5 GeV^-1)
      SnpAR = 9d0/4d0*(0.184d0*log(Q/(2d0*1.6d0))+0.201d0*(1d0+2d0*
     &(-0.129d0)*log(0.09d0*xxpdf/(0.009d0+xxpdf))))*bt**2
      end

      double precision function SnpCRsimp(bt)
      implicit none
      double precision bt
c     Ca/Cf=3/(4/3)=9/4 : : the fit of Snp was made for quarks, so we scale it by a factor Ca/Cf to apply it to gluons
c     Collins-Roger (to use with btmax = 1.5 GeV^-1) simplified as the btmax-evolution function is 0
c     At btmax = btmax0 = 1.5 GeV^-1 (and g0(btmax0)=0.3 roughly)
      SnpCRsimp = 9d0/4d0*(g0*(1d0-exp(-Nc*alphas(bt)*bt**2/(pi*g0*btmax
     &**2)))*log(Q/(2d0*1.6d0))+0.21d0*(1d0-0.6d0*log(100*xxpdf)))*bt**2
      end

      double precision function SnpARfit(bt)
      implicit none
      double precision bt
      SnpARfit = 0.05*(log(Q/(2d0*1.6d0))+1d0)*bt**2
      end

      double precision function SnpCRsimpfit(bt)
      implicit none
      double precision bt
      SnpCRsimpfit = 9d0/4d0*(0.01*g0*(1d0-exp(-Nc*alphas(bt)*bt**2/(pi
     &*g0*btmax**2)))*log(Q/(2d0*1.6d0))+0.01*0.21d0*(1d0-0.6d0*log(
     &100*xxpdf)))*bt**2
      end

      double precision function SnpCRsimpfit2(bt)
      implicit none
      double precision bt
      SnpCRsimpfit2 = 9d0/4d0*(0.01*g0*(1d0-exp(-Nc*alphas(bt)*bt**2/(pi
     &*g0*btmax**2)))*log(Q/(2d0*1.6d0))+(0.01)*0.21d0*(1d0-0.6d0*log(
     &100*xxpdf)))*bt**2
      end

      double precision function Snp2(bt)
      implicit none
      double precision bt
      Snp2 = (sqrt(bt**2+btmax**2)-btmax)*log(Q**2/Q0**2)
      end

      double precision function Snp2fit(bt)
      implicit none
      double precision bt
      Snp2fit = 0.25*(sqrt(bt**2+btmax**2)-btmax)*log(Q/Q0)
      end

      double precision function Snpfit(bt)
      implicit none
      double precision bt,C1,C2
c     Ca/Cf=3/(4/3)=9/4 : : the fit of Snp was made for quarks, so we scale it by a factor Ca/Cf to apply it to gluons
      Snpfit = 3d0*(C1*log(Q/(2d0*1.6d0))+C2)*bt**2
      end

      double precision function SnpGauss(bt)
      implicit none
      double precision bt
      SnpGauss = bt**2
      end
