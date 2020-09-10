c**********
c This file contains various definitions of the scales bt and mu used in the computations, as well as alphas (the strong coupling)
c**********

      double precision function bc(bt)
      implicit none
      double precision bt
      bc=sqrt(bt**2+(b0/Q)**2)
      end

      double precision function btbar(bt)   ! btbar = btstar(bc(bt))
      implicit none
      double precision bt
      btbar=sqrt((bt**2+(b0/Q)**2)/(1+(bt/btmax)**2+(b0/(Q*btmax))**2))
      end

      double precision function mubar(bt)
      implicit none
      double precision bt
      mubar = b0/btbar(bt)
      end


      double precision function alphasbar(bt)
c     technically is alphas(mubar(bt)), we take a shortcut and include bt->mubar(bt) inside alphas definition
      implicit none
      double precision bt
      alphasbar = 4d0*pi/((11d0-2d0/3d0*nf)*log(mubar(bt)**2
     &/Lambda**2))
      end

c**********

      double precision function btstar(bt)
      implicit none
      double precision bt
      btstar = bt/(sqrt(1d0+bt**2/btmax**2))
      end

      double precision function mub(bt)
      implicit none
      double precision bt
      mub = b0/btstar(bt)
      end

      double precision function mubfrozen(bt)
      implicit none
      double precision bt
      mubfrozen = Q*b0/sqrt(Q**2*(btstar(bt))**2+b0**2)
      end

      double precision function alphasfroz(bt)
c     technically is alphas(mubfrozen(bt)), we take a shortcut and include bt->mubfrozen(bt) inside alphas definition
      implicit none
      double precision bt
      alphasfroz = 4d0*pi/((11d0-2d0/3d0*nf)*log(mubfrozen(bt)**2
     &/Lambda**2))
      end

      double precision function alphasstar(bt)
      implicit none
      double precision bt
      alphasstar = 4d0*pi/((11d0-2d0/3d0*nf)*log(mub(bt)**2/Lambda**2))
      end
