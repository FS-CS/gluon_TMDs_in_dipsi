c     *********   cos(theta)-integral of hard-scattering coefficients  *********

      double precision function fftheta_int()
         implicit none
         real*8 err4
         integer ierr4

         fftheta_int=dcadredo4(F100func,ctmin,ctmax,0.d0
     &   ,precision/10.,err4,ierr4)

      end


c*********************************************

      double precision function hhtheta_int()
         implicit none
         real*8 err4
         integer ierr4

         hhtheta_int=dcadredo4(F20func,ctmin,ctmax,0.d0
     &   ,precision/10.,err4,ierr4)

      end


c*********************************************

      double precision function w3afhtheta_int()
         implicit none
         real*8 err4
         integer ierr4

         w3afhtheta_int=dcadredo4(F300afunc,ctmin,ctmax,0.d0
     &   ,precision/10.,err4,ierr4)

      end


c*********************************************

      double precision function w3bhftheta_int()
         implicit none
         real*8 err4
         integer ierr4

         w3bhftheta_int=dcadredo4(F300bfunc,ctmin,ctmax,0.d0
     &   ,precision/10.,err4,ierr4)

      end




c*********************************************

      double precision function w4hhtheta_int()
         implicit none
         real*8 err4
         integer ierr4

         w4hhtheta_int=dcadredo4(F40func,ctmin,ctmax,0.d0
     &   ,precision/10.,err4,ierr4)

      end
