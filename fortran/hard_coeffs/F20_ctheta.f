      double precision function F20func(ctheta)
         implicit none
         double precision, intent(in) :: ctheta
         double precision F20
         F20=        (8388608*M**2*Pi**6*
     -    (-2*ctheta**2*
     -       (108 + (64*M**6)/Q**6 + 
     -         (272*M**4)/Q**4 - 
     -         (504*M**2)/Q**2) + 
     -      ctheta**4*(756 + (16*M**4)/Q**4)*
     -       (-1 + (4*M**2)/Q**2)**2 + 
     -      324*ctheta**8*
     -       (-1 + (4*M**2)/Q**2)**4 + 
     -      36*ctheta**6*
     -       (-1 + (4*M**2)/Q**2)**3*
     -       (24 + (4*M**2)/Q**2) + 
     -      (16*M**4)/Q**4)*R0Jpsi**4)/
     -  (27.*(11 - (2*nf)/3.)**4*
     -    (1 + ctheta**2*
     -        (-1 + (4*M**2)/Q**2))**4*Q**8*
     -    Log(Q**2/Lambda**2)**4)
         F20func=F20
      end function F20func