      double precision function F300afunc(ctheta)
         implicit none
         double precision, intent(in) :: ctheta
         double precision F300a
         F300a=        (-4194304*(-1 + ctheta**2)*Pi**6*
     -    (-6*ctheta**2*
     -       (264 + (64*M**6)/Q**6 + 
     -         (416*M**4)/Q**4 - 
     -         (1164*M**2)/Q**2) + 
     -      1944*ctheta**8*
     -       (-1 + (4*M**2)/Q**2)**4 + 
     -      18*ctheta**6*
     -       (-1 + (4*M**2)/Q**2)**3*
     -       (304 + (36*M**2)/Q**2) + 
     -      ctheta**4*(-1 + (4*M**2)/Q**2)**2*
     -       (5112 + (48*M**4)/Q**4 + 
     -         (64*M**2)/Q**2) + 
     -      (48*M**4)/Q**4 - (64*M**2)/Q**2)*
     -    (-1 + (4*M**2)/Q**2)*R0Jpsi**4*
     -    Cos(2*phid))/
     -  (81.*(11 - (2*nf)/3.)**4*
     -    (1 + ctheta**2*
     -        (-1 + (4*M**2)/Q**2))**4*Q**6*
     -    Log(Q**2/Lambda**2)**4)
         F300afunc=F300a
      end function F300afunc