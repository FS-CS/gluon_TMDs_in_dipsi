      double precision function F40func(ctheta)
         implicit none
         double precision, intent(in) :: ctheta
         double precision F40
         F40=        (524288*(-1 + ctheta**2)**2*Pi**6*
     -    (256 + 3888*ctheta**8*
     -       (-1 + (4*M**2)/Q**2)**4 + 
     -      72*ctheta**6*
     -       (-1 + (4*M**2)/Q**2)**3*
     -       (160 + (12*M**2)/Q**2) + 
     -      ctheta**4*(-1 + (4*M**2)/Q**2)**2*
     -       (11632 + (48*M**4)/Q**4 + 
     -         (128*M**2)/Q**2) - 
     -      2*ctheta**2*
     -       (-2128 + 
     -         (12*M**2*(36 + (4*M**2)/Q**2))/
     -          Q**2)*(-1 + (2*M)/Q)*
     -       (1 + (2*M)/Q) + (48*M**4)/Q**4 - 
     -      (128*M**2)/Q**2)*
     -    (-1 + (4*M**2)/Q**2)**2*R0Jpsi**4*
     -    Cos(4*phid))/
     -  (81.*M**2*(11 - (2*nf)/3.)**4*
     -    (1 + ctheta**2*
     -        (-1 + (4*M**2)/Q**2))**4*Q**4*
     -    Log(Q**2/Lambda**2)**4)
         F40func=F40
      end function F40func