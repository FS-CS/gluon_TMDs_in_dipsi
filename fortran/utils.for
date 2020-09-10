c**********
c This module combines all the utility/helper functions and global-scope variables needed by 'main.for' to compute TMD spectra with evolution in quarkonium-pair production. 
c**********

      module utils

c********** Abstract interfaces **********

      abstract interface
c     Allows to define a procedure pointer that can refer to different definitions of the scale-dependent functions Sa, Snp or alphas

         double precision function scalepointer(bt)
c        The procedure pointer must have the same input parameters as the procedure it points towards 
            implicit none
            double precision bt
         end function scalepointer

      end interface

      procedure (scalepointer), pointer :: Snp=>null()
      procedure (scalepointer), pointer :: Sa=>null()
      procedure (scalepointer), pointer :: alphas=>null()


      abstract interface

         double precision function PDFscalepointer(xpdf,bt)
c        A new interface is defined for scale-dependent PDFs because they require a different set of arguments
         implicit none
         double precision xpdf,bt
         end function PDFscalepointer

      end interface

      procedure (PDFscalepointer), pointer :: gMSTW=>null()
      procedure (PDFscalepointer), pointer :: qMSTW=>null()

c********** Global-scope variables **********

      double precision x1,x2,xxpdf,Q,precision,ctmin,ctmax
      double precision nf,Nc,btmax,b0,Lambda,Q0,pi,g0,Tf,Cf
      double precision M,R0Jpsi,phid
      double precision qt,vbt,vy1
      character prefix*50


c********** External files containing helper functions to include in the main **********

      contains

      include 'tmd_scales.f'
      include 'tmd_sudakovs.f'
c      include 'tmd_convols.f'
      include 'tmd_convols_qg.f'
      include 'hard_scatter_integrals.f'
      include 'MSTW2008/mstw2008code/mstwpdf_MOD.f'

c     hard-scattering coefficients for S-wave vector quarkonium pair production
      include 'hard_coeffs/F100_ctheta.f'
      include 'hard_coeffs/F20_ctheta.f'
      include 'hard_coeffs/F300a_ctheta.f'
      include 'hard_coeffs/F300b_ctheta.f'
      include 'hard_coeffs/F40_ctheta.f'

c     The 'dcadredo' package is used to compute numerical integrals. Nested integrals require a different version of the package for each integration variable.
      include 'dcadredo/dcadredo.f'
      include 'dcadredo/dcadredo2.f'
      include 'dcadredo/dcadredo3.f'
      include 'dcadredo/dcadredo4.f'


      end module utils

