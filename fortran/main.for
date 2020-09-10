c TO RUN THE SCRIPT IN ONE GO WITH GFORTRAN, USE THE COMMAND : gfortran -c utils.for && gfortran main.for utils.o && ./a.out

c**********
c This code computes qt-spectra and Q-spectra of S-wave vector quarkonium-pair production from gluon fusion in the TMD framework with TMD evolution. The user can fix the physical constants and kinematical parameters in the first part of the 'main' file.
c It is then possible to choose a definition of the perturbative and nonpertubative Sudakov factors Sa and Snp, as well as the gluon PDF gMSTW and the strong coupling alphas, provided that the definition is present in 'tmd_sudakovs.f' (for Sudakov factors) or in 'tmd_scales.f' (for the other ones).
c The user can then call one of the four subroutines defined in 'main'. One can compute lists of values of the TMD convolutions and their ratios, or ratios of the hard-scattering coefficients times their corresponding convolutions, as functions of qt or Q. The results are saved in a file in the 'ratiodata' folder which name contains information about the parameters used for the computation (cf. definition of the subroutines).
c It is possible to vary the parameter values and the definitions several times in one execution of the main. One simply needs to write the desired configuration and call the subroutines before starting over. A few use examples are provided. Bear in mind that the customisation of the produced data file names is limited, therefore chaining computations while changing parameters not embedded in the title will result in overwriting previous data files with identical names.
c**********

c     Start of the program
      program main

c     Use the 'utils' module to load all necessary helper functions and global-scope variables
      use utils

c     Variable declaration
      implicit none
      double precision qtmin,dqt,qtmax,Qmin,dQ
      integer i,j
      double precision, allocatable :: Qlist(:), qtlist(:)
      character strSud*14,strSa*14,strx1*10,strx2*10,strct*7

      write(6,*) char(10),'---------------------------------------------
     &---------------------',char(10),char(9),'UNPOLARISED DOUBLE J/psi 
     &Sthint WITH EVOLVED TMDs',char(10),'------------------------------
     &------------------------------------',char(10)

c     Location of the pdf data file (here mstw2008lo for leading-order PDFs)
      prefix = "MSTW2008/Grids/mstw2008lo"

c     Physical constants
      Nc = 3d0
      nf = 5d0
      b0 = 1.12292d0
      Q0 = 1d0
      Lambda = 0.217d0
      pi = dacos(-1.d0)
      M=3.0969d0
      R0Jpsi=sqrt(0.81d0)
      phid=0d0
      g0=0.3d0

c     Precision of numerical integral computations
      precision=1d-2

c     List of values for the hard scale Q
      allocate(Qlist(3))
      Qlist=(/8.0,12.0,21.0/)

c     List of values for the transverse momentum qt
      allocate(qtlist(106))
      qtmin=0d0
      dqt=0.2d0
      do i = 1, size(qtlist)
         qtlist(i)=(qtmin+(i-1)*dqt)
      end do

c     Values of gluons x fractions
      x1=1d-2
      x2=1d-2
      write(6,*) 'x1        x2'
      write(6,"(A,1pe8.2,2x,1pe8.2,2x,A)") ' ',x1,x2,char(10)
      xxpdf=x1*x2
      write(strx1,'(1pe8.2,2x)') x1
      write(strx2,'(1pe8.2,2x)') x2

c     Values of the rapidity difference Delta_y of the produced quarkonium pair
c     Large theta_cs (small Delta_y)
      ctmin=0d0   ! cos(theta) min
      ctmax=0.25d0  ! cos(theta) max
      write(6,*) '0<cos(theta)<0.25 => central rapidities'
      write(strct,'(A)') 'central'

c     Small theta_cs (large Delta_y)
c      ctmin=0.25d0
c      ctmax=0.5d0
c      write(6,*) '0.25<cos(theta)<0.5 => forward rapidities'
c      write(strct,'(A)') 'forward'


c----------   qT spectra   ----------

c----------   convol & convol ratios   ----------

c     Select a function for the perturbative Sudakov factor Sa, the PDFs qMSTW and gMSTW and the strong coupling alphas
      write(strSa,"(A)") 'SaO2'
      Sa=>SaO2
      gMSTW=>gMSTWbar
      qMSTW=>qMSTWbar
      alphas=>alphasbar

c     Select a function for the nonperturbative Sudakov factor Snp (as well as the btmax value)
      btmax = 1.5d0

      write(6,*) char(10),'Snp=>Snpbc2 (becomes 0 at bt=2GeV-1 for Q=8Ge
     &V), Sa=>SaO2'
      Snp=>Snpbc2
       write(strSud,"(A)") 'Snpbc2'

c     Call the subroutine to compute convol ratios as functions of qt
      call ComputeRatiosqtlist()

c     Change Snp and repeat the computation
      write(6,*) char(10),'Snp=>Snpbc8 (becomes 0 at bt=8GeV-1 for Q=8Ge
     &V), Sa=>SaO2'
      Snp=>Snpbc8
       write(strSud,"(A)") 'Snpbc8'
      call ComputeRatiosqtlist()

c----------   total ratios (hard scattering coefficients x convols)    ----------

      call TotRatiosqtlist()

c----------   Q spectra   ----------

c----------   convol & convol ratios   ----------

      deallocate(Qlist,qtlist)
      allocate(qtlist(3),Qlist(30))
      qtlist=(/1.0,2.42,4.0/)

      Qmin=8d0
      dQ=0.5d0

      do i = 1, size(Qlist)
         Qlist(i)=(Qmin+(i-1)*dQ)
      end do

      write(6,*) char(10),'Snp=>Snpbc2 (becomes 0 at bt=2GeV-1 for Q=8Ge
     &V), Sa=>SaO2'
      Snp=>Snpbc2
       write(strSud,"(A)") 'Snpbc2'
      call ComputeRatiosQlist()

      write(6,*) char(10),'Snp=>Snpbc8 (becomes 0 at bt=8GeV-1 for Q=8Ge
     &V), Sa=>SaO2'
      Snp=>Snpbc8
       write(strSud,"(A)") 'Snpbc8'
      call ComputeRatiosQlist()

c----------   total ratios (hard scattering coefficients x convols)    ----------

      write(6,*) char(10),'Snp=>Snpbc2 (becomes 0 at bt=2GeV-1 for Q=8Ge
     &V), Sa=>SaO2'
      Snp=>Snpbc2
       write(strSud,"(A)") 'Snpbc2'
      call TotRatiosQlist()

      write(6,*) char(10),'Snp=>Snpbc8 (becomes 0 at bt=8GeV-1 for Q=8Ge
     &V), Sa=>SaO2'
      Snp=>Snpbc8
       write(strSud,"(A)") 'Snpbc8'
      call TotRatiosQlist()

c********** Subroutines **********

c This section contains the subroutines called to compute spectra of single convolutions or convolutions combined with hard-scattering coefficients, as functions of qt or Q

      contains

c Convolutions and their ratios only (no hard-scattering coeficients)

c     Compute convol & convol ratios as functions of qt for fixed values of Q
      subroutine ComputeRatiosqtlist()
      implicit none
      character strQ*10
      double precision Cval(106,5)
c     The routine repeats the computation for each value of Q in Qlist 
      do i = 1, size(Qlist)
         Q = Qlist(i)
c        The maximal considered value of qt is Q/2
         qtmax=Q/2d0
c        Write the results header in the terminal
         write(6,*) char(10),'---------------'
         write(6,"(A,F0.1,A)") ' Q = ',Q,' GeV',char(10)
         write (6,*)'#qt          Cff         Chh         Cw3afh      Cw
     &4hh       Chh/Cff     Cw3afh/Cff  Cw4hh/Cff'
c        Store the value of Q in a string
         write(strQ,'(F0.1)') Q

c        Create (or overwrite) a file containing the computed values of convolutions
c        The file title contains the value of Q and the Sa, Snp functions used for the computations
         open(unit=10,file='ratiodata/convols_Q'//strQ(1:lentrim(strQ))
     &//'_'//strSud(1:lentrim(strSud))//'_'//strSa(1:lentrim(strSa))//'_
     &091018.dat')
c        Create (or overwrite) a file containing the computed ratios of convolutions
         open(unit=11,file='ratiodata/ratios_Q'//strQ(1:lentrim(strQ))//
     &'_'//strSud(1:lentrim(strSud))//'_'//strSa(1:lentrim(strSa))//'_09
     &1018.dat')
c        Write the data headers for both files
         write(10,*) '#Evolved convols at Q='//strQ(1:lentrim(strQ))//'G
     &eV, x1='//strx1//', x2='//strx2//', Snp='//strSud(1:lentrim(strSud
     &))//' and Sa='//strSa(1:lentrim(strSa))
         write(11,*) '#Evolved convols ratios at Q='//strQ(1:lentrim(
     &strQ))//'GeV, x1='//strx1//', x2='//strx2//', Snp='//strSud(1:
     &lentrim(strSud))//' and Sa='//strSa(1:lentrim(strSa))
         write (10,*)'#  qt             Cff            Chh            C
     &w3afh         Cw4hh'
         write (11,*)'#qt          Chh/Cff     Cw3afh/Cff  Cw4hh/Cff'
c        Compute the convolutions and their ratios for each value of qt until reaching qtmax
         do j = 1,size(qtlist)
            qt=qtlist(j)
            if (qt.le.qtmax) then
               Cval(j,1)=Cff()
               Cval(j,2)=Chh()
               Cval(j,3)=Cw3afh()
               Cval(j,4)=Cw4hh()
               write(6,16) qt,Cval(j,1),Cval(j,2),Cval(j,3),Cval(j,4),
     &         Cval(j,2)/Cval(j,1),Cval(j,3)/Cval(j,1),
     &         Cval(j,4)/Cval(j,1)
               write(10,17) qt,Cval(j,1),Cval(j,2),Cval(j,3),Cval(j,4)
               write(11,17) qt,Cval(j,2)/Cval(j,1),Cval(j,3)/Cval(j,1),
     &         Cval(j,4)/Cval(j,1)
            else
               goto 12
            end if
 12      end do
c        Close the files
         close(10)
         close(11)

      end do
c     Formats used to write the data in the files
 16   format(8(1pe10.2,2x))
 17   format(8(1pe15.5,2x))

      end subroutine ComputeRatiosqtlist


c     Compute convol & convol ratios as functions of Q for fixed values of qt
      subroutine ComputeRatiosQlist()
      implicit none
      character strqt*10
      double precision Cval(106,5)

      do i = 1,size(qtlist)
         qt=qtlist(i)
         write(6,*) char(10),'---------------'
         write(6,"(A,F0.1,A)") ' qt = ',qt,' GeV',char(10)
         write (6,*)' Q           Cff         Chh         Cw3afh      Cw
     &4hh       Chh/Cff     Cw3afh/Cff  Cw4hh/Cff'
         write(strqt,'(F0.1)') qt

         open(unit=10,file='ratiodata/convols_qt'//strqt(1:lentrim(strqt
     &))//'_'//strSud(1:lentrim(strSud))//'_'//strSa(1:lentrim(strSa))//
     &'_091018.dat')
         open(unit=11,file='ratiodata/ratios_qt'//strqt(1:lentrim(strqt)
     &)//'_'//strSud(1:lentrim(strSud))//'_'//strSa(1:lentrim(strSa))//'
     &_091018.dat')
         write(10,*) '#Evolved convols at qt='//strqt(1:lentrim(strqt))
     &//'GeV, x1='//strx1//', x2='//strx2//', Snp='//strSud(1:lentrim(
     &strSud))//' and Sa='//strSa(1:lentrim(strSa))
         write(11,*) '#Evolved convols ratios at qt='//strqt(1:lentrim(
     &strqt))//'GeV, x1='//strx1//', x2='//strx2//', Snp='//strSud(1:
     &lentrim(strSud))//' and Sa='//strSa(1:lentrim(strSa))
         write (10,*)'#Q          Cff         Chh         Cw3afh      Cw
     &4hh'
         write (11,*)'#Q          Chh/Cff     Cw3afh/Cff  Cw4hh/Cff'

         do j = 1,size(Qlist)
            Q=Qlist(j)
            Cval(j,1)=Cff()
            Cval(j,2)=Chh()
            Cval(j,3)=Cw3afh()
            Cval(j,4)=Cw4hh()
            write(6,16) Q,Cval(j,1),Cval(j,2),Cval(j,3),Cval(j,4),
     &      Cval(j,2)/Cval(j,1),Cval(j,3)/Cval(j,1),
     &      Cval(j,4)/Cval(j,1)
            write(10,17) Q,Cval(j,1),Cval(j,2),Cval(j,3),Cval(j,4)
            write(11,17) Q,Cval(j,2)/Cval(j,1),Cval(j,3)/Cval(j,1),
     &      Cval(j,4)/Cval(j,1)
 12      end do

         close(10)
         close(11)

      end do

 16   format(8(1pe10.2,2x))
 17   format(8(1pe15.5,2x))

      end subroutine ComputeRatiosQlist

c Ratios including the hard-scattering coefficients

c     Compute total ratios as functions of qt for fixed values of Q
      subroutine TotRatiosqtlist()
      implicit none
      character strQ*10
      double precision Rval(106,3)

      do i = 1, size(Qlist)
         Q = Qlist(i)
         qtmax=Q/2d0
         write(6,*) char(10),'---------------'
         write(6,"(A,F0.1,A)") ' Q = ',Q,' GeV',char(10)
         write (6,*)'#qt          R0          R2          R4'
         write(strQ,'(F0.1)') Q

         open(unit=11,file='ratiodata/totratios_Q'//strQ(1:lentrim(strQ)
     &)//'_'//strSud(1:lentrim(strSud))//'_'//strSa(1:lentrim(strSa))//'
     &_'//strct//'_121118.dat')
         write(11,*) '#Evolved total ratios at Q='//strQ(1:lentrim(
     &strQ))//'GeV, x1='//strx1//', x2='//strx2//', Snp='//strSud(1:
     &lentrim(strSud))//' and Sa='//strSa(1:lentrim(strSa))//' at '//
     &strct//' rapidities'
         write (11,*)'#qt          R0     R2  R4'

         do j = 1,size(qtlist)
            qt=qtlist(j)
            if (qt.le.qtmax) then
               Rval(j,1)=hhtheta_int()/fftheta_int()*Chh()/Cff()
               Rval(j,2)=w3afhtheta_int()/fftheta_int()*Cw3afh()/Cff()
               Rval(j,3)=w4hhtheta_int()/fftheta_int()*Cw4hh()/Cff()
               write(6,16) qt,Rval(j,1),Rval(j,2),Rval(j,3)
               write(11,17) qt,Rval(j,1),Rval(j,2),Rval(j,3)
            else
               goto 12
            end if
 12      end do

         close(11)

      end do

 16   format(8(1pe10.2,2x))
 17   format(8(1pe15.5,2x))

      end subroutine TotRatiosqtlist


c     Compute total ratios as functions of Q for fixed values of qt
      subroutine TotRatiosQlist()
      implicit none
      character strqt*10
      double precision Rval(106,3)

      do i = 1,size(qtlist)
         qt=qtlist(i)
         write(6,*) char(10),'---------------'
         write(6,"(A,F0.1,A)") ' qt = ',qt,' GeV',char(10)
         write (6,*)' Q           R0          R2          R4'
         write(strqt,'(F0.1)') qt

         open(unit=11,file='ratiodata/totratios_qt'//strqt(1:lentrim(
     &strqt))//'_'//strSud(1:lentrim(strSud))//'_'//strSa(1:lentrim(
     &strSa))//'_'//strct//'_121118.dat')
         write(11,*) '#Evolved total ratios at qt='//strqt(1:lentrim(
     &strqt))//'GeV, x1='//strx1//', x2='//strx2//', Snp='//strSud(1:
     &lentrim(strSud))//' and Sa='//strSa(1:lentrim(strSa))//' at '//
     &strct//' rapidities'
         write (11,*)'#Q          R0     R2  R4'
         write(strqt,'(F0.1)') qt

         do j = 1,size(Qlist)
            Q=Qlist(j)
            Rval(j,1)=hhtheta_int()/fftheta_int()*Chh()/Cff()
            Rval(j,2)=w3afhtheta_int()/fftheta_int()*Cw3afh()/Cff()
            Rval(j,3)=w4hhtheta_int()/fftheta_int()*Cw4hh()/Cff()
            write(6,16) Q,Rval(j,1),Rval(j,2),Rval(j,3)
            write(11,17) Q,Rval(j,1),Rval(j,2),Rval(j,3)
 12      end do

         close(11)

      end do

 16   format(8(1pe10.2,2x))
 17   format(8(1pe15.5,2x))

      end subroutine TotRatiosQlist


      end

