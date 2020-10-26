

      subroutine mltdg (iout, iprint, v, lenv, norbs, x, alpang, betang,
     $     gamang, mode)

      implicit none


C +------------------------------------------------------------+
C |  mltdg  --   CAJH, 07.2011                                 |
C |                                                            |
C |                                     (based on PHFB mltdg)  |
C |                                                            |
C |                                                            |
C |  Perform matrix multiplications involving the Dg matrix    |
C |  of the spin projection operator in an ORTHONORMAL basis   |
C |  (the AO orthonormal basis).                               |
C |                                                            |
C |  The matrix Dg is given by the product                     |
C |                                                            |
C |    Dg = D_a . D_b . D_g                                    |
C |                                                            |
C |  where a == alpha, b == beta, g == gamma. Each of these    |
C |  matrices can be written in spin blocks as                 |
C |                                                            |
C |    D_a = ( exp(-i*alpha/2) . I         0           ),      |
C |          (        0             exp(i*alpha/2) . I )       |
C |                                                            |
C |    D_b = (   cos(beta/2) . I    -sin(beta/2) . I   ),      |
C |          (  -sin(beta/2) . I     cos(beta/2) . I   )       |
C |                                                            |
C |    D_g = ( exp(-i*gamma/2) . I         0           ),      |
C |          (        0             exp(i*gamma/2) . I )       |
C |                                                            |
C |  where I is the identity matrix of dimension norbs.        |
C |                                                            |
C |  This subroutine exploits the very simple structure of     |
C |  the Dg matrix in an ORTHONORMAL basis to compute          |
C |  products of Dg with a matrix X given by                   |
C |                                                            |
C |    mode = 1,  Dg . X                                       |
C |         = 2,  Dg* . X                                      |
C |         = 3,  transpose(Dg) . X                            |
C |         = 4,  Dg! . X                                      |
C |         = 5,  X . Dg                                       |
C |         = 6,  X . Dg*                                      |
C |         = 7,  X . transpose(Dg)                            |
C |         = 8,  X . Dg!                                      |
C |                                                            |
C |  On output, the matrix X is overwritten with the           |
C |  resulting transformed matrix.                             |
C |                                                            |
C |   ** NOTE: This subroutine expects the matrix X to be      |
C |            fed in full-square form, with dimension         |
C |            4*norbs^2.                                      |
C |                                                            |
C |  scratch space required:                                   |
C |                                                            |
C |    32*nosq,  in all cases                                  |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       norbs  - number of orbitals
C       x      - matrix to transform [ in / out ]
C       alpang - angle alpha
C       betang - angle beta
C       gamang - angle gamma
C       mode   - controls what operation should be performed

      real*8      v(*), x(*), alpang, betang, gamang
      integer     iout, iprint, lenv, norbs
      integer     mode

C     other variables

      integer     jxaa, jxab, jxba, jxbb
      integer     jyaa, jyab, jyba, jybb
      integer     jscr, lscr, jend, nosq
      complex*16  fac1, fac2, fac3, fac4, fac5, fac6
      complex*16  f1, f2, f3, f4, f5, f6, f7, f8


C     Construct useful variables.

      nosq = norbs * norbs


C     Memory allocation.

C       Allocate space to store:
C         - spin blocks of original matrix (jxaa, ...)
C         - spin blocks of transformed matrix (jyaa, ...)
C         - scratch space to call dmblck (jscr)

      lscr = 16*nosq

      jxaa = 1
      jxab = jxaa + 2*nosq
      jxba = jxab + 2*nosq
      jxbb = jxba + 2*nosq
      jyaa = jxbb + 2*nosq
      jyab = jyaa + 2*nosq
      jyba = jyab + 2*nosq
      jybb = jyba + 2*nosq
      jscr = jybb + 2*nosq
      jend = jscr + lscr

      call tstcor (jend-1, lenv, 'mltdg')


C     Compute overlap factors given the angles alpha, beta and gamma.

      fac1 = exp (cmplx (0.0d0, -alpang / 2.0d0))
      fac2 = exp (cmplx (0.0d0,  alpang / 2.0d0))
      fac3 = exp (cmplx (0.0d0, -gamang / 2.0d0))
      fac4 = exp (cmplx (0.0d0,  gamang / 2.0d0))
      fac5 = cmplx (cos (betang / 2.0d0), 0.0d0)
      fac6 = -cmplx (sin (betang / 2.0d0), 0.0d0)


C     Split original matrix into spin-blocks.

      call dmblck (iout, iprint, v(jscr), lscr, norbs, x, v(jxaa),
     $     v(jxab), v(jxba), v(jxbb), 2, -1)


C     Decide where to go depending on mode.

      if ( mode .ge. 1 .and. mode .lt. 5 ) goto 10
      if ( mode .ge. 5 .and. mode .lt. 9 ) goto 20
      call gauerr ('mode not supported in mltdg')



 10   continue

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Form the product  f(Dg) . X  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C       mode = 1,  X  =  Dg . X
C            = 2,  X  =  Dg* . X
C            = 3,  X  =  transpose(Dg) . X
C            = 4,  X  =  Dg! . X


      if ( mode .eq. 1 ) then
        f1 =  fac1*fac3*fac5
        f2 =  fac1*fac4*fac6
        f3 =  fac1*fac3*fac5
        f4 =  fac1*fac4*fac6
        f5 =  fac2*fac4*fac5
        f6 = -fac2*fac3*fac6
        f7 =  fac2*fac4*fac5
        f8 = -fac2*fac3*fac6

      elseif ( mode .eq. 2 ) then
        f1 =  fac2*fac4*fac5
        f2 =  fac2*fac3*fac6
        f3 =  fac2*fac4*fac5
        f4 =  fac2*fac3*fac6
        f5 =  fac1*fac3*fac5
        f6 = -fac1*fac4*fac6
        f7 =  fac1*fac3*fac5
        f8 = -fac1*fac4*fac6

      elseif ( mode .eq. 3 ) then
        f1 =  fac1*fac3*fac5
        f2 = -fac2*fac3*fac6
        f3 =  fac1*fac3*fac5
        f4 = -fac2*fac3*fac6
        f5 =  fac2*fac4*fac5
        f6 =  fac1*fac4*fac6
        f7 =  fac2*fac4*fac5
        f8 =  fac1*fac4*fac6

      elseif ( mode .eq. 4 ) then
        f1 =  fac2*fac4*fac5
        f2 = -fac1*fac4*fac6
        f3 =  fac2*fac4*fac5
        f4 = -fac1*fac4*fac6
        f5 =  fac1*fac3*fac5
        f6 =  fac2*fac3*fac6
        f7 =  fac1*fac3*fac5
        f8 =  fac2*fac3*fac6
      endif


C     aa block

      call amove (2*nosq, v(jxaa), v(jyaa))
      call amove (2*nosq, v(jxba), v(jscr))

      call zscal (nosq, f1, v(jyaa), 1)
      call zscal (nosq, f2, v(jscr), 1)

      call aadd (2*nosq, v(jyaa), v(jscr), v(jyaa))

C     ab block

      call amove (2*nosq, v(jxab), v(jyab))
      call amove (2*nosq, v(jxbb), v(jscr))
      
      call zscal (nosq, f3, v(jyab), 1)
      call zscal (nosq, f4, v(jscr), 1)

      call aadd (2*nosq, v(jyab), v(jscr), v(jyab))

C     ba block

      call amove (2*nosq, v(jxba), v(jyba))
      call amove (2*nosq, v(jxaa), v(jscr))

      call zscal (nosq, f5, v(jyba), 1)
      call zscal (nosq, f6, v(jscr), 1)
      
      call aadd (2*nosq, v(jyba), v(jscr), v(jyba))

C     bb block

      call amove (2*nosq, v(jxbb), v(jybb))
      call amove (2*nosq, v(jxab), v(jscr))

      call zscal (nosq, f7, v(jybb), 1)
      call zscal (nosq, f8, v(jscr), 1)

      call aadd (2*nosq, v(jybb), v(jscr), v(jybb))

      goto 100



 20   continue

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Form the product  X . f(Dg)  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C       mode = 5,  X  =  X . Dg
C            = 6,  X  =  X . Dg*
C            = 7,  X  =  X . transpose(Dg)
C            = 8,  X  =  X . Dg!


      if ( mode .eq. 5 ) then
        f1 =  fac1*fac3*fac5
        f2 = -fac2*fac3*fac6
        f3 =  fac2*fac4*fac5
        f4 =  fac1*fac4*fac6
        f5 =  fac1*fac3*fac5
        f6 = -fac2*fac3*fac6
        f7 =  fac2*fac4*fac5
        f8 =  fac1*fac4*fac6

      elseif ( mode .eq. 6 ) then
        f1 =  fac2*fac4*fac5
        f2 = -fac1*fac4*fac6
        f3 =  fac1*fac3*fac5
        f4 =  fac2*fac3*fac6
        f5 =  fac2*fac4*fac5
        f6 = -fac1*fac4*fac6
        f7 =  fac1*fac3*fac5
        f8 =  fac2*fac3*fac6

      elseif ( mode .eq. 7 ) then
        f1 =  fac1*fac3*fac5
        f2 =  fac1*fac4*fac6
        f3 =  fac2*fac4*fac5
        f4 = -fac2*fac3*fac6
        f5 =  fac1*fac3*fac5
        f6 =  fac1*fac4*fac6
        f7 =  fac2*fac4*fac5
        f8 = -fac2*fac3*fac6

      elseif ( mode .eq. 8 ) then
        f1 =  fac2*fac4*fac5
        f2 =  fac2*fac3*fac6
        f3 =  fac1*fac3*fac5
        f4 = -fac1*fac4*fac6
        f5 =  fac2*fac4*fac5
        f6 =  fac2*fac3*fac6
        f7 =  fac1*fac3*fac5
        f8 = -fac1*fac4*fac6
      endif


C     aa block

      call amove (2*nosq, v(jxaa), v(jyaa))
      call amove (2*nosq, v(jxab), v(jscr))

      call zscal (nosq, f1, v(jyaa), 1)
      call zscal (nosq, f2, v(jscr), 1)

      call aadd (2*nosq, v(jyaa), v(jscr), v(jyaa))

C     ab block

      call amove (2*nosq, v(jxab), v(jyab))
      call amove (2*nosq, v(jxaa), v(jscr))
      
      call zscal (nosq, f3, v(jyab), 1)
      call zscal (nosq, f4, v(jscr), 1)

      call aadd (2*nosq, v(jyab), v(jscr), v(jyab))

C     ba block

      call amove (2*nosq, v(jxba), v(jyba))
      call amove (2*nosq, v(jxbb), v(jscr))

      call zscal (nosq, f5, v(jyba), 1)
      call zscal (nosq, f6, v(jscr), 1)
      
      call aadd (2*nosq, v(jyba), v(jscr), v(jyba))

C     bb block

      call amove (2*nosq, v(jxbb), v(jybb))
      call amove (2*nosq, v(jxba), v(jscr))

      call zscal (nosq, f7, v(jybb), 1)
      call zscal (nosq, f8, v(jscr), 1)

      call aadd (2*nosq, v(jybb), v(jscr), v(jybb))

      goto 100



 100  continue

C     Combine spin-blocks of transformed matrix back into the original
C     location.

      call dmblck (iout, iprint, v(jscr), lscr, norbs, x, v(jyaa),
     $     v(jyab), v(jyba), v(jybb), 1, -1)


      return
      end


