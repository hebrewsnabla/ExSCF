

      subroutine evlygg (iout, iprint, itype, n1c, n2c, na, nb, norbs,
     $     detC, detC2, inobas, dodiag, piv, eig, scr1, scr2, lscr1,
     $     ngg, xgg)

      implicit none


C +------------------------------------------------------------+
C |  evlygg   --   CAJH, 08.2011                               |
C |                                                            |
C |                                                            |
C |  Evaluate the overlap-related weighting function x(g,g')   |
C |  in a similar way as evalsmt builds the function x(g).     |
C |  x(g,g') is then used to build y(g,g'), hence the name of  |
C |  the subroutine. Such integration, however, is carried     |
C |  out outside of here.                                      |
C |                                                            |
C |  The function x(g,g') is given by (possible cases)         |
C |                                                            |
C |    x(g,g'; i,j) =                                          |
C |             [D^S_{K'M} (Omega')] . [D^S_{MK} (Omega)]* .   |
C |                                                            |
C |     { <phi| I! R_{MK'} (Omega')! R_{MK} (Omega) I | phi>   |
C |     { <phi| I! R_{MK'} (Omega')! R_{MK} (Omega) K | phi>   |
C |     { <phi| K! R_{MK'} (Omega')! R_{MK} (Omega) I | phi>   |
C |     { <phi| K! R_{MK'} (Omega')! R_{MK} (Omega) K | phi>   |
C |                                                            |
C |  Here, R_{MK} (Omega) are rotation operators, g == Omega   |
C |  defines the space of the projection operator, and         |
C |  D^S_{MK} (Omega) is Wigner's rotation matrix. S and M     |
C |  denote spin quantum numbers of the projected state.       |
C |                                                            |
C |  The indices i,j are used to select the appropriate K,K'   |
C |  indices in the spin projection operator and the           |
C |  appropriate operators in the matrix element.              |
C |                                                            |
C |    Note that the number of CI states (and thus of x(g,g')  |
C |    functions) depends on whether spin projection and       |
C |    complex conjugation restoration are being performed     |
C |    or not.                                                 |
C |                                                            |
C |  The matrix elements are evaluated as                      |
C |                                                            |
C |    <phi| I! R(Omega')! R(Omega) I |phi>  =                 |
C |                       det ( C! R(Omega')! R(Omega) C )     |
C |    <phi| I! R(Omega')! R(Omega) K |phi>  =                 |
C |                       det ( C! R(Omega')! R(Omega) C* )    |
C |    <phi| K! R(Omega')! R(Omega) I |phi>  =                 |
C |                       det ( C^T R(Omega')! R(Omega) C )    |
C |    <phi| K! R(Omega')! R(Omega) K |phi>  =                 |
C |                       det ( C^T R(Omega')! R(Omega) C* ),  |
C |                                                            |
C |  where |phi> is the reference HF state, ^T indicates       |
C |  transpose and C is the matrix of occupied orbital         |
C |  coefficients (thus C is an M x N matrix).                 |
C |                                                            |
C |  We then evaluate the function y(g,g') from the values     |
C |  x(g,g') and the CI coefficients stored in civec.          |
C |  The weighting function y(g,g') is given by                |
C |                                                            |
C |                        c_i* c_j  x(g,g'; i,j)              |
C |  y(g,g'; i,j) =  ---------------------------------------,  |
C |                sum_{k,l} c_k* c_l  int dg  x(g,g'; k,l)    |
C |                                                            |
C |  where there is a parametric dependence on the i and j     |
C |  indices in the functions x(g,g') and y(g,g').             |
C |                                                            |
C |  Here, c_i are CI coefficients obtained from solving the   |
C |  CI problem described in solvci, that is, the              |
C |  diagonalization of the Hamiltonian matrix in the basis    |
C |  basis of the nci states.                                  |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  CAJH, 08.2011                                             |
C |                                                            |
C |  Major modification of the subroutine. It now computes     |
C |  the weighting factor x(g,g') for a single integration     |
C |  point. The integration over gauge angles and the          |
C |  multiplication times appropriate CI coefficients is       |
C |  done in the parent subroutine phfdm1.F.                   |
C |                                                            |
C |  dimension of scratch arrays:                              |
C |                                                            |
C |    itype = 1, 2,                                           |
C |                                                            |
C |      piv  - intowp (na)                                    |
C |      eig  - 2*na                                           |
C |      scr1 - max (2*na*na, 4*na)                            |
C |      scr2 - 2*na                                           |
C |                                                            |
C |    itype = 3,                                              |
C |                                                            |
C |      piv  - intowp (nel)                                   |
C |      eig  - 2*nel                                          |
C |      scr1 - max (2*nel*nel, 4*nel)                         |
C |      scr2 - 2*nel                                          |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  CAJH, 08.2011                                             |
C |                                                            |
C |  Use of the true NO basis has been enabled even for        |
C |  methods with complex conjugation.                         |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  CAJH, 12.2012                                             |
C |                                                            |
C |  Removed the evaluation of Wigner's matrix from here.      |
C |  A lot of input variables have thus been removed.          |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       itype  - type of PHF calculation (see phfdrv for details)
C       n?c    - indices for complex conjugation projection
C       na     - number of alpha electrons
C       nb     - number of beta electrons
C       norbs  - number of orbitals
C       detC   - determinant of occupied-occupied block of MOs
C       detC2  - determinant of occupied-occupied (C*) block of MOs
C       inobas - if .true., force the use of the true NO basis
C       dodiag - whether to compute determinant by diagonalization
C       piv    - scratch array
C       eig    - scratch array
C       scr?   - scratch arrays
C       lscr1  - length of scratch array 1, in real*8 words
C       ngg    - matrix N(g,g')
C       xgg    - weigthing factor x(g,g')  [ out ]

      integer     iout, iprint, itype, n1c, n2c, na, nb, norbs
      integer     lscr1
      real*8      piv(*), eig(*), scr1(*), scr2(*), ngg(*)
      complex*16  detC, detC2, xgg
      logical     inobas, dodiag

C     other variables

      real*8      fac1, fac2, fac3, fac4
      integer     nel, nosq, info
      complex*16  detN

C     constants

      complex*16  zero, one, im


      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)
      im   = cmplx (0.0d0, 1.0d0)


C     Construct useful variables.

      nel = na + nb
      nosq = norbs * norbs



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Evaluate the overlap matrix elements  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Evaluate the overlap matrix elements as described above...


C     Compute the quantity
C       det ( N(g,g') )

      if ( itype .eq. 1 ) then

        if ( .not. dodiag ) then
          call zgetrf (na, na, ngg, na, piv, info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evlygg')

          call detmat (na, ngg, piv, fac1, fac2)

        else
          call zgeev ('n', 'n', na, ngg, na, eig, zero, 1, zero,
     $         1, scr1, lscr1/2, scr2, info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evlygg')

          call detmt2 (na, eig, fac1, fac2)
        endif

        detN = cmplx (fac1) * exp (im * cmplx (fac2))
        detN = detN**2

      elseif ( itype .eq. 2 ) then

        if ( .not. dodiag ) then
          call zgetrf (na, na, ngg, na, piv, info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evlygg')

          call detmat (na, ngg, piv, fac1, fac2)

          if ( nb .ge. 1 ) then
          call zgetrf (nb, nb, ngg(1+2*na*na), nb, piv, info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evlygg')

          call detmat (nb, ngg(1+2*na*na), piv, fac3, fac4)
          else
            fac3 = 1.0d0
            fac4 = 0.0d0
          endif

        else
          call zgeev ('n', 'n', na, ngg, na, eig, zero, 1,
     $         zero, 1, scr1, lscr1/2, scr2, info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evlygg')

          call detmt2 (na, eig, fac1, fac2)

          if ( nb .ge. 1 ) then
          call zgeev ('n', 'n', nb, ngg(1+2*na*na), nb, eig, zero,
     $         1, zero, 1, scr1, lscr1/2, scr2, info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evlygg')

          call detmt2 (nb, eig, fac3, fac4)
          else
            fac3 = 1.0d0
            fac4 = 0.0d0
          endif
        endif

        detN = cmplx (fac1) * exp (im * cmplx (fac2))
        detN = detN * cmplx (fac3) * exp (im * cmplx (fac4))

      elseif ( itype .eq. 3 ) then

        if ( .not. dodiag ) then
          call zgetrf (nel, nel, ngg, nel, piv, info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evlygg')

          call detmat (nel, ngg, piv, fac1, fac2)

        else
          call zgeev ('n', 'n', nel, ngg, nel, eig, zero, 1,
     $         zero, 1, scr1, lscr1/2, scr2, info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evlygg')

          call detmt2 (nel, eig, fac1, fac2)
        endif

        detN = cmplx (fac1) * exp (im * cmplx (fac2))
      endif


C     Compute factor x(g,g').

      if ( .not. inobas ) then
        if ( n1c .eq. 1 .and. n2c .eq. 1 ) then
          xgg = 1.0d0 / (detN * detC * conjg (detC))

        elseif ( n1c .eq. 1 .and. n2c .eq. 2 ) then
          xgg = 1.0d0 / (detN * detC * detC)

        elseif ( n1c .eq. 2 .and. n2c .eq. 1 ) then
          xgg = 1.0d0 / (detN * conjg (detC) * conjg (detC))

        elseif ( n1c .eq. 2 .and. n2c .eq. 2 ) then
          xgg = 1.0d0 / (detN * conjg (detC) * detC)
        endif

      elseif ( inobas ) then
        if ( n1c .eq. 1 .and. n2c .eq. 1 ) then
          xgg = 1.0d0 / (detN * detC * conjg (detC))

        elseif ( n1c .eq. 1 .and. n2c .eq. 2 ) then
          xgg = 1.0d0 / (detN * detC * conjg (detC2))

        elseif ( n1c .eq. 2 .and. n2c .eq. 1 ) then
          xgg = 1.0d0 / (detN * detC2 * conjg (detC))

        elseif ( n1c .eq. 2 .and. n2c .eq. 2 ) then
          xgg = 1.0d0 / (detN * detC2 * conjg (detC2))
        endif
      endif


      return
      end


