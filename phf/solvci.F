

      subroutine solvci (iout, iprint, v, lenv, nci, hmt, smt, civec)

      implicit none


C +------------------------------------------------------------+
C |  solvci   --   CAJH, 07.2011                               |
C |                                                            |
C |                                    (based on PHFB solvci)  |
C |                                                            |
C |                                                            |
C |  Solve the CI problem associated with spin projection and  |
C |  / or complex conjugation restoration. The CI              |
C |  coefficients obtained by diagonalizing the Hamiltonian    |
C |  matrix are used to compute the projected energy of the    |
C |  model.                                                    |
C |                                                            |
C |  The coefficients corresponding to the lowest eigenvalue   |
C |  are returned in the vector civec.                         |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  The projected energy is obtained by diagonalizing the     |
C |  Hamiltonian matrix in the basis { P^S_{MK} I |0>,         |
C |  P^S_{MK} K | 0>} (if both spin projection and complex     |
C |  conjugation restoration are performed).                   |
C |                                                            |
C |  Here, M,K = -S, -S+1, ..., S-1, S, where S is the         |
C |  quantum number associated with the S^2 operator. Note     |
C |  that if the wavefunction is of UHF type, then one does    |
C |  not need to solve the CI problem, and M,K = S = Na - Nb,  |
C |  where Na and Nb are the number of alpha and beta          |
C |  electrons, respectively.                                  |
C |                                                            |
C |  Additionally, the operators {I, K} (identity and complex  |
C |  conjugation) define the basis of states necessary to      |
C |  perform complex conjugation restoration.                  |
C |                                                            |
C |  If ncis is the number of basis states required for spin   |
C |  projection, then the total number of CI basis states is   |
C |  given by nci = 2 * ncis if complex conjugation            |
C |  restoration is also being performed.                      |
C |                                                            |
C |  The Hamiltonian and overlap matrices that constitute the  |
C |  CI problem are thus given by                              |
C |                                                            |
C |    Hmt(i,j) =  < 0 | nu! H P_{MK} (Omega) mu | 0 >,        |
C |                                                            |
C |    Smt(i,j) =  < 0 | nu! P_{MK} (Omega) mu | 0 >,          |
C |                                                            |
C |  where mu,nu = I,K, and P^S_{MK} is the spin projection    |
C |  operator. The indices i,j loop over all possible          |
C |  combination of indices nu, mu, M, K.                      |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       nci    - number of basis states in CI problem
C       hmt    - Hamiltonian matrix in CI problem
C       smt    - overlap matrix in CI problem
C       civec  - vector of CI coefficients [ out ]

      real*8      v(*)
      integer     iout, iprint, lenv, nci
      complex*16  hmt(nci,nci), smt(nci,nci), civec(*)

C     other variables

      real*8      fac
      integer     ncit, lind, info
      integer     i, ind1, ind2
      integer     jhn, jsmt, jscr, jscr1, jscr2, jend
      integer     jsvec, jsval, jhvec, jhval

C     constants

      real*8      thresh
      complex*16  zero, one

      parameter   ( thresh = 1.0d-12 )



      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)


C     Print overlap and Hamiltonian matrices if requested.

      if ( iprint .gt. 1 ) then
        call ltoutx (iout, 1, 2, 1, 0, 'Smt', nci, nci, smt, smt, 0)
        call ltoutx (iout, 1, 2, 1, 0, 'Hmt', nci, nci, hmt, hmt, 0)
      endif


C     Return if the number of CI states is equal to 1.

      if ( nci .eq. 1 ) then
        civec(1) = one / sqrt(dble(smt(1,1)))
        return
      endif


C     Define useful quantities.

      ncit = nci * (nci+1) / 2


C     Memory allocation.

C       Allocate space for:
C         - transformed Hamiltonian matrix (jhn)
C         - overlap matrix in lower triangular form (jsmt)
C         - eigenvalues and vectors of overlap matrix (jsval, jsvec)
C         - eigenvalues and vectors of Hamiltonian matrix (jhval, jhvec)
C         - scratch array to perform various operations (jscr)
C         - scratch arrays needed for zhpev (jscr1, jscr2)

      jhn   = 1
      jsmt  = jhn   + 2*nci*nci
      jsval = jsmt  + 2*ncit
      jsvec = jsval + nci
      jhval = jsvec + 2*nci*nci
      jhvec = jhval + nci
      jscr  = jhvec + 2*nci*nci
      jscr1 = jscr  + 2*nci*nci
      jscr2 = jscr1 + 2*(2*nci-1)
      jend  = jscr2 + 3*nci-2

      call tstcor (jend-1, lenv, 'solvci')



C     Hermitize the overlap matrix.
C     Transform into lower triangular.

      call apadag (1, nci, 0.50d0, smt, v(jscr))

      call trspnc (.false., nci, nci, v(jscr))
      call zlinear (v(jscr), v(jsmt), nci, nci)

C     Apply a factor of -1 (so that larger eigenvalues appear first).

      call zdscal (ncit, -1.0d0, v(jsmt), 1)


C     Diagonalize the overlap matrix.

C       Recall that we store the lower triangle as an upper triangle. To
C       make it a true upper triangle we thus take the complex
C       conjugate. We then send to LaPack as a true upper triangle.

      call vconjg (ncit, v(jsmt))

      call zhpev ('v', 'u', nci, v(jsmt), v(jsval), v(jsvec), nci,
     $     v(jscr1), v(jscr2), info)

      if ( info .ne. 0 ) then
        call gauerr ('zhpev failed in solvci')
      endif


C     Form the transformation matrix

C       X = U . s^(-1/2),

C     where U is the matrix of eigenvectors of S, and s^(-1/2) is the
C     inverse of the sqrt of the diagonalized overlap matrix. The
C     transformation matrix is stored in jsvec.

C     lind counts the number of linearly independent vectors.

      lind = 0

      do 10 i = 1, nci
        ind1 = (i-1)*2*nci + jsvec
        ind2 = i-1 + jsval

        if ( abs (v(ind2)) .gt. thresh ) then
          fac = sqrt (-1.0d0 / v(ind2))
          lind = lind + 1
        else
          fac = 0.0d0
        endif

        call zdscal (nci, fac, v(ind1), 1)
 10   continue


C     Transform the Hamiltonian matrix to orthonormal basis.

C       Hn = X' . H . X.

      call zgemm ('c', 'n', lind, nci, nci, one, v(jsvec), nci, hmt,
     $     nci, zero, v(jscr), lind)
      call zgemm ('n', 'n', lind, lind, nci, one, v(jscr), lind,
     $     v(jsvec), nci, zero, v(jhn), lind)


C     Hermitize the Hamiltonian matrix.
C     Transform into lower triangular.

      call apadag (1, lind, 0.50d0, v(jhn), v(jscr))

      if ( lind .gt. 1 ) then
        call trspnc (.false., lind, lind, v(jscr))
        call zlinear (v(jscr), v(jhn), lind, lind)
      endif


C     Diagonalize the Hamiltonian matrix if the number of linearly
C     independent states is greater than 1.

      if ( lind .eq. 1 ) then
        v(jhvec) = 1.0d0
        v(jhvec+1) = 0.0d0

      else
        call vconjg (lind*lind, v(jhn))

        call zhpev ('v', 'u', lind, v(jhn), v(jhval), v(jhvec), lind,
     $       v(jscr1), v(jscr2), info)

        if ( info .ne. 0 ) then
          call gauerr ('zhpev failed in solvci')
        endif
      endif


C     Transform the Hamiltonian eigenvectors into the original basis.

      call zgemm ('n', 'n', nci, lind, lind, one, v(jsvec), nci,
     $     v(jhvec), lind, zero, v(jscr), nci)


C     Save CI coefficients for the lowest eigenvalue.

      call amove (2*nci, v(jscr), civec)


      return
      end


