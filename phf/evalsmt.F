

      subroutine evalsmt (iout, iprint, v, lenv, iphf, nci, ncis, ncik, 
     $     ngrdb, norbs, na, nb, irwnos, inobas, ngrdt, indarr, grda,
     $     wgta, grdb, wgtb, grdg, wgtg, dmt, dind, irwvec, irwng, yg, 
     $     smt)

      implicit none


C +------------------------------------------------------------+
C |  evalsmt   --   CAJH, 07.2011                              |
C |                                                            |
C |                                   (based on PHFB evalsmt)  |
C |                                                            |
C |                                                            |
C |  Evaluate the overlap matrix elements (Smt) needed to      |
C |  solve the CI problem associated with complex conjugation  |
C |  restoration and spin projection.                          |
C |                                                            |
C |  Additionally, this routine evaluates the weighting        |
C |  function x(g) appearing in the projected energy           |
C |  expression.                                               |
C |                                                            |
C |  The overlap matrix elements are given by                  |
C |                                                            |
C |    Smt(i,j) = int dg x(g; i,j),                            |
C |                                                            |
C |  where there is a parametric dependence on the i and j     |
C |  indices in the function x(g).                             |
C |                                                            |
C |  In other words, we have nci*nci different functions       |
C |  x(g), where nci is the number of CI states to consider.   |
C |                                                            |
C |  The function x(g) is given by (possible cases)            |
C |                                                            |
C |    x(g; i,j) =                                             |
C |                        { <phi| I! R_{MK} (Omega) I | phi>  |
C |    [D^S_{MK} (Omega)]* { <phi| I! R_{MK} (Omega) K | phi>  |
C |                        { <phi| K! R_{MK} (Omega) I | phi>  |
C |                        { <phi| K! R_{MK} (Omega) K | phi>  |
C |                                                            |
C |  Here, R_{MK} (Omega) are rotation operators, g == Omega   |
C |  defines the space of the projection operator, and         |
C |  D^S_{MK} (Omega) is Wigner's rotation matrix.             |
C |                                                            |
C |  The indices i,j are used to select the appropriate M,K    |
C |  indices in the spin projection operator and the           |
C |  appropriate operators in the matrix element.              |
C |                                                            |
C |    Note that the number of CI states (and thus of x(g)     |
C |    functions) depends on whether spin projection and       |
C |    complex conjugation restoration are being performed     |
C |    or not.                                                 |
C |                                                            |
C |  The matrix elements are evaluated as                      |
C |                                                            |
C |    <phi| I! R(Omega) I |phi>  =  det ( C! R(Omega) C )     |
C |    <phi| I! R(Omega) K |phi>  =  det ( C! R(Omega) C* )    |
C |    <phi| K! R(Omega) I |phi>  =  det ( C^T R(Omega) C )    |
C |    <phi| K! R(Omega) K |phi>  =  det ( C^T R(Omega) C* ),  |
C |                                                            |
C |  where |phi> is the reference HF state, ^T indicates       |
C |  transpose and C is the matrix of occupied orbital         |
C |  coefficients (thus C is an M x N matrix).                 |
C |                                                            |
C |  Details of how to turn these matrix elements into         |
C |  functionals of the density matrix, and how we compute     |
C |  them as such, will be hopefully published soon.           |
C |                                                            |
C |    Note that we could also evaluate the overlap matrix     |
C |    elements as stated above. We go through the more        |
C |    convoluted way to see if there is anything breaking     |
C |    down in our development.                                |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  NOTE: This routine evaluates the overlap factors x(g)     |
C |        and stores them in the vector yg. Other routines    |
C |        will convert the x(g) into y(g).                    |
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
C |  Reworked subroutine to handle the CI over spin            |
C |  internally.                                               |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       nci    - number of basis states in CI problem
C       ncis   - number of basis states for spin projection
C       ncik   - number of basis states for complex conj projection
C       ngrdb  - number of grid points in beta integration
C       norbs  - number of orbitals
C       na     - number of alpha electrons
C       nb     - number of beta electrons
C       irwnos - RWF file with transformation matrix from orthonormal AO
C                to NO basis
C       inobas - if .true., force the use of the true NO basis
C       ngrdt  - total number of grid points
C       indarr - mapping array igrdt => igrdX
C       grdX   - integration grid for each variable
C       wgtX   - array with quadrature weights for each variable
C       dmt    - Wigner's small d-matrix d^S_{MK} (beta)
C       dind   - vector with indices {M} for Wigner rotation matrices
C       irwvec - RWF file with orbital coefficients
C                (orthonormal AO basis)
C       irwng  - RWF file with matrices N(g)
C       yg     - vector of x(g) values [ out ]
C       smt    - overlap matrix in CI problem [ out ]

      real*8      grda(*), grdb(*), grdg(*)
      real*8      wgta(*), wgtb(*), wgtg(*)
      integer     iout, iprint, lenv, iphf, nci, ncis, ncik
      integer     ngrdb, na, nb
      integer     norbs, ngrdt, indarr(ngrdt,*)
      integer     irwnos, irwvec, irwng
      logical     inobas
      complex*16  yg(*), smt(nci,nci)

      real*8      v(*), dmt(ngrdb,ncis,ncis), dind(ncis)

C     other variables

      integer     nel, nosq, szn, szbfn, itype, ind
      integer     jorb, jorb2, jnos, jnost, jend
      integer     jzab, jz, jscrd, lscrd
      integer     ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8
      integer     info, ic, ica, icb
      integer     lpiv, leig, lscr1, lscr2
      integer     lenp, npmem, npmax, npuse, ipar
      integer     n1, n2, n1c, n2c, n1s, n2s, ip, iblock, nthrd
      logical     icmplx, ispuhf, ispghf, throk
      complex*16  detC, detC2, sums
      character   indc1*2, indc2*2, str*10

      real*8      facw, fac1, fac2, fac3, fac4
      integer     igrdt, iint, ibufn, ing, iend, inga, ingb
      integer     ipiv, ieig, iscr1, iscr2
      integer     indy, indi, igrda, igrdb, igrdg
      complex*16  facd, detN

C     functions

      integer     intowp, nproc, lappar

C     dodiag - whether to perform diagonalization to compute
C              determinants

      logical     dodiag
      save        dodiag
      data        dodiag / .true. /

C     constants

      complex*16  zero, one, im


      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)
      im   = cmplx (0.0d0, 1.0d0)



C     Construct useful variables.

      nel = na + nb
      nosq = norbs * norbs


C     Decipher iphf.

C       icmplx - whether complex conjugation projection is on
C       ispuhf - whether UHF-based spin projection is used
C       ispghf - whether GHF-based spin projection is used
C       itype  - decides which spin blocks of transition density
C                matrices are stored

C       itype  = 1,  closed shell [alpha-alpha block]
C              = 2,  open shell [alpha-alpha, beta-beta blocks]
C              = 3,  general [all spin blocks are active]

      icmplx = iphf .le. 3 .or. iphf .eq. 5 .or. iphf .eq. 7
      ispuhf = iphf .eq. 4 .or. iphf .eq. 5
      ispghf = iphf .eq. 6 .or. iphf .eq. 7

      itype = 1

      if ( iphf .eq. 2 ) itype = 2
      if ( iphf .ge. 3 ) itype = 3

C     Error checking.

      if ( iphf .lt. 1 .or. iphf .gt. 7 ) then
        call gauerr ('Unrecognized iphf in evalsmt.')
      endif


C     Define useful array dimensions:
C       - szn,  size of N(g) matrices

      if ( itype .eq. 1 ) then
        szn = 2*na*na
      elseif ( itype .eq. 2 ) then
        szn = 2*na*na + 2*nb*nb
      elseif ( itype .eq. 3 ) then
        szn = 2*nel*nel
      endif

C       .. buffer to load computed Ng ..

      szbfn = ncik*ncik * szn


C     Memory allocation.

C       Allocate space for:
C         - matrix of orbital coefficients in orthonormal AO basis (jorb)
C         - orthonormal AO => NO basis transformation matrix (jnos)
C         - scratch array to perform basis transformation (jnost)

      if ( itype .eq. 1 ) then
        jorb  = 1
        jnos  = jorb  + 2*nosq
        jnost = jnos  + 2*nosq
        jend  = jnost + 2*nosq

      elseif ( itype .eq. 2 ) then
        jorb  = 1
        jnos  = jorb  + 4*nosq
        jnost = jnos  + 4*nosq
        jend  = jnost + 2*nosq

      elseif ( itype .eq. 3 ) then
        jorb  = 1
        jnos  = jorb  + 8*nosq
        jnost = jnos  + 8*nosq
        jend  = jnost + 8*nosq
      endif

C       Allocate space for:
C         - second copy of matrix of orbital coefficients (jorb2)

      if ( inobas .and. itype .eq. 1 ) then
        jorb2 = jend
        jend  = jorb2 + 2*nosq

      elseif ( inobas .and. itype .eq. 2 ) then
        jorb2 = jend
        jend  = jorb2 + 4*nosq

      elseif ( inobas .and. itype .eq. 3 ) then
        jorb2 = jend
        jend  = jorb2 + 8*nosq
      endif


C     Allocate scratch space for:
C       - scratch arrays needed for ispuhf = .true. (see below)

      if ( ispuhf ) then
        lscrd = 16*nosq

        jzab  = jend
        jscrd = jzab  + 2*nosq
        jz    = jscrd + lscrd
        jend  = jz    + 8*nosq
      endif


C     Do some preparation for parallel processing:
C       - lenp,  dimension of scratch space per processor
C       - npuse, number of processors to use

      lpiv = 0
      leig = 0
      lscr1 = 0
      lscr2 = 0

      if ( itype .eq. 1 .or. itype .eq. 2 ) then
        if ( .not. dodiag ) lpiv = intowp (na)

        if ( dodiag ) leig  = 2*na
        if ( dodiag ) lscr1 = max (2*na*na, 4*na)
        if ( dodiag ) lscr2 = 2*na

      elseif ( itype .eq. 3 ) then
        if ( .not. dodiag ) lpiv = intowp (nel)

        if ( dodiag ) leig  = 2*nel
        if ( dodiag ) lscr1 = max (2*nel*nel, 4*nel)
        if ( dodiag ) lscr2 = 2*nel
      endif

      lenp = 2*nci*nci + szbfn + lpiv + leig + lscr1 + lscr2

      npmem = (lenv - jend + 1) / lenp

      if ( npmem .lt. 1 ) then
        call gauerr ('Not enough memory in evalsmt for a single CPU.')
      endif

C      write(iout, *) 'End mem in smt'
      npmax = nproc (0)
      npuse = min (npmem, npmax)
      
      ipar  = lappar (0)
      throk = ipar .eq. 1

      if ( ipar .eq. -1 ) npuse = 1

      call tstcor (jend + lenp*npuse - 1, lenv, 'evalsmt')


C      write(iout, *) 'End tstcor in smt'
C     Load matrix of orbital coefficients.

      if ( itype .eq. 1 ) then
        call fileio (2, -irwvec, 2*nosq, v(jorb), 0)

      elseif ( itype .eq. 2 ) then
        call fileio (2, -irwvec, 4*nosq, v(jorb), 0)

      elseif ( itype .eq. 3 ) then

        if ( .not. ispuhf ) then
          call fileio (2, -irwvec, 8*nosq, v(jorb), 0)
        else

C           For ispuhf = .true., we only have the alpha-alpha and
C           beta-beta blocks of the MO coefficients. We prepare the full
C           matrix here.

          call fileio (2, -irwvec, 4*nosq, v(jnost), 0)

          call aclear (2*nosq, v(jzab))

          call dmblck (iout, iprint, v(jscrd), lscrd, norbs, v(jz),
     $         v(jnost), v(jzab), v(jzab), v(jnost+2*nosq), 1, -1)
          write (iout, *) ' '
          call ltoutx (iout, 1, 2, 1, 0, 'C ortho', 2*norbs,
     $           2*norbs, v(jz), v(jz), 0)


C           Reorganize orbitals. See hfdmdg for details.

          ind1 = jz
          ind2 = jz + 4*nosq
          ind3 = jz + 4*norbs*na
          ind4 = jz + 4*nosq + 4*norbs*nb

          ind5 = jorb
          ind6 = jorb + 4*norbs*na
          ind7 = jorb + 4*norbs*nel
          ind8 = jorb + 4*nosq + 4*norbs*nb

          call amove (4*norbs*na, v(ind1), v(ind5))
          call amove (4*norbs*nb, v(ind2), v(ind6))
          call amove (4*norbs*(norbs-na), v(ind3), v(ind7))
          call amove (4*norbs*(norbs-nb), v(ind4), v(ind8))
          write (iout, *) ' '
          call ltoutx (iout, 1, 2, 1, 0, 'C ortho reorganized', 2*norbs,
     $           2*norbs, v(jorb), v(jorb), 0)
        endif
      endif

C        write(iout, *) 'End fileio in smt'

C     Create a copy of jorb into jorb2. Store orb* in orb2.

      if ( inobas .and. itype .eq. 1 ) then
        call amove (2*nosq, v(jorb), v(jorb2))
        call vconjg (nosq, v(jorb2))

      elseif ( inobas .and. itype .eq. 2 ) then
        call amove (4*nosq, v(jorb), v(jorb2))
        call vconjg (2*nosq, v(jorb2))

      elseif ( inobas .and. itype .eq. 3 ) then
        call amove (8*nosq, v(jorb), v(jorb2))
        call vconjg (4*nosq, v(jorb2))
      endif


C        write(iout, *) 'End copy in smt'
C     Load transformation matrix (orthonormal AO => NO basis).

      if ( itype .eq. 1 ) then
        call fileio (2, -irwnos, 2*nosq, v(jnos), 0)

      elseif ( itype .eq. 2 ) then
        call fileio (2, -irwnos, 4*nosq, v(jnos), 0)

      elseif ( itype .eq. 3 ) then
        call fileio (2, -irwnos, 8*nosq, v(jnos), 0)
      endif


C     Transform matrix of orbital coefficients to NO basis:
C       nos' * orb

      if ( itype .eq. 1 ) then

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jorb), norbs, zero, v(jnost), norbs)

        call amove (2*nosq, v(jnost), v(jorb))

      elseif ( itype .eq. 2 ) then

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jorb), norbs, zero, v(jnost), norbs)

        call amove (2*nosq, v(jnost), v(jorb))

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos+2*nosq),
     $       norbs, v(jorb+2*nosq), norbs, zero, v(jnost), norbs)

        call amove (2*nosq, v(jnost), v(jorb+2*nosq))

      elseif ( itype .eq. 3 ) then

        call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jnos),
     $       2*norbs, v(jorb), 2*norbs, zero, v(jnost), 2*norbs)

        call amove (8*nosq, v(jnost), v(jorb))
          write (iout, *) ' '
          call ltoutx (iout, 1, 2, 1, 0, 'C (NO)', 2*norbs,
     $           2*norbs, v(jnost), v(jnost), 0)
      endif


C     Transform matrix of orbital coefficients (orb*) to NO basis:
C       nos' * orb*

      if ( inobas .and. itype .eq. 1 ) then

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jorb2), norbs, zero, v(jnost), norbs)

        call amove (2*nosq, v(jnost), v(jorb2))

      elseif ( inobas .and. itype .eq. 2 ) then

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jorb2), norbs, zero, v(jnost), norbs)

        call amove (2*nosq, v(jnost), v(jorb2))

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos+2*nosq),
     $       norbs, v(jorb2+2*nosq), norbs, zero, v(jnost), norbs)

        call amove (2*nosq, v(jnost), v(jorb2+2*nosq))

      elseif ( inobas .and. itype .eq. 3 ) then

        call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jnos),
     $       2*norbs, v(jorb2), 2*norbs, zero, v(jnost), 2*norbs)

        call amove (8*nosq, v(jnost), v(jorb2))
      endif

c        write(iout, *) 'End transform in smt'


C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Evaluate the determinant of C  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     The procedure we use to evaluate the overlap matrix elements
C     relies on computing the inverse of the occupied-occupied block of
C     the matrix of orbital coefficients. To understand what we mean,
C     let us write the matrix of orbital coefficients, D, as

C       D =  (  D_oo  D_ov  )
C            (  D_vo  D_vv  )

C     Here, D is an M x M matrix, while the block D_oo is of size N x N,
C     where N is the number of electrons. Clearly, D_oo depends on the
C     basis in which we write down D. We use a basis in which D_oo has
C     an inverse. Let us now call C = D_oo.

C     The overlap matrix elements appearing in x(g) (see above) can then
C     be computed as

C       <phi| I! R(Omega) I |phi>  =  1 / ( det ( C . N(g) . C! ),
C       <phi| I! R(Omega) K |phi>  =  1 / ( det ( C . N(g) . C^T ),
C       <phi| K! R(Omega) I |phi>  =  1 / ( det ( C* . N(g) . C! ),
C       <phi| K! R(Omega) K |phi>  =  1 / ( det ( C* . N(g) . C^T ),

C     where the matrices N(g) are those given in formng.

C     We now proceed to evaluate the determinant of C...


C     Memory allocation.

      if ( .not. dodiag ) then
        ic   = jend + 2
        ipiv = ic   + szn
        iend = ipiv + lpiv

      else
        ic    = jend  + 2
        ieig  = ic    + szn
        iscr1 = ieig  + leig
        iscr2 = iscr1 + lscr1
        iend  = iscr2 + lscr2
      endif


C     Useful quantities for itype = 2.

      if ( itype .eq. 2 ) then
        ica = ic
        icb = ic + 2*na*na
      endif


C     Extract C (occupied-occupied block) from full matrix of orbital
C     coefficients.

      if ( itype .eq. 1 ) then
        call sptmat (2, norbs, norbs, v(jorb), na, na, v(ic), 1)

      elseif ( itype .eq. 2 ) then
        call sptmat (2, norbs, norbs, v(jorb), na, na, v(ica), 1)

        if ( nb .ge. 1 ) then
        call sptmat (2, norbs, norbs, v(jorb+2*nosq), nb, nb, v(icb), 1)
        endif

      elseif ( itype .eq. 3 ) then
        call sptmat (2, 2*norbs, 2*norbs, v(jorb), nel, nel, v(ic), 1)
      endif


c        write(iout, *) 'End extract C in smt'
C     Compute the quantity
C       detC = det ( C )

      if ( itype .eq. 1 ) then

        if ( .not. dodiag ) then
          call zgetrf (na, na, v(ic), na, v(ipiv), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

          call detmat (na, v(ic), v(ipiv), fac1, fac2)

        else
          call zgeev ('n', 'n', na, v(ic), na, v(ieig), zero, 1, zero,
     $         1, v(iscr1), lscr1/2, v(iscr2), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

          call detmt2 (na, v(ieig), fac1, fac2)
        endif

        detC = cmplx (fac1) * exp (im * cmplx (fac2))
        detC = detC**2

      elseif ( itype .eq. 2 ) then

        if ( .not. dodiag ) then
          call zgetrf (na, na, v(ica), na, v(ipiv), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

          call detmat (na, v(ica), v(ipiv), fac1, fac2)

          if ( nb .ge. 1 ) then
          call zgetrf (nb, nb, v(icb), nb, v(ipiv), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

          call detmat (nb, v(icb), v(ipiv), fac3, fac4)
          else
            fac3 = 1.0d0
            fac4 = 0.0d0
          endif

        else
          call zgeev ('n', 'n', na, v(ica), na, v(ieig), zero, 1, zero,
     $         1, v(iscr1), lscr1/2, v(iscr2), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

          call detmt2 (na, v(ieig), fac1, fac2)

          if ( nb .ge. 1 ) then
          call zgeev ('n', 'n', nb, v(icb), nb, v(ieig), zero, 1, zero,
     $         1, v(iscr1), lscr1/2, v(iscr2), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

          call detmt2 (nb, v(ieig), fac3, fac4)
          else
            fac3 = 1.0d0
            fac4 = 0.0d0
          endif
        endif

        detC = cmplx (fac1) * exp (im * cmplx (fac2))
        detC = detC * cmplx (fac3) * exp (im * cmplx (fac4))

      elseif ( itype .eq. 3 ) then

        if ( .not. dodiag ) then
          call zgetrf (nel, nel, v(ic), nel, v(ipiv), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

          call detmat (nel, v(ic), v(ipiv), fac1, fac2)

        else
C          write(iout, *) 'Start zgeev in smt'
          call zgeev ('n', 'n', nel, v(ic), nel, v(ieig), zero, 1, zero,
     $         1, v(iscr1), lscr1/2, v(iscr2), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

C          write(iout, *) 'Start detmt2 in smt'
          call detmt2 (nel, v(ieig), fac1, fac2)
        endif

        detC = cmplx (fac1) * exp (im * cmplx (fac2))
        write(iout, 66) detC
 66   Format("DetC = ", F10.5, "+", F10.5, "i")
      endif


C        write(iout, *) 'End det C in smt'
C     Extract C (occupied-occupied block) from full matrix of orbital
C     coefficients (orb*).

      if ( inobas .and. itype .eq. 1 ) then
        call sptmat (2, norbs, norbs, v(jorb2), na, na, v(ic), 1)

      elseif ( inobas .and. itype .eq. 2 ) then
      call sptmat (2, norbs, norbs, v(jorb2), na, na, v(ica), 1)

        if ( nb .ge. 1 ) then
      call sptmat (2, norbs, norbs, v(jorb2+2*nosq), nb, nb, v(icb), 1)
        endif

      elseif ( inobas .and. itype .eq. 3 ) then
        call sptmat (2, 2*norbs, 2*norbs, v(jorb2), nel, nel, v(ic), 1)
      endif


C     Compute the quantity
C       detC2 = det ( C* )

      if ( inobas .and. itype .eq. 1 ) then

        if ( .not. dodiag ) then
          call zgetrf (na, na, v(ic), na, v(ipiv), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

          call detmat (na, v(ic), v(ipiv), fac1, fac2)

        else
          call zgeev ('n', 'n', na, v(ic), na, v(ieig), zero, 1, zero,
     $         1, v(iscr1), lscr1/2, v(iscr2), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

          call detmt2 (na, v(ieig), fac1, fac2)
        endif

        detC2 = cmplx (fac1) * exp (im * cmplx (fac2))
        detC2 = detC2**2

      elseif ( inobas .and. itype .eq. 2 ) then

        if ( .not. dodiag ) then
          call zgetrf (na, na, v(ica), na, v(ipiv), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

          call detmat (na, v(ica), v(ipiv), fac1, fac2)

          if ( nb .ge. 1 ) then
          call zgetrf (nb, nb, v(icb), nb, v(ipiv), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

          call detmat (nb, v(icb), v(ipiv), fac3, fac4)
          else
            fac3 = 1.0d0
            fac4 = 0.0d0
          endif

        else
          call zgeev ('n', 'n', na, v(ica), na, v(ieig), zero, 1, zero,
     $         1, v(iscr1), lscr1/2, v(iscr2), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

          call detmt2 (na, v(ieig), fac1, fac2)

          if ( nb .ge. 1 ) then
          call zgeev ('n', 'n', nb, v(icb), nb, v(ieig), zero, 1, zero,
     $         1, v(iscr1), lscr1/2, v(iscr2), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

          call detmt2 (nb, v(ieig), fac3, fac4)
          else
            fac3 = 1.0d0
            fac4 = 0.0d0
          endif
        endif

        detC2 = cmplx (fac1) * exp (im * cmplx (fac2))
        detC2 = detC2 * cmplx (fac3) * exp (im * cmplx (fac4))

      elseif ( inobas .and. itype .eq. 3 ) then

        if ( .not. dodiag ) then
          call zgetrf (nel, nel, v(ic), nel, v(ipiv), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

          call detmat (nel, v(ic), v(ipiv), fac1, fac2)

        else
          call zgeev ('n', 'n', nel, v(ic), nel, v(ieig), zero, 1, zero,
     $         1, v(iscr1), lscr1/2, v(iscr2), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in evalsmt')

          call detmt2 (nel, v(ieig), fac1, fac2)
        endif

        detC2 = cmplx (fac1) * exp (im * cmplx (fac2))
C        write(iout, 66) detC2
C 166  Format("DetC2 = ", F10.5)
      endif


      write(iout, *) 'End C in smt'

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Evaluate the overlap matrix elements  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Evaluate the overlap matrix elements as described above...


C     Clear iint values for all processors.

      do 301 ip = 1, npuse
        iint = jend + (ip - 1) * lenp

        call aclear (2*nci*nci, v(iint))
 301  continue


C     Loop over grid points.

      do 103 iblock = 1, ngrdt, npuse
        nthrd = min (npuse, ngrdt - iblock + 1)


C     Load matrices. Loop not parallelized.

        do 104 ip = 1, nthrd

          igrdt = iblock + ip - 1

C     Memory allocation per processor:
C       iint  - integrated quantity x(g)
C       ibufn - N(g) buffer at each grid point

          iint  = jend  + (ip - 1) * lenp
          ibufn = iint  + 2*nci*nci
          iend  = ibufn + szbfn

          ind = (igrdt - 1)*szbfn

          call fileio (2, -irwng, szbfn, v(ibufn), ind)
 104    continue


C     Computation loop. Parallelized.

        throk = throk .and. npuse .gt. 1

C$omp   parallel do if (throk) schedule(static,1) default(shared)
C$omp+  private (igrdt, iint, ibufn, ing, iend)
C$omp+  private (ipiv, ieig, iscr1, iscr2, inga, ingb)
C$omp+  private (igrda, igrdb, igrdg)
C$omp+  private (info, facw, facd, fac1, fac2, fac3, fac4, detN)
C$omp+  private (n1c, n2c, n1s, n2s, indy, n1, n2, indi)

        do 105 ip = 1, nthrd

          igrdt = iblock + ip - 1

C     Memory allocation per processor:
C       iint  - integrated quantity x(g)
C       ibufn - N(g) buffer at each grid point

          iint  = jend  + (ip - 1) * lenp
          ibufn = iint  + 2*nci*nci
          iend  = ibufn + szbfn

C       ipiv - scratch array with pivoting indices
C       ieig - scratch array with eigenvalues
C       iscr1, iscr2 - extra scratch space required by lapack

          if ( .not. dodiag ) then
            ipiv = iend
            iend = ipiv + lpiv

          else
            ieig  = iend
            iscr1 = ieig  + leig
            iscr2 = iscr1 + lscr1
            iend  = iscr2 + lscr2
          endif

C     Retrieve some useful factors.

          igrda = indarr(igrdt,1)
          igrdb = indarr(igrdt,2)
          igrdg = indarr(igrdt,3)


C     Compute useful factor
C       facw = wgt(g)

          facw = wgta(igrda) * wgtb(igrdb) * wgtg(igrdg)
C          write(iout, 866) wgta(igrda), wgtb(igrdb), wgtg(igrdg) 
866   Format(F10.5, F10.5, F10.5)


C     Loop over CI matrix elements (complex conjugation).

          do 101 n1c = 1, ncik
            do 102 n2c = 1, ncik

              ing = (n1c-1)*ncik*szn + (n2c-1)*szn + ibufn
              indy = (igrdt-1)*ncik*ncik + (n1c-1)*ncik + n2c

C     Useful indices for itype = 2.

              if ( itype .eq. 2 ) then
                inga = ing
                ingb = ing + 2*na*na
              endif

C     Compute the quantity
C       det ( N(g) )

              if ( itype .eq. 1 ) then

                if ( .not. dodiag ) then
                  call zgetrf (na, na, v(ing), na, v(ipiv), info)
                  if ( info .ne. 0 )
     $              call gauerr ('info .ne. 0 in evalsmt')

                  call detmat (na, v(ing), v(ipiv), fac1, fac2)

                else
                  call zgeev ('n', 'n', na, v(ing), na, v(ieig), zero,
     $                 1, zero, 1, v(iscr1), lscr1/2, v(iscr2), info)
                  if ( info .ne. 0 )
     $              call gauerr ('info .ne. 0 in evalsmt')

                  call detmt2 (na, v(ieig), fac1, fac2)
                endif

                detN = cmplx (fac1) * exp (im * cmplx (fac2))
                detN = detN**2

              elseif ( itype .eq. 2 ) then

                if ( .not. dodiag ) then
                  call zgetrf (na, na, v(inga), na, v(ipiv), info)
                  if ( info .ne. 0 )
     $              call gauerr ('info .ne. 0 in evalsmt')

                  call detmat (na, v(inga), v(ipiv), fac1, fac2)

                  if ( nb .ge. 1 ) then
                  call zgetrf (nb, nb, v(ingb), nb, v(ipiv), info)
                  if ( info .ne. 0 )
     $              call gauerr ('info .ne. 0 in evalsmt')

                  call detmat (nb, v(ingb), v(ipiv), fac3, fac4)
                  else
                    fac3 = 1.0d0
                    fac4 = 0.0d0
                  endif

                else
                  call zgeev ('n', 'n', na, v(inga), na, v(ieig), zero,
     $                 1, zero, 1, v(iscr1), lscr1/2, v(iscr2), info)
                  if ( info .ne. 0 )
     $              call gauerr ('info .ne. 0 in evalsmt')

                  call detmt2 (na, v(ieig), fac1, fac2)

                  if ( nb .ge. 1 ) then
                  call zgeev ('n', 'n', nb, v(ingb), nb, v(ieig), zero,
     $                 1, zero, 1, v(iscr1), lscr1/2, v(iscr2), info)
                  if ( info .ne. 0 )
     $              call gauerr ('info .ne. 0 in evalsmt')

                  call detmt2 (nb, v(ieig), fac3, fac4)
                  else
                    fac3 = 1.0d0
                    fac4 = 0.0d0
                  endif
                endif

                detN = cmplx (fac1) * exp (im * cmplx (fac2))
                detN = detN * cmplx (fac3) * exp (im * cmplx (fac4))

              elseif ( itype .eq. 3 ) then

                if ( .not. dodiag ) then
                  call zgetrf (nel, nel, v(ing), nel, v(ipiv), info)
                  if ( info .ne. 0 )
     $              call gauerr ('info .ne. 0 in evalsmt')

                  call detmat (nel, v(ing), v(ipiv), fac1, fac2)

                else
                  call zgeev ('n', 'n', nel, v(ing), nel, v(ieig), zero,
     $                 1, zero, 1, v(iscr1), lscr1/2, v(iscr2), info)
                  if ( info .ne. 0 )
     $              call gauerr ('info .ne. 0 in evalsmt')

                  call detmt2 (nel, v(ieig), fac1, fac2)
                endif

                detN = cmplx (fac1) * exp (im * cmplx (fac2))
                write(iout, 266) detN
 266   Format("DetN = ", F10.5, " +", F10.5, "i")
              endif


C     Compute factor x(g).

              if ( .not. inobas ) then
                if ( n1c .eq. 1 .and. n2c .eq. 1 ) then
                  yg(indy) = 1.0d0/ (detN * detC * conjg (detC))

                elseif ( n1c .eq. 1 .and. n2c .eq. 2 ) then
                  yg(indy) = 1.0d0/ (detN * detC * detC)

                elseif ( n1c .eq. 2 .and. n2c .eq. 1 ) then
                  yg(indy) = 1.0d0/ (detN * conjg (detC) * conjg (detC))

                elseif ( n1c .eq. 2 .and. n2c .eq. 2 ) then
                  yg(indy) = 1.0d0/ (detN * conjg (detC) * detC)
                endif

              elseif ( inobas ) then
                if ( n1c .eq. 1 .and. n2c .eq. 1 ) then
                  yg(indy) = 1.0d0/ (detN * detC * conjg (detC))

                elseif ( n1c .eq. 1 .and. n2c .eq. 2 ) then
                  yg(indy) = 1.0d0/ (detN * detC * conjg (detC2))

                elseif ( n1c .eq. 2 .and. n2c .eq. 1 ) then
                  yg(indy) = 1.0d0/ (detN * detC2 * conjg (detC))

                elseif ( n1c .eq. 2 .and. n2c .eq. 2 ) then
                  yg(indy) = 1.0d0/ (detN * detC2 * conjg (detC2))
                endif
              endif

C     Loop over CI matrix elements (spin projection).

              do 201 n1s = 1, ncis
                do 202 n2s = 1, ncis

                  n1 = (n1c-1)*ncis + n1s
                  n2 = (n2c-1)*ncis + n2s
                  indi = (n1-1)*nci*2 + (n2-1)*2 + iint

C     Compute Wigner's rotation matrix factor
C       facd = ( D^S_{MK} (Omega) )*

C       Note that
C         ( D^S_{MK} (Omega) )*  =  exp (+i * alpha * M) *
C                                   exp (+i * gamma * K) *
C                                   d^S_{MK} (beta),

C       where d^S_{MK} (beta) is Wigner's small d-matrix.

                  if ( ispghf ) then
                    facd = exp (+im * cmplx (grda(igrda) * dind(n1s))) *
     $                     exp (+im * cmplx (grdg(igrdg) * dind(n2s))) *
     $                     cmplx (dmt(igrdb,n1s,n2s))

                  elseif ( ispuhf ) then
                    facd = cmplx (dmt(igrdb,1,1))

                  else
                    facd = one
                  endif

                  v(indi)   = v(indi) + facw * dble (yg(indy) * facd)
                  v(indi+1) = v(indi+1) + facw * aimag (yg(indy) * facd)
C                  write(iout, 1066) facw, dble(yg(indy)), facd
1066  Format('facw = ', F10.5, ', yg = ', F10.5, ', facd = ', F10.5)

 202            continue
 201          continue

 102        continue
 101      continue
 105    continue

C$omp   end parallel do

 103  continue

      write(iout, 966) (yg(indy), indy=1,8)
966   Format('yg(indy) = ', 8(F9.6, ' +', F5.3, 'i '))

C     Accumulate all integrals.

      do 401 n1 = 1, nci
        do 402 n2 = 1, nci
          sums = zero

          do 302 ip = 1, npuse
            iint = jend + (ip - 1) * lenp
            indi = (n1-1)*nci*2 + (n2-1)*2 + iint
            sums = sums + cmplx (v(indi), v(indi+1))
 302      continue

C     Save overlap matrix elements.

          smt(n1,n2) = sums
C          write(iout, 'I5, I5, F10.5, F10.5') n1, n2, sums
 402    continue
 401  continue


      return
      end


