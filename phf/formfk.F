

      subroutine formfk (iout, iprint, v, lenv, iphf, nci, ncis, ncik,
     $     ngrdb, nbasis, norbs, na, nb, irwh, irwx, irwnos, ngrdt,
     $     indarr, grda, wgta, grdb, wgtb, grdg, wgtg, dmt, dind, civec,
     $     irwng, irwpg, irwgg, irwyg, intxg, yg, pno, pstrno, inobas,
     $     hsp, hph, f, irwf)

      implicit none


C +------------------------------------------------------------+
C |  formfk   --   CAJH, 07.2011                               |
C |                                                            |
C |                                                            |
C |  Build the effective Fock matrix associated with the       |
C |  projected Hartree-Fock method. This effective Fock        |
C |  matrix is obtained by taking derivatives of the energy    |
C |  expression with respect to elements of the HF density     |
C |  matrix:                                                   |
C |                                                            |
C |    d (E) / d(P_{kl})  ==  F_{lk},                          |
C |                                                            |
C |  where P stands for the density matrix and F stands for    |
C |  the effective Fock matrix.                                |
C |                                                            |
C |  The energy expression for PHF is given by                 |
C |                                                            |
C |    E = int dg y(g; i,j) *                                  |
C |               { Hsp (g; i,j) + Hph (g; i,j) },             |
C |                                                            |
C |  where there is a parametric dependence on the i and j     |
C |  indices in the functions y(g), Hsp(g), Hph(g). The        |
C |  functions Hsp (g) and Hph (g) are given by                |
C |                                                            |
C |     Hsp (g; i,j) = Tr [ H P(g; i,j) ],                     |
C |     Hph (g; i,j) = 1/2 * Tr [ G(g; i,j) P(g; i,j) ].       |
C |                                                            |
C |  Here, H is the core Hamiltonian, P(g) are transition      |
C |  density matrices, and G(g) are the matrices resulting     |
C |  from the contraction of P(g) with electron-repulsion      |
C |  integrals.                                                |
C |                                                            |
C |                                                            |
C |  The matrix F can be written explicitly as                 |
C |                                                            |
C |     F  =  int dg y(g) *                                    |
C |             { [ Hsp(g) + Hph(g) ] . Y(g) +  Theta(g) }     |
C |                                                            |
C |  Here, the Y(g) matrices are those formed in the           |
C |  subroutine formyg.F. The explicit expression for the      |
C |  matrix Theta(g) depends on the rotation matrices Dg       |
C |  associated with the spin projection operator as well      |
C |  as on the operators nu, mu associated with complex        |
C |  conjugation restoration. They are given by                |
C |                                                            |
C |                                                            |
C |  if nu = I, mu = I  ( <0| I! R(Omega) I |0> )              |
C |                                                            |
C |    Theta(g) =                                              |
C |      N(g) . ( P_oo  P_vo ) . F(g) . (I - P(g)) . Dg        |
C |                                                            |
C |  +   (I - P(g)) . F(g) . Dg . ( P_oo ) . N(g)              |
C |                               ( P_vo )                     |
C |                                                            |
C |                                                            |
C |  if nu = K, mu = I  ( <0| I! R(Omega) K |0> )              |
C |                                                            |
C |    Theta(g) =                                              |
C |    [ N(g) . ( P_oo  P_vo ) . F(g) . (I - P(g)) . Dg ]^T    |
C |                                                            |
C |  +   (I - P(g)) . F(g) . Dg . ( P_oo  P_vo )^T . N(g)      |
C |                                                            |
C |                                                            |
C |  if nu = I, mu = K  ( <0| K! R(Omega) I |0> )              |
C |                                                            |
C |    Theta(g) =                                              |
C |      N(g) . ( P_oo )^T . F(g) . (I - P(g)) . Dg            |
C |             ( P_vo )                                       |
C |                                                            |
C |  + [ (I - P(g)) . F(g) . Dg . ( P_oo ) . N(g) ]^T          |
C |    [                          ( P_vo )        ]            |
C |                                                            |
C |                                                            |
C |  if nu = K, mu = K  ( <0| K! R(Omega) K |0> )              |
C |                                                            |
C |    Theta(g) =                                              |
C |    [ N(g) . ( P_oo )^T . F(g) . (I - P(g)) . Dg ]^T        |
C |    [        ( P_vo )                            ]          |
C |                                                            |
C |  + [ (I - P(g)) . F(g) . Dg . ( P_oo  P_vo )^T . N(g) ]^T  |
C |                                                            |
C |                                                            |
C |  where nu, mu = I, K are operators related to complex      |
C |  conjugation restoration. (See evalsmt for more details.)  |
C |  The matrices N(g) are those built in formng. We have      |
C |  also defined the Fock-type matrices F(g) such that        |
C |                                                            |
C |    F(g) = H + G(g),                                        |
C |                                                            |
C |  where H is the core Hamiltonian. Lastly, P is the         |
C |  density matrix associated with the HF state while P(g)    |
C |  are transition density matrices (see formpg for more      |
C |  details).                                                 |
C |                                                            |
C |  In the equations provided above, N(g) is of dimensions    |
C |  N x N. Dg, F(g), and P(g) are all of dimensions           |
C |  nbasis x nbasis. The matrices involving P are therefore   |
C |  rectangular matrices of the appropriate dimensions.       |
C |  Namely, the subscripts (o,v) denote the dimensions of     |
C |  the block.                                                |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  In this subroutine we compute the effective Fock matrix   |
C |  in the NO basis in full square form. At the end, we       |
C |  transform the matrix back to the orthonormal AO basis,    |
C |  Hermitize it, and store it in lower triangular form in    |
C |  the appropriate RWF file. In this way, everything is      |
C |  ready for fockdg to diagonalize the effective Fock        |
C |  matrix.                                                   |
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
C       nbasis - number of basis functions
C       norbs  - number of orbitals
C       na     - number of alpha electrons
C       nb     - number of beta electrons
C       irwh   - RWF file with core Hamiltonian matrix
C       irwx   - RWF file with transformation matrix (=S^(-1/2))
C       irwnos - RWF file with transformation matrix from orthonormal AO
C                to NO basis
C       ngrdt  - total number of grid points
C       indarr - mapping array igrdt => igrdX
C       grdX   - integration grid for each variable
C       wgtX   - array with quadrature weights for each variable
C       dmt    - Wigner's small d-matrix d^S_{MK} (beta)
C       dind   - vector with indices {M} for Wigner rotation matrices
C       civec  - vector of CI coefficients
C       irwng  - RWF file with matrices N(g)
C       irwpg  - RWF file with matrices P(g)
C       irwgg  - RWF file with matrices G(g)
C       irwyg  - RWF file with matrices Y(g)
C       intxg  - integrated matrix X(g)
C       yg     - vector of y(g) values
C       pno    - density matrix in NO basis
C       pstrno - density matrix (P*) in NO basis
C       inobas - if .true., force the use of the true NO basis
C       hsp    - single-particle energy Hsp  [ out ]
C       hph    - particle-hole energy Hph  [ out ]
C       f      - effective Fock matrix  [ out ]
C       irwf   - RWF file with effective Fock matrix  [ updated here ]

      real*8      grda(*), grdb(*), grdg(*)
      real*8      wgta(*), wgtb(*), wgtg(*)
      integer     iout, iprint, lenv, iphf, nci, ncis, ncik
      integer     ngrdb, na, nb
      integer     nbasis, norbs, ngrdt, indarr(ngrdt,*)
      integer     irwh, irwx, irwnos
      integer     irwng, irwpg, irwgg, irwyg, irwf
      logical     inobas

      real*8      v(*), dmt(ngrdb,ncis,ncis), dind(ncis)
      real*8      pno(*), pstrno(*), f(*), hsp, hph
      complex*16  yg(*), civec(*), intxg(*)

C     other variables

      integer     nel, ntt, ntto, ntt2o, nosq
      integer     itype, iopcl, nri, szn, szp, szbfn, szbfp, ind1, ind2
      integer     jh1, jh2, jx, jscr1, ju, jend
      integer     jnos, jnost, jhno, jhno2, jf, jf2, lscrdg
      integer     jpnt, jptn, jpant, jpatn, jpbnt, jpbtn
      integer     jp2nt, jp2tn, jp2ant, jp2atn, jp2bnt, jp2btn
      integer     jzaa, jzab, jzba, jzbb, jscr2, lscr2
      integer     lenp, npmem, npmax, npuse, ipar
      integer     n1, n2, n1c, n2c, n1s, n2s, ip, iblock, nthrd
      logical     icmplx, ispuhf, ispghf, throk
      complex*16  sumhsp, sumhph

      real*8      facw, ener
      integer     igrdt, ihsp, ihph, ifock, ix, iy, iend
      integer     ibufn, ibufp, ibufg, ibufy
      integer     ing, inga, ingb, ipg, igg, iyg
      integer     it1, it1a, it1b, it2, it2a, it2b
      integer     it3, it3a, it3b, it4, it4a, it4b
      integer     iz1f, iz2f, iscrdg, ix2, indy
      integer     igrda, igrdb, igrdg
      real*8      anga, angb, angg
      complex*16  facd, fac, trh, trg

C     functions

      integer     nproc, lappar
      complex*16  trcabc

C     constants

      complex*16  zero, one, two, im


      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)
      two  = cmplx (2.0d0, 0.0d0)
      im   = cmplx (0.0d0, 1.0d0)


C     Construct useful variables.

      nel = na + nb

      ntt = nbasis * (nbasis + 1) / 2
      ntto = norbs * (norbs + 1) / 2
      ntt2o = norbs * (2*norbs + 1)

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
        call gauerr ('Unrecognized iphf in formfk.')
      endif


C     Read iopcl from ILSW file.

C       iopcl = 0,  real RHF
C             = 1,  real UHF
C             = 2,  complex RHF
C             = 3,  complex UHF
C             = 6,  complex GHF (there is no real GHF)

      call ilsw (2, 1, iopcl)

      if ( iopcl .gt. 3 .and. iopcl .ne. 6 ) then
        call gauerr ('Incorrect iopcl in formfk.')
      endif

      nri = 2

      if ( iopcl .le. 1 ) nri = 1


C     Define useful array dimensions:
C       - szn,  size of N(g) matrices
C       - szp,  size of P(g) matrices

      if ( itype .eq. 1 ) then
        szn = 2*na*na
        szp = 2*nosq

      elseif ( itype .eq. 2 ) then
        szn = 2*na*na + 2*nb*nb
        szp = 4*nosq

      elseif ( itype .eq. 3 ) then
        szn = 2*nel*nel
        szp = 8*nosq
      endif

C       .. buffer to load computed Ng, Pg, Yg, Gg ..
C       .. buffer to save computed Yg ..

      szbfn = ncik*ncik * szn
      szbfp = ncik*ncik * szp


C     Memory allocation.

C       Allocate space for:
C         - core Hamiltonian matrix (AO basis) (jh1)
C         - core Hamiltonian matrix (orthonormal basis) (jh2)
C         - AO => orthonormal AO transformation matrix (jx)
C         - scratch array to transform core Hamiltonian (jscr1)

      jh1   = 1
      jh2   = jh1   + nri*nbasis*nbasis
      jx    = jh2   + 2*nosq
      jscr1 = jx    + nri*nbasis*norbs
      jend  = jscr1 + nri*nbasis*norbs

C       Allocate space for:
C         - orthonormal AO => NO basis transformation matrix (jnos)
C         - scratch array to perform basis transformation (jnost)
C         - core Hamiltonian matrix in NO basis (jhno, jhno2)

      if ( itype .eq. 1 ) then
        jnos  = jend
        jnost = jnos  + 2*nosq
        jhno  = jnost + 2*nosq
        jhno2 = jhno  + 2*nosq
        jend  = jhno2 + 2*nosq

      elseif ( itype .eq. 2 ) then
        jnos  = jend
        jnost = jnos  + 4*nosq
        jhno  = jnost + 2*nosq
        jhno2 = jhno  + 4*nosq
        jend  = jhno2 + 4*nosq

      elseif ( itype .eq. 3 ) then
        jnos  = jend
        jnost = jnos  + 8*nosq
        jhno  = jnost + 8*nosq
        jhno2 = jhno  + 8*nosq
        jend  = jhno2 + 8*nosq
      endif

C       Allocate space for:
C         - rectangular blocks of density matrix in NO basis (jpnt, jptn)
C         - integrated Fock matrix (square and LT) (jf, jf2)
C         - rectangular blocks of density matrix (P*) in NO basis (jp2*)

      if ( itype .eq. 1 ) then
        jpnt = jend
        jptn = jpnt + 2*na*norbs
        jf   = jptn + 2*na*norbs
        jend = jf   + 2*nosq

      elseif ( itype .eq. 2 ) then
        jpant = jend
        jpatn = jpant + 2*na*norbs
        jpbnt = jpatn + 2*na*norbs
        jpbtn = jpbnt + 2*nb*norbs
        jf    = jpbtn + 2*nb*norbs
        jend  = jf    + 4*nosq

      elseif ( itype .eq. 3 ) then
        jpnt = jend
        jptn = jpnt + 4*nel*norbs
        jf   = jptn + 4*nel*norbs
        jend = jf   + 8*nosq
      endif

      if ( itype .lt. 3 .or. ispuhf ) then
        jf2  = jend
        jend = jf2  + 2*nosq

      elseif ( itype .eq. 3 ) then
        jf2  = jend
        jend = jf2  + 8*nosq
      endif

      if ( inobas .and. itype .eq. 1 ) then
        jp2nt = jend
        jp2tn = jp2nt + 2*na*norbs
        jend  = jp2tn + 2*na*norbs

      elseif ( inobas .and. itype .eq. 2 ) then
        jp2ant = jend
        jp2atn = jp2ant + 2*na*norbs
        jp2bnt = jp2atn + 2*na*norbs
        jp2btn = jp2bnt + 2*nb*norbs
        jend   = jp2btn + 2*nb*norbs

      elseif ( inobas .and. itype .eq. 3 ) then
        jp2nt = jend
        jp2tn = jp2nt + 4*nel*norbs
        jend  = jp2tn + 4*nel*norbs
      endif

C       Allocate space for:
C         - transpose(nos) . nos product (ju)

      if ( inobas .and. itype .eq. 1 ) then
        ju   = jend
        jend = ju   + 2*nosq

      elseif ( inobas .and. itype .eq. 2 ) then
        ju   = jend
        jend = ju   + 4*nosq

      elseif ( inobas .and. itype .eq. 3 ) then
        ju   = jend
        jend = ju   + 8*nosq
      endif

C       Allocate some extra scratch arrays for itype = 3.

      if ( itype .eq. 3 ) then
        lscr2 = 16*nosq

        jzaa  = jend
        jzab  = jzaa  + 2*nosq
        jzba  = jzab  + 2*nosq
        jzbb  = jzba  + 2*nosq
        jscr2 = jzbb  + 2*nosq
        jend  = jscr2 + lscr2
      endif


C     Do some preparation for parallel processing:
C       - lenp,  dimension of scratch space per processor
C       - npuse, number of processors to use

      lscrdg = 0

      if ( ispuhf .or. ispghf ) then
        lscrdg = 32*nosq

        lenp = 4 + 3*szp + lscrdg + 8*nel*norbs + 32*nosq

      elseif ( itype .eq. 1 ) then
        lenp = 4 + 3*szp + 4*na*norbs + 4*nosq

      elseif ( itype .eq. 2 ) then
        lenp = 4 + 3*szp + 4*na*norbs + 4*nb*norbs + 8*nosq

      elseif ( itype .eq. 3 ) then
        lenp = 4 + 3*szp + 8*nel*norbs + 16*nosq
      endif

C       .. length of Ng buffer, Pg buffer ..

      lenp = lenp + szbfn + 3*szbfp

      if ( inobas ) then
        if ( itype .lt. 3 ) lenp = lenp + 2*nosq
        if ( itype .eq. 3 ) lenp = lenp + 8*nosq
      endif


      npmem = (lenv - jend + 1) / lenp

      if ( npmem .lt. 1 ) then
        call gauerr ('Not enough memory in formfk for a single CPU.')
      endif

      npmax = nproc (0)
      npuse = min (npmem, npmax)
      
      ipar  = lappar (0)
      throk = ipar .eq. 1

      if ( ipar .eq. -1 ) npuse = 1

      call tstcor (jend + lenp*npuse - 1, lenv, 'formfk')



C     Load core Hamiltonian matrix.
C     Transform into a square matrix.

      call fileio (2, -irwh, nri*ntt, v(jh1), 0)

      if ( nri .eq. 1 ) then
        call square (v(jh1), v(jh1), nbasis, nbasis, 0)
      elseif ( nri .eq. 2 ) then
        call zsquare (v(jh1), v(jh1), nbasis, nbasis, 0)
      endif


C     Load transformation matrix.

      call fileio (2, -irwx, nri*nbasis*norbs, v(jx), 0)


C     Load transformation matrix (orthonormal AO => NO basis).

      if ( itype .eq. 1 ) then
        call fileio (2, -irwnos, 2*nosq, v(jnos), 0)

      elseif ( itype .eq. 2 ) then
        call fileio (2, -irwnos, 4*nosq, v(jnos), 0)

      elseif ( itype .eq. 3 ) then
        call fileio (2, -irwnos, 8*nosq, v(jnos), 0)
      endif



C     Transform core Hamiltonian matrix into orthonormal AO basis by
C     computing the product X! . H . X.

      if ( nri .eq. 1 ) then
        call dgemm ('t', 'n', norbs, nbasis, nbasis, 1.0d0, v(jx),
     $       nbasis, v(jh1), nbasis, 0.0d0, v(jscr1), norbs)
        call dgemm ('n', 'n', norbs, norbs, nbasis, 1.0d0, v(jscr1),
     $       norbs, v(jx), nbasis, 0.0d0, v(jh1), norbs)

C       .. transform core Hamiltonian matrix into complex ..

        call zrecmp (3, v(jh2), v(jh1), nosq)

      elseif ( nri .eq. 2 ) then
        call zgemm ('c', 'n', norbs, nbasis, nbasis, one, v(jx),
     $       nbasis, v(jh1), nbasis, zero, v(jscr1), norbs)
        call zgemm ('n', 'n', norbs, norbs, nbasis, one, v(jscr1),
     $       norbs, v(jx), nbasis, zero, v(jh1), norbs)

        call amove (2*nosq, v(jh1), v(jh2))
      endif


C     debugging...
C       Print core Hamiltonian in orthonormal basis.

c$$$       write (iout, *) ' '
c$$$       call ltoutx (iout, 1, 2, 1, 0, 'Core Hamiltonian', norbs,
c$$$     $      norbs, v(jh2), v(jh2), 0)


C     Transform core Hamiltonian into NO basis:
C       nos! . H . nos

      if ( itype .eq. 1 ) then

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jh2), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos), norbs, zero, v(jhno), norbs)

      elseif ( itype .eq. 2 ) then

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jh2), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos), norbs, zero, v(jhno), norbs)

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos+2*nosq),
     $       norbs, v(jh2), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos+2*nosq), norbs, zero, v(jhno+2*nosq), norbs)

      elseif ( itype .eq. 3 ) then

C       .. form core Hamiltonian in spin orbital basis ..

        call aclear (2*nosq, v(jzab))
        call aclear (2*nosq, v(jzba))

        call dmblck (iout, iprint, v(jscr2), lscr2, norbs, v(jhno),
     $       v(jh2), v(jzab), v(jzba), v(jh2), 1, -1)

C       .. perform the transformation ..

        call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $       v(jnos), 2*norbs, v(jhno), 2*norbs, zero, v(jnost),
     $       2*norbs)
        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $       v(jnost), 2*norbs, v(jnos), 2*norbs, zero, v(jhno),
     $       2*norbs)
      endif


C     debugging...
C       Print core Hamiltonian in NO basis.

c$$$       if ( itype .eq. 1 ) then
c$$$         write (iout, *) ' '
c$$$         call ltoutx (iout, 1, 2, 1, 0, 'Core Hamiltonian (NO)', norbs,
c$$$     $        norbs, v(jhno), v(jhno), 0)
c$$$
c$$$       elseif ( itype .eq. 2 ) then
c$$$         write (iout, *) ' '
c$$$         call ltoutx (iout, 1, 3, 1, 0, 'Core Hamiltonian (NO)', norbs,
c$$$     $        norbs, v(jhno), v(jhno+2*nosq), 0)
c$$$
c$$$       elseif ( itype .eq. 3 ) then
         write (iout, *) ' '
         call ltoutx (iout, 1, 2, 1, 0, 'Core Hamiltonian (NO)',
     $        2*norbs, 2*norbs, v(jhno), v(jhno), 0)
c$$$       endif


C     Take complex conjugate of core Hamiltonian (to take traces
C     appropriately).

      if ( itype .eq. 1 ) then
        call amove (2*nosq, v(jhno), v(jhno2))
        call vconjg (nosq, v(jhno2))

      elseif ( itype .eq. 2 ) then
        call amove (2*nosq, v(jhno), v(jhno2))
        call amove (2*nosq, v(jhno+2*nosq), v(jhno2+2*nosq))

        call vconjg (nosq, v(jhno2))
        call vconjg (nosq, v(jhno2+2*nosq))

      elseif ( itype .eq. 3 ) then
        call amove (8*nosq, v(jhno), v(jhno2))
        call vconjg (4*nosq, v(jhno2))
      endif


C     Extract rectangular blocks of density matrix.

      if ( itype .eq. 1 ) then
        call sptmat (2, norbs, norbs, pno, na, norbs, v(jpnt), 1)
        call sptmat (2, norbs, norbs, pno, norbs, na, v(jptn), 1)

      elseif ( itype .eq. 2 ) then
        call sptmat (2, norbs, norbs, pno, na, norbs, v(jpant), 1)
        call sptmat (2, norbs, norbs, pno, norbs, na, v(jpatn), 1)

        if ( nb .ge. 1 ) then
        call sptmat (2, norbs, norbs, pno(1+2*nosq), nb, norbs,
     $       v(jpbnt), 1)
        call sptmat (2, norbs, norbs, pno(1+2*nosq), norbs, nb,
     $       v(jpbtn), 1)
        endif

      elseif ( itype .eq. 3 ) then
        call sptmat (2, 2*norbs, 2*norbs, pno, nel, 2*norbs, v(jpnt), 1)
        call sptmat (2, 2*norbs, 2*norbs, pno, 2*norbs, nel, v(jptn), 1)
      endif

      if ( inobas .and. itype .eq. 1 ) then
        call sptmat (2, norbs, norbs, pstrno, na, norbs, v(jp2nt), 1)
        call sptmat (2, norbs, norbs, pstrno, norbs, na, v(jp2tn), 1)

      elseif ( inobas .and. itype .eq. 2 ) then
        call sptmat (2, norbs, norbs, pstrno, na, norbs, v(jp2ant), 1)
        call sptmat (2, norbs, norbs, pstrno, norbs, na, v(jp2atn), 1)

        if ( nb .ge. 1 ) then
        call sptmat (2, norbs, norbs, pstrno(1+2*nosq), nb, norbs,
     $       v(jp2bnt), 1)
        call sptmat (2, norbs, norbs, pstrno(1+2*nosq), norbs, nb,
     $       v(jp2btn), 1)
        endif

      elseif ( inobas .and. itype .eq. 3 ) then
        call sptmat (2, 2*norbs, 2*norbs, pstrno, nel, 2*norbs,
     $       v(jp2nt), 1)
        call sptmat (2, 2*norbs, 2*norbs, pstrno, 2*norbs, nel,
     $       v(jp2tn), 1)
      endif


C     Compute product
C       U = transpose (nos) . nos

      if ( inobas .and. itype .eq. 1 ) then
        call zgemm ('t', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jnos), norbs, zero, v(ju), norbs)

      elseif ( inobas .and. itype .eq. 2 ) then
        call zgemm ('t', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jnos), norbs, zero, v(ju), norbs)
        call zgemm ('t', 'n', norbs, norbs, norbs, one, v(jnos+2*nosq),
     $       norbs, v(jnos+2*nosq), norbs, zero, v(ju+2*nosq), norbs)

      elseif ( inobas .and. itype .eq. 3 ) then
        call zgemm ('t', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jnos),
     $       2*norbs, v(jnos), 2*norbs, zero, v(ju), 2*norbs)
      endif



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Form effective Fock matrix  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Form effective Fock matrix as described above...
C     Compute the single-particle and the particle-hole energies.


C     Clear fock matrices and accumulated quantities for all processors.

      do 11 ip = 1, npuse
        ihsp  = jend  + (ip - 1) * lenp
        ihph  = ihsp  + 2
        ifock = ihph  + 2
        iend  = ifock + szp

        call aclear (4, v(ihsp))
        call aclear (szp, v(ifock))
 11   continue


C     Loop over grid points.

      do 103 iblock = 1, ngrdt, npuse
        nthrd = min (npuse, ngrdt - iblock + 1)


C     Load matrices. Loop not parallelized.

        do 104 ip = 1, nthrd

          igrdt = iblock + ip - 1

C     Memory allocation per processor:
C       ihsp  - integrated quantity - int dg y(g) Hsp(g)
C       ihph  - integrated quantity - int dg y(g) Hph(d)
C       ifock - integrated Fock matrix
C       ibufn - N(g) buffer at each grid point
C       ibufp - P(g) buffer at each grid point
C       ibufg - G(g) buffer at each grid point
C       ibufy - Y(g) buffer at each grid point

          ihsp  = jend  + (ip - 1) * lenp
          ihph  = ihsp  + 2
          ifock = ihph  + 2
          ibufn = ifock + szp
          ibufp = ibufn + szbfn
          ibufg = ibufp + szbfp
          ibufy = ibufg + szbfp
          iend  = ibufy + szbfp

          ind1 = (igrdt - 1)*szbfn
          ind2 = (igrdt - 1)*szbfp

          call fileio (2, -irwng, szbfn, v(ibufn), ind1)
          call fileio (2, -irwpg, szbfp, v(ibufp), ind2)
          call fileio (2, -irwgg, szbfp, v(ibufg), ind2)
          call fileio (2, -irwyg, szbfp, v(ibufy), ind2)
 104    continue


C     Computation loop. Parallelized.

        throk = throk .and. npuse .gt. 1

C$omp   parallel do if (throk) schedule(static,1) default(shared)
C$omp+  private (igrdt, ihsp, ihph, ifock, ix, iy, iend)
C$omp+  private (ibufn, ibufp, ibufg, ibufy)
C$omp+  private (ing, inga, ingb, ipg, igg, iyg)
C$omp+  private (it1, it1a, it1b, it2, it2a, it2b)
C$omp+  private (it3, it3a, it3b, it4, it4a, it4b)
C$omp+  private (iz1f, iz2f, iscrdg, ix2)
C$omp+  private (igrda, igrdb, igrdg)
C$omp+  private (anga, angb, angg, facw, facd, fac)
C$omp+  private (trh, trg)
C$omp+  private (n1c, n2c, n1s, n2s, indy, n1, n2)

        do 105 ip = 1, nthrd

          igrdt = iblock + ip - 1

C     Memory allocation per processor:
C       ihsp  - integrated quantity - int dg y(g) Hsp(g)
C       ihph  - integrated quantity - int dg y(g) Hph(g)
C       ifock - integrated Fock matrix
C       ibufn - N(g) buffer at each grid point
C       ibufp - P(g) buffer at each grid point
C       ibufg - G(g) buffer at each grid point
C       ibufy - Y(g) buffer at each grid point
C       ix,iy - scratch arrays used to build effective Fock

          ihsp  = jend  + (ip - 1) * lenp
          ihph  = ihsp  + 2
          ifock = ihph  + 2
          ibufn = ifock + szp
          ibufp = ibufn + szbfn
          ibufg = ibufp + szbfp
          ibufy = ibufg + szbfp
          ix    = ibufy + szbfp
          iy    = ix    + szp
          iend  = iy    + szp


C       it1, it2, it3, it4, iz1f, iz2f - scratch arrays

          if ( itype .eq. 1 ) then
            it1  = iend
            it2  = it1  + 2*nosq
            it3  = it2  + 2*nosq
            it4  = it3  + 2*na*norbs
            iend = it4  + 2*na*norbs

          elseif ( itype .eq. 2 ) then
            it1a = iend
            it1b = it1a + 2*nosq
            it2a = it1b + 2*nosq
            it2b = it2a + 2*nosq
            it3a = it2b + 2*nosq
            it3b = it3a + 2*na*norbs
            it4a = it3b + 2*nb*norbs
            it4b = it4a + 2*na*norbs
            iend = it4b + 2*nb*norbs

          elseif ( itype .eq. 3 ) then
            it1  = iend
            it2  = it1  + 8*nosq
            it3  = it2  + 8*nosq
            it4  = it3  + 4*nel*norbs
            iend = it4  + 4*nel*norbs
          endif

          if ( ispuhf .or. ispghf ) then
            iz1f = iend
            iz2f = iz1f + 8*nosq
            iend = iz2f + 8*nosq
          endif

C       ix2 - scratch array to fix basis if inobas = .true.

          if ( inobas .and. itype .lt. 3 ) then
            ix2  = iend
            iend = ix2 + 2*nosq

          elseif ( inobas .and. itype .eq. 3 ) then
            ix2  = iend
            iend = ix2 + 8*nosq
          endif

C       iscrdg - scratch array for multiplication of Dg matrices

          if ( ispuhf .or. ispghf ) then
            iscrdg = iend
            iend   = iscrdg + lscrdg
          endif

C     Retrieve some useful factors.

          igrda = indarr(igrdt,1)
          igrdb = indarr(igrdt,2)
          igrdg = indarr(igrdt,3)

          anga = grda(igrda)
          angb = grdb(igrdb)
          angg = grdg(igrdg)


C     Compute useful factor
C       facw = wgt(g)

          facw = wgta(igrda) * wgtb(igrdb) * wgtg(igrdg)


C     Transform Dg matrix to NO basis:
C       nos' . Dg . nos

C     This product is stored in the array iz2f.

          if ( ispuhf .or. ispghf ) then
            call amove (8*nosq, v(jnos), v(iz1f))

            call mltdg (iout, iprint, v(iscrdg), lscrdg, norbs,
     $           v(iz1f), anga, angb, angg, 1)

            call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $           v(jnos), 2*norbs, v(iz1f), 2*norbs, zero, v(iz2f),
     $           2*norbs)
          endif


C     Loop over CI matrix elements (complex conjugation).

          do 101 n1c = 1, ncik
            do 102 n2c = 1, ncik

              ing = (n1c-1)*ncik*szn + (n2c-1)*szn + ibufn
              ipg = (n1c-1)*ncik*szp + (n2c-1)*szp + ibufp
              igg = (n1c-1)*ncik*szp + (n2c-1)*szp + ibufg
              iyg = (n1c-1)*ncik*szp + (n2c-1)*szp + ibufy
              indy = (igrdt-1)*ncik*ncik + (n1c-1)*ncik + n2c


C     Useful indices for itype = 2.

              if ( itype .eq. 2 ) then
                inga = ing
                ingb = inga + 2*na*na
              endif


C     ------------------------------

C     Compute the trace
C       tr ( H . P(g) )
C     and add it to ihsp.

C       ** trcabc computes tr ( a . b* )

              if ( itype .eq. 1 ) then
                trh = trcabc (1, norbs, v(ipg), v(jhno2))
                trh = two * trh

              elseif ( itype .eq. 2 ) then
                trh = trcabc (1, norbs, v(ipg), v(jhno2))
                trh = trh + 
     $                trcabc (1, norbs, v(ipg+2*nosq), v(jhno2+2*nosq))

              elseif ( itype .eq. 3 ) then
                trh = trcabc (1, 2*norbs, v(ipg), v(jhno2))
              endif


C     Compute the trace
C       1/2 . tr ( G(g) . P(g) )
C     and add it to ihph.

              if ( itype .eq. 1 ) then
                call vconjg (nosq, v(igg))

                trg = trcabc (1, norbs, v(ipg), v(igg))

                call vconjg (nosq, v(igg))

              elseif ( itype .eq. 2 ) then
                call vconjg (nosq, v(igg))
                call vconjg (nosq, v(igg+2*nosq))

                trg = trcabc (1, norbs, v(ipg), v(igg))
                trg = trg + 
     $                trcabc (1, norbs, v(ipg+2*nosq), v(igg+2*nosq))
                trg = trg / two

                call vconjg (nosq, v(igg))
                call vconjg (nosq, v(igg+2*nosq))

              elseif ( itype .eq. 3 ) then
                call vconjg (4*nosq, v(igg))

                trg = trcabc (1, 2*norbs, v(ipg), v(igg))
                trg = trg / two

                call vconjg (4*nosq, v(igg))
              endif
              
              write(iout, '("trh", 2F10.5, "trg", 2F10.5)') trh, trg

C     ------------------------------

C     Form the quantity
C       ( I - P(g) )

C     and store it in ipg.

              if ( itype .eq. 1 ) then
                call aunitm (.false., 2, norbs, norbs, v(it1))
                call asub (2*nosq, v(it1), v(ipg), v(it1))

                call amove (2*nosq, v(it1), v(ipg))

              elseif ( itype .eq. 2 ) then
                call aunitm (.false., 2, norbs, norbs, v(it1a))
                call asub (2*nosq, v(it1a), v(ipg), v(it1a))

                call aunitm (.false., 2, norbs, norbs, v(it1b))
                call asub (2*nosq, v(it1b), v(ipg+2*nosq), v(it1b))

                call amove (2*nosq, v(it1a), v(ipg))
                call amove (2*nosq, v(it1b), v(ipg+2*nosq))

              elseif ( itype .eq. 3 ) then
                call aunitm (.false., 2, 2*norbs, 2*norbs, v(it1))
                call asub (8*nosq, v(it1), v(ipg), v(it1))

                call amove (8*nosq, v(it1), v(ipg))
              endif


C     Form the quantity
C       F(g) = H + G(g)

C     and store it in igg. It is sort of a gauge-dependent Fock matrix.

              if ( itype .eq. 1 ) then
                call aadd (2*nosq, v(igg), v(jhno), v(igg))

              elseif ( itype .eq. 2 ) then
                call aadd (2*nosq, v(igg), v(jhno), v(igg))
                call aadd (2*nosq, v(igg+2*nosq), v(jhno+2*nosq),
     $               v(igg+2*nosq))

              elseif ( itype .eq. 3 ) then
                call aadd (8*nosq, v(igg), v(jhno), v(igg))
              endif


C     ------------------------------

C     We are now ready to start building the fields...


C     First, we will focus on the block that goes like

C       N(g) . ( P_oo   P_ov ) . F(g) . (I - P(g)) . Dg

C     We start from right to left...


C     Compute the product
C       T1  =  ( I - P(g) ) . Dg

              if ( ispuhf .or. ispghf ) then
                call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $               v(ipg), 2*norbs, v(iz2f), 2*norbs, zero, v(it1),
     $               2*norbs)

              else
                if ( itype .eq. 1 ) then
                  call amove (2*nosq, v(ipg), v(it1))

                elseif ( itype .eq. 2 ) then
                  call amove (2*nosq, v(ipg), v(it1a))
                  call amove (2*nosq, v(ipg+2*nosq), v(it1b))

                elseif ( itype .eq. 3 ) then
                  call amove (8*nosq, v(ipg), v(it1))
                endif
              endif


C     Compute the product
C       T2  =  F(g) . T1

              if ( itype .eq. 1 ) then
                call zgemm ('n', 'n', norbs, norbs, norbs, one, v(igg),
     $               norbs, v(it1), norbs, zero, v(it2), norbs)

              elseif ( itype .eq. 2 ) then
                call zgemm ('n', 'n', norbs, norbs, norbs, one, v(igg),
     $               norbs, v(it1a), norbs, zero, v(it2a), norbs)

                call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $               v(igg+2*nosq), norbs, v(it1b), norbs, zero,
     $               v(it2b), norbs)

              elseif ( itype .eq. 3 ) then
                call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $               v(igg), 2*norbs, v(it1), 2*norbs, zero, v(it2),
     $               2*norbs)
              endif


C     Compute the product

C       ( T3_oo   T3_ov )  =  ( P_oo   P_ov ) . ( T2_oo   T2_ov )
C                                               ( T2_vo   T2_vv )

C     or

C       ( T3_oo   T3_ov )  =  ( P_oo )^T . ( T2_oo   T2_ov )
C                             ( P_vo )     ( T2_vo   T2_vv )

              if ( n1c .eq. 1 ) then
                if ( itype .eq. 1 ) then
                  call zgemm ('n', 'n', na, norbs, norbs, one, v(jpnt),
     $                 na, v(it2), norbs, zero, v(it3), na)

                elseif ( itype .eq. 2 ) then
                  call zgemm ('n', 'n', na, norbs, norbs, one, v(jpant),
     $                 na, v(it2a), norbs, zero, v(it3a), na)

                  if ( nb .ge. 1 ) then
                  call zgemm ('n', 'n', nb, norbs, norbs, one, v(jpbnt),
     $                 nb, v(it2b), norbs, zero, v(it3b), nb)
                  endif

                elseif ( itype .eq. 3 ) then
                  call zgemm ('n', 'n', nel, 2*norbs, 2*norbs, one,
     $                 v(jpnt), nel, v(it2), 2*norbs, zero, v(it3), nel)
                endif

              elseif ( n1c .eq. 2 .and. inobas ) then
                if ( itype .eq. 1 ) then
                  call zgemm ('n', 'n', na, norbs, norbs, one, v(jp2nt),
     $                 na, v(it2), norbs, zero, v(it3), na)

                elseif ( itype .eq. 2 ) then
                  call zgemm ('n', 'n', na, norbs, norbs, one,
     $                 v(jp2ant), na, v(it2a), norbs, zero, v(it3a), na)

                  if ( nb .ge. 1 ) then
                  call zgemm ('n', 'n', nb, norbs, norbs, one,
     $                 v(jp2bnt), nb, v(it2b), norbs, zero, v(it3b), nb)
                  endif

                elseif ( itype .eq. 3 ) then
                  call zgemm ('n', 'n', nel, 2*norbs, 2*norbs, one,
     $                 v(jp2nt), nel, v(it2), 2*norbs, zero, v(it3),
     $                 nel)
                endif

              elseif ( n1c .eq. 2 .and. .not. inobas ) then
                if ( itype .eq. 1 ) then
                  call zgemm ('t', 'n', na, norbs, norbs, one, v(jptn),
     $                 norbs, v(it2), norbs, zero, v(it3), na)

                elseif ( itype .eq. 2 ) then
                  call zgemm ('t', 'n', na, norbs, norbs, one, v(jpatn),
     $                 norbs, v(it2a), norbs, zero, v(it3a), na) 

                  if ( nb .ge. 1 ) then
                  call zgemm ('t', 'n', nb, norbs, norbs, one, v(jpbtn),
     $                 norbs, v(it2b), norbs, zero, v(it3b), nb)
                  endif

                elseif ( itype .eq. 3 ) then
                  call zgemm ('t', 'n', nel, 2*norbs, 2*norbs, one,
     $                 v(jptn), 2*norbs, v(it2), 2*norbs, zero, v(it3),
     $                 nel)
                endif
              endif


C     Compute the product

C       ( T4_oo   T4_ov )  =  N(g) . ( T3_oo   T3_ov )

              if ( itype .eq. 1 ) then
                call zgemm ('n', 'n', na, norbs, na, one, v(ing), na,
     $               v(it3), na, zero, v(it4), na)

              elseif ( itype .eq. 2 ) then
                call zgemm ('n', 'n', na, norbs, na, one, v(inga), na,
     $               v(it3a), na, zero, v(it4a), na)

                if ( nb .ge. 1 ) then
                call zgemm ('n', 'n', nb, norbs, nb, one, v(ingb), nb,
     $               v(it3b), nb, zero, v(it4b), nb)
                endif

              elseif ( itype .eq. 3 ) then
                call zgemm ('n', 'n', nel, 2*norbs, nel, one, v(ing),
     $               nel, v(it3), nel, zero, v(it4), nel)
              endif


C     Transform T4 to full norbs x norbs, so that we can later add it to
C     ifock. Note that for n2c = 2 we need to take the transpose.

              if ( n2c .eq. 1 ) then
                if ( itype .eq. 1 ) then
                  call sptmat (2, norbs, norbs, v(ix), na, norbs,
     $                 v(it4), 2)

                elseif ( itype .eq. 2 ) then
                  call sptmat (2, norbs, norbs, v(ix), na, norbs,
     $                 v(it4a), 2)

                  if ( nb .ge. 1 ) then
                  call sptmat (2, norbs, norbs, v(ix+2*nosq), nb, norbs,
     $                 v(it4b), 2)
                  else
                    call aclear (2*nosq, v(ix+2*nosq))
                  endif

                elseif ( itype .eq. 3 ) then
                  call sptmat (2, 2*norbs, 2*norbs, v(ix), nel, 2*norbs,
     $                 v(it4), 2)
c          write (iout, *) ' '
c          call ltoutx (iout, 1, 2, 1, 0, 'T4', 2*norbs,
c     $         2*norbs, v(ix), v(ix), 0)
                endif

              elseif ( n2c .eq. 2 ) then
                if ( itype .eq. 1 ) then
                  call trprct (2, na, norbs, v(it4), v(it3))

                  call sptmat (2, norbs, norbs, v(ix), norbs, na,
     $                 v(it3), 2)

                elseif ( itype .eq. 2 ) then
                  call trprct (2, na, norbs, v(it4a), v(it3a))

                  if ( nb .ge. 1 ) then
                  call trprct (2, nb, norbs, v(it4b), v(it3b))
                  endif

                  call sptmat (2, norbs, norbs, v(ix), norbs, na,
     $                 v(it3a), 2)

                  if ( nb .ge. 1 ) then
                  call sptmat (2, norbs, norbs, v(ix+2*nosq), norbs, nb,
     $                 v(it3b), 2)
                  else
                    call aclear (2*nosq, v(ix+2*nosq))
                  endif

                elseif ( itype .eq. 3 ) then
                  call trprct (2, nel, 2*norbs, v(it4), v(it3))

                  call sptmat (2, 2*norbs, 2*norbs, v(ix), 2*norbs, nel,
     $                 v(it3), 2)
                endif
              endif

C              .. fix basis for inobas = .true. ..

              if ( n2c .eq. 2 .and. inobas ) then
                if ( itype .eq. 1 ) then
                  call zgemm ('c', 'n', norbs, norbs, norbs, one,
     $                 v(ju), norbs, v(ix), norbs, zero, v(ix2), norbs)
                  call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $                 v(ix2), norbs, v(ju), norbs, zero, v(ix), norbs)

                elseif ( itype .eq. 2 ) then
                  call zgemm ('c', 'n', norbs, norbs, norbs, one,
     $                 v(ju), norbs, v(ix), norbs, zero, v(ix2), norbs)
                  call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $                 v(ix2), norbs, v(ju), norbs, zero, v(ix), norbs)

                  call zgemm ('c', 'n', norbs, norbs, norbs, one,
     $                 v(ju+2*nosq), norbs, v(ix+2*nosq), norbs, zero,
     $                 v(ix2), norbs)
                  call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $                 v(ix2), norbs, v(ju+2*nosq), norbs, zero,
     $                 v(ix+2*nosq), norbs)

                elseif ( itype .eq. 3 ) then
                  call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $                 v(ju), 2*norbs, v(ix), 2*norbs, zero, v(ix2),
     $                 2*norbs)
                  call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $                 v(ix2), 2*norbs, v(ju), 2*norbs, zero, v(ix),
     $                 2*norbs)
                endif
              endif

C     Save current contribution in iy.

              call amove (szp, v(ix), v(iy))


C     ------------------------------

C     Now, we will focus on the block that goes like

C       (I - P(g)) . F(g) . Dg . ( P_oo ) . N(g)
C                                ( P_vo )

C     We start from left to right...


C     Compute the product
C       T1  =  ( I - P(g) ) . F(g)

              if ( itype .eq. 1 ) then
                call zgemm ('n', 'n', norbs, norbs, norbs, one, v(ipg),
     $               norbs, v(igg), norbs, zero, v(it1), norbs)

              elseif ( itype .eq. 2 ) then
                call zgemm ('n', 'n', norbs, norbs, norbs, one, v(ipg),
     $               norbs, v(igg), norbs, zero, v(it1a), norbs)

                call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $               v(ipg+2*nosq), norbs, v(igg+2*nosq), norbs, zero,
     $               v(it1b), norbs)

              elseif ( itype .eq. 3 ) then
                call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $               v(ipg), 2*norbs, v(igg), 2*norbs, zero, v(it1),
     $               2*norbs)
              endif


C     Compute the product
C       T2  =  T1 . Dg

              if ( ispuhf .or. ispghf ) then
                call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $               v(it1), 2*norbs, v(iz2f), 2*norbs, zero, v(it2),
     $               2*norbs)

              else
                if ( itype .eq. 1 ) then
                  call amove (2*nosq, v(it1), v(it2))

                elseif ( itype .eq. 2 ) then
                  call amove (2*nosq, v(it1a), v(it2a))
                  call amove (2*nosq, v(it1b), v(it2b))

                elseif ( itype .eq. 3 ) then
                  call amove (8*nosq, v(it1), v(it2))
                endif
              endif


C     Compute the product

C       ( T3_oo )  =  ( T2_oo   T2_ov ) . ( P_oo )
C       ( T3_vo )     ( T2_vo   T2_vv )   ( P_vo )

C     or

C       ( T3_oo )  =  ( T2_oo   T2_ov ) . ( P_oo   P_ov )^T
C       ( T3_vo )     ( T2_vo   T2_vv )

              if ( n2c .eq. 1 ) then
                if ( itype .eq. 1 ) then
                  call zgemm ('n', 'n', norbs, na, norbs, one, v(it2),
     $                 norbs, v(jptn), norbs, zero, v(it3), norbs)

                elseif ( itype .eq. 2 ) then
                  call zgemm ('n', 'n', norbs, na, norbs, one, v(it2a),
     $                 norbs, v(jpatn), norbs, zero, v(it3a), norbs)

                  if ( nb .ge. 1 ) then
                  call zgemm ('n', 'n', norbs, nb, norbs, one, v(it2b),
     $                 norbs, v(jpbtn), norbs, zero, v(it3b), norbs)
                  endif

                elseif ( itype .eq. 3 ) then
                  call zgemm ('n', 'n', 2*norbs, nel, 2*norbs, one,
     $                 v(it2), 2*norbs, v(jptn), 2*norbs, zero, v(it3),
     $                 2*norbs)
                endif

              elseif ( n2c .eq. 2 .and. inobas ) then
                if ( itype .eq. 1 ) then
                  call zgemm ('n', 'n', norbs, na, norbs, one, v(it2),
     $                 norbs, v(jp2tn), norbs, zero, v(it3), norbs)

                elseif ( itype .eq. 2 ) then
                  call zgemm ('n', 'n', norbs, na, norbs, one, v(it2a),
     $                 norbs, v(jp2atn), norbs, zero, v(it3a), norbs)

                  if ( nb .ge. 1 ) then
                  call zgemm ('n', 'n', norbs, nb, norbs, one, v(it2b),
     $                 norbs, v(jp2btn), norbs, zero, v(it3b), norbs)
                  endif

                elseif ( itype .eq. 3 ) then
                  call zgemm ('n', 'n', 2*norbs, nel, 2*norbs, one,
     $                 v(it2), 2*norbs, v(jp2tn), 2*norbs, zero, v(it3),
     $                 2*norbs)
                endif

              elseif ( n2c .eq. 2 .and. .not. inobas ) then
                if ( itype .eq. 1 ) then
                  call zgemm ('n', 't', norbs, na, norbs, one, v(it2),
     $                 norbs, v(jpnt), na, zero, v(it3), norbs)

                elseif ( itype .eq. 2 ) then
                  call zgemm ('n', 't', norbs, na, norbs, one, v(it2a),
     $                 norbs, v(jpant), na, zero, v(it3a), norbs)

                  if ( nb .ge. 1 ) then
                  call zgemm ('n', 't', norbs, nb, norbs, one, v(it2b),
     $                 norbs, v(jpbnt), nb, zero, v(it3b), norbs)
                  endif

                elseif ( itype .eq. 3 ) then
                  call zgemm ('n', 't', 2*norbs, nel, 2*norbs, one,
     $                 v(it2), 2*norbs, v(jpnt), nel, zero, v(it3),
     $                 2*norbs)
                endif
              endif


C     Compute the product

C       ( T4_oo )  =  ( T3_oo ) . N(g)
C       ( T4_vo )     ( T3_vo )

              if ( itype .eq. 1 ) then
                call zgemm ('n', 'n', norbs, na, na, one, v(it3), norbs,
     $               v(ing), na, zero, v(it4), norbs)

              elseif ( itype .eq. 2 ) then
                call zgemm ('n', 'n', norbs, na, na, one, v(it3a),
     $               norbs, v(inga), na, zero, v(it4a), norbs)

                if ( nb .ge. 1 ) then
                call zgemm ('n', 'n', norbs, nb, nb, one, v(it3b),
     $               norbs, v(ingb), nb, zero, v(it4b), norbs)
                endif

              elseif ( itype .eq. 3 ) then
                call zgemm ('n', 'n', 2*norbs, nel, nel, one, v(it3),
     $               2*norbs, v(ing), nel, zero, v(it4), 2*norbs)
              endif


C     Transform T4 to full norbs x norbs, so that we can later add it to
C     ifock. Note that for n1c = 2 we need to take the transpose.

              if ( n1c .eq. 1 ) then
                if ( itype .eq. 1 ) then
                  call sptmat (2, norbs, norbs, v(ix), norbs, na,
     $                 v(it4), 2)

                elseif ( itype .eq. 2 ) then
                  call sptmat (2, norbs, norbs, v(ix), norbs, na,
     $                 v(it4a), 2)

                  if ( nb .ge. 1 ) then
                  call sptmat (2, norbs, norbs, v(ix+2*nosq), norbs, nb,
     $                 v(it4b), 2)
                  else
                    call aclear (2*nosq, v(ix+2*nosq))
                  endif

                elseif ( itype .eq. 3 ) then
                  call sptmat (2, 2*norbs, 2*norbs, v(ix), 2*norbs, nel,
     $                 v(it4), 2)
c          write (iout, *) ' '
c          call ltoutx (iout, 1, 2, 1, 0, 'T4p', 2*norbs,
c    $         2*norbs, v(ix), v(ix), 0)
                endif

              elseif ( n1c .eq. 2 ) then
                if ( itype .eq. 1 ) then
                  call trprct (2, norbs, na, v(it4), v(it3))

                  call sptmat (2, norbs, norbs, v(ix), na, norbs,
     $                 v(it3), 2)

                elseif ( itype .eq. 2 ) then
                  call trprct (2, norbs, na, v(it4a), v(it3a))

                  if ( nb .ge. 1 ) then
                  call trprct (2, norbs, nb, v(it4b), v(it3b))
                  endif

                  call sptmat (2, norbs, norbs, v(ix), na, norbs,
     $                 v(it3a), 2)

                  if ( nb .ge. 1 ) then
                  call sptmat (2, norbs, norbs, v(ix+2*nosq), nb, norbs,
     $                 v(it3b), 2)
                  else
                    call aclear (2*nosq, v(ix+2*nosq))
                  endif

                elseif ( itype .eq. 3 ) then
                  call trprct (2, 2*norbs, nel, v(it4), v(it3))

                  call sptmat (2, 2*norbs, 2*norbs, v(ix), nel, 2*norbs,
     $                 v(it3), 2)
                endif
              endif

C              .. fix basis for inobas = .true. ..

              if ( n1c .eq. 2 .and. inobas ) then
                if ( itype .eq. 1 ) then
                  call zgemm ('c', 'n', norbs, norbs, norbs, one,
     $                 v(ju), norbs, v(ix), norbs, zero, v(ix2), norbs)
                  call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $                 v(ix2), norbs, v(ju), norbs, zero, v(ix), norbs)

                elseif ( itype .eq. 2 ) then
                  call zgemm ('c', 'n', norbs, norbs, norbs, one,
     $                 v(ju), norbs, v(ix), norbs, zero, v(ix2), norbs)
                  call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $                 v(ix2), norbs, v(ju), norbs, zero, v(ix), norbs)

                  call zgemm ('c', 'n', norbs, norbs, norbs, one,
     $                 v(ju+2*nosq), norbs, v(ix+2*nosq), norbs, zero,
     $                 v(ix2), norbs)
                  call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $                 v(ix2), norbs, v(ju+2*nosq), norbs, zero,
     $                 v(ix+2*nosq), norbs)

                elseif ( itype .eq. 3 ) then
                  call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $                 v(ju), 2*norbs, v(ix), 2*norbs, zero, v(ix2),
     $                 2*norbs)
                  call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $                 v(ix2), 2*norbs, v(ju), 2*norbs, zero, v(ix),
     $                 2*norbs)
                endif
              endif


C     Add contribution from iy.

              call aadd (szp, v(ix), v(iy), v(ix))


C     Add contribution of current grid point to the integrated quantity.

C     Loop over CI matrix elements (spin projection).

              do 201 n1s = 1, ncis
                do 202 n2s = 1, ncis

                  n1 = (n1c-1)*ncis + n1s
                  n2 = (n2c-1)*ncis + n2s

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

                  fac = yg(indy) * facd *
     $                  conjg(civec(n1)) * civec(n2)

C     Contributions to the energy.

                  v(ihsp) = v(ihsp) + facw * dble (fac * trh)
                  v(ihsp+1) = v(ihsp+1) + facw * aimag (fac * trh)

                  v(ihph) = v(ihph) + facw * dble (fac * trg)
                  v(ihph+1) = v(ihph+1) + facw * aimag (fac * trg)

C     Contributions to Fock matrix.

                  call amove (szp, v(ix), v(iy))
                  call zscal (szp/2, facw*fac, v(iy), 1)
                  call aadd (szp, v(ifock), v(iy), v(ifock))


C     Add contributions
C       wgt(g) . y(g) . [ Hsp(g) + Hph(g) ] . Y(g)
C     to ifock.

c          write (iout, *) ' '
c          call ltoutx (iout, 1, 2, 1, 0, 'iyg', 2*norbs,
c     $         2*norbs, v(iyg), v(iyg), 0)
                  call amove (szp, v(iyg), v(iy))
                  call zscal (szp/2, facw*fac*(trh + trg), v(iy), 1)
                  call aadd (szp, v(ifock), v(iy), v(ifock))

 202            continue
 201          continue
 102        continue
 101      continue

 105    continue

C$omp   end parallel do

 103  continue


C     Merge results from all processors.

      sumhsp = zero
      sumhph = zero

      call aclear (szp, v(jf))

      do 21 ip = 1, npuse
        ihsp  = jend  + (ip - 1) * lenp
        ihph  = ihsp  + 2
        ifock = ihph  + 2
        iend  = ifock + szp

        sumhsp = sumhsp + cmplx (v(ihsp), v(ihsp+1))
        sumhph = sumhph + cmplx (v(ihph), v(ihph+1))

        call aadd (szp, v(jf), v(ifock), v(jf))
 21   continue

      hsp = dble (sumhsp)
      hph = dble (sumhph)

      ener = hsp + hph


C     Add contribution from intxg.

      call zdscal (szp/2, -ener, intxg(1), 1)
      call aadd (szp, v(jf), intxg(1), v(jf), 1)



          write (iout, *) ' '
          call ltoutx (iout, 1, 2, 1, 0, 'Feff', 2*norbs,
     $         2*norbs, v(jf), v(jf), 0)
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Save effective Fock to RWF file  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Transform effective Fock matrix into orthonormal AO basis:
C       nos . F . nos!

      if ( itype .eq. 1 ) then

        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jf), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'c', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos), norbs, zero, v(jf), norbs)

      elseif ( itype .eq. 2 ) then

        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jf), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'c', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos), norbs, zero, v(jf), norbs)

        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnos+2*nosq),
     $       norbs, v(jf+2*nosq), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'c', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos+2*nosq), norbs, zero, v(jf+2*nosq), norbs)

      elseif ( itype .eq. 3 ) then

        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jnos),
     $       2*norbs, v(jf), 2*norbs, zero, v(jnost), 2*norbs)
        call zgemm ('n', 'c', 2*norbs, 2*norbs, 2*norbs, one, v(jnost),
     $       2*norbs, v(jnos), 2*norbs, zero, v(jf), 2*norbs)
      endif


C     Zero-out alpha-beta and beta-alpha spin blocks for ispuhf = .true.
C     Transform into UHF-type Fock matrix (only alpha-alpha and beta
C     -beta blocks are stored).

      if ( ispuhf ) then
        call dmblck (iout, iprint, v(jscr2), lscr2, norbs, v(jf),
     $       v(jzaa), v(jzab), v(jzba), v(jzbb), 2, -1)

        call amove (2*nosq, v(jzaa), v(jf))
        call amove (2*nosq, v(jzbb), v(jf+2*nosq))
      endif


C     Hermitize the effective Fock matrix and store as lower triangular.

      if ( itype .eq. 1 ) then
        call apadag (1, norbs, 0.50d0, v(jf), v(jf2))

        call trspnc (.false., norbs, norbs, v(jf2))
        call zlinear (v(jf2), f, norbs, norbs)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call apadag (1, norbs, 0.50d0, v(jf), v(jf2))

        call trspnc (.false., norbs, norbs, v(jf2))
        call zlinear (v(jf2), f, norbs, norbs)

        call apadag (1, norbs, 0.50d0, v(jf+2*nosq), v(jf2))

        call trspnc (.false., norbs, norbs, v(jf2))
        call zlinear (v(jf2), f(1+2*ntto), norbs, norbs)

      elseif ( itype .eq. 3 ) then
        call apadag (1, 2*norbs, 0.50d0, v(jf), v(jf2))

        call trspnc (.false., 2*norbs, 2*norbs, v(jf2))
        call zlinear (v(jf2), f, 2*norbs, 2*norbs)
      endif


C     Save matrices in RWF file.

      if ( itype .eq. 1 ) then
        call fileio (1, -irwf, 2*ntto, f, 0)
      elseif ( itype .eq. 2 .or. ispuhf ) then
        call fileio (1, -irwf, 4*ntto, f, 0)
      elseif ( itype .eq. 3 ) then
        call fileio (1, -irwf, 2*ntt2o, f, 0)
      endif


C     Print effective Fock matrices.

      if ( iprint .gt. 3 ) then
        if ( itype .eq. 1 ) then
          write (iout, *) ' '
          call ltoutx (iout, 0, 2, 1, 0, 'Effective Fock', norbs,
     $         norbs, f, f, 0)

        elseif ( itype .eq. 2 .or. ispuhf ) then
          write (iout, *) ' '
          call ltoutx (iout, 0, 3, 1, 0, 'Effective Fock', norbs,
     $         norbs, f, f(1+2*ntto), 0)

        elseif ( itype .eq. 3 ) then
          write (iout, *) ' '
          call ltoutx (iout, 0, 2, 1, 0, 'Effective Fock', 2*norbs,
     $         2*norbs, f, f, 0)
        endif
      endif


      return
      end


