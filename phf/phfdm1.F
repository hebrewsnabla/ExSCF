

      subroutine phfdm1 (iout, iprint, v, lenv, iphf, nci, ncis, ncik,
     $     nbasis, norbs, na, nb, s, isint, irwh, irwx, irwvec, irwnos,
     $     ngrdt, indarr, grda, wgta, grdb, wgtb, grdg, wgtg, irw1dm,
     $     irw2dm, pno, pstrno, civec, inobas, savedm, prtdm, dotwom,
     $     useao)

      implicit none


C +------------------------------------------------------------+
C |  phfdm1   ---   CAJH, 08.2011                              |
C |                                                            |
C |                                                            |
C |  Build the one- and the two-particle density matrices      |
C |  associated with the PHF state using the transition        |
C |  density matrices P(g,g') built in bldpgg.                 |
C |                                                            |
C |  The one-particle density matrix, gamma, is built          |
C |  according to                                              |
C |                                                            |
C |    gamma_i^j  =  < PHF | a_i! a_j | PHF >                  |
C |               =  int dg dg'  y(g,g')  [P(g,g')]_i^j        |
C |                                                            |
C |  The two-particle density matrix, Gamma, is built          |
C |  according to                                              |
C |                                                            |
C |    Gamma_{ij}^{lk}  =  1/2 < | a_i! a_j! a_k a_l | >       |
C |                     =  int dg dg'  y(g,g') .               |
C |                        { [P(g,g')]_i^l [P(g,g')]_j^k -     |
C |                          [P(g,g')]_i^k [P(g,g')]_j^l }     |
C |                                                            |
C |  Such forms of the density matrices can be derived by      |
C |  application of Wick's theorem to the set of rotated       |
C |  states involved in the projection.                        |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  The logical arguments savedm, prtdm, dotwom, useao        |
C |  control what is actually evaluated in this subroutine.    |
C |  Please look at the description of each of them below.     |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  CAJH, 08.2011                                             |
C |                                                            |
C |  Use of the true NO basis has been enabled even for        |
C |  methods with complex conjugation.                         |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       nci    - number of basis states in CI problem
C       ncis   - number of basis states for spin projection
C       ncik   - number of basis states for complex conj projection
C       nbasis - number of basis functions
C       norbs  - number of orbitals
C       na     - number of alpha electrons
C       nb     - number of beta electrons
C       s      - quantum number s on which to project
C       isint  - .true. for integer s, .false. for half-integer s
C       irwh   - RWF file with core Hamiltonian matrix
C       irwx   - RWF file with transformation matrix (=S^(-1/2))
C       irwvec - RWF file with orbital coefficients
C                (orthonormal AO basis)
C       irwnos - RWF file with transformation matrix from orthonormal AO
C                to NO basis
C       ngrdt  - total number of grid points
C       indarr - mapping array igrdt => igrdX
C       grdX   - integration grid for each variable
C       wgtX   - array with quadrature weights for each variable
C       irw1dm - RWF file to store 1PDM [ updated here ]
C       irw2dm - RWF file to store 2PDM [ updated here ]
C       pno    - density matrix in NO basis
C       pstrno - density matrix (P*) in NO basis
C       civec  - vector of CI coefficients
C       inobas - if .true., force the use of the true NO basis
C       savedm - whether density matrices should be saved
C       prtdm  - whether to print density matrices to files
C       dotwom - whether the two-matrix needs to be computed
C       useao  - if true, then the regular AO basis is used to build the
C                density matrices

      real*8      v(*), pno(*), pstrno(*)
      real*8      grda(*), grdb(*), grdg(*)
      real*8      wgta(*), wgtb(*), wgtg(*)
      complex*16  civec(*)
      integer     iout, iprint, lenv, iphf, nci, ncis, ncik, na, nb, s
      integer     nbasis, norbs, ngrdt, indarr(ngrdt,*)
      integer     irwh, irwx, irwvec, irwnos, irw1dm, irw2dm

      logical     isint, inobas, savedm, prtdm, dotwom, useao


C     other variables

      real*8      junk
      real*8      fac1, fac2, fac3, fac4, facw
      real*8      sclfac, s0, s1, seq
      complex*16  detC, detC2, intxgg, fac, xgg
      integer     nel, isz, itype, itype1, key
      integer     ntt, ntto, ntt2o, nosq, nbsuse, nbsq, nbcb
      integer     szn, szp, szp1, sz1pdm, sz2pdm
      integer     jx, jorb, jnos, jnost
      integer     jp, jocc, jwrk1, jwrk2, jorb2
      integer     jpnt, jptn, jpant, jpatn, jpbnt, jpbtn
      integer     jp2nt, jp2tn, jp2ant, jp2atn, jp2bnt, jp2btn
      integer     lscrd, jzab, jscrd, jz
      integer     jend, mdv, ic, ica, icb, info
      integer     jend1, jdmno, j1pdm, j2pdm
      integer     lpiv, leig, lscr1, lscr2, lscr3, lscrdg
      integer     lenp, npmem, npmax, npuse, ipar
      integer     ip, iblock, nthrd
      integer     n1, n2, n1c, n2c, n1s, n2s, is1, is2
      integer     ind5, ind6, ind7, ind8, ind
      logical     icmplx, ispuhf, ispghf, throk

      real*8      anga, angb, angg
      real*8      cgf0, cgf1, cgf2
      real*8      cgf0_0, cgf1_0, cgf1_p1, cgf1_m1
      real*8      cgf2_0, cgf2_p1, cgf2_m1, cgf2_p2, cgf2_m2
      integer     igrdt, iint, idmno, i1pdm, i2pdm, ingg, ipgg
      integer     ix, iza, izb, izab, izba, izf
      integer     it1, it2, it1a, it1b, it2a, it2b
      integer     ipiv, ieig, iscr1, iscr2, iscr3
      integer     iz1f, iz2f, iscrdg, iend
      integer     j2w10, j2w20, j2w2p1, j2w2m1, j2w1c, j2w2c
      integer     j2s10, j2s20, j2s2p1, j2s2m1, j2s1c, j2s2c
      integer     j2pdm1, j2pdm2
      integer     j1aaaa, j1bbbb, j1abba, j1baab, j1abab, j1baba
      integer     j1aabb, j1bbaa, j1aaab, j1aaba, j1abaa, j1baaa
      integer     j1bbba, j1bbab, j1babb, j1abbb
      integer     j2t10, j2t20, j2t2p1, j2t2m1, j2t30, j2t3p1, j2t3m1
      integer     j2t40, j2t60, j2t6p1, j2t6m1, j2t6p2, j2t6m2
      integer     j2t50, j2t5p1, j2t5m1
      integer     j2t1c, j2t2c, j2t3c, j2t4c, j2t5c, j2t6c
      integer     igrda, igrdb, igrdg
      integer     i1, i2, i3, i4, indt, ind1, ind2, ind3, ind4
      complex*16  wgtf0, wgtfp1, wgtfm1, wgtfp2, wgtfm2

      integer     iopcl, nri
      integer     jh1, jh2, jhno, jdmtst, jscr4

      complex*16  trh, trp

C     functions

      real*8      clgord, wignerd
      integer     nproc, lappar, intowp
      complex*16  trcabc


C     dodiag - whether to perform diagonalization to compute
C              determinants

C     debug  - whether to enter the debugging part of the code

      logical     dodiag, debug
      save        dodiag, debug
      data        dodiag / .true. /, debug / .false. /

C     constants

      complex*16  zero, one, two, im


      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)
      two  = cmplx (2.0d0, 0.0d0)
      im   = cmplx (0.0d0, 1.0d0)

      junk = 0.0d0


C     Construct useful variables.

      nel = na + nb

      if ( isint ) then
        isz = (na - nb) / 2
      else
        isz = na - nb
      endif

      ntt   = nbasis * (nbasis + 1) / 2
      ntto  = norbs * (norbs + 1) / 2
      ntt2o = norbs * (2*norbs + 1)

      nosq = norbs * norbs

      if ( useao ) then
        nbsuse = nbasis
      else
        nbsuse = norbs
      endif

      nbsq = nbsuse * nbsuse
      nbcb = nbsq * nbsuse


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

      itype1 = itype
      if ( ispuhf .or. ispghf ) itype1 = 2


C     Error checking.

      if ( iphf .lt. 1 .or. iphf .gt. 7 ) then
        call gauerr ('Unrecognized iphf in phfdm1.')
      endif


C     Read iopcl from ILSW file.

C       iopcl = 0,  real RHF
C             = 1,  real UHF
C             = 2,  complex RHF
C             = 3,  complex UHF
C             = 6,  complex GHF (there is no real GHF)

      call ilsw (2, 1, iopcl)

      if ( iopcl .gt. 3 .and. iopcl .ne. 6 ) then
        call gauerr ('Incorrect iopcl in phfdm1.')
      endif

      nri = 2

      if ( iopcl .le. 1 ) nri = 1


C     key controls which elements need to be transformed in any given
C     call to the basis set transformation routine (aobstf). The   
C     documentation of this routine says:

C       key = 1,  all spin blocks are active
C           = 2,  [aa] and [bb] blocks are active
C           = 3,  [aa] block is active

C     Therefore, we adjust key accordingly.

      if ( itype .eq. 1 ) then
        key = 3
      elseif ( itype .eq. 2 ) then
        key = 2
      elseif ( itype .eq. 3 ) then
        key = 1
      endif



C     Define useful array dimensions:
C       - szn,  size of N(g,g') matrices
C       - szp,  size of P(g,g') matrices

      if ( itype .eq. 1 ) then
        szn    = 2*na*na
        szp    = 2*nosq
        sz1pdm = 2*nbsq
        sz2pdm = 4*nbsq*nbsq

      elseif ( itype .eq. 2 ) then
        szn    = 2*na*na + 2*nb*nb
        szp    = 4*nosq
        sz1pdm = 4*nbsq
        sz2pdm = 6*nbsq*nbsq

      elseif ( itype .eq. 3 ) then
        szn    = 2*nel*nel
        szp    = 8*nosq

        if ( ispuhf .or. ispghf ) then
          sz1pdm = 4*nbsq
          sz2pdm = 6*nbsq*nbsq
        else
          sz1pdm = 8*nbsq
          sz2pdm = 18*nbsq*nbsq
        end if
      endif

      szp1 = szp
      if ( ispuhf .or. ispghf ) szp1 = 4*nosq


C     Memory allocation.

C       Allocate space for:
C         - AO => orthonormal AO transformation matrix (jx)

      jx   = 1
      jend = jx + nri*nbasis*norbs


C       Allocate space for:
C         - matrix of orbital coefficients in orthonormal AO basis (jorb)
C         - orthonormal AO => NO basis transformation matrix (jnos)
C         - scratch array to perform basis transformation (jnost)
C         - density matrix in lower triangular form (jp)
C         - vector of natural orbital occupations (jocc)
C         - scratch arrays needed for diagonalization (jwrk1, jwrk2)

      if ( itype .eq. 1 ) then
        jorb  = jend
        jnos  = jorb  + 2*nosq
        jnost = jnos  + 2*nosq
        jend  = jnost + 2*nosq

      elseif ( itype .eq. 2 ) then
        jorb  = jend
        jnos  = jorb  + 4*nosq
        jnost = jnos  + 4*nosq
        jend  = jnost + 2*nosq

      elseif ( itype .eq. 3 ) then
        jorb  = jend
        jnos  = jorb  + 8*nosq
        jnost = jnos  + 8*nosq
        jend  = jnost + 8*nosq
      endif

      if ( itype1 .eq. 1 ) then
        jp    = jend
        jocc  = jp    + 2*ntto
        jwrk1 = jocc  + norbs
        jwrk2 = jwrk1 + 2*(2*norbs-1)
        jend  = jwrk2 + 3*norbs-2

      elseif ( itype1 .eq. 2 ) then
        jp    = jend
        jocc  = jp    + 4*ntto
        jwrk1 = jocc  + 2*norbs
        jwrk2 = jwrk1 + 2*(2*norbs-1)
        jend  = jwrk2 + 3*norbs-2

      elseif ( itype1 .eq. 3 ) then
        jp    = jend
        jocc  = jp    + 2*ntt2o
        jwrk1 = jocc  + 2*norbs
        jwrk2 = jwrk1 + 2*(4*norbs-1)
        jend  = jwrk2 + 6*norbs-2
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

C       Allocate space for:
C         - rectangular blocks of density matrix in NO basis (jpnt, jptn)
C         - rectangular blocks of density matrix (P*) in NO basis (jp2*)

      if ( itype .eq. 1 ) then
        jpnt = jend
        jptn = jpnt + 2*na*norbs
        jend = jptn + 2*na*norbs

      elseif ( itype .eq. 2 ) then
        jpant = jend
        jpatn = jpant + 2*na*norbs
        jpbnt = jpatn + 2*na*norbs
        jpbtn = jpbnt + 2*nb*norbs
        jend  = jpbtn + 2*nb*norbs

      elseif ( itype .eq. 3 ) then
        jpnt = jend
        jptn = jpnt + 4*nel*norbs
        jend = jptn + 4*nel*norbs
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

C     Allocate scratch space for:
C       - scratch arrays needed for ispuhf = .true. (see below)

      if ( ispuhf ) then
        lscrd = 16*nosq

        jzab  = jend
        jscrd = jzab  + 2*nosq
        jz    = jscrd + lscrd
        jend  = jz    + 8*nosq
      endif


      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'phfdm1')


C     Dimensioning of some arrays for parallel processing:

      lpiv = 0
      leig = 0
      lscr1 = 0
      lscr2 = 0
      lscr3 = 0
      lscrdg = 0

      if ( itype .eq. 1 .or. itype .eq. 2 ) then
        lpiv = intowp (na)

        if ( dodiag ) leig  = 2*na
        if ( dodiag ) lscr1 = max (2*na*na, 4*na)
        if ( dodiag ) lscr2 = 2*na

      elseif ( itype .eq. 3 ) then
        lpiv = intowp (nel)

        if ( dodiag ) leig  = 2*nel
        if ( dodiag ) lscr1 = max (2*nel*nel, 4*nel)
        if ( dodiag ) lscr2 = 2*nel
      endif

      if ( ispuhf .or. ispghf ) lscrdg = 32*nosq

      if ( ( useao .and. itype .eq. 3 )
     $     .or. ispuhf .or. ispghf ) then
        lscr3 = max (16*nbsq, 4*nbasis*norbs)
      elseif ( useao ) then
        lscr3 = 4*nbasis*norbs
      endif


C     Do some preparation for parallel processing:
C       - lenp,  dimension of scratch space per processor
C       - npuse, number of processors to use

      if ( itype .eq. 1 ) then
        lenp = 2 + szn + szp + szp1 + 2*nosq + 4*na*norbs

        if ( savedm ) lenp = lenp + 2*nbsq

      elseif ( itype .eq. 2 ) then
        lenp = 2 + szn + szp + szp1 + 2*nosq + 4*na*norbs + 4*nb*norbs

        if ( savedm ) lenp = lenp + 4*nbsq

      elseif ( itype .eq. 3 ) then
        lenp = 2 + szn + szp + szp1 + 8*nosq + 8*nel*norbs

        if ( savedm .or. ispuhf .or. ispghf ) lenp = lenp + 8*nbsq
        if ( useao .or. ispuhf .or. ispghf )  lenp = lenp + 8*nbsq
        if ( ispuhf .or. ispghf ) lenp = lenp + 16*nosq
      endif

      lenp = lenp + lscr1 + lscr2 + lscr3 + lpiv + leig + lscrdg

      if ( savedm ) lenp = lenp + sz1pdm
      if ( dotwom ) lenp = lenp + sz2pdm

      if ( ispuhf .or. ispghf ) then
        lenp = lenp + 12*nosq
        if ( savedm ) lenp = lenp + 12*nbsq
        if ( dotwom ) lenp = lenp + 58*nbsq*nbsq
      endif

      npmem = (lenv - jend + 1) / lenp

      if ( npmem .lt. 1 ) then
        call gauerr ('Not enough memory in phfdm1 for a single CPU.')
      endif

      npmax = nproc (0)
      npuse = min (npmem, npmax)
      
      ipar  = lappar (0)
      throk = ipar .eq. 1

      if ( ipar .eq. -1 ) npuse = 1

      call tstcor (jend + lenp*npuse - 1, lenv, 'phfdm1')



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Load transformation matrices  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Load orbital coefficients  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
        endif

      endif

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

C     The overlap matrix elements appearing in x(g,g') can then be
C     computed as

C       <phi| I! R(Omega')! R(Omega) I |phi>  =
C                                      1 / ( det ( C . N(g,g') . C! ),

C       <phi| I! R(Omega')! R(Omega) K |phi>  =
C                                      1 / ( det ( C . N(g,g') . C^T ),

C       <phi| K! R(Omega')! R(Omega) I |phi>  =
C                                      1 / ( det ( C* . N(g,g') . C! ),

C       <phi| K! R(Omega')! R(Omega) K |phi>  =
C                                      1 / ( det ( C* . N(g,g') . C^T ),

C     where the matrices N(g,g') are those given in bldngg.

C     We now proceed to evaluate the determinant of C...


C     Memory allocation.

      if ( .not. dodiag ) then
        ic   = jend
        ipiv = ic   + szn
        iend = ipiv + lpiv

      else
        ic    = jend
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


C     Compute the quantity
C       detC = det ( C )

      if ( itype .eq. 1 ) then

        if ( .not. dodiag ) then
          call zgetrf (na, na, v(ic), na, v(ipiv), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in phfdm1')

          call detmat (na, v(ic), v(ipiv), fac1, fac2)

        else
          call zgeev ('n', 'n', na, v(ic), na, v(ieig), zero, 1, zero,
     $         1, v(iscr1), lscr1/2, v(iscr2), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in phfdm1')

          call detmt2 (na, v(ieig), fac1, fac2)
        endif

        detC = cmplx (fac1) * exp (im * cmplx (fac2))
        detC = detC**2

      elseif ( itype .eq. 2 ) then

        if ( .not. dodiag ) then
          call zgetrf (na, na, v(ica), na, v(ipiv), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in phfdm1')

          call detmat (na, v(ica), v(ipiv), fac1, fac2)

          if ( nb .ge. 1 ) then
          call zgetrf (nb, nb, v(icb), nb, v(ipiv), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in phfdm1')

          call detmat (nb, v(icb), v(ipiv), fac3, fac4)
          else
            fac3 = 1.0d0
            fac4 = 0.0d0
          endif

        else
          call zgeev ('n', 'n', na, v(ica), na, v(ieig), zero, 1, zero,
     $         1, v(iscr1), lscr1/2, v(iscr2), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in phfdm1')

          call detmt2 (na, v(ieig), fac1, fac2)

          if ( nb .ge. 1 ) then
          call zgeev ('n', 'n', nb, v(icb), nb, v(ieig), zero, 1, zero,
     $         1, v(iscr1), lscr1/2, v(iscr2), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in phfdm1')

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
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in phfdm1')

          call detmat (nel, v(ic), v(ipiv), fac1, fac2)

        else
          call zgeev ('n', 'n', nel, v(ic), nel, v(ieig), zero, 1, zero,
     $         1, v(iscr1), lscr1/2, v(iscr2), info)
          if ( info .ne. 0 ) call gauerr ('info .ne. 0 in phfdm1')

          call detmt2 (nel, v(ieig), fac1, fac2)
        endif

        detC = cmplx (fac1) * exp (im * cmplx (fac2))
      endif


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
      endif




C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Prepare HF density matrix in NO basis  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Build the density matrices  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Build the one- and possibly the two-particle density matrix as
C     described above.


C     Clear up arrays for all processors.

      do 11 ip = 1, npuse
        iint  = jend  + (ip - 1) * lenp
        idmno = iint  + 2
        iend  = idmno + szp1

        v(iint) = 0.0d0
        v(iint+1) = 0.0d0

        call aclear (szp1, v(idmno))

        if ( savedm ) then
          i1pdm = iend
          iend  = i1pdm + sz1pdm

          call aclear (sz1pdm, v(i1pdm))
        endif

        if ( dotwom ) then
          i2pdm = iend
          iend  = i2pdm + sz2pdm

          call aclear (sz2pdm, v(i2pdm))
        endif
 11   continue


C     Some Clebsch-Gordan coefficients for spin projection.

      if ( ispuhf .or. ispghf ) then
        cgf2 = clgord (s, isz, isint, 2, 0, .true., s, isz, isint)
        cgf1 = clgord (s, isz, isint, 1, 0, .true., s, isz, isint)
        cgf0 = clgord (s, isz, isint, 0, 0, .true., s, isz, isint)
      end if


C     Loop over grid points.

      do 103 iblock = 1, ngrdt, npuse
        nthrd = min (npuse, ngrdt - iblock + 1)


C     Computation loop. Parallelized.

        throk = throk .and. npuse .gt. 1

C$omp   parallel do if (throk) schedule(static,1) default(shared)
C$omp+  private (igrdt, iint, idmno, i1pdm, i2pdm, ingg, ipgg)
C$omp+  private (ix, iza, izb, izab, izba, izf)
C$omp+  private (it1, it2, it1a, it1b, it2a, it2b)
C$omp+  private (ipiv, ieig, iscr1, iscr2, iscr3)
C$omp+  private (iz1f, iz2f, iscrdg, iend)
C$omp+  private (igrda, igrdb, igrdg)
C$omp+  private (anga, angb, angg)
C$omp+  private (i1, i2, i3, i4, ind1, ind2, ind3, ind4, indt)
C$omp+  private (fac, facw, xgg)
C$omp+  private (cgf0_0, cgf1_0, cgf1_p1, cgf1_m1)
C$omp+  private (cgf2_0, cgf2_p1, cgf2_m1, cgf2_p2, cgf2_m2)
C$omp+  private (wgtf0, wgtfp1, wgtfm1, wgtfp2, wgtfm2)
C$omp+  private (n1c, n2c, n1s, n2s, n1, n2, is1, is2)
C$omp+  private (j2w10, j2w20, j2w2p1, j2w2m1, j2w1c, j2w2c)
C$omp+  private (j2s10, j2s20, j2s2p1, j2s2m1, j2s1c, j2s2c)
C$omp+  private (j2pdm1, j2pdm2)
C$omp+  private (j1aaaa, j1bbbb, j1abba, j1baab, j1abab, j1baba)
C$omp+  private (j1aabb, j1bbaa, j1aaab, j1aaba, j1abaa, j1baaa)
C$omp+  private (j1bbba, j1bbab, j1babb, j1abbb)
C$omp+  private (j2t10, j2t20, j2t2p1, j2t2m1)
C$omp+  private (j2t40, j2t60, j2t6p1, j2t6m1, j2t6p2, j2t6m2)
C$omp+  private (j2t30, j2t3p1, j2t3m1, j2t50, j2t5p1, j2t5m1)
C$omp+  private (j2t1c, j2t2c, j2t3c, j2t4c, j2t6c, j2t5c)

        do 104 ip = 1, nthrd

          igrdt = iblock + ip - 1

C     Memory allocation per processor:
C       iint  - integrated quantity  c_i* . c_j . x(g)
C       idmno - integrated density matrix (NO basis)
C       i1pdm - integrated 1PDM in user-specified basis
C       i2pdm - integrated 2PDM in user-specified basis
C       ingg  - N(g) matrix at each grid point
C       ipgg  - P(g) matrix at each grid point

          iint  = jend  + (ip - 1) * lenp
          idmno = iint  + 2
          iend  = idmno + szp1

          if ( savedm ) then
            i1pdm = iend
            iend  = i1pdm + sz1pdm
          endif

          if ( dotwom ) then
            i2pdm = iend
            iend  = i2pdm + sz2pdm
          endif

          ingg = iend
          ipgg = ingg + szn
          iend = ipgg + szp

C       ix, it1, it2 - scratch arrays

          if ( itype .eq. 1 ) then
            ix   = iend
            it1  = ix   + 2*nosq
            it2  = it1  + 2*na*norbs
            iend = it2  + 2*na*norbs

          elseif ( itype .eq. 2 ) then
            ix   = iend
            it1a = ix   + 2*nosq
            it1b = it1a + 2*na*norbs
            it2a = it1b + 2*nb*norbs
            it2b = it2a + 2*na*norbs
            iend = it2b + 2*nb*norbs

          elseif ( itype .eq. 3 ) then
            ix   = iend
            it1  = ix   + 8*nosq
            it2  = it1  + 4*nel*norbs
            iend = it2  + 4*nel*norbs
          endif

C       iz1f, iz2f, iz3f, iscrdg - scratch arrays to store Dg matrices

          if ( ispuhf .or. ispghf ) then
            iz1f   = iend
            iz2f   = iz1f   + 8*nosq
            iscrdg = iz2f   + 8*nosq
            iend   = iscrdg + lscrdg
          endif

C       iza, izb, izab, izba, izf - scratch arrays

          if ( savedm .and. itype .eq. 1 ) then
            iza  = iend
            iend = iza  + 2*nbsq

          elseif ( savedm .and. itype .eq. 2 ) then
            iza  = iend
            izb  = iza  + 2*nbsq
            iend = izb  + 2*nbsq

          elseif ( ( savedm .and. itype .eq. 3 )
     $             .or. ispuhf .or. ispghf ) then
            izf  = iend
            iend = izf  + 8*nbsq
          endif

          if ( ( useao .and. itype .eq. 3 )
     $         .or. ispuhf .or. ispghf ) then
            iza  = iend
            izb  = iza  + 2*nbsq
            izab = izb  + 2*nbsq
            izba = izab + 2*nbsq
            iend = izba + 2*nbsq
          endif

C       ipiv, leig, iscr1, iscr2, iscr3 - scratch arrays

          ipiv  = iend
          ieig  = ipiv  + lpiv
          iscr1 = ieig  + leig
          iscr2 = iscr1 + lscr1
          iscr3 = iscr2 + lscr2
          iend  = iscr3 + lscr3

          if ( ispuhf .or. ispghf ) then
            j2w10  = iend
            j2w20  = j2w10  + 2*nosq
            j2w2p1 = j2w20  + 2*nosq
            j2w2m1 = j2w2p1 + 2*nosq
            j2w1c  = j2w2m1 + 2*nosq
            j2w2c  = j2w1c  + 2*nosq
            iend   = j2w2c  + 2*nosq
          endif

          if ( savedm .and. (ispuhf .or. ispghf) ) then
            j2s10  = iend
            j2s20  = j2s10  + 2*nbsq
            j2s2p1 = j2s20  + 2*nbsq
            j2s2m1 = j2s2p1 + 2*nbsq
            j2s1c  = j2s2m1 + 2*nbsq
            j2s2c  = j2s1c  + 2*nbsq
            iend   = j2s2c  + 2*nbsq
          endif

          if ( dotwom .and. (ispuhf .or. ispghf) ) then
            j2pdm1 = iend
            j2pdm2 = j2pdm1 + 32*nbsq*nbsq
            iend   = j2pdm2 + 26*nbsq*nbsq
c$$$            iend   = j2pdm2 + 32*nbsq*nbsq

            j1aaaa = j2pdm1
            j1bbbb = j1aaaa + 2*nbsq*nbsq
            j1abba = j1bbbb + 2*nbsq*nbsq
            j1baab = j1abba + 2*nbsq*nbsq
            j1abab = j1baab + 2*nbsq*nbsq
            j1baba = j1abab + 2*nbsq*nbsq
            j1aabb = j1baba + 2*nbsq*nbsq
            j1bbaa = j1aabb + 2*nbsq*nbsq
            j1aaab = j1bbaa + 2*nbsq*nbsq
            j1aaba = j1aaab + 2*nbsq*nbsq
            j1abaa = j1aaba + 2*nbsq*nbsq
            j1baaa = j1abaa + 2*nbsq*nbsq
            j1bbba = j1baaa + 2*nbsq*nbsq
            j1bbab = j1bbba + 2*nbsq*nbsq
            j1babb = j1bbab + 2*nbsq*nbsq
            j1abbb = j1babb + 2*nbsq*nbsq

            j2t10  = j2pdm2
            j2t20  = j2t10  + 2*nbsq*nbsq
            j2t2p1 = j2t20  + 2*nbsq*nbsq
            j2t2m1 = j2t2p1 + 2*nbsq*nbsq
            j2t30  = j2t2m1 + 2*nbsq*nbsq
            j2t3p1 = j2t30  + 2*nbsq*nbsq
            j2t3m1 = j2t3p1 + 2*nbsq*nbsq
            j2t40  = j2t3m1 + 2*nbsq*nbsq
            j2t60  = j2t40  + 2*nbsq*nbsq
            j2t6p1 = j2t60  + 2*nbsq*nbsq
            j2t6m1 = j2t6p1 + 2*nbsq*nbsq
            j2t6p2 = j2t6m1 + 2*nbsq*nbsq
            j2t6m2 = j2t6p2 + 2*nbsq*nbsq
c$$$            j2t50  = j2t6m2 + 2*nbsq*nbsq
c$$$            j2t5p1 = j2t50  + 2*nbsq*nbsq
c$$$            j2t5m1 = j2t5p1 + 2*nbsq*nbsq

            j2t1c  = j2pdm1
            j2t2c  = j2t1c + 2*nbsq*nbsq
            j2t3c  = j2t2c + 2*nbsq*nbsq
            j2t4c  = j2t3c + 2*nbsq*nbsq
            j2t6c  = j2t4c + 2*nbsq*nbsq
c$$$            j2t5c  = j2t6c + 2*nbsq*nbsq
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

C     Evaluate matrix N(g).

              if ( ispuhf .or. ispghf ) then
                call bldngg (iout, iprint, itype, ispuhf, ispghf, n1c,
     $               n2c, na, nb, norbs, inobas, v(jpnt), v(jptn), junk,
     $               junk, v(jp2nt), v(jp2tn), junk, junk, v(it1), junk,
     $               v(it2), junk, v(iz2f), v(ipiv), v(ingg))

              elseif ( itype .eq. 1 .or. itype .eq. 3 ) then
                call bldngg (iout, iprint, itype, ispuhf, ispghf, n1c,
     $               n2c, na, nb, norbs, inobas, v(jpnt), v(jptn), junk,
     $               junk, v(jp2nt), v(jp2tn), junk, junk, v(it1), junk,
     $               v(it2), junk, junk, v(ipiv), v(ingg))

              elseif ( itype .eq. 2 ) then
                call bldngg (iout, iprint, itype, ispuhf, ispghf, n1c,
     $               n2c, na, nb, norbs, inobas, v(jpant), v(jpatn),
     $               v(jpbnt), v(jpbtn), v(jp2ant), v(jp2atn),
     $               v(jp2bnt), v(jp2btn), v(it1a), v(it1b), v(it2a),
     $               v(it2b), junk, v(ipiv), v(ingg))
              endif


C     Evaluate matrix P(g).

              if ( ispuhf .or. ispghf ) then
                call bldpgg (iout, iprint, itype, ispuhf, ispghf, n1c,
     $               n2c, na, nb, norbs, inobas, v(jpnt), v(jptn), junk,
     $               junk, v(jp2nt), v(jp2tn), junk, junk, v(it1), junk,
     $               v(it2), junk, v(iz2f), v(ingg), v(ipgg))

              elseif ( itype .eq. 1 .or. itype .eq. 3 ) then
                call bldpgg (iout, iprint, itype, ispuhf, ispghf, n1c,
     $               n2c, na, nb, norbs, inobas, v(jpnt), v(jptn), junk,
     $               junk, v(jp2nt), v(jp2tn), junk, junk, v(it1), junk,
     $               v(it2), junk, junk, v(ingg), v(ipgg))

              elseif ( itype .eq. 2 ) then
                call bldpgg (iout, iprint, itype, ispuhf, ispghf, n1c,
     $               n2c, na, nb, norbs, inobas, v(jpant), v(jpatn),
     $               v(jpbnt), v(jpbtn), v(jp2ant), v(jp2atn),
     $               v(jp2bnt), v(jp2btn), v(it1a), v(it1b), v(it2a),
     $               v(it2b), junk, v(ingg), v(ipgg))
              endif


C     Evaluate factor x(g,g').

              call evlygg (iout, iprint, itype, n1c, n2c, na, nb, norbs,
     $             detC, detC2, inobas, dodiag, v(ipiv), v(ieig),
     $             v(iscr1), v(iscr2), lscr1, v(ingg), xgg)


C     Accumulate xgg factors.

              if ( .not. ispuhf .and. .not. ispghf ) then
                xgg = xgg * conjg(civec(n1c)) * civec(n2c)

                v(iint) = v(iint) + dble (xgg)
                v(iint+1) = v(iint+1) + aimag (xgg)
              endif

C     Accumulate density matrix in NO basis.

              if ( itype .eq. 1 ) then
                call amove (2*nosq, v(ipgg), v(ix))
                call zscal (nosq, xgg, v(ix), 1)

                call aadd (2*nosq, v(idmno), v(ix), v(idmno))

              elseif ( itype .eq. 2 ) then
                call amove (2*nosq, v(ipgg), v(ix))
                call zscal (nosq, xgg, v(ix), 1)

                call aadd (2*nosq, v(idmno), v(ix), v(idmno))

                call amove (2*nosq, v(ipgg+2*nosq), v(ix))
                call zscal (nosq, xgg, v(ix), 1)

                call aadd (2*nosq, v(idmno+2*nosq), v(ix),
     $               v(idmno+2*nosq))

              elseif ( itype .eq. 3 .and.
     $                 (.not. ispuhf) .and. (.not. ispghf) ) then
                call amove (8*nosq, v(ipgg), v(ix))
                call zscal (4*nosq, xgg, v(ix), 1)

                call aadd (8*nosq, v(idmno), v(ix), v(idmno))
              endif


C     Transform to orthonormal AO basis.
C       nos . P(g) . nos!

              if ( savedm .and. itype .eq. 1 ) then
                call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnos),
     $               norbs, v(ipgg), norbs, zero, v(ix), norbs)
                call zgemm ('n', 'c', norbs, norbs, norbs, one, v(ix),
     $               norbs, v(jnos), norbs, zero, v(iza), norbs)

              elseif ( savedm .and. itype .eq. 2 ) then
                call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnos),
     $               norbs, v(ipgg), norbs, zero, v(ix), norbs)
                call zgemm ('n', 'c', norbs, norbs, norbs, one, v(ix),
     $               norbs, v(jnos), norbs, zero, v(iza), norbs)

                call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $               v(jnos+2*nosq), norbs, v(ipgg+2*nosq), norbs, zero,
     $               v(ix), norbs)
                call zgemm ('n', 'c', norbs, norbs, norbs, one, v(ix),
     $               norbs, v(jnos+2*nosq), norbs, zero, v(izb), norbs)

              elseif ( ( savedm .and. itype .eq. 3 )
     $                 .or. ispuhf .or. ispghf ) then
                call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $               v(jnos), 2*norbs, v(ipgg), 2*norbs, zero, v(ix),
     $               2*norbs)
                call zgemm ('n', 'c', 2*norbs, 2*norbs, 2*norbs, one,
     $               v(ix), 2*norbs, v(jnos), 2*norbs, zero, v(izf),
     $               2*norbs)

C               .. split into spin blocks (in square format) ..

                if ( useao .or. ispuhf .or. ispghf ) then
                  call dmblck (iout, iprint, v(iscr3), lscr3, norbs,
     $                 v(izf), v(iza), v(izab), v(izba), v(izb), 2, -1)
                endif
              endif

              if ( ispuhf .or. ispghf ) then
                call amove (2*nosq, v(iza), v(j2w10))
                call aadd (2*nosq, v(j2w10), v(izb), v(j2w10))
                call zscal (nosq, 1.d0/sqrt(2.0d0)*xgg, v(j2w10), 1)

                call amove (2*nosq, v(iza), v(j2w20))
                call asub (2*nosq, v(j2w20), v(izb), v(j2w20))
                call zscal (nosq, 1.d0/sqrt(2.0d0)*xgg, v(j2w20), 1)

                call amove (2*nosq, v(izba), v(j2w2p1))
                call zscal (nosq, -xgg, v(j2w2p1), 1)

                call amove (2*nosq, v(izab), v(j2w2m1))
                call zscal (nosq, xgg, v(j2w2m1), 1)
              end if


C     Transform all the spin blocks from orthonormal AO to the regular
C     AO basis if the user requested the regular AO basis.

              if ( useao ) then
                call aobstf (iout, iprint, v(iscr3), lscr3, norbs,
     $               nbasis,0, v(jx), v(iza), v(izab), v(izba), v(izb),
     $               .false., key, 1)

C             .. combine spin blocks ..

                if ( itype .eq. 3 ) then
                  call dmblck (iout, iprint, v(iscr3), lscr3, nbasis,
     $                 v(izf), v(iza), v(izab), v(izba), v(izb), 1, -1)
                endif
              endif

C     Prepare the 1PDM.

              if ( savedm ) then
                if ( itype .eq. 1 ) then
                  call zaxpy (nbsq, xgg, v(iza), 1, v(i1pdm), 1)

                elseif ( itype .eq. 2 ) then
                  call zaxpy (nbsq, xgg, v(iza), 1, v(i1pdm), 1)

                  call zaxpy (nbsq, xgg, v(izb), 1, v(i1pdm+2*nbsq), 1)

                elseif ( itype .eq. 3 .and.
     $                   (.not. ispuhf) .and. (.not. ispghf) ) then

                  call zaxpy (4*nbsq, xgg, v(izf), 1, v(i1pdm), 1)
                endif
              endif

              if ( savedm .and. (ispuhf .or. ispghf) ) then
                call amove (2*nbsq, v(iza), v(j2s10))
                call aadd (2*nbsq, v(j2s10), v(izb), v(j2s10))
                call zscal (nbsq, 1.d0/sqrt(2.0d0)*xgg, v(j2s10), 1)

                call amove (2*nbsq, v(iza), v(j2s20))
                call asub (2*nbsq, v(j2s20), v(izb), v(j2s20))
                call zscal (nbsq, 1.d0/sqrt(2.0d0)*xgg, v(j2s20), 1)

                call amove (2*nbsq, v(izba), v(j2s2p1))
                call zscal (nbsq, -xgg, v(j2s2p1), 1)

                call amove (2*nbsq, v(izab), v(j2s2m1))
                call zscal (nbsq, xgg, v(j2s2m1), 1)
              end if


C     Prepare the 2PDM.

              if ( dotwom .and. itype .eq. 1 ) then

                do 1101 i1 = 1, nbsuse
                  do 1102 i2 = 1, nbsuse
                    do 1103 i3 = 1, nbsuse
                      do 1104 i4 = 1, nbsuse

C                       .. < a! a! a a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + iza
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + iza
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + iza
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + iza

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + i2pdm

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < a! b! b a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + iza
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + iza

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + i2pdm
     $                       + 2*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

 1104                 continue
 1103               continue
 1102             continue
 1101           continue

              elseif ( dotwom .and. itype .eq. 2 ) then

                do 1201 i1 = 1, nbsuse
                  do 1202 i2 = 1, nbsuse
                    do 1203 i3 = 1, nbsuse
                      do 1204 i4 = 1, nbsuse

C                       .. < a! a! a a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + iza
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + iza
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + iza
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + iza

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + i2pdm

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < b! b! b b >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izb
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izb
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izb
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izb

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + i2pdm
     $                       + 2*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < a! b! b a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + iza
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izb

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + i2pdm
     $                       + 4*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

 1204                 continue
 1203               continue
 1202             continue
 1201           continue

              elseif ( dotwom .and. itype .eq. 3 .and.
     $                 (.not. ispuhf) .and. (.not. ispghf) ) then

                do 1301 i1 = 1, nbsuse
                  do 1302 i2 = 1, nbsuse
                    do 1303 i3 = 1, nbsuse
                      do 1304 i4 = 1, nbsuse

C                       .. < a! a! a a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + iza
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + iza
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + iza
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + iza

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + i2pdm

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < b! b! b b >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izb
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izb
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izb
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izb

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + i2pdm
     $                       + 2*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < a! b! b a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + iza
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izb
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izba
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izab

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + i2pdm
     $                       + 4*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < a! a! b b >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izba
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izba
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izba
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izba

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + i2pdm
     $                       + 6*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < b! b! a a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izab
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izab
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izab
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izab

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + i2pdm
     $                       + 8*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < a! a! a b >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izba
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + iza
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + iza
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izba

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + i2pdm
     $                       + 10*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < a! b! a a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + iza
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izab
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + iza
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izab

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + i2pdm
     $                       + 12*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < b! b! b a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izab
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izb
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izb
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izab

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + i2pdm
     $                       + 14*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < b! a! b b >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izb
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izba
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izb
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izba

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + i2pdm
     $                       + 16*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

 1304                 continue
 1303               continue
 1302             continue
 1301           continue
              end if


              if ( dotwom .and. (ispuhf .or. ispghf) ) then

                call aclear (32*nbsq*nbsq, v(j2pdm1))

                do 1401 i1 = 1, nbsuse
                  do 1402 i2 = 1, nbsuse
                    do 1403 i3 = 1, nbsuse
                      do 1404 i4 = 1, nbsuse

C                       .. < a! a! a a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + iza
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + iza
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + iza
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + iza

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < b! b! b b >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izb
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izb
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izb
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izb

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1
     $                       + 2*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < a! b! b a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + iza
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izb
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izba
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izab

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1
     $                       + 4*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < b! a! a b >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izb
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + iza
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izab
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izba

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1
     $                       + 6*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < a! b! a b >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izba
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izab
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + iza
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izb

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1
     $                       + 8*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < b! a! b a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izab
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izba
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izb
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + iza

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1
     $                       + 10*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < a! a! b b >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izba
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izba
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izba
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izba

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1
     $                       + 12*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < b! b! a a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izab
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izab
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izab
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izab

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1
     $                       + 14*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < a! a! a b >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izba
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + iza
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + iza
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izba

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1
     $                       + 16*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < a! a! b a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + iza
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izba
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izba
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + iza

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1
     $                       + 18*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < a! b! a a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + iza
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izab
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + iza
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izab

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1
     $                       + 20*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < b! a! a a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izab
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + iza
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izab
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + iza

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1
     $                       + 22*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < b! b! b a >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izab
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izb
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izb
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izab

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1
     $                       + 24*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < b! b! a b >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izb
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izab
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izab
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izb

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1
     $                       + 26*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < b! a! b b >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izb
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izba
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izb
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izba

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1
     $                       + 28*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

C                       .. < a! b! b b >  ..

                        ind1 = (i1-1)*2*nbsuse + (i4-1)*2 + izba
                        ind2 = (i2-1)*2*nbsuse + (i3-1)*2 + izb
                        ind3 = (i1-1)*2*nbsuse + (i3-1)*2 + izba
                        ind4 = (i2-1)*2*nbsuse + (i4-1)*2 + izb

                        indt = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $                       + (i3-1)*2*nbsuse + (i4-1)*2 + j2pdm1
     $                       + 30*nbsq*nbsq

                        fac = xgg / two *
     $                        ( cmplx (v(ind1), v(ind1+1)) *
     $                          cmplx (v(ind2), v(ind2+1)) -
     $                          cmplx (v(ind3), v(ind3+1)) *
     $                          cmplx (v(ind4), v(ind4+1)) )

                        v(indt) = v(indt) + dble (fac)
                        v(indt+1) = v(indt+1) + aimag (fac)

 1404                 continue
 1403               continue
 1402             continue
 1401           continue

C               .. combine into irreducible tensor components ..

C               T1_0

                call amove (2*nbsq*nbsq, v(j1aaaa), v(j2t10))
                call zdscal (nbsq*nbsq, 0.5d0, v(j2t10), 1)
                call daxpy (2*nbsq*nbsq, 0.5d0, v(j1bbbb), 1,
     $               v(j2t10), 1)
                call daxpy (2*nbsq*nbsq, 0.5d0, v(j1abba), 1,
     $               v(j2t10), 1)
                call daxpy (2*nbsq*nbsq, 0.5d0, v(j1baab), 1,
     $               v(j2t10), 1)

C               T2_0

                call amove (2*nbsq*nbsq, v(j1aaaa), v(j2t20))
                call zdscal (nbsq*nbsq, 0.5d0, v(j2t20), 1)
                call daxpy (2*nbsq*nbsq, -0.5d0, v(j1bbbb), 1,
     $               v(j2t20), 1)
                call daxpy (2*nbsq*nbsq, -0.5d0, v(j1abba), 1,
     $               v(j2t20), 1)
                call daxpy (2*nbsq*nbsq, 0.5d0, v(j1baab), 1,
     $               v(j2t20), 1)

C               T2_+1

                call amove (2*nbsq*nbsq, v(j1aaba), v(j2t2p1))
                call zdscal (nbsq*nbsq, -1.d0/sqrt(2.d0), v(j2t2p1), 1)
                call daxpy (2*nbsq*nbsq, -1.d0/sqrt(2.d0), v(j1babb), 1,
     $               v(j2t2p1), 1)

C               T2_-1

                call amove (2*nbsq*nbsq, v(j1abaa), v(j2t2m1))
                call zdscal (nbsq*nbsq, 1.d0/sqrt(2.d0), v(j2t2m1), 1)
                call daxpy (2*nbsq*nbsq, 1.d0/sqrt(2.d0), v(j1bbab), 1,
     $               v(j2t2m1), 1)

C               T3_0

                call amove (2*nbsq*nbsq, v(j1aaaa), v(j2t30))
                call zdscal (nbsq*nbsq, 0.5d0, v(j2t30), 1)
                call daxpy (2*nbsq*nbsq, -0.5d0, v(j1bbbb), 1,
     $               v(j2t30), 1)
                call daxpy (2*nbsq*nbsq, 0.5d0, v(j1abba), 1,
     $               v(j2t30), 1)
                call daxpy (2*nbsq*nbsq, -0.5d0, v(j1baab), 1,
     $               v(j2t30), 1)

C               T3_+1

                call amove (2*nbsq*nbsq, v(j1aaab), v(j2t3p1))
                call zdscal (nbsq*nbsq, -1.d0/sqrt(2.d0), v(j2t3p1), 1)
                call daxpy (2*nbsq*nbsq, -1.d0/sqrt(2.d0), v(j1abbb), 1,
     $               v(j2t3p1), 1)

C               T3_-1

                call amove (2*nbsq*nbsq, v(j1baaa), v(j2t3m1))
                call zdscal (nbsq*nbsq, 1.d0/sqrt(2.d0), v(j2t3m1), 1)
                call daxpy (2*nbsq*nbsq, 1.d0/sqrt(2.d0), v(j1bbba), 1,
     $               v(j2t3m1), 1)

C               T4_0

                call amove (2*nbsq*nbsq, v(j1aaaa), v(j2t40))
                call zdscal (nbsq*nbsq, -1.d0/(2.d0*sqrt(3.d0)),
     $               v(j2t40), 1)
                call daxpy (2*nbsq*nbsq, -1.d0/(2.d0*sqrt(3.d0)),
     $               v(j1bbbb), 1, v(j2t40), 1)
                call daxpy (2*nbsq*nbsq, 1.d0/(2.d0*sqrt(3.d0)),
     $               v(j1abba), 1, v(j2t40), 1)
                call daxpy (2*nbsq*nbsq, 1.d0/(2.d0*sqrt(3.d0)),
     $               v(j1baab), 1, v(j2t40), 1)
                call daxpy (2*nbsq*nbsq, -1.d0/sqrt(3.d0),
     $               v(j1abab), 1, v(j2t40), 1)
                call daxpy (2*nbsq*nbsq, -1.d0/sqrt(3.d0),
     $               v(j1baba), 1, v(j2t40), 1)

C               T5_0

c$$$                call amove (2*nbsq*nbsq, v(j1abab), v(j2t50))
c$$$                call zdscal (nbsq*nbsq, -1.d0/sqrt(2.d0), v(j2t50), 1)
c$$$                call daxpy (2*nbsq*nbsq, 1.d0/sqrt(2.d0), v(j1baba), 1,
c$$$     $               v(j2t50), 1)

C               T5_p1

c$$$                call amove (2*nbsq*nbsq, v(j1aaab), v(j2t5p1))
c$$$                call zdscal (nbsq*nbsq, -0.5d0, v(j2t5p1), 1)
c$$$                call daxpy (2*nbsq*nbsq, 0.5d0, v(j1abbb), 1,
c$$$     $               v(j2t5p1), 1)
c$$$                call daxpy (2*nbsq*nbsq, 0.5d0, v(j1aaba), 1,
c$$$     $               v(j2t5p1), 1)
c$$$                call daxpy (2*nbsq*nbsq, -0.5d0, v(j1babb), 1,
c$$$     $               v(j2t5p1), 1)

C               T5_m1

c$$$                call amove (2*nbsq*nbsq, v(j1baaa), v(j2t5m1))
c$$$                call zdscal (nbsq*nbsq, -0.5d0, v(j2t5m1), 1)
c$$$                call daxpy (2*nbsq*nbsq, 0.5d0, v(j1bbba), 1,
c$$$     $               v(j2t5m1), 1)
c$$$                call daxpy (2*nbsq*nbsq, 0.5d0, v(j1abaa), 1,
c$$$     $               v(j2t5m1), 1)
c$$$                call daxpy (2*nbsq*nbsq, -0.5d0, v(j1bbab), 1,
c$$$     $               v(j2t5m1), 1)

C               T6_0

                call amove (2*nbsq*nbsq, v(j1aaaa), v(j2t60))
                call zdscal (nbsq*nbsq, 1.d0/sqrt(6.d0),
     $               v(j2t60), 1)
                call daxpy (2*nbsq*nbsq, 1.d0/sqrt(6.d0),
     $               v(j1bbbb), 1, v(j2t60), 1)
                call daxpy (2*nbsq*nbsq, -1.d0/sqrt(6.d0),
     $               v(j1abba), 1, v(j2t60), 1)
                call daxpy (2*nbsq*nbsq, -1.d0/sqrt(6.d0),
     $               v(j1baab), 1, v(j2t60), 1)
                call daxpy (2*nbsq*nbsq, -1.d0/sqrt(6.d0),
     $               v(j1abab), 1, v(j2t60), 1)
                call daxpy (2*nbsq*nbsq, -1.d0/sqrt(6.d0),
     $               v(j1baba), 1, v(j2t60), 1)

C               T6_p1

                call amove (2*nbsq*nbsq, v(j1aaab), v(j2t6p1))
                call zdscal (nbsq*nbsq, -0.5d0, v(j2t6p1), 1)
                call daxpy (2*nbsq*nbsq, 0.5d0, v(j1abbb), 1,
     $               v(j2t6p1), 1)
                call daxpy (2*nbsq*nbsq, -0.5d0, v(j1aaba), 1,
     $               v(j2t6p1), 1)
                call daxpy (2*nbsq*nbsq, 0.5d0, v(j1babb), 1,
     $               v(j2t6p1), 1)

C               T6_m1

                call amove (2*nbsq*nbsq, v(j1baaa), v(j2t6m1))
                call zdscal (nbsq*nbsq, 0.5d0, v(j2t6m1), 1)
                call daxpy (2*nbsq*nbsq, -0.5d0, v(j1bbba), 1,
     $               v(j2t6m1), 1)
                call daxpy (2*nbsq*nbsq, 0.5d0, v(j1abaa), 1,
     $               v(j2t6m1), 1)
                call daxpy (2*nbsq*nbsq, -0.5d0, v(j1bbab), 1,
     $               v(j2t6m1), 1)

C               T6_p2

                call amove (2*nbsq*nbsq, v(j1aabb), v(j2t6p2))

C               T6_m2

                call amove (2*nbsq*nbsq, v(j1bbaa), v(j2t6m2))

              endif


C     Loop over CI matrix elements (spin projection).

              if ( ispuhf .or. ispghf ) then

                if ( ispuhf ) is1 = isz
                if ( ispghf ) is1 = s

                do 201 n1s = 1, ncis

                  if ( ispuhf ) is2 = isz
                  if ( ispghf ) is2 = s

                  do 202 n2s = 1, ncis

                    n1 = (n1c-1)*ncis + n1s
                    n2 = (n2c-1)*ncis + n2s

                    fac = conjg(civec(n1)) * civec(n2) * facw

C     Wigner and Clebsch-Gordan coefficients

                    if ( isint ) then

                      wgtf0  = wignerd (isint, s, is1, is2, angb) *
     $                         exp (+im * anga * dble(is1)) *
     $                         exp (+im * angg * dble(is2))

                      cgf2_0 = clgord (s, is1, isint, 2, 0, .true., 
     $                                 s, is1, isint)
                      cgf1_0 = clgord (s, is1, isint, 1, 0, .true.,
     $                                 s, is1, isint)
                      cgf0_0 = clgord (s, is1, isint, 0, 0, .true.,
     $                                 s, is1, isint)

                      if ( abs(is1-1) .le. s ) then
                      wgtfp1  = wignerd (isint, s, is1-1, is2, angb) *
     $                          exp (+im * anga * dble(is1-1)) *
     $                          exp (+im * angg * dble(is2))

                      cgf2_p1 = clgord (s, is1, isint, 2, -1, .true., 
     $                                  s, is1-1, isint)
                      cgf1_p1 = clgord (s, is1, isint, 1, -1, .true.,
     $                                  s, is1-1, isint)

                      else
                        wgtfp1  = 0.0d0
                        cgf2_p1 = 0.0d0
                        cgf1_p1 = 0.0d0
                      end if

                      if ( abs(is1+1) .le. s ) then
                      wgtfm1  = wignerd (isint, s, is1+1, is2, angb) *
     $                          exp (+im * anga * dble(is1+1)) *
     $                          exp (+im * angg * dble(is2))

                      cgf2_m1 = clgord (s, is1, isint, 2, +1, .true., 
     $                                  s, is1+1, isint)
                      cgf1_m1 = clgord (s, is1, isint, 1, +1, .true.,
     $                                  s, is1+1, isint)

                      else
                        wgtfm1  = 0.0d0
                        cgf2_m1 = 0.0d0
                        cgf1_m1 = 0.0d0
                      end if

                      if ( abs(is1-2) .le. s ) then
                      wgtfp2  = wignerd (isint, s, is1-2, is2, angb) *
     $                          exp (+im * anga * dble(is1-2)) *
     $                          exp (+im * angg * dble(is2))

                      cgf2_p2 = clgord (s, is1, isint, 2, -2, .true., 
     $                                  s, is1-2, isint)

                      else
                        wgtfp2  = 0.0d0
                        cgf2_p2 = 0.0d0
                      end if

                      if ( abs(is1+2) .le. s ) then
                      wgtfm2  = wignerd (isint, s, is1+2, is2, angb) *
     $                          exp (+im * anga * dble(is1+2)) *
     $                          exp (+im * angg * dble(is2))

                      cgf2_m2 = clgord (s, is1, isint, 2, +2, .true., 
     $                                  s, is1+2, isint)

                      else
                        wgtfm2  = 0.0d0
                        cgf2_m2 = 0.0d0
                      end if

                    else

                      wgtf0  = wignerd (isint, s, is1, is2, angb) *
     $                         exp (+im * anga * dble(is1)/2.0d0) *
     $                         exp (+im * angg * dble(is2)/2.0d0)

                      cgf2_0 = clgord (s, is1, isint, 2, 0, .true., 
     $                                 s, is1, isint)
                      cgf1_0 = clgord (s, is1, isint, 1, 0, .true.,
     $                                 s, is1, isint)
                      cgf0_0 = clgord (s, is1, isint, 0, 0, .true.,
     $                                 s, is1, isint)

                      if ( abs(is1-2) .le. s ) then
                      wgtfp1  = wignerd (isint, s, is1-2, is2, angb) *
     $                          exp (+im * anga * dble(is1-2)/2.0d0) *
     $                          exp (+im * angg * dble(is2)/2.0d0)

                      cgf2_p1 = clgord (s, is1, isint, 2, -1, .true., 
     $                                  s, is1-2, isint)
                      cgf1_p1 = clgord (s, is1, isint, 1, -1, .true.,
     $                                  s, is1-2, isint)

                      else
                        wgtfp1  = 0.0d0
                        cgf2_p1 = 0.0d0
                        cgf1_p1 = 0.0d0
                      end if

                      if ( abs(is1+2) .le. s ) then
                      wgtfm1  = wignerd (isint, s, is1+2, is2, angb) *
     $                          exp (+im * anga * dble(is1+2)/2.0d0) *
     $                          exp (+im * angg * dble(is2)/2.0d0)

                      cgf2_m1 = clgord (s, is1, isint, 2, +1, .true., 
     $                                  s, is1+2, isint)
                      cgf1_m1 = clgord (s, is1, isint, 1, +1, .true.,
     $                                  s, is1+2, isint)

                      else
                        wgtfm1  = 0.0d0
                        cgf2_m1 = 0.0d0
                        cgf1_m1 = 0.0d0
                      end if

                      if ( abs(is1-4) .le. s ) then
                      wgtfp2  = wignerd (isint, s, is1-4, is2, angb) *
     $                          exp (+im * anga * dble(is1-4)/2.0d0) *
     $                          exp (+im * angg * dble(is2)/2.0d0)

                      cgf2_p2 = clgord (s, is1, isint, 2, -2, .true., 
     $                                  s, is1-4, isint)

                      else
                        wgtfp2  = 0.0d0
                        cgf2_p2 = 0.0d0
                      end if

                      if ( abs(is1+4) .le. s ) then
                      wgtfm2  = wignerd (isint, s, is1+4, is2, angb) *
     $                          exp (+im * anga * dble(is1+4)/2.0d0) *
     $                          exp (+im * angg * dble(is2)/2.0d0)

                      cgf2_m2 = clgord (s, is1, isint, 2, +2, .true., 
     $                                  s, is1+4, isint)

                      else
                        wgtfm2  = 0.0d0
                        cgf2_m2 = 0.0d0
                      end if

                    end if

                    v(iint) = v(iint) + dble (xgg*fac*wgtf0)
                    v(iint+1) = v(iint+1) + aimag (xgg*fac*wgtf0)

                    call amove (2*nosq, v(j2w10), v(j2w1c))
                    call zscal (nosq, wgtf0*cgf0_0, v(j2w1c), 1)

                    call amove (2*nosq, v(j2w20), v(j2w2c))
                    call zscal (nosq, wgtf0*cgf1_0, v(j2w2c), 1)
                    call zaxpy (nosq, -wgtfp1*cgf1_p1, v(j2w2p1), 1,
     $                   v(j2w2c), 1)
                    call zaxpy (nosq, -wgtfm1*cgf1_m1, v(j2w2m1), 1,
     $                   v(j2w2c), 1)

                    if ( savedm ) then
                      call amove (2*nbsq, v(j2s10), v(j2s1c))
                      call zscal (nbsq, wgtf0*cgf0_0, v(j2s1c), 1)

                      call amove (2*nbsq, v(j2s20), v(j2s2c))
                      call zscal (nbsq, wgtf0*cgf1_0, v(j2s2c), 1)
                      call zaxpy (nbsq, -wgtfp1*cgf1_p1, v(j2s2p1), 1,
     $                     v(j2s2c), 1)
                      call zaxpy (nbsq, -wgtfm1*cgf1_m1, v(j2s2m1), 1,
     $                     v(j2s2c), 1)
                    endif

                    if ( dotwom ) then
                      call amove (2*nbsq*nbsq, v(j2t10), v(j2t1c))
                      call zscal (nbsq*nbsq, wgtf0*cgf0_0, v(j2t1c), 1)

                      call amove (2*nbsq*nbsq, v(j2t20), v(j2t2c))
                      call zscal (nbsq*nbsq, wgtf0*cgf1_0, v(j2t2c), 1)
                      call zaxpy (nbsq*nbsq, -wgtfp1*cgf1_p1, v(j2t2p1),
     $                     1, v(j2t2c), 1)
                      call zaxpy (nbsq*nbsq, -wgtfm1*cgf1_m1, v(j2t2m1),
     $                     1, v(j2t2c), 1)

                      call amove (2*nbsq*nbsq, v(j2t30), v(j2t3c))
                      call zscal (nbsq*nbsq, wgtf0*cgf1_0, v(j2t3c), 1)
                      call zaxpy (nbsq*nbsq, -wgtfp1*cgf1_p1, v(j2t3p1),
     $                     1, v(j2t3c), 1)
                      call zaxpy (nbsq*nbsq, -wgtfm1*cgf1_m1, v(j2t3m1),
     $                     1, v(j2t3c), 1)

                      call amove (2*nbsq*nbsq, v(j2t40), v(j2t4c))
                      call zscal (nbsq*nbsq, wgtf0*cgf0_0, v(j2t4c), 1)
c$$$
c$$$                      call amove (2*nbsq*nbsq, v(j2t50), v(j2t5c))
c$$$                      call zscal (nbsq*nbsq, wgtf0*cgf1_0, v(j2t5c), 1)
c$$$                      call zaxpy (nbsq*nbsq, -wgtfp1*cgf1_p1, v(j2t5p1),
c$$$     $                     1, v(j2t5c), 1)
c$$$                      call zaxpy (nbsq*nbsq, -wgtfm1*cgf1_m1, v(j2t5m1),
c$$$     $                     1, v(j2t5c), 1)

                      call amove (2*nbsq*nbsq, v(j2t60), v(j2t6c))
                      call zscal (nbsq*nbsq, wgtf0*cgf2_0, v(j2t6c), 1)
                      call zaxpy (nbsq*nbsq, -wgtfp1*cgf2_p1, v(j2t6p1),
     $                     1, v(j2t6c), 1)
                      call zaxpy (nbsq*nbsq, -wgtfm1*cgf2_m1, v(j2t6m1),
     $                     1, v(j2t6c), 1)
                      call zaxpy (nbsq*nbsq, +wgtfp2*cgf2_p2, v(j2t6p2),
     $                     1, v(j2t6c), 1)
                      call zaxpy (nbsq*nbsq, +wgtfm2*cgf2_m2, v(j2t6m2),
     $                     1, v(j2t6c), 1)
                    endif

                    call zaxpy (nosq, +fac*cgf0*sqrt(2.d0)/2.d0,
     $                   v(j2w1c), 1, v(idmno), 1)
                    call zaxpy (nosq, +fac*cgf1*sqrt(2.d0)/2.d0,
     $                   v(j2w2c), 1, v(idmno), 1)

                    call zaxpy (nosq, +fac*cgf0*sqrt(2.d0)/2.d0,
     $                   v(j2w1c), 1, v(idmno+2*nosq), 1)
                    call zaxpy (nosq, -fac*cgf1*sqrt(2.d0)/2.d0,
     $                   v(j2w2c), 1, v(idmno+2*nosq), 1)

                    if ( savedm ) then
                      call zaxpy (nbsq, +fac*cgf0*sqrt(2.d0)/2.d0,
     $                     v(j2s1c), 1, v(i1pdm), 1)
                      call zaxpy (nbsq, +fac*cgf1*sqrt(2.d0)/2.d0,
     $                     v(j2s2c), 1, v(i1pdm), 1)

                      call zaxpy (nbsq, +fac*cgf0*sqrt(2.d0)/2.d0,
     $                     v(j2s1c), 1, v(i1pdm+2*nbsq), 1)
                      call zaxpy (nbsq, -fac*cgf1*sqrt(2.d0)/2.d0,
     $                     v(j2s2c), 1, v(i1pdm+2*nbsq), 1)
                    endif

                    if ( dotwom ) then
                    call zaxpy (nbsq*nbsq, +fac*cgf0/2.d0,
     $                   v(j2t1c), 1, v(i2pdm), 1)
                    call zaxpy (nbsq*nbsq, -fac*cgf0/(2.d0*sqrt(3.d0)),
     $                   v(j2t4c), 1, v(i2pdm), 1)
                    call zaxpy (nbsq*nbsq, +fac*cgf1/2.d0,
     $                   v(j2t2c), 1, v(i2pdm), 1)
                    call zaxpy (nbsq*nbsq, +fac*cgf1/2.d0,
     $                   v(j2t3c), 1, v(i2pdm), 1)
                    call zaxpy (nbsq*nbsq, +fac*cgf2/sqrt(6.d0),
     $                   v(j2t6c), 1, v(i2pdm), 1)

                    call zaxpy (nbsq*nbsq, +fac*cgf0/2.d0,
     $                   v(j2t1c), 1, v(i2pdm+2*nbsq*nbsq), 1)
                    call zaxpy (nbsq*nbsq, -fac*cgf0/(2.d0*sqrt(3.d0)),
     $                   v(j2t4c), 1, v(i2pdm+2*nbsq*nbsq), 1)
                    call zaxpy (nbsq*nbsq, -fac*cgf1/2.d0,
     $                   v(j2t2c), 1, v(i2pdm+2*nbsq*nbsq), 1)
                    call zaxpy (nbsq*nbsq, -fac*cgf1/2.d0,
     $                   v(j2t3c), 1, v(i2pdm+2*nbsq*nbsq), 1)
                    call zaxpy (nbsq*nbsq, +fac*cgf2/sqrt(6.d0),
     $                   v(j2t6c), 1, v(i2pdm+2*nbsq*nbsq), 1)

                    call zaxpy (nbsq*nbsq, +fac*cgf0/2.d0,
     $                   v(j2t1c), 1, v(i2pdm+4*nbsq*nbsq), 1)
                    call zaxpy (nbsq*nbsq, +fac*cgf0/(2.d0*sqrt(3.d0)),
     $                   v(j2t4c), 1, v(i2pdm+4*nbsq*nbsq), 1)
                    call zaxpy (nbsq*nbsq, -fac*cgf1/2.d0,
     $                   v(j2t2c), 1, v(i2pdm+4*nbsq*nbsq), 1)
                    call zaxpy (nbsq*nbsq, +fac*cgf1/2.d0,
     $                   v(j2t3c), 1, v(i2pdm+4*nbsq*nbsq), 1)
                    call zaxpy (nbsq*nbsq, -fac*cgf2/sqrt(6.d0),
     $                   v(j2t6c), 1, v(i2pdm+4*nbsq*nbsq), 1)
                    endif


                    if ( isint ) then
                      is2 = is2 - 1
                    else
                      is2 = is2 - 2
                    endif
 202              continue

                  if ( isint ) then
                    is1 = is1 - 1
                  else
                    is1 = is1 - 2
                  endif
 201            continue
              end if

 102        continue
 101      continue

 104    continue

C$omp   end parallel do

 103  continue



C     Merge results from all processors.

      jdmno = jend  + 2
      jend1 = jdmno + szp1

      if ( savedm ) then
        j1pdm = jend1
        jend1 = j1pdm + sz1pdm
      endif

      if ( dotwom ) then
        j2pdm = jend1
        jend1 = j2pdm + sz2pdm
      endif

      intxgg = cmplx (v(jend), v(jend+1))

      do 21 ip = 2, npuse
        iint  = jend  + (ip - 1) * lenp
        idmno = iint  + 2
        iend  = idmno + szp1

        intxgg = intxgg + cmplx (v(iint), v(iint+1))

        call aadd (szp1, v(jdmno), v(idmno), v(jdmno))

        if ( savedm ) then
          i1pdm = iend
          iend  = i1pdm + sz1pdm

          call aadd (sz1pdm, v(j1pdm), v(i1pdm), v(j1pdm))
        endif

        if ( dotwom ) then
          i2pdm = iend
          iend  = i2pdm + sz2pdm

          call aadd (sz2pdm, v(j2pdm), v(i2pdm), v(j2pdm))
        endif
 21   continue

      jend = jend1


C     Scale density matrices by overlap factor.

      call zscal (szp1/2, one/intxgg, v(jdmno), 1)

      if ( savedm ) then
        call zscal (sz1pdm/2, one/intxgg, v(j1pdm), 1)
      endif

      if ( dotwom ) then
        call zscal (sz2pdm/2, one/intxgg, v(j2pdm), 1)
      endif


C     Print matrices if requested.

       if ( savedm .and. iprint .gt. 2 ) then
        write (iout, *) ' '
        if ( itype1 .eq. 1 ) then
          call ltoutx (iout, 1, 2, 1, 0, 'PHF density matrix', nbsuse,
     $         nbsuse, v(j1pdm), v(j1pdm), 0)

        elseif ( itype1 .eq. 2 ) then
          call ltoutx (iout, 1, 3, 1, 0, 'PHF density matrix', nbsuse,
     $         nbsuse, v(j1pdm), v(j1pdm+2*nbsq), 0)

        elseif ( itype1 .eq. 3 ) then
          call ltoutx (iout, 1, 2, 1, 0, 'PHF density matrix',
     $         2*nbsuse, 2*nbsuse, v(j1pdm), v(j1pdm), 0)
        endif
      endif


C     Save density matrices in RWF files.

      if ( savedm ) then
        call conddf (irw1dm, sz1pdm)
        call fileio (1, -irw1dm, sz1pdm, v(j1pdm), 0)
      endif

      if ( dotwom ) then
        call conddf (irw2dm, sz2pdm)
        call fileio (1, -irw2dm, sz2pdm, v(j2pdm), 0)
      endif


C     Print density matrices to external file.

      if ( prtdm ) then
        call prtdmt (dotwom, nbsuse, itype1, v(j1pdm), v(j2pdm))
      endif



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Diagonalize one-particle density matrix  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Diagonalize one-particle density matrix to get natural orbital
C     occupations.


C     Hermitize density matrix and transform to lower triangular.
C     Apply a scale factor of -1 (so that larger eigenvalues appear first).

      if ( itype1 .eq. 1 ) then
        call apadag (1, norbs, 0.50d0, v(jdmno), v(jnost))

        call trspnc (.false., norbs, norbs, v(jnost))
        call zlinear (v(jnost), v(jp), norbs, norbs)

        call zdscal (ntto, -1.0d0, v(jp), 1)

      elseif ( itype1 .eq. 2 ) then
        call apadag (1, norbs, 0.50d0, v(jdmno), v(jnost))

        call trspnc (.false., norbs, norbs, v(jnost))
        call zlinear (v(jnost), v(jp), norbs, norbs)

        call zdscal (ntto, -1.0d0, v(jp), 1)

        call apadag (1, norbs, 0.50d0, v(jdmno+2*nosq), v(jnost))

        call trspnc (.false., norbs, norbs, v(jnost))
        call zlinear (v(jnost), v(jp+2*ntto), norbs, norbs)

        call zdscal (ntto, -1.0d0, v(jp+2*ntto), 1)

      elseif ( itype1 .eq. 3 ) then
        call apadag (1, 2*norbs, 0.50d0, v(jdmno), v(jnost))

        call trspnc (.false., 2*norbs, 2*norbs, v(jnost))
        call zlinear (v(jnost), v(jp), 2*norbs, 2*norbs)

        call zdscal (ntt2o, -1.0d0, v(jp), 1)
      endif


C     Diagonalize the density matrix.

      if ( itype1 .eq. 1 ) then
        call vconjg (ntto, v(jp))

        call zhpev ('n', 'u', norbs, v(jp), v(jocc), junk, 1, v(jwrk1),
     $       v(jwrk2), info)

        if ( info .ne. 0 ) then
          call gauerr ('zhpev failed in phfdm1')
        endif

      elseif ( itype1 .eq. 2 ) then
        call vconjg (ntto, v(jp))
        call vconjg (ntto, v(jp+2*ntto))

        call zhpev ('n', 'u', norbs, v(jp), v(jocc), junk, 1, v(jwrk1),
     $       v(jwrk2), info)

        if ( info .ne. 0 ) then
          call gauerr ('zhpev failed in phfdm1')
        endif

        call zhpev ('n', 'u', norbs, v(jp+2*ntto), v(jocc+norbs), junk,
     $       1, v(jwrk1), v(jwrk2), info)

        if ( info .ne. 0 ) then
          call gauerr ('zhpev failed in phfdm1')
        endif

      elseif ( itype1 .eq. 3 ) then
        call vconjg (ntt2o, v(jp))

        call zhpev ('n', 'u', 2*norbs, v(jp), v(jocc), junk, 1,
     $       v(jwrk1), v(jwrk2), info)

        if ( info .ne. 0 ) then
          call gauerr ('zhpev failed in phfdm1')
        endif
      endif


C     Print natural orbital occupations.

      if ( itype1 .eq. 1 ) then
        call ascale (norbs, -1.0d0, v(jocc), v(jocc))

        call prtocc (iout, 'Natural occupations', 'orbital',
     $       'occupation', norbs, v(jocc))

      elseif ( itype1 .eq. 2 ) then
        call ascale (2*norbs, -1.0d0, v(jocc), v(jocc))

        call prtocc (iout, 'Natural occupations (alpha)', 'orbital',
     $       'occupation', norbs, v(jocc))
        call prtocc (iout, 'Natural occupations (beta)', 'orbital',
     $       'occupation', norbs, v(jocc+norbs))

      elseif ( itype1 .eq. 3 ) then
        call ascale (2*norbs, -1.0d0, v(jocc), v(jocc))

        call prtocc (iout, 'Natural occupations', 'orbital',
     $       'occupation', 2*norbs, v(jocc))
      endif



C     %%%%%%%%%%%%%%%%%%%%%
C     %  Compute entropy  %
C     %%%%%%%%%%%%%%%%%%%%%


C     The entropy is computed according to

C       S = - sum_i n_i log n_i

C     Note that we actually report S / S_0, where S_0 is given
C     by the state with equal occupations among all the orbitals.


      s0 = 0.0d0
      s1 = 0.0d0

      seq = dble (nel) / dble (2*norbs)

      s0 = - dble (2*norbs) * seq * log (seq)
     $     - dble (2*norbs) * (1.0d0 - seq) * log (1.0d0 - seq)

      if ( itype1 .eq. 1 ) then
        do 401 i1 = 1, norbs
          ind1 = (i1-1) + jocc

          if ( v(ind1) .gt. 1.0d-12 ) then
            s1 = s1 - v(ind1) * log (v(ind1))
          elseif ( 1.0d0 - v(ind1) .gt. 1.0d-12 ) then
            s1 = s1 - (1.0d0 - v(ind1)) * log (1.0d0 - v(ind1))
          endif
 401    continue

        s1 = 2.0d0 * s1

      elseif ( itype1 .eq. 2 .or. itype1 .eq. 3 ) then
        do 402 i1 = 1, 2*norbs
          ind1 = (i1-1) + jocc
 
          if ( v(ind1) .gt. 1.0d-12 ) then
            s1 = s1 - v(ind1) * log (v(ind1))
          elseif ( 1.0d0 - v(ind1) .gt. 1.0d-12 ) then
            s1 = s1 - (1.0d0 - v(ind1)) * log (1.0d0 - v(ind1))
          endif
 402    continue
      endif


      write (iout, *) ' '
      write (iout, 403) s1/s0
 403  format (5X, ' S / S0 = ', 6X, F20.12)
      write (iout, *) ' '



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  DEBUGGING - Return if false  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      if ( .not. debug ) then
        return
      endif



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Check density matrices  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%



      write (iout, *)
      write (iout, *) ' debugging density matrices ... '
      write (iout, *)


C     Single-particle energy
C     ======================


C     Memory allocation

C       Allocate space for:
C         - core Hamiltonian matrix (AO basis) (jh1)
C         - core Hamiltonian matrix (orthonormal basis) (jh2)
C         - scratch array to transform core Hamiltonian (jscr4)

      jh1   = jend
      jh2   = jh1   + nri*nbasis*nbasis
      jscr4 = jh2   + 2*nosq
      jend  = jscr4 + nri*nbasis*norbs


C       Allocate space for:
C         - core Hamiltonian matrix in NO basis (jhno)

      if ( itype1 .eq. 1 ) then
        jhno = jend
        jend = jhno + 2*nosq

      elseif ( itype1 .eq. 2 ) then
        jhno = jend
        jend = jhno + 4*nosq

      elseif ( itype1 .eq. 3 ) then
        jhno = jend
        jend = jhno + 8*nosq
      endif

C       Allocate some extra scratch arrays for itype = 3.

      if ( itype1 .eq. 3 ) then
        lscrd = 16*nosq

        jzab  = jend
        jscrd = jzab  + 2*nosq
        jend  = jscrd + lscrd
      endif

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'phfdm1')


C     Load core Hamiltonian matrix.
C     Transform into a square matrix.

      call fileio (2, -irwh, nri*ntt, v(jh1), 0)

      if ( nri .eq. 1 ) then
        call square (v(jh1), v(jh1), nbasis, nbasis, 0)
      elseif ( nri .eq. 2 ) then
        call zsquare (v(jh1), v(jh1), nbasis, nbasis, 0)
      endif


C     Transform core Hamiltonian matrix into orthonormal AO basis by
C     computing the product X! . H . X.

      if ( nri .eq. 1 ) then
        call dgemm ('t', 'n', norbs, nbasis, nbasis, 1.0d0, v(jx),
     $       nbasis, v(jh1), nbasis, 0.0d0, v(jscr4), norbs)
        call dgemm ('n', 'n', norbs, norbs, nbasis, 1.0d0, v(jscr4),
     $       norbs, v(jx), nbasis, 0.0d0, v(jh1), norbs)

C       .. transform core Hamiltonian matrix into complex ..

        call zrecmp (3, v(jh2), v(jh1), nosq)

      elseif ( nri .eq. 2 ) then
        call zgemm ('c', 'n', norbs, nbasis, nbasis, one, v(jx),
     $       nbasis, v(jh1), nbasis, zero, v(jscr4), norbs)
        call zgemm ('n', 'n', norbs, norbs, nbasis, one, v(jscr4),
     $       norbs, v(jx), nbasis, zero, v(jh1), norbs)

        call amove (2*nosq, v(jh1), v(jh2))
      endif


C     Transform core Hamiltonian into NO basis:
C       nos! . H . nos

      if ( itype1 .eq. 1 ) then

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jh2), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos), norbs, zero, v(jhno), norbs)

      elseif ( itype1 .eq. 2 ) then

        if ( ispuhf .or. ispghf ) then
          call amove (2*nosq, v(jh2), v(jhno))
          call amove (2*nosq, v(jh2), v(jhno+2*nosq))

        else
        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jh2), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos), norbs, zero, v(jhno), norbs)

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos+2*nosq),
     $       norbs, v(jh2), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos+2*nosq), norbs, zero, v(jhno+2*nosq), norbs)
        endif

      elseif ( itype1 .eq. 3 ) then

C       .. form core Hamiltonian in spin orbital basis ..

        call aclear (2*nosq, v(jzab))

        call dmblck (iout, iprint, v(jscrd), lscrd, norbs, v(jhno),
     $       v(jh2), v(jzab), v(jzab), v(jh2), 1, -1)

C       .. perform the transformation ..

        call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jnos),
     $       2*norbs, v(jhno), 2*norbs, zero, v(jnost), 2*norbs)
        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jnost),
     $       2*norbs, v(jnos), 2*norbs, zero, v(jhno), 2*norbs)
      endif


C     Take complex conjugate of core Hamiltonian (to take traces
C     appropriately).

      if ( itype1 .eq. 1 ) then
        call vconjg (nosq, v(jhno))

      elseif ( itype1 .eq. 2 ) then
        call vconjg (nosq, v(jhno))
        call vconjg (nosq, v(jhno+2*nosq))

      elseif ( itype1 .eq. 3 ) then
        call vconjg (4*nosq, v(jhno))
      endif


C     Evaluate the single-particle energy.

      trh = zero

      if ( itype1 .eq. 1 ) then
        trh = trcabc (1, norbs, v(jdmno), v(jhno))
        trh = two * trh

      elseif ( itype1 .eq. 2 ) then
        trh = trcabc (1, norbs, v(jdmno), v(jhno))
        trh = trh + 
     $        trcabc (1, norbs, v(jdmno+2*nosq), v(jhno+2*nosq))

      elseif ( itype1 .eq. 3 ) then
        trh = trcabc (1, 2*norbs, v(jdmno), v(jhno))
      endif


C     Print single-particle energy.

      write (iout, 3001) dble (trh)
 3001 format (5X, ' Hsp = ', 9X, F20.12)



C     Trace of 1PDM
C     =============


C     Compute trace of 1PDM.

      if ( savedm .and. .not. useao ) then

        trp = zero

        if ( itype1 .eq. 1 ) then
          do 2101 i1 = 1, nbsuse
            ind = (i1-1)*2*nbsuse + (i1-1)*2 + j1pdm
            trp = trp + cmplx (v(ind), v(ind+1))
 2101     continue

          trp = two * trp

        elseif ( itype1 .eq. 2 ) then
          do 2102 i1 = 1, nbsuse
            ind = (i1-1)*2*nbsuse + (i1-1)*2 + j1pdm
            trp = trp + cmplx (v(ind), v(ind+1))
 2102     continue

          do 2103 i1 = 1, nbsuse
            ind = (i1-1)*2*nbsuse + (i1-1)*2 + j1pdm + 2*nbsq
            trp = trp + cmplx (v(ind), v(ind+1))
 2103     continue

        elseif ( itype1 .eq. 3 ) then
          do 2104 i1 = 1, 2*nbsuse
            ind = (i1-1)*4*nbsuse + (i1-1)*2 + j1pdm
            trp = trp + cmplx (v(ind), v(ind+1))
 2104     continue
        endif


C     Print trace of 1PDM.

        write (iout, 3002) dble (trp)
 3002   format (5X, ' trace (1PDM) = ', F20.12)

      endif



C     Trace of 2PDM
C     =============


C     Compute trace of 2PDM.

      if ( savedm .and. dotwom .and. .not. useao ) then

        trp = zero

        if ( itype1 .eq. 1 ) then
          do 2201 i1 = 1, nbsuse
            do 2202 i2 = 1, nbsuse
              ind = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $            + (i2-1)*2*nbsuse + (i1-1)*2 + j2pdm

              trp = trp + cmplx (v(ind), v(ind+1))

              ind = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $            + (i2-1)*2*nbsuse + (i1-1)*2 + j2pdm 
     $            + 2*nbsq*nbsq

              trp = trp + cmplx (v(ind), v(ind+1))
 2202       continue
 2201     continue

          trp = two * trp

        elseif ( itype1 .eq. 2 .or. itype1 .eq. 3 ) then
          do 2203 i1 = 1, nbsuse
            do 2204 i2 = 1, nbsuse
              ind = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $            + (i2-1)*2*nbsuse + (i1-1)*2 + j2pdm

              trp = trp + cmplx (v(ind), v(ind+1))

              ind = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $            + (i2-1)*2*nbsuse + (i1-1)*2 + j2pdm 
     $            + 2*nbsq*nbsq

              trp = trp + cmplx (v(ind), v(ind+1))

              ind = (i1-1)*2*nbcb + (i2-1)*2*nbsq
     $            + (i2-1)*2*nbsuse + (i1-1)*2 + j2pdm 
     $            + 4*nbsq*nbsq

              trp = trp + two * cmplx (v(ind), v(ind+1))
 2204       continue
 2203     continue
        endif


C     Print trace of 2PDM.

        write (iout, 3003) dble (trp)
 3003   format (5X, ' trace (2PDM) = ', F20.12)

      endif



C     Does the 2PDM trace to the 1PDM?
C     ================================


      if ( savedm .and. dotwom .and.
     $     iprint .gt. 2 .and. .not. useao ) then

        if ( itype1 .eq. 1 ) then
          jdmtst = jend
          jend   = jdmtst + 2*nbsq

          call aclear (2*nbsq, v(jdmtst))

        elseif ( itype1 .eq. 2 ) then
          jdmtst = jend
          jend   = jdmtst + 4*nbsq

          call aclear (4*nbsq, v(jdmtst))

        elseif ( itype1 .eq. 3 ) then
          jdmtst = jend
          jend   = jdmtst + 8*nbsq

          call aclear (8*nbsq, v(jdmtst))
        endif

        call tstcor (jend-1, lenv, 'phfdm1')


        if ( itype1 .eq. 1 ) then
          do 2301 i1 = 1, nbsuse
            do 2302 i2 = 1, nbsuse
              do 2303 i3 = 1, nbsuse
                ind1 = (i1-1)*2*nbsuse + (i2-1)*2 + jdmtst

                indt = (i1-1)*2*nbcb + (i3-1)*2*nbsq
     $               + (i3-1)*2*nbsuse + (i2-1)*2 + j2pdm

                v(ind1) = v(ind1) + v(indt)
                v(ind1+1) = v(ind1+1) + v(indt+1)

                indt = (i1-1)*2*nbcb + (i3-1)*2*nbsq
     $               + (i3-1)*2*nbsuse + (i2-1)*2 + j2pdm
     $               + 2*nbsq*nbsq

                v(ind1) = v(ind1) + v(indt)
                v(ind1+1) = v(ind1+1) + v(indt+1)
 2303         continue
 2302       continue
 2301     continue

        elseif ( itype1 .eq. 2 ) then
          do 2304 i1 = 1, nbsuse
            do 2305 i2 = 1, nbsuse
              do 2306 i3 = 1, nbsuse
                ind1 = (i1-1)*2*nbsuse + (i2-1)*2 + jdmtst
                ind2 = (i1-1)*2*nbsuse + (i2-1)*2 + jdmtst + 2*nbsq

                indt = (i1-1)*2*nbcb + (i3-1)*2*nbsq
     $               + (i3-1)*2*nbsuse + (i2-1)*2 + j2pdm

                v(ind1) = v(ind1) + v(indt)
                v(ind1+1) = v(ind1+1) + v(indt+1)

                indt = (i1-1)*2*nbcb + (i3-1)*2*nbsq
     $               + (i3-1)*2*nbsuse + (i2-1)*2 + j2pdm
     $               + 2*nbsq*nbsq

                v(ind2) = v(ind2) + v(indt)
                v(ind2+1) = v(ind2+1) + v(indt+1)

                indt = (i1-1)*2*nbcb + (i3-1)*2*nbsq
     $               + (i3-1)*2*nbsuse + (i2-1)*2 + j2pdm
     $               + 4*nbsq*nbsq

                v(ind1) = v(ind1) + v(indt)
                v(ind1+1) = v(ind1+1) + v(indt+1)

                indt = (i3-1)*2*nbcb + (i1-1)*2*nbsq
     $               + (i2-1)*2*nbsuse + (i3-1)*2 + j2pdm
     $               + 4*nbsq*nbsq

                v(ind2) = v(ind2) + v(indt)
                v(ind2+1) = v(ind2+1) + v(indt+1)
 2306         continue
 2305       continue
 2304     continue

        elseif ( itype1 .eq. 3 ) then
          do 2307 i1 = 1, nbsuse
            do 2308 i2 = 1, nbsuse
              do 2309 i3 = 1, nbsuse
                ind1 = (i1-1)*4*nbsuse + (i2-1)*2 + jdmtst
                ind2 = (i1+nbsuse-1)*4*nbsuse + (i2+nbsuse-1)*2 + jdmtst
                ind3 = (i1+nbsuse-1)*4*nbsuse + (i2-1)*2 + jdmtst
                ind4 = (i1-1)*4*nbsuse + (i2+nbsuse-1)*2 + jdmtst

                indt = (i1-1)*2*nbcb + (i3-1)*2*nbsq
     $               + (i3-1)*2*nbsuse + (i2-1)*2 + j2pdm

                v(ind1) = v(ind1) + v(indt)
                v(ind1+1) = v(ind1+1) + v(indt+1)

                indt = (i1-1)*2*nbcb + (i3-1)*2*nbsq
     $               + (i3-1)*2*nbsuse + (i2-1)*2 + j2pdm
     $               + 2*nbsq*nbsq

                v(ind2) = v(ind2) + v(indt)
                v(ind2+1) = v(ind2+1) + v(indt+1)

                indt = (i1-1)*2*nbcb + (i3-1)*2*nbsq
     $               + (i3-1)*2*nbsuse + (i2-1)*2 + j2pdm
     $               + 4*nbsq*nbsq

                v(ind1) = v(ind1) + v(indt)
                v(ind1+1) = v(ind1+1) + v(indt+1)

                indt = (i3-1)*2*nbcb + (i1-1)*2*nbsq
     $               + (i2-1)*2*nbsuse + (i3-1)*2 + j2pdm
     $               + 4*nbsq*nbsq

                v(ind2) = v(ind2) + v(indt)
                v(ind2+1) = v(ind2+1) + v(indt+1)

                indt = (i3-1)*2*nbcb + (i1-1)*2*nbsq
     $               + (i2-1)*2*nbsuse + (i3-1)*2 + j2pdm
     $               + 12*nbsq*nbsq

                v(ind3) = v(ind3) + v(indt)
                v(ind3+1) = v(ind3+1) + v(indt+1)

                indt = (i1-1)*2*nbcb + (i3-1)*2*nbsq
     $               + (i3-1)*2*nbsuse + (i2-1)*2 + j2pdm
     $               + 14*nbsq*nbsq

                v(ind3) = v(ind3) + v(indt)
                v(ind3+1) = v(ind3+1) + v(indt+1)

                indt = (i3-1)*2*nbcb + (i1-1)*2*nbsq
     $               + (i2-1)*2*nbsuse + (i3-1)*2 + j2pdm
     $               + 16*nbsq*nbsq

                v(ind4) = v(ind4) + v(indt)
                v(ind4+1) = v(ind4+1) + v(indt+1)

                indt = (i1-1)*2*nbcb + (i3-1)*2*nbsq
     $               + (i3-1)*2*nbsuse + (i2-1)*2 + j2pdm
     $               + 10*nbsq*nbsq

                v(ind4) = v(ind4) + v(indt)
                v(ind4+1) = v(ind4+1) + v(indt+1)
 2309         continue
 2308       continue
 2307     continue
        endif

        sclfac = 2.0d0 / dble (nel-1)

        if ( itype1 .eq. 1 ) then
          call zdscal (nosq, sclfac, v(jdmtst), 1)
        elseif ( itype1 .eq. 2 ) then
          call zdscal (2*nosq, sclfac, v(jdmtst), 1)
        elseif ( itype1 .eq. 2 ) then
          call zdscal (4*nosq, sclfac, v(jdmtst), 1)
        endif


C     Print 1PDM obtained from 2PDM.

        write (iout, *) ' '

        if ( itype1 .eq. 1 ) then
          call ltoutx (iout, 1, 2, 1, 0, 'density matrix test', nbsuse,
     $         nbsuse, v(jdmtst), v(jdmtst), 0)

        elseif ( itype1 .eq. 2 ) then
          call ltoutx (iout, 1, 3, 1, 0, 'density matrix test', nbsuse,
     $         nbsuse, v(jdmtst), v(jdmtst+2*nbsq), 0)

        elseif ( itype1 .eq. 3 ) then
          call ltoutx (iout, 1, 2, 1, 0, 'density matrix test',
     $         2*nbsuse, 2*nbsuse, v(jdmtst), v(jdmtst), 0)
        endif

      endif


      write (iout, *)


      return
      end


