

      subroutine hffock (iout, iprint, v, lenv, iphf, norbs, nbasis,
     $     nbas6d, irwp, irwf, f, irwh, irwx, energy, allowp, ipflag,
     $     jsym2e, nsymop, nop1, neqbas, neqshl, neqsh2, rotop, neqatm,
     $     accdes)

      implicit none


C +------------------------------------------------------------+
C |  hffock   --   CAJH, 08.2011                               |
C |                                                            |
C |                                                            |
C |  Build the Fock matrix associated with the HF state and    |
C |  compute the energy (sum of single- and two-particle       |
C |  energies) associated with such state.                     |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       norbs  - number of orbitals
C       nbasis - number of basis functions (pure)
C       nbas6d - number of basis functions (cartesian)
C       irwp   - RWF file with HF density matrix
C       irwf   - RWF file with saved HF Fock matrix  [ updated here ]
C       f      - HF Fock matrix  [ out ]
C       irwh   - RWF file with core Hamiltonian matrix
C       irwx   - RWF file with transformation matrix (= S^(-1/2))
C       energy - sum of single and two-particle energies  [ out ]

      real*8      v(*), f(*), energy, accdes
      integer     iout, iprint, lenv, iphf
      integer     norbs, nbasis, nbas6d
      integer     irwp, irwf, irwh, irwx

C       .. prism variables .. [ input variables ]

      real*8      rotop(*)
      integer     ipflag, nsymop, nop1, jsym2e
      integer     neqbas(nbasis,*), neqshl(*), neqsh2(*), neqatm(*)
      logical     allowp(50)

C     other variables

      real*8      junk
      complex*16  enerc
      integer     ntt, ntt6d, ntto, ntt2o, nbsq, nosq, itype
      integer     iprtf, key, iopcl, nri
      integer     nmaxs, nmaxt, nmaxa
      integer     cntas, cntat, cntaa, cntbs, cntbt, cntba
      integer     jstrt, jend, mdv
      integer     jh1, jh2, jx, jscr1, jf, jf2, jz1, jz2, jz3
      integer     jzaa, jzab, jzba, jzbb, jzf
      integer     jscrrs, jscrra, jscris, jscria
      integer     jscr, lscr, lscr1, lscr2, lscr3, lscr4
      integer     jxas, jxat, jxaa, jyas, jyat, jyaa
      integer     jxbs, jxbt, jxba, jybs, jybt, jyba
      integer     ind1, ind2, ind3, ind4
      logical     icmplx, ispuhf, ispghf

C     functions

      complex*16  trcabc

C     constants

      complex*16  zero, one


      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)

      junk = 0.0d0


C     Printing level for PRISM...

C      iprtf = iprint - 2
C      if ( iprtf .lt. 0 ) iprtf = 0
      iprtf = 1  ! wsr

C     Construct useful variables.

      ntt   = nbasis * (nbasis + 1) / 2
      ntt6d = nbas6d * (nbas6d + 1) / 2

      ntto = norbs * (norbs + 1) / 2
      ntt2o = norbs * (2*norbs + 1)

      nbsq = nbasis * nbasis
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
        call gauerr ('Unrecognized iphf in hffock.')
      endif


C     Read iopcl from ILSW file.

C       iopcl = 0,  real RHF
C             = 1,  real UHF
C             = 2,  complex RHF
C             = 3,  complex UHF
C             = 6,  complex GHF (there is no real GHF)

      call ilsw (2, 1, iopcl)

      if ( iopcl .gt. 3 .and. iopcl .ne. 6 ) then
        call gauerr ('Incorrect iopcl in hffock.')
      endif

      nri = 2

      if ( iopcl .le. 1 ) nri = 1


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
C         - Fock matrix (jf)
C         - scratch array to transform to LT (jf2)

      if ( itype .eq. 1 ) then
        jf   = jend
        jf2  = jf   + 2*nosq
        jend = jf2  + 2*nosq

      elseif ( itype .eq. 2 .or. ispuhf ) then
        jf   = jend
        jf2  = jf   + 4*nosq
        jend = jf2  + 2*nosq

      elseif ( itype .eq. 3 ) then
        jf   = jend
        jf2  = jf   + 8*nosq
        jend = jf2  + 8*nosq
      endif

C       Allocate space to store:
C         - scratch arrays for lt to square transformation (jz1, jz2, jz3)

      if ( itype .eq. 1 ) then
        jz1  = jend
        jz2  = jz1  + 2*ntto
        jz3  = jz2  + nosq
        jend = jz3  + nosq

      elseif ( itype .eq. 2 .or. ispuhf ) then
        jz1  = jend
        jz2  = jz1  + 2*ntto
        jz3  = jz2  + nosq
        jend = jz3  + nosq

      elseif ( itype .eq. 3 ) then
        jz1  = jend
        jz2  = jz1  + 2*ntt2o
        jz3  = jz2  + 4*nosq
        jend = jz3  + 4*nosq
      endif

C       Allocate space for:
C         - scratch square spin blocks (jzaa, jzab, ...)
C         - scratch array for the full matrix (for itype = 3) (jzf)

      if ( itype .eq. 1 ) then
        jzaa = jend
        jend = jzaa + 2*nbsq

      elseif ( itype .eq. 2 .or. ispuhf ) then
        jzaa = jend
        jzbb = jzaa + 2*nbsq
        jend = jzbb + 2*nbsq

      elseif ( itype .eq. 3 ) then
        jzaa = jend
        jzab = jzaa + 2*nbsq
        jzba = jzab + 2*nbsq
        jzbb = jzba + 2*nbsq
        jzf  = jzbb + 2*nbsq
        jend = jzf  + 8*nbsq
      endif

C       Allocate space for:
C         - four scratch lower triangular arrays
C             (real-symmetric, real-antisymmetric,
C              imag-symmetric, imag-antisymmetric)

      jscrrs = jend
      jscrra = jscrrs + ntt
      jscris = jscrra + ntt
      jscria = jscris + ntt
      jend   = jscria + ntt

C       Allocate space for:
C         - scratch array to call aobstf, dmblck, and sptblk

C         - lscr1, length of scratch space required for aobstf
C         - lscr2, length of scratch space required for sptblk, mode = 1
C         - lscr3, length of scratch space required for sptblk, mode = 2
C         - lscr4, length of scratch space required for dmblck

      lscr1 = 4*nbasis*norbs
      lscr2 = 2*nbsq + 4*ntt
      lscr3 = 3*nbsq

      lscr2 = max (lscr2, lscr3)
      lscr1 = max (lscr1, lscr2)

      if ( itype .eq. 3 .and. .not. ispuhf ) then
        lscr4 = 16*nosq
        lscr = max (lscr1, lscr4)
      else
        lscr = lscr1
      endif

      jscr = jend
      jend = jscr + lscr

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'hffock')



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Add core-Hamiltonian to Fock matrix  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

       write (iout, *) ' '
       call ltoutx (iout, 1, 2, 1, 0, 'Core Hamiltonian', norbs,
     $      norbs, v(jh2), v(jh2), 0)


C     Add contributions of H to Fock matrix.

      if ( itype .eq. 1 ) then
        call aclear (2*nosq, v(jf))

        call aadd (2*nosq, v(jf), v(jh2), v(jf))

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call aclear (4*nosq, v(jf))

        call aadd (2*nosq, v(jf), v(jh2), v(jf))
        call aadd (2*nosq, v(jf+2*nosq), v(jh2), v(jf+2*nosq))

      elseif ( itype .eq. 3 ) then

        call aclear (2*nosq, v(jzab))
        call aclear (2*nosq, v(jzba))

C         .. combine spin blocks ..

        call dmblck (iout, iprint, v(jscr), lscr, norbs, v(jzf),
     $       v(jh2), v(jzab), v(jzba), v(jh2), 1, -1)

        call aclear (8*nosq, v(jf))

        call aadd (8*nosq, v(jf), v(jzf), v(jf))
      endif



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Do some preparation to call PRISM later  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Specify iopcl.

C       iopcl = 0,  when only alpha matrices are sent to PRISM
C             = 1,  when both alpha and beta matrices are sent to PRISM

      if ( itype .eq. 1 ) then
        iopcl = 0
      elseif ( itype .eq. 2 .or. itype .eq. 3 ) then
        iopcl = 1
      endif


C     key controls which elements need to be transformed in any given
C     call to the basis set transformation routine (aobstf). The   
C     documentation of this routine says:

C       key = 1,  all spin blocks are active
C           = 2,  [aa] and [bb] blocks are active
C           = 3,  [aa] block is active

C     Therefore, we adjust key accordingly.

      if ( itype .eq. 1 ) then
        key = 3
      elseif ( itype .eq. 2 .or. ispuhf ) then
        key = 2
      elseif ( itype .eq. 3 ) then
        key = 1
      endif


C       Specify the maximum number of singlet, triplet, and antisymmetric
C       blocks for each itype case.

C       This can be determined by noting that the spin blocks of P
C       can be decomposed into:

C              spin block     # singlet   # triplet   # antisymm
C              -------------------------------------------------
C                [ aa ]           1                        1
C                [ ab ]                        2           2
C                [ ba ]                        2           2
C                [ bb ]           1                        1

C         ** The [aa] and [bb] blocks are 'singlets' because they need to
C            be contracted with Coulomb and exchange integrals, while the
C            [ab] and [ba] blocks only need to be contracted with exchange
C            integrals. All antisymmetric matrices are considered as
C            antisymmetric.

C         ** We exploit the fact that P_ab = P_ba!

C       Furthermore, for itype = 2, 3, we will send both alpha and beta
C       matrices, while for itype = 1 we will only send alpha matrices:

C         [ aa ] block is sent as alpha
C         [ ab ] block is sent as alpha, beta
C         [ bb ] block is sent as beta
C         [ ba ] block is NOT contracted

      if ( itype .eq. 1 .or. itype .eq. 2 .or. ispuhf ) then
        nmaxs = 1
        nmaxt = 0
        nmaxa = 1

      elseif ( itype .eq. 3 ) then
        nmaxs = 1
        nmaxt = 1
        nmaxa = 2
      endif


C       Allocate space for:
C         - arrays storing all the matrices to send to PRISM per batch:
C             x - density matrices (sent to PRISM)
C             y - Fock matrices (returned from PRISM)

C           Naming convention:
C             jxas - jx, alpha - singlet
C             jxat - jx, alpha - triplet
C             jxaa - jx, alpha - antisymmetric
C             ...

C       ** Note that PRISM requires that all matrices be allocated ntt6d,
C          as opposed to ntt. However, all matrices (singlet, triplet,
C          antisymmetric) are sent as a single block. That's why
C          allocation may look tricky.

      jstrt = jend

      if ( itype .lt. 4 ) then
        jxas = jstrt
        jxat = jxas + ntt * nmaxs
        jxaa = jxat + ntt * nmaxt
        jyas = jxas + ntt6d * (nmaxs + nmaxt + nmaxa)
        jyat = jyas + ntt * nmaxs
        jyaa = jyat + ntt * nmaxt
        jend = jyas + ntt6d * (nmaxs + nmaxt + nmaxa)
      endif

      if ( itype .eq. 2 .or. itype .eq. 3 ) then
        jxbs = jend
        jxbt = jxbs + ntt * nmaxs
        jxba = jxbt + ntt * nmaxt
        jybs = jxbs + ntt6d * (nmaxs + nmaxt + nmaxa)
        jybt = jybs + ntt * nmaxs
        jyba = jybt + ntt * nmaxt
        jend = jybs + ntt6d * (nmaxs + nmaxt + nmaxa)
      endif

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'hffock')



C     %%%%%%%%%%%%
C     %  Form G  %
C     %%%%%%%%%%%%

C     Contract P with two-electron integrals...


C     Recover all blocks of P and organize them in preparation for
C     erictr, which actually calls FoFDir to perform the contraction.

C       cnt?s - keeps track of number of singlet matrices
C       cnt?t - keeps track of number of triplet matrices
C       cnt?a - keeps track of number of antisymmetric matrices

      cntas = 0
      cntat = 0
      cntaa = 0

      cntbs = 0
      cntbt = 0
      cntba = 0


      if ( itype .eq. 1 ) then
        call fileio (2, -irwp, 2*ntto, v(jzaa), 0)

        call zrecmp (1, v(jzaa), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jzaa), v(jz2), nosq)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call fileio (2, -irwp, 2*ntto, v(jzaa), 0)
        call fileio (2, -irwp, 2*ntto, v(jzbb), 2*ntto)

        call zrecmp (1, v(jzaa), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jzaa), v(jz2), nosq)

        call zrecmp (1, v(jzbb), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jzbb), v(jz2), nosq)

      elseif ( itype .eq. 3 ) then
        call fileio (2, -irwp, 2*ntt2o, v(jzf), 0)

        call zrecmp (1, v(jzf), v(jz1), ntt2o)
        call square (v(jz1), v(jz2), 2*norbs, 2*norbs, 0)
        call square (v(jz1+ntt2o), v(jz3), 2*norbs, 2*norbs, -1)
        call zrecmp (2, v(jzf), v(jz2), 4*nosq)

C         .. split into spin blocks (in square format) ..

        call dmblck (iout, iprint, v(jscr), lscr, norbs, v(jzf),
     $       v(jzaa), v(jzab), v(jzba), v(jzbb), 2, -1)
      endif


C     Transform all the spin blocks from orthonormal AO to the regular
C     AO basis.

      call aobstf (iout, iprint, v(jscr), lscr, norbs, nbasis, irwx,
     $     junk, v(jzaa), v(jzab), v(jzba), v(jzbb), .true., key, 1)


C     Split each block into real-imaginary, symmetric-antisymmetric.

C       .. [alpha, alpha] block

      if ( itype .lt. 4 ) then
        ind1 = jxas + cntas * ntt
        ind2 = jxaa + cntaa * ntt

        call sptblk (iout, iprint, v(jscr), lscr, nbasis, v(jzaa),
     $       v(jscrrs), v(jscrra), v(jscris), v(jscria), 1)

        call amove (ntt, v(jscrrs), v(ind1))
        call amove (ntt, v(jscria), v(ind2))
          
        cntas = cntas + 1
        cntaa = cntaa + 1
      endif

C       .. [beta, beta] block

      if ( itype .gt. 1 ) then
        ind1 = jxbs + cntbs * ntt
        ind2 = jxba + cntba * ntt

        call sptblk (iout, iprint, v(jscr), lscr, nbasis, v(jzbb),
     $       v(jscrrs), v(jscrra), v(jscris), v(jscria), 1)

        call amove (ntt, v(jscrrs), v(ind1))
        call amove (ntt, v(jscria), v(ind2))
          
        cntbs = cntbs + 1
        cntba = cntba + 1
      endif

C       .. [alpha, beta] block

      if ( itype .eq. 3 .and. .not. ispuhf ) then
        ind1 = jxat + cntat * ntt
        ind2 = jxaa + cntaa * ntt
        ind3 = jxbt + cntbt * ntt
        ind4 = jxba + cntba * ntt

        call sptblk (iout, iprint, v(jscr), lscr, nbasis, v(jzab),
     $       v(jscrrs), v(jscrra), v(jscris), v(jscria), 1)

        call amove (ntt, v(jscrrs), v(ind1))
        call amove (ntt, v(jscrra), v(ind4))
        call amove (ntt, v(jscris), v(ind3))
        call amove (ntt, v(jscria), v(ind2))
          
        cntat = cntat + 1
        cntaa = cntaa + 1
        cntbt = cntbt + 1
        cntba = cntba + 1
      endif


C     Contract P against two-electron integrals in erictr.

      if ( iprtf .gt. 0 ) write (iout, *)
      call erictr (iout, iprtf, v(jend), mdv, nbasis, iopcl, cntas,
     $     cntat, cntaa, v(jxas), v(jxbs), v(jyas), v(jybs), ipflag,
     $     allowp, jsym2e, nsymop, nop1, neqbas, neqshl, neqsh2, rotop,
     $     neqatm, accdes)


C     Collect all matrices into appropriate places in G.

      cntas = 0
      cntat = 0
      cntaa = 0

      cntbs = 0
      cntbt = 0
      cntba = 0

C     Merge the real-imaginary, symmetric-antisymmetric blocks into a
C     single spin-block.

C       .. [alpha, alpha] block

      if ( itype .lt. 4 ) then
        ind1 = jyas + cntas * ntt
        ind2 = jyaa + cntaa * ntt

C         .. scale the antisymmetric blocks by a factor of -1 ..

        call ascale (ntt, -1.0d0, v(ind2), v(ind2))
    
        call amove (ntt, v(ind1), v(jscrrs))
        call amove (ntt, v(ind2), v(jscria))

        call aclear (ntt, v(jscrra))
        call aclear (ntt, v(jscris))

        call sptblk (iout, iprint, v(jscr), lscr, nbasis, v(jzaa),
     $       v(jscrrs), v(jscrra), v(jscris), v(jscria), 2)

        cntas = cntas + 1
        cntaa = cntaa + 1
      endif

C       .. [beta, beta] block

      if ( itype .gt. 1 ) then
        ind1 = jybs + cntbs * ntt
        ind2 = jyba + cntba * ntt

C         .. scale the antisymmetric blocks by a factor of -1 ..

        call ascale (ntt, -1.0d0, v(ind2), v(ind2))

        call amove (ntt, v(ind1), v(jscrrs))
        call amove (ntt, v(ind2), v(jscria))

        call aclear (ntt, v(jscrra))
        call aclear (ntt, v(jscris))

        call sptblk (iout, iprint, v(jscr), lscr, nbasis, v(jzbb),
     $       v(jscrrs), v(jscrra), v(jscris), v(jscria), 2)
          
        cntbs = cntbs + 1
        cntba = cntba + 1
      endif

C       .. [alpha, beta] block

      if ( itype .eq. 3 .and. .not. ispuhf ) then
        ind1 = jyat + cntat * ntt
        ind2 = jyaa + cntaa * ntt
        ind3 = jybt + cntbt * ntt
        ind4 = jyba + cntba * ntt

C         .. scale the antisymmetric blocks by a factor of -1 ..

        call ascale (ntt, -1.0d0, v(ind2), v(ind2))
        call ascale (ntt, -1.0d0, v(ind4), v(ind4))

        call amove (ntt, v(ind1), v(jscrrs))
        call amove (ntt, v(ind4), v(jscrra))
        call amove (ntt, v(ind3), v(jscris))
        call amove (ntt, v(ind2), v(jscria))
          
        call sptblk (iout, iprint, v(jscr), lscr, nbasis, v(jzab),
     $       v(jscrrs), v(jscrra), v(jscris), v(jscria), 2)

        cntat = cntat + 1
        cntaa = cntaa + 1
        cntbt = cntbt + 1
        cntba = cntba + 1
      endif

C       .. [beta, alpha] block

C          ** Generated as the dagger of [alpha, beta] block.

      if ( itype .eq. 3 .and. .not. ispuhf ) then
        call amove (2*nbsq, v(jzab), v(jzba))
        call trspnc (.true., nbasis, nbasis, v(jzba))
      endif


C     Transform all the spin blocks to the orthonormal AO basis.

      call aobstf (iout, iprint, v(jscr), lscr, norbs, nbasis, irwx,
     $     junk, v(jzaa), v(jzab), v(jzba), v(jzbb), .true., key, 2)


C     Add contributions of G to Fock matrix.

      if ( itype .eq. 1 ) then
        call aadd (2*nosq, v(jf), v(jzaa), v(jf))

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call aadd (2*nosq, v(jf), v(jzaa), v(jf))
        call aadd (2*nosq, v(jf+2*nosq), v(jzbb), v(jf+2*nosq))
C        write (iout, *) ' '
C        call ltoutx (iout, 0, 2, 1, 0, 'Veff matrix', norbs, norbs,
C     $         f, f, 0)

      elseif ( itype .eq. 3 ) then

C         .. combine spin blocks ..

        call dmblck (iout, iprint, v(jscr), lscr, norbs, v(jzf),
     $       v(jzaa), v(jzab), v(jzba), v(jzbb), 1, -1)

        call aadd (8*nosq, v(jf), v(jzf), v(jf))
      endif



C     %%%%%%%%%%%%%%%%%%%%%%
C     %  Save Fock matrix  %
C     %%%%%%%%%%%%%%%%%%%%%%


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


C     Save Fock matrix into appropriate RWF file.

      if ( itype .eq. 1 ) then
        call fileio (1, -irwf, 2*ntto, f, 0)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call fileio (1, -irwf, 4*ntto, f, 0)

      elseif ( itype .eq. 3 ) then
        call fileio (1, -irwf, 2*ntt2o, f, 0)
      endif


C     Print HF Fock matrix.

      if ( iprint .gt. 3 ) then
        if ( itype .eq. 1 ) then
          write (iout, *) ' '
          call ltoutx (iout, 0, 2, 1, 0, 'HF Fock matrix', norbs, norbs,
     $         f, f, 0)

        elseif ( itype .eq. 2 .or. ispuhf ) then
          write (iout, *) ' '
          call ltoutx (iout, 0, 3, 1, 0, 'HF Fock matrix', norbs, norbs,
     $         f, f(1+2*ntto), 0)

        elseif ( itype .eq. 3 ) then
          write (iout, *) ' '
          call ltoutx (iout, 0, 2, 1, 0, 'HF Fock matrix', 2*norbs,
     $         2*norbs, f, f, 0)
        endif
      endif



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Compute energy of HF state  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Load density matrix of HF state.

      if ( itype .eq. 1 ) then
        call fileio (2, -irwp, 2*ntto, v(jzaa), 0)

        call zrecmp (1, v(jzaa), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jzaa), v(jz2), nosq)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call fileio (2, -irwp, 2*ntto, v(jzaa), 0)
        call fileio (2, -irwp, 2*ntto, v(jzbb), 2*ntto)

        call zrecmp (1, v(jzaa), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jzaa), v(jz2), nosq)

        call zrecmp (1, v(jzbb), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jzbb), v(jz2), nosq)

      elseif ( itype .eq. 3 ) then
        call fileio (2, -irwp, 2*ntt2o, v(jzf), 0)

        call zrecmp (1, v(jzf), v(jz1), ntt2o)
        call square (v(jz1), v(jz2), 2*norbs, 2*norbs, 0)
        call square (v(jz1+ntt2o), v(jz3), 2*norbs, 2*norbs, -1)
        call zrecmp (2, v(jzf), v(jz2), 4*nosq)

C         .. split into spin blocks (in square format) ..

        call dmblck (iout, iprint, v(jscr), lscr, norbs, v(jzf),
     $       v(jzaa), v(jzab), v(jzba), v(jzbb), 2, -1)
      endif


C     Compute energy.

      if ( itype .eq. 1 ) then
        call vconjg (nosq, v(jzaa))

        enerc = trcabc (1, norbs, v(jh2), v(jzaa))
     $        + trcabc (1, norbs, v(jf), v(jzaa))
        energy = dble (enerc)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call vconjg (nosq, v(jzaa))
        call vconjg (nosq, v(jzbb))

        enerc = trcabc (1, norbs, v(jh2), v(jzaa))
     $        + trcabc (1, norbs, v(jh2), v(jzbb))
     $        + trcabc (1, norbs, v(jf), v(jzaa))
     $        + trcabc (1, norbs, v(jf+2*nosq), v(jzbb))
        energy = dble (enerc) / 2.0d0

      elseif ( itype .eq. 3 ) then
        call vconjg (nosq, v(jzaa))
        call vconjg (nosq, v(jzbb))
        call vconjg (4*nosq, v(jzf))

        enerc = trcabc (1, norbs, v(jh2), v(jzaa))
     $        + trcabc (1, norbs, v(jh2), v(jzbb))
     $        + trcabc (1, 2*norbs, v(jf), v(jzf))
        energy = dble (enerc) / 2.0d0
      endif


      return
      end


