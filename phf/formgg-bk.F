      

      subroutine formgg (iout, iprint, v, lenv, iphf, nci, norbs,
     $     nbasis, nbas6d, ngrdt, irwpg, irwgg, irwx, irwnos, allowp,
     $     ipflag, jsym2e, nsymop, nop1, neqbas, neqshl, neqsh2, rotop,
     $     neqatm, accdes)

      implicit none


C +------------------------------------------------------------+
C |  formgg   --   CAJH, 07.2011                               |
C |                                                            |
C |                                   (based on PHFB form2gg)  |
C |                                                            |
C |                                                            |
C |  Form the matrices G(g) by contracting P(g) with           |
C |  electron-repulsion integrals (Coulomb and exchange).      |
C |                                                            |
C |------------------------------------------------------------|
C |                                                            |
C |  On input, P(g) are stored in some orthonormal basis       |
C |  (might be the NO basis). Before sending them to PRISM,    |
C |  we need to go back to the regular AO basis. After all     |
C |  G(g) matrices are formed, we transform them back to the   |
C |  original basis.                                           |
C |                                                            |
C |    ( The transformation matrix from the orthonormal AO     |
C |    basis to the current basis should be stored in          |
C |    irwnos. We call that basis the NO basis even though     |
C |    it might be a different basis, as explained in          |
C |    evalsmt. )                                              |
C |                                                            |
C |  The routine is currently organized so that we call PRISM  |
C |  (FoFDir) only a few times. We loop over the number of CI  |
C |  matrix elements externally (this is hopefully small).     |
C |  Then, for each CI matrix element, we send the ngrdt       |
C |  matrices P(g) in batches. The number of density matrices  |
C |  sent per batch is decided here based on the amount of     |
C |  memory available (see below).                             |
C |                                                            |
C |  Therefore, we call PRISM at least nCI*nCI times...        |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  Comment: This routine is currently ugly in the sense      |
C |           that it uses a lot of logic. This has to do      |
C |           with the fact that PRISM accepts                 |
C |           lower-triangular density matrices in the         |
C |           spatial orbital basis. Therefore, we have to     |
C |           perform a lot of logical steps to transform      |
C |           the density matrices into such a form.           |
C |                                                            |
C |           Moreover, the contraction in PRISM is done       |
C |           depending on the type of matrix (singlet,        |
C |           triplet, antisymmetric). This introduces         |
C |           another degree of complexity in the logic.       |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  CAJH, 12.2012                                             |
C |                                                            |
C |  Modified to avoid loop over nCI.                          |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       nci    - number of basis states in CI problem
C       norbs  - number of orbitals
C       nbasis - number of basis functions (pure)
C       nbas6d - number of basis functions (cartesian)
C       ngrdt  - total number of grid points
C       ngrid  - number of grid points
C       irwpg  - RWF file with matrices P(g)
C       irwgg  - RWF file with matrices G(g) [ updated here ]
C       irwx   - RWF file with transformation matrix (= S^(-1/2))
C       irwnos - RWF file with transformation matrix from orthonormal AO
C                to NO basis

      real*8      v(*), accdes
      integer     iout, iprint, lenv, iphf, nci
      integer     norbs, nbasis, nbas6d, ngrdt
      integer     irwpg, irwgg, irwx, irwnos

C       .. prism variables .. [ input variables ]

      real*8      rotop(*)
      integer     ipflag, nsymop, nop1, jsym2e
      integer     neqbas(nbasis,*), neqshl(*), neqsh2(*), neqatm(*)
      logical     allowp(50)

C     other variables

      real*8      junk
      integer     ntt, ntt6d, nbsq, nosq, itype
      integer     igrdt, iprtf, key, iopcl
      integer     ngrdtt, n, igrd
      integer     nmaxs, nmaxt, nmaxa
      integer     cntas, cntat, cntaa, cntbs, cntbt, cntba
      integer     jstrt, jend, mdv, ngot, jnos, jnost
      integer     iblock, nblock, ip, szg, nbatch
      integer     jscr, lscr, lscr1, lscr2, lscr3, lscr4
      integer     jzaa, jzab, jzba, jzbb, jzf
      integer     jscrrs, jscrra, jscris, jscria
      integer     jxas, jxat, jxaa, jyas, jyat, jyaa
      integer     jxbs, jxbt, jxba, jybs, jybt, jyba
      integer     indx, indxa, indxb, ind1, ind2, ind3, ind4
      logical     icmplx, ispuhf, ispghf
      character   indc1*5, indc2*3, str*12

C     constants

      complex*16  zero, one


      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)

      junk = 0.0d0


C     Printing level for PRISM...

      iprtf = iprint - 2
      if ( iprtf .lt. 0 ) iprtf = 0


C     Construct useful variables.

      ntt   = nbasis * (nbasis + 1) / 2
      ntt6d = nbas6d * (nbas6d + 1) / 2

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
        call gauerr ('Unrecognized iphf in formgg.')
      endif


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
      elseif ( itype .eq. 2 ) then
        key = 2
      elseif ( itype .eq. 3 ) then
        key = 1
      endif



C     Memory allocation.

C       Allocate space for:
C         - orthonormal AO => NO basis transformation matrix (jnos)
C         - scratch array to perform basis transformation (jnost)

      if ( itype .eq. 1 ) then
        jnos  = 1
        jnost = jnos  + 2*nosq
        jend  = jnost + 2*nosq

      elseif ( itype .eq. 2 ) then
        jnos  = 1
        jnost = jnos  + 4*nosq
        jend  = jnost + 2*nosq

      elseif ( itype .eq. 3 ) then
        jnos  = 1
        jnost = jnos  + 8*nosq
        jend  = jnost + 8*nosq
      endif


C       Allocate space for:
C         - scratch square spin blocks (jzaa, jzab, ...)
C         - scratch array for the full matrix (for itype = 3) (jzf)

      if ( itype .eq. 1 ) then
        jzaa = jend
        jend = jzaa + 2*nbsq

      elseif ( itype .eq. 2 ) then
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

      if ( itype .eq. 3 ) then
        lscr4 = 16*nosq
        lscr = max (lscr1, lscr4)
      else
        lscr = lscr1
      endif

      jscr = jend
      jend = jscr + lscr


C       Specify the maximum number of singlet, triplet, and antisymmetric
C       blocks (per grid point) for each itype case.

C       This can be determined by noting that the spin blocks of P(g)
C       can be decomposed into:

C              spin block     # singlet   # triplet   # antisymm
C              -------------------------------------------------
C                [ aa ]           2                        2
C                [ ab ]                        2           2
C                [ ba ]                        2           2
C                [ bb ]           2                        2

C         ** The [aa] and [bb] blocks are 'singlets' because they need to
C            be contracted with Coulomb and exchange integrals, while the
C            [ab] and [ba] blocks only need to be contracted with exchange
C            integrals. All antisymmetric matrices are considered as
C            antisymmetric.

C       Furthermore, for itype = 2, 3, we will send both alpha and beta
C       matrices, while for itype = 1 we will only send alpha matrices:

C         [ aa ] block is sent as alpha
C         [ ab ] block is sent as alpha
C         [ bb ] block is sent as beta
C         [ ba ] block is sent as beta


      if ( itype .eq. 1 .or. itype .eq. 2 ) then
        nmaxs = 2
        nmaxt = 0
        nmaxa = 2

      elseif ( itype .eq. 3 ) then
        nmaxs = 2
        nmaxt = 2
        nmaxa = 4
      endif


C       Determine the batch size (number of grid points sent in a single
C       batch). Our current strategy is rather naive: we strive to use
C       only half the available memory in storing the density and Fock
C       matrices to be sent to PRISM. This strategy could be improved
C       once we know how to properly estimate the memory requirements in
C       a call to FoFDir.

      ngot = lenv - jend + 1

      if ( itype .eq. 1 ) then
        szg = 2*ntt6d * (nmaxs + nmaxt + nmaxa)

      elseif ( itype .eq. 2 .or. itype .eq. 3 ) then
        szg = 4*ntt6d * (nmaxs + nmaxt + nmaxa)
      endif

      ngrdtt = ngrdt*nci*nci

      nbatch = ngot / (2*szg)
      nbatch = min (nbatch, ngrdtt)


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
        jxat = jxas + ntt * nbatch * nmaxs
        jxaa = jxat + ntt * nbatch * nmaxt
        jyas = jxas + ntt6d * nbatch * (nmaxs + nmaxt + nmaxa)
        jyat = jyas + ntt * nbatch * nmaxs
        jyaa = jyat + ntt * nbatch * nmaxt
        jend = jyas + ntt6d * nbatch * (nmaxs + nmaxt + nmaxa)
      endif

      if ( itype .eq. 2 .or. itype .eq. 3 ) then
        jxbs = jend
        jxbt = jxbs + ntt * nbatch * nmaxs
        jxba = jxbt + ntt * nbatch * nmaxt
        jybs = jxbs + ntt6d * nbatch * (nmaxs + nmaxt + nmaxa)
        jybt = jybs + ntt * nbatch * nmaxs
        jyba = jybt + ntt * nbatch * nmaxt
        jend = jybs + ntt6d * nbatch * (nmaxs + nmaxt + nmaxa)
      endif

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'formgg')


C     Load transformation matrix (orthonormal AO => NO basis).

      if ( itype .eq. 1 ) then
        call fileio (2, -irwnos, 2*nosq, v(jnos), 0)

      elseif ( itype .eq. 2 ) then
        call fileio (2, -irwnos, 4*nosq, v(jnos), 0)

      elseif ( itype .eq. 3 ) then
        call fileio (2, -irwnos, 8*nosq, v(jnos), 0)
      endif



C     %%%%%%%%%%%%%%%
C     %  Form G(g)  %
C     %%%%%%%%%%%%%%%

C     Contract P(g) with two-electron integrals...

C     Loop over batches.

      do 101 iblock = 1, ngrdtt, nbatch
        nblock = min (nbatch, ngrdtt - iblock + 1)

C     Memory allocation... (same code as before)

        if ( itype .lt. 4 ) then
          jxas = jstrt
          jxat = jxas + ntt * nblock * nmaxs
          jxaa = jxat + ntt * nblock * nmaxt
          jyas = jxas + ntt6d * nblock * (nmaxs + nmaxt + nmaxa)
          jyat = jyas + ntt * nblock * nmaxs
          jyaa = jyat + ntt * nblock * nmaxt
          jend = jyas + ntt6d * nblock * (nmaxs + nmaxt + nmaxa)
        endif

        if ( itype .eq. 2 .or. itype .eq. 3 ) then
          jxbs = jend
          jxbt = jxbs + ntt * nblock * nmaxs
          jxba = jxbt + ntt * nblock * nmaxt
          jybs = jxbs + ntt6d * nblock * (nmaxs + nmaxt + nmaxa)
          jybt = jybs + ntt * nblock * nmaxs
          jyba = jybt + ntt * nblock * nmaxt
          jend = jybs + ntt6d * nblock * (nmaxs + nmaxt + nmaxa)
        endif


C     Recover all P(g) in a given batch and organize them in preparation
C     for erictr, which actually calls FoFDir to perform the
C     contraction.

C       cnt?s - keeps track of number of singlet matrices
C       cnt?t - keeps track of number of triplet matrices
C       cnt?a - keeps track of number of antisymmetric matrices

        cntas = 0
        cntat = 0
        cntaa = 0

        cntbs = 0
        cntbt = 0
        cntba = 0

        do 102 ip = 1, nblock

          igrdt = iblock + ip - 1

C     Recover P(g) at current grid point. Split it into spin-blocks.

C     Transform from NO basis to orthonormal AO basis:
C       nos . P(g) . nos!

          if ( itype .eq. 1 ) then
            indx = (igrdt - 1) * 2*nosq

            call fileio (2, -irwpg, 2*nosq, v(jzaa), indx)

            call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnos),
     $           norbs, v(jzaa), norbs, zero, v(jnost), norbs)
            call zgemm ('n', 'c', norbs, norbs, norbs, one, v(jnost),
     $           norbs, v(jnos), norbs, zero, v(jzaa), norbs)

          elseif ( itype .eq. 2 ) then
            indxa = (igrdt - 1) * 4*nosq
            indxb = (igrdt - 1) * 4*nosq + 2*nosq

            call fileio (2, -irwpg, 2*nosq, v(jzaa), indxa)
            call fileio (2, -irwpg, 2*nosq, v(jzbb), indxb)

            call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnos),
     $           norbs, v(jzaa), norbs, zero, v(jnost), norbs)
            call zgemm ('n', 'c', norbs, norbs, norbs, one, v(jnost),
     $           norbs, v(jnos), norbs, zero, v(jzaa), norbs)

            call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $           v(jnos+2*nosq), norbs, v(jzbb), norbs, zero,
     $           v(jnost), norbs)
            call zgemm ('n', 'c', norbs, norbs, norbs, one, v(jnost),
     $           norbs, v(jnos+2*nosq), norbs, zero, v(jzbb), norbs)

          elseif ( itype .eq. 3 ) then
            indx = (igrdt - 1) * 8*nosq

            call fileio (2, -irwpg, 8*nosq, v(jzf), indx)

            call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $           v(jnos), 2*norbs, v(jzf), 2*norbs, zero, v(jnost),
     $           2*norbs)
            call zgemm ('n', 'c', 2*norbs, 2*norbs, 2*norbs, one,
     $           v(jnost), 2*norbs, v(jnos), 2*norbs, zero, v(jzf),
     $           2*norbs)

C             .. split into spin blocks (in square format) ..

            call dmblck (iout, iprint, v(jscr), lscr, norbs, v(jzf),
     $           v(jzaa), v(jzab), v(jzba), v(jzbb), 2, -1)
          endif


C     Transform all the spin blocks from orthonormal AO to the regular
C     AO basis.

          call aobstf (iout, iprint, v(jscr), lscr, norbs, nbasis,
     $         irwx, junk, v(jzaa), v(jzab), v(jzba), v(jzbb), .true.,
     $         key, 1)

C     Split each block into real-imaginary, symmetric-antisymmetric.

C           .. [alpha, alpha] block

          if ( itype .lt. 4 ) then
            ind1 = jxas + cntas * ntt
            ind2 = jxaa + cntaa * ntt
            ind3 = ind1 + ntt
            ind4 = ind2 + ntt

            call sptblk (iout, iprint, v(jscr), lscr, nbasis, v(jzaa),
     $           v(jscrrs), v(jscrra), v(jscris), v(jscria), 1)

            call amove (ntt, v(jscrrs), v(ind1))
            call amove (ntt, v(jscrra), v(ind2))
            call amove (ntt, v(jscris), v(ind3))
            call amove (ntt, v(jscria), v(ind4))
        
            cntas = cntas + 2
            cntaa = cntaa + 2
          endif

C           .. [beta, beta] block

          if ( itype .gt. 1 ) then
            ind1 = jxbs + cntbs * ntt
            ind2 = jxba + cntba * ntt
            ind3 = ind1 + ntt
            ind4 = ind2 + ntt

            call sptblk (iout, iprint, v(jscr), lscr, nbasis, v(jzbb),
     $           v(jscrrs), v(jscrra), v(jscris), v(jscria), 1)

            call amove (ntt, v(jscrrs), v(ind1))
            call amove (ntt, v(jscrra), v(ind2))
            call amove (ntt, v(jscris), v(ind3))
            call amove (ntt, v(jscria), v(ind4))
        
            cntbs = cntbs + 2
            cntba = cntba + 2
          endif

C           .. [alpha, beta] block

          if ( itype .eq. 3 ) then
            ind1 = jxat + cntat * ntt
            ind2 = jxaa + cntaa * ntt
            ind3 = ind1 + ntt
            ind4 = ind2 + ntt

            call sptblk (iout, iprint, v(jscr), lscr, nbasis, v(jzab),
     $           v(jscrrs), v(jscrra), v(jscris), v(jscria), 1)

            call amove (ntt, v(jscrrs), v(ind1))
            call amove (ntt, v(jscrra), v(ind2))
            call amove (ntt, v(jscris), v(ind3))
            call amove (ntt, v(jscria), v(ind4))
        
            cntat = cntat + 2
            cntaa = cntaa + 2
          endif

C           .. [beta, alpha] block

          if ( itype .eq. 3 ) then
            ind1 = jxbt + cntbt * ntt
            ind2 = jxba + cntba * ntt
            ind3 = ind1 + ntt
            ind4 = ind2 + ntt

            call sptblk (iout, iprint, v(jscr), lscr, nbasis, v(jzba),
     $           v(jscrrs), v(jscrra), v(jscris), v(jscria), 1)

            call amove (ntt, v(jscrrs), v(ind1))
            call amove (ntt, v(jscrra), v(ind2))
            call amove (ntt, v(jscris), v(ind3))
            call amove (ntt, v(jscria), v(ind4))
        
            cntbt = cntbt + 2
            cntba = cntba + 2
          endif

 102    continue

C     Contract P(g) against two-electron integrals in erictr.

        if ( iprtf .gt. 0 ) write (iout, *)
        call erictr (iout, iprtf, v(jend), mdv, nbasis, iopcl, cntas,
     $       cntat, cntaa, v(jxas), v(jxbs), v(jyas), v(jybs), ipflag,
     $       allowp, jsym2e, nsymop, nop1, neqbas, neqshl, neqsh2,
     $       rotop, neqatm, accdes)


C     Collect all matrices into appropriate places in G(g).

        cntas = 0
        cntat = 0
        cntaa = 0

        cntbs = 0
        cntbt = 0
        cntba = 0

        do 103 ip = 1, nblock

          igrdt = iblock + ip - 1

          n    = mod(igrdt-1,nci*nci) + 1
          igrd = (igrdt-1)/(nci*nci) + 1


C     Merge the real-imaginary, symmetric-antisymmetric blocks into a
C     single spin-block.

C           .. [alpha, alpha] block

          if ( itype .lt. 4 ) then
            ind1 = jyas + cntas * ntt
            ind2 = jyaa + cntaa * ntt
            ind3 = ind1 + ntt
            ind4 = ind2 + ntt

C             .. scale the antisymmetric blocks by a factor of -1 ..

            call ascale (ntt, -1.0d0, v(ind2), v(ind2))
            call ascale (ntt, -1.0d0, v(ind4), v(ind4))

            call amove (ntt, v(ind1), v(jscrrs))
            call amove (ntt, v(ind2), v(jscrra))
            call amove (ntt, v(ind3), v(jscris))
            call amove (ntt, v(ind4), v(jscria))

            call sptblk (iout, iprint, v(jscr), lscr, nbasis, v(jzaa),
     $           v(jscrrs), v(jscrra), v(jscris), v(jscria), 2)

            cntas = cntas + 2
            cntaa = cntaa + 2
          endif

C           .. [beta, beta] block

          if ( itype .gt. 1 ) then
            ind1 = jybs + cntbs * ntt
            ind2 = jyba + cntba * ntt
            ind3 = ind1 + ntt
            ind4 = ind2 + ntt

C             .. scale the antisymmetric blocks by a factor of -1 ..

            call ascale (ntt, -1.0d0, v(ind2), v(ind2))
            call ascale (ntt, -1.0d0, v(ind4), v(ind4))

            call amove (ntt, v(ind1), v(jscrrs))
            call amove (ntt, v(ind2), v(jscrra))
            call amove (ntt, v(ind3), v(jscris))
            call amove (ntt, v(ind4), v(jscria))

            call sptblk (iout, iprint, v(jscr), lscr, nbasis, v(jzbb),
     $           v(jscrrs), v(jscrra), v(jscris), v(jscria), 2)
        
            cntbs = cntbs + 2
            cntba = cntba + 2
          endif

C           .. [alpha, beta] block

          if ( itype .eq. 3 ) then
            ind1 = jyat + cntat * ntt
            ind2 = jyaa + cntaa * ntt
            ind3 = ind1 + ntt
            ind4 = ind2 + ntt

C             .. scale the antisymmetric blocks by a factor of -1 ..

            call ascale (ntt, -1.0d0, v(ind2), v(ind2))
            call ascale (ntt, -1.0d0, v(ind4), v(ind4))

            call amove (ntt, v(ind1), v(jscrrs))
            call amove (ntt, v(ind2), v(jscrra))
            call amove (ntt, v(ind3), v(jscris))
            call amove (ntt, v(ind4), v(jscria))
        
            call sptblk (iout, iprint, v(jscr), lscr, nbasis, v(jzab),
     $           v(jscrrs), v(jscrra), v(jscris), v(jscria), 2)

            cntat = cntat + 2
            cntaa = cntaa + 2
          endif

C           .. [beta, alpha] block

          if ( itype .eq. 3 ) then
            ind1 = jybt + cntbt * ntt
            ind2 = jyba + cntba * ntt
            ind3 = ind1 + ntt
            ind4 = ind2 + ntt

C             .. scale the antisymmetric blocks by a factor of -1 ..

            call ascale (ntt, -1.0d0, v(ind2), v(ind2))
            call ascale (ntt, -1.0d0, v(ind4), v(ind4))

            call amove (ntt, v(ind1), v(jscrrs))
            call amove (ntt, v(ind2), v(jscrra))
            call amove (ntt, v(ind3), v(jscris))
            call amove (ntt, v(ind4), v(jscria))
        
            call sptblk (iout, iprint, v(jscr), lscr, nbasis, v(jzba),
     $           v(jscrrs), v(jscrra), v(jscris), v(jscria), 2)

            cntbt = cntbt + 2
            cntba = cntba + 2
          endif


C     Transform all the spin blocks to the orthonormal AO basis.

          call aobstf (iout, iprint, v(jscr), lscr, norbs, nbasis,
     $         irwx, junk, v(jzaa), v(jzab), v(jzba), v(jzbb), .true.,
     $         key, 2)


C     Collect all spin blocks into appropriate places in G(g).

C     Transform back into the original (NO) basis:
C       nos! . G(g) . nos

          if ( itype .eq. 1 ) then
            indx = (igrdt - 1) * 2*nosq

            call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos),
     $           norbs, v(jzaa), norbs, zero, v(jnost), norbs)
            call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost),
     $           norbs, v(jnos), norbs, zero, v(jzaa), norbs)

            call fileio (1, -irwgg, 2*nosq, v(jzaa), indx)

          elseif ( itype .eq. 2 ) then
            indxa = (igrdt - 1) * 4*nosq
            indxb = (igrdt - 1) * 4*nosq + 2*nosq

            call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos),
     $           norbs, v(jzaa), norbs, zero, v(jnost), norbs)
            call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost),
     $           norbs, v(jnos), norbs, zero, v(jzaa), norbs)

            call zgemm ('c', 'n', norbs, norbs, norbs, one,
     $           v(jnos+2*nosq), norbs, v(jzbb), norbs, zero,
     $           v(jnost), norbs)
            call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost),
     $           norbs, v(jnos+2*nosq), norbs, zero, v(jzbb), norbs)

            call fileio (1, -irwgg, 2*nosq, v(jzaa), indxa)
            call fileio (1, -irwgg, 2*nosq, v(jzbb), indxb)

          elseif ( itype .eq. 3 ) then
            indx = (igrdt - 1) * 8*nosq

C             .. combine spin blocks ..

            call dmblck (iout, iprint, v(jscr), lscr, norbs, v(jzf),
     $           v(jzaa), v(jzab), v(jzba), v(jzbb), 1, -1)

            call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $           v(jnos), 2*norbs, v(jzf), 2*norbs, zero, v(jnost),
     $           2*norbs)
            call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $           v(jnost), 2*norbs, v(jnos), 2*norbs, zero, v(jzf),
     $           2*norbs)

            call fileio (1, -irwgg, 8*nosq, v(jzf), indx)
          endif


C     Print matrices if requested.

          if ( iprint .gt. 4 ) then
            write (indc1, '(I5)') igrd
            write (indc2, '(I3)') n

            str = ' (' // indc1 // ',' // indc2 // ')'

            if ( itype .eq. 1 ) then
              write (iout, *) ' '
              call ltoutx (iout, 1, 2, 1, 0, 'G(g)' // str, norbs,
     $             norbs, v(jzaa), v(jzaa), 0)

            elseif ( itype .eq. 2 ) then
              write (iout, *) ' '
              call ltoutx (iout, 1, 3, 1, 0, 'G(g)' // str, norbs,
     $             norbs, v(jzaa), v(jzbb), 0)

            elseif ( itype .eq. 3 ) then
              write (iout, *) ' '
              call ltoutx (iout, 1, 2, 1, 0, 'G(g)' // str, 2*norbs,
     $             2*norbs, v(jzf), v(jzf), 0)
            endif
          endif

 103    continue
 101  continue


      return
      end


