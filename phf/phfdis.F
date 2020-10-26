

      subroutine phfdis (iout, iprint, v, lenv, iphf, norbs, it, idisgo,
     $     istrat, rises, maxcom, irwpa, irwpb, irwfa, irwfb, errmax,
     $     nersav, indfok, energy, encyc, errcyc, nsvdis, iobdis,
     $     iopadi, iopbdi, iofadi, iofbdi, ioerdi)

      implicit none


C +------------------------------------------------------------+
C |  phfdis   --   CAJH, 07.2011                               |
C |                                                            |
C |                                   (based on PHFB phfbdis)  |
C |                                                            |
C |                                                            |
C |  Calls the routine SCFDIS to perform the extrapolation     |
C |  of the effective Fock matrix F using the direct           |
C |  inversion of iterative subspace (DIIS) algorithm.         |
C |                                                            |
C |  Convergence of the PHF equations is achieved when         |
C |                                                            |
C |    [ F, P ] = 0,                                           |
C |                                                            |
C |  where F is the effective Fock matrix and P is the         |
C |  density matrix associated with the HF state.              |
C |                                                            |
C |  DIIS extrapolates the effective Fock matrix F to reduce   |
C |  the error in the commutator. (SCFDIS should always be     |
C |  asked to do commutator DIIS here!)                        |
C |                                                            |
C +------------------------------------------------------------+


C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       norbs  - number of orbitals
C       it     - current SCF iteration
C       idisgo - SCF iteration at which DIIS should start
C       istrat - strategy to follow in SCF procedure
C       rises  - whether current energy is higher than previous cycle
C       maxcom - number of DIIS vectors
C       irwpX  - RWF file with HF density matrix (X spin)
C       irwfX  - RWF file with effective Fock matrix (X spin)
C       errmax - array with maximum DIIS error
C       nersav - number of errors saved
C       indfok - pointers to saved Fock matrices on disk
C       energy - energy computed for current HF density matrix
C       encyc  - queue of saved energies
C       errcyc - queue of saved errors
C       nsvdis - number of Fock matrices currently stored (do not touch)
C       iobdis, iopadi, iopbdi, iofadi, iofbdi, ioerdi - DIIS controls

      integer     iout, iprint, lenv, iphf
      integer     norbs, it, idisgo, istrat, maxcom
      integer     irwpa, irwpb, irwfa, irwfb
      integer     nersav, indfok(*), nsvdis
      integer     iobdis, iopadi, iopbdi, iofadi, iofbdi, ioerdi
      logical     rises

      real*8      v(*), errmax(0:nersav), encyc(*), errcyc(*), energy

C     other variables

      integer     itype, n, jend, mdv
      integer     lscr1, lscr2, lscr3, lscr4, lscr5, lscr6
      integer     jscr1, jscr2, jscr3, jscr4, jscr5, jscr6
      integer     ifsel, niter, numfs, idityp, iopcl
      integer     irwpxa, irwpxb, irwfxa, irwfxb
      integer     irws, irwv, irwerr, iprtf
      logical     icmplx, ispuhf, ispghf


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
        call gauerr ('Unrecognized iphf in phfdis.')
      endif


C     Memory allocation.

C       Specify dimensions of scratch arrays (see SCFDIS for details).
C       Be careful with the fact that we are sending complex matrices!

      if ( itype .lt. 3 .or. ispuhf ) then
        n = max (norbs, maxcom+1)
      elseif ( itype .eq. 3 ) then
        n = max (2*norbs, maxcom+1)
      endif

      lscr1 = 2*n*(n+1)
      lscr2 = 2*n*n
      lscr4 = (maxcom+1)**2
      lscr5 = 7*(maxcom+1)
      lscr6 = 2*n*n

      if ( itype .lt. 3 .or. ispuhf ) then
        lscr3 = max (norbs, (maxcom+1)**2)
      elseif ( itype .eq. 3 ) then
        lscr3 = max (2*norbs, (maxcom+1)**2)
      endif


C       Allocate space for scratch arrays.

      jscr1 = 1
      jscr2 = jscr1 + lscr1
      jscr3 = jscr2 + lscr2
      jscr4 = jscr3 + lscr3
      jscr5 = jscr4 + lscr4
      jscr6 = jscr5 + lscr5
      jend  = jscr6 + lscr6


      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'phfdis')


C     Set variable ifsel (see SCFDIS for details).

      if ( istrat .ne. 3 ) then
        if ( it .ge. idisgo ) then
          ifsel = 0
        else
          ifsel = 3
        endif

      elseif ( istrat .eq. 3 ) then
        if ( mod (it,maxcom) .eq. 0 ) then
          ifsel = 0
        else
          ifsel = 3
        endif
      endif

      if ( istrat .eq. 2 .and. ifsel .eq. 0 .and. rises ) ifsel = 5


C     Set variables in preparation for SCFDIS.

      niter = it
      numfs = 1
      idityp = 1

      if ( itype .eq. 1 ) then
        iopcl = 2
      elseif ( itype .eq. 2 .or. ispuhf ) then
        iopcl = 3
      elseif ( itype .eq. 3 ) then
        iopcl = 6
      endif

      irwpxa = 0
      irwpxb = 0
      irwfxa = 0
      irwfxb = 0

      irws = 0
      irwv = 0
      irwerr = 0

      iprtf = iprint
      if ( iprtf .lt. 0 ) iprtf = 0


      call scfdis (iout, iprtf, idityp, ifsel, maxcom, niter, numfs,
     $     iopcl, norbs, norbs, irwpa, irwpxa, irwfa, irwfxa, irwpb,
     $     irwpxb, irwfb, irwfxb, irws, irwv, irwerr, nsvdis, iobdis,
     $     iopadi, iopbdi, iofadi, iofbdi, ioerdi, 0.0d0, v(jscr1),
     $     v(jscr2), v(jscr3), v(jscr4), v(jscr5), errmax, indfok,
     $     energy, encyc, errcyc, v(jscr6), lscr6)


      return
      end


