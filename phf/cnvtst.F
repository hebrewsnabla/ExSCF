

      subroutine cnvtst (iout, iprint, v, lenv, iphf, itst, acurcy,
     $     norbs, irwp, irwp1, rmsdp, maxp, delen, done)

      implicit none


C +------------------------------------------------------------+
C |  cnvtst   ---   CAJH, 07.2011                              |
C |                                                            |
C |                                    (based on PHFB cnvtst)  |
C |                                                            |
C |                                                            |
C |  Test whether convergence in the SCF cycle has been        |
C |  achieved. The following quantities are computed:          |
C |                                                            |
C |    rmsdp - RMSD of HF density matrix elements              |
C |    maxp  - maximum change in HF density matrix             |
C |                                                            |
C |  The quantity itst controls which variable should be       |
C |  used to determine convergence:                            |
C |                                                            |
C |    itst = 1,  check for convergence in energy              |
C |         = 2,  check for convergence in rmsdp               |
C |         = 3,  check for convergence in maxp                |
C |                                                            |
C |  This quantity is set by iop(13). See phfdrv for details.  |
C |                                                            |
C |  The convergence criterion is stored in the variable       |
C |  acurcy, which is in turn set by iop(14).                  |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       itst   - determines what should be used to test for convergence
C       acurcy - accuracy desired in the SCF procedure
C       norbs  - number of orbitals
C       irwp   - RWF file with new HF density matrix
C       irwp1  - RWF file with old HF density matrix
C       rmsdp  - RMSD of density matrix [ out ]
C       maxp   - maximum change in density matrix [ out ]
C       delen  - change in energy
C       done   - whether convergence has been achieved [ out ]

      real*8      v(*), acurcy
      real*8      rmsdp, maxp, delen
      integer     iout, iprint, lenv, iphf, itst, norbs
      integer     irwp, irwp1
      logical     done

C     other variables

      real*8      trp
      integer     ntto, ntt2o, nosq, itype, szp
      integer     jpnew, jpold, mdv, jend
      integer     jscr1, jscr2, jscr3, jpsq1, jpsq2
      logical     icmplx, ispuhf, ispghf

C     functions

      complex*16  trcabc
      real*8      armaxc


C     Construct useful variables.

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
        call gauerr ('Unrecognized iphf in cnvtst.')
      endif


C     Define size of density matrix.

      if ( itype .eq. 1 ) then
        szp = 2*ntto
      elseif ( itype .eq. 2 .or. ispuhf ) then
        szp = 4*ntto
      elseif ( itype .eq. 3 ) then
        szp = 2*ntt2o
      endif


C     Memory allocation.

C       Allocate space for old and new density matrices.

      jpnew = 1
      jpold = jpnew + szp
      jend  = jpold + szp

C       Allocate space for:
C         - one scratch array of dimension lower triangular complex
C         - two scratch arrays of dimension square

      if ( itype .lt. 3 .or. ispuhf ) then
        jscr1 = jend
        jscr2 = jscr1 + 2*ntto
        jscr3 = jscr2 + nosq
        jend  = jscr3 + nosq

      elseif ( itype .eq. 3 ) then
        jscr1 = jend
        jscr2 = jscr1 + 2*ntt2o
        jscr3 = jscr2 + 4*nosq
        jend  = jscr3 + 4*nosq
      endif

C       Allocate space for two square copies of the density matrix.

      if ( itype .eq. 1 ) then
        jpsq1 = jend
        jpsq2 = jpsq1 + 2*nosq
        jend  = jpsq2 + 2*nosq

      elseif ( itype .eq. 2 .or. ispuhf ) then
        jpsq1 = jend
        jpsq2 = jpsq1 + 4*nosq
        jend  = jpsq2 + 4*nosq

      elseif ( itype .eq. 3 ) then
        jpsq1 = jend
        jpsq2 = jpsq1 + 8*nosq
        jend  = jpsq2 + 8*nosq
      endif

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'cnvtst')



C     Load old and new density matrices.

      call fileio (2, -irwp, szp, v(jpnew), 0)
      call fileio (2, -irwp1, szp, v(jpold), 0)


C     Subtract old density matrix from the new one.

      call asub (szp, v(jpnew), v(jpold), v(jpnew))


C     Determine maximum change in density matrix.

      maxp = armaxc (2, v(jpnew), szp/2)


C     Expand density matrix to square.

      if ( itype .eq. 1 ) then
        call zrecmp (1, v(jpnew), v(jscr1), ntto)
        call square (v(jscr1), v(jscr2), norbs, norbs, 0)
        call square (v(jscr1+ntto), v(jscr3), norbs, norbs, -1)
        call zrecmp (2, v(jpsq1), v(jscr2), nosq)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call zrecmp (1, v(jpnew), v(jscr1), ntto)
        call square (v(jscr1), v(jscr2), norbs, norbs, 0)
        call square (v(jscr1+ntto), v(jscr3), norbs, norbs, -1)
        call zrecmp (2, v(jpsq1), v(jscr2), nosq)

        call zrecmp (1, v(jpnew+2*ntto), v(jscr1), ntto)
        call square (v(jscr1), v(jscr2), norbs, norbs, 0)
        call square (v(jscr1+ntto), v(jscr3), norbs, norbs, -1)
        call zrecmp (2, v(jpsq1+2*nosq), v(jscr2), nosq)

      elseif ( itype .eq. 3 ) then
        call zrecmp (1, v(jpnew), v(jscr1), ntt2o)
        call square (v(jscr1), v(jscr2), 2*norbs, 2*norbs, 0)
        call square (v(jscr1+ntt2o), v(jscr3), 2*norbs, 2*norbs, -1)
        call zrecmp (2, v(jpsq1), v(jscr2), 4*nosq)
      endif


C     Copy density matrix to a second array.
C     Tranpose the second array.

      if ( itype .eq. 1 ) then
        call amove (2*nosq, v(jpsq1), v(jpsq2))

        call trspnc (.false., norbs, norbs, v(jpsq2))

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call amove (4*nosq, v(jpsq1), v(jpsq2))

        call trspnc (.false., norbs, norbs, v(jpsq2))
        call trspnc (.false., norbs, norbs, v(jpsq2+2*nosq))

      elseif ( itype .eq. 3 ) then
        call amove (8*nosq, v(jpsq1), v(jpsq2))

        call trspnc (.false., 2*norbs, 2*norbs, v(jpsq2))
      endif


C     Compute the sum:

C       trp = trace [ (Pnew - Pold) * (Pnew - Pold)' ]

      if ( itype .eq. 1 ) then
        trp = dble (trcabc (1, norbs, v(jpsq1), v(jpsq2)))

      elseif ( itype .eq. 2 .or. ispuhf ) then
        trp = dble (trcabc (1, norbs, v(jpsq1), v(jpsq2)))
        trp = trp +
     $    dble (trcabc (1, norbs, v(jpsq1+2*nosq), v(jpsq2+2*nosq)))

      elseif ( itype .eq. 3 ) then
        trp = dble (trcabc (1, 2*norbs, v(jpsq1), v(jpsq2)))
      endif


C     Compute rmsdp.

      if ( itype .eq. 1 ) then
        rmsdp = sqrt (trp / dble (nosq))

      elseif ( itype .eq. 2 .or. ispuhf ) then
        rmsdp = sqrt (trp / dble (2*nosq))

      elseif ( itype .eq. 3 ) then
        rmsdp = sqrt (trp / dble (4*nosq))
      endif


C     Decide whether convergence has been achieved.

      done = .false.

      if ( itst .eq. 1 ) then
        if ( delen .le. 0.0d0 .and. -delen .lt. acurcy ) done = .true.

      elseif ( itst .eq. 2 ) then
        if ( rmsdp .lt. acurcy ) done = .true.

      elseif ( itst .eq. 3 ) then
        if ( maxp .lt. acurcy ) done = .true.
      endif


      return
      end


