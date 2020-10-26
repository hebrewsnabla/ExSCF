

      subroutine fockdg (iout, iprint, v, lenv, iphf, norbs, na, nb,
     $     ivshft, vshft, gap, f, p, irwvec, irwval, irwp)

      implicit none


C +------------------------------------------------------------+
C |  fockdg   --   CAJH, 07.2011                               |
C |                                                            |
C |                                                            |
C |  Diagonalize the effective Fock matrix associated with     |
C |  the PHF formalism to obtain a new set of molecular        |
C |  orbitals (that define the new reference HF). We then      |
C |  construct the density matrix associated with the          |
C |  reference HF state.                                       |
C |                                                            |
C |  On input, F is expected in lower triangular form, ready   |
C |  for diagonalization. On output, we update P with the      |
C |  new density matrix and save the eigenvectors in irwvec,   |
C |  the Fock eigenvalues in irwval, and the new density       |
C |  matrix in irwp.                                           |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       norbs  - number of orbitals
C       na     - number of alpha electrons
C       nb     - number of beta electrons
C       ivshft - level shifting control
C       vshft  - level shifting applied
C       gap    - HOMO-LUMO difference in Fock eigenvalues
C                (without level shifting)  [ out ]
C       f      - effective Fock matrix
C       p      - HF density matrix  [ in / out ]
C       irwvec - RWF file with orbital coefficients
C                (orthonormal AO basis)  [ updated here ]
C       irwval - RWF file with Fock eigenvalues  [ updated here ]
C       irwp   - RWF file with HF density matrix  [ updated here ]

      real*8      v(*), f(*), p(*), gap
      real*8      vshft, gap
      integer     iout, iprint, lenv, iphf, norbs, na, nb
      integer     ivshft, irwvec, irwval, irwp

C     other variables

      real*8      thresh, gapa, gapb
      integer     nel, ntto, ntt2o, nosq, itype, info
      integer     jf, jp, jvec, jval, jscr1, jscr2, jend
      logical     icmplx, ispuhf, ispghf

      parameter   ( thresh = 5.0d-15 )



C     Construct useful variables.

      nel = na + nb

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
        call gauerr ('Unrecognized iphf in fockdg.')
      endif


C     Memory allocation.

C       Allocate space for:
C         - scratch copy of Fock matrix (jf)
C         - scratch array for density matrix, in square format (jp)
C         - eigenvectors of Fock matrix (jvec)
C         - eigenvalues of Fock matrix (jval)
C         - scratch arrays needed for zhpev (jscr1, jscr2)

      if ( itype .eq. 1 ) then
        jf    = 1
        jp    = jf    + 2*ntto
        jvec  = jp    + 2*nosq
        jval  = jvec  + 2*nosq
        jscr1 = jval  + norbs
        jscr2 = jscr1 + 2*(2*norbs - 1)
        jend  = jscr2 + 3*norbs - 2

      elseif ( itype .eq. 2 .or. ispuhf ) then
        jf    = 1
        jp    = jf    + 4*ntto
        jvec  = jp    + 4*nosq
        jval  = jvec  + 4*nosq
        jscr1 = jval  + 2*norbs
        jscr2 = jscr1 + 2*(2*norbs - 1)
        jend  = jscr2 + 3*norbs - 2

      elseif ( itype .eq. 3 ) then
        jf    = 1
        jp    = jf    + 2*ntt2o
        jvec  = jp    + 8*nosq
        jval  = jvec  + 8*nosq
        jscr1 = jval  + 2*norbs
        jscr2 = jscr1 + 2*(4*norbs - 1)
        jend  = jscr2 + 6*norbs - 2
      endif
      
      call tstcor (jend-1, lenv, 'solvci')
      

C     Move Fock matrix to scratch array.

      if ( itype .eq. 1 ) then
        call amove (2*ntto, f, v(jf))
      elseif ( itype .eq. 2 .or. ispuhf ) then
        call amove (4*ntto, f, v(jf))
      elseif ( itype .eq. 3 ) then
        call amove (2*ntt2o, f, v(jf))
      endif


C     Diagonalize Fock matrix.
C     Take the complex conjugate of the Fock matrix (so that we are
C     properly sending the upper triangle to LaPack).

      if ( itype .eq. 1 ) then
        call vconjg (ntto, v(jf))

        call zhpev ('v', 'u', norbs, v(jf), v(jval), v(jvec), norbs,
     $       v(jscr1), v(jscr2), info)

        if ( info .ne. 0 ) then
          call gauerr ('zhpev failed in fockdg.')
        endif

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call vconjg (2*ntto, v(jf))

        call zhpev ('v', 'u', norbs, v(jf), v(jval), v(jvec), norbs,
     $       v(jscr1), v(jscr2), info)

        if ( info .ne. 0 ) then
          call gauerr ('zhpev failed in fockdg.')
        endif

        if ( nb .ge. 1 ) then
        call zhpev ('v', 'u', norbs, v(jf+2*ntto), v(jval+norbs),
     $       v(jvec+2*nosq), norbs, v(jscr1), v(jscr2), info)

        if ( info .ne. 0 ) then
          call gauerr ('zhpev failed in fockdg.')
        endif
        else
          call aclear (2*nosq, v(jvec+2*nosq))
        endif

      elseif ( itype .eq. 3 ) then
        call vconjg (ntt2o, v(jf))

        call zhpev ('v', 'u', 2*norbs, v(jf), v(jval), v(jvec), 2*norbs,
     $       v(jscr1), v(jscr2), info)

        if ( info .ne. 0 ) then
          call gauerr ('zhpev failed in fockdg.')
        endif
      endif


C     Generate new density matrix.

      if ( itype .eq. 1 ) then
        call zherk ('u', 'n', norbs, na, 1.0d0, v(jvec), norbs, 0.0d0,
     $       v(jp), norbs)

        call zlinear (v(jp), p, norbs, norbs)
        call vconjg (ntto, p)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call zherk ('u', 'n', norbs, na, 1.0d0, v(jvec), norbs, 0.0d0,
     $       v(jp), norbs)

        if ( nb .ge. 1 ) then
        call zherk ('u', 'n', norbs, nb, 1.0d0, v(jvec+2*nosq), norbs,
     $       0.0d0, v(jp+2*nosq), norbs)
        else
          call aclear (2*nosq, v(jp+2*nosq))
        endif

        call zlinear (v(jp), p, norbs, norbs)
        call zlinear (v(jp+2*nosq), p(1+2*ntto), norbs, norbs)
        call vconjg (2*ntto, p)
c        write (iout, *) ' '
c        call ltoutx (iout, 1, 3, 1, 0, 'New P',
c     $        norbs, norbs, v(jp), v(jp+2*nosq), 0)
        write (iout, *) ' '
        call ltoutx (iout, 0, 3, 1, 0, 'New P',
     $        norbs, norbs, p, p(1+2*ntto), 0)

      elseif ( itype .eq. 3 ) then
        call zherk ('u', 'n', 2*norbs, nel, 1.0d0, v(jvec), 2*norbs,
     $       0.0d0, v(jp), 2*norbs)

        call zlinear (v(jp), p, 2*norbs, 2*norbs)
        call vconjg (ntt2o, p)
      endif


C     Clean-up P by eliminating very small elements.

      if ( itype .eq. 1 ) then
        call rndvec (2*ntto, p, thresh)
      elseif ( itype .eq. 2 .or. ispuhf ) then
        call rndvec (4*ntto, p, thresh)
      elseif ( itype .eq. 3 ) then
        call rndvec (2*ntt2o, p, thresh)
      endif


C     Subtract level shifting from occupied orbitals.

c$$$      if ( ivshft .ne. -1 ) then
c$$$        if ( itype .eq. 1 ) then
c$$$          call ainc (na, vshft, v(jval), v(jval))
c$$$
c$$$        elseif ( itype .eq. 2 .or. ispuhf ) then
c$$$          call ainc (na, vshft, v(jval), v(jval))
c$$$          call ainc (nb, vshft, v(jval+norbs), v(jval+norbs))
c$$$
c$$$        elseif ( itype .eq. 3 ) then
c$$$          call ainc (nel, vshft, v(jval), v(jval))
c$$$        endif
c$$$      endif


C     Determine gap.

      if ( itype .eq. 1 ) then
        gap = v(jval+na) - v(jval+na-1)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        gapa = v(jval+na) - v(jval+na-1)
        gapb = v(jval+norbs+nb) - v(jval+norbs+nb-1)

        gap = min (gapa, gapb)

      elseif ( itype .eq. 3 ) then
        gap = v(jval+nel) - v(jval+nel-1)
      endif


C     Print Fock eigenvalues.

      if ( iprint .gt. 1 ) then
        if ( itype .eq. 1 ) then
          call prtocc (iout, 'Fock eigenvalues', 'orbital',
     $         'energy', norbs, v(jval))

        elseif ( itype .eq. 2 .or. ispuhf ) then
          call prtocc (iout, 'Fock eigenvalues (alpha)', 'orbital',
     $         'energy', norbs, v(jval))
          call prtocc (iout, 'Fock eigenvalues (beta)', 'orbital',
     $         'energy', norbs, v(jval+norbs))

        elseif ( itype .eq. 3 ) then
          call prtocc (iout, 'Fock eigenvalues', 'orbital',
     $         'energy', 2*norbs, v(jval))
        endif
      endif


C     Save arrays into appropriate RWF files.

      if ( itype .eq. 1 ) then
        call fileio (1, -irwvec, 2*nosq, v(jvec), 0)
        call fileio (1, -irwval, norbs, v(jval), 0)
        call fileio (1, -irwp, 2*ntto, p, 0)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call fileio (1, -irwvec, 4*nosq, v(jvec), 0)
        call fileio (1, -irwval, 2*norbs, v(jval), 0)
        call fileio (1, -irwp, 4*ntto, p, 0)

      elseif ( itype .eq. 3 ) then
        call fileio (1, -irwvec, 8*nosq, v(jvec), 0)
        call fileio (1, -irwval, 2*norbs, v(jval), 0)
        call fileio (1, -irwp, 2*ntt2o, p, 0)
      endif


      return
      end


