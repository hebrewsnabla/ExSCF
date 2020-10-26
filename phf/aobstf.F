

      subroutine aobstf (iout, iprint, v, lenv, norbs, nbasis, irwx,
     $     x, paa, pab, pba, pbb, loadx, key, mode)

      implicit none


C +------------------------------------------------------------+
C |  aobstf   --   CAJH, 07.2011                               |
C |                                                            |
C |                                    (based on PHFB aobstf)  |
C |                                                            |
C |                                                            |
C |  Perform either of two basis set transformations           |
C |  (depends on mode) between the regular AO basis and an     |
C |  orthonormal AO basis:                                     |
C |                                                            |
C |    mode = 1,   X * P * X!                                  |
C |         = 2,   X! * P * X                                  |
C |                                                            |
C |  where P is the matrix to transform and X is the           |
C |  transformation matrix (=S^(-1/2)) stored in irwx.         |
C |                                                            |
C |  key controls which spin blocks are allocated and should   |
C |  be transformed:                                           |
C |                                                            |
C |    key = 1,  all spin blocks are active                    |
C |        = 2,  [aa] and [bb] blocks are active               |
C |        = 3,  [aa] block is active                          |
C |                                                            |
C |  In all cases (regardless of mode), all spin blocks        |
C |  should be allocated with size nbasis * nbasis.            |
C |                                                            |
C |  scratch space required:                                   |
C |                                                            |
C |    4*nbasis * norbs,  in all cases                         |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  Note that if loadx = .false., the matrix loaded into the  |
C |  x array has to have the same form as that from the RWF    |
C |  file (DO NOT MANIPULATE).                                 |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       norbs  - number of orbitals
C       nbasis - number of basis functions
C       irwx   - RWF file with transformation matrix (= S^(-1/2))
C                ( not used when loadx = .false. )
C       x      - transformation matrix (= S^(-1/2))
C                ( not used when loadx = .true. )
C       paa    - [alpha, alpha] block of P matrix [ in / out ]
C       pab    - [alpha, beta]  block of P matrix [ in / out ]
C       pba    - [beta, alpha]  block of P matrix [ in / out ]
C       pbb    - [beta, beta]   block of P matrix [ in / out ]
C       loadx  - whether to load transformation matrix from RWF file
C       key    - determines which blocks are non-zero
C       mode   - controls which transformation to perform:
C                  [ mode = 1,  X * P * X! ]
C                  [ mode = 2,  X! * P * X ]

      real*8      v(*), x(*), paa(*), pab(*), pba(*), pbb(*)
      integer     iout, iprint, lenv
      integer     norbs, nbasis, irwx, key, mode
      logical     loadx

C     other variables

      integer     jx, jscr, jend, mdv
      integer     iopcl, nri
      logical     doaa, doab, doba, dobb

C     constants

      complex*16  zero, one


      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)


C     Read iopcl from ILSW file.

C       iopcl = 0,  real RHF
C             = 1,  real UHF
C             = 2,  complex RHF
C             = 3,  complex UHF
C             = 6,  complex GHF (there is no real GHF)

      call ilsw (2, 1, iopcl)

      if ( iopcl .gt. 3 .and. iopcl .ne. 6 ) then
        call gauerr ('Incorrect iopcl in aobstf.')
      endif

      nri = 2

      if ( iopcl .le. 1 ) nri = 1


C     Memory allocation.

C       Allocate space to store:
C         - transformation matrix (jx)
C         - a scratch array of dimension nbasis * norbs

      jx   = 1
      jscr = jx   + 2*nbasis * norbs
      jend = jscr + 2*nbasis * norbs

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'aobstf')
      

C     Decipher key.

      doaa = .false.
      doab = .false.
      doba = .false.
      dobb = .false.

      if ( key .eq. 1 ) then
        doaa = .true.
        doab = .true.
        doba = .true.
        dobb = .true.

      elseif ( key .eq. 2 ) then
        doaa = .true.
        dobb = .true.

      elseif ( key .eq. 3 ) then
        doaa = .true.
      endif


C     Load transformation matrix.

      if ( loadx ) then
        call fileio (2, -irwx, nri*nbasis*norbs, v(jx), 0)
      else
        call amove (nri*nbasis*norbs, x(1), v(jx))
      endif

      if ( nri .eq. 1 ) then
        call zrecmp (3, v(jscr), v(jx), nbasis*norbs)
        call amove (2*nbasis*norbs, v(jscr), v(jx))
      endif


C     Perform transformation.

      if ( mode .eq. 1 ) then

        if ( doaa ) then
          call zgemm ('n', 'n', nbasis, norbs, norbs, one, v(jx),
     $         nbasis, paa, norbs, zero, v(jscr), nbasis)
          call zgemm ('n', 'c', nbasis, nbasis, norbs, one, v(jscr),
     $         nbasis, v(jx), nbasis, zero, paa, nbasis)
        endif

        if ( doab ) then
          call zgemm ('n', 'n', nbasis, norbs, norbs, one, v(jx),
     $         nbasis, pab, norbs, zero, v(jscr), nbasis)
          call zgemm ('n', 'c', nbasis, nbasis, norbs, one, v(jscr),
     $         nbasis, v(jx), nbasis, zero, pab, nbasis)
        endif

        if ( doba ) then
          call zgemm ('n', 'n', nbasis, norbs, norbs, one, v(jx),
     $         nbasis, pba, norbs, zero, v(jscr), nbasis)
          call zgemm ('n', 'c', nbasis, nbasis, norbs, one, v(jscr),
     $         nbasis, v(jx), nbasis, zero, pba, nbasis)
        endif

        if ( dobb ) then
          call zgemm ('n', 'n', nbasis, norbs, norbs, one, v(jx),
     $         nbasis, pbb, norbs, zero, v(jscr), nbasis)
          call zgemm ('n', 'c', nbasis, nbasis, norbs, one, v(jscr),
     $         nbasis, v(jx), nbasis, zero, pbb, nbasis)
        endif

      elseif ( mode .eq. 2 ) then

        if ( doaa ) then
          call zgemm ('n', 'n', nbasis, norbs, nbasis, one, paa,
     $         nbasis, v(jx), nbasis, zero, v(jscr), nbasis)
          call zgemm ('c', 'n', norbs, norbs, nbasis, one, v(jx),
     $         nbasis, v(jscr), nbasis, zero, paa, norbs)
        endif

        if ( doab ) then
          call zgemm ('n', 'n', nbasis, norbs, nbasis, one, pab,
     $         nbasis, v(jx), nbasis, zero, v(jscr), nbasis)
          call zgemm ('c', 'n', norbs, norbs, nbasis, one, v(jx),
     $         nbasis, v(jscr), nbasis, zero, pab, norbs)
        endif

        if ( doba ) then
          call zgemm ('n', 'n', nbasis, norbs, nbasis, one, pba,
     $         nbasis, v(jx), nbasis, zero, v(jscr), nbasis)
          call zgemm ('c', 'n', norbs, norbs, nbasis, one, v(jx),
     $         nbasis, v(jscr), nbasis, zero, pba, norbs)
        endif

        if ( dobb ) then
          call zgemm ('n', 'n', nbasis, norbs, nbasis, one, pbb,
     $         nbasis, v(jx), nbasis, zero, v(jscr), nbasis)
          call zgemm ('c', 'n', norbs, norbs, nbasis, one, v(jx),
     $         nbasis, v(jscr), nbasis, zero, pbb, norbs)
        endif

      endif


      return
      end


