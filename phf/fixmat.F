

      subroutine fixmat (iout, iprint, v, lenv, iopcl, nbasis, norbs,
     $     irws, irwx, irwh, irws1, irwx1, irwh1)

      implicit none


C +------------------------------------------------------------+
C |  fixmat   --   CAJH, 07.2011                               |
C |                                                            |
C |                                                            |
C |  Simple routine that creates working versions (for the     |
C |  PHF code) of the overlap, transformation (=S^(-1/2)),     |
C |  and core Hamiltonian matrices.                            |
C |                                                            |
C |  In particular, the PHF subroutines all assume that the    |
C |  RWF files storing this matrices keep only the alpha-      |
C |  alpha block (there is no need to store the full matrix).  |
C |  However, when the user requests for a GHF calculation,    |
C |  Gaussian stores these matrices in full form (with all     |
C |  spin blocks present).                                     |
C |                                                            |
C |  fixmat thus creates copies of these matrices and creates  |
C |  RWF files (irws1, irwx1, irwh1) where such matrices are   |
C |  stored in the expected form.                              |
C |                                                            |
C |  Note that if the user requested an RHF or a UHF           |
C |  wavefunction, then there is no need to do anything.       |
C |                                                            |
C |  The variable iopcl is used to determine how Gaussian is   |
C |  storing these matrices.                                   |
C |                                                            |
C |     iopcl = 0,  real RHF                                   |
C |           = 1,  real UHF                                   |
C |           = 2,  complex RHF                                |
C |           = 3,  complex UHF                                |
C |           = 6,  complex GHF (there is no real GHF)         |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iopcl  - describes type of HF wavefunction
C       nbasis - number of basis functions
C       norbs  - number of orbitals
C       irws   - RWF file with overlap matrix
C       irwx   - RWF file with transformation matrix (=S^(-1/2))
C       irwh   - RWF file with core Hamiltonian matrix
C       irws1  - working version of irws [ out ]
C       irwx1  - working version of irwx [ out ]
C       irwh1  - working version of irwh [ out ]

      real*8      v(*)
      integer     iout, iprint, lenv, iopcl, nbasis, norbs
      integer     irws, irwx, irwh, irws1, irwx1, irwh1

C     other variables

      integer     ntt, ntt2, i, j, ind1, ind2
      integer     js1, js2, js3, jh1, jh2, jh3
      integer     jx1, jx2, jscr, jend

C     functions

      integer     igetrw



C     Return if iopcl <= 3.

      if ( iopcl .le. 3 ) then
        irws1 = irws
        irwx1 = irwx
        irwh1 = irwh

        return
      endif


C     Construct useful variables.

      ntt = nbasis * (nbasis + 1) / 2
      ntt2 = nbasis * (2*nbasis + 1)


C     Create scratch files for matrices of the correct dimension.

      irws1 = igetrw (2*ntt, 0)
      irwx1 = igetrw (2*nbasis*norbs, 0)
      irwh1 = igetrw (2*ntt, 0)


C     Memory allocation.

C       Allocate space for:
C         - overlap matrix (full, small, lt: js1, js2, js3)
C         - core Hamiltonian matrix (full, small, lt: jh1, jh2, jh3)
C         - transformation matrix (full, small: jx1, jx2)
C         - scratch array to reorganize columns (jscr)

      js1  = 1
      js2  = js1  + 8*nbasis*nbasis
      js3  = js2  + 2*nbasis*nbasis
      jh1  = js3  + 2*ntt
      jh2  = jh1  + 8*nbasis*nbasis
      jh3  = jh2  + 2*nbasis*nbasis
      jx1  = jh3  + 2*ntt
      jx2  = jx1  + 8*nbasis*norbs
      jscr = jx2  + 2*nbasis*norbs
      jend = jscr + 8*nbasis*nbasis

      call tstcor (jend-1, lenv, 'fixmat')


C     Extract single spin-block from full overlap matrix.

C       .. recover overlap matrix from file

      call fileio (2, -irws, 2*ntt2, v(js1), 0)

C       .. expand to square matrix

      call zsquare (v(js1), v(js1), 2*nbasis, 2*nbasis, 0)


C     Gaussian stores matrices by ordering basis functions as
C       1a, 1b, 2a, 2b, ...

C     On the other hand, we want to have them organized as
C       1a, 2a, ..., 1b, 2b, ...

C     This is so that our matrices have a spin-block structure:

C          S = ( S_aa  S_ab )
C              ( S_ba  S_bb )

C     We thus need to reorganize columns and rows.

C       .. first reorganize columns

      do 10 i = 1, 2*nbasis
        if ( mod (i,2) .eq. 1 ) then
          j = (i-1)/2 + 1
        elseif ( mod (i,2) .eq. 0 ) then
          j = i/2 + nbasis
        endif

        ind1 = (i - 1) * 4*nbasis + js1
        ind2 = (j - 1) * 4*nbasis + jscr

        call amove (4*nbasis, v(ind1), v(ind2))
 10   continue

C       .. now reorganize rows

      call trspnc (.false., 2*nbasis, 2*nbasis, v(jscr))

      do 11 i = 1, 2*nbasis
        if ( mod (i,2) .eq. 1 ) then
          j = (i-1)/2 + 1
        elseif ( mod (i,2) .eq. 0 ) then
          j = i/2 + nbasis
        endif

        ind1 = (i - 1) * 4*nbasis + jscr
        ind2 = (j - 1) * 4*nbasis + js1

        call amove (4*nbasis, v(ind1), v(ind2))
 11   continue

      call trspnc (.false., 2*nbasis, 2*nbasis, v(js1))

C       .. obtain the alpha-alpha block of the matrix

      call sptmat (2, 2*nbasis, 2*nbasis, v(js1), nbasis, nbasis,
     $     v(js2), 1)

C       .. transform to lower triangular

      call zlinear (v(js2), v(js3), nbasis, nbasis)

C       .. save to scratch file

      call fileio (1, -irws1, 2*ntt, v(js3), 0)


C     debugging...
C       Print overlap matrix.

      write (iout, *) ' '
      call ltoutx (iout, 1, 2, 1, 0, 'Overlap (full)', 2*nbasis,
     $     2*nbasis, v(js1), v(js1), 0)

      write (iout, *) ' '
      call ltoutx (iout, 0, 2, 1, 0, 'Overlap (small)', nbasis, nbasis,
     $     v(js3), v(js3), 0)



C     Extract single spin-block from full core Hamiltonian matrix in the
C     same way as was done for the overlap matrix.

      call fileio (2, -irwh, 2*ntt2, v(jh1), 0)

      call zsquare (v(jh1), v(jh1), 2*nbasis, 2*nbasis, 0)

      do 20 i = 1, 2*nbasis
        if ( mod (i,2) .eq. 1 ) then
          j = (i-1)/2 + 1
        elseif ( mod (i,2) .eq. 0 ) then
          j = i/2 + nbasis
        endif

        ind1 = (i - 1) * 4*nbasis + jh1
        ind2 = (j - 1) * 4*nbasis + jscr

        call amove (4*nbasis, v(ind1), v(ind2))
 20   continue

      call trspnc (.false., 2*nbasis, 2*nbasis, v(jscr))

      do 21 i = 1, 2*nbasis
        if ( mod (i,2) .eq. 1 ) then
          j = (i-1)/2 + 1
        elseif ( mod (i,2) .eq. 0 ) then
          j = i/2 + nbasis
        endif

        ind1 = (i - 1) * 4*nbasis + jscr
        ind2 = (j - 1) * 4*nbasis + jh1

        call amove (4*nbasis, v(ind1), v(ind2))
 21   continue

      call trspnc (.false., 2*nbasis, 2*nbasis, v(jh1))

      call sptmat (2, 2*nbasis, 2*nbasis, v(jh1), nbasis, nbasis,
     $     v(jh2), 1)

      call zlinear (v(jh2), v(jh3), nbasis, nbasis)

      call fileio (1, -irwh1, 2*ntt, v(jh3), 0)


C     debugging...
C       Print core Hamiltonian matrix.

      write (iout, *) 'Start itoutx '
      call ltoutx (iout, 1, 2, 1, 1, 'Core Hamiltonian (full)',
     $     2*nbasis, 2*nbasis, v(jh1), v(jh1), 0)

      write (iout, *) ' '
      call ltoutx (iout, 0, 2, 1, 1, 'Core Hamiltonian (small)',
     $     nbasis, nbasis, v(jh3), v(jh3), 0)



C     Extract single spin-block from full transformation matrix in the
C     same way as was done for the overlap matrix.

      call fileio (2, -irwx, 8*nbasis*norbs, v(jx1), 0)

      do 30 i = 1, 2*norbs
        if ( mod (i,2) .eq. 1 ) then
          j = (i-1)/2 + 1
        elseif ( mod (i,2) .eq. 0 ) then
          j = i/2 + norbs
        endif

        ind1 = (i - 1) * 4*nbasis + jx1
        ind2 = (j - 1) * 4*nbasis + jscr

        call amove (4*nbasis, v(ind1), v(ind2))
 30   continue

      call trspnc (.false., 2*nbasis, 2*norbs, v(jscr))

      do 31 i = 1, 2*nbasis
        if ( mod (i,2) .eq. 1 ) then
          j = (i-1)/2 + 1
        elseif ( mod (i,2) .eq. 0 ) then
          j = i/2 + nbasis
        endif

        ind1 = (i - 1) * 4*norbs + jscr
        ind2 = (j - 1) * 4*norbs + jx1

        call amove (4*norbs, v(ind1), v(ind2))
 31   continue

      call trspnc (.false., 2*norbs, 2*nbasis, v(jx1))

      call sptmat (2, 2*nbasis, 2*norbs, v(jx1), nbasis, norbs,
     $     v(jx2), 1)

      call fileio (1, -irwx1, 2*nbasis*norbs, v(jx2), 0)


C     debugging...
C       Print transformation matrix.

      write (iout, *) ' '
      call ltoutx (iout, 3, 2, 1, 0, 'Transformation (full)', 2*nbasis,
     $     2*norbs, v(jx1), v(jx1), 0)

      write (iout, *) ' '
      call ltoutx (iout, 3, 2, 1, 0, 'Transformation (small)', nbasis,
     $     norbs, v(jx2), v(jx2), 0)

      write (iout, *) 'End fixmat '

      return
      end


