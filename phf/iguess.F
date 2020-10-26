

      subroutine iguess (iout, iprint, v, lenv, iphf, iopcl, na, nb,
     $     nbasis, norbs, irwx, irws, irwca, irwcb, irwvec, irwp, p,
     $     igues1)

      implicit none


C +------------------------------------------------------------+
C |  iguess   --   CAJH, 07.2011                               |
C |                                                            |
C |                                                            |
C |  Prepare initial guess for PHF calculations. The igues1    |
C |  variable controls how this initialization takes place.    |
C |                                                            |
C |    igues1 = 1,  use GDV initial guess                      |
C |           = 2,  use GDV initial guess +                    |
C |                 apply complex factor to some elements      |
C |                 of HOMO to break c.c. symmetry **          |
C |                                                            |
C |    ** for UHF initial guesses, the phase factor is         |
C |       applied to both the alpha HOMO and the beta HOMO     |
C |                                                            |
C |                                                            |
C |  igues1 reads GDV RWF files with MO coefficients (which    |
C |  are stored in the regular AO basis). It transforms these  |
C |  sets of coefficients into an orthonormal AO basis.        |
C |                                                            |
C |  It then transforms the orbitals according to igues1,      |
C |  and saves the resulting set of MO coefficients in the     |
C |  appropriate RWF file, irwvec. Using these MO              |
C |  coefficients, the density matrix is constructed in the    |
C |  orthonormal AO basis and saved both in the array p and    |
C |  in the RWF file irwp.                                     |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       iopcl  - describes type of HF wavefunction
C       na     - number of alpha electrons
C       nb     - number of beta electrons
C       nbasis - number of basis functions
C       norbs  - number of orbitals
C       irwx   - RWF file with transformation matrix (=S^(-1/2))
C       irws   - RWF file with overlap matrix
C       irwca  - Gaussian RWF file for MO coefficients (alpha)
C       irwcb  - Gaussian RWF file for MO coefficients (beta)
C       irwvec - RWF file with orbital coefficients
C                (orthonormal AO basis) [ updated here ]
C       irwp   - RWF file with density matrix (orthonormal AO basis)
C                [ updated here ]
C       p      - HF density matrix in orthonormal AO basis [ out ]

      real*8      v(*), p(*)
      integer     iout, iprint, lenv, iphf, iopcl
      integer     na, nb, nbasis, norbs, igues1
      integer     irwx, irws, irwca, irwcb
      integer     irwvec, irwp

C     other variables

      integer     nel, ntt, ntto, ntt2o, nbsq, nosq, nri, itype
      integer     i, j, ind1, ind2
      integer     mdv, jend, jscrx
      integer     jx, js, jt, jorb, jp
      integer     jzaa, jzab, jzba, jzbb, jscr, jscr1, lscr1
      integer     jtst1
      logical     icmplx, ispuhf, ispghf
      logical     irhf, iuhf, ighf
      complex*16  fac

C     constants

      real*8      pi
      complex*16  zero, one, im


      pi = 4.0d0 * atan (1.0d0)

      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)
      im   = cmplx (0.0d0, 1.0d0)


C     Construct useful variables.

      nel = na + nb

      ntt = nbasis * (nbasis + 1) / 2

      ntto = norbs * (norbs + 1) / 2
      ntt2o = norbs * (2*norbs + 1)

      nbsq = nbasis * nbasis
      nosq = norbs * norbs


C     Determine nri.

      nri = 2

      if ( iopcl .le. 1 ) nri = 1


C     Check iopcl and decipher iopcl.

      if ( iopcl .gt. 3 .and. iopcl .ne. 6 ) then
        call gauerr ('Incorrect iopcl in iguess.')
      endif

      irhf = iopcl .eq. 0 .or. iopcl .eq. 2
      iuhf = iopcl .eq. 1 .or. iopcl .eq. 3
      ighf = iopcl .eq. 6


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
        call gauerr ('Unrecognized iphf in iguess.')
      endif


C     Make sure that iopcl and itype are consistent.
C     Kill the calculation if there are inconsistencies.

      if ( itype .eq. 1 .and. .not. irhf ) then
        call gauerr ('Inconsistent iphf and initial guess orbitals.')
      endif

      if ( ( itype .eq. 2 .or. ispuhf ) .and. .not. iuhf ) then
        call gauerr ('Inconsistent iphf and initial guess orbitals.')
      endif

      if ( ( itype .eq. 3 .and. .not. ispuhf ) .and. .not. ighf ) then
        call gauerr ('Inconsistent iphf and initial guess orbitals.')
      endif


C     Memory allocation.

C       Allocate space for:
C         - transformation matrix (jx)
C         - overlap matrix (js)
C         - scratch array for matrix multiplications (jt)
C         - orbital coefficients (jorb)
C         - density matrix, square format (jp)

      jx   = 1
      js   = jx + nri*nbasis*norbs
      jt   = js + nri*nbsq
      jend = jt + nri*nbasis*norbs

      if ( itype .eq. 1 ) then
        jorb = jend
        jp   = jorb + 2*nbsq
        jend = jp   + 2*nosq
      
      elseif ( itype .eq. 2 .or. ispuhf ) then
        jorb = jend
        jp   = jorb + 4*nbsq
        jend = jp   + 2*nosq

      elseif ( itype .eq. 3 ) then
        jorb = jend
        jp   = jorb + 8*nbsq
        jend = jp   + 8*nosq
      endif

C       Allocate space for:
C         - real to complex transformation of orbitals (jscrx)

      if ( nri .eq. 1 ) then
        jscrx = jend
        jend  = jscrx + nosq
      endif

C       Allocate space for:
C         - scratch array for matrix reorganization if itype = 3 (jscr)
C         - scratch spin blocks (jzaa, jzab, jzba, jzbb)
C         - scratch space for calls to dmblck (jscr1)

      if ( itype .eq. 3 ) then
        lscr1 = 16*nbsq

        jscr  = jend
        jzaa  = jscr  + 8*nbsq
        jzab  = jzaa  + 2*nbsq
        jzba  = jzab  + 2*nbsq
        jzbb  = jzba  + 2*nbsq
        jscr1 = jzbb  + 2*nbsq
        jend  = jscr1 + lscr1
      endif

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'iguess')



C     %%%%%%%%%%%%%%%%%%%
C     %  Load matrices  %
C     %%%%%%%%%%%%%%%%%%%


C     Load transformation matrix.

      call fileio (2, -irwx, nri*nbasis*norbs, v(jx), 0)
      write (iout, *) ' '
      call ltoutx (iout, 3, 0, 1, 0, 'Transform Matrix', nbasis, nbasis,
     $       v(jx), v(jx), 0)


C     Load overlap matrix.

      call fileio (2, -irws, nri*ntt, v(js), 0)
C        write (iout, *) ' '
C        call ltoutx (iout, 0, 0, 1, 0, 'Overlap Matrix', nbasis, nbasis,
C     $       v(js), v(js), 0)

C     Expand the overlap to square matrix.

      if ( nri .eq. 1 ) then
        call square (v(js), v(js), nbasis, nbasis, 0)

      elseif ( nri .eq. 2 ) then
        call zsquare (v(js), v(js), nbasis, nbasis, 0)
      endif


C     Clear arrays of orbital coefficients.

      if ( itype .eq. 1 ) then
        call aclear (2*nbsq, v(jorb))

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call aclear (4*nbsq, v(jorb))

      elseif ( itype .eq. 3 ) then
        call aclear (8*nbsq, v(jorb))
      endif


C     Read gaussian orbitals from file.

      if ( irhf ) then
        call fileio (2, -irwca, nri*nbasis*norbs, v(jorb), 0)

      elseif ( iuhf ) then
        call fileio (2, -irwca, nri*nbasis*norbs, v(jorb), 0)
        call fileio (2, -irwcb, nri*nbasis*norbs, v(jorb+2*nbsq), 0)

      elseif ( ighf ) then
        call fileio (2, -irwca, 8*nbasis*norbs, v(jorb), 0)
      endif


C     debugging...
C       print molecular orbitals

      if ( itype .eq. 1 .and. nri .eq. 1 ) then
        write (iout, *) ' '
        call ltoutx (iout, 3, 0, 1, 0, 'MO coefficients', nbasis, norbs,
     $       v(jorb), v(jorb), 0)

      elseif ( ( itype .eq. 2 .or. ispuhf ) .and. nri .eq. 1 ) then
        write (iout, *) ' '
        call ltoutx (iout, 3, 1, 1, 0, 'MO coefficients', nbasis, norbs,
     $       v(jorb), v(jorb+2*nbsq), 0)

      elseif ( itype .eq. 3 ) then
        write (iout, *) ' '
        call ltoutx (iout, 3, 2, 1, 0, 'MO coefficients', 2*nbasis,
     $       2*norbs, v(jorb), v(jorb), 0)
      endif


C     Reorganize the matrix for GHF-type guesses (see fixmat for more
C     details).

      if ( ighf ) then

C         .. DO NOT reorganize columns

c$$$        do 10 i = 1, 2*nbasis
c$$$          if ( mod (i,2) .eq. 1 ) then
c$$$            j = (i-1)/2 + 1
c$$$          elseif ( mod (i,2) .eq. 0 ) then
c$$$            j = i/2 + nbasis
c$$$          endif
c$$$
c$$$          ind1 = (i - 1) * 4*nbasis + jorb
c$$$          ind2 = (j - 1) * 4*nbasis + jscr
c$$$
c$$$          call amove (4*nbasis, v(ind1), v(ind2))
c$$$ 10     continue

        call amove (8*nbsq, v(jorb), v(jscr))

C         .. now reorganize rows

        call trspnc (.false., 2*nbasis, 2*nbasis, v(jscr))

        do 11 i = 1, 2*nbasis
          if ( mod (i,2) .eq. 1 ) then
            j = (i-1)/2 + 1
          elseif ( mod (i,2) .eq. 0 ) then
            j = i/2 + nbasis
          endif

          ind1 = (i - 1) * 4*nbasis + jscr
          ind2 = (j - 1) * 4*nbasis + jorb

          call amove (4*nbasis, v(ind1), v(ind2))
 11     continue

        call trspnc (.false., 2*nbasis, 2*nbasis, v(jorb))
      endif


C     debugging...
C       print molecular orbitals

      if ( itype .eq. 3 .and. .not. ispuhf ) then
        write (iout, *) ' '
        call ltoutx (iout, 3, 2, 1, 0, 'MO coefficients', 2*nbasis,
     $       2*norbs, v(jorb), v(jorb), 0)
      endif



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Transform orbitals to orthonormal AO basis  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Transform orbitals to orthonormal AO basis:
C       X! . S . orb

C       This transformation can be understood as follows. In the AO basis,
C       we have the following orthonormality condition:
C         orb! . S . orb = I

C       Evidently, we can split it as
C         [ orb! . S^(1/2) ] . [ S^(1/2) . orb ] = I

C       Thus our orthonormal orbitals should be built from
C         S^(1/2) . orb

C       We compute S^(1/2) quickly by multiplying S and S^(-1/2), both of
C       which are available already: X . S.


      if ( irhf ) then

        if ( nri .eq. 1 ) then
          call dgemm ('n', 'n', nbasis, norbs, nbasis, 1.0d0, v(js),
     $         nbasis, v(jorb), nbasis, 0.0d0, v(jt), nbasis)
          call dgemm ('t', 'n', norbs, norbs, nbasis, 1.0d0, v(jx),
     $         nbasis, v(jt), nbasis, 0.0d0, v(jorb), norbs)

        elseif ( nri .eq. 2 ) then
          call zgemm ('n', 'n', nbasis, norbs, nbasis, one, v(js),
     $         nbasis, v(jorb), nbasis, zero, v(jt), nbasis)
          call zgemm ('c', 'n', norbs, norbs, nbasis, one, v(jx),
     $         nbasis, v(jt), nbasis, zero, v(jorb), norbs)
        endif

      elseif ( iuhf ) then

        if ( nri .eq. 1 ) then
          call dgemm ('n', 'n', nbasis, norbs, nbasis, 1.0d0, v(js),
     $         nbasis, v(jorb), nbasis, 0.0d0, v(jt), nbasis)
          call dgemm ('t', 'n', norbs, norbs, nbasis, 1.0d0, v(jx),
     $         nbasis, v(jt), nbasis, 0.0d0, v(jorb), norbs)

          call dgemm ('n', 'n', nbasis, norbs, nbasis, 1.0d0, v(js),
     $         nbasis, v(jorb+2*nbsq), nbasis, 0.0d0, v(jt), nbasis)
          call dgemm ('t', 'n', norbs, norbs, nbasis, 1.0d0, v(jx),
     $         nbasis, v(jt), nbasis, 0.0d0, v(jorb+2*nbsq), norbs)

        elseif ( nri .eq. 2 ) then
          call zgemm ('n', 'n', nbasis, norbs, nbasis, one, v(js),
     $         nbasis, v(jorb), nbasis, zero, v(jt), nbasis)
          call zgemm ('c', 'n', norbs, norbs, nbasis, one, v(jx),
     $         nbasis, v(jt), nbasis, zero, v(jorb), norbs)

          call zgemm ('n', 'n', nbasis, norbs, nbasis, one, v(js),
     $         nbasis, v(jorb+2*nbsq), nbasis, zero, v(jt), nbasis)
          call zgemm ('c', 'n', norbs, norbs, nbasis, one, v(jx),
     $         nbasis, v(jt), nbasis, zero, v(jorb+2*nbsq), norbs)
        endif

      elseif ( ighf ) then

C         .. split into spin blocks ..

        call dmblck (iout, iprint, v(jscr1), lscr1, nbasis, v(jorb),
     $       v(jzaa), v(jzab), v(jzba), v(jzbb), 2, -1)

C         .. transform [aa] block

        call zgemm ('n', 'n', nbasis, norbs, nbasis, one, v(js), nbasis,
     $       v(jzaa), nbasis, zero, v(jt), nbasis)
        call zgemm ('c', 'n', norbs, norbs, nbasis, one, v(jx), nbasis,
     $       v(jt), nbasis, zero, v(jzaa), norbs)

C         .. transform [ab] block

        call zgemm ('n', 'n', nbasis, norbs, nbasis, one, v(js), nbasis,
     $       v(jzab), nbasis, zero, v(jt), nbasis)
        call zgemm ('c', 'n', norbs, norbs, nbasis, one, v(jx), nbasis,
     $       v(jt), nbasis, zero, v(jzab), norbs)

C         .. transform [ba] block

        call zgemm ('n', 'n', nbasis, norbs, nbasis, one, v(js), nbasis,
     $       v(jzba), nbasis, zero, v(jt), nbasis)
        call zgemm ('c', 'n', norbs, norbs, nbasis, one, v(jx), nbasis,
     $       v(jt), nbasis, zero, v(jzba), norbs)

C         .. transform [bb] block

        call zgemm ('n', 'n', nbasis, norbs, nbasis, one, v(js), nbasis,
     $       v(jzbb), nbasis, zero, v(jt), nbasis)
        call zgemm ('c', 'n', norbs, norbs, nbasis, one, v(jx), nbasis,
     $       v(jt), nbasis, zero, v(jzbb), norbs)

C         .. combine spin blocks ..

        call dmblck (iout, iprint, v(jscr1), lscr1, norbs, v(jorb),
     $       v(jzaa), v(jzab), v(jzba), v(jzbb), 1, -1)
      endif


C     Make the orbitals complex if they are still real.

      if ( irhf .and. nri .eq. 1 ) then
        call amove (nosq, v(jorb), v(jscrx))
        call zrecmp (3, v(jorb), v(jscrx), nosq)

      elseif ( iuhf .and. nri .eq. 1 ) then
        call amove (nosq, v(jorb), v(jscrx))
        call zrecmp (3, v(jorb), v(jscrx), nosq)

        call amove (nosq, v(jorb+2*nbsq), v(jscrx))
        call zrecmp (3, v(jorb+2*nbsq), v(jscrx), nosq)
      endif


C     debugging...
C       print molecular orbitals
      write (iout, *) 'After transformation '

      if ( itype .eq. 1 ) then
        write (iout, *) ' '
        call ltoutx (iout, 1, 2, 1, 0, 'MO coefficients', norbs, norbs,
     $       v(jorb), v(jorb), 0)
        write (iout, *) ' '
        call ltoutx (iout, 0, 2, 1, 0, 'Overlap Matrix', nbasis, nbasis,
     $       v(js), v(js), 0)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        write (iout, *) ' '
        call ltoutx (iout, 1, 3, 1, 0, 'MO coefficients', norbs, norbs,
     $       v(jorb), v(jorb+2*nbsq), 0)
        write (iout, *) ' '
        call ltoutx (iout, 1, 0, 1, 0, 'Overlap Matrix', nbasis, nbasis,
     $       v(js), v(js), 0)
        write (iout, *) ' '
        call ltoutx (iout, 1, 0, 1, 0,'Transform Matrix', nbasis,nbasis,
     $       v(jx), v(jx), 0)
        call dgemm ('t', 'n', nbasis, nbasis, nbasis, 1.0d0, v(jx),
     $         nbasis, v(js), nbasis, 0.0d0, v(jt), nbasis)
        write (iout, *) ' '
        call ltoutx (iout, 1, 0, 1, 0,'S^(1/2) Matrix', nbasis,nbasis,
     $       v(jt), v(jt), 0)

      elseif ( itype .eq. 3 ) then
        write (iout, *) ' '
        call ltoutx (iout, 1, 2, 1, 0, 'MO coefficients', 2*norbs,
     $       2*norbs, v(jorb), v(jorb), 0)
      endif


C     debugging...
C       Is C! . C = 1 ?.

c$$$      if ( itype .eq. 1 ) then
c$$$        jtst1 = jend
c$$$        jend  = jtst1 + 2*nosq
c$$$
c$$$        call tstcor (jend-1, lenv, 'iguess')
c$$$
c$$$        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jorb), norbs,
c$$$     $       v(jorb), norbs, zero, v(jtst1), norbs)
c$$$
c$$$        write (iout, *) ' '
c$$$        call ltoutx (iout, 1, 2, 1, 0, 'C! . C', norbs, norbs, v(jtst1),
c$$$     $       v(jtst1), 0)
c$$$
c$$$      elseif ( itype .eq. 2 .or. ispuhf ) then
c$$$        jtst1 = jend
c$$$        jend  = jtst1 + 4*nosq
c$$$
c$$$        call tstcor (jend-1, lenv, 'iguess')
c$$$
c$$$        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jorb), norbs,
c$$$     $       v(jorb), norbs, zero, v(jtst1), norbs)
c$$$        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jorb+2*nbsq),
c$$$     $       norbs, v(jorb+2*nbsq), norbs, zero, v(jtst1+2*nosq), norbs)
c$$$
c$$$        write (iout, *) ' '
c$$$        call ltoutx (iout, 1, 3, 1, 0, 'C! . C', norbs, norbs, v(jtst1),
c$$$     $       v(jtst1+2*nosq), 0)
c$$$
c$$$      elseif ( itype .eq. 3 ) then
c$$$        jtst1 = jend
c$$$        jend  = jtst1 + 8*nosq
c$$$
c$$$        call tstcor (jend-1, lenv, 'iguess')
c$$$
c$$$        call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jorb),
c$$$     $       2*norbs, v(jorb), 2*norbs, zero, v(jtst1), 2*norbs)
c$$$
c$$$        write (iout, *) ' '
c$$$        call ltoutx (iout, 1, 2, 1, 0, 'C! . C', 2*norbs, 2*norbs,
c$$$     $       v(jtst1), v(jtst1), 0)
c$$$      endif


C     Apply a complex factor to some of the elements of the HOMO.
C       ( for UHF type wavefunctions, the factor is applied in
C         both the alpha HOMO and the beta HOMO )

      if ( igues1 .eq. 2 ) then

        fac = exp (im * cmplx (pi) / cmplx (3.0d0)) ! used to be 8.0d0

        if ( itype .eq. 1 ) then
          ind1 = (na - 1) * 2*norbs + jorb

          call zscal (na, fac, v(ind1), 1)

        elseif ( itype .eq. 2 .or. ispuhf ) then
          ind1 = (na - 1) * 2*norbs + jorb
          ind2 = (nb - 1) * 2*norbs + jorb + 2*nbsq

C         does scaling just beta work better?

          call zscal (na, fac, v(ind1), 1)
c$$$          call zscal (nb, fac, v(ind2), 1)

        elseif ( itype .eq. 3 ) then
          ind1 = (nel - 1) * 4*norbs + jorb

          call zscal (nel, fac, v(ind1), 1)
        endif
      endif


C     Save orbitals to RWF file.

      if ( itype .eq. 1 ) then
        call fileio (1, -irwvec, 2*nosq, v(jorb), 0)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call fileio (1, -irwvec, 2*nosq, v(jorb), 0)
        call fileio (1, -irwvec, 2*nosq, v(jorb+2*nbsq), 2*nosq)

      elseif ( itype .eq. 3 ) then
        call fileio (1, -irwvec, 8*nosq, v(jorb), 0)
      endif



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Form density matrix in orthonormal AO basis  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Form density matrix in orthonormal basis.

C       If the matrix of orbital coefficients is C, then the density
C       matrix is formed by
C         P = C . N . C!,

C       where N is a diagonal matrix with canonical orbital occupations.

C       Here, we use the LaPack routine zherk, which forms the rank-k
C       product of C with C!. See the zherk documentation for more
C       details.

      if ( itype .eq. 1 ) then
        call zherk ('u', 'n', norbs, na, 1.0d0, v(jorb), norbs, 0.0d0,
     $       v(jp), norbs)

        call zlinear (v(jp), p, norbs, norbs)
        call vconjg (ntto, p)

        call fileio (1, -irwp, 2*ntto, p, 0)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call zherk ('u', 'n', norbs, na, 1.0d0, v(jorb), norbs, 0.0d0,
     $       v(jp), norbs)

        call zlinear (v(jp), p, norbs, norbs)
        call vconjg (ntto, p)

        call fileio (1, -irwp, 2*ntto, p, 0)

        if ( nb .ge. 1 ) then
        call zherk ('u', 'n', norbs, nb, 1.0d0, v(jorb+2*nbsq), norbs,
     $       0.0d0, v(jp), norbs)
        else
          call aclear (2*nosq, v(jp))
        endif

        call zlinear (v(jp), p(1+2*ntto), norbs, norbs)
        call vconjg (ntto, p(1+2*ntto))

        call fileio (1, -irwp, 2*ntto, p(1+2*ntto), 2*ntto)

      elseif ( itype .eq. 3 ) then
        call zherk ('u', 'n', 2*norbs, nel, 1.0d0, v(jorb), 2*norbs,
     $       0.0d0, v(jp), 2*norbs)

        call zlinear (v(jp), p, 2*norbs, 2*norbs)
        call vconjg (ntt2o, p)

        call fileio (1, -irwp, 2*ntt2o, p, 0)
      endif


C     debugging...
C       print density matrix

      if ( itype .eq. 1 ) then
        write (iout, *) ' '
        call ltoutx (iout, 0, 2, 1, 0, 'density matrix', norbs, norbs,
     $       p, p, 0)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        write (iout, *) ' '
        call ltoutx (iout, 0, 3, 1, 0, 'density matrix', norbs, norbs,
     $       p, p(1+2*ntto), 0)

      elseif ( itype .eq. 3 ) then
        write (iout, *) ' '
        call ltoutx (iout, 0, 2, 1, 0, 'density matrix', 2*norbs,
     $       2*norbs, p, p, 0)
      endif


      return
      end


