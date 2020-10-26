

      subroutine gaufmt (iout, iprint, v, lenv, iphf, nbasis, norbs,
     $     irwx, irwvec, irwp, irwca, irwcb, irwpa, irwpb)

      implicit none


C +------------------------------------------------------------+
C |  gaufmt   --   CAJH, 07.2011                               |
C |                                                            |
C |                                                            |
C |  This routine retrieves the converged density matrix and   |
C |  MO coefficients from a PHF calculation (which are stored  |
C |  stored in an orthonormal AO basis). It then transforms    |
C |  them into the regular AO basis and stores them in the     |
C |  RWF files that Gaussian uses for MO coefficients          |
C |  (irwca, irwcb) and density matrices (irwpa, irwpb).       |
C |                                                            |
C |  The transformation of GHF-type matrices deserves a few    |
C |  words. Throughout the PHF code, we handle GHF matrices    |
C |  in terms of spin blocks. That is, a general matrix S      |
C |  will be stored as                                         |
C |                                                            |
C |    S  =  ( S_aa  S_ab )                                    |
C |          ( S_ba  S_bb )                                    |
C |                                                            |
C |  On the other hand, Gaussian stores GHF matrices in a      |
C |  different order. The basis functions are ordered          |
C |  according to                                              |
C |                                                            |
C |    PHF:       1a, 2a, ..., 1b, 2b, ...                     |
C |    Gaussian:  1a, 1b, 2a, 2b, ...                          |
C |                                                            |
C |  Thus, for GHF matrices, we need to reorganize rows and    |
C |  columns to bring the matrices back to Gaussian standard   |
C |  format.                                                   |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       nbasis - number of basis functions
C       norbs  - number of orbitals
C       irwx   - RWF file with transformation matrix (=S^(-1/2))
C       irwvec - RWF file with orbital coefficients
C                (orthonormal AO basis)
C       irwp   - RWF file with density matrix (orthonormal AO basis)
C       irwca  - Gaussian RWF file for MO coefficients (alpha)
C       irwcb  - Gaussian RWF file for MO coefficients (beta)
C       irwpa  - Gaussian RWF file for density matrix  (alpha)
C       irwpb  - Gaussian RWF file for density matrix  (beta)

C         ** irwca, irwcb, irwpa, irwpb are all updated here

      real*8      v(*)
      integer     iout, iprint, lenv, iphf, nbasis, norbs
      integer     irwx, irwvec, irwp, irwca, irwcb, irwpa, irwpb

C     other variables

      real*8      junk
      integer     ntt, ntt2, ntto, ntt2o, nosq, nbsq
      integer     itype, iopcl, nri
      integer     i, j, ind1, ind2
      integer     mdv, jend
      integer     jx, jxt, jp, jp2, jorb, jscr
      integer     jz1, jz2, jz3, jscrx
      integer     jzaa, jzab, jzba, jzbb
      integer     lscr1, lscr2, jscr1, jscr2
      logical     icmplx, ispuhf, ispghf

C     constants

      complex*16  zero, one



      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)

      junk = 0.0d0


C     Construct useful variables.

      ntt = nbasis * (nbasis + 1) / 2
      ntt2 = nbasis * (2*nbasis + 1)

      ntto = norbs * (norbs + 1) / 2
      ntt2o = norbs * (2*norbs + 1)

      nosq = norbs * norbs
      nbsq = nbasis * nbasis


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
        call gauerr ('Unrecognized iphf in gaufmt.')
      endif


C     Read iopcl from ILSW file.

C       iopcl = 0,  real RHF
C             = 1,  real UHF
C             = 2,  complex RHF
C             = 3,  complex UHF
C             = 6,  complex GHF (there is no real GHF)

      call ilsw (2, 1, iopcl)

      if ( iopcl .gt. 3 .and. iopcl .ne. 6 ) then
        call gauerr ('Incorrect iopcl in gaufmt.')
      endif

      nri = 2

      if ( iopcl .le. 1 ) nri = 1


C     Check that the Gaussian RWF files have the correct dimension.

      if ( itype .eq. 1 ) then
        call conddf (irwca, 2*nbasis*norbs)
        call conddf (irwpa, 2*ntt)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call conddf (irwca, 2*nbasis*norbs)
        call conddf (irwpa, 2*ntt)

        call conddf (irwcb, 2*nbasis*norbs)
        call conddf (irwpb, 2*ntt)

      elseif ( itype .eq. 3 ) then
        call conddf (irwca, 8*nbasis*norbs)
        call conddf (irwpa, 2*ntt2)
      endif


C     Memory allocation.

C       Allocate space for:
C         - transformation matrix (jx)
C         - matrix of orbital coefficients (jorb)
C         - scratch space for MO transformation (jxt)
C         - scratch space for matrix reorganization (jscr) (itype = 3)
C         - density matrix (jp)
C         - scratch arrays for lt to square transformation (jz1, jz2, jz3)
C         - density matrix in lower triangular form (jp2)

      if ( itype .eq. 1 ) then
        jx   = 1
        jorb = jx   + 2*nbasis*norbs
        jxt  = jorb + 2*nbasis*norbs
        jp   = jxt  + 2*nbasis*norbs
        jz1  = jp   + 2*nbsq
        jz2  = jz1  + 2*ntto
        jz3  = jz2  + nosq
        jp2  = jz3  + nosq
        jend = jp2  + 2*ntt

      elseif ( itype .eq. 2 .or. ispuhf ) then
        jx   = 1
        jorb = jx   + 2*nbasis*norbs
        jxt  = jorb + 4*nbsq
        jp   = jxt  + 2*nbasis*norbs
        jz1  = jp   + 4*nbsq
        jz2  = jz1  + 2*ntto
        jz3  = jz2  + nosq
        jp2  = jz3  + nosq
        jend = jp2  + 4*ntt

      elseif ( itype .eq. 3 ) then
        jx   = 1
        jorb = jx   + 2*nbasis*norbs
        jxt  = jorb + 8*nbsq
        jscr = jxt  + 2*nbasis*norbs
        jp   = jscr + 8*nbsq
        jz1  = jp   + 8*nbsq
        jz2  = jz1  + 2*ntt2o
        jz3  = jz2  + 4*nosq
        jp2  = jz3  + 4*nosq
        jend = jp2  + 2*ntt2

      endif

C       Allocate scratch space if transformation matrix is real.

      if ( nri .eq. 1 ) then
        jscrx = jend
        jend  = jscrx + 2*nbasis*norbs
      endif

C       Allocate extra scratch arrays for itype = 3.
C         - scratch spin blocks (jzaa, jzab, jzba, jzbb)
C         - scratch space for calls to dmblck (jscr1)

      if ( itype .eq. 3 .and. .not. ispuhf ) then
        lscr1 = 16*nosq

        jzaa  = jend
        jzab  = jzaa  + 2*nbsq
        jzba  = jzab  + 2*nbsq
        jzbb  = jzba  + 2*nbsq
        jscr1 = jzbb  + 2*nbsq
        jend  = jscr1 + lscr1
      endif

C       Allocate space for call to aobstf.

      lscr2 = 4*nbasis*norbs

      jscr2 = jend
      jend  = jscr2 + lscr2


      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'gaufmt')



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Transform MO coefficients  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Load transformation matrix.

      call fileio (2, -irwx, nri*nbasis*norbs, v(jx), 0)

      if ( nri .eq. 1 ) then
        call zrecmp (3, v(jscrx), v(jx), nbasis*norbs)
        call amove (2*nbasis*norbs, v(jscrx), v(jx))
      endif


C     Load matrix of orbital coefficients in orthonormal AO basis.

      if ( itype .eq. 1 ) then
        call fileio (2, -irwvec, 2*nosq, v(jorb), 0)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call fileio (2, -irwvec, 2*nosq, v(jorb), 0)
        call fileio (2, -irwvec, 2*nosq, v(jorb+2*nbsq), 2*nosq)

      elseif ( itype .eq. 3 ) then
        call fileio (2, -irwvec, 8*nosq, v(jorb), 0)

C         .. split into spin blocks ..

        call dmblck (iout, iprint, v(jscr1), lscr1, norbs, v(jorb),
     $       v(jzaa), v(jzab), v(jzba), v(jzbb), 2, -1)
      endif


C     Transform orbital coefficients into AO basis:
C       X . C,

C     where X is the transformation matrix and C is the set of orbital
C     coefficients.

      if ( itype .eq. 1 ) then
        call zgemm ('n', 'n', nbasis, norbs, norbs, one, v(jx), nbasis,
     $       v(jorb), norbs, zero, v(jxt), nbasis)

        call amove (2*nbasis*norbs, v(jxt), v(jorb))

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call zgemm ('n', 'n', nbasis, norbs, norbs, one, v(jx), nbasis,
     $       v(jorb), norbs, zero, v(jxt), nbasis)

        call amove (2*nbasis*norbs, v(jxt), v(jorb))

        call zgemm ('n', 'n', nbasis, norbs, norbs, one, v(jx), nbasis,
     $       v(jorb+2*nbsq), norbs, zero, v(jxt), nbasis)

        call amove (2*nbasis*norbs, v(jxt), v(jorb+2*nbsq))

      elseif ( itype .eq. 3 ) then

C         .. [aa] block

        call zgemm ('n', 'n', nbasis, norbs, norbs, one, v(jx), nbasis,
     $       v(jzaa), norbs, zero, v(jxt), nbasis)

        call aclear (2*nbsq, v(jzaa))
        call amove (2*nbasis*norbs, v(jxt), v(jzaa))

C         .. [ab] block

        call zgemm ('n', 'n', nbasis, norbs, norbs, one, v(jx), nbasis,
     $       v(jzab), norbs, zero, v(jxt), nbasis)

        call aclear (2*nbsq, v(jzab))
        call amove (2*nbasis*norbs, v(jxt), v(jzab))

C         .. [ba] block

        call zgemm ('n', 'n', nbasis, norbs, norbs, one, v(jx), nbasis,
     $       v(jzba), norbs, zero, v(jxt), nbasis)

        call aclear (2*nbsq, v(jzba))
        call amove (2*nbasis*norbs, v(jxt), v(jzba))

C         .. [bb] block

        call zgemm ('n', 'n', nbasis, norbs, norbs, one, v(jx), nbasis,
     $       v(jzbb), norbs, zero, v(jxt), nbasis)

        call aclear (2*nbsq, v(jzbb))
        call amove (2*nbasis*norbs, v(jxt), v(jzbb))

C         .. combine spin blocks ..

        call dmblck (iout, iprint, v(jscr1), lscr1, nbasis, v(jorb),
     $       v(jzaa), v(jzab), v(jzba), v(jzbb), 1, -1)
      endif

C     debugging...
C       print molecular orbitals

C$$$      if ( itype .eq. 1 ) then
C$$$        write (iout, *) ' '
C$$$        call ltoutx (iout, 3, 2, 1, 0, 'MO coefficients', nbasis, norbs,
C$$$     $       v(jorb), v(jorb), 0)
C$$$
C$$$      elseif ( itype .eq. 2 .or. ispuhf ) then
C$$$        write (iout, *) ' '
C$$$        call ltoutx (iout, 3, 3, 1, 0, 'MO coefficients', nbasis, norbs,
C$$$     $       v(jorb), v(jorb+2*nosq), 0)
C$$$
C$$$      elseif ( itype .eq. 3 ) then
C$$$        write (iout, *) ' '
C$$$        call ltoutx (iout, 3, 2, 1, 0, 'MO coefficients', 2*nbasis,
C$$$     $       2*norbs, v(jorb), v(jorb), 0)
C$$$      endif

C     Reorder rows and columns for GHF-type MO coefficients.

C       Gaussian stores matrices by ordering basis functions as
C         1a, 1b, 2a, 2b, ...

C       On the other hand, we have them organized as
C         1a, 2a, ..., 1b, 2b, ...

      if ( itype .eq. 3 .and. .not. ispuhf ) then

C         .. DO NOT reorganize columns

c$$$        do 10 i = 1, 2*nbasis
c$$$          if ( i .le. nbasis ) then
c$$$            j = 2*i - 1
c$$$          elseif ( i .gt. nbasis ) then
c$$$            j = 2*i - 2*nbasis
c$$$          endif
c$$$
c$$$          ind1 = (i - 1) * 4*nbasis + jorb
c$$$          ind2 = (j - 1) * 4*nbasis + jscr
c$$$ 10     continue

        call amove (8*nbsq, v(jorb), v(jscr))

C         .. now reorganize rows

        call trspnc (.false., 2*nbasis, 2*nbasis, v(jscr))

        do 11 i = 1, 2*nbasis
          if ( i .le. nbasis ) then
            j = 2*i - 1
          elseif ( i .gt. nbasis ) then
            j = 2*i - 2*nbasis
          endif

          ind1 = (i - 1) * 4*nbasis + jscr
          ind2 = (j - 1) * 4*nbasis + jorb

          call amove (4*nbasis, v(ind1), v(ind2))
 11     continue

        call trspnc (.false., 2*nbasis, 2*nbasis, v(jorb))
      endif


C     debugging...
C       print molecular orbitals

C$$$      if ( itype .eq. 1 ) then
C$$$        write (iout, *) ' '
C$$$        call ltoutx (iout, 3, 2, 1, 0, 'MO coefficients', nbasis, norbs,
C$$$     $       v(jorb), v(jorb), 0)
C$$$
C$$$      elseif ( itype .eq. 2 .or. ispuhf ) then
C$$$        write (iout, *) ' '
C$$$        call ltoutx (iout, 3, 3, 1, 0, 'MO coefficients', nbasis, norbs,
C$$$     $       v(jorb), v(jorb+2*nosq), 0)
C$$$
C$$$      elseif ( itype .eq. 3 ) then
C$$$        write (iout, *) ' '
C$$$        call ltoutx (iout, 3, 2, 1, 0, 'MO coefficients', 2*nbasis,
C$$$     $       2*norbs, v(jorb), v(jorb), 0)
C$$$      endif



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Transform density matrix  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Load density matrix in orthonormal AO basis.

      if ( itype .eq. 1 ) then
        call fileio (2, -irwp, 2*ntto, v(jp), 0)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call fileio (2, -irwp, 2*ntto, v(jp), 0)
        call fileio (2, -irwp, 2*ntto, v(jp+2*nbsq), 2*ntto)

      elseif ( itype .eq. 3 ) then
        call fileio (2, -irwp, 2*ntt2o, v(jp), 0)
      endif


C     Transform density matrix to square format.

      if ( itype .eq. 1 ) then
        call zrecmp (1, v(jp), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jp), v(jz2), nosq)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call zrecmp (1, v(jp), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jp), v(jz2), nosq)

        call zrecmp (1, v(jp+2*nbsq), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jp+2*nbsq), v(jz2), nosq)

      elseif ( itype .eq. 3 ) then
        call zrecmp (1, v(jp), v(jz1), ntt2o)
        call square (v(jz1), v(jz2), 2*norbs, 2*norbs, 0)
        call square (v(jz1+ntt2o), v(jz3), 2*norbs, 2*norbs, -1)
        call zrecmp (2, v(jp), v(jz2), 4*nosq)
      endif


C     Transform density matrix to AO basis.
C     We use the aobstf subroutine to do such transformation.

      if ( itype .eq. 1 ) then
        call aobstf (iout, iprint, v(jscr2), lscr2, norbs, nbasis, irwx,
     $       junk, v(jp), 0.0d0, 0.0d0, 0.0d0, .true., 3, 1)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call aobstf (iout, iprint, v(jscr2), lscr2, norbs, nbasis, irwx,
     $       junk, v(jp), 0.0d0, 0.0d0, v(jp+2*nbsq), .true., 2, 1)

      elseif ( itype .eq. 3 ) then

C         .. split into spin blocks ..

        call dmblck (iout, iprint, v(jscr1), lscr1, norbs, v(jp),
     $       v(jzaa), v(jzab), v(jzba), v(jzbb), 2, -1)

C         .. perform transformation ..

        call aobstf (iout, iprint, v(jscr2), lscr2, norbs, nbasis, irwx,
     $       junk, v(jzaa), v(jzab), v(jzba), v(jzbb), .true., 1, 1)

C         .. combine spin blocks ..

        call dmblck (iout, iprint, v(jend), mdv, nbasis, v(jp),
     $       v(jzaa), v(jzab), v(jzba), v(jzbb), 1, -1)
      endif

C     Reorder rows and columns for GHF-type density matrix.

      if ( itype .eq. 3 .and. .not. ispuhf ) then

C       .. reorganize columns

        do 20 i = 1, 2*nbasis
          if ( i .le. nbasis ) then
            j = 2*i - 1
          elseif ( i .gt. nbasis ) then
            j = 2*i - 2*nbasis
          endif

          ind1 = (i - 1) * 4*nbasis + jp
          ind2 = (j - 1) * 4*nbasis + jscr

          call amove (4*nbasis, v(ind1), v(ind2))
 20     continue

C       .. now reorganize rows

        call trspnc (.false., 2*nbasis, 2*nbasis, v(jscr))

        do 21 i = 1, 2*nbasis
          if ( i .le. nbasis ) then
            j = 2*i - 1
          elseif ( i .gt. nbasis ) then
            j = 2*i - 2*nbasis
          endif

          ind1 = (i - 1) * 4*nbasis + jscr
          ind2 = (j - 1) * 4*nbasis + jp

          call amove (4*nbasis, v(ind1), v(ind2))
 21     continue

        call trspnc (.false., 2*nbasis, 2*nbasis, v(jp))
      endif


C     Bring density matrix to lower triangular form.

      if ( itype .eq. 1 ) then
        call zlinear (v(jp), v(jp2), nbasis, nbasis)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call zlinear (v(jp), v(jp2), nbasis, nbasis)
        call zlinear (v(jp+2*nbsq), v(jp2+2*ntt), nbasis, nbasis)

      elseif ( itype .eq. 3 ) then
        call zlinear (v(jp), v(jp2), 2*nbasis, 2*nbasis)
      endif



C     %%%%%%%%%%%%%%%%
C     %  Save files  %
C     %%%%%%%%%%%%%%%%


C     Save Gaussian files.

      if ( itype .eq. 1 ) then
        call fileio (1, -irwca, 2*nbasis*norbs, v(jorb), 0)
        call fileio (1, -irwpa, 2*ntt, v(jp2), 0)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call fileio (1, -irwca, 2*nbasis*norbs, v(jorb), 0)
        call fileio (1, -irwpa, 2*ntt, v(jp2), 0)

        call fileio (1, -irwcb, 2*nbasis*norbs, v(jorb+2*nbsq), 0)
        call fileio (1, -irwpb, 2*ntt, v(jp2+2*ntt), 0)

      elseif ( itype .eq. 3 ) then
        call fileio (1, -irwca, 8*nbasis*norbs, v(jorb), 0)
        call fileio (1, -irwpa, 2*ntt2, v(jp2), 0)
      endif


      return
      end


