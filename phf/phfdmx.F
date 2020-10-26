      

      subroutine phfdmx (iout, iprint, v, lenv, iphf, norbs, irw1dm,
     $     irwp)

      implicit none


C +------------------------------------------------------------+
C |  phfdmx   --   CAJH, 08.2011                               |
C |                                                            |
C |                                                            |
C |  Simple routine that retrieves the PHF density matrix      |
C |  stored in the RWF file irw1dm, Hermitizes it, and         |
C |  stores it in lower triangular form in irwp. This will     |
C |  then allow it to be converted to the expected Gaussian    |
C |  format in the subroutine gaufmt.F.                        |
C |                                                            |
C |  Note that for ispuhf = .true., only the alpha-alpha       |
C |  and beta-beta blocks are stored in irwp. We assume here   |
C |  that the other blocks exactly vanish.                     |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       norbs  - number of orbitals
C       irw1dm - RWF file with PHF 1PDM
C       irwp   - RWF file with Hermitized PHF 1PDM  [ updated here ]

      real*8      v(*)
      integer     iout, iprint, lenv, iphf, norbs
      integer     irw1dm, irwp

C     other variables

      integer     nosq, ntto, ntt2o, itype, itype1
      integer     jp, jx, jzf, jzaa, jzbb, jzab, jzba, jscr, jend
      logical     icmplx, ispuhf, ispghf


C     Construct useful variables.

      nosq  = norbs * norbs
      ntto  = norbs * (norbs + 1) / 2
      ntt2o = norbs * (2*norbs + 1)


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
        call gauerr ('Unrecognized iphf in phfdmx.')
      endif


C     Memory allocation.

C       Allocate space to store:
C         - PHF density matrix (jp)
C         - scratch array (jx)

      if ( itype .eq. 1 ) then
        jp   = 1
        jx   = jp + 2*nosq
        jend = jx + 2*nosq

      elseif ( itype .eq. 2 .or. ispuhf ) then
        jp   = 1
        jx   = jp + 4*nosq
        jend = jx + 2*nosq

      elseif ( itype .eq. 3 ) then
        jp   = 1
        jx   = jp + 8*nosq
        jend = jx + 8*nosq
      endif

C       Allocate space to store:
C         - extra scratch arrays for ispghf = .true. (jzf, jzaa, ...)

      if ( ispghf ) then
        jzf  = jend
        jzaa = jzf  + 8*nosq
        jzbb = jzaa + 2*nosq
        jzab = jzbb + 2*nosq
        jzba = jzab + 2*nosq
        jscr = jzba + 2*nosq
        jend = jscr + 16*nosq
      endif

      call tstcor (jend-1, lenv, 'phfdmx')


C     Read density matrix from 1PDM RWF file.

      if ( itype1 .eq. 1 ) then
        call fileio (2, -irw1dm, 2*nosq, v(jp), 0)

      elseif ( itype1 .eq. 2 ) then
        call fileio (2, -irw1dm, 4*nosq, v(jp), 0)

      elseif ( itype1 .eq. 3 ) then
        call fileio (2, -irw1dm, 8*nosq, v(jp), 0)
      endif


C     Fill with zeros for ispghf.

      if ( ispghf ) then
        call aclear (2*nosq, v(jzab))
        call aclear (2*nosq, v(jzba))
        call amove (2*nosq, v(jp), v(jzaa))
        call amove (2*nosq, v(jp+2*nosq), v(jzbb))

        call dmblck (iout, iprint, v(jscr), 16*nosq, norbs, v(jzf),
     $       v(jzaa), v(jzab), v(jzba), v(jzbb), 1, -1)

        call amove (8*nosq, v(jzf), v(jp))
      endif


C     Hermitize density matrix and transform to lower triangular.

      if ( itype .eq. 1 ) then
        call apadag (1, norbs, 0.50d0, v(jp), v(jx))

        call trspnc (.false., norbs, norbs, v(jx))
        call zlinear (v(jx), v(jp), norbs, norbs)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call apadag (1, norbs, 0.50d0, v(jp), v(jx))

        call trspnc (.false., norbs, norbs, v(jx))
        call zlinear (v(jx), v(jp), norbs, norbs)

        call apadag (1, norbs, 0.50d0, v(jp+2*nosq), v(jx))

        call trspnc (.false., norbs, norbs, v(jx))
        call zlinear (v(jx), v(jp+2*ntto), norbs, norbs)

      elseif ( itype .eq. 3 ) then
        call apadag (1, 2*norbs, 0.50d0, v(jp), v(jx))

        call trspnc (.false., 2*norbs, 2*norbs, v(jx))
        call zlinear (v(jx), v(jp), 2*norbs, 2*norbs)
      endif


C     Save density matrix in output RWF file.

      if ( itype .eq. 1 ) then
        call fileio (1, -irwp, 2*ntto, v(jp), 0)
      elseif ( itype .eq. 2 .or. ispuhf ) then
        call fileio (1, -irwp, 4*ntto, v(jp), 0)
      elseif ( itype .eq. 3 ) then
        call fileio (1, -irwp, 2*ntt2o, v(jp), 0)
      endif


      return
      end


