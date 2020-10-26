

      subroutine fixmt2 (iout, iprint, v, lenv, iopcl, nbasis, norbs,
     $     irws, irwx, irwh)

      implicit none


C +------------------------------------------------------------+
C |  fixmt2   --   CAJH, 08.2011                               |
C |                                                            |
C |                                                            |
C |  Simple routine that transforms the overlap, core          |
C |  Hamiltonian, and transformation (=S^(-1/2)) matrices to   |
C |  complex if they are currently real. That way we can       |
C |  avoid later problems in Gaussian.                         |
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

      real*8      v(*)
      integer     iout, iprint, lenv, iopcl, nbasis, norbs
      integer     irws, irwx, irwh

C     other variables

      integer     ntt, js, js2, jx, jx2, jh, jh2, jend



C     Return if iopcl is already complex.

      if ( iopcl .eq. 2 .or. iopcl .eq. 3 .or. iopcl .eq. 6 ) then
        return
      endif


C     Construct useful variables.

      ntt = nbasis * (nbasis + 1) / 2


C     Memory allocation.

C       Allocate space for:
C         - overlap matrix (js, js2)
C         - transformation matrix (jx, jx2)
C         - core Hamiltonian matrix (jh, jh2)

      js   = 1
      js2  = js   + 2*ntt
      jx   = js2  + ntt
      jx2  = jx   + 2*nbasis*norbs
      jh   = jx2  + nbasis*norbs
      jh2  = jh   + 2*ntt
      jend = jh2  + ntt

      call tstcor (jend-1, lenv, 'fixmt2')


C     Load matrices from RWF files.

      call fileio (2, -irws, ntt, v(js2), 0)
      call fileio (2, -irwx, nbasis*norbs, v(jx2), 0)
      call fileio (2, -irwh, ntt, v(jh2), 0)


C     Transform to complex.

      call zrecmp (3, v(js), v(js2), ntt)
      call zrecmp (3, v(jx), v(jx2), nbasis*norbs)
      call zrecmp (3, v(jh), v(jh2), ntt)


C     Fix size of RWF files.

      call conddf (irws, 2*ntt)
      call conddf (irwx, 2*nbasis*norbs)
      call conddf (irwh, 8*ntt)


C     Save complex matrices.

      call fileio (1, -irws, 2*ntt, v(js), 0)
      call fileio (1, -irwx, 2*nbasis*norbs, v(jx), 0)
      call fileio (1, -irwh, 2*ntt, v(jh), 0)
      call fileio (1, -irwh, 2*ntt, v(jh), 2*ntt)


      return
      end


