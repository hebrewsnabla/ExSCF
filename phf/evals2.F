

      subroutine evals2 (iout, iprint, v, lenv, iphf, nci, ncis, ncik,
     $     ngrdb, norbs, ngrdt, indarr, grda, wgta, grdb, wgtb, grdg,
     $     wgtg, dmt, dind, yg, irwnos, irwpg, civec, s2)

      implicit none


C +------------------------------------------------------------+
C |  evals2   --   CAJH, 07.2011                               |
C |                                                            |
C |                                    (based on PHFB evals2)  |
C |                                                            |
C |                                                            |
C |  Compute the expectation value <S**2> for a PHF            |
C |  wavefunction.                                             |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  Before giving explicit expressions for the expectation    |
C |  values, let us recall what the PHF 1PDM and 2PDM look     |
C |  like:                                                     |
C |                                                            |
C |    [1PDM]^j_i = int dg y(g) [P(g)]^j_i                     |
C |                                                            |
C |    [2PDM]^{kl}_{ij} = 1/2 * int dg y(g)                    |
C |                           { [P(g)]^k_i [P(g)]^l_j -        |
C |                             [P(g)]^l_i [P(g)]^k_j }        |
C |                                                            |
C |  where P(g) are transition density matrices (built in the  |
C |  subroutine formpg.F).                                     |
C |                                                            |
C |  This form of the 1PDM and 2PDM can be easily verified     |
C |  by observing that it is exactly these quantities that     |
C |  are traced with one- and two-electron integrals in the    |
C |  evaluation of the PHF energy.                             |
C |                                                            |
C |  ** The provided expression for the 2PDM is only           |
C |     appropriate for scalar operators that commute with     |
C |     the projection operators considered.                   |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  The spin blocks of the density matrix can be written as   |
C |                                                            |
C |    P  =  ( Pc + Mz      Mx - i*My )                        |
C |          ( Mx + i*My     Pc - Mz  ),                       |
C |                                                            |
C |  with                                                      |
C |                                                            |
C |    Pc =  1/2 * (Paa + Pbb)                                 |
C |    Mz =  1/2 * (Paa - Pbb)                                 |
C |    My =  i/2 * (Pab - Pba)                                 |
C |    Mx =  1/2 * (Pab + Pba)                                 |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  Here, I provide the appropriate expressions for <S**2>.   |
C |  All the M, Pc, and K matrices are assumed to be           |
C |  gauge-dependent and thus a single integration over the    |
C |  gauge angle needs to be performed at the end to obtain    |
C |  the actual value.                                         |
C |                                                            |
C |  The expectation value of <S**2> is given by:              |
C |                                                            |
C |    <S**2> = trace^2 (Mx) + trace^2 (My) + trace^2 (Mz) +   |
C |             3/2 * trace (Pc - Pc^2) +                      |
C |             1/2 * trace (Mx^2 + My^2 + Mz^2)               |
C |                                                            |
C |  A derivation of the above expression should appear in     |
C |  the literature in the near future...                      |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  CAJH, 12.2012                                             |
C |                                                            |
C |  Reworked subroutine to handle the CI over spin            |
C |  internally.                                               |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       nci    - number of basis states in CI problem
C       ncis   - number of basis states for spin projection
C       ncik   - number of basis states for complex conj projection
C       ngrdb  - number of grid points in beta integration
C       norbs  - number of orbitals
C       ngrdt  - total number of grid points
C       indarr - mapping array igrdt => igrdX
C       grdX   - integration grid for each variable
C       wgtX   - array with quadrature weights for each variable
C       dmt    - Wigner's small d-matrix d^S_{MK} (beta)
C       dind   - vector with indices {M} for Wigner rotation matrices
C       yg     - vector of y(g) values
C       irwnos - RWF file with transformation matrix from orthonormal AO
C                to NO basis
C       irwpg  - RWF file with matrices P(g)
C       civec  - vector of CI coefficients
C       s2     - expectation value < S**2 > [ out ]

      real*8      grda(*), grdb(*), grdg(*)
      real*8      wgta(*), wgtb(*), wgtg(*)
      integer     iout, iprint, lenv, iphf, nci, ncis, ncik, ngrdb
      integer     norbs, ngrdt, indarr(ngrdt,*)
      integer     irwnos, irwpg
      complex*16  yg(*), civec(*)

      real*8      v(*), dmt(ngrdb,ncis,ncis), dind(ncis), s2

C     other variables

      integer     nosq, itype, szp, szbfp, ind
      integer     jnos, jend, lscrx
      integer     lenp, npmem, npmax, npuse, ipar
      integer     n1, n2, n1c, n2c, n1s, n2s, ip, iblock, nthrd
      logical     icmplx, ispuhf, ispghf, throk
      complex*16  sums2

      real*8      facw
      integer     igrdt, ibufp, ipg, is2, iend
      integer     ipa, ipb, ipab, ipba
      integer     ip0, imx, imy, imz
      integer     iscrx, it, j, ind1, ind2, ind3, ind4
      integer     indy, igrda, igrdb, igrdg
      complex*16  fac, facd, facs2
      complex*16  trp0, trmx, trmy, trmz
      complex*16  trp02, trmx2, trmy2, trmz2

C     functions

      integer     nproc, lappar
      complex*16  trcabc

C     constants

      complex*16  zero, one, f1o2, f3o2, im


      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)
      f1o2 = cmplx (0.5d0, 0.0d0)
      f3o2 = cmplx (1.5d0, 0.0d0)
      im   = cmplx (0.0d0, 1.0d0)


C     Construct useful variables.

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
        call gauerr ('Unrecognized iphf in evals2.')
      endif


C     Define useful array dimensions:
C       - szp,  size of P(g)

      if ( itype .eq. 1 ) then
        szp = 2*nosq
      elseif ( itype .eq. 2 ) then
        szp = 4*nosq
      elseif ( itype .eq. 3 ) then
        szp = 8*nosq
      endif

C       .. buffer to load computed Pg ..

      szbfp = ncik*ncik * szp


C     Memory allocation.

C       Allocate space for:
C         - orthonormal AO => NO basis transformation matrix (jnos)

      if ( itype .eq. 1 ) then
        jnos = 1
        jend = jnos + 2*nosq

      elseif ( itype .eq. 2 ) then
        jnos = 1
        jend = jnos + 4*nosq

      elseif ( itype .eq. 3 ) then
        jnos = 1
        jend = jnos + 8*nosq
      endif


C     Do some preparation for parallel processing:
C       - lenp,  dimension of scratch space per processor
C       - npuse, number of processors to use

      lscrx = 0

      if ( itype .eq. 1 ) then
        lenp = 2 + 6*nosq

      elseif ( itype .eq. 2 ) then
        lenp = 2 + 10*nosq

      elseif ( itype .eq. 3 ) then
        lscrx = 16*nosq

        lenp = 2 + lscrx + 24*nosq
      endif

      lenp = lenp + szbfp

      npmem = (lenv - jend + 1) / lenp

      if ( npmem .lt. 1 ) then
        call gauerr ('Not enough memory in evals2 for a single CPU.')
      endif

      npmax = nproc (0)
      npuse = min (npmem, npmax)
      
      ipar  = lappar (0)
      throk = ipar .eq. 1

      if ( ipar .eq. -1 ) npuse = 1

      call tstcor (jend + lenp*npuse - 1, lenv, 'evals2')


C     Load transformation matrix (orthonormal AO => NO basis).

      if ( itype .eq. 1 ) then
        call fileio (2, -irwnos, 2*nosq, v(jnos), 0)

      elseif ( itype .eq. 2 ) then
        call fileio (2, -irwnos, 4*nosq, v(jnos), 0)

      elseif ( itype .eq. 3 ) then
        call fileio (2, -irwnos, 8*nosq, v(jnos), 0)
      endif



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Compute expectation values  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Evaluate the expectation value of S^2 as described above.


C     Clear accumulated quantities for all processors.

      do 11 ip = 1, npuse
        is2  = jend + (ip - 1) * lenp

        call aclear (2, v(is2))
 11   continue


C     Loop over grid points.

      do 103 iblock = 1, ngrdt, npuse
        nthrd = min (npuse, ngrdt - iblock + 1)


C     Load matrices. Loop not parallelized.

        do 104 ip = 1, nthrd

          igrdt = iblock + ip - 1

C     Memory allocation per processor:
C       is2   - accumulated quantity <S**2>
C       ibufp - P(g) buffer at each grid point

          is2   = jend  + (ip - 1) * lenp
          ibufp = is2   + 2
          iend  = ibufp + szbfp

          ind = (igrdt - 1)*szbfp

          call fileio (2, -irwpg, szbfp, v(ibufp), ind)
 104    continue


C     Computation loop. Parallelized.

        throk = throk .and. npuse .gt. 1

C$omp   parallel do if (throk) schedule(static,1) default(shared)
C$omp+  private (igrdt, ibufp, ipg, iend, is2)
C$omp+  private (ipa, ipb, ipab, ipba)
C$omp+  private (ip0, imx, imy, imz)
C$omp+  private (iscrx, it, fac, facw, facd, facs2)
C$omp+  private (j, ind1, ind2, ind3, ind4)
C$omp+  private (igrda, igrdb, igrdg)
C$omp+  private (trp0, trmx, trmy, trmz)
C$omp+  private (trp02, trmx2, trmy2, trmz2)
C$omp+  private (n1c, n2c, n1s, n2s, indy, n1, n2)

        do 105 ip = 1, nthrd

          igrdt = iblock + ip - 1

C     Memory allocation per processor:
C       is2   - accumulated quantity <S**2>
C       ibufp - P(g) buffer at each grid point

          is2   = jend  + (ip - 1) * lenp
          ibufp = is2   + 2
          iend  = ibufp + szbfp


C       ipa, ipb, ... - spin blocks of ipg
C       ip0, imz, ... - charge and spin density matrices

          if ( itype .eq. 1 ) then
            ipa  = iend
            ip0  = ipa  + 2*nosq
            iend = ip0  + 2*nosq

          elseif ( itype .eq. 2 ) then
            ipa  = iend
            ipb  = ipa  + 2*nosq
            ip0  = ipb  + 2*nosq
            imz  = ip0  + 2*nosq
            iend = imz  + 2*nosq

          elseif ( itype .eq. 3 ) then
            ipa   = iend
            ipb   = ipa  + 2*nosq
            ipab  = ipb  + 2*nosq
            ipba  = ipab + 2*nosq
            ip0   = ipba + 2*nosq
            imz   = ip0  + 2*nosq
            imy   = imz  + 2*nosq
            imx   = imy  + 2*nosq
            iend  = imx  + 2*nosq
          endif

C       it    - scratch array to transform P(g) to NO basis
C       iscrx - scratch array to split P(g) into spin blocks

          if ( itype .eq. 1 ) then
            it   = iend
            iend = it   + 2*nosq

          elseif ( itype .eq. 2 ) then
            it   = iend
            iend = it   + 2*nosq

          elseif ( itype .eq. 3 ) then
            it    = iend
            iscrx = it    + 8*nosq
            iend  = iscrx + lscrx
          endif

C     Retrieve some useful factors.

          igrda = indarr(igrdt,1)
          igrdb = indarr(igrdt,2)
          igrdg = indarr(igrdt,3)


C     Compute useful factor
C       facw = wgt(g)

          facw = wgta(igrda) * wgtb(igrdb) * wgtg(igrdg)


C     Loop over CI matrix elements (complex conjugation).

          do 101 n1c = 1, ncik
            do 102 n2c = 1, ncik

              ipg = (n1c-1)*ncik*szp + (n2c-1)*szp + ibufp
              indy = (igrdt-1)*ncik*ncik + (n1c-1)*ncik + n2c


C     Transform P(g) from NO basis to orthonormal AO basis:
C       nos . P(g) . nos!

              if ( itype .eq. 1 ) then
                call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnos),
     $               norbs, v(ipg), norbs, zero, v(it), norbs)
                call zgemm ('n', 'c', norbs, norbs, norbs, one, v(it),
     $               norbs, v(jnos), norbs, zero, v(ipg), norbs)

              elseif ( itype .eq. 2 ) then
                call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnos),
     $               norbs, v(ipg), norbs, zero, v(it), norbs)
                call zgemm ('n', 'c', norbs, norbs, norbs, one, v(it),
     $               norbs, v(jnos), norbs, zero, v(ipg), norbs)

                call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $               v(jnos+2*nosq), norbs, v(ipg+2*nosq), norbs, zero,
     $               v(it), norbs)
                call zgemm ('n', 'c', norbs, norbs, norbs, one, v(it),
     $               norbs, v(jnos+2*nosq), norbs, zero, v(ipg+2*nosq),
     $               norbs)

              elseif ( itype .eq. 3 ) then
                call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $               v(jnos), 2*norbs, v(ipg), 2*norbs, zero, v(it),
     $               2*norbs)
                call zgemm ('n', 'c', 2*norbs, 2*norbs, 2*norbs, one,
     $               v(it), 2*norbs, v(jnos), 2*norbs, zero, v(ipg),
     $               2*norbs)
              endif


C     Split density matrix into spin blocks.

              if ( itype .eq. 1 ) then
                call amove (2*nosq, v(ipg), v(ipa))

              elseif ( itype .eq. 2 ) then
                call amove (2*nosq, v(ipg), v(ipa))
                call amove (2*nosq, v(ipg+2*nosq), v(ipb))

              elseif ( itype .eq. 3 ) then
                call dmblck (iout, iprint, v(iscrx), lscrx, norbs,
     $               v(ipg), v(ipa), v(ipab), v(ipba), v(ipb), 2, -1)
              endif


C     Form charge density P0 and spin density matrices (Mx, My, Mz).

              if ( itype .eq. 1 ) then
                call amove (2*nosq, v(ipa), v(ip0))

              elseif ( itype .eq. 2 ) then
                call aadd (2*nosq, v(ipa), v(ipb), v(ip0))
                call asub (2*nosq, v(ipa), v(ipb), v(imz))

                call zdscal (nosq, 0.50d0, v(ip0), 1)
                call zdscal (nosq, 0.50d0, v(imz), 1)

              elseif ( itype .eq. 3 ) then
                call aadd (2*nosq, v(ipa), v(ipb), v(ip0))
                call asub (2*nosq, v(ipa), v(ipb), v(imz))
                call asub (2*nosq, v(ipab), v(ipba), v(imy))
                call aadd (2*nosq, v(ipab), v(ipba), v(imx))

                call zdscal (nosq, 0.50d0, v(ip0), 1)
                call zdscal (nosq, 0.50d0, v(imz), 1)
                call zdscal (nosq, 0.50d0, v(imy), 1)
                call zdscal (nosq, 0.50d0, v(imx), 1)

                call zscal (nosq, im, v(imy), 1)
              endif


C     Compute trace of P0, Mx, My, Mz.

              trp0 = zero
              trmx = zero
              trmy = zero
              trmz = zero

              do 200 j = 1, norbs
                ind1 = ip0 + (j-1) * 2*norbs + (j-1) * 2

                if ( itype .eq. 2 ) then
                  ind2 = imz + (j-1) * 2*norbs + (j-1) * 2

                elseif ( itype .eq. 3 ) then
                  ind2 = imz + (j-1) * 2*norbs + (j-1) * 2
                  ind3 = imy + (j-1) * 2*norbs + (j-1) * 2
                  ind4 = imx + (j-1) * 2*norbs + (j-1) * 2
                endif

                trp0 = trp0 + cmplx (v(ind1), v(ind1+1))

                if ( itype .eq. 2 ) then
                  trmz = trmz + cmplx (v(ind2), v(ind2+1))

                elseif ( itype .eq. 3 ) then
                  trmz = trmz + cmplx (v(ind2), v(ind2+1))
                  trmy = trmy + cmplx (v(ind3), v(ind3+1))
                  trmx = trmx + cmplx (v(ind4), v(ind4+1))
                endif
 200          continue
              write(iout, '8(F12.6)') trp0, trmx, trmy, trmz


C     Compute
C       trp02 = tr ( P0 . P0 ),
C       trmx2 = tr ( Mx . Mx ),
C       trmy2 = tr ( My . My ),
C       trmz2 = tr ( Mz . Mz ).

              call amove (2*nosq, v(ip0), v(ipa))
              call vconjg (nosq, v(ipa))

              trp02 = trcabc (1, norbs, v(ip0), v(ipa))

              if ( itype .eq. 2 ) then
                call amove (2*nosq, v(imz), v(ipb))
                call vconjg (nosq, v(ipb))

                trmz2 = trcabc (1, norbs, v(imz), v(ipb))

              elseif ( itype .eq. 3 ) then
                call amove (2*nosq, v(imz), v(ipb))
                call amove (2*nosq, v(imy), v(ipab))
                call amove (2*nosq, v(imx), v(ipba))
                call vconjg (nosq, v(ipb))
                call vconjg (nosq, v(ipab))
                call vconjg (nosq, v(ipba))

                trmz2 = trcabc (1, norbs, v(imz), v(ipb))
                trmy2 = trcabc (1, norbs, v(imy), v(ipab))
                trmx2 = trcabc (1, norbs, v(imx), v(ipba))
              
              write(iout, '8(F10.5)') trp02, trmx2, trmy2, trmz2
              endif


C     Add contributions to <S**2>.

              if ( itype .eq. 1 ) then
                facs2 = f3o2 * (trp0 - trp02)

              elseif ( itype .eq. 2 ) then
                facs2 = trmz**2 +
     $                  f1o2 * trmz2 +
     $                  f3o2 * (trp0 - trp02)

              elseif ( itype .eq. 3 ) then
                facs2 = trmx**2 + trmy**2 + trmz**2 +
     $                  f1o2 * (trmx2 + trmy2 + trmz2) +
     $                  f3o2 * (trp0 - trp02)
              endif


C     Loop over CI matrix elements (spin projection).

              do 201 n1s = 1, ncis
                do 202 n2s = 1, ncis

                  n1 = (n1c-1)*ncis + n1s
                  n2 = (n2c-1)*ncis + n2s


C     Compute Wigner's rotation matrix factor
C       facd = ( D^S_{MK} (Omega) )*

C       Note that
C         ( D^S_{MK} (Omega) )*  =  exp (+i * alpha * M) *
C                                   exp (+i * gamma * K) *
C                                   d^S_{MK} (beta),

C       where d^S_{MK} (beta) is Wigner's small d-matrix.

                  if ( ispghf ) then
                    facd = exp (+im * cmplx (grda(igrda) * dind(n1s))) *
     $                     exp (+im * cmplx (grdg(igrdg) * dind(n2s))) *
     $                     cmplx (dmt(igrdb,n1s,n2s))

                  elseif ( ispuhf ) then
                    facd = cmplx (dmt(igrdb,1,1))

                  else
                    facd = one
                  endif

                  fac = yg(indy) * facd *
     $                  conjg(civec(n1)) * civec(n2)

                  v(is2) = v(is2) + facw * dble (facs2 * fac)
                  v(is2+1) = v(is2+1) + facw * aimag (facs2 * fac)

 202            continue
 201          continue

 102        continue
 101      continue

 105    continue

C$omp   end parallel do

 103  continue


C     Merge results from all processors.

      sums2 = zero

      do 21 ip = 1, npuse
        is2  = jend + (ip - 1) * lenp

        sums2 = sums2 + cmplx (v(is2), v(is2+1))
 21   continue


C     Convert to real numbers.

      s2 = dble (sums2)


      return
      end


