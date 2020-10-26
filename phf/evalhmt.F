

      subroutine evalhmt (iout, iprint, v, lenv, iphf, nci, ncis, ncik, 
     $     ngrdb, nbasis, norbs, irwh, irwx, irwnos, ngrdt, indarr,
     $     grda, wgta, grdb, wgtb, grdg, wgtg, dmt, dind, yg, irwgg,
     $     irwpg, hmt)

      implicit none


C +------------------------------------------------------------+
C |  evalhmt   --   CAJH, 07.2011                              |
C |                                                            |
C |                                   (based on PHFB evalhmt)  |
C |                                                            |
C |                                                            |
C |  Evaluate the Hamiltonian matrix elements (Hmt) needed to  |
C |  solve the CI problem associated with complex conjugation  |
C |  restoration and spin projection.                          |
C |                                                            |
C |  All of these matrix elements can be written as            |
C |                                                            |
C |    Hmt(i,j) = int dg x(g; i,j) *                           |
C |                   { Hsp (g; i,j) + Hph (g; i,j) },         |
C |                                                            |
C |  where there is a parametric dependence on the i and j     |
C |  indices in the functions x(g), Hsp(g), Hph(g). The        |
C |  functions Hsp (g) and Hph (g) are given by                |
C |                                                            |
C |     Hsp (g; i,j) = Tr [ H P(g; i,j) ],                     |
C |     Hph (g; i,j) = 1/2 * Tr [ G(g; i,j) P(g; i,j) ].       |
C |                                                            |
C |  Here, H is the core Hamiltonian, P(g) are transition      |
C |  density matrices, and G(g) are the matrices resulting     |
C |  from the contraction of P(g) with electron-repulsion      |
C |  integrals.                                                |
C |                                                            |
C |    ( Everything is evaluated in what we call the NO basis  |
C |    throughout. This does not actually need to be the NO    |
C |    basis, but just some orthonormal basis as explained     |
C |    in evalsmt. )                                           |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  NOTE: This routine should be called right after evalsmt.  |
C |        The vector yg is loaded with the x(g) factors that  |
C |        we need here.                                       |
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
C       nbasis - number of basis functions
C       norbs  - number of orbitals
C       irwh   - RWF file with core Hamiltonian matrix
C       irwx   - RWF file with transformation matrix (=S^(-1/2))
C       irwnos - RWF file with transformation matrix from orthonormal AO
C                to NO basis
C       ngrdt  - total number of grid points
C       indarr - mapping array igrdt => igrdX
C       grdX   - integration grid for each variable
C       wgtX   - array with quadrature weights for each variable
C       dmt    - Wigner's small d-matrix d^S_{MK} (beta)
C       dind   - vector with indices {M} for Wigner rotation matrices
C       yg     - vector of x(g) values
C       irwgg  - RWF file with matrices G(g)
C       irwpg  - RWF file with matrices P(g)
C       hmt    - Hamiltonian matrix in CI problem [ out ]

      real*8      grda(*), grdb(*), grdg(*)
      real*8      wgta(*), wgtb(*), wgtg(*)
      integer     iout, iprint, lenv, iphf, nci, ncis, ncik, ngrdb
      integer     nbasis, norbs, ngrdt, indarr(ngrdt,*)
      integer     irwh, irwx, irwnos
      integer     irwgg, irwpg
      complex*16  yg(*), hmt(nci,nci)

      real*8      v(*), dmt(ngrdb,ncis,ncis), dind(ncis)

C     other variables

      integer     ntt, nosq, szp, szbfp, itype, iopcl, nri, ind
      integer     jh1, jh2, jx, jscr1, jend
      integer     jnos, jnost, jhno, jzab, jzba, jscr2, lscr2
      integer     lenp, npmem, npmax, npuse, ipar
      integer     n1, n2, n1c, n2c, n1s, n2s, ip, iblock, nthrd
      logical     icmplx, ispuhf, ispghf, throk
      complex*16  sumh

      real*8      facw
      integer     igrdt, iint, ibufg, ibufp, igg, ipg, iend
      integer     indy, indi, igrda, igrdb, igrdg
      complex*16  facd, trh, trg, zval

C     functions

      integer     nproc, lappar
      complex*16  trcabc

C     constants

      complex*16  zero, one, im, two


      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)
      im   = cmplx (0.0d0, 1.0d0)
      two  = cmplx (2.0d0, 0.0d0)



C     Construct useful variables.

      ntt  = nbasis * (nbasis + 1) / 2
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
        call gauerr ('Unrecognized iphf in evalhmt.')
      endif


C     Read iopcl from ILSW file.

C       iopcl = 0,  real RHF
C             = 1,  real UHF
C             = 2,  complex RHF
C             = 3,  complex UHF
C             = 6,  complex GHF (there is no real GHF)

      call ilsw (2, 1, iopcl)

      if ( iopcl .gt. 3 .and. iopcl .ne. 6 ) then
        call gauerr ('Incorrect iopcl in evalhmt.')
      endif

      nri = 2

      if ( iopcl .le. 1 ) nri = 1


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
C         - core Hamiltonian matrix (AO basis) (jh1)
C         - core Hamiltonian matrix (orthonormal basis) (jh2)
C         - AO => orthonormal AO transformation matrix (jx)
C         - scratch array to transform core Hamiltonian (jscr1)

      jh1   = 1
      jh2   = jh1   + nri*nbasis*nbasis
      jx    = jh2   + 2*nosq
      jscr1 = jx    + nri*nbasis*norbs
      jend  = jscr1 + nri*nbasis*norbs


C       Allocate space for:
C         - orthonormal AO => NO basis transformation matrix (jnos)
C         - scratch array to perform basis transformation (jnost)
C         - core Hamiltonian matrix in NO basis (jhno)

      if ( itype .eq. 1 ) then
        jnos  = jend
        jnost = jnos  + 2*nosq
        jhno  = jnost + 2*nosq
        jend  = jhno  + 2*nosq

      elseif ( itype .eq. 2 ) then
        jnos  = jend
        jnost = jnos  + 4*nosq
        jhno  = jnost + 2*nosq
        jend  = jhno  + 4*nosq

      elseif ( itype .eq. 3 ) then
        jnos  = jend
        jnost = jnos  + 8*nosq
        jhno  = jnost + 8*nosq
        jend  = jhno  + 8*nosq
      endif

C       Allocate some extra scratch arrays for itype = 3.

      if ( itype .eq. 3 ) then
        lscr2 = 16*nosq

        jzab  = jend
        jzba  = jzab  + 2*nosq
        jscr2 = jzba  + 2*nosq
        jend  = jscr2 + lscr2
      endif


C     Do some preparation for parallel processing:
C       - lenp,  dimension of scratch space per processor
C       - npuse, number of processors to use

      lenp = 2*nci*nci + 2*szbfp

      npmem = (lenv - jend + 1) / lenp

      if ( npmem .lt. 1 ) then
        call gauerr ('Not enough memory in evalhmt for a single CPU.')
      endif

      npmax = nproc (0)
      npuse = min (npmem, npmax)
      
      ipar  = lappar (0)
      throk = ipar .eq. 1

      if ( ipar .eq. -1 ) npuse = 1

      call tstcor (jend + lenp*npuse - 1, lenv, 'evalhmt')



C     Load core Hamiltonian matrix.
C     Transform into a square matrix.

      call fileio (2, -irwh, nri*ntt, v(jh1), 0)

      if ( nri .eq. 1 ) then
        call square (v(jh1), v(jh1), nbasis, nbasis, 0)
      elseif ( nri .eq. 2 ) then
        call zsquare (v(jh1), v(jh1), nbasis, nbasis, 0)
      endif


C     Load transformation matrix.

      call fileio (2, -irwx, nri*nbasis*norbs, v(jx), 0)


C     Load transformation matrix (orthonormal AO => NO basis).

      if ( itype .eq. 1 ) then
        call fileio (2, -irwnos, 2*nosq, v(jnos), 0)

      elseif ( itype .eq. 2 ) then
        call fileio (2, -irwnos, 4*nosq, v(jnos), 0)

      elseif ( itype .eq. 3 ) then
        call fileio (2, -irwnos, 8*nosq, v(jnos), 0)
      endif



C     Transform core Hamiltonian matrix into orthonormal AO basis by
C     computing the product X! . H . X.

      if ( nri .eq. 1 ) then
        call dgemm ('t', 'n', norbs, nbasis, nbasis, 1.0d0, v(jx),
     $       nbasis, v(jh1), nbasis, 0.0d0, v(jscr1), norbs)
        call dgemm ('n', 'n', norbs, norbs, nbasis, 1.0d0, v(jscr1),
     $       norbs, v(jx), nbasis, 0.0d0, v(jh1), norbs)

C       .. transform core Hamiltonian matrix into complex ..

        call zrecmp (3, v(jh2), v(jh1), nosq)

      elseif ( nri .eq. 2 ) then
        call zgemm ('c', 'n', norbs, nbasis, nbasis, one, v(jx),
     $       nbasis, v(jh1), nbasis, zero, v(jscr1), norbs)
        call zgemm ('n', 'n', norbs, norbs, nbasis, one, v(jscr1),
     $       norbs, v(jx), nbasis, zero, v(jh1), norbs)

        call amove (2*nosq, v(jh1), v(jh2))
      endif


C     debugging...
C       Print core Hamiltonian in orthonormal basis.

       write (iout, *) ' '
       call ltoutx (iout, 1, 2, 1, 0, 'Core Hamiltonian (ortho)', norbs,
     $      norbs, v(jh2), v(jh2), 0)


C     Transform core Hamiltonian into NO basis:
C       nos! . H . nos

      if ( itype .eq. 1 ) then

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jh2), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos), norbs, zero, v(jhno), norbs)

      elseif ( itype .eq. 2 ) then

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jh2), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos), norbs, zero, v(jhno), norbs)

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos+2*nosq),
     $       norbs, v(jh2), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos+2*nosq), norbs, zero, v(jhno+2*nosq), norbs)

      elseif ( itype .eq. 3 ) then

C       .. form core Hamiltonian in spin orbital basis ..

        call aclear (2*nosq, v(jzab))
        call aclear (2*nosq, v(jzba))

        call dmblck (iout, iprint, v(jscr2), lscr2, norbs, v(jhno),
     $       v(jh2), v(jzab), v(jzba), v(jh2), 1, -1)

C       .. perform the transformation ..

        call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $       v(jnos), 2*norbs, v(jhno), 2*norbs, zero, v(jnost),
     $       2*norbs)
        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $       v(jnost), 2*norbs, v(jnos), 2*norbs, zero, v(jhno),
     $       2*norbs)
      endif


C     debugging...
C       Print core Hamiltonian in NO basis.

c$$$       if ( itype .eq. 1 ) then
c$$$         write (iout, *) ' '
c$$$         call ltoutx (iout, 1, 2, 1, 0, 'Core Hamiltonian (NO)', norbs,
c$$$     $        norbs, v(jhno), v(jhno), 0)
c$$$
c$$$       elseif ( itype .eq. 2 ) then
c$$$         write (iout, *) ' '
c$$$         call ltoutx (iout, 1, 3, 1, 0, 'Core Hamiltonian (NO)', norbs,
c$$$     $        norbs, v(jhno), v(jhno+2*nosq), 0)
c$$$
c$$$       elseif ( itype .eq. 3 ) then
         write (iout, *) ' '
         call ltoutx (iout, 1, 2, 1, 0, 'Core Hamiltonian (NO)',
     $        2*norbs, 2*norbs, v(jhno), v(jhno), 0)
c$$$       endif


C     Take complex conjugate of core Hamiltonian (to take traces
C     appropriately).

      if ( itype .eq. 1 ) then
        call vconjg (nosq, v(jhno))

      elseif ( itype .eq. 2 ) then
        call vconjg (nosq, v(jhno))
        call vconjg (nosq, v(jhno+2*nosq))

      elseif ( itype .eq. 3 ) then
        call vconjg (4*nosq, v(jhno))
      endif



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Evaluate the Hamiltonian matrix elements  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Evaluate the Hamiltonian matrix elements as described above...


C     Clear iint values for all processors.

      do 301 ip = 1, npuse
        iint = jend + (ip - 1) * lenp

        call aclear (2*nci*nci, v(iint))
 301  continue


C     Loop over grid points.

      do 103 iblock = 1, ngrdt, npuse
        nthrd = min (npuse, ngrdt - iblock + 1)


C     Load matrices. Loop not parallelized.

        do 104 ip = 1, nthrd

          igrdt = iblock + ip - 1

C     Memory allocation per processor:
C       iint  - integrated quantity x(g) * { Hsp (g) + Hph (g) }
C       ibufg - G(g) buffer at each grid point
C       ibufp - P(g) buffer at each grid point

          iint  = jend  + (ip - 1) * lenp
          ibufg = iint  + 2*nci*nci
          ibufp = ibufg + szbfp
          iend  = ibufp + szbfp

          ind = (igrdt - 1)*szbfp

          call fileio (2, -irwgg, szbfp, v(ibufg), ind)
          call fileio (2, -irwpg, szbfp, v(ibufp), ind)
 104    continue


C     Computation loop. Parallelized.

        throk = throk .and. npuse .gt. 1

C$omp   parallel do if (throk) schedule(static,1) default(shared)
C$omp+  private (igrdt, iint, ibufg, ibufp, igg, ipg, iend)
C$omp+  private (igrda, igrdb, igrdg)
C$omp+  private (facw, facd, trh, trg, zval)
C$omp+  private (n1c, n2c, n1s, n2s, indy, n1, n2, indi)

        do 105 ip = 1, nthrd

          igrdt = iblock + ip - 1

C     Memory allocation per processor:
C       iint  - integrated quantity x(g) * { Hsp (g) + Hph (g) }
C       ibufg - G(g) buffer at each grid point
C       ibufp - P(g) buffer at each grid point

          iint  = jend  + (ip - 1) * lenp
          ibufg = iint  + 2*nci*nci
          ibufp = ibufg + szbfp
          iend  = ibufp + szbfp

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

              igg = (n1c-1)*ncik*szp + (n2c-1)*szp + ibufg
              ipg = (n1c-1)*ncik*szp + (n2c-1)*szp + ibufp
              indy = (igrdt-1)*ncik*ncik + (n1c-1)*ncik + n2c


C     Compute the trace
C       tr ( H . P(g) )
C     and add it to iint.

C       ** trcabc computes tr ( a . b* )

              if ( itype .eq. 1 ) then
                trh = trcabc (1, norbs, v(ipg), v(jhno))
                trh = two * trh

              elseif ( itype .eq. 2 ) then
                trh = trcabc (1, norbs, v(ipg), v(jhno))
                trh = trh + 
     $                trcabc (1, norbs, v(ipg+2*nosq), v(jhno+2*nosq))

              elseif ( itype .eq. 3 ) then
                trh = trcabc (1, 2*norbs, v(ipg), v(jhno))
              endif

              zval = trh * yg(indy)


C     Compute the trace
C       1/2 . tr ( G(g) . P(g) )
C     and add it to iint.

              if ( itype .eq. 1 ) then
                call vconjg (nosq, v(igg))

                trg = trcabc (1, norbs, v(ipg), v(igg))

              elseif ( itype .eq. 2 ) then
                call vconjg (nosq, v(igg))
                call vconjg (nosq, v(igg+2*nosq))

                trg = trcabc (1, norbs, v(ipg), v(igg))
                trg = trg + 
     $                trcabc (1, norbs, v(ipg+2*nosq), v(igg+2*nosq))
                trg = trg / two

              elseif ( itype .eq. 3 ) then
                call vconjg (4*nosq, v(igg))

                trg = trcabc (1, 2*norbs, v(ipg), v(igg))
                trg = trg / two
              endif

              zval = zval + trg * yg(indy)
              write(iout, 666) zval
666   Format('zval = ', F10.5, F10.5)


C     Loop over CI matrix elements (spin projection).

              do 201 n1s = 1, ncis
                do 202 n2s = 1, ncis

                  n1 = (n1c-1)*ncis + n1s
                  n2 = (n2c-1)*ncis + n2s
                  indi = (n1-1)*nci*2 + (n2-1)*2 + iint

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

                  v(indi)   = v(indi) + facw * dble (zval * facd)
                  v(indi+1) = v(indi+1) + facw * aimag (zval * facd)

 202            continue
 201          continue

 102        continue
 101      continue
 105    continue

C$omp   end parallel do

 103  continue


C     Accumulate all integrals.

      do 401 n1 = 1, nci
        do 402 n2 = 1, nci
          sumh = zero

          do 302 ip = 1, npuse
            iint = jend + (ip - 1) * lenp
            indi = (n1-1)*nci*2 + (n2-1)*2 + iint
            sumh = sumh + cmplx (v(indi), v(indi+1))
 302      continue

C     Save overlap matrix elements.

          hmt(n1,n2) = sumh
 402    continue
 401  continue



      return
      end


