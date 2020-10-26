

      subroutine formyg (iout, iprint, v, lenv, iphf, nci, ncis, ncik,
     $     ngrdb, norbs, na, nb, ngrdt, indarr, grda, wgta, grdb, wgtb,
     $     grdg, wgtg, dmt, dind, civec, irwnos, irwng, irwyg, intxg,
     $     yg, pno, pstrno, inobas)

      implicit none


C +------------------------------------------------------------+
C |  formyg   --   CAJH, 07.2011                               |
C |                                                            |
C |                                                            |
C |  Build the matrices Y(g), related to the derivatives of    |
C |  the weighting function y(g) with respect to elements of   |
C |  the HF density matrix:                                    |
C |                                                            |
C |    d (y(g)) / d(P_{kl})  ==  y(g) [Y(g)]_{lk}              |
C |                                                            |
C |  Before giving actual expressions for Y(g), it is          |
C |  convenient to introduce the matrices X(g), defined by     |
C |                                                            |
C |    d (x(g)) / d(P_{kl})  ==  x(g) [X(g)]_{lk},             |
C |                                                            |
C |  where x(g) are the overlap-related quantities (see        |
C |  evalsmt for details).                                     |
C |                                                            |
C |  Given X(g), the matrices Y(g) can be constructed as       |
C |                                                            |
C |    Y(g) = X(g) - int dg y(g) X(g).                         |
C |                                                            |
C |                                                            |
C |  We now provide the explicit forms of the matrices X(g).   |
C |  The Dg matrix introduced is a rotation matrix related     |
C |  to the spin projection operator with an angle Omega.      |
C |  The X(g) matrices are:                                    |
C |                                                            |
C |                                                            |
C |  if nu = I, mu = I  ( <0| I! R(Omega) I |0> )              |
C |                                                            |
C |    X(g)  =  Dg . ( P_oo ) . N(g)  +                        |
C |                  ( P_vo )                                  |
C |                                                            |
C |             N(g) . ( P_oo   P_ov ) . Dg                    |
C |                                                            |
C |                                                            |
C |  if nu = K, mu = I  ( <0| I! R(Omega) K |0> )              |
C |                                                            |
C |    X(g)  =  Dg . ( P_oo  P_ov )^T . N(g)  +                |
C |                                                            |
C |             [ N(g) . ( P_oo   P_ov ) . Dg ]^T              |
C |                                                            |
C |                                                            |
C |  if nu = I, mu = K  ( <0| K! R(Omega) I |0> )              |
C |                                                            |
C |    X(g)  =  [ Dg . ( P_oo ) . N(g) ]^T  +                  |
C |             [      ( P_vo )        ]                       |
C |                                                            |
C |             N(g) . ( P_oo )^T . Dg                         |
C |                    ( P_vo )                                |
C |                                                            |
C |                                                            |
C |  if nu = K, mu = K  ( <0| K! R(Omega) K |0> )              |
C |                                                            |
C |    X(g)  =  [ Dg . ( P_oo  P_ov )^T . N(g) ]^T  +          |
C |                                                            |
C |             [ N(g) . ( P_oo )^T . Dg ]^T                   |
C |             [        ( P_vo )        ]                     |
C |                                                            |
C |                                                            |
C |  where nu, mu = I, K are operators related to complex      |
C |  conjugation restoration. (See evalsmt for more details.)  |
C |  The matrices N(g) are those built in formng. Lastly, P    |
C |  is the density matrix corresponding to the deformed HF    |
C |  state.                                                    |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  CAJH, 08.2011                                             |
C |                                                            |
C |  Use of the true NO basis has been enabled even for        |
C |  methods with complex conjugation.                         |
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
C       na     - number of alpha electrons
C       nb     - number of beta electrons
C       ngrdt  - total number of grid points
C       indarr - mapping array igrdt => igrdX
C       grdX   - integration grid for each variable
C       wgtX   - array with quadrature weights for each variable
C       dmt    - Wigner's small d-matrix d^S_{MK} (beta)
C       dind   - vector with indices {M} for Wigner rotation matrices
C       civec  - vector of CI coefficients
C       irwnos - RWF file with transformation matrix from orthonormal AO
C                to NO basis
C       irwng  - RWF file with matrices N(g)
C       irwyg  - RWF file with matrices Y(g) [ updated here ]
C       intxg  - integrated matrix X(g) [ updated here ]
C       yg     - vector of y(g) values
C       pno    - density matrix in NO basis
C       pstrno - density matrix (P*) in NO basis
C       inobas - if .true., force the use of the true NO basis

      real*8      grda(*), grdb(*), grdg(*)
      real*8      wgta(*), wgtb(*), wgtg(*)
      integer     iout, iprint, lenv, iphf, nci, ncis, ncik
      integer     ngrdb, na, nb
      integer     norbs, ngrdt, indarr(ngrdt,*)
      integer     irwnos, irwng, irwyg
      logical     inobas

      real*8      v(*), dmt(ngrdb,ncis,ncis), dind(ncis)
      real*8      pno(*), pstrno(*)
      complex*16  yg(*), civec(*), intxg(*)

C     other variables

      integer     nel, nosq, itype, szn, szp, szbfn, szbfp, ind
      integer     jpnt, jptn, jpant, jpatn, jpbnt, jpbtn
      integer     jp2nt, jp2tn, jp2ant, jp2atn, jp2bnt, jp2btn
      integer     lscrdg, jnos, jy, ju, jend
      integer     lenp, npmem, npmax, npuse, ipar
      integer     n1, n2, n1c, n2c, n1s, n2s, ip, iblock, nthrd
      logical     icmplx, ispuhf, ispghf, throk
      character   indc1*5, indc2*2, indc3*2, str*14

      real*8      facw
      integer     igrdt, ibufn, ibufp, ing, inga, ingb
      integer     iy, iyg, ix, iend
      integer     it1, it1a, it1b, it2, it2a, it2b
      integer     iz1f, iz2f, iscrdg, ix2, indy
      integer     igrda, igrdb, igrdg
      real*8      anga, angb, angg
      complex*16  facd, fac

C     functions

      integer     nproc, lappar

C     constants

      complex*16  zero, one, im


      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)
      im   = cmplx (0.0d0, 1.0d0)


C     Construct useful variables.

      nel = na + nb
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
        call gauerr ('Unrecognized iphf in formyg.')
      endif


C     Define useful array dimensions:
C       - szn,  size of N(g) matrices
C       - szp,  size of Y(g) matrices

      if ( itype .eq. 1 ) then
        szn = 2*na*na
        szp = 2*nosq

      elseif ( itype .eq. 2 ) then
        szn = 2*na*na + 2*nb*nb
        szp = 4*nosq

      elseif ( itype .eq. 3 ) then
        szn = 2*nel*nel
        szp = 8*nosq
      endif

C       .. buffer to load computed Ng ..
C       .. buffer to save computed Yg ..

      szbfn = ncik*ncik * szn
      szbfp = ncik*ncik * szp


C     Memory allocation.

C       Allocate space for:
C         - orthonormal AO => NO basis transformation matrix (jnos)
C         - rectangular blocks of density matrix in NO basis (jpnt, jptn)
C         - integrated quantity: int dg X(g) (jy)
C         - rectangular blocks of density matrix (P*) in NO basis (jp2*)

      if ( itype .eq. 1 ) then
        jpnt = 1
        jptn = jpnt + 2*na*norbs
        jy   = jptn + 2*na*norbs
        jend = jy   + 2*nosq

      elseif ( itype .eq. 2 ) then
        jpant = 1
        jpatn = jpant + 2*na*norbs
        jpbnt = jpatn + 2*na*norbs
        jpbtn = jpbnt + 2*nb*norbs
        jy    = jpbtn + 2*nb*norbs
        jend  = jy    + 4*nosq

      elseif ( itype .eq. 3 ) then
        jpnt = 1
        jptn = jpnt + 4*nel*norbs
        jy   = jptn + 4*nel*norbs
        jend = jy   + 8*nosq
      endif

      if ( inobas .and. itype .eq. 1 ) then
        jp2nt = jend
        jp2tn = jp2nt + 2*na*norbs
        jend  = jp2tn + 2*na*norbs

      elseif ( inobas .and. itype .eq. 2 ) then
        jp2ant = jend
        jp2atn = jp2ant + 2*na*norbs
        jp2bnt = jp2atn + 2*na*norbs
        jp2btn = jp2bnt + 2*nb*norbs
        jend   = jp2btn + 2*nb*norbs

      elseif ( inobas .and. itype .eq. 3 ) then
        jp2nt = jend
        jp2tn = jp2nt + 4*nel*norbs
        jend  = jp2tn + 4*nel*norbs
      endif

      if ( inobas .and. itype .eq. 1 ) then
        jnos = jend
        jend = jnos + 2*nosq

      elseif ( inobas .and. itype .eq. 2 ) then
        jnos = jend
        jend = jnos + 4*nosq

      elseif ( inobas .and. itype .eq. 3 .or.
     $         ispuhf .or. ispghf ) then
        jnos = jend
        jend = jnos + 8*nosq
      endif

C       Allocate space for:
C         - transpose(nos) . nos product (ju)

      if ( inobas .and. itype .eq. 1 ) then
        ju   = jend
        jend = ju   + 2*nosq

      elseif ( inobas .and. itype .eq. 2 ) then
        ju   = jend
        jend = ju   + 4*nosq

      elseif ( inobas .and. itype .eq. 3 ) then
        ju   = jend
        jend = ju   + 8*nosq
      endif


C     Do some preparation for parallel processing:
C       - lenp,  dimension of scratch space per processor
C       - npuse, number of processors to use

      lscrdg = 0

      if ( ispuhf .or. ispghf ) then
        lscrdg = 32*nosq

        lenp = 2*szp + lscrdg + 8*nel*norbs + 16*nosq

      elseif ( itype .eq. 1 ) then
        lenp = 2*szp + 4*na*norbs

      elseif ( itype .eq. 2 ) then
        lenp = 2*szp + 4*na*norbs + 4*nb*norbs

      elseif ( itype .eq. 3 ) then
        lenp = 2*szp + 8*nel*norbs
      endif

C       .. length of Ng buffer, Yg buffer ..

      lenp = lenp + szbfn + szbfp

      if ( inobas ) then
        if ( itype .lt. 3 ) lenp = lenp + 2*nosq
        if ( itype .eq. 3 ) lenp = lenp + 8*nosq
      endif


      npmem = (lenv - jend + 1) / lenp

      if ( npmem .lt. 1 ) then
        call gauerr ('Not enough memory in formyg for a single CPU.')
      endif

      npmax = nproc (0)
      npuse = min (npmem, npmax)
      
      ipar  = lappar (0)
      throk = ipar .eq. 1

      if ( ipar .eq. -1 ) npuse = 1

      call tstcor (jend + lenp*npuse - 1, lenv, 'formyg')


C     Load transformation matrix (orthonormal AO => NO basis).

      if ( inobas .and. itype .eq. 1 ) then
        call fileio (2, -irwnos, 2*nosq, v(jnos), 0)

      elseif ( inobas .and. itype .eq. 2 ) then
        call fileio (2, -irwnos, 4*nosq, v(jnos), 0)

      elseif ( inobas .and. itype .eq. 3 .or.
     $         ispuhf .or. ispghf ) then
        call fileio (2, -irwnos, 8*nosq, v(jnos), 0)
      endif


C     Extract rectangular blocks of density matrix.

      if ( itype .eq. 1 ) then
        call sptmat (2, norbs, norbs, pno, na, norbs, v(jpnt), 1)
        call sptmat (2, norbs, norbs, pno, norbs, na, v(jptn), 1)

      elseif ( itype .eq. 2 ) then
        call sptmat (2, norbs, norbs, pno, na, norbs, v(jpant), 1)
        call sptmat (2, norbs, norbs, pno, norbs, na, v(jpatn), 1)

        if ( nb .ge. 1 ) then
        call sptmat (2, norbs, norbs, pno(1+2*nosq), nb, norbs,
     $       v(jpbnt), 1)
        call sptmat (2, norbs, norbs, pno(1+2*nosq), norbs, nb,
     $       v(jpbtn), 1)
        endif

      elseif ( itype .eq. 3 ) then
        call sptmat (2, 2*norbs, 2*norbs, pno, nel, 2*norbs, v(jpnt), 1)
        call sptmat (2, 2*norbs, 2*norbs, pno, 2*norbs, nel, v(jptn), 1)
c        call ltoutx (iout, 1, 2, 1, 0, 'jp', 2*norbs, 2*norbs,
c     $                  pno, pno, 0)
c        call ltoutx (iout, 3, 2, 1, 0, 'jpnt', nel, 2*norbs,
c     $                  v(jpnt), v(jpnt), 0)
c        call ltoutx (iout, 3, 2, 1, 0, 'jptn', 2*norbs, nel,
c     $                  v(jptn), v(jptn), 0)
      endif

      if ( inobas .and. itype .eq. 1 ) then
        call sptmat (2, norbs, norbs, pstrno, na, norbs, v(jp2nt), 1)
        call sptmat (2, norbs, norbs, pstrno, norbs, na, v(jp2tn), 1)

      elseif ( inobas .and. itype .eq. 2 ) then
        call sptmat (2, norbs, norbs, pstrno, na, norbs, v(jp2ant), 1)
        call sptmat (2, norbs, norbs, pstrno, norbs, na, v(jp2atn), 1)

        if ( nb .ge. 1 ) then
        call sptmat (2, norbs, norbs, pstrno(1+2*nosq), nb, norbs,
     $       v(jp2bnt), 1)
        call sptmat (2, norbs, norbs, pstrno(1+2*nosq), norbs, nb,
     $       v(jp2btn), 1)
        endif

      elseif ( inobas .and. itype .eq. 3 ) then
        call sptmat (2, 2*norbs, 2*norbs, pstrno, nel, 2*norbs,
     $       v(jp2nt), 1)
        call sptmat (2, 2*norbs, 2*norbs, pstrno, 2*norbs, nel,
     $       v(jp2tn), 1)
      endif


C     Compute product
C       U = transpose (nos) . nos

      if ( inobas .and. itype .eq. 1 ) then
        call zgemm ('t', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jnos), norbs, zero, v(ju), norbs)

      elseif ( inobas .and. itype .eq. 2 ) then
        call zgemm ('t', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jnos), norbs, zero, v(ju), norbs)
        call zgemm ('t', 'n', norbs, norbs, norbs, one, v(jnos+2*nosq),
     $       norbs, v(jnos+2*nosq), norbs, zero, v(ju+2*nosq), norbs)

      elseif ( inobas .and. itype .eq. 3 ) then
        call zgemm ('t', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jnos),
     $       2*norbs, v(jnos), 2*norbs, zero, v(ju), 2*norbs)
      endif



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Construct the matrices X(g)  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Construct the matrices X(g) as described above...
C     Additionally, we will compute the integrated quantity

C       int dg y(g) X(g),

C     stored in the array v(jy), and useful for computing Y(g).


C     Clear y matrices for all processors.

      do 11 ip = 1, npuse
        iy = jend + (ip - 1) * lenp

        call aclear (szp, v(iy))
 11   continue


C     Loop over grid points.

      do 103 iblock = 1, ngrdt, npuse
        nthrd = min (npuse, ngrdt - iblock + 1)


C     Load matrices. Loop not parallelized.

        do 104 ip = 1, nthrd

          igrdt = iblock + ip - 1

C     Memory allocation per processor:
C       iy    - integrated matrix:  int dg y(g) X(g)
C       ibufn - N(g) buffer at each grid point
C       ibufp - Y(g) buffer at each grid point

          iy    = jend  + (ip - 1) * lenp
          ibufn = iy    + szp
          ibufp = ibufn + szbfn
          iend  = ibufp + szbfp

          ind = (igrdt - 1)*szbfn

          call fileio (2, -irwng, szbfn, v(ibufn), ind)
 104    continue


C     Computation loop. Parallelized.

        throk = throk .and. npuse .gt. 1

C$omp   parallel do if (throk) schedule(static,1) default(shared)
C$omp+  private (igrdt, ibufn, ibufp, ing, inga, ingb)
C$omp+  private (iy, iyg, ix, iend)
C$omp+  private (it1, it1a, it1b, it2, it2a, it2b)
C$omp+  private (iz1f, iz2f, iscrdg, ix2)
C$omp+  private (igrda, igrdb, igrdg)
C$omp+  private (anga, angb, angg, facw, facd, fac)
C$omp+  private (n1c, n2c, n1s, n2s, indy, n1, n2)

        do 105 ip = 1, nthrd

          igrdt = iblock + ip - 1

C     Memory allocation per processor:
C       iy    - integrated matrix:  int dg y(g) X(g)
C       ibufn - N(g) buffer at each grid point
C       ibufp - Y(g) buffer at each grid point
C       ix    - scratch array used to build X(g) (ix)

          iy    = jend  + (ip - 1) * lenp
          ibufn = iy    + szp
          ibufp = ibufn + szbfn
          ix    = ibufp + szbfp
          iend  = ix    + szp

C       it1, it2, iz1f, iz2f - scratch arrays

          if ( itype .eq. 1 ) then
            it1  = iend
            it2  = it1  + 2*na*norbs
            iend = it2  + 2*na*norbs

          elseif ( itype .eq. 2 ) then
            it1a = iend
            it1b = it1a + 2*na*norbs
            it2a = it1b + 2*nb*norbs
            it2b = it2a + 2*na*norbs
            iend = it2b + 2*nb*norbs

          elseif ( itype .eq. 3 ) then
            it1  = iend
            it2  = it1  + 4*nel*norbs
            iend = it2  + 4*nel*norbs
          endif

          if ( ispuhf .or. ispghf ) then
            iz1f = iend
            iz2f = iz1f + 8*nosq
            iend = iz2f + 8*nosq
          endif

C       ix2 - scratch array to fix basis if inobas = .true.

          if ( inobas .and. itype .lt. 3 ) then
            ix2  = iend
            iend = ix2 + 2*nosq

          elseif ( inobas .and. itype .eq. 3 ) then
            ix2  = iend
            iend = ix2 + 8*nosq
          endif

C       iscrdg - scratch array for multiplication of Dg matrices

          if ( ispuhf .or. ispghf ) then
            iscrdg = iend
            iend   = iscrdg + lscrdg
          endif

C     Retrieve some useful factors.

          igrda = indarr(igrdt,1)
          igrdb = indarr(igrdt,2)
          igrdg = indarr(igrdt,3)

          anga = grda(igrda)
          angb = grdb(igrdb)
          angg = grdg(igrdg)


C     Compute useful factor
C       facw = wgt(g)

          facw = wgta(igrda) * wgtb(igrdb) * wgtg(igrdg)


C     Transform Dg matrix to NO basis:
C       nos' . Dg . nos

          if ( ispuhf .or. ispghf ) then
            call amove (8*nosq, v(jnos), v(iz1f))

            call mltdg (iout, iprint, v(iscrdg), lscrdg, norbs,
     $           v(iz1f), anga, angb, angg, 1)

            call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $           v(jnos), 2*norbs, v(iz1f), 2*norbs, zero, v(iz2f),
     $           2*norbs)
          endif
c          call ltoutx (iout, 1, 2, 1, 0, 'D(g) (NO)', 2*norbs,
c     $                 2*norbs, v(iz2f), v(iz2f), 0)


C     Loop over CI matrix elements (complex conjugation).

          do 101 n1c = 1, ncik
            do 102 n2c = 1, ncik

              ing = (n1c-1)*ncik*szn + (n2c-1)*szn + ibufn
              iyg = (n1c-1)*ncik*szp + (n2c-1)*szp + ibufp
              indy = (igrdt-1)*ncik*ncik + (n1c-1)*ncik + n2c


C     Useful indices for itype = 2.

              if ( itype .eq. 2 ) then
                inga = ing
                ingb = inga + 2*na*na
              endif


C     Clear yg matrices.

              call aclear (szp, v(iyg))


C     Compute the product:

C       ( T1_oo )  =  ( Dg_oo   Dg_ov ) . ( P_oo )
C       ( T1_vo )     ( Dg_vo   Dg_vv )   ( P_vo )

C     or

C       ( T1_oo )  =  ( Dg_oo   Dg_ov ) . ( P_oo   P_ov )^T
C       ( T1_vo )     ( Dg_vo   Dg_vv )

              if ( ispuhf .or. ispghf ) then
                if ( n2c .eq. 1 ) then
                  call zgemm ('n', 'n', 2*norbs, nel, 2*norbs, one,
     $                 v(iz2f), 2*norbs, v(jptn), 2*norbs, zero, v(it1),
     $                 2*norbs)

                elseif ( n2c .eq. 2 .and. inobas ) then
                  call zgemm ('n', 'n', 2*norbs, nel, 2*norbs, one,
     $                 v(iz2f), 2*norbs, v(jp2tn), 2*norbs, zero,
     $                 v(it1), 2*norbs)

                elseif ( n2c .eq. 2 .and. .not. inobas ) then
                  call zgemm ('n', 't', 2*norbs, nel, 2*norbs, one,
     $                 v(iz2f), 2*norbs, v(jpnt), nel, zero, v(it1),
     $                 2*norbs)
                endif

              else
                if ( n2c .eq. 1 ) then
                  if ( itype .eq. 1 ) then
                    call amove (2*na*norbs, v(jptn), v(it1))

                  elseif ( itype .eq. 2 ) then
                    call amove (2*na*norbs, v(jpatn), v(it1a))

                    if ( nb .ge. 1 ) then
                    call amove (2*nb*norbs, v(jpbtn), v(it1b))
                    endif

                  elseif ( itype .eq. 3 ) then
                    call amove (4*nel*norbs, v(jptn), v(it1))
                  endif

                elseif ( n2c .eq. 2 .and. inobas ) then
                  if ( itype .eq. 1 ) then
                    call amove (2*na*norbs, v(jp2tn), v(it1))

                  elseif ( itype .eq. 2 ) then
                    call amove (2*na*norbs, v(jp2atn), v(it1a))

                    if ( nb .ge. 1 ) then
                    call amove (2*nb*norbs, v(jp2btn), v(it1b))
                    endif

                  elseif ( itype .eq. 3 ) then
                    call amove (4*nel*norbs, v(jp2tn), v(it1))
                  endif

                elseif ( n2c .eq. 2 .and. .not. inobas ) then
                  if ( itype .eq. 1 ) then
                    call amove (2*na*norbs, v(jpnt), v(it2))

                    call trprct (2, na, norbs, v(it2), v(it1))

                  elseif ( itype .eq. 2 ) then
                    call amove (2*na*norbs, v(jpant), v(it2a))

                    if ( nb .ge. 1 ) then
                    call amove (2*nb*norbs, v(jpbnt), v(it2b))
                    endif

                    call trprct (2, na, norbs, v(it2a), v(it1a))

                    if ( nb .ge. 1 ) then
                    call trprct (2, nb, norbs, v(it2b), v(it1b))
                    endif

                  elseif ( itype .eq. 3 ) then
                    call amove (4*nel*norbs, v(jpnt), v(it2))

                    call trprct (2, nel, 2*norbs, v(it2), v(it1))
                  endif
                endif
              endif


C     Compute the product:

C       ( T2_oo )  =  ( T1_oo ) . Ng
C       ( T2_vo )     ( T1_vo )

              if ( itype .eq. 1 ) then
                call zgemm ('n', 'n', norbs, na, na, one, v(it1), norbs,
     $               v(ing), na, zero, v(it2), norbs)

              elseif ( itype .eq. 2 ) then
                call zgemm ('n', 'n', norbs, na, na, one, v(it1a),
     $               norbs, v(inga), na, zero, v(it2a), norbs)

                if ( nb .ge. 1 ) then
                call zgemm ('n', 'n', norbs, nb, nb, one, v(it1b),
     $               norbs, v(ingb), nb, zero, v(it2b), norbs)
                endif

              elseif ( itype .eq. 3 ) then
                call zgemm ('n', 'n', 2*norbs, nel, nel, one, v(it1),
     $               2*norbs, v(ing), nel, zero, v(it2), 2*norbs)
c                call ltoutx (iout, 3, 2, 1, 0, 'X(g)2', 2*norbs,
c     $                 nel, v(it2), v(it2), 0)
              endif


C     Add T2 to X(g). Note that for n1c = 2 we need to take the transpose.

              if ( n1c .eq. 1 ) then
                if ( itype .eq. 1 ) then
                  call sptmat (2, norbs, norbs, v(ix), norbs, na,
     $                 v(it2), 2)

                elseif ( itype .eq. 2 ) then
                  call sptmat (2, norbs, norbs, v(ix), norbs, na,
     $                 v(it2a), 2)

                  if ( nb .ge. 1 ) then
                  call sptmat (2, norbs, norbs, v(ix+2*nosq), norbs, nb,
     $                 v(it2b), 2)
                  else
                    call aclear (2*nosq, v(ix+2*nosq))
                  endif

                elseif ( itype .eq. 3 ) then
                  call sptmat (2, 2*norbs, 2*norbs, v(ix), 2*norbs, nel,
     $                 v(it2), 2)
                endif

              elseif ( n1c .eq. 2 ) then
                if ( itype .eq. 1 ) then
                  call trprct (2, norbs, na, v(it2), v(it1))

                  call sptmat (2, norbs, norbs, v(ix), na, norbs,
     $                 v(it1), 2)

                elseif ( itype .eq. 2 ) then
                  call trprct (2, norbs, na, v(it2a), v(it1a))

                  if ( nb .ge. 1 ) then
                  call trprct (2, norbs, nb, v(it2b), v(it1b))
                  endif

                  call sptmat (2, norbs, norbs, v(ix), na, norbs,
     $                 v(it1a), 2)

                  if ( nb .ge. 1 ) then
                  call sptmat (2, norbs, norbs, v(ix+2*nosq), nb, norbs,
     $                 v(it1b), 2)
                  else
                    call aclear (2*nosq, v(ix+2*nosq))
                  endif

                elseif ( itype .eq. 3 ) then
                  call trprct (2, 2*norbs, nel, v(it2), v(it1))

                  call sptmat (2, 2*norbs, 2*norbs, v(ix), nel, 2*norbs,
     $                 v(it1), 2)
                endif
              endif

C              .. fix basis for inobas = .true. ..

              if ( n1c .eq. 2 .and. inobas ) then
                if ( itype .eq. 1 ) then
                  call zgemm ('c', 'n', norbs, norbs, norbs, one,
     $                 v(ju), norbs, v(ix), norbs, zero, v(ix2), norbs)
                  call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $                 v(ix2), norbs, v(ju), norbs, zero, v(ix), norbs)

                elseif ( itype .eq. 2 ) then
                  call zgemm ('c', 'n', norbs, norbs, norbs, one,
     $                 v(ju), norbs, v(ix), norbs, zero, v(ix2), norbs)
                  call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $                 v(ix2), norbs, v(ju), norbs, zero, v(ix), norbs)

                  call zgemm ('c', 'n', norbs, norbs, norbs, one,
     $                 v(ju+2*nosq), norbs, v(ix+2*nosq), norbs, zero,
     $                 v(ix2), norbs)
                  call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $                 v(ix2), norbs, v(ju+2*nosq), norbs, zero,
     $                 v(ix+2*nosq), norbs)

                elseif ( itype .eq. 3 ) then
                  call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $                 v(ju), 2*norbs, v(ix), 2*norbs, zero, v(ix2),
     $                 2*norbs)
                  call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $                 v(ix2), 2*norbs, v(ju), 2*norbs, zero, v(ix),
     $                 2*norbs)
                endif
              endif

              call aadd (szp, v(iyg), v(ix), v(iyg))



C     Compute the product:

C       ( T1_oo   T1_ov )  =  ( P_oo   P_ov ) . ( Dg_oo   Dg_ov )
C                                               ( Dg_vo   Dg_vv )

C     or

C       ( T1_oo   T1_ov )  =  ( P_oo )^T . ( Dg_oo   Dg_ov )
C                             ( P_vo )     ( Dg_vo   Dg_vv )

              if ( ispuhf .or. ispghf ) then
                if ( n1c .eq. 1 ) then
                  call zgemm ('n', 'n', nel, 2*norbs, 2*norbs, one,
     $                 v(jpnt), nel, v(iz2f), 2*norbs, zero, v(it1),
     $                 nel)

                elseif ( n1c .eq. 2 .and. inobas ) then
                  call zgemm ('n', 'n', nel, 2*norbs, 2*norbs, one,
     $                 v(jp2nt), nel, v(iz2f), 2*norbs, zero, v(it1),
     $                 nel)

                elseif ( n1c .eq. 2 .and. .not. inobas ) then
                  call zgemm ('t', 'n', nel, 2*norbs, 2*norbs, one,
     $                 v(jptn), 2*norbs, v(iz2f), 2*norbs, zero, v(it1),
     $                 nel)
                endif

              else
                if ( n1c .eq. 1 ) then
                  if ( itype .eq. 1 ) then
                    call amove (2*na*norbs, v(jpnt), v(it1))

                  elseif ( itype .eq. 2 ) then
                    call amove (2*na*norbs, v(jpant), v(it1a))

                    if ( nb .ge. 1 ) then
                    call amove (2*nb*norbs, v(jpbnt), v(it1b))
                    endif

                  elseif ( itype .eq. 3 ) then
                    call amove (4*nel*norbs, v(jpnt), v(it1))
                  endif

                elseif ( n1c .eq. 2 .and. inobas ) then
                  if ( itype .eq. 1 ) then
                    call amove (2*na*norbs, v(jp2nt), v(it1))

                  elseif ( itype .eq. 2 ) then
                    call amove (2*na*norbs, v(jp2ant), v(it1a))

                    if ( nb .ge. 1 ) then
                    call amove (2*nb*norbs, v(jp2bnt), v(it1b))
                    endif

                  elseif ( itype .eq. 3 ) then
                    call amove (4*nel*norbs, v(jp2nt), v(it1))
                  endif

                elseif ( n1c .eq. 2 .and. .not. inobas ) then
                  if ( itype .eq. 1 ) then
                    call amove (2*na*norbs, v(jptn), v(it2))

                    call trprct (2, norbs, na, v(it2), v(it1))

                  elseif ( itype .eq. 2 ) then
                    call amove (2*na*norbs, v(jpatn), v(it2a))

                    if ( nb .ge. 1 ) then
                    call amove (2*nb*norbs, v(jpbtn), v(it2b))
                    endif

                    call trprct (2, norbs, na, v(it2a), v(it1a))

                    if ( nb .ge. 1 ) then
                    call trprct (2, norbs, nb, v(it2b), v(it1b))
                    endif

                  elseif ( itype .eq. 3 ) then
                    call amove (4*nel*norbs, v(jptn), v(it2))

                    call trprct (2, 2*norbs, nel, v(it2), v(it1))
                  endif
                endif
              endif


C     Compute the product:

C       ( T2_oo   T2_ov )  =  Ng . ( T1_oo   T1_ov )

              if ( itype .eq. 1 ) then
                call zgemm ('n', 'n', na, norbs, na, one, v(ing), na,
     $               v(it1), na, zero, v(it2), na)

              elseif ( itype .eq. 2 ) then
                call zgemm ('n', 'n', na, norbs, na, one, v(inga), na,
     $               v(it1a), na, zero, v(it2a), na)

                if ( nb .ge. 1 ) then
                call zgemm ('n', 'n', nb, norbs, nb, one, v(ingb), nb,
     $               v(it1b), nb, zero, v(it2b), nb)
                endif

              elseif ( itype .eq. 3 ) then
                call zgemm ('n', 'n', nel, 2*norbs, nel, one, v(ing),
     $               nel, v(it1), nel, zero, v(it2), nel)
c                call ltoutx (iout, 3, 2, 1, 0, 'X(g)2', nel, 2*norbs,
c     $                 v(it2), v(it2), 0)
              endif


C     Add T2 to X(g). Note that for n2c = 2 we need to take the transpose.

              if ( n2c .eq. 1 ) then
                if ( itype .eq. 1 ) then
                  call sptmat (2, norbs, norbs, v(ix), na, norbs,
     $                 v(it2), 2)

                elseif ( itype .eq. 2 ) then
                  call sptmat (2, norbs, norbs, v(ix), na, norbs,
     $                 v(it2a), 2)

                  if ( nb .ge. 1 ) then
                  call sptmat (2, norbs, norbs, v(ix+2*nosq), nb, norbs,
     $                 v(it2b), 2)
                  else
                    call aclear (2*nosq, v(ix+2*nosq))
                  endif

                elseif ( itype .eq. 3 ) then
                  call sptmat (2, 2*norbs, 2*norbs, v(ix), nel, 2*norbs,
     $                 v(it2), 2)
                endif

              elseif ( n2c .eq. 2 ) then
                if ( itype .eq. 1 ) then
                  call trprct (2, na, norbs, v(it2), v(it1))

                  call sptmat (2, norbs, norbs, v(ix), norbs, na,
     $                 v(it1), 2)

                elseif ( itype .eq. 2 ) then
                  call trprct (2, na, norbs, v(it2a), v(it1a))

                  if ( nb .ge. 1 ) then
                  call trprct (2, nb, norbs, v(it2b), v(it1b))
                  endif

                  call sptmat (2, norbs, norbs, v(ix), norbs, na,
     $                 v(it1a), 2)

                  if ( nb .ge. 1 ) then
                  call sptmat (2, norbs, norbs, v(ix+2*nosq), norbs, nb,
     $                 v(it1b), 2)
                  else
                    call aclear (2*nosq, v(ix+2*nosq))
                  endif

                elseif ( itype .eq. 3 ) then
                  call trprct (2, nel, 2*norbs, v(it2), v(it1))

                  call sptmat (2, 2*norbs, 2*norbs, v(ix), 2*norbs, nel,
     $                 v(it1), 2)
                endif
              endif

C              .. fix basis for inobas = .true. ..

              if ( n2c .eq. 2 .and. inobas ) then
                if ( itype .eq. 1 ) then
                  call zgemm ('c', 'n', norbs, norbs, norbs, one,
     $                 v(ju), norbs, v(ix), norbs, zero, v(ix2), norbs)
                  call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $                 v(ix2), norbs, v(ju), norbs, zero, v(ix), norbs)

                elseif ( itype .eq. 2 ) then
                  call zgemm ('c', 'n', norbs, norbs, norbs, one,
     $                 v(ju), norbs, v(ix), norbs, zero, v(ix2), norbs)
                  call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $                 v(ix2), norbs, v(ju), norbs, zero, v(ix), norbs)

                  call zgemm ('c', 'n', norbs, norbs, norbs, one,
     $                 v(ju+2*nosq), norbs, v(ix+2*nosq), norbs, zero,
     $                 v(ix2), norbs)
                  call zgemm ('n', 'n', norbs, norbs, norbs, one,
     $                 v(ix2), norbs, v(ju+2*nosq), norbs, zero,
     $                 v(ix+2*nosq), norbs)

                elseif ( itype .eq. 3 ) then
                  call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $                 v(ju), 2*norbs, v(ix), 2*norbs, zero, v(ix2),
     $                 2*norbs)
                  call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one,
     $                 v(ix2), 2*norbs, v(ju), 2*norbs, zero, v(ix),
     $                 2*norbs)
                endif
              endif

              call aadd (szp, v(iyg), v(ix), v(iyg))


C     Add contribution of current grid point to the integrated quantity.

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
               write(iout, 1166) indy, yg(indy), facd , civec(n1), facw
1166  Format('indy = ', I5, ', yg = ', F10.5, F10.5 , ', facd = ', 
     $  F6.3, F6.3 , ', civec = ', F10.5, F10.5, ', facw = ', F10.5
     $ )

                  call amove (szp, v(iyg), v(ix))
                  call zscal (szp/2, facw*fac, v(ix), 1)

                  call aadd (szp, v(iy), v(ix), v(iy))
 202            continue
 201          continue
 102        continue
 101      continue

 105    continue

C$omp   end parallel do


C     Save matrices. Loop not parallelized.

        do 106 ip = 1, nthrd

          igrdt = iblock + ip - 1

C     Memory allocation per processor:
C       iy    - integrated matrix:  int dg y(g) X(g)
C       ibufn - N(g) buffer at each grid point
C       ibufp - Y(g) buffer at each grid point

          iy    = jend  + (ip - 1) * lenp
          ibufn = iy    + szp
          ibufp = ibufn + szbfn
          iend  = ibufp + szbfp

          ind = (igrdt - 1)*szbfp

C     Save X(g) matrices into RWF file.

          call fileio (1, -irwyg, szbfp, v(ibufp), ind)


C     Print matrices if requested.

          do 301 n1 = 1, ncik
            do 302 n2 = 1, ncik

              iyg = (n1-1)*ncik*szp + (n2-1)*szp + ibufp

              if ( iprint .gt. 4 ) then
                write (indc1, '(I5)') igrdt
                write (indc2, '(I2)') n1
                write (indc3, '(I2)') n2
                str = ' (' // indc1 // ',' // indc2 //
     $                ',' // indc3 // ')'

                write (iout, *) ' '
                if ( itype .eq. 1 ) then
                  call ltoutx (iout, 1, 2, 1, 0, 'X(g)' // str, norbs,
     $                 norbs, v(iyg), v(iyg), 0)

                elseif ( itype .eq. 2 ) then
                  call ltoutx (iout, 1, 3, 1, 0, 'X(g)' // str, norbs,
     $                 norbs, v(iyg), v(iyg+2*nosq), 0)

                elseif ( itype .eq. 3 ) then
                  call ltoutx (iout, 1, 2, 1, 0, 'X(g)' // str, 2*norbs,
     $                 2*norbs, v(iyg), v(iyg), 0)
                endif
              endif

 302        continue
 301      continue

 106    continue
 103  continue


C     Merge results from all processors.

      call aclear (szp, v(jy))

      do 21 ip = 1, npuse
        iy = jend + (ip - 1) * lenp

        call aadd (szp, v(jy), v(iy), v(jy))
 21   continue


C     Save integrated X(g) matrix.

      call amove (szp, v(jy), intxg)
      call ltoutx (iout, 1, 2, 1, 0, 'X(g) int' , 2*norbs,
     $                 2*norbs, v(jy), v(jy), 0)


      return
      end


