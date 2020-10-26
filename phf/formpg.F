      

      subroutine formpg (iout, iprint, v, lenv, iphf, ncik, norbs,
     $     na, nb, ngrdt, indarr, grda, grdb, grdg, irwnos, irwng,
     $     irwpg, pno, pstrno, inobas)

      implicit none


C +------------------------------------------------------------+
C |  formpg   --   CAJH, 07.2011                               |
C |                                                            |
C |                                                            |
C |  Build the transition density matrices, P(g), which are    |
C |  central to our projected HF scheme.                       |
C |                                                            |
C |  Given a rotation matrix Dg, related to the spin           |
C |  projection operator with an angle Omega, the matrices     |
C |  P(g) can take the form                                    |
C |                                                            |
C |  if nu = I, mu = I  ( <0| I! R(Omega) I |0> )              |
C |                                                            |
C |    P(g)  =  Dg . ( P_oo ) . N(g) . ( P_oo  P_ov )          |
C |                  ( P_vo )                                  |
C |                                                            |
C |  if nu = K, mu = I  ( <0| I! R(Omega) K |0> )              |
C |                                                            |
C |    P(g)  =  Dg . ( P_oo  P_ov )^T . N(g) . ( P_oo  P_ov )  |
C |                                                            |
C |  if nu = I, mu = K  ( <0| K! R(Omega) I |0> )              |
C |                                                            |
C |    P(g)  =  Dg . ( P_oo ) . N(g) . ( P_oo )^T              |
C |                  ( P_vo )          ( P_vo )                |
C |                                                            |
C |  if nu = K, mu = K  ( <0| K! R(Omega) K |0> )              |
C |                                                            |
C |    P(g)  =  Dg . ( P_oo  P_ov )^T . N(g) . ( P_oo )^T      |
C |                                            ( P_vo )        |
C |                                                            |
C |  where nu, mu = I, K are operators related to complex      |
C |  conjugation restoration. (See evalsmt for more details.)  |
C |  The matrices N(g) are those built in formng.              |
C |                                                            |
C |  Other than that, it is worth recalling that P is the      |
C |  density matrix of the unprojected HF state. The indices   |
C |  o and v define the dimensions (occupied, virtual) of the  |
C |  rectangular matrices. The Dg matrices are always          |
C |  of dimension nbasis x nbasis, whereas the N(g) matrices   |
C |  are of dimension N x N, where N is the number of          |
C |  electrons in the system.                                  |
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
C |  Reworked subroutine to avoid looping over CI matrix       |
C |  elements from spin.                                       |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       ncik   - number of basis states for complex conj projection
C       norbs  - number of orbitals
C       na     - number of alpha electrons
C       nb     - number of beta electrons
C       ngrdt  - total number of grid points
C       indarr - mapping array igrdt => igrdX
C       grdX   - integration grid for each variable
C       irwnos - RWF file with transformation matrix from orthonormal AO
C                to NO basis
C       irwng  - RWF file with matrices N(g)
C       irwpg  - RWF file with matrices P(g) [ updated here ]
C       pno    - density matrix in NO basis
C       pstrno - density matrix (P*) in NO basis
C       inobas - if .true., force the use of the true NO basis

      real*8      grda(*), grdb(*), grdg(*)
      integer     iout, iprint, lenv, iphf, ncik, na, nb
      integer     norbs, ngrdt, indarr(ngrdt,*)
      integer     irwnos, irwng, irwpg
      logical     inobas

      real*8      v(*), pno(*), pstrno(*)

C     other variables

      integer     nel, nosq, itype, szn, szp, szbfn, szbfp, ind
      integer     jnos, jpnt, jptn, jpant, jpatn, jpbnt, jpbtn, jend
      integer     jp2nt, jp2tn, jp2ant, jp2atn, jp2bnt, jp2btn
      integer     lscrdg
      integer     lenp, npmem, npmax, npuse, ipar
      integer     n1, n2, n1c, n2c, n1s, n2s, ip, iblock, nthrd
      logical     icmplx, ispuhf, ispghf, throk
      character   indc1*5, indc2*2, indc3*2, str*14

      integer     igrdt, ibufn, ibufp, ing, inga, ingb, ipg, iend
      integer     it1, it1a, it1b, it2, it2a, it2b
      integer     iz1f, iz2f, iscrdg
      integer     igrda, igrdb, igrdg
      real*8      anga, angb, angg

C     functions

      integer     nproc, lappar

C     constants

      complex*16  zero, one


      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)


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
        call gauerr ('Unrecognized iphf in formpg.')
      endif


C     Define useful array dimensions:
C       - szn,  size of N(g) matrices
C       - szp,  size of P(g) matrices

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

C       .. buffer to save computed Pg, load computed Ng ..

      szbfn = ncik*ncik * szn
      szbfp = ncik*ncik * szp


C     Memory allocation.

C       Allocate space for:
C         - orthonormal AO => NO basis transformation matrix (jnos)
C         - rectangular blocks of density matrix in NO basis (jpnt, jptn)
C         - rectangular blocks of density matrix (P*) in NO basis (jp2*)

      if ( itype .eq. 1 ) then
        jpnt = 1
        jptn = jpnt + 2*na*norbs
        jend = jptn + 2*na*norbs

      elseif ( itype .eq. 2 ) then
        jpant = 1
        jpatn = jpant + 2*na*norbs
        jpbnt = jpatn + 2*na*norbs
        jpbtn = jpbnt + 2*nb*norbs
        jend  = jpbtn + 2*nb*norbs

      elseif ( itype .eq. 3 ) then
        jpnt = 1
        jptn = jpnt + 4*nel*norbs
        jend = jptn + 4*nel*norbs
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

      if ( ispuhf .or. ispghf ) then
        jnos = jend
        jend = jnos + 8*nosq
      endif


C     Do some preparation for parallel processing:
C       - lenp,  dimension of scratch space per processor
C       - npuse, number of processors to use

      lscrdg = 0

      if ( ispuhf .or. ispghf ) then
        lscrdg = 32*nosq

        lenp = lscrdg + 8*nel*norbs + 16*nosq

      elseif ( itype .eq. 1 ) then
        lenp = 4*na*norbs

      elseif ( itype .eq. 2 ) then
        lenp = 4*na*norbs + 4*nb*norbs

      elseif ( itype .eq. 3 ) then
        lenp = 8*nel*norbs
      endif

C       .. length of Ng buffer, Pg buffer ..

      lenp = lenp + szbfn + szbfp

      npmem = (lenv - jend + 1) / lenp

      if ( npmem .lt. 1 ) then
        call gauerr ('Not enough memory in formpg for a single CPU.')
      endif

      npmax = nproc (0)
      npuse = min (npmem, npmax)
      
      ipar  = lappar (0)
      throk = ipar .eq. 1

      if ( ipar .eq. -1 ) npuse = 1

      call tstcor (jend + lenp*npuse - 1, lenv, 'formpg')


C     Load transformation matrix (orthonormal AO => NO basis).

      if ( ispuhf .or. ispghf ) then
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



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Construct the matrices P(g)  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Construct the matrices P(g) as described above...

C     Loop over grid points.

      do 103 iblock = 1, ngrdt, npuse
        nthrd = min (npuse, ngrdt - iblock + 1)


C     Load matrices. Loop not parallelized.

        do 104 ip = 1, nthrd

          igrdt = iblock + ip - 1

C     Memory allocation per processor:
C       ibufn - N(g) buffer at each grid point

          ibufn = jend  + (ip - 1) * lenp
          iend  = ibufn + szbfn

          ind = (igrdt - 1)*szbfn

          call fileio (2, -irwng, szbfn, v(ibufn), ind)
 104    continue


C     Computation loop. Parallelized.

        throk = throk .and. npuse .gt. 1

C$omp   parallel do if (throk) schedule(static,1) default(shared)
C$omp+  private (igrdt, ibufn, ibufp, ing, inga, ingb, ipg, iend)
C$omp+  private (it1, it1a, it1b, it2, it2a, it2b)
C$omp+  private (iz1f, iz2f, iscrdg)
C$omp+  private (igrda, igrdb, igrdg)
C$omp+  private (anga, angb, angg)
C$omp+  private (n1, n2)

        do 105 ip = 1, nthrd

          igrdt = iblock + ip - 1

C     Memory allocation per processor:
C       ibufn - N(g) buffer at each grid point
C       ibufp - P(g) buffer at each grid point

          ibufn = jend  + (ip - 1) * lenp
          ibufp = ibufn + szbfn
          iend  = ibufp + szbfp

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


C     Loop over CI matrix elements (complex conjugation).

          do 101 n1 = 1, ncik
            do 102 n2 = 1, ncik

              ing = (n1-1)*ncik*szn + (n2-1)*szn + ibufn
              ipg = (n1-1)*ncik*szp + (n2-1)*szp + ibufp

C     Useful indices for itype = 2.

              if ( itype .eq. 2 ) then
                inga = ing
                ingb = inga + 2*na*na
              endif


C     Compute the product:

C       ( T1_oo )  =  ( Dg_oo   Dg_ov ) . ( P_oo )
C       ( T1_vo )     ( Dg_vo   Dg_vv )   ( P_vo )

C     or

C       ( T1_oo )  =  ( Dg_oo   Dg_ov ) . ( P_oo   P_ov )^T
C       ( T1_vo )     ( Dg_vo   Dg_vv )

              if ( ispuhf .or. ispghf ) then
                if ( n2 .eq. 1 ) then
                  call zgemm ('n', 'n', 2*norbs, nel, 2*norbs, one,
     $                 v(iz2f), 2*norbs, v(jptn), 2*norbs, zero, v(it1),
     $                 2*norbs)

                elseif ( n2 .eq. 2 .and. inobas ) then
                  call zgemm ('n', 'n', 2*norbs, nel, 2*norbs, one,
     $                 v(iz2f), 2*norbs, v(jp2tn), 2*norbs, zero,
     $                 v(it1), 2*norbs)

                elseif ( n2 .eq. 2 .and. .not. inobas ) then
                  call zgemm ('n', 't', 2*norbs, nel, 2*norbs, one,
     $                 v(iz2f), 2*norbs, v(jpnt), nel, zero, v(it1),
     $                 2*norbs)
                endif

              else
                if ( n2 .eq. 1 ) then
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

                elseif ( n2 .eq. 2 .and. inobas ) then
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

                elseif ( n2 .eq. 2 .and. .not. inobas ) then
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

C       ( T2_oo )  =  ( T1_oo ) . N(g)
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
              endif


C     Compute the product:

C       ( P(g)_oo   P(g)_ov )  =  ( T2_oo ) . ( P_oo   P_ov )
C       ( P(g)_vo   P(g)_vv )     ( T2_vo )

C     or

C       ( P(g)_oo   P(g)_ov )  =  ( T2_oo ) . ( P_oo )^T
C       ( P(g)_vo   P(g)_vv )     ( T2_vo )   ( P_vo )

              if ( n1 .eq. 1 ) then
                if ( itype .eq. 1 ) then
                  call zgemm ('n', 'n', norbs, norbs, na, one, v(it2),
     $                 norbs, v(jpnt), na, zero, v(ipg), norbs)

                elseif ( itype .eq. 2 ) then
                  call zgemm ('n', 'n', norbs, norbs, na, one, v(it2a),
     $                 norbs, v(jpant), na, zero, v(ipg), norbs)

                  if ( nb .ge. 1 ) then
                  call zgemm ('n', 'n', norbs, norbs, nb, one, v(it2b),
     $                 norbs, v(jpbnt), nb, zero, v(ipg+2*nosq), norbs)
                  else
                    call aclear (2*nosq, v(ipg+2*nosq))
                  endif

                elseif ( itype .eq. 3 ) then
                  call zgemm ('n', 'n', 2*norbs, 2*norbs, nel, one,
     $                 v(it2), 2*norbs, v(jpnt), nel, zero, v(ipg),
     $                 2*norbs)
                endif

              elseif ( n1 .eq. 2 .and. inobas ) then
                if ( itype .eq. 1 ) then
                  call zgemm ('n', 'n', norbs, norbs, na, one, v(it2),
     $                 norbs, v(jp2nt), na, zero, v(ipg), norbs)

                elseif ( itype .eq. 2 ) then
                  call zgemm ('n', 'n', norbs, norbs, na, one, v(it2a),
     $                 norbs, v(jp2ant), na, zero, v(ipg), norbs)

                  if ( nb .ge. 1 ) then
                  call zgemm ('n', 'n', norbs, norbs, nb, one, v(it2b),
     $                 norbs, v(jp2bnt), nb, zero, v(ipg+2*nosq), norbs)
                  else
                    call aclear (2*nosq, v(ipg+2*nosq))
                  endif

                elseif ( itype .eq. 3 ) then
                  call zgemm ('n', 'n', 2*norbs, 2*norbs, nel, one,
     $                 v(it2), 2*norbs, v(jp2nt), nel, zero, v(ipg),
     $                 2*norbs)
                endif

              elseif ( n1 .eq. 2 .and. .not. inobas ) then
                if ( itype .eq. 1 ) then
                  call zgemm ('n', 't', norbs, norbs, na, one, v(it2),
     $                 norbs, v(jptn), norbs, zero, v(ipg), norbs)

                elseif ( itype .eq. 2 ) then
                  call zgemm ('n', 't', norbs, norbs, na, one, v(it2a),
     $                 norbs, v(jpatn), norbs, zero, v(ipg), norbs)

                  if ( nb .ge. 1 ) then
                  call zgemm ('n', 't', norbs, norbs, nb, one, v(it2b),
     $                 norbs, v(jpbtn), norbs, zero, v(ipg+2*nosq),
     $                 norbs)
                  else
                    call aclear (2*nosq, v(ipg+2*nosq))
                  endif

                elseif ( itype .eq. 3 ) then
                  call zgemm ('n', 't', 2*norbs, 2*norbs, nel, one,
     $                 v(it2), 2*norbs, v(jptn), 2*norbs, zero, v(ipg),
     $                 2*norbs)
                endif
              endif

 102        continue
 101      continue
 105    continue

C$omp   end parallel do


C     Save matrices. Loop not parallelized.

        do 106 ip = 1, nthrd

          igrdt = iblock + ip - 1

C     Memory allocation per processor:
C       ibufn - N(g) buffer at each grid point
C       ibufp - P(g) buffer at each grid point

          ibufn = jend  + (ip - 1) * lenp
          ibufp = ibufn + szbfn
          iend  = ibufp + szbfp

          ind = (igrdt - 1)*szbfp

C     Save P(g) matrices into RWF file.

          call fileio (1, -irwpg, szbfp, v(ibufp), ind)


C     Print matrices if requested.

          do 201 n1 = 1, ncik
            do 202 n2 = 1, ncik

              ipg = (n1-1)*ncik*szp + (n2-1)*szp + ibufp

              if ( iprint .gt. 4 ) then
                write (indc1, '(I5)') igrdt
                write (indc2, '(I2)') n1
                write (indc3, '(I2)') n2
                str = ' (' // indc1 // ',' // indc2 //
     $                ',' // indc3 // ')'

                write (iout, *) ' '
                if ( itype .eq. 1 ) then
                  call ltoutx (iout, 1, 2, 1, 0, 'P(g)' // str, norbs,
     $                 norbs, v(ipg), v(ipg), 0)

                elseif ( itype .eq. 2 ) then
                  call ltoutx (iout, 1, 3, 1, 0, 'P(g)' // str, norbs,
     $                 norbs, v(ipg), v(ipg+2*nosq), 0)

                elseif ( itype .eq. 3 ) then
                  call ltoutx (iout, 1, 2, 1, 0, 'P(g)' // str, 2*norbs,
     $                 2*norbs, v(ipg), v(ipg), 0)
                endif
              endif

 202        continue
 201      continue

 106    continue
 103  continue


      return
      end


