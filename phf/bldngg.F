

      subroutine bldngg (iout, iprint, itype, ispuhf, ispghf, n1c, n2c,
     $     na, nb, norbs, inobas, pnt, ptn, pbnt, pbtn, p2nt, p2tn,
     $     p2bnt, p2btn, scr1, scr1b, scr2, scr2b, dg1, piv, ngg)

      implicit none


C +------------------------------------------------------------+
C |  bldngg   --   CAJH, 08.2011                               |
C |                                                            |
C |                                                            |
C |  Build the overlap related matrices, N(g,g'), in a         |
C |  similar way as formng builds the matrices N(g).           |
C |                                                            |
C |  Given a rotation matrix Dg, related to the spin           |
C |  projection operator with an angle Omega, the matrices     |
C |  N(g,g') can take the form:                                |
C |                                                            |
C |  if nu = I, mu = I  ( <0| I! R(Omega')! R(Omega) I |0> )   |
C |                                                            |
C |    [N(g,g')]^{-1}  =                                       |
C |       ( P_oo  P_ov ) . Dg'! . Dg . ( P_oo )                |
C |                                    ( P_vo )                |
C |                                                            |
C |  if nu = K, mu = I  ( <0| I! R(Omega')! R(Omega) K |0> )   |
C |                                                            |
C |    [N(g,g')]^{-1}  =                                       |
C |       ( P_oo  P_ov ) . Dg'! . Dg . ( P_oo  P_ov )^T        |
C |                                                            |
C |  if nu = I, mu = K  ( <0| K! R(Omega')! R(Omega) I |0> )   |
C |                                                            |
C |    [N(g,g')]^{-1}  =                                       |
C |       ( P_oo )^T  . Dg'! . Dg . ( P_oo )                   |
C |       ( P_vo )                  ( P_vo )                   |
C |                                                            |
C |  if nu = K, mu = K  ( <0| K! R(Omega')! R(Omega) K |0> )   |
C |                                                            |
C |    [N(g,g')]^{-1}  =                                       |
C |       ( P_oo )^T . Dg'! . Dg . ( P_oo  P_ov )^T            |
C |       ( P_vo )                                             |
C |                                                            |
C |  where nu, mu = I, K are operators related to complex      |
C |  conjugation restoration. (See evalsmt for more details.)  |
C |                                                            |
C |  Here, P is the density matrix (in the NO basis) and T     |
C |  denotes the transpose of a matrix. The indices o and v    |
C |  define the dimensions (occupied, virtual) of the          |
C |  rectangular matrices. The Dg matrices are always          |
C |  nbasis x nbasis matrices.                                 |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  CAJH, 08.2011                                             |
C |                                                            |
C |  Major modification of the subroutine. It now computes     |
C |  the matrix N(g,g') for a single integration point g,g'.   |
C |  The loop over grid points is performed in phfdm1.F.       |
C |                                                            |
C |  Array assignment and dimensions                           |
C |                                                            |
C |    itype = 1,                                              |
C |                                                            |
C |      pnt  -  2*na*norbs     [ in ]                         |
C |      ptn  -  2*na*norbs     [ in ]                         |
C |      scr1 -  2*na*norbs                                    |
C |      scr2 -  2*na*norbs                                    |
C |      piv  -  intowp (na)                                   |
C |      ngg  -  2*na*na        [ out ]                        |
C |                                                            |
C |    itype = 2,                                              |
C |                                                            |
C |      pnt   -  2*na*norbs    [ in ]                         |
C |      ptn   -  2*na*norbs    [ in ]                         |
C |      pbnt  -  2*nb*norbs    [ in ]                         |
C |      pbtn  -  2*nb*norbs    [ in ]                         |
C |      scr1  -  2*na*norbs                                   |
C |      scr1b -  2*nb*norbs                                   |
C |      scr2  -  2*na*norbs                                   |
C |      scr2b -  2*nb*norbs                                   |
C |      piv   -  intowp (na)                                  |
C |      ngg   -  2*na*na + 2*nb*nb   [ out ]                  |
C |                                                            |
C |    itype = 3,                                              |
C |                                                            |
C |      pnt  -  2*nel*norbs    [ in ]                         |
C |      ptn  -  2*nel*norbs    [ in ]                         |
C |      scr1 -  2*nel*norbs                                   |
C |      scr2 -  2*nel*norbs                                   |
C |      piv  -  intowp (nel)                                  |
C |      ngg  -  2*nel*nel      [ out ]                        |
C |                                                            |
C |    if ( ispuhf .or. ispghf )                               |
C |                                                            |
C |      dg1 -  8*nosq                                         |
C |      dg2 -  8*nosq                                         |
C |                                                            |
C |  Here,                                                     |
C |                                                            |
C |    p?nt - N rows    of density matrix in NO basis          |
C |    p?tn - N columns of density matrix in NO basis          |
C |    ngg  - output matrix N(g,g')                            |
C |    dg1  - rotation matrix for angle g                      |
C |    dg2  - rotation matrix for angle g'                     |
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
C |  Moved to single integration. dg2 is no longer declared.   |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       itype  - type of PHF calculation (see phfdrv for details)
C       ispuhf - whether spin projection on UHF wavefunction is on
C       ispghf - whether spin projection on GHF wavefunction is on
C       n?c    - indices for complex conjugation projection
C       na     - number of alpha electrons
C       nb     - number of beta electrons
C       norbs  - number of orbitals
C       inobas - if .true., force the use of the true NO basis
C       p?nt   - N rows of density matrix in NO basis
C       p?tn   - N columns of density matrix in NO basis
C       p2?nt  - N rows of density matrix in NO basis (P*)
C       p2?tn  - N columns of density matrix in NO basis (P*)
C       scr1?  - scratch array
C       scr2?  - scratch array
C       dg1    - rotation matrix for angle g
C       piv    - scratch array
C       ngg    - matrix N(g,g')  [ out ]

      integer     iout, iprint, itype, n1c, n2c, na, nb, norbs
      real*8      pnt(*), ptn(*), pbnt(*), pbtn(*)
      real*8      p2nt(*), p2tn(*), p2bnt(*), p2btn(*)
      real*8      scr1(*), scr2(*), scr1b(*), scr2b(*)
      real*8      dg1(*), piv(*), ngg(*)
      logical     ispuhf, ispghf, inobas

C     other variables

      integer     nel, nosq, info

C     constants

      complex*16  zero, one


      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)


C     Construct useful variables.

      nel = na + nb
      nosq = norbs * norbs



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Construct the matrices N(g,g')  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Construct the matrices N(g,g') as described above...


C     Compute the product:

C       ( T1_oo  T1_ov )  =
C         ( P_oo   P_ov ) . ( Dg'_oo   Dg'_ov )! .  (Dg_oo   Dg_ov),
C                           ( Dg'_vo   Dg'_vv )     (Dg_vo   Dg_vv)

C     or

C       ( T1_oo  T1_ov )  =
C         ( P_oo )^T . ( Dg'_oo   Dg'_ov )! . ( Dg_oo   Dg_ov ).
C         ( P_vo )     ( Dg'_vo   Dg'_vv )    ( Dg_vo   Dg_vv )

      if ( ispuhf .or. ispghf ) then
        if ( n1c .eq. 1 ) then
          call zgemm ('n', 'n', nel, 2*norbs, 2*norbs, one, pnt, nel,
     $         dg1, 2*norbs, zero, scr1, nel)

        elseif ( n1c .eq. 2 .and. inobas ) then
          call zgemm ('n', 'n', nel, 2*norbs, 2*norbs, one, p2nt, nel,
     $         dg1, 2*norbs, zero, scr1, nel)

        elseif ( n1c .eq. 2 .and. .not. inobas ) then
          call zgemm ('t', 'n', nel, 2*norbs, 2*norbs, one, ptn,
     $         2*norbs, dg1, 2*norbs, zero, scr1, nel)
        endif

      else
        if ( n1c .eq. 1 ) then
          if ( itype .eq. 1 ) then
            call amove (2*na*norbs, pnt, scr1)

          elseif ( itype .eq. 2 ) then
            call amove (2*na*norbs, pnt, scr1)

            if ( nb .ge. 1 ) then
            call amove (2*nb*norbs, pbnt, scr1b)
            endif

          elseif ( itype .eq. 3 ) then
            call amove (4*nel*norbs, pnt, scr1)
          endif

        elseif ( n1c .eq. 2 .and. inobas ) then
          if ( itype .eq. 1 ) then
            call amove (2*na*norbs, p2nt, scr1)

          elseif ( itype .eq. 2 ) then
            call amove (2*na*norbs, p2nt, scr1)

            if ( nb .ge. 1 ) then
            call amove (2*nb*norbs, p2bnt, scr1b)
            endif

          elseif ( itype .eq. 3 ) then
            call amove (4*nel*norbs, p2nt, scr1)
          endif

        elseif ( n1c .eq. 2 .and. .not. inobas ) then
          if ( itype .eq. 1 ) then
            call amove (2*na*norbs, ptn, scr2)

            call trprct (2, norbs, na, scr2, scr1)

          elseif ( itype .eq. 2 ) then
            call amove (2*na*norbs, ptn, scr2)

            if ( nb .ge. 1 ) then
            call amove (2*nb*norbs, pbtn, scr2b)
            endif

            call trprct (2, norbs, na, scr2, scr1)

            if ( nb .ge. 1 ) then
            call trprct (2, norbs, nb, scr2b, scr1b)
            endif

          elseif ( itype .eq. 3 ) then
            call amove (4*nel*norbs, ptn, scr2)

            call trprct (2, 2*norbs, nel, scr2, scr1)
          endif
        endif
      endif


C     Compute the product:

C       ( T2_oo )  =  ( T1_oo   T1_ov ) . ( P_oo ),
C                                         ( P_vo )

C     or

C       ( T2_oo )  =  ( T1_oo   T1_ov ) . ( P_oo   P_ov )^T

C     Note that the product T2 is stored in N(g,g'), in preparation for
C     inversion...

      if ( n2c .eq. 1 ) then
        if ( itype .eq. 1 ) then
          call zgemm ('n', 'n', na, na, norbs, one, scr1, na, ptn,
     $         norbs, zero, ngg, na)

        elseif ( itype .eq. 2 ) then
          call zgemm ('n', 'n', na, na, norbs, one, scr1, na, ptn,
     $         norbs, zero, ngg, na)

          if ( nb .ge. 1 ) then
          call zgemm ('n', 'n', nb, nb, norbs, one, scr1b, nb, pbtn,
     $         norbs, zero, ngg(1+2*na*na), nb)
          endif

        elseif ( itype .eq. 3 ) then
          call zgemm ('n', 'n', nel, nel, 2*norbs, one, scr1, nel, ptn,
     $         2*norbs, zero, ngg, nel)
        endif

      elseif ( n2c .eq. 2 .and. inobas ) then
        if ( itype .eq. 1 ) then
          call zgemm ('n', 'n', na, na, norbs, one, scr1, na, p2tn,
     $         norbs, zero, ngg, na)

        elseif ( itype .eq. 2 ) then
          call zgemm ('n', 'n', na, na, norbs, one, scr1, na, p2tn,
     $         norbs, zero, ngg, na)

          if ( nb .ge. 1 ) then
          call zgemm ('n', 'n', nb, nb, norbs, one, scr1b, nb, p2btn,
     $         norbs, zero, ngg(1+2*na*na), nb)
          endif

        elseif ( itype .eq. 3 ) then
          call zgemm ('n', 'n', nel, nel, 2*norbs, one, scr1, nel, p2tn,
     $         2*norbs, zero, ngg, nel)
        endif

      elseif ( n2c .eq. 2 .and. .not. inobas ) then
        if ( itype .eq. 1 ) then
          call zgemm ('n', 't', na, na, norbs, one, scr1, na, pnt, na,
     $         zero, ngg, na)

        elseif ( itype .eq. 2 ) then
          call zgemm ('n', 't', na, na, norbs, one, scr1, na, pnt, na,
     $         zero, ngg, na)

          if ( nb .ge. 1 ) then
          call zgemm ('n', 't', nb, nb, norbs, one, scr1b, nb, pbnt, nb,
     $         zero, ngg(1+2*na*na), nb)
          endif

        elseif ( itype .eq. 3 ) then
          call zgemm ('n', 't', nel, nel, 2*norbs, one, scr1, nel, pnt,
     $         nel, zero, ngg, nel)
        endif
      endif


C     Evaluate the matrix N(g,g'):
C       N(g,g') = ( T2 )^{-1}.

      if ( itype .eq. 1 ) then
        call zgetrf (na, na, ngg, na, piv, info)
        if ( info .ne. 0 ) call gauerr ('info .ne. 0 in bldngg')

        call zgetri (na, ngg, na, piv, scr1, na*na, info)
        if ( info .ne. 0 ) call gauerr ('info .ne. 0 in bldngg')

      elseif ( itype .eq. 2 ) then
        call zgetrf (na, na, ngg, na, piv, info)
        if ( info .ne. 0 ) call gauerr ('info .ne. 0 in bldngg')

        call zgetri (na, ngg, na, piv, scr1, na*na, info)
        if ( info .ne. 0 ) call gauerr ('info .ne. 0 in bldngg')

        if ( nb .ge. 1 ) then
        call zgetrf (nb, nb, ngg(1+2*na*na), nb, piv, info)
        if ( info .ne. 0 ) call gauerr ('info .ne. 0 in bldngg')

        call zgetri (nb, ngg(1+2*na*na), nb, piv, scr1b, nb*nb, info)
        if ( info .ne. 0 ) call gauerr ('info .ne. 0 in bldngg')
        endif

      elseif ( itype .eq. 3 ) then
        call zgetrf (nel, nel, ngg, nel, piv, info)
        if ( info .ne. 0 ) call gauerr ('info .ne. 0 in bldngg')

        call zgetri (nel, ngg, nel, piv, scr1, nel*nel, info)
        if ( info .ne. 0 ) call gauerr ('info .ne. 0 in bldngg')
      endif


      return
      end


