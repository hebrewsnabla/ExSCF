      

      subroutine bldpgg (iout, iprint, itype, ispuhf, ispghf, n1c, n2c,
     $     na, nb, norbs, inobas, pnt, ptn, pbnt, pbtn, p2nt, p2tn,
     $     p2bnt, p2btn, scr1, scr1b, scr2, scr2b, dg1, ngg, pgg)

      implicit none


C +------------------------------------------------------------+
C |  bldpgg   --   CAJH, 08.2011                               |
C |                                                            |
C |                                                            |
C |  Build the transition density matrices, P(g,g'), in the    |
C |  same way as formpg forms the matrices P(g).               |
C |                                                            |
C |  Given a rotation matrix Dg, related to the spin           |
C |  projection operator with an angle Omega, the matrices     |
C |  P(g,g') can take the form                                 |
C |                                                            |
C |  if nu = I, mu = I  ( <0| I! R(Omega')! R(Omega) I |0> )   |
C |                                                            |
C |    P(g,g')  =                                              |
C |      Dg . ( P_oo ) . N(g,g') . ( P_oo  P_ov ) . Dg'!       |
C |           ( P_vo )                                         |
C |                                                            |
C |  if nu = K, mu = I  ( <0| I! R(Omega')! R(Omega) K |0> )   |
C |                                                            |
C |    P(g,g')  =                                              |
C |  Dg . ( P_oo  P_ov )^T . N(g,g') . ( P_oo  P_ov ) . Dg'!   |
C |                                                            |
C |  if nu = I, mu = K  ( <0| K! R(Omega')! R(Omega) I |0> )   |
C |                                                            |
C |    P(g,g')  =                                              |
C |      Dg . ( P_oo ) . N(g,g') . ( P_oo )^T . Dg'!           |
C |           ( P_vo )             ( P_vo )                    |
C |                                                            |
C |  if nu = K, mu = K  ( <0| K! R(Omega')! R(Omega) K |0> )   |
C |                                                            |
C |    P(g,g')  =                                              |
C |      Dg . ( P_oo  P_ov )^T . N(g,g') . ( P_oo )^T . Dg'!   |
C |                                        ( P_vo )            |
C |                                                            |
C |  where nu, mu = I, K are operators related to complex      |
C |  conjugation restoration. (See evalsmt for more details.)  |
C |  The matrices N(g,g') are those built in bldngg.           |
C |                                                            |
C |  Other than that, it is worth recalling that P is the      |
C |  density matrix of the unprojected HF state. The indices   |
C |  o and v define the dimensions (occupied, virtual) of the  |
C |  rectangular matrices. The Dg matrices are always          |
C |  of dimension nbasis x nbasis, whereas the N(g,g')         |
C |  matrices are of dimension N x N, where N is the number    |
C |  of electrons in the system.                               |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  CAJH, 08.2011                                             |
C |                                                            |
C |  Major modification of the subroutine. It now computes     |
C |  the matrix P(g,g') for a single integration point g,g'.   |
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
C |      ngg  -  2*na*na        [ in ]                         |
C |      pgg  -  2*nosq         [ out ]                        |
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
C |      ngg   -  2*na*na + 2*nb*nb   [ in ]                   |
C |      pgg   -  4*nosq        [ out ]                        |
C |                                                            |
C |    itype = 3,                                              |
C |                                                            |
C |      pnt  -  2*nel*norbs    [ in ]                         |
C |      ptn  -  2*nel*norbs    [ in ]                         |
C |      scr1 -  2*nel*norbs                                   |
C |      scr2 -  2*nel*norbs                                   |
C |      ngg  -  2*nel*nel      [ in ]                         |
C |      pgg  -  8*nosq         [ out ]                        |
C |                                                            |
C |    if ( ispuhf .or. ispghf )                               |
C |                                                            |
C |      dg1   -  8*nosq                                       |
C |      dg2   -  8*nosq                                       |
C |      scrdg -  8*nosq                                       |
C |                                                            |
C |  Here,                                                     |
C |                                                            |
C |    p?nt - N rows    of density matrix in NO basis          |
C |    p?tn - N columns of density matrix in NO basis          |
C |    ngg  - input matrix N(g,g')                             |
C |    pgg  - output matrix P(g,g')                            |
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
C |  scrdg is no longer declared.                              |
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
C       ngg    - matrix N(g,g')
C       pgg    - matrix P(g,g')  [ out ]

      integer     iout, iprint, itype, n1c, n2c, na, nb, norbs
      real*8      pnt(*), ptn(*), pbnt(*), pbtn(*)
      real*8      p2nt(*), p2tn(*), p2bnt(*), p2btn(*)
      real*8      scr1(*), scr2(*), scr1b(*), scr2b(*)
      real*8      dg1(*), ngg(*), pgg(*)
      logical     ispuhf, ispghf, inobas

C     other variables

      integer     nel, nosq

C     constants

      complex*16  zero, one


      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)


C     Construct useful variables.

      nel = na + nb
      nosq = norbs * norbs



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Construct the matrices P(g,g')  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Construct the matrices P(g,g') as described above...


C     Compute the product:

C       ( T1_oo )  =  ( Dg_oo   Dg_ov ) . ( P_oo )
C       ( T1_vo )     ( Dg_vo   Dg_vv )   ( P_vo )

C     or

C       ( T1_oo )  =  ( Dg_oo   Dg_ov ) . ( P_oo   P_ov )^T
C       ( T1_vo )     ( Dg_vo   Dg_vv )

      if ( ispuhf .or. ispghf ) then
        if ( n2c .eq. 1 ) then
          call zgemm ('n', 'n', 2*norbs, nel, 2*norbs, one, dg1,
     $         2*norbs, ptn, 2*norbs, zero, scr1, 2*norbs)

        elseif ( n2c .eq. 2 .and. inobas ) then
          call zgemm ('n', 'n', 2*norbs, nel, 2*norbs, one, dg1,
     $         2*norbs, p2tn, 2*norbs, zero, scr1, 2*norbs)

        elseif ( n2c .eq. 2 .and. .not. inobas ) then
          call zgemm ('n', 't', 2*norbs, nel, 2*norbs, one, dg1,
     $         2*norbs, pnt, nel, zero, scr1, 2*norbs)
        endif

      else
        if ( n2c .eq. 1 ) then
          if ( itype .eq. 1 ) then
            call amove (2*na*norbs, ptn, scr1)

          elseif ( itype .eq. 2 ) then
            call amove (2*na*norbs, ptn, scr1)

            if ( nb .ge. 1 ) then
            call amove (2*nb*norbs, pbtn, scr1b)
            endif

          elseif ( itype .eq. 3 ) then
            call amove (4*nel*norbs, ptn, scr1)
          endif

        elseif ( n2c .eq. 2 .and. inobas ) then
          if ( itype .eq. 1 ) then
            call amove (2*na*norbs, p2tn, scr1)

          elseif ( itype .eq. 2 ) then
            call amove (2*na*norbs, p2tn, scr1)

            if ( nb .ge. 1 ) then
            call amove (2*nb*norbs, p2btn, scr1b)
            endif

          elseif ( itype .eq. 3 ) then
            call amove (4*nel*norbs, p2tn, scr1)
          endif

        elseif ( n2c .eq. 2 .and. .not. inobas ) then
          if ( itype .eq. 1 ) then
            call amove (2*na*norbs, pnt, scr2)

            call trprct (2, na, norbs, scr2, scr1)

          elseif ( itype .eq. 2 ) then
            call amove (2*na*norbs, pnt, scr2)

            if ( nb .ge. 1 ) then
            call amove (2*nb*norbs, pbnt, scr2b)
            endif

            call trprct (2, na, norbs, scr2, scr1)

            if ( nb .ge. 1 ) then
            call trprct (2, nb, norbs, scr2b, scr1b)
            endif

          elseif ( itype .eq. 3 ) then
            call amove (4*nel*norbs, pnt, scr2)

            call trprct (2, nel, 2*norbs, scr2, scr1)
          endif
        endif
      endif


C     Compute the product:

C       ( T2_oo )  =  ( T1_oo ) . N(g,g')
C       ( T2_vo )     ( T1_vo )

      if ( itype .eq. 1 ) then
        call zgemm ('n', 'n', norbs, na, na, one, scr1, norbs,
     $       ngg, na, zero, scr2, norbs)

      elseif ( itype .eq. 2 ) then
        call zgemm ('n', 'n', norbs, na, na, one, scr1, norbs,
     $       ngg, na, zero, scr2, norbs)

        if ( nb .ge. 1 ) then
        call zgemm ('n', 'n', norbs, nb, nb, one, scr1b, norbs,
     $       ngg(1+2*na*na), nb, zero, scr2b, norbs)
        endif

      elseif ( itype .eq. 3 ) then
        call zgemm ('n', 'n', 2*norbs, nel, nel, one, scr1, 2*norbs,
     $       ngg, nel, zero, scr2, 2*norbs)
      endif


C     Compute the product:

C       ( T3_oo   T3_ov )  =  ( T2_oo ) . ( P_oo   P_ov )
C       ( T3_vo   T3_vv )     ( T2_vo )

C     or

C       ( T3_oo   T3_ov )  =  ( T2_oo ) . ( P_oo )^T
C       ( T3_vo   T3_vv )     ( T2_vo )   ( P_vo )

C     ** Note that this product is already stored in ipgg.

      if ( n1c .eq. 1 ) then
        if ( itype .eq. 1 ) then
          call zgemm ('n', 'n', norbs, norbs, na, one, scr2, norbs,
     $         pnt, na, zero, pgg, norbs)

        elseif ( itype .eq. 2 ) then
          call zgemm ('n', 'n', norbs, norbs, na, one, scr2, norbs,
     $         pnt, na, zero, pgg, norbs)

          if ( nb .ge. 1 ) then
          call zgemm ('n', 'n', norbs, norbs, nb, one, scr2b, norbs,
     $         pbnt, nb, zero, pgg(1+2*nosq), norbs)
          else
            call aclear (2*nosq, pgg(1+2*nosq))
          endif

        elseif ( itype .eq. 3 ) then
          call zgemm ('n', 'n', 2*norbs, 2*norbs, nel, one, scr2,
     $         2*norbs, pnt, nel, zero, pgg, 2*norbs)
        endif

      elseif ( n1c .eq. 2 .and. inobas ) then
        if ( itype .eq. 1 ) then
          call zgemm ('n', 'n', norbs, norbs, na, one, scr2, norbs,
     $         p2nt, na, zero, pgg, norbs)

        elseif ( itype .eq. 2 ) then
          call zgemm ('n', 'n', norbs, norbs, na, one, scr2, norbs,
     $         p2nt, na, zero, pgg, norbs)

          if ( nb .ge. 1 ) then
          call zgemm ('n', 'n', norbs, norbs, nb, one, scr2b, norbs,
     $         p2bnt, nb, zero, pgg(1+2*nosq), norbs)
          else
            call aclear (2*nosq, pgg(1+2*nosq))
          endif

        elseif ( itype .eq. 3 ) then
          call zgemm ('n', 'n', 2*norbs, 2*norbs, nel, one, scr2,
     $         2*norbs, p2nt, nel, zero, pgg, 2*norbs)
        endif

      elseif ( n1c .eq. 2 .and. .not. inobas ) then
        if ( itype .eq. 1 ) then
          call zgemm ('n', 't', norbs, norbs, na, one, scr2, norbs,
     $         ptn, norbs, zero, pgg, norbs)

        elseif ( itype .eq. 2 ) then
          call zgemm ('n', 't', norbs, norbs, na, one, scr2, norbs,
     $         ptn, norbs, zero, pgg, norbs)

          if ( nb .ge. 1 ) then
          call zgemm ('n', 't', norbs, norbs, nb, one, scr2b, norbs,
     $         pbtn, norbs, zero, pgg(1+2*nosq), norbs)
          else
            call aclear (2*nosq, pgg(1+2*nosq))
          endif

        elseif ( itype .eq. 3 ) then
          call zgemm ('n', 't', 2*norbs, 2*norbs, nel, one, scr2,
     $         2*norbs, ptn, 2*norbs, zero, pgg, 2*norbs)
        endif
      endif


      return
      end


