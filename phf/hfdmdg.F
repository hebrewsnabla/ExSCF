

      subroutine hfdmdg (iout, iprint, v, lenv, iphf, na, nb, norbs,
     $     irwnos, p, pno, pstrno, dodiag, inobas)

      implicit none


C +------------------------------------------------------------+
C |  hfdmdg   --   CAJH, 07.2011                               |
C |                                                            |
C |                                                            |
C |  This routine diagonalizes the density matrix to find the  |
C |  NO basis. It stores the orthonormal AO => NO basis        |
C |  transformation matrix in the RWF file irwnos. It also     |
C |  stores the density matrix in the NO basis in square       |
C |  format in the array pno.                                  |
C |                                                            |
C |  IMPORTANT NOTE: For icmplx = .true., we do not find       |
C |  the NO basis. Instead, we diagonalize real(P), such       |
C |  that the "NO" basis is an orthogonal (as opposed to a     |
C |  unitary matrix). This simplifies the basis set            |
C |  transformations for other subroutines.                    |
C |                                                            |
C |    In particular, in the equations involving complex       |
C |    conjugation restoration, P and P* appear frequently.    |
C |    In the orthonormal AO basis where they are prepared,    |
C |    they are related by                                     |
C |                                                            |
C |      P* = P^T,   because  P = P!                           |
C |                                                            |
C |    However, in a different basis, this relationship is     |
C |    no longer valid. To see this, let U be a unitary        |
C |    matrix. We then have                                    |
C |                                                            |
C |      (P_NO)  = U! . P . U                                  |
C |      (P_NO)* = U^T . P* . U*                               |
C |                                                            |
C |    However, PNO* does not correspond to P* being           |
C |    transformed from one basis to another. Such             |
C |    transformation should be given by                       |
C |                                                            |
C |      (P_NO)  = U! . P . U                                  |
C |      (P*_NO) = U! . P* . U                                 |
C |                                                            |
C |  In summary, one needs to be careful when using unitary    |
C |  matrices for basis transformations. On the other hand,    |
C |  the relationships above DO hold for orthogonal            |
C |  transformations. We therefore prefer to deal with         |
C |  an orthogonal matrix that simplifies later                |
C |  developments.                                             |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  CAJH, 08.2011                                             |
C |                                                            |
C |  Added the input argument dodiag, which is set to .false.  |
C |  can avoid the diagonalization of the density matrix.      |
C |  Instead, it uses the set of eigenvectors already stored   |
C |  in the appropriate RWF file.                              |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  CAJH, 08.2011                                             |
C |                                                            |
C |  The subroutine has now the capability to force the use    |
C |  of the actual NO basis (as opposed to the real part       |
C |  of the NO basis). Note that we now have to output both    |
C |  P and P* in the NO basis.                                 |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       na     - number of alpha electrons
C       nb     - number of beta electrons
C       norbs  - number of orbitals
C       irwnos - RWF file with transformation matrix from orthonormal AO
C                to NO basis [ updated here ]
C       p      - density matrix in orthonormal AO basis
C       pno    - density matrix in NO basis [ out ]
C       pstrno - density matrix (P*) in NO basis [ out ]
C       dodiag - whether to diagonalize the density matrix (.true.) or
C                just use the set of eigenvectors already stored in irwnos
C       inobas - if .true., force the use of the true NO basis

      real*8      v(*), p(*), pno(*), pstrno(*)
      integer     iout, iprint, lenv, iphf, na, nb, norbs, irwnos
      logical     dodiag, inobas

C     other variables

      integer     itype, nel, ntto, ntt2o, nosq, info
      integer     jend, mdv
      integer     jp, jpstr, jvec, jval, jt
      integer     jscr1, jscr2, jscrx, jscrd, lscrd, jzab, jz
      integer     jz1, jz2, jz3
      integer     ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8
      logical     icmplx, ispuhf, ispghf, useno
      integer     ijvec, jjvec

C     constants

      complex*16  zero, one


      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)


C     Construct useful variables.

      nel = na + nb

      ntto = norbs * (norbs + 1) / 2
      ntt2o = norbs * (2*norbs + 1)

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
        call gauerr ('Unrecognized iphf in hfdmdg.')
      endif



C     Skip the section of the code that diagonalizes the density matrix
C     if it is not needed.

      if ( .not. dodiag ) then
        goto 200
      endif



C     Decide whether we will use the true NO basis (useno = .true.) or
C     the basis that diagonalizes the real part of P.

      useno = inobas .or. .not. icmplx


C     Memory allocation.

C       Allocate space to store:
C         - density matrix as square matrix (jp)
C         - eigenvalues of density matrix (jval)
C         - eigenvectors of density matrix (jvec)
C         - scratch space for matrix multiplications (jt)

      if ( itype .eq. 1 ) then
        jp   = 1
        jval = jp   + 2*nosq
        jvec = jval + norbs
        jt   = jvec + 2*nosq
        jend = jt   + 2*nosq

      elseif ( itype .eq. 2 ) then
        jp   = 1
        jval = jp   + 4*nosq
        jvec = jval + 2*norbs
        jt   = jvec + 4*nosq
        jend = jt   + 2*nosq

      elseif ( itype .eq. 3 ) then
        jp   = 1
        jval = jp   + 8*nosq
        jvec = jval + 2*norbs
        jt   = jvec + 8*nosq
        jend = jt   + 8*nosq
      endif

C       Allocate space to store:
C         - scratch arrays for calls to zhpev, dspev (jscr1, jscr2)
C         - scratch array to perform real to complex transformation of
C           eigenvectors (jscrx)

      if ( useno ) then
        if ( itype .lt. 3 .or. ispuhf ) then
          jscr1 = jend
          jscr2 = jscr1 + 2*(2*norbs - 1)
          jend  = jscr2 + 3*norbs - 2

        elseif ( itype .eq. 3 ) then
          jscr1 = jend
          jscr2 = jscr1 + 2*(4*norbs - 1)
          jend  = jscr2 + 6*norbs - 2
        endif

      else
        if ( itype .lt. 3 .or. ispuhf ) then
          jscr1 = jend
          jscrx = jscr1 + 3*norbs
          jend  = jscrx + nosq

        elseif ( itype .eq. 3 ) then
          jscr1 = jend
          jscrx = jscr1 + 6*norbs
          jend  = jscrx + 4*nosq
        endif
      endif

C     Allocate space for some ispuhf specific transformations:
C       - scratch array for calls to dmblck (jscrd, jzab)
C       - scratch array to reorganize eigenvectors (jz)

      if ( ispuhf ) then
        lscrd = 16*nosq

        jzab  = jend
        jscrd = jzab  + 2*nosq
        jz    = jscrd + lscrd
        jend  = jz    + 8*nosq
      endif

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'hfdmdg')

C      write(iout, *) 'End mem alloc in hfdmdg'

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Diagonalize density matrix  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Copy density matrix into array p.
C     Scale the density matrix by a factor of -1 (so that the largest
C     eigenvalues come first after diagonalization).

C       If useno = .true., the full density matrix is copied.
C       If useno = .false., then only the real part is copied.

      if ( useno ) then
        if ( itype .eq. 1 ) then
          call amove (2*ntto, p, v(jp))
          call zdscal (ntto, -1.0d0, v(jp), 1)

        elseif ( itype .eq. 2 .or. ispuhf ) then
          call amove (4*ntto, p, v(jp))
          call zdscal (2*ntto, -1.0d0, v(jp), 1)

        elseif ( itype .eq. 3 ) then
          call amove (2*ntt2o, p, v(jp))
          call zdscal (ntt2o, -1.0d0, v(jp), 1)
        endif

      else
        if ( itype .eq. 1 ) then
          call zrecmp (4, p, v(jp), ntto)
          call ascale (ntto, -1.0d0, v(jp), v(jp))

        elseif ( itype .eq. 2 .or. ispuhf ) then
          call zrecmp (4, p, v(jp), 2*ntto)
          call ascale (2*ntto, -1.0d0, v(jp), v(jp))

        elseif ( itype .eq. 3 ) then
          call zrecmp (4, p, v(jp), ntt2o)
          call ascale (ntt2o, -1.0d0, v(jp), v(jp))
        endif
      endif
C      write(iout, *) 'End copy in hfdmdg'

C     debugging...
C       print density matrix to be diagonalized

c$$$      if ( useno ) then
c$$$        if ( itype .eq. 1 ) then
c$$$          write (iout, *) ' '
c$$$          call ltoutx (iout, 0, 2, 1, 0, 'density diag', norbs, norbs,
c$$$     $         v(jp), v(jp), 0)
c$$$
c$$$        elseif ( itype .eq. 2 .or. ispuhf ) then
c          DO ijvec = 0, 3
c           WRITE(*,665) v(jp)
c           WRITE(*,665) v(jp+2), v(jp+8)
c           WRITE(*,665) v(jp+4), v(jp+10), v(jp+14)
c           WRITE(*,665) v(jp+6), v(jp+12), v(jp+16), v(jp+18)
c          END DO
          write (iout, *) ' '
          call ltoutx (iout, 0, 3, 1, 0, 'dm to be diaged', norbs,norbs,
     $         v(jp), v(jp+2*ntto), 0)
c          write (iout, *) ' dm to be diaged'
c          call outmtc (iout, 3, 0, 2, 1, 2, ! 'dm to be diaged', norbs,norbs,
c     $         v(jp) ,!v(jp+2*ntto), 0)
c     $         norbs, norbs, norbs, norbs)
c$$$
c$$$        elseif ( itype .eq. 3 ) then
c$$$          write (iout, *) ' '
c$$$          call ltoutx (iout, 0, 2, 1, 0, 'density diag', 2*norbs,
c$$$     $         2*norbs, v(jp), v(jp), 0) 
c$$$        endif
c$$$
c$$$      else
c$$$        if ( itype .eq. 1 ) then
c$$$          write (iout, *) ' '
c$$$          call ltoutx (iout, 0, 0, 1, 0, 'density diag', norbs, norbs,
c$$$     $         v(jp), v(jp), 0)
c$$$
c$$$        elseif ( itype .eq. 2 .or. ispuhf ) then
c$$$          write (iout, *) ' '
c$$$          call ltoutx (iout, 0, 1, 1, 0, 'density diag', norbs, norbs,
c$$$     $         v(jp), v(jp+ntto), 0)
c$$$
c$$$        elseif ( itype .eq. 3 ) then
c$$$          write (iout, *) ' '
c$$$          call ltoutx (iout, 0, 0, 1, 0, 'density diag', 2*norbs,
c$$$     $         2*norbs, v(jp), v(jp), 0)
c$$$        endif
c$$$      endif


C     Diagonalize the density matrix or the real part of the density
C     matrix.

C       .. the set of eigenvalues are stored in jval
C       .. the set of eigenvectors are stored in jvec

      if ( useno ) then

        if ( itype .eq. 1 ) then
          call vconjg (ntto, v(jp))
          call zhpev ('v', 'u', norbs, v(jp), v(jval), v(jvec), norbs,
     $         v(jscr1), v(jscr2), info)

          if ( info .ne. 0 ) then
            call gauerr ('zhpev failed in hfdmdg.')
          endif

        elseif ( itype .eq. 2 .or. ispuhf ) then
C          write(iout, *) 'Enter itype=2'
          call vconjg (ntto, v(jp))
C          write(iout, *) 'End vconjg1'
          call zhpev ('v', 'u', norbs, v(jp), v(jval), v(jvec), norbs,
     $         v(jscr1), v(jscr2), info)
C          write(iout, *) 'End zhpev1'

          if ( info .ne. 0 ) then
            call gauerr ('zhpev failed in hfdmdg.')
          endif

          call vconjg (ntto, v(jp+2*ntto))
C          write(iout, *) 'End vconjg2'
          call zhpev ('v', 'u', norbs, v(jp+2*ntto), v(jval+norbs),
     $         v(jvec+2*nosq), norbs, v(jscr1), v(jscr2), info)
C          write(iout, *) 'End zhpev2 '

          if ( info .ne. 0 ) then
            call gauerr ('zhpev failed in hfdmdg.')
          endif

        elseif ( itype .eq. 3 ) then
          call vconjg (ntt2o, v(jp))
          call zhpev ('v', 'u', 2*norbs, v(jp), v(jval), v(jvec),
     $         2*norbs, v(jscr1), v(jscr2), info)

          if ( info .ne. 0 ) then
            call gauerr ('zhpev failed in hfdmdg.')
          endif
        endif

      else

        if ( itype .eq. 1 ) then
          call dspev ('v', 'u', norbs, v(jp), v(jval), v(jvec), norbs,
     $         v(jscr1), info)

          if ( info .ne. 0 ) then
            call gauerr ('dspev failed in hfdmdg.')
          endif

        elseif ( itype .eq. 2 .or. ispuhf ) then
          call dspev ('v', 'u', norbs, v(jp), v(jval), v(jvec), norbs,
     $         v(jscr1), info)

          if ( info .ne. 0 ) then
            call gauerr ('dspev failed in hfdmdg.')
          endif

          call dspev ('v', 'u', norbs, v(jp+ntto), v(jval+norbs),
     $         v(jvec+2*nosq), norbs, v(jscr1), info)

          if ( info .ne. 0 ) then
            call gauerr ('dspev failed in hfdmdg.')
          endif

        elseif ( itype .eq. 3 ) then
          call dspev ('v', 'u', 2*norbs, v(jp), v(jval), v(jvec),
     $         2*norbs, v(jscr1), info)

          if ( info .ne. 0 ) then
            call gauerr ('dspev failed in hfdmdg.')
          endif
        endif
      endif

C      write(iout, *) 'End diag1 in hfdmdg'

C     debugging...
C       print eigenvectors of density matrix

      if ( useno ) then
        if ( itype .eq. 1 ) then
          write (iout, *) ' '
          call ltoutx (iout, 1, 2, 1, 0, 'eigenvectors', norbs, norbs,
     $         v(jvec), v(jvec), 0)

        elseif ( itype .eq. 2 .or. ispuhf ) then
          write (iout, *) ' '
          call ltoutx (iout, 1, 3, 1, 0, 'eigenvectors', norbs, norbs,
     $         v(jvec), v(jvec+2*nosq), 0)
c          DO ijvec = 0, 3
c           WRITE(*,665) ( v(jvec + ijvec*8 + jjvec ), jjvec = 0, 7, 2 )
c          END DO
665   FORMAT( 11(:,1X,E23.15) )
          write(iout, 667) v(jval), v(jval+1), v(jval+2), v(jval+3)
          write(iout, 667) v(jval+4), v(jval+5), v(jval+6), v(jval+7)
666   format('eigenvalues ', F12.6, F12.6, F12.6, F12.6) 
667   format('eigenvalues ', 4E23.15) 
C          write (iout, *) ' '
C          call ltoutx (iout, 0, 3, 1, 0, 'dm', norbs, norbs,
C     $         v(jp), v(jp+2*ntto), 0)
        elseif ( itype .eq. 3 ) then
          write (iout, *) ' '
          call ltoutx (iout, 1, 2, 1, 0, 'eigenvectors', 2*norbs,
     $         2*norbs, v(jvec), v(jvec), 0) 
        endif

      else
        if ( itype .eq. 1 ) then
          write (iout, *) ' '
          call ltoutx (iout, 1, 0, 1, 0, 'eigenvectors', norbs, norbs,
     $         v(jvec), v(jvec), 0)

        elseif ( itype .eq. 2 .or. ispuhf ) then
          write (iout, *) ' '
          call ltoutx (iout, 1, 1, 1, 0, 'eigenvectors', norbs, norbs,
     $         v(jvec), v(jvec+2*nosq), 0)

        elseif ( itype .eq. 3 ) then
          write (iout, *) ' '
          call ltoutx (iout, 1, 0, 1, 0, 'eigenvectors', 2*norbs,
     $         2*norbs, v(jvec), v(jvec), 0)
        endif
      endif


C     Transform the eigenvectors to complex if they are real.

      if ( .not. useno ) then
        if ( itype .eq. 1 ) then
          call amove (nosq, v(jvec), v(jscrx))
          call zrecmp (3, v(jvec), v(jscrx), nosq)

        elseif ( itype .eq. 2 .or. ispuhf ) then
          call amove (nosq, v(jvec), v(jscrx))
          call zrecmp (3, v(jvec), v(jscrx), nosq)

          call amove (nosq, v(jvec+2*nosq), v(jscrx))
          call zrecmp (3, v(jvec+2*nosq), v(jscrx), nosq)

        elseif ( itype .eq. 3 ) then
          call amove (4*nosq, v(jvec), v(jscrx))
          call zrecmp (3, v(jvec), v(jscrx), 4*nosq)
        endif
      endif


C     For ispuhf = .true., there is some reorganization to do. First, we
C     need to expand the set of eigenvectors to the full spin-orbital
C     basis. (We currently only have the alpha-alpha and the beta-beta
C     blocks of the set of eigenvectors stored separately.) After that,
C     we need to order the eigenvectors so that the eigenvalues of full
C     P are approximately in descending order.

C       We currently take a simplistic approach, and organize the
C       eigenvectors as follows (in terms of eigenvalues):

C         1 .. na               - largest na alpha eigenvalues
C         na+1 .. nel           - largest nb beta eigenvalues
C         nel+1 .. norbs+nb     - remaining norbs-na alpha eigenvalues
C         norbs+nb+1 .. 2*norbs - remaining norbs-nb beta eigenvalues

      if ( ispuhf ) then
        call aclear (2*nosq, v(jzab))

        call dmblck (iout, iprint, v(jscrd), lscrd, norbs, v(jz),
     $       v(jvec), v(jzab), v(jzab), v(jvec+2*nosq), 1, -1)

        ind1 = jz
        ind2 = jz + 4*nosq
        ind3 = jz + 4*norbs*na
        ind4 = jz + 4*nosq + 4*norbs*nb

        ind5 = jvec
        ind6 = jvec + 4*norbs*na
        ind7 = jvec + 4*norbs*nel
        ind8 = jvec + 4*nosq + 4*norbs*nb

        call amove (4*norbs*na, v(ind1), v(ind5))
        call amove (4*norbs*nb, v(ind2), v(ind6))
        call amove (4*norbs*(norbs-na), v(ind3), v(ind7))
        call amove (4*norbs*(norbs-nb), v(ind4), v(ind8))
      endif


C     Save orthonormal AO => NO basis transformation matrix, i.e., the
C     eigenvectors that we just computed.

      if ( itype .eq. 1 ) then
        call fileio (1, -irwnos, 2*nosq, v(jvec), 0)

      elseif ( itype .eq. 2 ) then
        call fileio (1, -irwnos, 4*nosq, v(jvec), 0)

      elseif ( itype .eq. 3 ) then
        call fileio (1, -irwnos, 8*nosq, v(jvec), 0)
        write (iout, *) ' '
        call ltoutx (iout, 3, 2, 1, 0, 'expanded eigenvectors', 2*norbs,
     $ 2*norbs, v(jvec), v(jvec+8*nosq), 0)
      endif

C      write(iout, *) 'End diag in hfdmdg'
      goto 200



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Transform density matrix to NO basis  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 200  continue


C     Memory allocation.

C       Allocate space to store:
C         - density matrix as square matrix (jp)
C         - eigenvectors of density matrix (jvec)
C         - scratch space for matrix multiplications (jt)

      if ( itype .eq. 1 ) then
        jp   = 1
        jvec = jp   + 2*nosq
        jt   = jvec + 2*nosq
        jend = jt   + 2*nosq

      elseif ( itype .eq. 2 ) then
        jp   = 1
        jvec = jp   + 4*nosq
        jt   = jvec + 4*nosq
        jend = jt   + 2*nosq

      elseif ( itype .eq. 3 ) then
        jp   = 1
        jvec = jp   + 8*nosq
        jt   = jvec + 8*nosq
        jend = jt   + 8*nosq
      endif

C       Allocate space to store:
C         - density matrix P* as square matrix (jp)

      if ( inobas .and. itype .eq. 1 ) then
        jpstr = jend
        jend  = jpstr + 2*nosq

      elseif ( inobas .and. itype .eq. 2 ) then
        jpstr = jend
        jend  = jpstr + 4*nosq

      elseif ( inobas .and. itype .eq. 3 ) then
        jpstr = jend
        jend  = jpstr + 8*nosq
      endif


C       Allocate space to store:
C         - scratch arrays for lt to square transformation (jz1, jz2, jz3)

      if ( itype .eq. 1 ) then
        jz1  = jend
        jz2  = jz1  + 2*ntto
        jz3  = jz2  + nosq
        jend = jz3  + nosq

      elseif ( itype .eq. 2 .or. ispuhf ) then
        jz1  = jend
        jz2  = jz1  + 2*ntto
        jz3  = jz2  + nosq
        jend = jz3  + nosq

      elseif ( itype .eq. 3 ) then
        jz1  = jend
        jz2  = jz1  + 2*ntt2o
        jz3  = jz2  + 4*nosq
        jend = jz3  + 4*nosq
      endif

C     Allocate space for some ispuhf specific transformations:
C       - scratch array for calls to dmblck (jscrd, jzab)
C       - scratch array to form full density matrix (jz)

      if ( ispuhf ) then
        lscrd = 16*nosq

        jzab  = jend
        jscrd = jzab  + 2*nosq
        jz    = jscrd + lscrd
        jend  = jz    + 8*nosq
      endif

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'hfdmdg')


C     Load orthonormal AO => NO basis transformation matrix, i.e., the
C     eigenvectors that we just computed.

      if ( itype .eq. 1 ) then
        call fileio (2, -irwnos, 2*nosq, v(jvec), 0)

      elseif ( itype .eq. 2 ) then
        call fileio (2, -irwnos, 4*nosq, v(jvec), 0)

      elseif ( itype .eq. 3 ) then
        call fileio (2, -irwnos, 8*nosq, v(jvec), 0)
      endif


C     Transform the density matrix (orthonormal AO basis) into square
C     format.

      if ( itype .eq. 1 ) then
        call zrecmp (1, p, v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jp), v(jz2), nosq)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call zrecmp (1, p, v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jp), v(jz2), nosq)

        call zrecmp (1, p(1+2*ntto), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jp+2*nosq), v(jz2), nosq)

      elseif ( itype .eq. 3 ) then
        call zrecmp (1, p, v(jz1), ntt2o)
        call square (v(jz1), v(jz2), 2*norbs, 2*norbs, 0)
        call square (v(jz1+ntt2o), v(jz3), 2*norbs, 2*norbs, -1)
        call zrecmp (2, v(jp), v(jz2), 4*nosq)
      endif


C     For ispuhf = .true., we need to expand the density matrix to full
C     spin-orbital basis.

      if ( ispuhf ) then
        call aclear (2*nosq, v(jzab))

        call dmblck (iout, iprint, v(jscrd), lscrd, norbs, v(jz), v(jp),
     $       v(jzab), v(jzab), v(jp+2*nosq), 1, -1)

        call amove (8*nosq, v(jz), v(jp))
        write (iout, *) ' '
        call ltoutx (iout, 3, 2, 1, 0, 'expanded density matrix P', 
     $ 2*norbs, 2*norbs, v(jp), v(jp+8*nosq), 0)
      endif


C     Form P*.

      if ( inobas ) then
        if ( itype .eq. 1 ) then
          call amove (2*nosq, v(jp), v(jpstr))
          call vconjg (nosq, v(jpstr))

        elseif ( itype .eq. 2 ) then
          call amove (4*nosq, v(jp), v(jpstr))
          call vconjg (2*nosq, v(jpstr))

        elseif ( itype .eq. 3 ) then
          call amove (8*nosq, v(jp), v(jpstr))
          call vconjg (4*nosq, v(jpstr))
C          write (iout, *) ' '
C          call ltoutx (iout, 3, 2, 1, 0, 'expanded density matrix P*', 
C     $ 2*norbs, 2*norbs, v(jpstr), v(jpstr+8*nosq), 0)
        endif
      endif


C     Transform the density matrix into NO basis:
C       vec' . P . vec

      if ( itype .eq. 1 ) then
        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jvec), norbs,
     $       v(jp), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jvec), norbs, zero, v(jp), norbs)

      elseif ( itype .eq. 2 ) then
        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jvec), norbs,
     $       v(jp), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jvec), norbs, zero, v(jp), norbs)

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jvec+2*nosq),
     $       norbs, v(jp+2*nosq), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jvec+2*nosq), norbs, zero, v(jp+2*nosq), norbs)

      elseif ( itype .eq. 3 ) then
        call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jvec),
     $       2*norbs, v(jp), 2*norbs, zero, v(jt), 2*norbs)
        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jt),
     $       2*norbs, v(jvec), 2*norbs, zero, v(jp), 2*norbs)
      endif


C     Transform the density matrix into NO basis:
C       vec' . P* . vec

      if ( inobas .and. itype .eq. 1 ) then
        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jvec), norbs,
     $       v(jpstr), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jvec), norbs, zero, v(jpstr), norbs)

      elseif ( inobas .and. itype .eq. 2 ) then
        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jvec), norbs,
     $       v(jpstr), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jvec), norbs, zero, v(jpstr), norbs)

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jvec+2*nosq),
     $       norbs, v(jpstr+2*nosq), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jvec+2*nosq), norbs, zero, v(jpstr+2*nosq), norbs)

      elseif ( inobas .and. itype .eq. 3 ) then
        call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jvec),
     $       2*norbs, v(jpstr), 2*norbs, zero, v(jt), 2*norbs)
        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jt),
     $       2*norbs, v(jvec), 2*norbs, zero, v(jpstr), 2*norbs)
      endif



C     Save density matrix in the NO basis (in square format) in the
C     array pno.

      if ( itype .eq. 1 ) then
        call amove (2*nosq, v(jp), pno)

      elseif ( itype .eq. 2 ) then
        call amove (4*nosq, v(jp), pno)

      elseif ( itype .eq. 3 ) then
        call amove (8*nosq, v(jp), pno)
      endif

      if ( inobas .and. itype .eq. 1 ) then
        call amove (2*nosq, v(jpstr), pstrno)

      elseif ( inobas .and. itype .eq. 2 ) then
        call amove (4*nosq, v(jpstr), pstrno)

      elseif ( inobas .and. itype .eq. 3 ) then
        call amove (8*nosq, v(jpstr), pstrno)
      endif


C     Print density matrices in the NO basis.

      if ( iprint .gt. 4 ) then
        if ( itype .eq. 1 ) then
          write (iout, *) ' '
          call ltoutx (iout, 1, 2, 1, 0, 'density matrix (NO)', norbs,
     $         norbs, pno, pno, 0)

        elseif ( itype .eq. 2 ) then
          write (iout, *) ' '
          call ltoutx (iout, 1, 3, 1, 0, 'density matrix (NO)', norbs,
     $         norbs, pno, pno(1+2*nosq), 0)

        elseif ( itype .eq. 3 ) then
          write (iout, *) ' '
          call ltoutx (iout, 1, 2, 1, 0, 'density matrix (NO)', 2*norbs,
     $         2*norbs, pno, pno, 0)
        endif
      endif


      return
      end


