

      subroutine phfdm (iout, iprint, v, lenv, iphf, na, nb, s, isint,
     $     nbasis, norbs, ngrda, ngrdb, ngrdg, p, irwp, irwvec, irwnos,
     $     inobas, irwciv, irwh, irwx, irw1dm, irw2dm, savedm, prtdm,
     $     dotwom, useao)

      implicit none


C +------------------------------------------------------------+
C |  phfdm   ---   CAJH, 08.2011                               |
C |                                                            |
C |                                                            |
C |  Driver routine for the evaluation of density matrices     |
C |  associated with the PHF state. This is independent of     |
C |  phfcyc because double integration over the Euler angles   |
C |  is required in the case of spin projection.               |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  CAJH, 12.2012                                             |
C |                                                            |
C |  Switched to single integration.                           |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables (described elsewhere)

      real*8      v(*), p(*)
      integer     iout, iprint, lenv, iphf, na, nb, s
      integer     nbasis, norbs
      integer     ngrda, ngrdb, ngrdg
      integer     irwp, irwvec, irwnos, irwciv, irwh, irwx
      integer     irw1dm, irw2dm
      logical     isint, inobas, savedm, prtdm, dotwom, useao

C     memory allocation variables

      integer     indarr, jpno, jpstno, jend, mdv
      integer     jgrda, jgrdb, jgrdg
      integer     jwgta, jwgtb, jwgtg
      integer     jciv

C     other variables

      real*8      pi
      integer     nel, ntto, ntt2o, nosq, ngrdt, itype
      integer     nci, ncis, ncik, szp2
      integer     igrdb, ind1, ind2
      logical     icmplx, ispuhf, ispghf

C     functions

      integer     intowp, igetrw



C     Useful quantities.

      nel = na + nb

      ntto  = norbs * (norbs + 1) / 2
      ntt2o = norbs * (2*norbs + 1)
      nosq  = norbs * norbs


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
        call gauerr ('Unrecognized iphf in phfdm.')
      endif


C     Adjust number of grid points for ispuhf.
C     Determine total number of grid points.

c$$$      if ( ispuhf ) then
c$$$        ngrda = ngrdb + 1
c$$$      endif

c$$$      ngrdt = ngrda**2 * ngrdb**2 * ngrdg**2

      ngrdt = ngrda*ngrdb*ngrdg


C     Determine size of CI problem.

C       .. number of CI vectors for spin projection (>1 only if GHF)

      if ( ispghf ) then
        if ( isint ) then
          ncis = 2*s + 1
        else
          ncis = s + 1
        endif
      else
        ncis = 1
      endif

      nci = ncis

C       .. double the number of CI states if complex conjugation
C          restoration will be performed

      ncik = 1
      if ( icmplx ) ncik = 2

      nci = nci*ncik


C     Determine size for total density matrix P (szp2).

      if ( itype .eq. 1 ) then
        szp2 = 2*ntto

      elseif ( itype .eq. 2 ) then
        szp2 = 4*ntto

      elseif ( itype .eq. 3 ) then
        if ( .not. ispuhf ) szp2 = 2*ntt2o
        if ( ispuhf ) szp2 = 4*ntto
      endif


C     Memory allocation.

C       Allocate space for integration grid.

      indarr = 1
      jgrda  = indarr + intowp (6*ngrdt)
      jwgta  = jgrda  + ngrda
      jgrdb  = jwgta  + ngrda
      jwgtb  = jgrdb  + ngrdb
      jgrdg  = jwgtb  + ngrdb
      jwgtg  = jgrdg  + ngrdg
      jend   = jwgtg  + ngrdg

C       Allocate space for vector of CI coefficients.

      jciv = jend
      jend = jciv + 2*nci

C       Allocate space for density matrix in NO basis.

      if ( itype .eq. 1 ) then
        jpno = jend
        jend = jpno + 2*nosq

      elseif ( itype .eq. 2 ) then
        jpno = jend
        jend = jpno + 4*nosq

      elseif ( itype .eq. 3 ) then
        jpno = jend
        jend = jpno + 8*nosq
      endif

      jpstno = jend

      if ( inobas .and. itype .eq. 1 ) then
        jend = jpstno + 2*nosq
      elseif ( inobas .and. itype .eq. 2 ) then
        jend = jpstno + 4*nosq
      elseif ( inobas .and. itype .eq. 3 ) then
        jend = jpstno + 8*nosq
      endif


C       Check whether there is enough memory available.

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'phfcyc')



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Prepare integration grid  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Prepare grid associated with spin projection.

      call indgrd (iout, iprint, v(indarr), ngrdt, ngrda, ngrdb, ngrdg,
     $     1)

      pi = 4.0d0 * atan (1.0d0)

      if ( ngrda .gt. 1 ) then
        call trpgrd (0.0d0, 2*pi, v(jgrda), v(jwgta), ngrda)
      else
        v(jgrda) = 0.0d0
        v(jwgta) = 1.0d0
      endif

      if ( ngrdb .gt. 1 ) then
        call gauleg (0.0d0, pi, v(jgrdb), v(jwgtb), ngrdb)

        do 10 igrdb = 1, ngrdb
          ind1 = jgrdb + igrdb - 1
          ind2 = jwgtb + igrdb - 1
          v(ind2) = v(ind2) * sin (v(ind1))
 10     continue
      else
        v(jgrdb) = 0.0d0
        v(jwgtb) = 1.0d0
      endif

      if ( ngrdg .gt. 1 ) then
        call trpgrd (0.0d0, 2*pi, v(jgrdg), v(jwgtg), ngrdg)
      else
        v(jgrdg) = 0.0d0
        v(jwgtg) = 1.0d0
      endif


C     Load CI vector.

      call fileio (2, -irwciv, 2*nci, v(jciv), 0)



C     %%%%%%%%%%%%%%%%%%%
C     %  Find NO basis  %
C     %%%%%%%%%%%%%%%%%%%


C     Load HF density matrix.

      call fileio (2, -irwp, szp2, p, 0)


C     Diagonalize the density matrix to find the NO basis (or some other
C     orthonormal basis).

C     The density matrix in such basis is stored in jpno.

C       NOTE: This has been changed. Now we use the eigenvectors already
C             stored in irwnos. This avoids the need to recompute the
C             CI coefficients.

      call hfdmdg (iout, iprint, v(jend), mdv, iphf, na, nb, norbs,
     $     irwnos, p, v(jpno), v(jpstno), .false., inobas)



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Build PHF density matrices  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Build 1PDM and possibly 2PDM associated with the PHF state.
C     Diagonalize to recover the natural orbital occupations.


      call phfdm1 (iout, iprint, v(jend), mdv, iphf, nci, ncis, ncik,
     $     nbasis, norbs, na, nb, s, isint, irwh, irwx, irwvec, irwnos,
     $     ngrdt, v(indarr), v(jgrda), v(jwgta), v(jgrdb), v(jwgtb),
     $     v(jgrdg), v(jwgtg), irw1dm, irw2dm, v(jpno), v(jpstno),
     $     v(jciv), inobas, savedm, prtdm, dotwom, useao)


      return
      end


