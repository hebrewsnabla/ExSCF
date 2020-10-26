

      subroutine phfcyc (iout, iprint, v, lenv, iphf, na, nb, s, isint,
     $     nbasis, norbs, nbas6d, ngrda, ngrdb, ngrdg, p, f, irwp, irwf,
     $     irwvec, irwval, irwnos, inobas, irwciv, irwh, irwx, maxcyc,
     $     itst, acurcy, ifsel, ivshft, vsgoal, diis, mxscfc, idisgo,
     $     istrat, allowp, ipflag, jsym2e, nsymop, nop1, neqbas, neqshl,
     $     neqsh2, rotop, neqatm, accdes, enr, hsp, hph, energy, s2,
     $     convg)

      implicit none


C +------------------------------------------------------------+
C |  phfcyc   ---   CAJH, 07.2011                              |
C |                                                            |
C |                                   (based on PHFB phfb2cy)  |
C |                                                            |
C |                                                            |
C |  Working routine for phfdrv. This routine solves the       |
C |  projected (complex conjugation and / or spin projection)  |
C |  HF equations by repeated diagonalization (SCF) starting   |
C |  from an initial guess prepared in phfdrv. Convergence     |
C |  is reached when the convergence criterion requested is    |
C |  met (see phfdrv for details).                             |
C |                                                            |
C |  The algorithm is simple: given the initial guess of the   |
C |  density matrix, one prepares the effective Fock matrix.   |
C |  Diagonalization of such effective Fock matrix yields      |
C |  new orbitals from which a new density matrix is           |
C |  obtained. This is used to prepare a new effective Fock    |
C |  matrix, and the procedure is repeated. When the density   |
C |  matrix stops changing from iteration to iteration (up     |
C |  to some predetermined threshold), then we declare the     |
C |  equations to be converged.                                |
C |                                                            |
C |  This subroutine works mostly in two basis: the            |
C |  orthonormal AO basis in which the total density matrix    |
C |  is stored, and some form of NO basis (see dmdiag for      |
C |  details) in which the effective Fock matrix is initially  |
C |  prepared.                                                 |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  Details about the strategy followed:                      |
C |                                                            |
C |  1. The standard procedure performs repeated               |
C |     diagonalization without convergence acceleration.      |
C |                                                            |
C |  2. Level-shifting is controlled by ivshft. See phfdrv     |
C |     for details.                                           |
C |                                                            |
C |  3. DIIS (direct inversion of iterative subspace)          |
C |     extrapolation is used if diis is set to true. The      |
C |     number of DIIS vectors is set in mxscfc. The           |
C |     iteration at which DIIS starts in is set in idisgo.    |
C |                                                            |
C |  4. istrat controls the strategy to follow to improve      |
C |     convergence (see phfdrv for details).                  |
C |                                                            |
C |  5. ifsel controls the structure of the Fock matrix to     |
C |     be diagonalized every iteration (see phfdrv for        |
C |     details).                                              |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables (described elsewhere)

      real*8      v(*), p(*), f(*)
      real*8      acurcy, vsgoal, accdes, convg
      real*8      enr, hsp, hph, energy, s2
      integer     iout, iprint, lenv, iphf, na, nb, s
      integer     nbasis, norbs, nbas6d
      integer     ngrda, ngrdb, ngrdg
      integer     irwp, irwf, irwvec, irwval, irwnos, irwciv
      integer     irwh, irwx
      integer     maxcyc, itst
      integer     ifsel, ivshft, mxscfc, idisgo, istrat
      logical     isint, inobas, diis

C       .. prism variables .. [ input variables ]

      real*8      rotop(*)
      integer     ipflag, nsymop, nop1, jsym2e
      integer     neqbas(nbasis,*), neqshl(*), neqsh2(*), neqatm(*)
      logical     allowp(50)

C     DIIS variables
      
      integer     nersav, iencyc, iercyc, iindfk, mxscf1
      integer     nsvdis, iobdis, iopadi, iopbdi, iofadi, iofbdi, ioerdi

      parameter   ( nersav = 2 )

      real*8      errmax(0:nersav)

C     memory allocation variables

      integer     indarr, jvec, jfhf, jpno, jpstno, jx, jend, jend1, mdv
      integer     jgrda, jgrdb, jgrdg
      integer     jwgta, jwgtb, jwgtg
      integer     jdmt, jdind, jyg, jintxg, jsmt, jhmt, jciv

C     scratch RWF files

      integer     irwp1, irwv1, irwfhf, irwng, irwpg, irwgg, irwyg
      integer     irwpa, irwpb, irwfa, irwfb

C     other variables

      real*8      pi, gap, vshft
      real*8      enold, delen, hfen, rmsdp, maxp, thrris
      integer     nel, ntto, ntt2o, nosq, ngrdt, itype, it
      integer     nci, ncis, ncik, szp, szn, szp2, szvec
      integer     igrdb, ind1, ind2
      logical     icmplx, ispuhf, ispghf, done, rises

      parameter   ( thrris = 1.0d-12 )

C     functions

      integer     intowp, igetrw


C     format statements

 1110 format (1X, '***', ' cycle = ', I4, ' ***')
 1120 format (5X, '   E = ', F20.12)
 1121 format (5X, '   E = ', F20.12, 3X, '   DE = ', F20.12)
 1122 format (5X, ' Hsp = ', F20.12, 3X, '  Hph = ', F20.12)
 1123 format (3X, '<S**2> = ', F20.12)
 1124 format (4X, 'RMSDP = ', 2X, 1P, D14.4, 8X, 'MaxP = ', 2X, D14.4)
 1125 format (3X, 'E (HF) = ', F20.12)
 1130 format (5X, 'WARNING: Energy rises - reducing DIIS space')
 1131 format (5X, 'WARNING: Energy rises - increasing level shift')
 1132 format (5X, 'WARNING: Energy rises - continue as usual')
 1200 format (1X, '********************')



C     Useful quantities.

      nel = na + nb

      ntto  = norbs * (norbs + 1) / 2
      ntt2o = norbs * (2*norbs + 1)
      nosq  = norbs * norbs

      ngrdt = ngrda * ngrdb * ngrdg


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
        call gauerr ('Unrecognized iphf in phfcyc.')
      endif


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


C     Determine sizes for P(g) and N(g) matrices (szp, szn).
C     Determine size for total density matrix P (szp2).
C     Determine size for orbital coefficients (szvec).

      if ( itype .eq. 1 ) then
        szp   = 2*nosq
        szn   = 2*na*na
        szp2  = 2*ntto
        szvec = 2*nosq

      elseif ( itype .eq. 2 ) then
        szp   = 4*nosq
        szn   = 2*na*na + 2*nb*nb
        szp2  = 4*ntto
        szvec = 4*nosq

      elseif ( itype .eq. 3 ) then
        szp  = 8*nosq
        szn  = 2*nel*nel

        if ( .not. ispuhf ) then
          szp2  = 2*ntt2o
          szvec = 8*nosq
        else
          szp2  = 4*ntto
          szvec = 4*nosq
        endif
      endif


C     Memory allocation.

C       Allocate space for integration grid.

      indarr = 1
      jgrda  = indarr + intowp (3*ngrdt)
      jwgta  = jgrda  + ngrda
      jgrdb  = jwgta  + ngrda
      jwgtb  = jgrdb  + ngrdb
      jgrdg  = jwgtb  + ngrdb
      jwgtg  = jgrdg  + ngrdg
      jend   = jwgtg  + ngrdg

C       Allocate space for Wigner's small d-matrix and vector with
C       indices for full Wigner rotation matrices.

      jdmt  = jend
      jdind = jdmt  + ngrdb * ncis * ncis
      jend  = jdind + ncis

C       Allocate space for yg vector, the overlap and Hamiltonian
C       matrices associated with the CI problem, and the CI
C       coefficients.

c$$$      jyg  = jend
c$$$      jsmt = jyg  + 2*ngrdt * nci*nci
      jyg  = jend
      jsmt = jyg  + 2*ngrdt * ncik*ncik
      jhmt = jsmt + 2*nci*nci
      jciv = jhmt + 2*nci*nci
      jend = jciv + 2*nci

C       Allocate space for integrated matrix intxg (see formyg).

      jintxg = jend
      jend   = jintxg + szp

C       Allocate space for set of orbital coefficients.
C       Allocate space for Fock matrix of HF state.

      jvec = jend
      jfhf = jvec + szvec
      jend = jfhf + szp2

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


C       Allocate space for DIIS buffers.

      if ( diis ) then
        mxscf1 = mxscfc + 1
      else
        mxscf1 = 0
      endif

      iencyc = jend
      iercyc = iencyc + mxscf1
      iindfk = iercyc + mxscf1
      jend   = iindfk + intowp (mxscf1)


C       Check whether there is enough memory available.

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'phfcyc')



C     Define scratch file for old density matrix and old orbitals.

      irwp1 = igetrw (szp2, 0)
      irwv1 = igetrw (szvec, 0)

C     Define scratch file for HF Fock matrix.

      irwfhf = igetrw (szp2, 0)


C     Define scratch files needed for DIIS.

      if ( ( itype .eq. 2 .or. ispuhf ) .and. diis ) then
        irwpa = igetrw (2*ntto, 0)
        irwpb = igetrw (2*ntto, 0)
        irwfa = igetrw (2*ntto, 0)
        irwfb = igetrw (2*ntto, 0)

      else
        irwpa = irwp
        irwpb = 0
        irwfa = irwf
        irwfb = 0
      endif


C     Define scratch array for matrices:
C       - irwng - N(g) matrices
C       - irwpg - P(g) matrices
C       - irwgg - G(g) matrices
C       - irwyg - Y(g) matrices

c$$$      irwng = igetrw (nci*nci * szn * ngrdt, 0)
c$$$      irwpg = igetrw (nci*nci * szp * ngrdt, 0)
c$$$      irwgg = igetrw (nci*nci * szp * ngrdt, 0)
c$$$      irwyg = igetrw (nci*nci * szp * ngrdt, 0)
      irwng = igetrw (ncik*ncik * szn * ngrdt, 0)
      irwpg = igetrw (ncik*ncik * szp * ngrdt, 0)
      irwgg = igetrw (ncik*ncik * szp * ngrdt, 0)
      irwyg = igetrw (ncik*ncik * szp * ngrdt, 0)



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Prepare integration grid  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
        write(iout, 1000) (v(jgrdb + igrdb) , igrdb=0,7)
        write(iout, 1001) (v(jwgtb + igrdb) , igrdb=0,7)
 1000   Format('beta in GauLeg', / , 'grd=', 8(F10.4))
 1001   Format('wgt=', 8(F10.4))
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


C     Prepare arrays associated with Wigner's rotation matrix.

      call wigarr (iout, iprint, s, isint, na, nb, ncis, ngrdb,
     $     v(jgrdb), v(jdmt), v(jdind))



C     Initialize DIIS variables.

      nsvdis = 0
      iobdis = 0
      iopadi = 0
      iopbdi = 0
      iofadi = 0
      iofbdi = 0
      ioerdi = 0


C     Initialize some variables.
C     Set convergence to false.

      gap = 0.0d0
      vshft = 0.0d0

      done = .false.



C     %%%%%%%%%%%%%%%%%%%%
C     %  SCF iterations  %
C     %%%%%%%%%%%%%%%%%%%%


      do 100 it = 1, maxcyc

        call tstamp (1, 'Start of iteration.')


C       Load and print total density matrix.

        call fileio (2, -irwp, szp2, p, 0)

        if ( iprint .gt. 3 ) then
          if ( itype .eq. 1 ) then
            write (iout, *) ' '
            call ltoutx (iout, 0, 2, 1, 0, 'HF density matrix', norbs,
     $           norbs, p, p, 0)

          elseif ( itype .eq. 2 .or. ispuhf ) then
            write (iout, *) ' '
            call ltoutx (iout, 0, 3, 1, 0, 'HF density matrix', norbs,
     $           norbs, p, p(1+2*ntto), 0)

          elseif ( itype .eq. 3 ) then
            write (iout, *) ' '
            call ltoutx (iout, 0, 2, 1, 0, 'HF density matrix', 2*norbs,
     $           2*norbs, p, p, 0)
          endif
        endif


C       Build Fock matrix associated with the HF state.
        write(iout, *) '\nStart hffock'

        call hffock (iout, iprint, v(jend), mdv, iphf, norbs, nbasis,
     $       nbas6d, irwp, irwfhf, v(jfhf), irwh, irwx, hfen, allowp,
     $       ipflag, jsym2e, nsymop, nop1, neqbas, neqshl, neqsh2,
     $       rotop, neqatm, accdes)

        hfen = enr + hfen



C       %%%%%%%%%%%%%%%%%%%
C       %  Find NO basis  %
C       %%%%%%%%%%%%%%%%%%%


C       Diagonalize the density matrix to find the NO basis (or some
C       other orthonormal basis).

C       The density matrix in such basis is stored in jpno.

        write(iout, *) '\nStart hfdmdg: find NO basis'
        call hfdmdg (iout, iprint, v(jend), mdv, iphf, na, nb, norbs,
     $       irwnos, p, v(jpno), v(jpstno), .true., inobas)

C        write(iout, *) 'End hfdmdg'


C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C       %  Build transition density matrices  %
C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        write(iout, *) '\nStart formng/pg'

C       Form overlap-related matrices N(g).

        call formng (iout, iprint, v(jend), mdv, iphf, ncik, norbs,
     $       na, nb, ngrdt, v(indarr), v(jgrda), v(jgrdb), v(jgrdg),
     $       irwnos, irwng, v(jpno), v(jpstno), inobas)


C       Build transition density matrices P(g).

        call formpg (iout, iprint, v(jend), mdv, iphf, ncik, norbs,
     $       na, nb, ngrdt, v(indarr), v(jgrda), v(jgrdb), v(jgrdg),
     $       irwnos, irwng, irwpg, v(jpno), v(jpstno), inobas)



C       %%%%%%%%%%%%%%%%
C       %  Build G(g)  %
C       %%%%%%%%%%%%%%%%


C       Contract P(g) with electron repulsion integrals to form the
C       matrices G(g).
        write(iout, *) '\nStart formgg'

        call formgg (iout, iprint, v(jend), mdv, iphf, ncik, norbs,
     $       nbasis, nbas6d, ngrdt, irwpg, irwgg, irwx, irwnos, allowp,
     $       ipflag, jsym2e, nsymop, nop1, neqbas, neqshl, neqsh2,
     $       rotop, neqatm, accdes)



C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C       %  Build CI matrix elements  %
C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C       Evaluate the coefficients x(g) and build Smt.
        write(iout, *) '\nStart evalsmt/hmt'
c        write(iout, 1001) (v(jwgtb + igrdb) , igrdb=0,7)

        call evalsmt (iout, iprint, v(jend), mdv, iphf, nci, ncis, ncik, 
     $       ngrdb, norbs, na, nb, irwnos, inobas, ngrdt, v(indarr),
     $       v(jgrda), v(jwgta), v(jgrdb), v(jwgtb), v(jgrdg), v(jwgtg),
     $       v(jdmt), v(jdind), irwvec, irwng, v(jyg), v(jsmt))


C       Evaluate the Hamiltonian matrix elements and build Hmt.

        call evalhmt (iout, iprint, v(jend), mdv, iphf, nci, ncis, ncik,
     $       ngrdb, nbasis, norbs, irwh, irwx, irwnos, ngrdt, v(indarr),
     $       v(jgrda), v(jwgta), v(jgrdb), v(jwgtb), v(jgrdg), v(jwgtg),
     $       v(jdmt), v(jdind), v(jyg), irwgg, irwpg, v(jhmt))


c        write(iout, *) 'End build CI'

C       %%%%%%%%%%%%%%%%%%%%%%
C       %  Solve CI problem  %
C       %%%%%%%%%%%%%%%%%%%%%%


C       Solve CI problem associated with spin projection and / or
C       complex conjugation restoration.

        call solvci (iout, iprint, v(jend), mdv, nci, v(jhmt), v(jsmt),
     $       v(jciv))


C       Evaluate the coefficients y(g) using x(g) and the CI
C       coefficients obtained above.

C       CAJH 12.2012 - subroutine evalyg has been removed.

c$$$
c$$$        call evalyg (iout, iprint, nci, ngrdt, v(jyg), v(jsmt), v(jciv))



C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C       %  Evaluate spin expectation values  %
C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C       Evaluate the expectation value of S^2.

        call evals2 (iout, iprint, v(jend), mdv, iphf, nci, ncis, ncik,
     $       ngrdb, norbs, ngrdt, v(indarr), v(jgrda), v(jwgta),
     $       v(jgrdb), v(jwgtb), v(jgrdg), v(jwgtg), v(jdmt), v(jdind),
     $       v(jyg), irwnos, irwpg, v(jciv), s2)



C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C       %  Construct effective Fock matrix  %
C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C       Form the matrices Y(g) (related to derivatives of the weighting
C       function y(g) with respect to elements of the density matrix).

        call formyg (iout, iprint, v(jend), mdv, iphf, nci, ncis, ncik,
     $       ngrdb, norbs, na, nb, ngrdt, v(indarr), v(jgrda), v(jwgta),
     $       v(jgrdb), v(jwgtb), v(jgrdg), v(jwgtg), v(jdmt), v(jdind),
     $       v(jciv), irwnos, irwng, irwyg, v(jintxg), v(jyg), v(jpno),
     $       v(jpstno), inobas)


C       Build the effective Fock matrix and store it in appropriate RWF
C       value. formfk also computes the single-particle (Hsp) and the
C       particle-hole energies (Hph).

        call formfk (iout, iprint, v(jend), mdv, iphf, nci, ncis, ncik,
     $       ngrdb, nbasis, norbs, na, nb, irwh, irwx, irwnos, ngrdt,
     $       v(indarr), v(jgrda), v(jwgta), v(jgrdb), v(jwgtb),
     $       v(jgrdg), v(jwgtg), v(jdmt), v(jdind), v(jciv), irwng,
     $       irwpg, irwgg, irwyg, v(jintxg), v(jyg), v(jpno), v(jpstno),
     $       inobas, hsp, hph, f, irwf)

        call focksel (iout, iprint, v(jend), mdv, iphf, na, nb, norbs,
     $       irwf, irwfhf, irwnos, f, p, ivshft, vsgoal, vshft, gap,
     $       ifsel)



C       %%%%%%%%%%%%%%%%%%%%%%%%%%
C       %  Compute total energy  %
C       %%%%%%%%%%%%%%%%%%%%%%%%%%


        if ( it .gt. 1 ) then
          enold = energy
        else
          enold = 0.0d0
        endif

        energy = enr + hsp + hph
        delen  = energy - enold


C       Increase level shift if energy is above previous energy and diis
C       is turned off (if requested).

        rises = delen .gt. thrris

        if ( rises .and. (.not. diis) .and. (istrat .eq. 2) ) then
          ivshft = ivshft + ivshft / 2
          vsgoal = 1.5d0 * vsgoal
        endif



C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C       %  Diagonalize effective Fock matrix  %
C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C       Save old density matrix.

        call fileio (1, -irwp1, szp2, p, 0)


C       Apply direct inversion of iterative subspace (DIIS)
C       extrapolation if requested.

        if ( diis ) then

C          Load PA, PB and FA, FB matrices if itype = 2.

           if ( itype .eq. 2 .or. ispuhf ) then
             jx    = jend
             jend1 = jx   + szp2

             call tstcor (jend1-1, lenv, 'phfcyc')

             call fileio (2, -irwp, szp2, v(jx), 0)
             call fileio (1, -irwpa, 2*ntto, v(jx), 0)
             call fileio (1, -irwpb, 2*ntto, v(jx+2*ntto), 0)

             call fileio (2, -irwf, szp2, v(jx), 0)
             call fileio (1, -irwfa, 2*ntto, v(jx), 0)
             call fileio (1, -irwfb, 2*ntto, v(jx+2*ntto), 0)
           endif


           call phfdis (iout, iprint, v(jend), mdv, iphf, norbs, it,
     $          idisgo, istrat, rises, mxscfc, irwpa, irwpb, irwfa,
     $          irwfb, errmax, nersav, v(iindfk), energy, v(iencyc),
     $          v(iercyc), nsvdis, iobdis, iopadi, iopbdi, iofadi,
     $          iofbdi, ioerdi)


C          Recover FA, FB matrices if itype = 2.

           if ( itype .eq. 2 .or. ispuhf ) then
             jx    = jend
             jend1 = jx   + szp2

             call tstcor (jend1-1, lenv, 'phfcyc')

             call fileio (2, -irwfa, 2*ntto, v(jx), 0)
             call fileio (2, -irwfb, 2*ntto, v(jx+2*ntto), 0)
             call fileio (1, -irwf, szp2, v(jx), 0)
           endif

        endif


C       Save old set of orbital coefficients.

        call fileio (2, -irwvec, szvec, v(jvec), 0)
        call fileio (1, -irwv1, szvec, v(jvec), 0)


C       Diagonalize the effective Fock matrix to get a new set of
C       molecular orbitals. Build the density matrix associated with
C       such orbitals.

C         Here, level shifting is applied to the Fock matrix if it was
C         requested before diagonalization takes place.

        call fileio (2, -irwf, szp2, f, 0)

        call fockdg (iout, iprint, v(jend), mdv, iphf, norbs, na, nb,
     $       ivshft, vshft, gap, f, p, irwvec, irwval, irwp)



C       %%%%%%%%%%%%%%%%%%%%%%%%%%
C       %  Test for convergence  %
C       %%%%%%%%%%%%%%%%%%%%%%%%%%


C       Test for convergence of the energy or the elements of the
C       density matrix.

        call cnvtst (iout, iprint, v(jend), mdv, iphf, itst, acurcy,
     $       norbs, irwp, irwp1, rmsdp, maxp, delen, done)



C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C       %  Print out cycle-related information  %
C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ( iprint .gt. 0 ) then
          write (iout, *)
          write (iout, 1110) it
          write (iout, *)

          if ( it .eq. 1 ) then
            write (iout, 1120) energy
          else
            write (iout, 1121) energy, delen
          endif

          write (iout, 1122) hsp, hph
          write (iout, 1123) s2
          write (iout, 1124) rmsdp, maxp

          write (iout, *)
          write (iout, 1125) hfen

          if ( rises ) then
            write (iout, *)

            if ( diis .and. istrat .eq. 2 ) then
              write (iout, 1130)
            elseif ( ivshft .ne. -1 .and. istrat .eq. 2 ) then
              write (iout, 1131)
            else
              write (iout, 1132)
            endif
          endif

          write (iout, *)
          write (iout, 1200)
          write (iout, *)
        endif


        call tstamp (1, 'End of iteration.')

C       Exit loop if we are done.

        if ( done ) goto 101

 100  continue
 101  continue


C     Go back to old density matrix and old set of orbital coefficients.

      call fileio (2, -irwp1, szp2, p, 0)
      call fileio (1, -irwp, szp2, p, 0)

      call fileio (2, -irwv1, szvec, v(jvec), 0)
      call fileio (1, -irwvec, szvec, v(jvec), 0)


C     Return some things to phfdrv.

      if ( it .lt. maxcyc ) maxcyc = it

      if ( itst .eq. 1 ) then
        convg = delen
      elseif ( itst .eq. 2 ) then
        convg = rmsdp
      elseif ( itst .eq. 3 ) then
        convg = maxp
      endif


C     Save CI vector.

      call conddf (irwciv, 2*nci)
      call fileio (1, -irwciv, 2*nci, v(jciv), 0)


C     Free scratch RWF files.

      call freerw (irwp1)
      call freerw (irwv1)

      if ( ( itype .eq. 2 .or. ispuhf ) .and. diis ) then
        call freerw (irwpa)
        call freerw (irwpb)
        call freerw (irwfa)
        call freerw (irwfb)
      endif

      call freerw (irwng)
      call freerw (irwpg)
      call freerw (irwgg)
      call freerw (irwyg)


      return
      end


