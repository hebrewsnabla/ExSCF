

      subroutine focksel (iout, iprint, v, lenv, iphf, na, nb, norbs,
     $     irwf, irwfhf, irwnos, f, p, ivshft, vsgoal, vshft, gap,
     $     ifsel)

      implicit none


C +------------------------------------------------------------+
C |  focksel   --   CAJH, 08.2011                              |
C |                                                            |
C |                                                            |
C |  Construct the Fock matrix to be diagonalized according    |
C |  to two different methods:                                 |
C |                                                            |
C |    ifsel = 1,  the Fock matrix used corresponds to the     |
C |                effective Fock matrix of PHF                |
C |                                                            |
C |          = 2,  the occupied-occupied and virtual-virtual   |
C |                block of the Fock matrix are taken from     |
C |                the regular HF Fock matrix; the ov and vo   |
C |                blocks correspond to the effective Fock     |
C |                matrix of PHF                               |
C |                                                            |
C |  Level-shifting is also applied to the Fock matrix to      |
C |  be diagonalized (see below).                              |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  Level-shifting.                                           |
C |                                                            |
C |  Level-shifting is used to shift the occupied orbital      |
C |  energies down (or the virtuals up) in order to reduce     |
C |  occupied-virtual mixing. This improves convergence in     |
C |  most cases.                                               |
C |                                                            |
C |  Specifically, we perform level shifting by letting        |
C |                                                            |
C |    F = F - vshft * P                                       |
C |                                                            |
C |  where F is the effective Fock matrix, P is the HF         |
C |  density matrix, and vshft is a constant (positive) shift  |
C |  applied.                                                  |
C |                                                            |
C |  Note that, at convergence, F and P commute, and thus      |
C |  the level-shifting truly only shifts the occupied         |
C |  energies down, without modifying the eigenvectors.        |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       na     - number of alpha electrons
C       nb     - number of beta electrons
C       norbs  - number of orbitals
C       irwf   - RWF file with effective Fock matrix  [ updated here ]
C       irwfhf - Fock matrix associated with HF state
C       irwnos - RWF file with transformation matrix from orthonormal AO
C                to NO basis
C       f      - effective Fock matrix  [ in / out ]
C       p      - HF density matrix
C       ivshft - level shifting control
C       vsgoal - goal gap in dynamic level shifting
C       vshft  - level shifting applied  [ out ]
C       gap    - HOMO-LUMO difference in Fock eigenvalues
C                (without level shifting)
C       ifsel  - determines how to construct Fock matrix (see above)

      real*8      vsgoal, vshft, gap
      integer     iout, iprint, lenv, iphf, na, nb, norbs, ivshft
      integer     irwf, irwfhf, irwnos, ifsel

      real*8      v(*), f(*), p(*)

C     other variables

      integer     ntto, ntt2o, nosq, nel, itype
      integer     jnos, jnost
      integer     jp, jq, jt, jf, jf2, jfhf
      integer     jfoo, jfov, jfvo, jfvv
      integer     jz1, jz2, jz3
      integer     jzaa, jzab, jzba, jzbb, jscr, lscr
      integer     mdv, jend
      integer     k1, k2, ind1, ind2
      logical     icmplx, ispuhf, ispghf

C     constants

      complex*16  zero, one

C     useproj - whether to use projection operators for ifsel = 2
C               or do everything manually

      logical     useproj
      save        useproj
      data        useproj / .true. /



      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)


C     Construct useful variables.

      ntto = norbs * (norbs + 1) / 2
      ntt2o = norbs * (2*norbs + 1)

      nosq = norbs * norbs

      nel = na + nb


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
        call gauerr ('Unrecognized iphf in focksel.')
      endif


C     Split the code according to the parameter selected.
      write(iout, '("ifsel = ", I5, ", useproj", L5)') ifsel, useproj
      if ( ifsel .eq. 1 ) goto 300
      if ( ifsel .eq. 2 .and. useproj ) goto 200
      if ( ifsel .eq. 2 .and. .not. useproj ) goto 100



C     %%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Modified Fock matrix  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%


C     Prepare a modified version of the Fock matrix, where the blocks
C     are constructed as:

C       occupied - occupied block   --  HF Fock matrix
C       occupied - virtual  block   --  PHF effective Fock matrix
C       virtual  - occupied block   --  PHF effective Fock matrix
C       virtual  - virtual  block   --  HF Fock matrix

C     We use projection operators in order to extract the appropriate
C     blocks of the HF or PHF Fock matrices.


C     ** Note that for icmplx = .true., the definition of occupied or
C        virtual is determined from the real part of the density matrix,
C        as opposed to the full matrix. See hfdmdg for more details.



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Manually build Fock matrix  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 100  continue


C     Memory allocation.

C       Allocate space for:
C         - orthonormal AO => NO basis transformation matrix (jnos)
C         - scratch array to perform basis transformation (jnost)
C         - scratch copy of PHF Fock matrix (jf)
C         - scratch copy of Fock matrix for LT transformation (jf2)
C         - scratch copy of HF Fock matrix (jfhf)

      if ( itype .eq. 1 ) then
        jnos  = 1
        jnost = jnos  + 2*nosq
        jf    = jnost + 2*nosq
        jfhf  = jf    + 2*nosq
        jend  = jfhf  + 2*nosq

      elseif ( itype .eq. 2 ) then
        jnos  = 1
        jnost = jnos  + 4*nosq
        jf    = jnost + 2*nosq
        jfhf  = jf    + 4*nosq
        jend  = jfhf  + 4*nosq

      elseif ( itype .eq. 3 ) then
        jnos  = 1
        jnost = jnos  + 8*nosq
        jf    = jnost + 8*nosq
        jfhf  = jf    + 8*nosq
        jend  = jfhf  + 8*nosq
      endif

      if ( itype .lt. 3 .or. ispuhf ) then
        jf2  = jend
        jend = jf2  + 2*nosq

      elseif ( itype .eq. 3 ) then
        jf2  = jend
        jend = jf2  + 8*nosq
      endif

C       Allocate some extra scratch arrays for ispuhf.

      if ( ispuhf ) then
        lscr  = 16*nosq

        jzaa = jend
        jzab = jzaa + 2*nosq
        jzba = jzab + 2*nosq
        jzbb = jzba + 2*nosq
        jscr = jzbb + 2*nosq
        jend = jscr + lscr
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

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'focksel')


C     Load transformation matrix (orthonormal AO => NO basis).

      if ( itype .eq. 1 ) then
        call fileio (2, -irwnos, 2*nosq, v(jnos), 0)

      elseif ( itype .eq. 2 ) then
        call fileio (2, -irwnos, 4*nosq, v(jnos), 0)

      elseif ( itype .eq. 3 ) then
        call fileio (2, -irwnos, 8*nosq, v(jnos), 0)
      endif


C     Load PHF and HF Fock matrices from RWF files.

      if ( itype .eq. 1 ) then
        call fileio (2, -irwf, 2*ntto, v(jf), 0)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call fileio (2, -irwf, 2*ntto, v(jf), 0)
        call fileio (2, -irwf, 2*ntto, v(jf+2*nosq), 2*ntto)

      elseif ( itype .eq. 3 ) then
        call fileio (2, -irwf, 2*ntt2o, v(jf), 0)
      endif

      if ( itype .eq. 1 ) then
        call fileio (2, -irwfhf, 2*ntto, v(jfhf), 0)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call fileio (2, -irwfhf, 2*ntto, v(jfhf), 0)
        call fileio (2, -irwfhf, 2*ntto, v(jfhf+2*nosq), 2*ntto)

      elseif ( itype .eq. 3 ) then
        call fileio (2, -irwfhf, 2*ntt2o, v(jfhf), 0)
      endif



C     Transform the PHF Fock matrix (orthonormal AO basis) into square
C     format.

      if ( itype .eq. 1 ) then
        call zrecmp (1, v(jf), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jf), v(jz2), nosq)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call zrecmp (1, v(jf), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jf), v(jz2), nosq)

        call zrecmp (1, v(jf+2*nosq), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jf+2*nosq), v(jz2), nosq)

      elseif ( itype .eq. 3 ) then
        call zrecmp (1, v(jf), v(jz1), ntt2o)
        call square (v(jz1), v(jz2), 2*norbs, 2*norbs, 0)
        call square (v(jz1+ntt2o), v(jz3), 2*norbs, 2*norbs, -1)
        call zrecmp (2, v(jf), v(jz2), 4*nosq)
      endif


C     Transform the HF Fock matrix (orthonormal AO basis) into square
C     format.

      if ( itype .eq. 1 ) then
        call zrecmp (1, v(jfhf), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jfhf), v(jz2), nosq)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call zrecmp (1, v(jfhf), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jfhf), v(jz2), nosq)

        call zrecmp (1, v(jfhf+2*nosq), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jfhf+2*nosq), v(jz2), nosq)

      elseif ( itype .eq. 3 ) then
        call zrecmp (1, v(jfhf), v(jz1), ntt2o)
        call square (v(jz1), v(jz2), 2*norbs, 2*norbs, 0)
        call square (v(jz1+ntt2o), v(jz3), 2*norbs, 2*norbs, -1)
        call zrecmp (2, v(jfhf), v(jz2), 4*nosq)
      endif


C     Build full spin-orbital matrices for ispuhf = .true.

      if ( ispuhf ) then
        call amove (2*nosq, v(jf), v(jzaa))
        call amove (2*nosq, v(jf+2*nosq), v(jzbb))

        call aclear (2*nosq, v(jzab))
        call aclear (2*nosq, v(jzba))

        call dmblck (iout, iprint, v(jscr), lscr, norbs, v(jf),
     $       v(jzaa), v(jzab), v(jzba), v(jzbb), 1, -1)

        call amove (2*nosq, v(jfhf), v(jzaa))
        call amove (2*nosq, v(jfhf+2*nosq), v(jzbb))

        call aclear (2*nosq, v(jzab))
        call aclear (2*nosq, v(jzba))

        call dmblck (iout, iprint, v(jscr), lscr, norbs, v(jfhf),
     $       v(jzaa), v(jzab), v(jzba), v(jzbb), 1, -1)
      endif


C     Transform Fock matrices into NO basis:
C       nos! . F . nos

      if ( itype .eq. 1 ) then

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jf), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos), norbs, zero, v(jf), norbs)

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jfhf), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos), norbs, zero, v(jfhf), norbs)

      elseif ( itype .eq. 2 ) then

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jf), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos), norbs, zero, v(jf), norbs)

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos+2*nosq),
     $       norbs, v(jf+2*nosq), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos+2*nosq), norbs, zero, v(jf+2*nosq), norbs)

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jfhf), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos), norbs, zero, v(jfhf), norbs)

        call zgemm ('c', 'n', norbs, norbs, norbs, one, v(jnos+2*nosq),
     $       norbs, v(jfhf+2*nosq), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos+2*nosq), norbs, zero, v(jfhf+2*nosq), norbs)

      elseif ( itype .eq. 3 ) then

        call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jnos),
     $       2*norbs, v(jf), 2*norbs, zero, v(jnost), 2*norbs)
        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jnost),
     $       2*norbs, v(jnos), 2*norbs, zero, v(jf), 2*norbs)

        call zgemm ('c', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jnos),
     $       2*norbs, v(jfhf), 2*norbs, zero, v(jnost), 2*norbs)
        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jnost),
     $       2*norbs, v(jnos), 2*norbs, zero, v(jfhf), 2*norbs)
      endif


C     Use oo and vv blocks from regular HF Fock matrix in effective Fock
C     matrix.

      if ( itype .eq. 1 ) then

        do 1101 k1 = 1, na
          do 1102 k2 = 1, na
            ind1 = (k1-1)*2*norbs + (k2-1)*2 + jf
            ind2 = (k1-1)*2*norbs + (k2-1)*2 + jfhf

            v(ind1) = v(ind2)
            v(ind1+1) = v(ind2+1)
 1102     continue
 1101   continue

        do 1103 k1 = na+1, norbs
          do 1104 k2 = na+1, norbs
            ind1 = (k1-1)*2*norbs + (k2-1)*2 + jf
            ind2 = (k1-1)*2*norbs + (k2-1)*2 + jfhf

            v(ind1) = v(ind2)
            v(ind1+1) = v(ind2+1)
 1104     continue
 1103   continue

      elseif ( itype .eq. 2 ) then

        do 1201 k1 = 1, na
          do 1202 k2 = 1, na
            ind1 = (k1-1)*2*norbs + (k2-1)*2 + jf
            ind2 = (k1-1)*2*norbs + (k2-1)*2 + jfhf

            v(ind1) = v(ind2)
            v(ind1+1) = v(ind2+1)
 1202     continue
 1201   continue

        do 1203 k1 = na+1, norbs
          do 1204 k2 = na+1, norbs
            ind1 = (k1-1)*2*norbs + (k2-1)*2 + jf
            ind2 = (k1-1)*2*norbs + (k2-1)*2 + jfhf

            v(ind1) = v(ind2)
            v(ind1+1) = v(ind2+1)
 1204     continue
 1203   continue

        do 1205 k1 = 1, nb
          do 1206 k2 = 1, nb
            ind1 = (k1-1)*2*norbs + (k2-1)*2 + jf + 2*nosq
            ind2 = (k1-1)*2*norbs + (k2-1)*2 + jfhf + 2*nosq

            v(ind1) = v(ind2)
            v(ind1+1) = v(ind2+1)
 1206     continue
 1205   continue

        do 1207 k1 = nb+1, norbs
          do 1208 k2 = nb+1, norbs
            ind1 = (k1-1)*2*norbs + (k2-1)*2 + jf + 2*nosq
            ind2 = (k1-1)*2*norbs + (k2-1)*2 + jfhf + 2*nosq

            v(ind1) = v(ind2)
            v(ind1+1) = v(ind2+1)
 1208     continue
 1207   continue

      elseif ( itype .eq. 3 ) then

        do 1301 k1 = 1, nel
          do 1302 k2 = 1, nel
            ind1 = (k1-1)*4*norbs + (k2-1)*2 + jf
            ind2 = (k1-1)*4*norbs + (k2-1)*2 + jfhf

            v(ind1) = v(ind2)
            v(ind1+1) = v(ind2+1)
 1302     continue
 1301   continue

        do 1303 k1 = nel+1, 2*norbs
          do 1304 k2 = nel+1, 2*norbs
            ind1 = (k1-1)*4*norbs + (k2-1)*2 + jf
            ind2 = (k1-1)*4*norbs + (k2-1)*2 + jfhf

            v(ind1) = v(ind2)
            v(ind1+1) = v(ind2+1)
 1304     continue
 1303   continue

      endif


C     Transform Fock matrix back into orthonomal AO basis:
C       nos . F . nos!

      if ( itype .eq. 1 ) then

        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jf), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'c', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos), norbs, zero, v(jf), norbs)

      elseif ( itype .eq. 2 ) then

        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnos), norbs,
     $       v(jf), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'c', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos), norbs, zero, v(jf), norbs)

        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jnos+2*nosq),
     $       norbs, v(jf+2*nosq), norbs, zero, v(jnost), norbs)
        call zgemm ('n', 'c', norbs, norbs, norbs, one, v(jnost), norbs,
     $       v(jnos+2*nosq), norbs, zero, v(jf+2*nosq), norbs)

      elseif ( itype .eq. 3 ) then

        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jnos),
     $       2*norbs, v(jf), 2*norbs, zero, v(jnost), 2*norbs)
        call zgemm ('n', 'c', 2*norbs, 2*norbs, 2*norbs, one, v(jnost),
     $       2*norbs, v(jnos), 2*norbs, zero, v(jf), 2*norbs)
      endif


C     Zero-out alpha-beta and beta-alpha spin blocks for ispuhf = .true.
C     Transform into UHF-type Fock matrix (only alpha-alpha and beta
C     -beta blocks are stored).

      if ( ispuhf ) then
        call dmblck (iout, iprint, v(jscr), lscr, norbs, v(jf),
     $       v(jzaa), v(jzab), v(jzba), v(jzbb), 2, -1)

        call amove (2*nosq, v(jzaa), v(jf))
        call amove (2*nosq, v(jzbb), v(jf+2*nosq))
      endif


C     Hermitize the effective Fock matrix and store as lower triangular.

      if ( itype .eq. 1 ) then
        call apadag (1, norbs, 0.50d0, v(jf), v(jf2))

        call trspnc (.false., norbs, norbs, v(jf2))
        call zlinear (v(jf2), f, norbs, norbs)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call apadag (1, norbs, 0.50d0, v(jf), v(jf2))

        call trspnc (.false., norbs, norbs, v(jf2))
        call zlinear (v(jf2), f, norbs, norbs)

        call apadag (1, norbs, 0.50d0, v(jf+2*nosq), v(jf2))

        call trspnc (.false., norbs, norbs, v(jf2))
        call zlinear (v(jf2), f(1+2*ntto), norbs, norbs)

      elseif ( itype .eq. 3 ) then
        call apadag (1, 2*norbs, 0.50d0, v(jf), v(jf2))

        call trspnc (.false., 2*norbs, 2*norbs, v(jf2))
        call zlinear (v(jf2), f, 2*norbs, 2*norbs)
      endif

      goto 300



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Use projection operators  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 200  continue


C     Memory allocation.

C       Allocate space for:
C         - scratch copy of HF density matrix (jp)
C         - scratch copy of (I-P(HF)) (jq)
C         - scratch array for matrix multiplication (jt)
C         - scratch copy of PHF Fock matrix (jf)
C         - scratch copy of Fock matrix for LT transformation (jf2)
C         - scratch copy of HF Fock matrix (jfhf)
C         - blocks of Fock matrix (jfoo, jfov, jfvo, jfvv)

      if ( itype .eq. 1 ) then
        jp   = 1
        jq   = jp   + 2*nosq
        jt   = jq   + 2*nosq
        jf   = jt   + 2*nosq
        jf2  = jf   + 2*nosq
        jfhf = jf2  + 2*nosq
        jfoo = jfhf + 2*nosq
        jfov = jfoo + 2*nosq
        jfvo = jfov + 2*nosq
        jfvv = jfvo + 2*nosq
        jend = jfvv + 2*nosq

      elseif ( itype .eq. 2 .or. ispuhf ) then
        jp   = 1
        jq   = jp   + 4*nosq
        jt   = jq   + 4*nosq
        jf   = jt   + 2*nosq
        jf2  = jf   + 4*nosq
        jfhf = jf2  + 2*nosq
        jfoo = jfhf + 4*nosq
        jfov = jfoo + 4*nosq
        jfvo = jfov + 4*nosq
        jfvv = jfvo + 4*nosq
        jend = jfvv + 4*nosq

      elseif ( itype .eq. 3 ) then
        jp   = 1
        jq   = jp   + 8*nosq
        jt   = jq   + 8*nosq
        jf   = jt   + 8*nosq
        jf2  = jf   + 8*nosq
        jfhf = jf2  + 8*nosq
        jfoo = jfhf + 8*nosq
        jfov = jfoo + 8*nosq
        jfvo = jfov + 8*nosq
        jfvv = jfvo + 8*nosq
        jend = jfvv + 8*nosq
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


      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'focksel')



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


C     Load PHF and HF Fock matrices from RWF files.

      if ( itype .eq. 1 ) then
        call fileio (2, -irwf, 2*ntto, v(jf), 0)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call fileio (2, -irwf, 2*ntto, v(jf), 0)
        call fileio (2, -irwf, 2*ntto, v(jf+2*nosq), 2*ntto)

      elseif ( itype .eq. 3 ) then
        call fileio (2, -irwf, 2*ntt2o, v(jf), 0)
      endif

      if ( itype .eq. 1 ) then
        call fileio (2, -irwfhf, 2*ntto, v(jfhf), 0)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call fileio (2, -irwfhf, 2*ntto, v(jfhf), 0)
        call fileio (2, -irwfhf, 2*ntto, v(jfhf+2*nosq), 2*ntto)

      elseif ( itype .eq. 3 ) then
        call fileio (2, -irwfhf, 2*ntt2o, v(jfhf), 0)
      endif



C     Transform the PHF Fock matrix (orthonormal AO basis) into square
C     format.

      if ( itype .eq. 1 ) then
        call zrecmp (1, v(jf), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jf), v(jz2), nosq)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call zrecmp (1, v(jf), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jf), v(jz2), nosq)

        call zrecmp (1, v(jf+2*nosq), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jf+2*nosq), v(jz2), nosq)

      elseif ( itype .eq. 3 ) then
        call zrecmp (1, v(jf), v(jz1), ntt2o)
        call square (v(jz1), v(jz2), 2*norbs, 2*norbs, 0)
        call square (v(jz1+ntt2o), v(jz3), 2*norbs, 2*norbs, -1)
        call zrecmp (2, v(jf), v(jz2), 4*nosq)
      endif


C     Transform the HF Fock matrix (orthonormal AO basis) into square
C     format.

      if ( itype .eq. 1 ) then
        call zrecmp (1, v(jfhf), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jfhf), v(jz2), nosq)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call zrecmp (1, v(jfhf), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jfhf), v(jz2), nosq)

        call zrecmp (1, v(jfhf+2*nosq), v(jz1), ntto)
        call square (v(jz1), v(jz2), norbs, norbs, 0)
        call square (v(jz1+ntto), v(jz3), norbs, norbs, -1)
        call zrecmp (2, v(jfhf+2*nosq), v(jz2), nosq)

      elseif ( itype .eq. 3 ) then
        call zrecmp (1, v(jfhf), v(jz1), ntt2o)
        call square (v(jz1), v(jz2), 2*norbs, 2*norbs, 0)
        call square (v(jz1+ntt2o), v(jz3), 2*norbs, 2*norbs, -1)
        call zrecmp (2, v(jfhf), v(jz2), 4*nosq)
      endif


C     Form the quantity Q = I-P.

      if ( itype .eq. 1 ) then
        call aunitm (.false., 2, norbs, norbs, v(jq))
        call asub (2*nosq, v(jq), v(jp), v(jq))

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call aunitm (.false., 2, norbs, norbs, v(jq))
        call asub (2*nosq, v(jq), v(jp), v(jq))

        call aunitm (.false., 2, norbs, norbs, v(jq+2*nosq))
        call asub (2*nosq, v(jq+2*nosq), v(jp+2*nosq), v(jq+2*nosq))

      elseif ( itype .eq. 3 ) then
        call aunitm (.false., 2, 2*norbs, 2*norbs, v(jq))
        call asub (8*nosq, v(jq), v(jp), v(jq))
      endif


C     Form occupied-occupied block as
C       F_oo  =  P . F(HF) . P

      if ( itype .eq. 1 ) then
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jp), norbs,
     $       v(jfhf), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jp), norbs, zero, v(jfoo), norbs)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jp), norbs,
     $       v(jfhf), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jp), norbs, zero, v(jfoo), norbs)

        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jp+2*nosq),
     $       norbs, v(jfhf+2*nosq), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jp+2*nosq), norbs, zero, v(jfoo+2*nosq), norbs)

      elseif ( itype .eq. 3 ) then
        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jp),
     $       2*norbs, v(jfhf), 2*norbs, zero, v(jt), 2*norbs)
        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jt),
     $       2*norbs, v(jp), 2*norbs, zero, v(jfoo), 2*norbs)
      endif


C     Form virtual-virtual block as
C       F_vv  =  Q . F(HF) . Q

      if ( itype .eq. 1 ) then
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jq), norbs,
     $       v(jfhf), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jq), norbs, zero, v(jfvv), norbs)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jq), norbs,
     $       v(jfhf), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jq), norbs, zero, v(jfvv), norbs)

        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jq+2*nosq),
     $       norbs, v(jfhf+2*nosq), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jq+2*nosq), norbs, zero, v(jfvv+2*nosq), norbs)

      elseif ( itype .eq. 3 ) then
        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jq),
     $       2*norbs, v(jfhf), 2*norbs, zero, v(jt), 2*norbs)
        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jt),
     $       2*norbs, v(jq), 2*norbs, zero, v(jfvv), 2*norbs)
      endif


C     Form occupied-virtual block as
C       F_ov  =  P . F(PHF) . Q

      if ( itype .eq. 1 ) then
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jp), norbs,
     $       v(jf), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jq), norbs, zero, v(jfov), norbs)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jp), norbs,
     $       v(jf), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jq), norbs, zero, v(jfov), norbs)

        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jp+2*nosq),
     $       norbs, v(jf+2*nosq), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jq+2*nosq), norbs, zero, v(jfov+2*nosq), norbs)

      elseif ( itype .eq. 3 ) then
        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jp),
     $       2*norbs, v(jf), 2*norbs, zero, v(jt), 2*norbs)
        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jt),
     $       2*norbs, v(jq), 2*norbs, zero, v(jfov), 2*norbs)
      endif


C     Form virtual-occupied block as
C       F_vo  =  Q . F(PHF) . P

      if ( itype .eq. 1 ) then
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jq), norbs,
     $       v(jf), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jp), norbs, zero, v(jfvo), norbs)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jq), norbs,
     $       v(jf), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jp), norbs, zero, v(jfvo), norbs)

        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jq+2*nosq),
     $       norbs, v(jf+2*nosq), norbs, zero, v(jt), norbs)
        call zgemm ('n', 'n', norbs, norbs, norbs, one, v(jt), norbs,
     $       v(jp+2*nosq), norbs, zero, v(jfvo+2*nosq), norbs)

      elseif ( itype .eq. 3 ) then
        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jq),
     $       2*norbs, v(jf), 2*norbs, zero, v(jt), 2*norbs)
        call zgemm ('n', 'n', 2*norbs, 2*norbs, 2*norbs, one, v(jt),
     $       2*norbs, v(jp), 2*norbs, zero, v(jfvo), 2*norbs)
      endif


C     Add all contributions to Fock.

      if ( itype .eq. 1 ) then
        call aclear (2*nosq, v(jf))

        call aadd (2*nosq, v(jf), v(jfoo), v(jf))
        call aadd (2*nosq, v(jf), v(jfov), v(jf))
        call aadd (2*nosq, v(jf), v(jfvo), v(jf))
        call aadd (2*nosq, v(jf), v(jfvv), v(jf))

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call aclear (2*nosq, v(jf))
        call aclear (2*nosq, v(jf+2*nosq))

        call aadd (2*nosq, v(jf), v(jfoo), v(jf))
        call aadd (2*nosq, v(jf), v(jfov), v(jf))
        call aadd (2*nosq, v(jf), v(jfvo), v(jf))
        call aadd (2*nosq, v(jf), v(jfvv), v(jf))

        call aadd (2*nosq, v(jf+2*nosq), v(jfoo+2*nosq), v(jf+2*nosq))
        call aadd (2*nosq, v(jf+2*nosq), v(jfov+2*nosq), v(jf+2*nosq))
        call aadd (2*nosq, v(jf+2*nosq), v(jfvo+2*nosq), v(jf+2*nosq))
        call aadd (2*nosq, v(jf+2*nosq), v(jfvv+2*nosq), v(jf+2*nosq))

      elseif ( itype .eq. 3 ) then
        call aclear (8*nosq, v(jf))

        call aadd (8*nosq, v(jf), v(jfoo), v(jf))
        call aadd (8*nosq, v(jf), v(jfov), v(jf))
        call aadd (8*nosq, v(jf), v(jfvo), v(jf))
        call aadd (8*nosq, v(jf), v(jfvv), v(jf))
      endif


C     Hermitize the effective Fock matrix and store as lower triangular.

      if ( itype .eq. 1 ) then
        call apadag (1, norbs, 0.50d0, v(jf), v(jf2))

        call trspnc (.false., norbs, norbs, v(jf2))
        call zlinear (v(jf2), f, norbs, norbs)

      elseif ( itype .eq. 2 .or. ispuhf ) then
        call apadag (1, norbs, 0.50d0, v(jf), v(jf2))

        call trspnc (.false., norbs, norbs, v(jf2))
        call zlinear (v(jf2), f, norbs, norbs)

        call apadag (1, norbs, 0.50d0, v(jf+2*nosq), v(jf2))

        call trspnc (.false., norbs, norbs, v(jf2))
        call zlinear (v(jf2), f(1+2*ntto), norbs, norbs)

      elseif ( itype .eq. 3 ) then
        call apadag (1, 2*norbs, 0.50d0, v(jf), v(jf2))

        call trspnc (.false., 2*norbs, 2*norbs, v(jf2))
        call zlinear (v(jf2), f, 2*norbs, 2*norbs)
      endif

      goto 300



C     %%%%%%%%%%%%%%%%%%%%
C     %  Level-shifting  %
C     %%%%%%%%%%%%%%%%%%%%


 300   continue


C     Memory allocation.

C       Allocate space for:
C         - scratch copy of HF density matrix (jp)

      if ( itype .eq. 1 ) then
        jp   = 1
        jend = jp + 2*ntto

      elseif ( itype .eq. 2 ) then
        jp   = 1
        jend = jp + 4*ntto

      elseif ( itype .eq. 3 ) then
        jp   = 1
        jend = jp + 2*ntt2o
      endif

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'focksel')


C     Apply level shifting.

C       ivshft =  0,  constant level shifting
C              =  1,  dynamically adjust to achieve a gap of at
C                     least vsgoal

      if ( ivshft .gt. 0 ) then
        vshft = abs (dble (ivshft) / dble (1000))

      elseif ( ivshft .lt. -1 ) then
        if ( gap .le. vsgoal ) then
          vshft = vsgoal - gap
        else
          vshft = 0.0d0
        endif
      endif

      if ( ivshft .ne. -1 ) then
        if ( itype .eq. 1 ) then
          call amove (2*ntto, p, v(jp))
          call zdscal (ntto, -vshft, v(jp), 1)
          call aadd (2*ntto, f, v(jp), f)

        elseif ( itype .eq. 2 .or. ispuhf ) then
          call amove (4*ntto, p, v(jp))
          call zdscal (2*ntto, -vshft, v(jp), 1)
          call aadd (4*ntto, f, v(jp), f)

        elseif ( itype .eq. 3 ) then
          call amove (2*ntt2o, p, v(jp))
          call zdscal (ntt2o, -vshft, v(jp), 1)
          call aadd (2*ntt2o, f, v(jp), f)
        endif
      endif


C     Save matrices in RWF file.

      if ( itype .eq. 1 ) then
        call fileio (1, -irwf, 2*ntto, f, 0)
      elseif ( itype .eq. 2 .or. ispuhf ) then
        call fileio (1, -irwf, 4*ntto, f, 0)
      elseif ( itype .eq. 3 ) then
        call fileio (1, -irwf, 2*ntt2o, f, 0)
      endif


C     Print effective Fock matrices.

      if ( iprint .gt. 3 ) then
        if ( itype .eq. 1 ) then
          write (iout, *) ' '
          call ltoutx (iout, 0, 2, 1, 0, 'Effective Fock', norbs,
     $         norbs, f, f, 0)

        elseif ( itype .eq. 2 .or. ispuhf ) then
          write (iout, *) ' '
          call ltoutx (iout, 0, 3, 1, 0, 'Effective Fock', norbs,
     $         norbs, f, f(1+2*ntto), 0)

        elseif ( itype .eq. 3 ) then
          write (iout, *) ' '
          call ltoutx (iout, 0, 2, 1, 0, 'Effective Fock', 2*norbs,
     $         2*norbs, f, f, 0)
        endif
      endif



      return
      end


