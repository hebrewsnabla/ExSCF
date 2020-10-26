

      subroutine dmblck (iout, iprint, v, lenv, norbs, xfull, xaa, xab,
     $     xba, xbb, mode, key)

      implicit none


C +------------------------------------------------------------+
C |  dmblck   --   CAJH, 07.2011                               |
C |                                                            |
C |                                    (based on PHFB dmblck)  |
C |                                                            |
C |                                                            |
C |  Construct the full X matrix from its spin-blocks, or      |
C |  split the full X matrix into its spin-blocks.             |
C |                                                            |
C |  ( X and its spin blocks are always assumed complex. )     |
C |                                                            |
C |    mode - controls what is being provided and what should  |
C |           be returned as output                            |
C |                                                            |
C |    key  - specifies the characteristics of X and how       |
C |           results should be returned                       |
C |                                                            |
C |  A summary of implemented combinations is shown next.      |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  mode = 1  ::  build full matrix                           |
C |                                                            |
C |    key =  1, X is a general matrix, all spin blocks are    |
C |              provided in square form                       |
C |                ( Xfull is returned as square )             |
C |                                                            |
C |    key = -1, same as key = 1                               |
C |                                                            |
C |    key =  2, X is a Hermitian matrix, spin blocks are      |
C |              provided in compact form                      |
C |                ( Xfull is returned lower triangular )      |
C |                                                            |
C |    key = -2, same as key = 2                               |
C |                ( return Xfull as square regardless )       |
C |                                                            |
C |    key =  3, X is an antisymmetric matrix, spin blocks     |
C |              are provided in compact form                  |
C |                ( Xfull is returned lower triangular )      |
C |                                                            |
C |    key = -3, same as key = 3                               |
C |                ( return Xfull as square regardless )       |
C |                                                            |
C |  scratch space required:                                   |
C |                                                            |
C |    key =  1,  16*nosq                                      |
C |        =  2,  14*nosq + 2*ntto + 8*nosq                    |
C |        =  3,  14*nosq + 2*ntto + 8*nosq                    |
C |        = -1,  16*nosq                                      |
C |        = -2,  14*nosq + 2*ntto                             |
C |        = -3,  14*nosq + 2*ntto                             |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  mode = 2  ::  split full matrix into spin-blocks          |
C |                                                            |
C |    key =  1, X is a general square matrix ***              |
C |                ( spin blocks returned square )             |
C |                                                            |
C |    key = -1, same as key = 1 ***                           |
C |                                                            |
C |    key =  2, X is a Hermitian matrix, stored in lower      |
C |              triangular form                               |
C |                ( spin blocks are returned compact )        |
C |                                                            |
C |    key = -2, same as key = 2 ***                           |
C |                ( return spin blocks square regardless )    |
C |                                                            |
C |    key =  3, X is an antisymmetric matrix, stored in       |
C |              lower triangular form                         |
C |                ( spin blocks are returned compact )        |
C |                                                            |
C |    key = -3, same as key = 3 ***                           |
C |                ( return spin blocks square regardless )    |
C |                                                            |
C |  *** all spin blocks should be allocated square            |
C |                                                            |
C |  scratch space required:                                   |
C |                                                            |
C |    key =  1,  16*nosq                                      |
C |        =  2,  14*nosq + 2*ntt2o                            |
C |        =  3,  14*nosq + 2*ntt2o                            |
C |        = -1,  16*nosq                                      |
C |        = -2,  14*nosq + 2*ntt2o + 2*nosq                   |
C |        = -3,  14*nosq + 2*ntt2o + 2*nosq                   |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  Lastly, a word on what we understand as compact storing   |
C |  for spin blocks (for Hermitian or antisymmetric           |
C |  matrices):                                                |
C |                                                            |
C |    [alpha, alpha] block - lower triangular                 |
C |    [alpha, beta]  block - full square                      |
C |    [beta, alpha]  block - not allocated                    |
C |    [beta, beta]   block - lower triangular                 |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v     - scratch array
C       lenv  - size of scratch array
C       norbs - number of orbitals
C       xfull - full X matrix [ in / out ]
C       xaa   - [alpha, alpha] block of X matrix [ in / out ]
C       xab   - [alpha, beta]  block of X matrix [ in / out ]
C       xba   - [beta, alpha]  block of X matrix [ in / out ]
C       xbb   - [beta, beta]   block of X matrix [ in / out ]
C       mode  - whether to build full X matrix or build spin blocks
C       key   - specifies characteristics of X and controls how
C               to return results

      real*8      v(*), xfull(*), xaa(*), xab(*), xba(*), xbb(*)
      integer     iout, iprint, lenv, norbs
      integer     mode, key

C     other variables

      integer     ntto, ntt2o, nosq, jend, mdv
      integer     szaa, szab, szba, szbb
      integer     jxaar, jxaai, jxabr, jxabi
      integer     jxbar, jxbai, jxbbr, jxbbi
      integer     jxfr, jxfi, jxf, jscr1, jscr2
      integer     j, k, indj, indk, indr1, indr2, inds1, inds2
      logical     outsq, inpsq, inpht, inpas


C     Construct useful variables.

      ntto  = norbs * (norbs + 1) / 2
      ntt2o = norbs * (2*norbs + 1)

      nosq = norbs * norbs


C     Split the code according to mode.

      if ( mode .eq. 1 ) goto 10
      if ( mode .eq. 2 ) goto 20
      call gauerr ('mode not supported in dmblck')



 10   continue

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Build Xfull from spin blocks  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Decipher key.

C       outsq - whether to output as full square
C       inpsq - whether input is full square
C       inpht - whether input is lower triangular, Hermitian
C       inpas - whether input is lower triangular, antisymmetric

      outsq = .false.
      inpsq = .false.
      inpht = .false.
      inpas = .false.

      if ( key .lt. 0 ) then
        outsq = .true.
      endif

      if ( abs (key) .eq. 1 ) then
        inpsq = .true.
      elseif ( abs (key) .eq. 2 ) then
        inpht = .true.
      elseif ( abs (key) .eq. 3 ) then
        inpas = .true.
      else
        call gauerr ('key not supported in dmblck')
      endif


C     Memory allocation.

C       Allocate space to store:
C         - real and imaginary parts of each spin block
C         - real part of Xfull (jxfr)
C         - imag part of Xfull (jxfi)
C         - full Xfull (jxf) - only needed for Hermitian or
C                              antisymmetric outputs
C         - two scratch lower triangular arrays (jscr1, jscr2)
C             for Hermitian or antisymmetric inputs

      if ( inpsq ) then
        szba = nosq
      else
        szba = 0
      endif

      jxaar = 1
      jxaai = jxaar + nosq
      jxabr = jxaai + nosq
      jxabi = jxabr + nosq
      jxbar = jxabi + nosq
      jxbai = jxbar + szba
      jxbbr = jxbai + szba
      jxbbi = jxbbr + nosq
      jxfr  = jxbbi + nosq
      jxfi  = jxfr  + 4*nosq
      jend  = jxfi  + 4*nosq

      if ( .not. outsq ) then
        jxf  = jend
        jend = jxf  + 8*nosq
      endif

      if ( inpht .or. inpas ) then
        jscr1 = jend
        jscr2 = jscr1 + ntto
        jend  = jscr2 + ntto
      endif

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'dmblck')


C     Expand spin blocks into square forms, with separate real and
C     imaginary parts.

      if ( inpsq ) then
        call zrecmp (1, xaa, v(jxaar), nosq)
        call zrecmp (1, xab, v(jxabr), nosq)
        call zrecmp (1, xba, v(jxbar), nosq)
        call zrecmp (1, xbb, v(jxbbr), nosq)

      elseif ( inpht ) then
        call zrecmp (1, xaa, v(jscr1), ntto)
        call square (v(jscr1), v(jxaar), norbs, norbs, 0)
        call square (v(jscr2), v(jxaai), norbs, norbs, -1)

        call zrecmp (1, xbb, v(jscr1), ntto)
        call square (v(jscr1), v(jxbbr), norbs, norbs, 0)
        call square (v(jscr2), v(jxbbi), norbs, norbs, -1)

        call zrecmp (1, xab, v(jxabr), nosq)

      elseif ( inpas ) then
        call zrecmp (1, xaa, v(jscr1), ntto)
        call square (v(jscr1), v(jxaar), norbs, norbs, -1)
        call square (v(jscr2), v(jxaai), norbs, norbs, -1)

        call zrecmp (1, xbb, v(jscr1), ntto)
        call square (v(jscr1), v(jxbbr), norbs, norbs, -1)
        call square (v(jscr2), v(jxbbi), norbs, norbs, -1)

        call zrecmp (1, xab, v(jxabr), nosq)
      endif


C     Clear spaces for real and imaginary Xfull.

      call aclear (4*nosq, v(jxfr))
      call aclear (4*nosq, v(jxfi))


C     Form Xfull
C     ==========

C     The full matrix should be constructed from its spin blocks as

C       Xfull =  [ Xaa   Xab ]
C                [ Xba   Xbb ]

C       if Xfull is Hermitian, then Xba = Xab'
C       if Xfull is antisymmetric, then Xba = -transpose(Xab)


C       .. [alpha, alpha] block

      do 101 j = 1, norbs
        do 102 k = 1, norbs
          indj = j
          indk = k

          indr1 = jxfr + (k - 1) * 2*norbs + j - 1
          indr2 = jxfi + (k - 1) * 2*norbs + j - 1

          inds1 = jxaar + (indk - 1) * norbs + indj - 1
          inds2 = jxaai + (indk - 1) * norbs + indj - 1

          v(indr1) = v(indr1) + v(inds1)
          v(indr2) = v(indr2) + v(inds2)
 102    continue
 101  continue

C       .. [alpha, beta] block

      do 103 j = 1, norbs
        do 104 k = norbs+1, 2*norbs
          indj = j
          indk = k - norbs

          indr1 = jxfr + (k - 1) * 2*norbs + j - 1
          indr2 = jxfi + (k - 1) * 2*norbs + j - 1

          inds1 = jxabr + (indk - 1) * norbs + indj - 1
          inds2 = jxabi + (indk - 1) * norbs + indj - 1

          v(indr1) = v(indr1) + v(inds1)
          v(indr2) = v(indr2) + v(inds2)
 104    continue
 103  continue

C       .. [beta, alpha] block

      do 105 j = norbs+1, 2*norbs
        do 106 k = 1, norbs
          indj = j - norbs
          indk = k

          indr1 = jxfr + (k - 1) * 2*norbs + j - 1
          indr2 = jxfi + (k - 1) * 2*norbs + j - 1
          
          if ( inpsq ) then
            inds1 = jxbar + (indk - 1) * norbs + indj - 1
            inds2 = jxbai + (indk - 1) * norbs + indj - 1

          elseif ( inpht .or. inpas ) then
            inds1 = jxabr + (indj - 1) * norbs + indk - 1
            inds2 = jxabi + (indj - 1) * norbs + indk - 1
          endif

          if ( inpsq ) then
            v(indr1) = v(indr1) + v(inds1)
            v(indr2) = v(indr2) + v(inds2)

          elseif ( inpht ) then
            v(indr1) = v(indr1) + v(inds1)
            v(indr2) = v(indr2) - v(inds2)

          elseif ( inpas ) then
            v(indr1) = v(indr1) - v(inds1)
            v(indr2) = v(indr2) - v(inds2)
          endif
 106    continue
 105  continue

C       .. [beta, beta] block

      do 107 j = norbs+1, 2*norbs
        do 108 k = norbs+1, 2*norbs
          indj = j - norbs
          indk = k - norbs

          indr1 = jxfr + (k - 1) * 2*norbs + j - 1
          indr2 = jxfi + (k - 1) * 2*norbs + j - 1

          inds1 = jxbbr + (indk - 1) * norbs + indj - 1
          inds2 = jxbbi + (indk - 1) * norbs + indj - 1

          v(indr1) = v(indr1) + v(inds1)
          v(indr2) = v(indr2) + v(inds2)
 108    continue
 107  continue


C     Merge real and imaginary parts of Xfull.
C     Transform Xfull to lower triangular if output is not square.
      

      if ( outsq ) then
        call zrecmp (2, xfull, v(jxfr), 4*nosq)

      else
        call zrecmp (2, v(jxf), v(jxfr), 4*nosq)

        call trspnc (.false., 2*norbs, 2*norbs, v(jxf))
        call zlinear (v(jxf), xfull, 2*norbs, 2*norbs)
      endif


      return




 20   continue

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Split Xfull into spin blocks  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C     Decipher key.

C       outsq - whether to output as square
C       inpsq - whether input is full square
C       inpht - whether input is lower triangular, Hermitian
C       inpas - whether input is lower triangular, antisymmetric

      outsq = .false.
      inpsq = .false.
      inpht = .false.
      inpas = .false.

      if ( key .lt. 0 ) then
        outsq = .true.
      endif

      if ( abs (key) .eq. 1 ) then
        inpsq = .true.
      elseif ( abs (key) .eq. 2 ) then
        inpht = .true.
      elseif ( abs (key) .eq. 3 ) then
        inpas = .true.
      else
        call gauerr ('key not supported in dmblck')
      endif


C     Memory allocation.

C       Allocate space to store:
C         - real and imaginary parts of each spin block
C         - real part of Xfull (jxfr)
C         - imag part of Xfull (jxfi)
C         - two scratch lower triangular arrays (jscr1, jscr2)
C             for Hermitian or antisymmetric inputs

      if ( outsq .or. inpsq ) then
        szba = nosq
      else
        szba = 0
      endif

      jxaar = 1
      jxaai = jxaar + nosq
      jxabr = jxaai + nosq
      jxabi = jxabr + nosq
      jxbar = jxabi + nosq
      jxbai = jxbar + szba
      jxbbr = jxbai + szba
      jxbbi = jxbbr + nosq
      jxfr  = jxbbi + nosq
      jxfi  = jxfr  + 4*nosq
      jend  = jxfi  + 4*nosq

      if ( inpht .or. inpas ) then
        jscr1 = jend
        jscr2 = jscr1 + ntt2o
        jend  = jscr2 + ntt2o
      endif

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'dmblck')


C     Expand full matrix into square form, with separate real and
C     imaginary parts.

      if ( inpsq ) then
        call zrecmp (1, xfull, v(jxfr), 4*nosq)

      elseif ( inpht ) then
        call zrecmp (1, xfull, v(jscr1), ntt2o)
        call square (v(jscr1), v(jxfr), 2*norbs, 2*norbs, 0)
        call square (v(jscr2), v(jxfi), 2*norbs, 2*norbs, -1)

      elseif ( inpas ) then
        call zrecmp (1, xfull, v(jscr1), ntt2o)
        call square (v(jscr1), v(jxfr), 2*norbs, 2*norbs, -1)
        call square (v(jscr2), v(jxfi), 2*norbs, 2*norbs, -1)
      endif


C     Clear spaces for real and imaginary spin blocks.

      call aclear (nosq, v(jxaar))
      call aclear (nosq, v(jxaai))
      call aclear (nosq, v(jxabr))
      call aclear (nosq, v(jxabi))
      call aclear (nosq, v(jxbbr))
      call aclear (nosq, v(jxbbi))

      if ( outsq .or. inpsq ) then
        call aclear (nosq, v(jxbar))
        call aclear (nosq, v(jxbai))
      endif


C     Split Xfull into spin blocks
C     ============================

C     The full matrix can be split into its spin blocks like

C       Xfull =  [ Xaa   Xab ]
C                [ Xba   Xbb ]

C       if Xfull is Hermitian, then Xba = Xab'
C       if Xfull is antisymmetric, then Xba = -transpose(Xab)


C       .. [alpha, alpha] block

      do 201 j = 1, norbs
        do 202 k = 1, norbs
          indj = j
          indk = k

          indr1 = jxfr + (k - 1) * 2*norbs + j - 1
          indr2 = jxfi + (k - 1) * 2*norbs + j - 1

          inds1 = jxaar + (indk - 1) * norbs + indj - 1
          inds2 = jxaai + (indk - 1) * norbs + indj - 1

          v(inds1) = v(inds1) + v(indr1)
          v(inds2) = v(inds2) + v(indr2)
 202    continue
 201  continue

C       .. [alpha, beta] block

      do 203 j = 1, norbs
        do 204 k = norbs+1, 2*norbs
          indj = j
          indk = k - norbs

          indr1 = jxfr + (k - 1) * 2*norbs + j - 1
          indr2 = jxfi + (k - 1) * 2*norbs + j - 1

          inds1 = jxabr + (indk - 1) * norbs + indj - 1
          inds2 = jxabi + (indk - 1) * norbs + indj - 1

          v(inds1) = v(inds1) + v(indr1)
          v(inds2) = v(inds2) + v(indr2)
 204    continue
 203  continue

C       .. [beta, alpha] block

      if ( outsq .or. inpsq ) then
        do 205 j = norbs+1, 2*norbs
          do 206 k = 1, norbs
            indj = j - norbs
            indk = k

            indr1 = jxfr + (k - 1) * 2*norbs + j - 1
            indr2 = jxfi + (k - 1) * 2*norbs + j - 1
          
            inds1 = jxbar + (indk - 1) * norbs + indj - 1
            inds2 = jxbai + (indk - 1) * norbs + indj - 1

            v(inds1) = v(inds1) + v(indr1)
            v(inds2) = v(inds2) + v(indr2)
 206      continue
 205    continue
      endif

C       .. [beta, beta] block

      do 207 j = norbs+1, 2*norbs
        do 208 k = norbs+1, 2*norbs
          indj = j - norbs
          indk = k - norbs

          indr1 = jxfr + (k - 1) * 2*norbs + j - 1
          indr2 = jxfi + (k - 1) * 2*norbs + j - 1

          inds1 = jxbbr + (indk - 1) * norbs + indj - 1
          inds2 = jxbbi + (indk - 1) * norbs + indj - 1

          v(inds1) = v(inds1) + v(indr1)
          v(inds2) = v(inds2) + v(indr2)
 208    continue
 207  continue


C     Merge real and imaginary parts of each spin block.
C     Transform [alpha, alpha] and [beta, beta] spin blocks to lower
C     triangular for Hermitian or antisymmetric outputs.
      

      if ( outsq .or. inpsq ) then
        call zrecmp (2, xaa, v(jxaar), nosq)
        call zrecmp (2, xab, v(jxabr), nosq)
        call zrecmp (2, xba, v(jxbar), nosq)
        call zrecmp (2, xbb, v(jxbbr), nosq)

      else
        call zrecmp (2, xab, v(jxabr), nosq)

C       .. we use the array jxabr as scratch space here ..

        call zrecmp (2, v(jxabr), v(jxaar), nosq)
        call trspnc (.false., norbs, norbs, v(jxabr))
        call zlinear (v(jxabr), xaa, norbs, norbs)

        call zrecmp (2, v(jxabr), v(jxbbr), nosq)
        call trspnc (.false., norbs, norbs, v(jxabr))
        call zlinear (v(jxabr), xbb, norbs, norbs)
      endif


      return
      end

