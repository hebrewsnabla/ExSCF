      

      subroutine sptblk (iout, iprint, v, lenv, nbasis, x, xrs, xra,
     $     xis, xia, mode)

      implicit none


C +------------------------------------------------------------+
C |  sptblk   --   CAJH, 07.2011                               |
C |                                                            |
C |                                    (based on PHFB sptblk)  |
C |                                                            |
C |                                                            |
C |  Split a matrix into:                                      |
C |                                                            |
C |    real-symmetric block (xrs)                              |
C |    real-antisymmetric block (xra)                          |
C |    imag-symmetric block (xis)                              |
C |    imag-antisymmetric block (xia)                          |
C |                                                            |
C |  or do the opposite (combine). Namely,                     |
C |                                                            |
C |    mode = 1,  split the matrix x                           |
C |         = 2,  combine xrs, xra, xis, xia into x            |
C |                                                            |
C |  In all cases, x should be allocated square, while each    |
C |  block is stored lower triangular.                         |
C |                                                            |
C |  scratch space required:                                   |
C |                                                            |
C |    mode = 1,  2*nbsq + 4*ntt                               |
C |         = 2,  3*nbsq                                       |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       nbasis - number of basis functions
C       x      - full matrix [ in / out ]
C       xrs    - real, symmetric part of x     [ in / out ]
C       xra    - real, antisymmetric part of x [ in / out ]
C       xis    - imag, symmetric part of x     [ in / out ]
C       xia    - imag, antisymmetric part of x [ in / out ]
C       key    - determines which blocks are non-zero
C       mode   - decides whether to split or combine the matrices

      real*8      v(*), x(*), xrs(*), xra(*), xis(*), xia(*)
      integer     iout, iprint, lenv
      integer     nbasis, mode

C     other variables

      integer     ntt, nbsq
      integer     jscr1, jscr2, jscr3, jend, mdv


C     Construct useful variables.

      ntt = nbasis * (nbasis + 1) / 2
      nbsq = nbasis * nbasis


C     Memory allocation.

      if ( mode .eq. 1 ) then
        jscr1 = 1
        jscr2 = jscr1 + 2*nbsq
        jscr3 = jscr2 + 2*ntt
        jend  = jscr3 + 2*ntt

      elseif ( mode .eq. 2 ) then
        jscr1 = 1
        jscr2 = jscr1 + 2*nbsq
        jend  = jscr2 + nbsq
      endif

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'sptblk')



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Split matrix into blocks  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      if ( mode .eq. 1 ) then

C       Form x + x!.

        call amove (2*nbsq, x, v(jscr1))
        call trspnc (.true., nbasis, nbasis, v(jscr1))

        call aadd (2*nbsq, x, v(jscr1), v(jscr1))

C       .. recover xrs and xia ..

        call trspnc (.false., nbasis, nbasis, v(jscr1))
        call zlinear (v(jscr1), v(jscr2), nbasis, nbasis)
        call zrecmp (1, v(jscr2), v(jscr3), ntt)
        call ascale (2*ntt, 0.50d0, v(jscr3), v(jscr3))

        call amove (ntt, v(jscr3), xrs)
        call amove (ntt, v(jscr3+ntt), xia)

C       Form x - x!.

        call amove (2*nbsq, x, v(jscr1))
        call trspnc (.true., nbasis, nbasis, v(jscr1))

        call asub (2*nbsq, x, v(jscr1), v(jscr1))

C       .. recover xra and xis ..

        call trspnc (.false., nbasis, nbasis, v(jscr1))
        call zlinear (v(jscr1), v(jscr2), nbasis, nbasis)
        call zrecmp (1, v(jscr2), v(jscr3), ntt)
        call ascale (2*ntt, 0.50d0, v(jscr3), v(jscr3))

        call amove (ntt, v(jscr3), xra)
        call amove (ntt, v(jscr3+ntt), xis)

      endif



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Combine blocks into full matrix  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      if ( mode .eq. 2 ) then

C       Clear matrix x.

        call aclear (2*nbsq, x)
        call aclear (2*nbsq, v(jscr1))

C       Form the Hermitian part of x.

        call square (xrs, v(jscr2), nbasis, nbasis, 0)
        call aadd (nbsq, v(jscr2), v(jscr1), v(jscr1))

        call square (xia, v(jscr2), nbasis, nbasis, -1)
        call aadd (nbsq, v(jscr2), v(jscr1+nbsq), v(jscr1+nbsq))

C       Now form the anti-Hermitian part of x.

        call square (xra, v(jscr2), nbasis, nbasis, -1)
        call aadd (nbsq, v(jscr2), v(jscr1), v(jscr1))

        call square (xis, v(jscr2), nbasis, nbasis, 0)
        call aadd (nbsq, v(jscr2), v(jscr1+nbsq), v(jscr1+nbsq))

        call zrecmp (2, x, v(jscr1), nbsq)

      endif


      return
      end


