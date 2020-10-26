

      subroutine detmt2 (n, veig, adet, pdet)

      implicit none


C +------------------------------------------------------------+
C |  detmt2   --   CAJH, 07.2011                               |
C |                                                            |
C |                                    (based on PHFB detmt2)  |
C |                                                            |
C |                                                            |
C |  Compute the absolute value and the phase of the           |
C |  determinant of a SQUARE, COMPLEX matrix A from its        |
C |  eigenvalues (veig).                                       |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       n    - dimension of the vector of eigenvalues
C       veig - vector with eigenvalues of matrix whose determinant needs
C              to be computed
C       adet - absolute value of determinant [ out ]
C       pdet - phase of determinant [ out ]

      real*8      adet, pdet
      integer     n
      complex*16  veig(*)

C     other variables

      integer     i
      real*8      xr, xi
      real*8      sphs, sabs, phsx, absx, logx


C     Initialize values of phase and the product of absolute values.

      sphs = 0.0d0
      sabs = 0.0d0

C     Loop over vector of eigenvalues.

      do 20 i = 1, n

C     Recover real and imaginary part of a(i,i).

        xr = dble (veig(i))
        xi = aimag (veig(i))

C     Compute phase and absolute value.

        phsx = atan2 (xi, xr)
        absx = sqrt (xr**2 + xi**2)
        logx = log (absx)

C     Add to overall quantities.

        sphs = sphs + phsx
        sabs = sabs + logx
 20   continue


C     Prepare quantities to return:
C       adet - absolute value of the determinant
C       pdet - phase of the determinant

      adet = exp (sabs)
      pdet = sphs


      return
      end


