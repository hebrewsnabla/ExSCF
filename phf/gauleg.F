

      subroutine gauleg (x1, x2, grd, wgt, n)

      implicit none


C +------------------------------------------------------------+
C |  gauleg   --   CAJH, 07.2011                               |
C |                                                            |
C |                                    (based on PHFB gauleg)  |
C |                                                            |
C |                                                            |
C |  Subroutine to determine the abscissas and weights of the  |
C |  Gauss-Legendre n-point quadrature formula. Copied from    |
C |  Numerical Recipes.                                        |
C |                                                            |
C |  Given the lower and upper limits of integration x1 and    |
C |  x2, and given n, this routine returns arrays grd(1:n)     |
C |  and wgt(1:n) of length n, containing the abscissas and    |
C |  weights of the Gauss-Legendre n-point quadrature formula. |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       n      - number of grid points
C       x1, x2 - lower and upper limits of integration
C       grd    - array of grid positions [ out ]
C       wgt    - array of grid weights [ out ]

      integer     n
      real*8      x1, x2, grd(n), wgt(n)

C     other variables

      real*8      eps, pi, xm, xl
      real*8      z, z1, p1, p2, p3, pp
      integer     m, i, j

      parameter   ( eps = 3.0d-14 )


      pi = 4.0d0 * atan (1.0d0)


C     Only half of the roots need to be found given that the polynomial
C     is symmetric about the center of the interval.

      m = (n + 1) / 2

      xm = 0.5d0 * (x2 + x1)
      xl = 0.5d0 * (x2 - x1)

      do 10 i = 1, m

C       Approximation to the i-th root.

        z = cos (pi * (dble(i) - 0.25d0) / (dble (n) + 0.50d0))

C       Refine the polynomial root by Newton's method.

 20     continue

        p1 = 1.0d0
        p2 = 0.0d0

C       Loop over the recurrence relation to evaluate the Legendre
C       polynomial at position z.

        do 30 j = 1, n
          p3 = p2
          p2 = p1
          p1 = ((2.0d0 * dble (j) - 1.0d0) * z * p2 -
     $          (dble (j) - 1.0d0) * p3) / dble (j)
 30     continue

C       p1 is now the desired Legendre polynomial. We next compute pp,
C       its derivative, by a standard relation involving also p2, the
C       polynomial of one lower order.

        pp = dble (n) * (z * p1 - p2) / (z*z - 1.0d0)
        z1 = z

C       Newton's method ...

        z = z1 - p1 / pp

        if ( abs (z-z1) .gt. eps ) goto 20

C       Scale the root to the desired interval, and put in its symmetric
C       counterpart. Compute the weight and its symmetric counterpart.

        grd(i) = xm - xl * z
        grd(n+1-i) = xm + xl * z

        wgt(i) = 2.0d0 * xl / ((1.0d0 - z*z) * pp * pp)
        wgt(n+1-i) = wgt(i)
        

 10   continue

      return
      end

