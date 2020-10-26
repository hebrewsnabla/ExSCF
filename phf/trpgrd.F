

      subroutine trpgrd (x1, x2, grd, wgt, n)

      implicit none


C +------------------------------------------------------------+
C |  trpgrd   --   CAJH, 07.2011                               |
C |                                                            |
C |                                    (based on PHFB trpgrd)  |
C |                                                            |
C |                                                            |
C |  Fills the grid positions and weights to perform a         |
C |  trapezoidal-based integration over an integral defined    |
C |  by the x1 and x2 limits using n points.                   |
C |                                                            |
C |  Note that the weights returned are normalized, such that  |
C |                                                            |
C |    int_{x1}^{x2} dx f(x)  ==  1                            |
C |                                                            |
C |  when f(x) = 1.                                            |
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

      real*8      spc, xl
      integer     i


C     Determine initial position and spacing.

      spc = (x2 - x1) / n
      xl = x1 + spc / 2


C     Fill-in array positions and weights.

      do 10 i = 1, n
        grd(i) = xl + dble (i-1) * spc
        wgt(i) = 1.0d0 / dble (n)
 10   continue


      return
      end


