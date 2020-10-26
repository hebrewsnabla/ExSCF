

      double precision function binom (j, k)

      implicit none


C +------------------------------------------------------------+
C |  binom   --   CAJH, 07.2011                                |
C |                                                            |
C |                                                            |
C |  Compute the binomial coefficient:                         |
C |                                                            |
C |    binom (j, k) =  j! / [ k! (j-k)! ],                     |
C |                                                            |
C |  for integer, non-negative j and k.                        |
C |                                                            |
C |                                                            |
C |  Note that this function is not suitable for large j, k.   |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

      integer     j, k

C     other variables

      real*8      z
      integer     maxk, mink, r


C     Error checking.

      if ( k .gt. j ) then
        call gauerr ('k .gt. j in binom.')
      elseif ( j .lt. 0 .or. k .lt. 0 ) then
        call gauerr ('Incorrect j, k in binom.')
      endif


C     A few trivial cases, where binom (j, k) = 1.

      if ( j .eq. k .or. j .le. 1 .or. k .eq. 0 ) then
        binom = 1.0d0
        return
      endif


C     Determine the maximum and minimum between k and j-k.

      maxk = max (k, j-k)
      mink = min (k, j-k)


C     Compute z = j! / maxk!.

      z = 1.0d0

      do 10 r = j, maxk+1, -1
        z = z * dble (r)
 10   continue


C     Now compute z = z / mink!.

      do 20 r = 2, mink
        z = z / dble (r)
 20   continue

      binom = z


      return
      end


