

      double precision function rtfact (j, k)

      implicit none


C +------------------------------------------------------------+
C |  rtfact   --   CAJH, 07.2011                               |
C |                                                            |
C |                                                            |
C |  Compute the ratio of two factorial numbers:               |
C |                                                            |
C |    rtfact (j, k) =  j! / k!,                               |
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
      integer     r


C     Error checking.

      if ( j .lt. 0 .or. k .lt. 0 ) then
        call gauerr ('Incorrect j, k in rtfact.')
      endif


      if ( j .eq. 0 ) j = 1
      if ( k .eq. 0 ) k = 1


C     If j == k, then the ratio is trivially 1.

      if ( j .eq. k ) then
        rtfact = 1.0d0
        return
      endif


C     We consider two cases: j > k, and j < k.
C     We handle each case independently, to avoid overflow or underflow
C     as much as possible.

      z = 1.0d0

      if ( j .gt. k ) then
        do 10 r = j, k+1, -1
          z = z * dble (r)
 10     continue

      else
        do 20 r = k, j+1, -1
          z = z / dble (r)
 20     continue
      endif

      rtfact = z


      return
      end



