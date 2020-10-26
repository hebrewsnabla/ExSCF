

      double precision function clgord (j1i, m1i, lint1, j2i, m2i,
     $     lint2, ji, mi, lint)

      implicit none


C +------------------------------------------------------------+
C |  clgord   --   CAJH, 12.2012                               |
C |                                                            |
C |                                                            |
C |  Compute Clebsch-Gordan coefficients:                      |
C |                                                            |
C |    <j1 m1, j2 m2 | j m >                                   |
C |                                                            |
C |  for integer or half-integer values of j,m. Only physical  |
C |  values of j-m pairs are allowed; j should be non-negative |
C |  and j >= abs(m).                                          |
C |                                                            |
C |  The logical parameters lint1, lint2, lint are used to     |
C |  determine whether the j-m pair should be interpreted as   |
C |  integer (.true.) or half-integer (.false.).               |
C |                                                            |
C |  Note that if lint = .false., then both j and m should be  |
C |  odd numbers.                                              |
C |                                                            |
C |  NOTE: This function should be used only for small values  |
C |    of angular momenta.                                     |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  The evaluation of the CG coefficient follows the formula  |
C |  provided in Wikipedia:                                    |
C |                                                            |
C |    <ji m1, j2 m2 | j m >  =  d_{m,m1+m2} *                 |
C |      sqrt[ (2*j+1) * (j+j1-j2)!*(j-j1+j2)!*                |
C |                      (j1+j2-j)!/(j1+j2+j+1)! ] *           |
C |      sqrt[ (j+m)!*(j-m)!*(j1+m1)!*(j1-m1)!*                |
C |                          (j2+m2)!*(j2-m2)! ] *             |
C |      sum_k (-1)^k *                                        |
C |         { (k)!*(j1+j2-j-k)!*(j1-m1-k)!*(j2+m2-k)!*         |
C |                (j-j2+m1+k)!*(j-j1-m2+k)! }.                |
C |                                                            |
C |  Here, the sum over k is extended over all integral k      |
C |  for which the argument of every factorial is nonnegative. |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  CAJH, 12.2012                                             |
C |                                                            |
C |  Switched to single integration.                           |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

      integer     j1i, m1i, j2i, m2i, ji, mi
      logical     lint1, lint2, lint

C     other variables

      real*8      fac1, fac2, sumf
      integer     jmp, jmm, j1m1p, j1m1m, j2m2p, j2m2m
      integer     ik, mink, maxk, sgnf
      integer     f1, f2, f3, f4, f5, f6
      integer     k1, k2, k3, k4, k5, k6
      integer     sf1, sf2, sf3, sf4, sf5, sf6
      integer     ssf1, ssf2, ssf3, ssf4, ssf5, ssf6

C     functions

      integer     cgfact


      clgord = 0.0d0

C     Look for negative j's.

      if ( (ji .lt. 0) .or. (j1i .lt. 0) .or. (j2i .lt. 0) ) then
        call gauerr ('Negative values of j in clgord.')
      end if

C     Check physical j-m relations.

      if ( abs(mi) .gt. ji ) then
        call gauerr ('Unphysical j,m values in clgord.')
      end if

      if ( abs(m1i) .gt. j1i ) then
        call gauerr ('Unphysical j1,m1 values in clgord.')
      end if

      if ( abs(m2i) .gt. j2i ) then
        call gauerr ('Unphysical j2,m2 values in clgord.')
      end if

C     Check integer / non-integer conditions.

      if ( .not. lint ) then
        if ( mod (abs(ji),2) .eq. 0 .or.
     $       mod (abs(mi),2) .eq. 0 ) then
          call gauerr ('Inconsistent lint with j,m in clgord.')
        end if
      end if

      if ( .not. lint1 ) then
        if ( mod (abs(j1i),2) .eq. 0 .or.
     $       mod (abs(m1i),2) .eq. 0 ) then
          call gauerr ('Inconsistent lint with j1,m1 in clgord.')
        end if
      end if

      if ( .not. lint2 ) then
        if ( mod (abs(j2i),2) .eq. 0 .or.
     $       mod (abs(m2i),2) .eq. 0 ) then
          call gauerr ('Inconsistent lint with j2,m2 in clgord.')
        end if
      end if

      if ( lint ) then
        if ( (lint1 .and. (.not. lint2)) .or.
     $       ((.not. lint1) .and. lint2) ) then
          call gauerr ('Inconsistent lint values in clgord.')
        end if

      else
        if ( (lint1 .and. lint2) .or.
     $       ((.not. lint1) .and. (.not. lint2)) ) then
          call gauerr ('Inconsistent lint values in clgord.')
        end if
      end if

C     Return if m1 + m2 - m /= 0: clgord = 0.

      if ( lint ) then
        if ( lint1 ) then
          k2 = m1i + m2i - mi
        else
          k2 = (m1i + m2i)/2 - mi
        end if

      else
        if ( lint1 ) then
          k2 = m1i + (m2i - mi)/2
        else
          k2 = m2i + (m1i - mi)/2
        end if
      end if

      if ( k2 .ne. 0 ) return


C     Generate a few quantities.

C       jmp = j + m;  jmm = j - m

      if ( lint ) then
        jmp = ji + mi
        jmm = ji - mi
      else
        jmp = (ji + mi)/2
        jmm = (ji - mi)/2
      end if

C       j1m1p = j1 + m1;  j1m1m = j1 - m1

      if ( lint1 ) then
        j1m1p = j1i + m1i
        j1m1m = j1i - m1i
      else
        j1m1p = (j1i + m1i)/2
        j1m1m = (j1i - m1i)/2
      end if

C       j2m2p = j2 + m2;  j2m2m = j2 - m2

      if ( lint2 ) then
        j2m2p = j2i + m2i
        j2m2m = j2i - m2i
      else
        j2m2p = (j2i + m2i)/2
        j2m2m = (j2i - m2i)/2
      end if

C       k1 = j1 + j2 + j
C       k2 = j1 + j2 - j
C       k3 = j + j2 - j1
C       k4 = j + j1 - j2
C       k5 = j - m2 - j1
C       k6 = j + m1 - j2

      if ( lint ) then
        if ( lint1 ) then
          k1 = j1i + j2i + ji
          k2 = j1i + j2i - ji
          k3 = ji + j2i - j1i
          k4 = ji + j1i - j2i
          k5 = ji - m2i - j1i
          k6 = ji + m1i - j2i
        else
          k1 = (j1i + j2i)/2 + ji
          k2 = (j1i + j2i)/2 - ji
          k3 = ji + (j2i - j1i)/2
          k4 = ji + (j1i - j2i)/2
          k5 = ji - (m2i + j1i)/2
          k6 = ji + (m1i - j2i)/2
        end if

      else
        if ( lint1 ) then
          k1 = j1i + (j2i + ji)/2
          k2 = j1i + (j2i - ji)/2
          k3 = (ji + j2i)/2 - j1i
          k4 = (ji - j2i)/2 + j1i
          k5 = (ji - m2i)/2 - j1i
          k6 = (ji - j2i)/2 + m1i
        else
          k1 = j2i + (j1i + ji)/2
          k2 = j2i + (j1i - ji)/2
          k3 = (ji - j1i)/2 + j2i
          k4 = (ji + j1i)/2 - j2i
          k5 = (ji - j1i)/2 - m2i
          k6 = (ji + m1i)/2 - j2i
        end if
      end if

C     Return if k2, k3 or k4 < 0: clgord = 0.

      if ( (k2 .lt. 0) .or. (k3 .lt. 0) .or. (k4 .lt. 0) ) return


C     Compute first prefactor.

      f1 = cgfact(k2)
      f2 = cgfact(k3)
      f3 = cgfact(k4)
      f4 = cgfact(k1+1)

      if ( lint ) then
        f5 = 2*ji + 1
      else
        f5 = ji + 1
      end if

      fac1 = sqrt (dble(f5)*dble(f1)*dble(f2)*dble(f3)/dble(f4))

C     Second prefactor.

      f1 = cgfact(jmp)
      f2 = cgfact(jmm)
      f3 = cgfact(j1m1p)
      f4 = cgfact(j1m1m)
      f5 = cgfact(j2m2p)
      f6 = cgfact(j2m2m)

      fac2 = sqrt(dble(f1)*dble(f2)*dble(f3)*dble(f4)*dble(f5)*dble(f6))

C     Decide lower and upper bounds in loop.

      mink = max (0, -k6, -k5)
      maxk = min (k2, j1m1m, j2m2p)

      sgnf = 1
      if ( mod(mink,2) .eq. 1 ) sgnf = -1


      sumf = 0.0d0

      do 100 ik = mink, maxk

        sf1 = ik
        sf2 = k2 - ik
        sf3 = j1m1m - ik
        sf4 = j2m2p - ik
        sf5 = k6 + ik
        sf6 = k5 + ik

        ssf1 = cgfact (sf1)
        ssf2 = cgfact (sf2)
        ssf3 = cgfact (sf3)
        ssf4 = cgfact (sf4)
        ssf5 = cgfact (sf5)
        ssf6 = cgfact (sf6)

        sumf = sumf +
     $    1.0d0*dble(sgnf)/dble(ssf1)/dble(ssf2)/dble(ssf3)/
     $                     dble(ssf4)/dble(ssf5)/dble(ssf6)

        sgnf = -sgnf
 100  continue

      clgord = sumf * fac1 * fac2


      return
      end



      integer function cgfact (j)

      implicit none

C +------------------------------------------------------------+
C |  cgfact   --   CAJH, 12.2012                               |
C |                                                            |
C |                                                            |
C |  Compute the factorial of j: j!.                           |
C |  Only small values of j (j <= 12) are allowed.             |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

      integer     j

C     other variables

      integer     k, z


C     Error checking.

      if ( j .lt. 0 ) then
        call gauerr ('Incorrect j in cgfact.')
      endif

      if ( j .gt. 12 ) then
        call gauerr ('j is too large in cgfact.')
      end if


C     A few trivial cases.

      if ( j .eq. 0 .or. j .eq. 1 ) then
        cgfact = 1
        return
      end if

C     Compute factorial.

      z = 1

      do 10 k = j, 2, -1
        z = z*k
 10   continue

      cgfact = z


      return
      end




