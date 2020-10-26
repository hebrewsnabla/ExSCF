

      double precision function wignerd (isint, jinp, minp, ninp, beta)

      implicit none


C +------------------------------------------------------------+
C |  wignerd   --   CAJH, 07.2011                              |
C |                                                            |
C |                                                            |
C |  Compute matrix elements of Wigner's d-matrix:             |
C |                                                            |
C |    wignerd (j, m, n, beta)  ==  d^j_{m,n} (beta)           |
C |                                                            |
C |  for integer or half-integer j, m, and n. j is required    |
C |  to be non-negative. m and n should satisfy                |
C |    -j <= m, n <= +j                                        |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  The logical parameter isint should be used to tell        |
C |  whether the indices j, m, and n, are integer or           |
C |  half-integer. That is,                                    |
C |                                                            |
C |    wignerd (.true., 3, 1, -1, beta)                        |
C |        == d^3_{1,-1} (beta),                               |
C |                                                            |
C |    wignerd (.false., 3, 1, -1, beta)                       |
C |        == d^(3/2)_{1/2,-1/2} (beta).                       |
C |                                                            |
C |  Note that if isint = .false., then all indices (j,m,n)    |
C |  should be odd numbers.                                    |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  The evaluation of the matrix elements is done using       |
C |  explicit expressions for j <= 2. For j > 2, we use        |
C |  the following expression:                                 |
C |       [ from Lai, et al., J. Math. Chem. 19, 131 (1996) ]  |
C |                                                            |
C |    d^j_{m,n} (beta) =                                      |
C |        [ (j+m)! (j-m)! / { (j+n)! (j-n)! } ]^(1/2) *       |
C |        sum_s { (-1)^(m-n+s) * binom (j+n,s) *              |
C |                binom (j-n,j-m-s) *                         |
C |                cos(beta/2)^(2j+n-m-2s) *                   |
C |                sin(beta/2)^(m-n+2s) },                     |
C |                                                            |
C |  where binom (j, k) is the binomial coefficient and the    |
C |  sum over s is restricted to the argument of any           |
C |  factorial which is non-negative.                          |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       isint - whether the j coefficient is integer (.true.) or
C               half-integer (.false.)

      real*8      beta
      logical     isint
      integer     jinp, minp, ninp

C     other variables

      real*8      z, rt1, rt2, sums, cosbt, sinbt, fac1, fac2, sgnfac
      integer     j, m, n
      integer     j1, j2, maxs, s, exp1, exp2, difmn
      logical     tstm, tstn, ltest

C     functions

      real*8      rtfact, binom

C     constants

      real*8      f12

      parameter   ( f12 = 1.0d0 / 2.0d0 )


      j = jinp
      m = minp
      n = ninp

C     Error checking.
C     Verify that the values of j, m, n provided are physical.

      if ( j .lt. 0 ) then
        call gauerr ('Wrong value of j provided to wignerd.')
      endif

      tstm = abs (m) .le. j
      tstn = abs (n) .le. j

      if ( .not. tstm .or. .not. tstn ) then
        call gauerr ('Unphysical values for m, n in wignerd.')
      endif

      if ( .not. isint ) then
        if ( mod (abs(j),2) .eq. 0 .or.
     $       mod (abs(m),2) .eq. 0 .or.
     $       mod (abs(n),2) .eq. 0 ) then
          call gauerr ('Inconsistent isint with j,m,n in wignerd.')
        end if
      end if


      ltest = ( isint .and. mod(n-m,2) .eq. 1 ) .or.
     $        ( .not. isint .and. mod((n-m)/2,2) .eq. 1 )


C     We do not need to code all the matrix elements. We make use of the
C     following relations between matrix elements:

C       d^j_{m,n} (beta)  =  (-1)^(m-n) * d^j_{n,m} (beta)
C                         =  (-1)^(m-n) * d^j_{-m,-n} (beta)

C     Thus, we only need to code the case m >= 0, m >= n.

C       Here, we define the variable sgnfac that keeps track of the sign.

      sgnfac = 1.0d0

      if ( m .ge. 0 ) then

        if ( m .ge. abs(n) ) then

C         nothing to do

        else if ( m .lt. n ) then
          m = ninp
          n = minp

          if ( ltest ) sgnfac = -sgnfac

        else if ( m .gt. n ) then
          m = -ninp
          n = -minp
        end if

      else

        if ( -m .ge. abs(n) ) then
          m = -minp
          n = -ninp

          if ( ltest ) sgnfac = -sgnfac

        else if ( m .lt. n ) then
          m = ninp
          n = minp

          if ( ltest ) sgnfac = -sgnfac

        else if ( m .gt. n ) then
          m = -ninp
          n = -minp
        end if

      end if


C     We use explicit expressions for the simplest cases (j <= 2). Such
C     expressions are available in, for instance, Wikipedia.

C     Here and in the following, we adopt the notation
C       d (j, m, n)  =  d^j_{m,n} (beta).

      z = 1.0d0


C     case j = 1/2

C       d (1/2, 1/2,  1/2) =  cos(beta/2)
C       d (1/2, 1/2, -1/2) = -sin(beta/2)

      if ( j .eq. 1 .and. .not. isint ) then
        if ( n .eq. 1 ) then
          z = cos (beta / 2.0d0)
        elseif ( n .eq. -1 ) then
          z = -sin (beta / 2.0d0)
        endif
      endif


C     case j = 1

C       d (1, 1,  1) =  (1 + cos(beta)) / 2
C       d (1, 1,  0) = -sin(beta) / sqrt(2)
C       d (1, 1, -1) =  (1 - cos(beta)) / 2
C       d (1, 0,  0) =  cos(beta)

      if ( j .eq. 1 .and. isint ) then
        if ( m .eq. 1 ) then
          if ( n .eq. 1 ) then
            z = f12 * (1.0d0 + cos (beta))
          elseif ( n .eq. 0 ) then
            z = -sin (beta) / sqrt (2.0d0)
          elseif ( n .eq. -1 ) then
            z = f12 * (1.0d0 - cos (beta))
          endif

        elseif ( m .eq. 0 ) then
          z = cos (beta)
        endif
      endif


C     case j = 3/2

C       d (3/2, 3/2,  3/2) =  1/2 * (1 + cos(beta)) * cos(beta/2)
C       d (3/2, 3/2,  1/2) = -sqrt(3)/2 * (1 + cos(beta)) * sin(beta/2)
C       d (3/2, 3/2, -1/2) =  sqrt(3)/2 * (1 - cos(beta)) * cos(beta/2)
C       d (3/2, 3/2, -3/2) = -1/2 * (1 - cos(beta)) * sin(beta/2)
C       d (3/2, 1/2,  1/2) =  1/2 * (3*cos(beta) - 1) * cos(beta/2)
C       d (3/2, 1/2, -1/2) = -1/2 * (3*cos(beta) + 1) * sin(beta/2)

      if ( j .eq. 3 .and. .not. isint ) then
        if ( m .eq. 3 ) then
          if ( n .eq. 3 ) then
            z = f12 * (1.0d0 + cos (beta)) * cos (beta / 2.0d0)
          elseif ( n .eq. 1 ) then
            z = -f12 * sqrt(3.0d0) * (1.0d0 + cos (beta)) *
     $                               sin (beta / 2.0d0)
          elseif ( n .eq. -1 ) then
            z = f12 * sqrt(3.0d0) * (1.0d0 - cos (beta)) *
     $                              cos (beta / 2.0d0)
          elseif ( n .eq. -3 ) then
            z = -f12 * (1.0d0 - cos (beta)) * sin (beta / 2.0d0)
          endif

        elseif ( m .eq. 1 ) then
          if ( n .eq. 1 ) then
            z = f12 * (3.0d0 * cos (beta) - 1.0d0) *
     $                cos (beta / 2.0d0)
          elseif ( n .eq. -1 ) then
            z = -f12 * (3.0d0 * cos (beta) + 1.0d0) *
     $                 sin (beta / 2.0d0)
          endif
        endif
      endif


C     case j = 2

C       d (2, 2,  2) =  ((1 + cos(beta)) / 2)^2
C       d (2, 2,  1) = -1/2 * (1 + cos(beta)) * sin(beta)
C       d (2, 2,  0) =  sqrt(6)/4 * (sin(beta))^2
C       d (2, 2, -1) = -1/2 * (1 - cos(beta)) * sin(beta)
C       d (2, 2, -2) =  ((1 - cos(beta)) / 2)^2
C       d (2, 1,  1) =  1/2 * (1 + cos(beta)) * (2*cos(beta) - 1)
C       d (2, 1,  0) = -sqrt(3/2) * sin(beta) * cos(beta)
C       d (2, 1, -1) =  1/2 * (1 - cos(beta)) * (2*cos(beta) + 1)
C       d (2, 0,  0) =  1/2 * (3*(cos(beta))^2 - 1)

      if ( j .eq. 2 .and. isint ) then
        if ( m .eq. 2 ) then
          if ( n .eq. 2 ) then
            z = (1.0d0 + cos (beta))**2 / 4.0d0
          elseif ( n .eq. 1 ) then
            z = -f12 * (1.0d0 + cos (beta)) * sin (beta)
          elseif ( n .eq. 0 ) then
            z = sqrt (6.0d0) / 4.0d0 * (sin (beta))**2
          elseif ( n .eq. -1 ) then
            z = -f12 * (1.0d0 - cos (beta)) * sin (beta)
          elseif ( n .eq. -2 ) then
            z = (1.0d0 - cos (beta))**2 / 4.0d0
          endif

        elseif ( m .eq. 1 ) then
          if ( n .eq. 1 ) then
            z = f12 * (1.0d0 + cos (beta)) *
     $                (2.0d0 * cos (beta) - 1.0d0)
          elseif ( n .eq. 0 ) then
            z = -sqrt (3.0d0 / 2.0d0) * sin (beta) * cos (beta)
          elseif ( n .eq. -1 ) then
            z = f12 * (1.0d0 - cos (beta)) *
     $                (2.0d0 * cos (beta) + 1.0d0)
          endif

        elseif ( m .eq. 0 ) then
          z = f12 * (3.0d0 * (cos (beta))**2 - 1.0d0)
        endif
      endif


C     all other cases, j > 2

C     In this case we compute the matrix element using the following
C     expression (from Lai, et al., J. Math. Chem. 19, 131 (1996))

C       d (j,m,n,beta) =  [ (j+m)! (j-m)! / { (j+n)! (j-n)! } ]^(1/2) *
C                         sum_s { (-1)^(m-n+s) * binom (j+n,s) *
C                                   binom (j-n,j-m-s) *
C                                   cos(beta/2)^(2j+n-m-2s) *
C                                   sin(beta/2)^(m-n+2s) },

C     where binom (j, k) is the binomial coefficient and the sum over s
C     is restricted to the argument of any factorial which is
C     non-negative.


      if ( ( j .gt. 3 .and. .not. isint ) .or.
     $     ( j .gt. 2 .and. isint ) ) then


C     Compute prefactor

C       [ (j+m)! (j-m)! / { (j+n)! (j-n)! } ]^(1/2).

C     Recall that m >= 0, so we know already that j+m > j-m. The same
C     does not need to be true for n.

        if ( n .ge. 0 ) then
          j1 = j+n
          j2 = j-n
        else
          j1 = j-n
          j2 = j+n
        endif

C       Compute the ratios of factorials
C         (j+m)! / j1!,  (j-m)! / j2!

        if ( isint ) then
          rt1 = rtfact (j+m, j1)
          rt2 = rtfact (j-m, j2)
        else
          rt1 = rtfact ((j+m)/2, j1/2)
          rt2 = rtfact ((j-m)/2, j2/2)
        endif

        z = sqrt (rt1 * rt2)

C       The prefactor is ready and loaded in z. Now we compute the sum
C       over s.


C     Loop over s to compute

C       sum  =  sum_s  (-1)^(m-n+s) *
C                           binom (j+n,s) * binom(j-n,j-m-s) *
C                           cos(beta/2)^(2j+n-m-2s) *
C                           sin(beta/2)^(m-n+2s),

C     where the sum over s is restricted to the argument of any
C     factorial which is non-negative.


C       Before getting into the actual computation, let us decide what
C       is the maximal value of s. Four conditions should be met:

C         s >= 0
C         s <= j+n
C         s <= j-m
C         n <= m+s   ==>   s >= n-m

C       These lead to  0 <= s <= min (j+n,j-m).

        sums = 0.0d0

        cosbt = cos (beta / 2.0d0)
        sinbt = sin (beta / 2.0d0)

        maxs = min (j+n, j-m)
        if ( .not. isint ) maxs = maxs / 2

        do 10 s = 0, maxs

          if ( isint ) then
            fac1 = binom (j+n, s)
            fac2 = binom (j-n, j-m-s)
          else
            fac1 = binom ((j+n)/2, s)
            fac2 = binom ((j-n)/2, (j-m)/2-s)
          endif

          exp1 = 2*j + n - m
          exp2 = m - n
          
          if ( .not. isint ) then
            exp1 = exp1 / 2
            exp2 = exp2 / 2
          endif

          if ( mod (exp2 + s, 2) .eq. 0 ) then
            sums = sums + fac1 * fac2 *
     $                    cosbt**(exp1 - 2*s) * sinbt**(exp2 + 2*s)
          else
            sums = sums - fac1 * fac2 *
     $                    cosbt**(exp1 - 2*s) * sinbt**(exp2 + 2*s)
          endif
 10     continue


C     Finish building z.

        z = z * sums
      endif


C     Account for the factor of (-1)^(m-n) if needed...

      z = sgnfac * z


C     Final result...

      wignerd = z


      return
      end


