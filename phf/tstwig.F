

      program tstwig

      implicit none

C     tstwig - A simple program to test the subroutine wignerd.
C              See wignerd for more details...

C     CAJH, 07.03.2011


      real*8      beta, z, pi, wignerd
      integer     j, m, n, l
      logical     isint

 21   format (2X, 'd (', 3(I4,','), F7.2, ') = ', 1P, E24.12E4)
 22   format (2X, 'd (', 3(I4,'/2,'), F7.2, ') = ', 1P, E24.12E4)


      pi = 4.0d0 * atan (1.0d0)


      write (*,*)
      write (*,*) 'Testing funtion wignerd'
      write (*,*)

      write (*,*) 'Given three integers, j, m, and n, the program will'
      write (*,*) 'print a table with values of the function'
      write (*,*) '  d (j, m, n, beta)'
      write (*,*) 'for beta = 0, 10, ..., ,180.'
      write (*,*)

      write (*,*) 'Please enter the value of j:'
      write (*,*) '  (use positive numbers for integer j, '
      write (*,*) '   negative numbers for half-integer j)'
      read (*,*) j

      if ( j .lt. 0 .and. mod (j,2) .eq. 0 ) then
        write (*,*) 'ERROR: incorrect value of j.'
        stop
      elseif ( abs (j) .gt. 15 ) then
        write (*,*) 'ERROR: selected value of j is too large.'
        stop
      endif

      isint = .true.

      if ( j .lt. 0 ) then
        j = -j
        isint = .false.
      endif


      write (*,*) 'Enter the value of m:'
      read (*,*) m

      write (*,*) 'Enter the value of n:'
      read (*,*) n

      if ( abs (m) .gt. j .or. abs (n) .gt. j ) then
        write (*,*) 'ERROR: incorrect value of m or n.'
        stop
      endif


      write (*,*)

      do 10 l = 0, 18, 1
        beta = dble (l*10) * pi / 180.0d0

        z = wignerd (isint, j, m, n, beta)

        if ( isint ) then
          write (*,21) j, m, n, dble (l*10), z
        else
          write (*,22) j, m, n, dble (l*10), z
        endif
 10   continue


      write (*,*)


      stop
      end


