

      subroutine prtocc (iout, name1, name2, name3, n, vec)

      implicit none


C +------------------------------------------------------------+
C |  prtocc   --   CAJH, 07.2011                               |
C |                                                            |
C |                                    (based on PHFB prtocc)  |
C |                                                            |
C |                                                            |
C |  Print a vector of real eigenvalues in a nice format.      |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       name1 - name of matrix whose eigenvalues are being printed
C       name2 - typically it should be just 'orbital'
C       name3 - name of the eigenvalues, i.e. 'occupation' or
C               'coefficient'

      integer     iout, n
      character   name1*(*), name2*(*), name3*(*)
      real*8      vec(n)


C     other variables

      real*8      prtvec(4), thresh
      integer     i, j

      parameter   ( thresh = 1.0d-14 )


 101  format (' ',A12,' -- ',1(14X,I3))
 102  format (' ',A12,' -- ',2(14X,I3))
 103  format (' ',A12,' -- ',3(14X,I3))
 104  format (' ',A12,' -- ',4(14X,I3))
 111  format (' ',A12,' -- ',1(3X,1PD14.7))
 112  format (' ',A12,' -- ',2(3X,1PD14.7))
 113  format (' ',A12,' -- ',3(3X,1PD14.7))
 114  format (' ',A12,' -- ',4(3X,1PD14.7))


      write (iout, *)
      write (iout, *) name1, ':'


      do 100 i = 1, n, 4

        do 300 j = 1, 4
          prtvec(j) = 0.0d0

          if ( ( i + j - 1 .le. n ) .and.
     $         ( abs (vec(i+j-1)) .ge. thresh ) ) then
            prtvec(j) = vec(i+j-1)
          endif
 300    continue

        if ( i + 3 .le. n ) then

          write (iout, 104) name2, i, i+1, i+2, i+3
          write (iout, 114) name3, prtvec(1), prtvec(2), prtvec(3),
     $                      prtvec(4)

          if ( i + 3 .eq. n ) goto 200

        elseif ( i + 2 .le. n ) then

          write (iout, 103) name2, i, i+1, i+2
          write (iout, 113) name3, prtvec(1), prtvec(2), prtvec(3)

          if ( i + 2 .eq. n ) goto 200

        elseif ( i + 1 .le. n ) then

          write (iout, 102) name2, i, i+1
          write (iout, 112) name3, prtvec(1), prtvec(2)

          if ( i + 1 .eq. n ) goto 200

        else

          write (iout, 101) name2, i
          write (iout, 111) name3, prtvec(1)

          goto 200
        endif

 100  continue
 200  continue


      return
      end


