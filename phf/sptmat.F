

      subroutine sptmat (nri, m1, m2, a, n1, n2, b, key)


C +------------------------------------------------------------+
C |  sptmat   --   CAJH, 07.2011                               |
C |                                                            |
C |                                                            |
C |  Extract an N1 x N2 matrix B from a larger M1 x M2 matrix  |
C |  A or do the opposite: paste a small N1 x N2 matrix B      |
C |  into a large M1 x M2 matrix A (with zeros everywhere      |
C |  else).                                                    |
C |                                                            |
C |    key = 1,  extract B from A                              |
C |    key = 2,  paste B into A                                |
C |                                                            |
C |  In particular, we always cut or paste B from the top-     |
C |  left corner of A. For instance, if A is given by          |
C |                                                            |
C |    A =  (  a11  a12  )  N1                                 |
C |         (  a21  a22  )  M1-N1                              |
C |            N2   M2-N2                                      |
C |                                                            |
C |  where we have specified the dimensions of the sub-        |
C |  matrices, then B = a11 if key = 1 is used.                |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       nri    - (= 1 for real matrices, = 2 for complex matrices)
C       m1, m2 - dimensions of a: a(m1,m2)
C       n1, n2 - dimensions of b: b(n1,n2)
C       key    - whether to extract from matrix or build the full matrix

      integer     nri, m1, m2, n1, n2
      real*8      a(nri,m1,m2), b(nri,n1,n2)

C     other variables

      integer     i, j, z


C     Error checking.

      if ( n1 .gt. m1 .or. n2 .gt. m2 ) then
        call gauerr ('Wrong dimensions in sptmat.')
      endif


C     Split code.

      if ( key .eq. 1 ) goto 10
      if ( key .eq. 2 ) goto 20
      call gauerr ('key not supported in sptmat')



C     %%%%%%%%%%%%%%%%%%%%%%
C     %  Extract B from A  %
C     %%%%%%%%%%%%%%%%%%%%%%


 10   continue

      call aclear (nri*n1*n2, b)

      do 101 i = 1, n2
        do 102 j = 1, n1
          do 103 z = 1, nri
            b(z,j,i) = a(z,j,i)
 103      continue
 102    continue
 101  continue

      return



C     %%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Build full A from B  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%


 20   continue

      call aclear (nri*m1*m2, a)

      do 201 i = 1, n2
        do 202 j = 1, n1
          do 203 z = 1, nri
            a(z,j,i) = b(z,j,i)
 203      continue
 202    continue
 201  continue


      return
      end


