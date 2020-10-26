

      subroutine trprct (nri, n, m, a, b)

      implicit none


C +------------------------------------------------------------+
C |  trprct   --   CAJH, 07.2011                               |
C |                                                            |
C |                                                            |
C |  Perform the matrix operation                              |
C |                                                            |
C |    B  =  A^T                                               |
C |                                                            |
C |  for a matrix A of dimensions N x M, possibly complex.     |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       nri - = 1 for real matrices, = 2 for complex
C       n   - leading dimension of matrix a
C       m   - leading dimension of matrix b
C       a   - input matrix
C       b   - output matrix  [ out ]

      integer     nri, n, m
      real*8      a(nri,n,m), b(nri,m,n)

C     other variables

      integer     i, j, k


      do 11 i = 1, n
        do 12 j = 1, m
          do 13 k = 1, nri
            b(k,j,i) = a(k,i,j)
 13       continue
 12     continue
 11   continue


      return
      end
      
