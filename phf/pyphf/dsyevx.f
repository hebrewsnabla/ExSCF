!===============================================================================
! Copyright 2009-2019 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

*
*  DSYEV Example.
*  ==============
*
*  Program computes all eigenvalues and eigenvectors of a real symmetric
*  matrix A:
*
*    1.96  -6.49  -0.47  -7.20  -0.65
*   -6.49   3.80  -6.39   1.50  -6.34
*   -0.47  -6.39   4.17  -1.51   2.67
*   -7.20   1.50  -1.51   5.70   1.80
*   -0.65  -6.34   2.67   1.80  -7.10
*
*  Description.
*  ============
*
*  The routine computes all eigenvalues and, optionally, eigenvectors of an
*  n-by-n real symmetric matrix A. The eigenvector v(j) of A satisfies
*
*  A*v(j) = lambda(j)*v(j)
*
*  where lambda(j) is its eigenvalue. The computed eigenvectors are
*  orthonormal.
*
*  Example Program Results.
*  ========================
*
* DSYEV Example Program Results
* 
* Eigenvalues
* -11.07  -6.23   0.86   8.87  16.09
* 
* Eigenvectors (stored columnwise)
*  -0.30  -0.61   0.40  -0.37   0.49
*  -0.51  -0.29  -0.41  -0.36  -0.61
*  -0.08  -0.38  -0.66   0.50   0.40
*   0.00  -0.45   0.46   0.62  -0.46
*  -0.80   0.45   0.17   0.31   0.16
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 4 )
      INTEGER          LDA
      PARAMETER        ( LDA = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK
*
*     .. Local Arrays ..
      DOUBLE PRECISION A( LDA, N ), W( N ), WORK( LWMAX )
      DATA             A/
     $ 0.6427421496861863, -0.4658213331496202, 
     $-0.0663840448532566,  0.0907113787628938,
     $-0.4658213331496202,  0.3375996338862058, 
     $ 0.0481112126978231, -0.0657422193453504,
     $-0.0663840448532566,  0.0481112126978231, 
     $ 0.0068563130848518, -0.0093689020386119,
     $ 0.0907113787628938, -0.0657422193453504, 
     $-0.0093689020386119,  0.0128022633043168
c     $ 0.64274215, -0.46582133, -0.06638404,  0.09071138,
c     $-0.46582133,  0.33759963,  0.04811121, -0.06574222,
c     $-0.06638404,  0.04811121,  0.00685631, -0.0093689 ,
c     $ 0.09071138, -0.06574222, -0.0093689 ,  0.01280226
c     $  1.96, 0.00, 0.00, 0.00, 0.00,
c     $ -6.49, 3.80, 0.00, 0.00, 0.00,
c     $ -0.47,-6.39, 4.17, 0.00, 0.00,
c     $ -7.20, 1.50,-1.51, 5.70, 0.00,
c     $ -0.65,-6.34, 2.67, 1.80,-7.10
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         DSYEV
      EXTERNAL         PRINT_MATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'DSYEV Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL DSYEV( 'Vectors', 'Upper', N, A, LDA, W, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve eigenproblem.
*
      CALL DSYEV( 'Vectors', 'Upper', N, A, LDA, W, WORK, LWORK, INFO )
*
*     Check for convergence.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF
*
*     Print eigenvalues.
*
      CALL PRINT_MATRIX( 'Eigenvalues', 1, N, W, 1 )
*
*     Print eigenvectors.
*
      CALL PRINT_MATRIX( 'Eigenvectors (stored columnwise)', N, N, A, 
     $                   LDA )
      STOP
      END
*
*     End of DSYEV Example.
*






*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      DOUBLE PRECISION A( LDA, * )
*
      INTEGER          I, J
*
*      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( 11(:,1X,E20.10) )
      RETURN
      END
