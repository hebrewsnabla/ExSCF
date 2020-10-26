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
*  DSYEVD Example.
*  ==============
*
*  Program computes all eigenvalues and eigenvectors of a real symmetric
*  matrix A using divide and conquer algorithm, where A is:
*
*    6.39   0.13  -8.23   5.71  -3.18
*    0.13   8.37  -4.46  -6.10   7.21
*   -8.23  -4.46  -9.58  -9.25  -7.42
*    5.71  -6.10  -9.25   3.72   8.54
*   -3.18   7.21  -7.42   8.54   2.51
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
*  If the eigenvectors are requested, then this routine uses a divide and
*  conquer algorithm to compute eigenvalues and eigenvectors.
*
*  Example Program Results.
*  ========================
*
* DSYEVD Example Program Results
* 
* Eigenvalues
* -17.44 -11.96   6.72  14.25  19.84
* 
* Eigenvectors (stored columnwise)
*  -0.26   0.31  -0.74   0.33   0.42
*  -0.17  -0.39  -0.38  -0.80   0.16
*  -0.89   0.04   0.09   0.03  -0.45
*  -0.29  -0.59   0.34   0.31   0.60
*  -0.19   0.63   0.44  -0.38   0.48
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
      INTEGER          INFO, LWORK, LIWORK
*
*     .. Local Arrays ..
      INTEGER          IWORK( LWMAX )
      DOUBLE PRECISION A( LDA, N ), W( N ), WORK( LWMAX )
      DATA             A/
     $ 0.642742149686186, 0.0d0, 0.0d0, 0.0d0,
     $-0.465821333149620,  0.337599633886206, 0.0d0, 0.0d0,
     $-0.066384044853257,  0.0481112126978231, 
     $0.0068563130848518, 0.0d0,
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
      EXTERNAL         DSYEVD
      EXTERNAL         PRINT_MATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'DSYEVD Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      LIWORK = -1
      CALL DSYEVD( 'Vectors', 'Upper', N, A, LDA, W, WORK, LWORK,
     $             IWORK, LIWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      LIWORK = MIN( LWMAX, IWORK( 1 ) )
*
*     Solve eigenproblem.
*
      CALL DSYEVD( 'Vectors', 'Upper', N, A, LDA, W, WORK, LWORK,
     $             IWORK, LIWORK, INFO )
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
*     End of DSYEVD Example.
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
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( 11(:,1X,E20.10) )
      RETURN
      END
