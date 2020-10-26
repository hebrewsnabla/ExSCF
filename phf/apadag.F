*Deck APADag
      Subroutine APADag(ISign,N,Alpha,A,B)
      Implicit Real*8(A-H,O-Z)
C
C     B = alpha * ( A +/- A^H )
C     ISign =  1 --> +
C           = -1 --> -
C     A is complex Hermitian.
C     If ISign =  1
C        B is complex Hermitian.
C     If ISign = -1
C        B is complex anit-Hermitian.
C
      Real*8 A(2,N,N), B(2,N,N)
      Save Zero
      Data Zero/0.0d0/
C
      TAlpha = Alpha + Alpha
      If(ISign.eq.1) then
        Do 10 I = 1, N
          B(1,I,I) = TAlpha*A(1,I,I)
          B(2,I,I) = Zero
          Do 10 J = 1, (I-1)
            B(1,I,J) = Alpha*(A(1,I,J)+A(1,J,I))
            B(2,I,J) = Alpha*(A(2,I,J)-A(2,J,I))
            B(1,J,I) = B(1,I,J)
   10       B(2,J,I) = -B(2,I,J)
      else if(ISign.eq.-1) then
        Do 20 I = 1, N
          B(1,I,I) = Zero
          B(2,I,I) = TAlpha*A(2,I,I)
          Do 20 J = 1, (I-1)
            B(1,I,J) = Alpha*(A(1,I,J)-A(1,J,I))
            B(2,I,J) = Alpha*(A(2,I,J)+A(2,J,I))
            B(1,J,I) = -B(1,I,J)
   20       B(2,J,I) = B(2,I,J)
      else
        Call GauErr('Incorrect Choice of ISign in APADag')
        endIf
      Return
      End
