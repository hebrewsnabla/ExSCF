*Deck VConjg
      Subroutine VConjg(N,A)
      Implicit Real*8(A-H,O-Z)
C
C     Complex conjugate of a vector.
C
      Dimension A(2,N)
C
      Do 10 I = 1, N
   10   A(2,I) = -A(2,I)
      Return
      End
