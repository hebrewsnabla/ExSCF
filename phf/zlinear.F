*Deck ZLinear
      Subroutine ZLinear(A,B,Max,N)
      Implicit Real*8(A-H,O-Z)
C
C     Places symmetric square array in linear form.
C
      Complex*16 A(Max,Max), B(*)
C
      Do 10 J = 1, N
        JJ = (J*(J-1))/2
        Do 10 I = 1, J
   10     B(JJ+I) = A(I,J)
      Return
      End
