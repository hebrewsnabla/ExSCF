*Deck SymAdS
      Subroutine SymAdS(N,NB,IndI,MaxCmp,NSABF,S,B)
      Implicit Real*8(A-H,O-Z)
C
C     Transform S into one symmetry block using integer symmetry vectors.
C
      Dimension NSABF(N,*), S(*), B(NB,NB)
      Save One
      Data One/1.0d0/
C
      RNSq = One / GFloat(MaxCmp)
      Call AClear(NB**2,B)
C$OMP Parallel Do Default(Shared)
C$OMP+ Private(I,IC,IX,SignI,J,JC,JX,SignIJ,IJX)
      Do 30 I = 1, NB
        Do 30 IC = 1, MaxCmp
          IX = NSABF(I+IndI,IC)
          If(IX.gt.0) then
            SignI = RNSq
          else
            IX = -IX
            SignI = -RNSq
            endIf
          If(IX.ne.0) then
            Do 20 J = 1, I
              Do 20 JC = 1, MaxCmp
                JX = NSABF(J+IndI,JC)
                If(JX.gt.0) then
                  SignIJ = SignI
                else
                  JX = -JX
                  SignIJ = -SignI
                  endIf
                If(IX.ge.JX) then
                  IJX = (IX*(IX-1))/2 + JX
                else
                  IJX = (JX*(JX-1))/2 + IX
                  endIf
   20           B(J,I) = B(J,I) + S(IJX)*SignIJ
            endIf
   30     Continue
      Do 40 J = 1, (NB-1)
        Do 40 I = (J+1), NB
   40     B(I,J) = B(J,I)
      Return
      End
