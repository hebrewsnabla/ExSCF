*Deck TrcABC
      Complex*16 Function TrcABC(NMat,N,A,B)
      Implicit Real*8(A-H,O-Z)
C
C     Trace(A*B) with complex square matrices.
C
      Complex*16 GCmplx
      Dimension A(2,N,N,NMat), B(2,N,N,NMat)
      Save Zero
      Data Zero/0.0d0/
C
      SR = Zero
      SI = Zero
      Do 10 IM = 1, NMat
        Do 10 I = 1, N
          Do 10 J = 1, N
            SR = SR + A(1,I,J,IM)*B(1,J,I,IM)
     $        + A(2,I,J,IM)*B(2,J,I,IM)
   10       SI = SI + A(2,I,J,IM)*B(1,J,I,IM)
     $        - A(1,I,J,IM)*B(2,J,I,IM)
      TrcABC = GCmplx(SR,SI)
      Return
      End
