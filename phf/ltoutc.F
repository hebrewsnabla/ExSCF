*Deck LTOutC
      Subroutine LTOutC(IOut,IOpt,NRI,IRI,NSpin,N,A,Key)
      Implicit Real*8(A-H,O-Z)
C
C     Working precision routine to print out the lower triangular part
C     of a matrix:
C
C     IOpt .. Output format:
C             1 ... D13.6
C             2 ... F13.3
C     NRI ... Number of components (1 for real, 2 for complex).
C     IRI ... Which component to print here (1=real, 2=imaginary).
C             Negative for anti-Hermetian matrices.
C     NSpin . -1 for just large component, 1 for just one spin component,
C             -2 for l and s components, 2 for a and b components.
C     N   ... Dimension of matrix.
C     A   ... Array to be printed.
C     Key ... <0: Suppress elements with absolute values less than
C                 10**(-6+Key)
C              1: Print complete matrix.
C
      Parameter (NumCol=5)
      Character*1 LabX(NumCol)
      Dimension A(NRI,*), S(NumCol), IRX(NumCol)
      Save Zero
      Data Zero/0.0d0/
 1010 Format(4X,10(7X,I6,A1))
 1020 Format(I7,A,D23.15,9D23.15)
 1021 Format(I7,A,D13.6,9D14.6)
 1030 Format(I7,A,F13.3,9F14.3)
C
      If(IOpt.lt.1.or.IOpt.gt.2) Call GauErr('Illegal IOpt in LTOutC.')
      Call PrtThr(Key,Thresh)
      NDo = IAbs(NSpin)*N
      IRIA = IAbs(IRI)
      Do 100 IStart = 1, NDo, NumCol
        IEnd = Min(IStart+NumCol-1,NDo)
        NCol = IEnd - IStart + 1
        Call LabSpn(NSpin,IStart-1,NCol,IRX,LabX)
        Write(IOut,1010) (IRX(I),LabX(I),I=1,NCol)
        Do 100 IRow = IStart, NDo
          ILim = Min(IRow-IStart+1,NumCol)
          L = (IRow*(IRow-1))/2 + IStart
          Do 10 I = 1, ILim
            S(I) = A(IRIA,L+I-1)
            If(Key.ne.1.and.Abs(S(I)).lt.Thresh) S(I) = Zero
   10       Continue
          If(NRI.eq.2.and.(IRI.eq.-1.or.IRI.eq.2)) then
            ILim1 = Min(ILim,IRow-IStart)
            Do 20 I = 1, ILim1
   20         S(I) = -S(I)
            endIf
          Call LabSpn(NSpin,IRow-1,1,IRX,LabX)
          If(IOpt.eq.1) then
            Write(IOut,1020) IRX(1), LabX(1), (S(I),I=1,ILim)
          else
            Write(IOut,1030) IRX(1), LabX(1), (S(I),I=1,ILim)
            endIf
  100     Continue
      Return
      End
