*Deck OutMtC
      Subroutine OutMtC(IOut,IOpt,NSpinr,Key,NRI,IRI,NSpin,X,M,N,MM,NN)
      Implicit Real*8(A-H,O-Z)
C
C     Print matrix X.  ISmall is >0 to print all elements, <=0 to
C     print elements greater than 10**(Key-6).  NRI is the number
C     of components (1 for real, 2 for complex), IRI which component
C     to print (1=real, 2=imaginary), NSpin is the number of spinor
C     terms (-1=just large, 1=just one spin, -2=l and s, 2=a and b).
C
C     IOpt .. Output format:
C             1 ... D14.6
C             2 ... F14.3
C             3 ... F20.8
C
      Parameter (NumCol=5)
      Logical NSpinr
      Character*1 LabX(NumCol)
      Dimension X(NRI,M,N), Y(NumCol), IRX(NumCol)
      Save Zero
      Data Zero/0.0d0/
 1000 Format(5X,8(7X,I6,A1))
 1010 Format(I7,A,8D23.15)
 1011 Format(I7,A,8D14.6)
 1020 Format(I7,A,8F14.3)
 1030 Format(I7,A,5F20.8)
 1040 Format(1X,8(13X,I6,A1))
C
      If(IOpt.lt.1.or.IOpt.gt.3) Call GauErr('Illegal IOpt in OutMtC.')
      Call PrtThr(Key,Small)
      MS = MM*IAbs(NSpin)
      If(NSpinr) then
        NS = NN*IAbs(NSpin)
        NSpinX = NSpin
      else
        NS = NN
        NSpinX = 1
        endIf
      Do 100 ILower = 1, NS, NumCol
        IUpper = Min(ILower+NumCol-1,NS)
        Num = IUpper - ILower + 1
        Call LabSpn(NSpinX,ILower-1,Num,IRX,LabX)
        If(IOpt.le.2) then
          Write(IOut,1000) (IRX(I),LabX(I),I=1,Num)
        else
          Write(IOut,1040) (IRX(I),LabX(I),I=1,Num)
          endIf
        Do 90 I = 1, MS
          Call LabSpn(NSpin,I-1,1,IRX,LabX)
          Do 10 J = ILower, IUpper
   10       Y(J-ILower+1) = X(IRI,I,J)
          Do 20 J = 1, Num
            If(Key.ne.1.and.Abs(Y(J)).lt.Small) Y(J) = Zero
   20       Continue
          If(IOpt.eq.1) then
            Write(IOut,1010) IRX(1), LabX(1), (Y(J),J=1,Num)
          else if(IOpt.eq.2) then
            Write(IOut,1020) IRX(1), LabX(1), (Y(J),J=1,Num)
          else
            Write(IOut,1030) IRX(1), LabX(1), (Y(J),J=1,Num)
            endIf
   90     Continue
  100   Continue
      Return
      End
