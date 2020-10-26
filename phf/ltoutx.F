*Deck LTOutX
      Subroutine LTOutX(IOut,IType,ISCF,NMat,IOfMat,Str,N,M,XA,XB,Key)
      Implicit Real*8(A-H,O-Z)
C
C     Print matrices which may be rectangular or lower triangular or over DBFs,
C     may be blocked by spin or run over spinor components:
C
C     IType = 0 ... lower triangular Hermetian
C             1 ... square
C             2 ... over dbfs (vectors of length N)
C             3 ... rectangular
C             4 ... lower triangular anti-Hermetian
C     ISCF  = Type of SCF.
C     Str   = String for labelling matrices, blank for none.
C     NMat  = number of matrices
C     N     = Basic dimension, column length per spinor component.
C     M     = Number of rows for rectangular matrices.
C     XA    = Matrices, or alpha matrices if blocked.
C     XB    = Beta spin matrices for RUHF/CUHF.
C
      Logical MSpinr, CBfn, LTMat, Spinor
      Character Str*(*), Label*256
      Dimension XA(*), XB(*)
 1000 Format(1X,A,':')
C
      Call DecoSC(ISCF,NSpBlk,NRI,NDimBl,CBfn,NSpBlX,Spinor)
      NDimBS = NDimBl
      If(CBfn) NDimBS = -NDimBl
      LTMat = IType.eq.0.or.IType.eq.4
      MSpinr = IType.le.1.or.IType.eq.4
      If(MSpinr) then
        ND = N*NDimBl
        MD = M*NDimBl
        MU = M
        NBlock = NSpBlk
      else if(IType.eq.2) then
        NRI = 1
        ND = N
        MD = 1
        MU = 1
        NBlock = 1
      else if(IType.eq.3) then
        ND = N*NDimBl
        MD = M
        MU = M
        NBlock = NSpBlk
      else
        Call GauErr('Unrecognized IType in LTOutX.')
        endIf
      If(IType.eq.0) then
        MNDim = NRI*(ND*(ND+1))/2
      else
        MNDim = NRI*ND*MD
        endIf
      Do 100 IMat = 1, NMat
        IndIM = 1 + (IMat-1)*MNDim
        Do 10 IRI = 1, NRI
          IRIX = IRI
          If(IType.eq.4) IRIX = -IRI
          Call PrtLab(Str,NRI,IRI,NBlock,1,NMat,IMat,IOfMat,Label,
     $      LenLab)
          If(LenLab.gt.0) Write(IOut,1000) Label(1:LenLab)
          If(LTMat) then
            Call LTOutC(IOut,1,NRI,IRIX,NDimBS,N,XA(IndIM),Key)
          else
            Call OutMtC(IOut,1,MSpinr,Key,NRI,IRI,NDimBS,XA(IndIM),ND,
     $        MD,N,MU)
            endIf
   10     Continue
        If(NBlock.eq.2) then
          Do 20 IRI = 1, NRI
            IRIX = IRI
            If(IType.eq.4) IRIX = -IRI
            Call PrtLab(Str,NRI,IRI,NBlock,2,NMat,IMat,IOfMat,Label,
     $        LenLab)
            If(LenLab.gt.0) Write(IOut,1000) Label(1:LenLab)
            If(LTMat) then
              Call LTOutC(IOut,1,NRI,IRIX,NDimBS,N,XB(IndIM),Key)
            else
              Call OutMtC(IOut,1,MSpinr,Key,NRI,IRI,NDimBS,XB(IndIM),ND,
     $          MD,N,MU)
              endIf
   20       Continue
          endIf
  100   Continue
      Return
      End
