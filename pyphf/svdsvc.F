*Deck SVDSVc
      Subroutine SVDSVc(IOut,IPrint,IOpt,NBlks,NBF,NBFU,MaxCmp,NBasis,
     $  NSABF,ThrOvl,S,A,B,E,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Perform SVD on a symmetry-blocked matrix.
C     IOpt = 1 ... Return an orthonormal set of vectors in A.  NBFU is
C                  loaded with the number of vectors in each irrep.
C            2 ... Return the lower triangular generalized inverse in S.
C                  A and NBFU are not used.
C     S ... Input lower-triangular matrix.  Receives inverse for IOpt=2.
C     B ... Storage for blocked matrices over irreps (Sum(I) NBF(I)**2).
C     A ... Scratch storage.  N^2 for IOpt=1, the size of the largest block
C           for IOpt=2.
C     E ... Scratch vector as long as the longest block.
C
C     For IOpt=1, S and A can be the same N^2 matrix.
C         IOpt=2, if NBlks=1 then S can be the same as either A or B but not both.
C                 if NBlks>1 then S, A, and B must be distinct.
C
      Parameter (MinPrt=1,MinMPr=2)
      Dimension NBF(*), NBFU(*), NSABF(NBasis,MaxCmp), S(*), A(*), B(*),
     $  E(*), V(MDV), XX(1)
      Save Zero, One, XX
      Data Zero/0.0d0/, One/1.0d0/, XX/0.0d0/
 1000 Format(' SVDSVc IR=',I1,/,(' SVDSVc V=',10(1PD9.2)))
 1010 Format(' Bad symmetry data in SVDSVc:  NBasis=',I5,' NT=',I5,
     $  ' NBlks=',I1,' MaxCmp=',I1,' NBF=',8I5)
C
C     First form all blocks of the symmetry-adapted matrix.
C
      Call TStamp(1,'Top of SVDSVc.')
      NT = IArSum(1,NBlks,NBF)
      If(NT.ne.NBasis) then
        Write(IOut,1010) NBasis, NT, NBlks, MaxCmp, (NBF(KK),KK=1,NBlks)
        Call GauErr('Bad data in SVDSVc.')
        endIf
      If(NBlks.eq.1) then
        Call Square(S,B,NBasis,NBasis,0)
        RNorm = One
      else
        RNorm = Sqrt(One/GFloat(MaxCmp))
        IndI = 0
        IndB = 1
        Do 10 IR = 1, NBlks
          Call SymAdS(NBasis,NBF(IR),IndI,MaxCmp,NSABF,S,B(IndB))
          IndI = IndI + NBF(IR)
   10     IndB = IndB + NBF(IR)**2
        endIf
      NTT = (NBasis*(NBasis+1))/2
      If(IOpt.eq.2) Call AClear(NTT,S)
      IndI = 0
      IndA = 1
      IndB = 1
      Do 100 IR = 1, NBlks
        If(NBF(IR).gt.0) then
          If(IPrint.ge.MinMPr) Call OutMtS(IOut,'SSO for IR=',IR,0,
     $      B(IndB),NBF(IR),NBF(IR),NBF(IR),NBF(IR))
          Call GGESVD('O','N',NBF(IR),NBF(IR),B(IndB),NBF(IR),E,XX,1,
     $      XX,1,V,MDV,Info)
          If(Info.ne.0) Call GauErr('SVD failed in SVDSVc.')
          If(IPrint.ge.MinPrt) Write(IOut,1000) IR, (E(I),I=1,NBF(IR))
          Do 20 I = NBF(IR), 1, -1
            If((E(I)*RNorm).ge.ThrOvl) goto 30
   20       Continue
          I = 0
   30     NBFU(IR) = I
          Call FixPha(1,NBF(IR),NBFU(IR),0,B(IndB))
          Do 40 I = 1, NBFU(IR)
   40       E(I) = RNorm / Sqrt(E(I))
          Call ScalDg(NBF(IR),NBF(IR),NBFU(IR),B(IndB),B(IndB),E,
     $      .False.,.False.,IErr)
          If(IOpt.eq.1) then
            Call UnSymB(NBasis,NBF(IR),NBFU(IR),IndI,MaxCmp,NSABF,
     $        B(IndB),A(IndA))
          else
            Call XGEMM(1,'N','T',NBF(IR),NBF(IR),NBFU(IR),One,B(IndB),
     $        NBF(IR),B(IndB),NBF(IR),Zero,A,NBF(IR))
            If(NBlks.eq.1) then
              Call Linear(A,S,NBasis,NBasis)
            else
              Call UnSySI(NBasis,NBF(IR),IndI,MaxCmp,NSABF,A,S)
              endIf
            endIf
        else
          NBFU(IR) = 0
          endIf
        IndI = IndI + NBF(IR)
        IndA = IndA + NBasis*NBFU(IR)
  100   IndB = IndB + NBF(IR)**2
      Call TStamp(1,'Bot of SVDSVc.')
      Return
      End
