*Deck GLowdi
      Subroutine GLowdi(IOut,Ptmp,P,XL,SSq,PS,Tmp1,Tmp2,NTT,IrWP,
     $  IOPCl,DLwdX,DLwdY,DLwdZ,Map,ChgLwX,ChgLwY,ChgLwZ,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Computes Lodwin for GHF
CJEP
      Common /PhyCon/ PhyCon(30)
      Parameter(IRwS=514,LMax=4,NCompM=((LMax+1)*(LMax+2)*(LMax+3))/6)
#include "commonmol.inc"
#include "commonb.inc"
      Dimension XL(NTT,4,2), P(NTT,4,2),PTmp(NTT,4,2)
      Integer IndJD(NCompM,3),Map(*),Evt
      Real*8 D(NCompM),ChgLwX(*),ChgLwY(*),
     $ ChgLwZ(*),DLwDX(*),DLwDY(*),
     $ DLwDZ(*),Tmp1(*),Tmp2(*),PS(NTT),SSq(NBasis,*),V(*)
 1005 Format(' Lowdin Atomic Magnetic Moments for ',2A)
 1010 Format(' Sum of Lowdin charges  =',F10.5)
 1015 Format(' Lowdin atomic charges with hydrogens summed into heavy',
     $  ' atoms:')
C
C     Define linear indexing of low-triangular matrix.
C
      LInd(I,J) = (Max(I,J)*(Max(I,J)-1))/2 + Min(I,J)
      ToAng  = PhyCon(1)
      ToE    = PhyCon(3)
      NBSq   = NBasis*NBasis
C
C  Prepare Density matrix to do the job
C
C  First, read and copy the density matrix
C
      Call FileIO(11,IRWP,LenP,RJunk,0)
      Call FileIO(2,-IRwP,LenP,P,0)
      Call AMove(LenP,P,Ptmp)
      Call RIComp(.True.,.False.,.False.,.True.,IOpCl,NBasis,1,1,1,0,
     $  NMat,NMatS,NMatT,Ptmp,Ptmp(1,1,2),V,MDV)
C
C Read Overlap    
C
      Call FileIO(11,IRWS,LenS,RJunk,0)
      Call FileIO(2,-IRwS,LenS,XL,0)
      Call RIComp(.True.,.False.,.False.,.True.,IOpCl,NBasis,1,1,1,0,
     $  NMat,NMatS,NMatT,XL,XL(1,1,2),V,MDV)
C
C     Convert to square format
C
      Call Square(XL(1,1,1),SSq,NBasis,NBasis,0)
C
C     Get S^1/2
C
      Call RootMt(SSq,PS,Tmp1,Tmp2,NBasis,NBasis,0)
C
C     Get the (S^1/2)(P)(S^1/2) matrix
C
C  X component
C
      Call AMove(NTT,Ptmp(1,2,1),PS)
      Call AScale(NTT,2.0D0,PS,PS)
c
      Call VFVMul(.False.,.False.,.False.,0,NBasis,NBasis,NBasis,
     $    PS,SSq,Tmp1,DLwdX,Tmp2,NBsq)
C
C  Y component
C
      Call AMove(NTT,Ptmp(1,2,2),PS)
      Call AScale(NTT,2.0D0,PS,PS)
C
      Call VFVMul(.False.,.False.,.False.,0,NBasis,NBasis,NBasis,
     $    PS,SSq,Tmp1,DLwdY,Tmp2,NBsq)
C
C  Z component
C
C      ZNa = SCFTrc(XL(1,1,1),Ptmp(1,1,1),NBasis,1)
C      ZNb = SCFTrc(XL(1,1,1),Ptmp(1,1,2),NBasis,1)
C      Write(IOut,*)'** >>>>> Na = ',ZNa
C      Write(IOut,*)'** >>>>> Nb = ',ZNb
C
C 
      Call ASub(NTT,Ptmp(1,1,1),Ptmp(1,1,2),PS)
      Call VFVMul(.False.,.False.,.False.,0,NBasis,NBasis,NBasis,
     $    PS,SSq,Tmp1,DLwdZ,Tmp2,NBsq)

C
C     Set map between AOs and atoms and compute Lowdin charges.
C
      Call RdMap(Map)
      Call AClear(NAtoms,ChgLwX)
      Call AClear(NAtoms,ChgLwY)
      Call AClear(NAtoms,ChgLwZ)
C 
C
      Do 10 I = 1, NBasis
        ChgLwX(Map(I)) = ChgLwX(Map(I)) + DLwdX(LInd(I,I))
        ChgLwY(Map(I)) = ChgLwY(Map(I)) + DLwdY(LInd(I,I))
   10   ChgLwZ(Map(I)) = ChgLwZ(Map(I)) + DLwdZ(LInd(I,I))
C
C     Print Lowdin charges.
C
      Write(IOut,1005),'Sx'
C      Call AtomPr(IOut,.true.,NAtoms,IAn,ChgLwX,NAtoms,1)
      Call AtomPr(IOut,NAtoms,IAn,ChgLwX,NAtoms,1)
      Write(IOut,1005),'Sy'
C      Call AtomPr(IOut,.true.,NAtoms,IAn,ChgLwY,NAtoms,1)
      Call AtomPr(IOut,NAtoms,IAn,ChgLwY,NAtoms,1)
      Write(IOut,1005),'Sz'
C      Call AtomPr(IOut,.true.,NAtoms,IAn,ChgLwZ,NAtoms,1)
      Call AtomPr(IOut,NAtoms,IAn,ChgLwZ,NAtoms,1)
C
C
      Return
      End
