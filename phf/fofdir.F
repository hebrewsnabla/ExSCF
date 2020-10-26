*Deck FoFDir
      Subroutine FoFDir(IOut,IPrint,IHMeth,IOpCl,ICntrl,IRafI,IPFlag,
     $  AllowP,ICnBeg,ICnEnd,AddH,InitF,DoPurF,AccDes,ScaHFX,NFrqR1,
     $  NMat,NMatS,NMatT,NMatD,NMatDS,NMatDT,NBasis,ISym2E,NSymOp,NOpAb,
     $  NEqAtm,NEqShl,RotOp,NEqBas,HA,HB,PZA,PZB,PA,PB,CIDenA,CIDenB,FA,
     $  FB,FDA,FDB,NAtoms,IAn,C,UseAtT,IAtTyp,Frozen,FXYZ,FFXYZ,LenSav,
     $  NSaved,R03,R1,R2,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Copyright (c) 1988,1990,1992,1993,1995,1998,2003,2009,2011,
C                Gaussian, Inc.  All Rights Reserved.
C
C     This is part of the Gaussian(R) 09 program.  It is based on
C     the Gaussian(R) 03 system (copyright 2003, Gaussian, Inc.),
C     the Gaussian(R) 98 system (copyright 1998, Gaussian, Inc.),
C     the Gaussian(R) 94 system (copyright 1995, Gaussian, Inc.),
C     the Gaussian 92(TM) system (copyright 1992, Gaussian, Inc.),
C     the Gaussian 90(TM) system (copyright 1990, Gaussian, Inc.),
C     the Gaussian 88(TM) system (copyright 1988, Gaussian, Inc.),
C     the Gaussian 86(TM) system (copyright 1986, Carnegie Mellon
C     University), and the Gaussian 82(TM) system (copyright 1983,
C     Carnegie Mellon University). Gaussian is a federally registered
C     trademark of Gaussian, Inc.
C
C     This software contains proprietary and confidential information,
C     including trade secrets, belonging to Gaussian, Inc.
C
C     This software is provided under written license and may be
C     used, copied, transmitted, or stored only in accord with that
C     written license.
C
C     The following legend is applicable only to US Government contracts
C     under DFARS:
C
C                        RESTRICTED RIGHTS LEGEND
C
C     Use, duplication or disclosure by the US Government is subject to
C     restrictions as set forth in subparagraph (c)(1)(ii) of the Rights
C     in Technical Data and Computer Software clause at DFARS
C     252.227-7013.
C
C     Gaussian, Inc., 340 Quinnipiac St., Bldg. 40, Wallingford CT 06492
C
C     The following legend is applicable only to US Government contracts
C     under FAR:
C
C                        RESTRICTED RIGHTS LEGEND
C
C     Use, reproduction and disclosure by the US Government is subject
C     to restrictions as set forth in subparagraph (c) of the Commercial
C     Computer Software - Restricted Rights clause at FAR 52.227-19.
C
C     Gaussian, Inc., 340 Quinnipiac St., Bldg. 40, Wallingford CT 06492
C
      Parameter (MinMtS=6,IRwSEP=768)
C
C     Form Fock matrices by computing and summing the necessary
C     integrals.  Arguments:
C
C     IOut   ... Output unit.
C     IPrint ... Debug print.  -1 means no print, 0 means no print
C                except when called for the first time, 1 means a
C                little print, 2 means more ... 5 prints integrals.
C     IHMeth ... Method:
C                0 ... Default (1).
C                1 ... Ab initio
C                2 ... CNDO/2
C                3 ... INDO/2
C                4 ... ZINDO/1 (NYI)
C                5 ... ZINDO/S
C                6 ... MINDO/3 (NYI)
C                7 ... MNDO
C                8 ... AM1
C                9 ... PM3
C               10 ... PM3MM
C            11,12 ... Harris functional.  All the terms for this
C                      are done in CalDSu, so this routine just
C                      initializes empty matrices.
C            13,14 ... Approximate DFT, but with normal Coulomb.
C                      Treated as ab initio here.
C               15 ... Regular SCF with separate K, for testing,
C                      Treated as ab initio here.
C               16 ... Regular SCF with separate NDDO K,
C                      Treated as ab initio here.
C               17 ... K part of model 16, used internally.
C               18 ... DFT-SCTB using interpolated values.
C               19 ... DFT-SCTB using analytic values.
C               20 ... EHT-SC.
C               21 ... 2e terms forced to be zero.
C            22-39 ... Reserved for non-zdo semi-empirical models.
C               40 ... PM6.
C               41 ... PDDG/PM3.
C            42-99 ... Reserved for zdo semi-empirical models.
C              101 ... Molecular Mechanics (not used here, but
C                      relevent for PCM in FoFCou).
C     IOpCl  ... 0/>0 Closed/Open shell, or type of integrals to
C                return.
C     ICntrl ... Algorithm control:
C                 0 ... Default (1227).
C                 1 ... Obsolete
C                 2 ... Obsolete
C                 3 ... Obsolete
C                 4 ... Prism up to L=1.
C                 5 ... Obsolete.
C                 6 ... Prism df (no sp done here at all).
C                 7 ... Default.  Do everything here.
C                 8 ... Prism up to L=2.
C                10 ... No cutoffs.
C                20 ... Cutoffs for 10**-10 accuracy.
C                30 ... Cutoffs for high accuracy.
C                40 ... Sleazy (10**-6) Cutoffs.
C               100 ... Do not compute operator matrices.
C               200 ... Compute SCF Fock matrices:
C                       (DA,DB)*(R1,R2)==>FA,FB (densities 1...NMatS)
C                       (DA,DB)*R2     ==>FA,FB
C                         (densities NMatS+1...NMatS+NMatT)
C                       (DA,DB)*R3     ==>FA,FB
C                         (densities NMatS+NMatT+1...NMat)
C               300 ... Compute CIS operators.  These are HF-style
C                       Fock matrices but done with square, non-symmetric
C                       densities passed in and square Fock matrices
C                       returned.
C               400 ... Compute WTilda terms.  These are exchange
C                       matrices computed with NMatS symmetric and
C                       NMat-NMatS antisymmetric lower triangular
C                       densities, returned lower triangular.  NMatT
C                       should be zero.
C               500 ... DFT: DA==>FJ (in FA) (IOpCl=0)
C                       GVB: DA==>FJ,FK (in FA,FB) (IOpCl=1)
C                       K:   DA==>FK (in FA) (IOpCl=2)
C                       K:   DA,DB==>FKA,FKB (in FA,FB) (IOpCl=3)
C               600 ... Compute regular integrals and load into R03,
C                       in canonical form if iopcl=0, in square form
C                       if iopcl=1, and 1c only lower triangualar if
C                       iopcl=2.  NMatS is used as the dimension of R0
C                       if iopcl=1, and should be the number of 1c basis
C                       pairs if iopcl=2.  For IOpCl=2, IRafI-1 indicates
C                       the maximum L value for retained shell pairs, or
C                       IRafI=0 requests all 1c pairs.  For IOpCl=3, the
C                       Coulomb integrals over pairs of DBFs given in /B2/
C                       are returned instead.
C               700 ... Compute raffenetti integrals:
C                       IOpCl=0 ... Load R1.
C                       IOpCl=1 ... Load R1 and R2.
C                       IOpCl=2 ... Load R1 and R03.
C                       IOpCl=3 ... Load R1, R2, and R03.
C                       IOpCl=4 ... Load R2 and R03.
C              1000 ... Do not compute forces.
C              2000 ... Compute forces.
C              3000 ... Make derivative Fock/Coulomb matrices.  Third order
C                       property derivatives PiGx(Pj) will also be done if
C                       NFrqR1>0.  The matrices are FA/FB for ICntr2<=0,
C                       or J/K if ICntr2=5.
C              4000 ... Make Gm(Pe) matrices over GIAOs.
C              5000 ... Compute forces using including CIS 2PDM terms.
C              6000 ... Compute dF/dB over GIAO's. X,Y and Z matrices
C                       are returned in FDA (and FDB, if IOpCl.eq.1)
C              7000 ... Make V and pVp integrals for DKH.  NYI here
C              8000 ... Do Forces and Fock matrices.
C              9000 ... Make V, pVp, and 1e SO ints for DKH.  NYI here.
C            100000 ... Compute second derivatives.
C                       Note  that SCF Fock matrices, Fock derivatives,
C                       and forces can be combined, but not J and K
C                       operators and derivatives or second derivatives.
C                       Also, square input and output matrices should be
C                       in arrays allocated 2*NTT6D*NMat rather than just
C                       NBas6D**2*NMat.
C            200000 ... Compute magnetic field second derivatives over
C                       GIAOs.  The 9 components are returned in FFXYZ.
C            300000 ... Make Gmx(P) matrices. (set ICntrl=306127)
C           N000000 ... IR12Ty=N.  Only N=1 and 2 are implemented.
C     IRafI  ... Whether to form Raffenetti integral combinations.
C                Only works if 1st derivatives are not computed;
C                Use of Raffenetti is only sensible if NMat is large
C                 0 ... Default (let the program decide the partition
C                       between regular and Raffenetti ERI's). The
C                       Raffenetti combinations are only made if NMat
C                       is large.
C                 n ... Do ERI's with degree of contraction greater
C                       than or equal to n as regular ERI's, those
C                       less than as Raffenetti.  A negative n means
C                       if NMat is at least -n, then all ERI's are made
C                       as Raffenetti, otherwise all are regular.
C                       So IRafI=1 means all integrals are done as
C                       regular, while IRafI=-1 means all are done as
C                       Raffenetti.
C     IPFlag ... Flags for Prism, PrfRaf, CalDSu, etc.:
C                1 ... Use AllowP(50) to turn off paths.
C                2 ... If set, force front and back end memory allocations
C                      to match.  Otherwise, allocate independently.
C                4 ... Run a script to execute Prism, PrmRaf, or CalDFT.
C                8 ... Use unformatted I/O on the data for the script process.
C               16 ... Do extra work to use cutoffs better in DFTMat.
C               32 ... Reverse normal choice of diagonal/canonical sampling
C                      in Prism and PrmRaf.  The default is diagonal only
C                      on vector machines.
C               64 ... Trace input and output using Linda/subprocess.
C              128 ... Force single matrix code in CPKS.
C              256 ... Force all near field in FMM
C              512 ... Turn off vFMM.
C             1024 ... Force square loops (currently only in PrismC).
C             2048 ... Force use of FoFCou, even if not doing FMM.
C     AllowP ... Logical(50) array for PRISM if IPFlag is set.
C     ICnBeg ... First center to include in Fock derivative formation.
C     ICnEnd ... Last center to include in Fock derivative formation.
C                The Fock derivative arrays should be dimensioned
C                3*(NAtoms+1)*NTT, or if not all atoms are being handled
C                then (NMtFkD+3)*NTT, where NMtFkD=(ICnEnd-ICnBeg+1),
C                unless ICnBeg=ICnEnd=0 in which case only 3*NTT are
C                required.
C     AddH   ... Whether to add HA/HB into FA/FB.
C     InitF  ... Whether to initialize Fock matrices or Fock derivatives
C                to 0.  If .False., F and/or FD are assumed to contain
C                contributions to the the matrices over cartesians.
C     DoPurF ... Whether to convert the Fock matrices or Fock derivatives
C                back to pure functions or to leave them in terms of
C                cartesians.
C     AccDes ... Desired accuracy.
C     NFrqR1 ... Number of different frequencies for the NMatD-NMat matrices
C                used for third-order properties.
C     NMat   ... Number of Fock matrices to form.
C     NMatS  ... Number of the NMat density matrices which are symmetric
C                singlets.
C     NMatT  ... Number of the NMat density matrices which are symmetric
C                triplets.  The next NMat-NMatS-NMatT densities are
C                assumed to be antisymmetric.
C     NMatD  ... The total number of density matrices passed in PA and PB.
C                This may be greater than NMat if polarizability derivatives
C                are being done along with F(x), in which case the NMatD-NMat
C                densities at the end of PA and PB are assumed to be the field
C                derivative densities.
C     NMatDS ... Number of the NMatD-NMat matrices for polar derivs which are
C                symmetric singlets.
C     NMatDT ... Number of the NMatD-NMat matrices for polar derivs which are
C                symmetric triplets.  The remaining NMatD-NMat-NMatDS-NMatDT
C                are assumed to be antisymmetric.
C     NBasis ... Number of basis functions.
C     ISym2E ... 1 for use of symmetry in integral evaluation by
C                symmetrization of Fock matrices, 2 for use of symmetry
C                by replicating integrals.  -1 to compute only unique
C                integrals but skip symmetrization here.
C     NSymOp ... Order of the point group to be used for symmetrization
C                of results (JSym2E=1).
C     NOpAb  ... Order of the abelian subgroup to use if integrals are
C                being replicated (JSym2E=2) or if the Abelian only
C                flag is set.
C     NEqAtm ... Used for derivatives only when symmetry is turned on.
C     NEqShl ... (nshell*nsymop), used if symmetry is turned on.  Should be
C                the array for /B2/ if doing (DBF|DBF) integrals.
C     RotOp  ... Rotation matrices for ISym2E=1.
C     NEqBas ... used for replication if ISym2E=2.
C     PZA    ... Alpha generalized density, used only if ICntr3 is 2 or
C                5. Should be allocated NTT6D.  Should be 1/2 the SCF
C                density for the SCF gradient case.
C     PZB    ... Beta generalized density, used only if ICntrl is 2 or
C                5. Should be allocated NTT6D, and must be separate from PZA
C                for the CIS case even for closed-shell. Should be 1/2 the
C                SCF density for the SCF gradient case.
C     PA     ... Input vector of alpha density matrices (NTT,NMatD),
C                where NMatD is NMat normally, but NMat+3 if the
C                polarizability derivative "forces" are computed. This
C                is the alpha density or 1/2 (PA+PB) for Coulomb.  May be
C                destroyed here.  Should be allocated NTT6D.
C     PB     ... Input vector of beta density matrices (NTT,NMatD) (only
C                used if IOpCl.eq.1).  May be destroyed here.
C                Should be allocated NTT6D.
C     CIDenA ... Alpha CIS 1PDM for CIS gradients, used only if ICntr3
C                is 5. Should be allocated NBSq6D.
C     CIDenB ... Beta CIS 1PDM for CIS gradients, used only if ICntr3
C                is 5. Should be allocated NBSq6D. Must be separate from
C                CIDenA for the CIS case even for closed-shell.
C     HA     ... Alpha Core Hamiltonian, only used if AddH is .True.
C     HB     ... Beta Core Hamiltonian, only used if AddH is .True.
C     FA     ... Output vector of alpha Fock matrices (NTT,NMat).
C                Should be allocated NTT6D.
C     FB     ... Output vector of beta Fock matrices (NTT,NMat) (only
C                used if IOpCl.ne.0). Should be allocated NTT6D.
C     FDA    ... Output alpha derivative Fock matrices... only used if
C                ICntrl is appropriately set.  Should be allocated NTT6D.
C     FDB    ... Output beta derivative Fock matrices  Should be
C                allocated NTT6D.
C     NAtoms ... Number of atoms -- can be 0 if forces are not requested
C                and the calculation is ab initio.
C     IAn    ... Atomic numbers, only needed for semi-empirical methods
C     C      ... Atomic coordinates.  Only used if forces are requested
C                and for semi-empirical methods.
C     IAtTyp ... Flags for active atoms.  Only used for nuclear F(x) and
C                for semi-empirical Hamiltonians.
C     Frozen ... Logical array indicating frozen atoms, only used if
C                forces are requested.
C     FXYZ   ... Returned forces, if requested.  This is incremented.
C     FFXYZ  ... Returned second derivatives, if requested.  This is
C                incremented.
C     LenSav ... >0 to stash the first LenSav computed integrals in R03;
C                <0 to re-use up to |LenSav| previously computed
C                integrals in R03.  Turns off density-based cutoffs.
C     NSaved ... Returned number of saved integrals if LenSav>0.
C     R03    ... For in-core calculations, regular or Raffenetti 3
C                integrals are returned here.  If LenSav > 0 then
C                computed integrals in the order formed are saved here.
C     R1     ... Recevies Raff 1 integrals for in-core.
C     R2     ... Recevies Raff 2 integrals for in-core.
C     MDV    ... Length of V.
C     V      ... Used for dynamic allocation.
C
C     Fock and density matrices will be lower triangular or square
C     depending on the type of job. Note that PA,PB,FA,FB must be
C     allocated as (NTT6D or NBSq6D,NMat) even if pure functions are in
C     use.  If Raffenetti ERI's are used with CIS then the dimension
C     must be NMat*NBas6D*(NBas6D+1).
C
#include "commonb2.inc"
      Real*8 MDCutO
      Common/SymInf/IDum1,IDum2,JTrans(3,8),T(3,3),TrVec(3)
      Dimension PA(*),PB(*),FA(*),FB(*),FDA(*),FDB(*),PZA(*),PZB(*),
     $  CIDenA(*),CIDenB(*),NEqShl(*),NEqBas(NBASIS,*),HA(*),HB(*),
     $  FFXYZ(*),C(3,*),FXYZ(3,*),NEqAtm(*),R03(*),R1(*),R2(*),RotOp(*),
     $  IAtTyp(*),IAn(*),V(*),JJA(1),XX(1)
      Logical Called, AddH, DoFock, DoForc, SqFock, Frozen(1), SEInt,
     $  DoRegI, DoRafI, DoFokD, InGaus, DoInt, JustIn, DoCISG, UseZ,
     $  DoFreq, HaveDB, HaveFB, DoP, InitF, DoPurF, AllowP(*), AM1Ham,
     $  LJJ(1), AbSym, NAbSym, AbOnly, UseAtT, DoPuFL, DoAdj, SEHam,
     $  HarHam, NoGath, LBit, DBFInt, DoChi, DoPolD, BraDBF, KetDBF,
     $  BrDBCon, BraDBO, DoGmx, DoGmPe, PeGmx, UseFSE
      Save Called, TenM6, TenM10, JJA, XX, LJJ
      Data TenM6/1.0d-6/, Called/.False./, TenM10/1.0d-10/, JJA/0/,
     $  XX/0.0d0/, LJJ/.False./
 1000 Format(' Enter FofDir, NBasis=',I4,' IHMeth=',I2,'.')
 1010 Format(' FoFDir skips integrals for Harris functional.')
 1060 Format(' Symmetry not used in FoFDir.')
 1070 Format(' Petite list used in FoFDir.')
 1080 Format(' Integrals replicated using symmetry in FoFDir.')
 1100 Format(' Unrecognized ICntrl=',I10,' in FoFDir.')
 1110 Format(' MinBra=',I2,' MaxBra=',I2,' Meth=',I2,'.')
 1140 Format(' Unrecognized ICntr2=',I10,' in FoFDir.')
 1150 Format(' Unrecognized ICntr3=',I10,' in FoFDir.')
 1160 Format(' Unrecognized ISym2E=',I10,' in FoFDir.')
 1180 Format(' IRaf=',I8,' NMat=',I4,' IRICut=',I8,' DoRegI=',L1,
     $       ' DoRafI=',L1,' ISym2E=',I2,' JSym2E=',I1,'.')
 1220 Format(' ICntr1=',I12,' error.')
 1230 Format(' Raf error: ICntr2=',I12,' IOpCl=',I1,'.')
 1240 Format(' Raf error: DoRegI=',L1,'.')
C
      IPrt1 = IPrint
      If(IPrt1.eq.0.and..not.Called) IPrt1 = 1
      If(IPrt1.lt.0) IPrt1 = 0
      If(IPrt1.ge.2) Write(IOut,1000) NBasis, IHMeth
      NoGath = LBit(IPFlag,15)
      If(IHMeth.le.1.or.IHMeth.eq.13) then
        Meth = 1
      else
        Meth = IHMeth
        endIf
      SEInt = SEHam(Meth)
      UseFSE = SEInt.and..not.(AM1Ham(IHMeth).and.ITqry(IRwSEP).gt.0)
      If(SEInt) then
        JSym2E = 0
        IPurAO = 0
      else
        JSym2E = IAbs(ISym2E)
        IPurAO = 1
        endIf
      If(JSym2E.eq.0) then
        NOp = 1
        If(IPrt1.ge.1) Write(IOut,1060)
      else if(JSym2E.eq.1) then
        NOp = NSymOp
        If(IPrt1.ge.1) Write(IOut,1070)
      else if(ISym2E.eq.2) then
        NOp = NOpAb
        If(IPrt1.ge.1) Write(IOut,1080)
      else
        Write(IOut,1160) ISym2E
        Call Lnk1E(0)
        endIf
      If(AbOnly(0)) NOp = Min(NOp,NOpAb)
      NOpAbU = Min(NOp,NOpAb)
      If(ISym2E.eq.2) then
        If(InGaus(0)) Call GetNOp(I,J)
      else
        JTrans(1,1) = 1
        JTrans(2,1) = 1
        JTrans(3,1) = 1
        endIf
      MxDCon = MxCont(0)**4
C
C     Determine partition between PRISM (BraKet), HGP, and Rys methods.
C
      Call DecICn(IOpCl,Meth,ICntrl,NFrqR1,ICntr0,ICntr1,ICntr2,ICntr3,
     $  ICntr4,ICnBeg,ICnEnd,NAtoms,3,NMat,NMatS,NMatT,NMatD,NMatDS,
     $  NMatDT,1,0,DoFock,DoForc,DoFokD,DoFreq,DoPolD,DoChi,HaveDB,
     $  HaveFB,NFSets,IFkDrv,NMtFkD,IR12Op,NMatF,NMatSF,NMatTF,NMtPBF,
     $  BraDBF,KetDBF,BrDBCon,BraDBO,DoGmx,DoGmPe,PeGmx)
      IDens = Max(ICntr2-1,0)
      DBFInt = IDens.eq.5.and.IOpCl.eq.3
      If(DBFInt) then
        NDimEq = NShelB
      else
        NDimEq = NShell
        endIf
      DoAdj = .not.SEInt
      ISEAll = 1
      IV = ISEAll
      Call TstCor(ISEAll,MDV,'FoFDir-SEAll')
      Call IClear(IntPWP(0),V(ISEAll))
      If(SEInt) then
        IR12Op = 3
        NBas6D = NBasis
        If(.not.UseFSE) Call GetSEP(IOut,IPrint,IHMeth,0,.False.,ISEAll,
     $    MaxASE,NSEAtI,NSEAtR,NCorPA,MaxCPA,NSEElI,NSEElR,jRSESc,jISEP,
     $    jRSEP,jSECP,jIElP,jRElP,IV,V,MDV)
      else if(DBFInt) then
        Call GetMBN(1,1,NShelB,ShelTB,ShelCB,NBas6D)
      else
        Call GetNB6(NBas6D)
        endIf
      LSEAll = IV - ISEAll
      If(LSEAll.eq.0) IV = ISEAll + 1
      If(UseFSE.or.(HarHam(Meth).and..not.DBFInt)) then
        MinBra = -1
        MaxBra = -1
      else if(DBFInt) then
        MinBra = 0
        MaxBra = MaxTyB
      else if(ICntr0.eq.4) then
        MinBra = 0
        MaxBra = 1
      else if(ICntr0.eq.6) then
        MinBra = 2
        MaxBra = 3
      else if(ICntr0.eq.7) then
        MinBra = 0
        MaxBra = MaxTyp
      else if(ICntr0.eq.8) then
        MinBra = 0
        MaxBra = 2
      else
        Write(IOut,1100) ICntrl
        Call Lnk1E(0)
        endIf
      If(ICntr2.lt.1.or.ICntr2.gt.7) then
        Write(IOut,1140) ICntr2
        Call Lnk1E(0)
        endIf
      If(ICntr3.lt.1.or.ICntr3.eq.7.or.ICntr3.gt.8) then
        Write(IOut,1150) ICntr3
        Call Lnk1E(0)
        endIf
C
C     Set flags for Fock matrix & force construction, & whether the
C     Fock & density matrices are square or lower triangular
C
      JustIn = ICntr2.eq.6.or.ICntr2.eq.7
      SqFock = ICntr2.eq.3
      DoInt = DoFock.or.JustIn
      DoCISG = ICntr3.eq.5
      UseZ = ICntr3.eq.2.or.ICntr3.eq.5.or.(DoFreq.and..not.DoFokD)
      If(ICntr2.eq.6.and.IOpCl.eq.2) then
        LMaxSP = IRafI
        IRaf = 0
      else
        LMaxSP = -1
        IRaf = IRafI
        endIf
      If(IPrt1.ge.1) Write(IOut,1110) MinBra, MaxBra, Meth
C
C     Decide what to do about Raffenetti integrals (decode IRaf)
C     In the event that the input is zero, we only do Raffenetti if
C     there are more than 15 matrices, and then we make the partition
C     between regular and Raffenetti at a degree of contraction of
C     NMat.  Symmetry in remaining regular integrals is done via
C     replication.  For IRaf=-N, if there are at least N matrices,
C     IRICut is set to a value greater than the maximum degree of
C     contraction, forcing all Raffenetti.  If -NMat<IRaf<0, IRICut
C     is set to 1, forcing all regular integrals.
C
      IRICut = ICutRI(IOut,.False.,SEInt,IOpCl,ICntr2,ICntr3,ICntr4,
     $  IR12Op,IRaf,NMat,NMatS,NMatT,MxDCon)
      DoRegI = IRICut.le.MxDCon
      DoRafI = IRICut.gt.1
      If(IPrt1.gt.0) Write(IOut,1180) IRaf,NMat,IRICut,DoRegI,DoRafI,
     $  ISym2E,JSym2E
      AbSym = (JSym2E.eq.1.and.NOp.eq.NOpAbU).or.JSym2E.eq.2
      NAbSym = NOp.gt.NOpAbU.and.JSym2E.eq.1
      If(AbSym.and.NAbSym) Call GauErr('Logic error #1 in FoFDir.')
      NFSDim = Max(1,NFSets)
C
C     Allocate memory
C
      NAt3 = 3*NAtoms
      If(IFkDrv.eq.2.or.DoChi) then
        NVDiff = 3
      else
        NVDiff = NAt3 + 3
        endIf
      IFTemp = IV
      If(PeGmx) then
        LenFT = NAt3*3*NFSets
      else if(DoForc) then
        LenFT = NVDiff*NFSets
      else
        LenFT = 0
        endIf
      If(DoChi) then
        LenFFT = NVDiff**2
      else if(DoFreq) then
        LenFFT = (NVDiff*(NVDiff+1))/2
      else
        LenFFT = 0
        endIf
      IFFTmp = IFTemp + LenFT
      JBeg = IFFTmp + LenFFT
      DoP = .not.JustIn
      If(ICntr2.eq.4) then
        NMatST = NMatS
      else
        NMatST = NMatS + NMatT
        endIf
      If(DoGmPe) then
        NMatSO = NMat
      else if(ICntr3.eq.6) then
        NMatSO = 0
      else
        NMatSO = NMatST
        endIf
      NMatE = NMatD - NMat
      NTT = (NBasis*(NBasis+1))/2
      NTT6D = (NBas6D*(NBas6D+1))/2
      IndE = NMat*NTT6D
      If(DBFInt) then
        IPureD = I5DB2
        IPureF = I7FB2
        Call GenAOS(1,1,NShelB,ShelTB,ShelCB,Junk,AOSB)
      else if(.not.SEInt) then
        If(DoP.and.NMatE.gt.0) then
          Call AMove1(NMatE*NTT,NMat*NTT,IndE,PA)
          If(HaveDB) Call AMove1(NMatE*NTT,NMat*NTT,IndE,PB)
          Call P2Cart(IOut,IPrint,.False.,DoP,.False.,.False.,.False.,
     $      .False.,SqFock,HaveDB,.True.,NBasis,NBas6D,NOpAbU,NMatE,
     $      NMatDS+NMatDT,0,IPureD,IPureF,PA(1+IndE),PB(1+IndE),XX,XX,
     $      XX,XX,NEqBas,JJA,XX,LJJ,LJJ,V(JBeg),MDV-JBeg+1)
          endIf
        Call P2Cart(IOut,IPrint,.False.,DoP,UseZ,DoCISG,.False.,.False.,
     $    SqFock,HaveDB,.True.,NBasis,NBas6D,NOpAbU,NMat,NMatST,0,
     $    IPureD,IPureF,PA,PB,PZA,PZB,CIDenA,CIDenB,NEqBas,JJA,XX,LJJ,
     $    LJJ,V(JBeg),MDV-JBeg+1)
        endIf
      If(SqFock) then
        NDim6D = NBas6D*NBas6D
      else
        NDim6D = NTT6D
        endIf
      If(DoFock.and.(NMatS+NMatT).gt.NMat)
     $  Call GauErr('Illegal combination of NMat, NMatS, and NMatT.')
      Call AClear(LenFT,V(IFTemp))
      If(InitF) then
        If(DoFock) then
          LenClr = NDim6D*NMat
          Call AClear(LenClr,FA)
          If(HaveFB) Call AClear(LenClr,FB)
          endIf
        If(DoFokD) then
          LenClr = NDim6D*NMtFkD*NMat
          Call AClear(LenClr,FDA)
          If(HaveFB) Call AClear(LenClr,FDB)
          endIf
      else if(DoAdj) then
        If(DoFock) then
          Call AdDF1N(SqFock,.False.,NDim6D,NMat,FA)
          If(HaveFB) Call AdDF1N(SqFock,.False.,NDim6D,NMat,FB)
          endIf
        If(DoFokD) then
          Call AdDF1N(SqFock,.False.,NDim6D,NMat*NMtFkD,FDA)
          If(HaveFB) Call AdDF1N(SqFock,.False.,NDim6D,NMat*NMtFkD,FDB)
          endIf
        endIf
C
C     Do integrals.
C
      If(ICntr1.eq.1) then
        IAcrcy = 201
      else if(ICntr1.eq.2) then
        If(AccDes.lt.TenM6) then
          IAcrcy = 10
        else
          IAcrcy = 12
          endIf
      else if(ICntr1.eq.3) then
        IAcrcy = 11
      else if(ICntr1.eq.4) then
        IAcrcy = 12
      else
        Write(IOut,1220) ICntr1
        Call Lnk1E(0)
        endIf
      If(ICntr1.eq.1) then
        CutWrt = MDCutO(0)
      else if(ICntr1.eq.2) then
        CutWrt = TenM10
      else
        CutWrt = TenM6
        endIf
      CutWrt = Max(AccDes,CutWrt)
      If(IFkDrv.ne.0.and.Mod(IAcrcy,100).ge.10) IAcrcy = IAcrcy - 10
      If(NAbSym.and.SqFock) Call GauErr('Illegal non-Abelian symmetry.')
      If(MaxBra.ge.0) then
        JLeft = MDV - JBeg + 1
        ILeft = JLeft*IntPWP(0)
        IWhat = 0
        If(DoInt) IWhat = IWhat + 1
        If(DoForc.or.DoFokD) IWhat = IWhat + 10
        If(DoFreq) IWhat = IWhat + 100
        If(DoChi) IWhat = IWhat + 200
        If(DoGmx.or.PeGmx) IWhat = 200
C       If NMat is large enough, we disable gather/scattering of Fock and
C       density matrices.  Note that CIS gather/scatter has not yet been debugged.
        MinSpl = Max(MinMtS,MDVLen(-1))
        If(NoGath) MinSpl = 1
        If(IDens.eq.3.or.JustIn) then
          IGathr = 1
        else if(IWhat.ge.10) then
          IGathr = 2
        else if(NMat.ge.MinSpl) then
          IGathr = 1
        else
          IGathr = 2
          endIf
        If(LenSav.eq.0) then
          If(ICntr2.eq.6.and.IOpCl.eq.0) then
            LenClr = (NTT6D*(NTT6D+1))/2
            Call AClear(LenClr,R03)
          else if(ICntr2.eq.6.and.IOpCl.eq.1) then
            Call AClear(NMatS*NTT6D,R03)
          else if(DBFInt) then
            Call AClear(NTT6D,R03)
          else if(ICntr2.eq.6.and.IOpCl.ge.2) then
            LenClr = (NMatS*(NMatS+1))/2
            Call AClear(LenClr,R03)
          else if(ICntr2.eq.7) then
            LenClr = (NTT6D*(NTT6D+1))/2
            If(IOpCl.le.3) Call AClear(LenClr,R1)
            If(IOpCl.eq.1.or.IOpCl.eq.3.or.IOpCl.eq.4)
     $        Call AClear(LenClr,R2)
            If(IOpCl.eq.2.or.IOpCl.eq.3.or.IOpCl.eq.4)
     $        Call AClear(LenClr,R03)
            endIf
          endIf
C       Do regular integrals
        Call AClear(LenFFT,V(IFFTmp))
        If(DoRegI) then
          If(IHMeth.ne.21) Call FoFBra(IOut,IPrint,IHMeth,IR12Op,IWhat,
     $      InitF,IDens,IGathr,IAcrcy,IPrt1,MinBra,MaxBra,IRaf,IRICut,
     $      DoRafI,IFkDrv,ICnBeg,ICnEnd,IOpCl,NBas6D,NFrqR1,NMat,NMatS,
     $      NMatT,NMatD,NMatDS,NMatDT,NMtFkD,NFSDim,NDim6D,NAtoms,
     $      JSym2E,NOp,NOpAbU,NDimEq,NEqShl,JTrans,PZA,PZB,DoCISG,
     $      CIDenA,CIDenB,PA,PB,FA,FB,FDA,FDB,IAn,C,IAtTyp,Frozen,
     $      V(IFTemp),V(IFFTmp),R03,R1,R2,AccDes,ScaHFX,IPFlag,AllowP,
     $      LSEAll,V(ISEAll),DoFock,DoForc,DoFokD,DoPolD,DoFreq,DoGmx,
     $      DoGmPe,PeGmx,NFSets,LenSav,NSaved,LMaxSP,V(JBeg),JLeft,
     $      V(JBeg),ILeft)
          If(LenSav.eq.0) then
            If(ICntr2.eq.4) then
              Call FixWT(NBas6D,NTT6D,NMat,NMatS,FA)
            else if(ICntr2.eq.6.and.IOpCl.eq.0) then
              Call FixRaf(NShell,ShellT,ShellC,3,NBas6D,R03)
            else if(ICntr2.eq.6.and.IOpCl.eq.1) then
              Call FixReg(NShell,ShellT,ShellC,NBas6D,NMatS,R03)
            else if(ICntr2.eq.6.and.IOpCl.eq.2) then
              IP1C2I = jBeg
              IP1C2J = IP1C2I + NMatS
              IMu2Sh = IP1C2J + NMatS
              Call TstCor(IMu2Sh+NShell,MDV,'FoFDir-FixR1C')
              Call FixR1C(IPurAO,NMatS,LMaxSP,NShell,ShellT,ShellC,JAn,
     $          R03,V(IP1C2I),V(IP1C2J),V(IMu2Sh))
            else if(DBFInt) then
C             Call FixDBI(NMatS,R03)
            else if(ICntr2.eq.7) then
              If(IOpCl.le.3) Call FixRaf(0,jja,jja,3,NBas6D,R1)
              If(IOpCl.eq.1.or.IOpCl.eq.3.or.IOpCl.eq.4)
     $          Call FixRaf(0,jja,jja,4,NBas6D,R2)
              If(IOpCl.eq.2.or.IOpCl.eq.3.or.IOpCl.eq.4)
     $          Call FixRaf(0,jja,jja,4,NBas6D,R03)
              endIf
            endIf
          endIf
        If(DoRafI) then
          NOp = Min(NOp,NOpAbU)
          AbSym = (JSym2E.eq.1.and.NOp.eq.NOpAbU).or.JSym2E.eq.2
          NAbSym = .False.
C         If we are doing CI singles by Raffenetti, the length of the
C         Fock & density matrices must be NBas6D*(NBas6D+1)
          If(SqFock) then
            LenFP = NBas6D*(NBas6D+1)/2
            NFPMat = 2*NMat
          else
            LenFP  = NDim6D
            NFPMat = NMat
            endIf
          If(ICntr2.eq.6.or.(ICntr2.eq.5.and.IOpCl.eq.0)) then
            Write(IOut,1230) ICntr2, IOpCl
            Call Lnk1E(0)
            endIf
          If(IHMeth.ne.21) Call FoFRaf(IOut,IPrint,IHMeth,IR12Op,IWhat,
     $      IDens,IAcrcy,MinBra,MaxBra,IRICut,CutWrt,IOpCl,0,NBas6D,
     $      NFPMat,NMatS,NMatT,NMatD,LenFP,NAtoms,JSym2E,NOpAbU,NDimEq,
     $      NEqShl,JTrans,PA,PB,FA,FB,C,R1,R2,AccDes,ScaHFX,IPFlag,
     $      AllowP,LMaxSP,V(JBeg),JLeft,V(JBeg),ILeft)
          If(ICntr2.eq.7) then
            If(DoRegI) then
              Write(IOut,1240) DoRegI
              Call Lnk1E(0)
              endIf
            If(IOpCl.le.3) Call FixRaf(0,jja,jja,1,NBas6D,R1)
            If(IOpCl.eq.1.or.IOpCl.eq.3.or.IOpCl.eq.4)
     $        Call FixRaf(0,jja,jja,2,NBas6D,R2)
            If(IOpCl.eq.2.or.IOpCl.eq.3.or.IOpCl.eq.4)
     $        Call FixRaf(0,jja,jja,2,NBas6D,R03)
            endIf
          endIf
      else if(HarHam(Meth)) then
        If(IPrt1.gt.0) Write(IOut,1010)
      else if(SEInt) then
        Call FoFSE(IOut,IPrint,Meth,IOpCl,IDens,DoFock,DoForc,DoFokD,
     $    SqFock,DoFreq,NDim6D,NMat,NMatS,NMatT,NAtoms,ICnBeg,ICnEnd,
     $    NBasis,NBas6D,AddH,InitF,IAn,IAtTyp,C,HA,HB,PA,PB,FA,FB,FDA,
     $    FDB,V(IFTemp+3),V(IFFTmp),V(JBeg),MDV-JBeg+1)
      else
        Call GauErr('Logic error #2 in FoFDir.')
        endIf
      DoPuFL = DoPurF.and..not.SEInt
      If(DoChi) then
        NOpEnd = 1
      else
        NOpEnd = NOp
        endIf
      Call EndFoF(IOut,IPrint,'FoFDir',DoFock,.False.,DoForc,DoFokD,
     $  DoFreq,DoChi,AbSym,NAbSym,DoPuFL,HaveFB,SqFock,AddH,DoAdj,
     $  ISym2E,NBasis,NBas6D,NAtoms,ICnEnd,1,NMat,NMatSO,NMtFkD,.False.,
     $  NOpEnd,NSymOp,NShell,LenFT,LenFFT,NDim6D,IPureD,IPureF,JTrans,
     $  NEqAtm,RotOp,NEqShl,NEqBas,ShellT,ShellC,0,0,JJA,0,JJA,JJA,HA,
     $  HB,V(IFTemp),V(IFFTmp),FA,FB,FDA,FDB,FXYZ,FFXYZ,PeGmx,V,MDV)
      If(DBFInt) then
        Call GenAOS(I5DB2,I7FB2,NShelB,ShelTB,ShelCB,Junk,AOSB)
      else if(.not.SEInt) then
        If(DoP) then
          Call F2Pure(IOut,IPrint,IPureD,IPureF,NMat,NMatST,0,NBasis,
     $      NBas6D,.True.,.True.,.False.,.False.,HaveDB,.True.,PA,PB,XX,
     $      XX,V,MDV)
          If(NMatE.gt.0) then
            Call F2Pure(IOut,IPrint,IPureD,IPureF,NMatE,NMatDS+NMatDT,0,
     $        NBasis,NBas6D,.True.,.True.,.False.,.False.,HaveDB,.True.,
     $        PA(1+IndE),PB(1+IndE),XX,XX,V,MDV)
            Call AMove1(NMatE*NTT,IndE,NMat*NTT,PA)
            If(HaveDB) Call AMove1(NMatE*NTT,IndE,NMat*NTT,PB)
            endIf
          endIf
        If(UseZ) Call F2Pure(IOut,IPrint,IPureD,IPureF,NMat,NMatS,0,
     $    NBasis,NBas6D,.True.,.True.,.False.,.False.,HaveDB,.True.,PZA,
     $    PZB,XX,XX,V,MDV)
        If(DoCISG) Call F2Pure(IOut,IPrint,IPureD,IPureF,NMat,NMatS,0,
     $    NBasis,NBas6D,.True.,.True.,.False.,.True.,.True.,.True.,
     $    CIDenA,CIDenB,XX,XX,V,MDV)
        endIf
      Called = .True.
      Return
      End
