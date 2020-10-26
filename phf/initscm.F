*Deck InitSCM
      Subroutine Initscm(IArg,MinSp,MaxSp,Work,IOff,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Allocate space using Work as a "pointer" and returning
C     the offset in IOff.  MinSp and MaxSp if non-zero specify
C     the minimum and maximum amounts to request.  If IArg is
C     non-zero, the IArg'th argument on the command line is
C     used to determine the basic value;  otherwise a default is used.
C     IOff receives the offset into Work for the start of the
C     allocated array and MDV receives the length in WP words.
C
C     This routine calls routine GAlloc, which has the
C     following calling sequence:
C
C     GAlloc (NWords, WordSize, IClear, RefArray, OffSet)
C
C     NWords   = (integer) # of words to allocate.
C     WordSize = (integer) size of words (in bytes).
C     IClear   = (integer) 1 to clear allocated memory.
C     RefArray = Array to reference created memory, of elements
C                of size WordSize (i.e., Real*8 if WordSize=8).
C     OffSet   = (integer) offset such that RefArray(OffSet+1)
C                is at Address.  Note that OffSet seems to actually
C                be returned such that RefArray(OffSet) is at
C                Address under 4.3bsd, but this routine and GetSCM
C                will work either way.
C
C     An adequate GAlloc (which works as documented, not as
C     implemented!) is provided in mdutil.c.
C
      Parameter (MemDef=(32*1024*1024))
      Logical GetEnI
      Integer GetCmN
      Dimension Work(1)
      Common /MaxMem/ MaxMem(2)
#ifdef P64
      Integer*8 IOff
#endif
C
      MaxMem(2) = 0
      If(IArg.gt.0) MaxMem(2) = GetCmN(IArg)
      If(MaxMem(2).eq.0) then
        If(GetEnI('GAUSS_MEMDEF',MemEnv)) then
          MaxMem(1) = MemEnv
        else
          MaxMem(1) = MemDef
          endIf
        If(MinSp.ne.0.and.MaxMem(1).lt.MinSp) MaxMem(1) = MinSp
        If(MaxSp.ne.0.and.MaxMem(1).gt.MaxSp) MaxMem(1) = MaxSp
      else
        MaxMem(1) = IAbs(MaxMem(2))
        endIf
      If(MaxMem(2).eq.0) MaxMem(2) = MaxMem(1)
      Call GAlloc(MaxMem(1), 8, 0, Work, IOff)
#ifdef INIT_SCM
      Call SetNaN(Work(IOff+1),MaxMem(1))
#else
      If(GetEnI('GAUSS_INITVAL',Junk))
     $  Call SetNaN(Work(IOff+1),MaxMem(1))
#endif
      MDV = MaxMem(1)
      Return
      End
