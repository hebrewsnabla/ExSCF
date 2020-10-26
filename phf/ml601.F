*Deck ML601
#ifdef GAUSS_LINDA
      Subroutine Real_Main
#else
      Program ML601
#endif
      Implicit Real*8(A-H,O-Z)
      Dimension Work(1)
      Common /IO/ In, IOut, IPunch
      Integer GSetJmp
#ifdef _ALIGN_CORE_
      Common /GXWork/ Work
#endif
#ifdef P64
      Integer*8 IOff
#endif
C
      In = 5
      IOut = 6
      IPunch = 7
      IStat = GSetJmp(0)
      If(IStat.eq.0) then
        Call InitSCM(1,0,0,Work,IOff,MDV)
        Call muldrv(Work(IOff+1),MDV)
      else
        Call PrtStat(IStat,'ML601')
        endIf
      End
