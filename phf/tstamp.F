*Deck TStamp
      Subroutine TStamp(IOper,String)
      Implicit Real*8(A-H,O-Z)
C
C     Print String with a time stamp.  Operations:
C     1 ... Print string and reset clock.
C     2 ... Print string with date and current cpu but don't reset clock.
C     3 ... Print string only.
C     4 ... Reset timer to unset; do not print anything.
C     5 ... Reset timer to current time; do not print anything.
C     6 ... Reset timer to current time only if using OpenMP; do not print.
C
      Parameter (IBTSt=14,IBPTst=24)
      Logical LBit1
      Common /IO/ In, IOut, IPunch
      Character DayTim*24, String*(*), Temp*20
      Common /TStampC/ CPULst
C$OMP ThreadPrivate (/TStampC/)
      Save Zero, Small
      Data Zero/0.0d0/, Small/1.d-10/
 1000 Format(1X,A,1X,A)
 1010 Format(1X,A,1X,A,F12.2,' cpu seconds.')
 1015 Format(1X,A,1X,A,F12.2,' cpu seconds (not reset).')
 1020 Format(1X,A)
 1030 Format(' Illegal IOper=',I12,' in TStamp.')
C
      If(LBit1(IBTSt)) then
        IThNum = MDPart(-1)
        If(IThNum.ge.0.and..not.LBit1(IBPTSt)) Return
        LS = Max(LinEnd(String),1)
        T = CPUTim(0)
        TT = T - CPULst
        Call GDate(DayTim)
        If(IOper.eq.1) then
          Temp = String(1:LS)
          If(CPULst.eq.Zero) then
            Write(IOut,1000) Temp, DayTim
          else
            Write(IOut,1010) Temp, DayTim, TT
            endIf
          CPULst = Max(T,Small)
        else if(IOper.eq.2) then
          If(CPULst.eq.Zero) then
            Write(IOut,1000) String(1:LS), DayTim
          else
            Write(IOut,1015) String(1:LS), DayTim, TT
            endIf
        else if(IOper.eq.3) then
          Write(IOut,1020) String(1:LS)
        else if(IOper.eq.4) then
          CPULst = Zero
        else if(IOper.eq.5) then
          CPULst = Max(T,Small)
        else if(IOper.eq.6) then
          If(IThNum.gt.0) CPULst = Max(T,Small)
        else
          Write(IOut,1030) IOper
          Call Lnk1E(0)
          endIf
        endIf
      Return
      End
