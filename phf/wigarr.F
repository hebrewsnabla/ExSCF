      

      subroutine wigarr (iout, iprint, s, isint, na, nb, ncis, ngrdb,
     $     grdb, dmt, dind)

      implicit none


C +------------------------------------------------------------+
C |  wigarr   --   CAJH, 07.2011                               |
C |                                                            |
C |                                                            |
C |  Prepare the arrays dmt and dind needed by evalsmt in      |
C |  PHF calculations.                                         |
C |                                                            |
C |  The array dind is the collection of {M} values that       |
C |  define the CI problem associated with spin restoration.   |
C |  Note that if the number of CI vectors = 1, then there     |
C |  array is of dimension 1.                                  |
C |                                                            |
C |    Let's consider an example. Suppose we are projecting    |
C |    onto s = 1/2, and we want to do full spin projection.   |
C |    Then, dind should store                                 |
C |                                                            |
C |      dind(1) =  0.50,                                      |
C |      dind(2) = -0.50.                                      |
C |                                                            |
C |  The array dmt collects the values of Wigner's small       |
C |  rotation matrix d^s_{MK} (beta). It is thus a three-      |
C |  dimensional array. In particular, for every combination   |
C |  of M and K it stores the value of d^s_{MK} (beta) by      |
C |  looping over all grid points associated with integration  |
C |  over beta.                                                |
C |                                                            |
C |  The array dmt is filled using the function wignerd,       |
C |  which yields explicitly d^s_{MK} (beta) for given s,      |
C |  M, K, and beta.                                           |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       s     - quantum number s on which to project
C       isint - .true. for integer s, .false. for half-integer s
C       na    - number of alpha electrons
C       nb    - number of beta electrons
C       ncis  - number of basis states for spin projection
C       ngrdb - number of grid points in beta-integration array
C       grdb  - integration grid for variable beta
C       dmt   - Wigner's small d-matrix d^S_{MK} (beta)  [ out ]
C       dind  - vector with indices {M} for Wigner matrices  [ out ]

      integer     iout, iprint, s, na, nb, ncis, ngrdb
      real*8      grdb(*), dmt(ngrdb,ncis,ncis), dind(ncis)
      logical     isint

C     other variables

      integer     n1, n2, i1, i2, j, sz
      real*8      wignerd


C     Prepare array dind.

      if ( ncis .eq. 1 ) then
        dind(1) = (dble (na) - dble (nb)) / 2.0d0
      else

        n1 = s

        do 10 i1 = 1, ncis
          if ( isint ) then
            dind(i1) = dble (n1)
            n1 = n1 - 1
          else
            dind(i1) = dble (n1) / 2.0d0
            n1 = n1 - 2
          endif
 10     continue
      endif


C     Now prepare the array dmt containing d^s_{MK} (beta).

C       ** Note that when ispuhf = .true., we only consider
C          d^s_{N0,N0} (beta), where N0 = (Na - Nb) / 2.

      if ( ncis .eq. 1 ) then

        if ( isint ) then
          sz = ( na - nb ) / 2
        else
          sz = ( na - nb )
        endif

        do 20 j = 1, ngrdb
          dmt(j,1,1) = wignerd (isint, s, sz, sz, grdb(j))
 20     continue

      else

        n1 = s

        do 21 i1 = 1, ncis
          n2 = s

          do 22 i2 = 1, ncis
            do 23 j = 1, ngrdb
              dmt(j,i1,i2) = wignerd (isint, s, n1, n2, grdb(j))
 23         continue

            if ( isint ) then
              n2 = n2 - 1
            else
              n2 = n2 - 2
            endif
 22       continue

          if ( isint ) then
            n1 = n1 - 1
          else
            n1 = n1 - 2
          endif
 21     continue
      endif

      write(iout, 1000) s, sz, isint, ncis, dind
      write(iout, 1001) dmt
1000  format('WigArr: s=', I5, ', sz=', I5, ', isint=', L5, ', ncis=', I5, 
     $  ', dind=', F10.4)
1001  format('dmt=', 8F10.4)

      return
      end


