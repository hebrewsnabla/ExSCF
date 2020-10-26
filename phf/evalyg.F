

      subroutine evalyg (iout, iprint, nci, ngrdt, yg, smt, civec)

      implicit none


C +------------------------------------------------------------+
C |  evalyg   --   CAJH, 07.2011                               |
C |                                                            |
C |                                    (based on PHFB evalyg)  |
C |                                                            |
C |                                                            |
C |  Evaluate the function y(g) from the values x(g) and the   |
C |  CI coefficients stored in civec. The weighting function   |
C |  y(g), appearing ubiquitously in the PHF expressions, is   |
C |  given by                                                  |
C |                                                            |
C |                          c_i* c_j  x(g; i,j)               |
C |    y(g; i,j) =  ---------------------------------------,   |
C |                  sum_{k,l} c_k* c_l  int dg  x(g; k,l)     |
C |                                                            |
C |  where there is a parametric dependence on the i and j     |
C |  indices in the functions x(g) and y(g).                   |
C |                                                            |
C |  Here, c_i are CI coefficients obtained from solving the   |
C |  CI problem described in solvci, that is, the              |
C |  diagonalization of the Hamiltonian matrix in the basis    |
C |  basis of the nci states.                                  |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  In this routine, we use the fact that the overlap matrix  |
C |  elements are given by                                     |
C |                                                            |
C |    Smt(i,j) = int dg x(g; i,j)                             |
C |                                                            |
C |  to compute the denominator appearing in the expression    |
C |  defining y(g).                                            |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       nci    - number of basis states in CI problem
C       ngrdt  - total number of grid points
C       yg     - vector of x(g) / y(g) values [ in / out ]
C       smt    - overlap matrix in CI problem
C       civec  - vector of CI coefficients

      integer     iout, iprint, nci, ngrdt
      complex*16  yg(*), smt(nci,nci), civec(*)

C     other variables

      integer     i, j, n1, n2, igrdt, ioffy, indy
      complex*16  intxg
      character   indc1*2, indc2*2, str*10



C     Compute the quantity

C       intxg = sum_{i,j} c_i* c_j  int dg x(g; i,j).

      intxg = cmplx (0.0d0, 0.0d0)

      do 101 i = 1, nci
        do 102 j = 1, nci
          intxg = intxg + smt(i,j) * conjg (civec(i)) * civec(j)
 102    continue
 101  continue



C     %%%%%%%%%%%%%%%%%%%
C     %  Evaluate y(g)  %
C     %%%%%%%%%%%%%%%%%%%


C     Compute y(g; i,j) from x(g; i,j) as described above...


C     Loop over CI matrix elements.

      do 201 n1 = 1, nci
        do 202 n2 = 1, nci

          ioffy = ((n1-1)*nci + n2-1) * ngrdt

C     Evaluate y(g).

          do 203 igrdt = 1, ngrdt
            indy = ioffy + igrdt

            yg(indy) = yg(indy) * conjg (civec(n1)) * civec(n2) / intxg
 203      continue


C     Print vector of y(g) if requested.

          if ( iprint .gt. 3 ) then
            write (indc1, '(I2)') n1
            write (indc2, '(I2)') n2
            str = ' (g,' // indc1 // ',' // indc2 // ')'

            write (iout, *) ' '
            call ltoutx (iout, 3, 2, 1, 0, 'y' // str, ngrdt, 1,
     $           yg(ioffy+1), yg(ioffy+1), 0)
          endif
 202    continue
 201  continue


      return
      end


