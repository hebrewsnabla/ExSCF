

      subroutine erictr (iout, iprint, v, lenv, nbasis, iopcl, nmts,
     $     nmtt, nmta, xa, xb, ya, yb, ipflag, allowp, jsym2e, nsymop,
     $     nop1, neqbas, neqshl, neqsh2, rotop, neqatm, accdes)

      implicit none


C +------------------------------------------------------------+
C |  erictr   --   CAJH, 07.2011                               |
C |                                                            |
C |                                   (based on PHFB version)  |
C |                                                            |
C |                                                            |
C |  This is the routine that calls PRISM (FoFDir), and        |
C |  should be called by formgg. It contracts P(g) with        |
C |  electron repulsion integrals.                             |
C |                                                            |
C |  The xa and xb arrays contain, on input, the alpha and     |
C |  beta density matrices to be contracted (they are          |
C |  DESTROYED). The ya and yb arrays are returned on output   |
C |  with the two-electron matrices. ALL MATRICES ARE STORED   |
C |  IN LOWER TRIANGULAR FORM.                                 |
C |                                                            |
C |  The integers nmts, nmtt, and nmta control the number of   |
C |  singlet, triplet, and antisymmetrices matrices,           |
C |  respectively. (Singlet matrices are contracted with       |
C |  Coulomb and exchange integrals, whereas triplet           |
C |  matrices are contracted only with exchange integrals).    |
C |                                                            |
C |  Note, however, that all matrices should be packed into    |
C |  the arrays xa, xb, with the following structure:          |
C |                                                            |
C |    . 1 ... nmts matrices should be singlet                 |
C |    . nmts+1 ... nmts+nmtt should be triplet                |
C |    . nmts+nmtt+1 ... nmts+nmtt+nmta should be antisymm     |
C |                                                            |
C |  All the logic of converting the output arrays ya and      |
C |  yb into G(g) should be done elsewhere.                    |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  A final cautionary note. The xa, xb, ya, and yb matrices  |
C |  are packed as a series of (nmts+nmtt+nmta) matrices,      |
C |  each of size NTT = nbasis*(nbasis+1)/2. However, the      |
C |  arrays xa, xb, ya, and yb should be allocated with size   |
C |  (nmts+nmtt+nmta)*NTT6D, where                             |
C |                                                            |
C |    NTT6D = nbas6d*(nbas6d+1)/2.                            |
C |                                                            |
C |  This is a PRISM requirement...                            |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       nbasis - number of basis functions (pure)
C       iopcl  - IOpCl (see FoFDir for details)
C       nmts   - number of singlet matrices to contract
C       nmts   - number of triplet matrices to contract
C       nmta   - number of antisymmetric matrices to contract
C       xa     - set of alpha density matrices
C       xb     - set of beta density matrices
C       ya     - set of alpha Hamiltonian matrices [ out ]
C       yb     - set of beta Hamiltonian matrices [ out ]

      real*8      v(*), accdes
      real*8      xa(*), xb(*), ya(*), yb(*)
      integer     iout, iprint, lenv
      integer     nbasis, iopcl, nmts, nmtt, nmta

C       .. prism variables .. [ input variables ]

      real*8      rotop(*)
      integer     ipflag, nsymop, nop1, jsym2e
      integer     neqbas(nbasis,*), neqshl(*), neqsh2(*), neqatm(*)
      logical     allowp(50)

C     other variables

      real*8      junk
      integer     icntrl
      integer     nmat, nmats, nmatt, nmatd
      logical     initf, dopurf


      junk = 0.0d0


C     Set up some options before calling FoFDir.

      icntrl = 0
      initf  = .true.
      dopurf = .true.

      nmat  = nmts + nmtt + nmta
      nmats = nmts
      nmatt = nmtt
      nmatd = nmat


C     Call FoFDir.

      call fofdir (iout, iprint, 1, iopcl, icntrl, 1, ipflag, allowp, 0,
     $     0, .false., initf, dopurf, accdes, 1.0d0, 0, nmat, nmats,
     $     nmatt, nmatd, 0, 0, nbasis, jsym2e, nsymop, nop1, neqatm,
     $     neqshl, rotop, neqbas, junk, junk, junk, junk, xa, xb, junk,
     $     junk, ya, yb, junk, junk, 0, 0, junk, .false., 0, .false.,
     $     junk, junk, 0, 0, junk, junk, junk, v, lenv)


      return
      end


