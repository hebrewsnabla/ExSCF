*Deck MulDrv
      Subroutine MulDrv(V,NGot)
      Implicit Real*8(A-H,O-Z)
C
C     Copyright (c) 1988,1990,1992,1993,1995,1998,2003,2009, Gaussian, Inc.
C                   All Rights Reserved.
C
C     This is part of the private, development version of the
C     the Gaussian(R) DV system of programs. It is based on
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
C1OV6
C
C     Overlay 6 is concerned with the computation of properties using
C     the wavefunctions.
C
C     Link 601:  Mulliken population analysis.
C
C     Link 602:  Determination and mapping of the electrostatic
C                potential, electric field, electric field gradient,
C                and charge multipole moments.
C
C     Link 603:  Determination of density basis expansion based on
C                electrostatic properties.
C
C     Link 604:  Evaluation of the density, density gradient, and
C                density divergence over a 2d or 3d grid of points.
C
C     Link 605:  Prepares input for RPAC.
C
C     Link 606:  RPAC.
C
C     Link 607:  Natural population and bond orbital analysis (NBO).
C
C     Link 608:  Density functional energies using densities from
C                traditional methods.
C
C2IOp(5)
C     IOp(5) ... No longer used.
C
C2IOp(6)
C     IOp(6) ... No longer used.
C
C2IOp(7)
C     IOp(7) ... Printing of MOs.
C          0 ... Default:  1 for molecules, 2 for PBC.
C          1 ... Print the occupied and first 5 virtual MOs.
C          2 ... Do not print any MOs.
C          3 ... Print all MOs.
C         10 ... Biorthogonalize unrestricted MOs.
C
C2IOp(8)
C2IOp(9)
C2IOp(10)
C2IOp(11)
C2IOp(12)
C          These options are print/no-print options.  The
C     possible values are:
C
C     0 ... DEFAULT.
C     1 ... PRINT THE NORMAL AMOUNT.
C     2 ... DO NOT PRINT.
C     3 ... PRINT VERBOSELY.
C
C     IOp(8) ... DENSITY MATRIX.  DEFAULT:  NO-PRINT.
C
C     IOp(9) ... FULL POPULATION ANALYSIS.  DEFAULT:  PRINT.
C
C     IOp(10) ... Gross orbital charges.  Default:  Print.
C
C     IOp(11) ... GROSS ORBITAL TYPE CHARGES.  DEFAULT:  NO-PRINT.
C
C     IOp(12) ... CONDENSED TO ATOMS.  DEFAULT:  PRINT.
C
C2IOp(13)
C     IOp(13) ... Whether to save computed electric field on disk
C                 for use in Tomasi RF calculations.
C           0 ... Default (No).
C           1 ... Yes.
C           2 ... No.
C
C2IOp(14)
C     IOp(14) ... SPECIFICATION OF OTHER PROPERTIES TO BE CALCULATED.
C     0 ... DEFAULT.  EVALUATE THE ELECTRIC POTENTIAL, THE ELECTRIC
C              FIELD, AND THE ELECTRIC FIELD GRADIENT AT EACH CENTER.
C     1 ... EVALUATE THE ELECTRIC POTENTIAL, THE ELECTRIC FIELD, AND
C              THE ELECTRIC FIELD GRADIENT AT EACH CENTER.
C     2 ... EVALUATE THE POTENTIAL AND THE ELECTRIC FIELD AT EACH
C              CENTER.
C     3 ... EVALUATE ONLY THE POTENTIAL AT EACH CENTER.
C     4 ... EVALUATE NONE.
C
C2IOp(15)
C     IOp(15) ... SPECIFICATION OF ADDITIONAL CENTERS.  IF MORE THAN
C                 ONE OF THESE IS REQUESTED, THE LISTS ARE IN SEPARATE
C                 INPUT SECTIONS IN THE ORDER LISTED BELOW.
C     0 ... NO ADDITIONAL CENTERS.  EVALUATE THE PROPERTIES ONLY AT
C              EACH ATOMIC CENTER.
C     1 ... READ ADDITIONAL CENTERS.  ONE CARD PER CENTER WITH THE
C              X, Y AND Z COORDINATES IN ANGSTROMS (FREE FORMAT).
C     2 ... READ IN COORDINATES AS FOR 1.  STARTING AT EACH POINT,
C              LOCATED THE NEAREST STATIONARY POINT IN THE ELECTRIC.
C              POTENTIAL.
C     4 ... Read in a set of cards specifying a grid of points at which
C              the electric potential will be computed.  Two forms of
C              specifications are allowed:
C
C          A.  Evenly spaced rectangular grid.
C              Three cards are required:
C              KTape,XO,YO,ZO   --- output unit and coordinates of one
C                                   corner of grid.  If KTape is 0, it
C                                   defaults to 51.
C              N1,X1,Y1,Z1      --- number of increments and vector.
C              N2,X2,Y2,Z2      --- number of increments and vector.
C              N1 records will be written to unit KTape, with N2
C              values in each record.
C
C          B.  An arbitrary list of points.
C              Only one card is needed:
C                     N,NEFG,LTape,KTape
C              The coordinates of N points in Angstroms will be read
C              unit LTape in format (3F20.12).  The potential (NEFG=3),
C              potential and field (NEFG=2), or potential, field, and
C              field gradient (NEFG=1) will be computed and written
C              along with the coordinates to unit KTape in format
C              (4F20.12).  Thus if NEFG=3 for each point there will be
C              4 cards written per point, containing:
C                     X-coord,Y-coord,Z-coord,Potential
C                     X-field,Y-field,Z-field,XX-EFG
C                     YY-EFG,ZZ-EFG,XY-EFG,XZ-EFG
C                     YZ-EFG
C
C          Note that either form of grid should be specified with respect
C          to the standard orientation of the molecule.
C
C     8 ... Do potential-derived charges.
C    16 ... Constrain the dipole in fitting charges.
C    32 ... Read in centers at which to evaluate the potential from the
C           rwf.
C   128 ... Read grid; do not default cube.
C
C2IOp(16)
C     IOp(16) ... Cutoffs in L602.
C     0 ... Use full accuracy in calculations at specific points,
C           but use sleazy cutoffs in mapping a grid of points.
C     1 ... Do all points to full accuracy.
C
C2IOp(17)
C     IOp(17) ... DEBUGGING CONTROL (L602).
C     0 ... COMPUTE ALL CONTRIBUTIONS TO SELECTED PROPERTIES.
C     1 ... COMPUTE ONLY THE NUCLEAR CONTRIBUTION.
C     2 ... COMPUTE ONLY THE ELECTRONIC CONTRIBUTION.
C    -N ... COMPUTE ONLY THE CONTRIBUTION OF SHELL N.
C
C2IOp(18)
C
C     Whether to update dipole rwf 0/1 yes/no.
C
C2IOp(19)
C
C     Whether to rotate exact polarizability before comparing
C     with approximate (which will be calculated in the standard
C     orientation).  This is like IOp(9) in L9999.
C
C     0 ... Default, same as 1.
C     1 ... Exact is still in standard orientation; use as-is.
C     2 ... Exact is already in z-matrix orientation, so rotate.
C
C2IOp(20)
C
C     How to do electrostatic-potential derived charges:
C
C     0 ... Default (1).
C    -1 ... Read a list of points at which to fit, one per line.
C     1 ... Merz-Kollman point selection
C     2 ... CHELP point selection.
C     3 ... CHELPG point selection.
C    00 ... Default radii are those defined with the selected method.
C    10 ... Force Merz-Kollman radii.
C    20 ... Force CHELP (Francl) recommended radii.
C    30 ... Force CHELPG (Breneman) recommended radii.
C    40 ... Force 2xUFF Radii.
C   100 ... Read in replacement radii for selected atom types as pairs
C           (IAn,Rad) or (Symbol,Rad), terminated by a blank line.
C   200 ... Read in replacment radii for selected atoms as pairs
C           (I,Rad), terminated by a blank line.
C  1000 ... Fit united atoms (heavy atoms only) rather than all atoms.
C 10000 ... Use only active atoms in the fit.
C
C2IOp(21)
C
C     Operation of L603:
C     0 ... Default (same as 2).
C     1 ... Read in density basis functions and compute populations.
C     2 ... Optimize density basis set.
C
C2IOp(22)
C
C     Selection of density matrix (currently only in L601, L602, L604):
C
C    -1x ... Read density matrices from .chk file.
C    +1x ... Read density matrices from .chk file.
C     -5 ... All available transition densities.
C     -4 ... Transition density between the states given by IOp(29) and
C            IOp(30).
C     -3 ... Density for the excited state given by IOp(29).
C     -2 ... Use all available density matrices.
C     -1 ... Use the density matrix for the current method, or the HF
C            density if the one for the current method is not available.
C N.ge.0 ... Use the density matrix for method N (see Link 1 for the
C            numbering scheme).
C
C2IOp(23)
C
C     Density values to evaluate over grid in L604:
C     0 ... Default (same as 3).
C     1 ... Density values.
C     2 ... Density values and gradients.
C     3 ... Density values, gradients and divergence.
C
C2IOp(24)
C
C     IOp(24) ... Frozen core:
C          -N ... Freeze N orbitals.
C           0 ... Default (Yes).
C           1 ... Yes.
C           2 ... No.
C
C2IOp(25)
C
C     IOp(25) ... Whether to compute coulomb self-energy in L601:
C           0 ... No
C           1 ... Yes, classically (including self terms -- requires
C                 2e integrals, O(N**4)).
C           2 ... Yes, quantum mechanically (no self terms -- requires
C                 2e integrals,  and only available for HF.  O(N**5)).
C
C2IOp(26)
C     IOp(26) ... Which density to use in L602 and L604:
C           0 ... Default (same as 1)
C           1 ... Total
C           2 ... Alpha
C           3 ... Beta
C           4 ... Spin
C
C2IOp(27)
C     IOp(27) ... Choice of population analysis:
C           0 ... Default (12).
C           1 ... Don't do Mulliken populations
C           2 ... Do Mulliken populations
C          10 ... Don't do bonding Mulliken Populations.
C          20 ... Do bonding Mulliken Populations.
C         100 ... Do Minimal population analysis.
C        1000 ... Read in weightings for atoms pairs for unequally
C                 split Mulliken.
C
C2IOp(28)
C     IOp(28) ... Mark SCF density as current density.
C           0 ... No: save SCF density, but do not mark.
C           1 ... Yes: mark as well.
C
C2IOp(29)
C     IOp(29) ... Excited state to use if requested by IOp(22).
C
C2IOp(30)
C     IOp(30) ... 2nd excited state for transition density:
C           0 ... Transition denstiy between state IOp(29) and g.s.
C           N ... Transition denstiy between state IOp(29) and state N.
C
C2IOp(31)
C     IOp(31) ... Whether to determine natural orbitals from densities:
C           0 ... No.
C           1 ... Yes, using total density.
C           2 ... Yes, using alpha and beta separately for UHF.
C           3 ... Store only alpha NOs.
C           4 ... Store only beta NOs.
C           5 ... Use spin density.
C
C2IOp(32)
C     IOp(32) ... CONTROL PARAMETERS FOR COVBON in L609 (NOT TO BE
C                 CHANGED UNDER MOST CIRCUMSTANCES):
C           10000*MItLoc+1000*ITlLoc+100*IDcInt+IPrLoc, where
C           MItLoc ... MItLoc*NOrb*(NOrb-1)/2 is the maximum number
C                      of iterations in localization of (spin)orbitals
C                      (1...9, default 6),
C           ITlLoc ... 10.**(-ITlLoc) is the convergence criterion
C                      for (spin)orbital localization
C                      (1...9, default 9),
C           IDcInt ... Localized (spin)orbitals with atomic occupancies
C                      less than 0.01*IDcInt are interpreted as lone pair
C                      MOs rather than bond MOs (1...99, default 10),
C           IPrLoc ... 0: Print the atomic occupancies of localized
C                         (spin)orbitals (default),
C                      1: Do not print the atomic occupancies.
C
C     IOp(32) ... L605, L606: naming of RPAC interface file.
C           0 ... Make this a scratch file.
C           1 ... Name this file 'rpac.11'
C
C2IOp(35)
C     IOp(35) ... WHAT TO DO:
C           0 ... Determine attractors, attractor interaction lines,
C                 ring points, and cage points.
C           1 ... Determine zero-flux surfaces (IDoZrF).
C           2 ... Compute charges of AIMs (IDoAtC).
C           4 ... Compute kinetic energies and multipole moments
C                 of AIMs (IDoPrp).
C          10 ... Compute energies of electrostatic interactions
C                 between AIMs (IDoPot).  This precludes calculations
C                 of atomic property derivatives with respect to nuclear
C                 displacements.
C          20 ... ??? currently unused ???.
C          40 ... ??? currently unused ???.
C         100 ... Compute atomic overlap matrices (IDoAOM).
C         200 ... Compute other atomic matrix elements (IDoAMa).
C         400 ... Include zero-flux surface relaxation terms in all atomic
C                 matrix elements (IDoSRe)
C        1000 ... Compute derivatives of atomic properties with respect
C                 to electric field (IDoSeP).  Note that IDoSRe should be
C                 set to 1 in order to obtain correct results!  Also note
C                 that analytical polarizabilities have to be available
C                 but force constants have to be absent!
C        2000 ... Compute derivatives of atomic properties with respect
C                 to nuclear displacements as well (IDoNuD).  Note that
C                 analytical force constants have to be available!
C        4000 ... ??? currently unused ???.
C       10000 ... Compute localized orbitals and bond orders (IDoLoc).
C       20000 ... Compute atomic orbitals in molecule (IDoAOs).
C       40000 ... ??? currently unused ???.
C      100000 ... If necessary, augment valence electron densities
C                 with relativistic core contributions, which is a default
C                 anyway (IHwAug=0).
C      200000 ... If necessary, augment valence electron densities
C                 with nonrelativistic core contributions (IHwAug=1).
C      400000 ... Abort if pseudopotentials have been used (IHwAug=3).
C     1000000 ... Reduce accuracy so atomic charges can be computed more
C                 rapidly (IQuick).  No other properties can be calculated.
C                 This option sets IPrNDe=5, IPrNAt=5, and IEpsIn=100.
C     2000000 ... Use numerical instead of analtyic integration.
C     3000000 ... Use numerical instead of analtyic integration and use reduced
C                 cutoffs.
C
C2IOp(36)
C     IOp(36) ... CONTROL PARAMETERS FOR NEGLECT OF ORBITALS AND PRIMITIVES
C                 in L609:
C           10000*INoZer+100*IPrNDe+IPrNAt, where
C           INoZer ... 0: Ignore (spin)orbitals with zero occupancies
C                      (default),
C                      1: Do not ignore (spin)orbitals with zero occupancies,
C           IPrNDe ... Neglect primitive contributions below 10.**(-IPrNDe)
C                      in evaluations of electron density and its derivatives
C                      (0...99, default 7),
C           IPrNAt ... Neglect primitive contributions below 10.**(-IPrNAt)
C                      in integrations over atomic basins (0...99, default 7).
C
C2IOp(37)
C     IOp(37) ... CONTROL PARAMETERS FOR ATINLI, RNGPNT, AND CAGPNT in L609
C                 (NOT TO BE CHANGED UNDER MOST CIRCUMSTANCES):
C           1000000*MxBpIt+100000*SBpMax+1000*NGrd+LookUp, where
C           MxBpIt ... Maximum number of iterations in trial path
C                      determination (1...99, default 10),
C           SBpMax ... Maximum value of the control sum (1...9, default 2),
C           NGrd ..... Length of Fourier expansion for the trial path
C                      (1...99, default 20),
C           LookUp ... Number of grid points in critical point search
C                      (1...999, default 100).
C
C2IOp(38)
C     IOp(38) ... CONTROL PARAMETERS FOR ZRFLUX AND OIGAPI in L609 (NOT TO BE
C                 CHANGED UNDER MOST CIRCUMSTANCES):
C           100000*INStRK+10000*IHowFa+1000*IGueDi
C           +100*IPraIn+10*IRScal+IRtFSe
C           INStRK ... 10*INStRK is the number of steps in the Runge-Kutta
C                      integrations along gradient paths (1...9, default 2),
C           IHowFa ... IHowFa is the maximum distance in the Runge-Kutta
C                      integrations along gradient paths (1...9, default 5),
C           IGueDi ... 10.**(-IGueDi) is the initial displacement from the
C                      critical point in the Runge-Kutta integrations
C                      (1...9, default 6),
C           IPraIn ... 10.*IPraIn is the cut-off for zero-flux surfaces
C                      (1...9, default 2),
C           IRScal ... IRScal is the scaling factor in the nonlinear
C                      transformation used in the intersection search
C                      (1...9, default 2),
C           IRtFSe ... 10.*IRtFSe is the safety factor used in the
C                      intersection search (1...9, default 2),
C2IOp(39)
C     IOp(39) ... More CONTROL PARAMETERS FOR ZRFLUX AND OIGAPI in L609
C                 (NOT TO BE CHANGED UNDER MOST CIRCUMSTANCES):
C           1000000*IToler+100000*INInGr+10000*INInCh+1000*IEpsSf
C           +10*IEpsIn+INTrig
C           IToler ... 10.**(-5-IToler) is the tolerance for the
C                      intersection search (1...9, default 5),
C           INInGr ... 10*INInGr is the initial number of grid points
C                      in theta and phi in the adaptive integration
C                      subroutine (1...9, default 2),
C           INInCh ... 5+INInCh is the initial number of sampling points
C                      in the intersection search (1...9, default 2),
C           IEpsSf ... IEpsSf is the safety factor used for patches with
C                      surface faults in the adaptive integration
C                      subroutine (1...9, default 6),
C           IEpsIn ... 0.0001*IEpsIn is the target for integration error
C                      (1...99, default 2),
C           INTrig ... 10*INTrig is the number of sine and cosine functions
C                      in the trial function for surface sheets (1...9,
C                      default 2).
C
C2IOp(40)
C
C     IOp(40) ... Control of link 607.
C           0 ... Default NBO analysis -- don't read input.
C           1 ... Read input data to control NBO analysis.
C           2 ... Delete selected elements of NBO Fock matrix
C                 and form a new density, whose energy can then
C                 be computed by one of the SCF links.  This link
C                 must have been invoked with IOp(40) = 0 or 1
C                 prior to invoking it with IOp(40)=2.
C           3 ... Read the deletion energy produced by a previous
C                 run with IOp(40)=2 and print it.
C
C2IOp(41)
C
C     IOp(41) ... Number of layers in esp charge fit.
C           0 ... Default (4).
C           N ... N layers, must be >=4.
C
C2IOp(42)
C
C     IOp(42) ... Density of points per unit area in esp fit.
C           0 ... Default (1).
C           N ... points per unit area.
C
C2IOp(43)
C
C     IOp(43) ... Increment between layers in MK charge fit.
C           0 ... Default (0.4/Sqrt(#layers))
C           N ... 0.01*N.
C
C2IOp(44)
C
C     IOp(44) ... Type of calculation in L604:
C           0 ... Default, same as 2.
C           1 ... Compute the molar volume
C           2 ... Evaluate the density over a cube of points
C           3 ... Evaluate MO's over a cube of points
C          10 ... Skip header information in cube file.
C
C2IOp(45)
C
C     IOp(45) ... Number of points per bohr**3 for Monte-Carlo
C                 calaulation of molar volume
C          -1 ... Read from input
C           0 ... Default (20)
C           N ... N points -- for tight accuracy, 50 is recommended.
C2IOp(46)
C
C     IOp(46) ... Threshold for molecular volume integration.
C           0 ... Default -- 10**-3
C          -1 ... Read from input.
C           N ... N*10**-4.
C
C2IOp(47)
C
C     IOp(47) ... Scale factor to apply to van der Waals radii for the
C                 box size during volume integration:
C           0 ... Default.
C           N ... N*0.01 -- for debugging.
C
C2IOp(48)
C
C     IOp(48) ... Use of cutoffs
C           0 ... Default (10**-6 accuracy for cubes, 1 digit better
C                 than desired acuracy for volumes).
C           N ... 10**-N
C
C2IOp(49)
C     IOp(49) ... Approximate number of points per side in cube in l602/l604:
C           0 ... Default (80)
C           N ... N points
C          -1 ... Read from cards.
C          -2 ... Coarse grid, 3 points/Bohr.
C          -3 ... Medium grid, 6 points/Bohr.
C          -4 ... Fine grid, 12 points/Bohr.
C        -N>4 ... Grid using 1000 / N points/Bohr.
C
C2IOp(51)
C     IOp(51) ... Whether to apply Extended Koopman's Theorem (EKT):
C           0 ... Default (No).
C           N ... Yes, on non-SCF densities, up to N IPs and EAs.
C          -1 ... Yes, on non-SCF densities, all possible IPs and EAs.
C          -2 ... No.
C
C2IOp(52)
C     IOp(52) ... Number of radial integration points in L609:
C           0 ... Default (100).
C           N ... N.
C
C2IOp(53)
C     IOp(53) ... Distribution of radial points in L609:
C           0 ... Default (cubic)
C           N ... Polynomial of order N.
C
C2IOp(54)
C     IOp(54) ... Maximum number of domains.
C           0 ... Default (100000).
C           N ... N.
C
C2IOp(55)
C     IOp(55) ... Number of inner angular points in numerical integration in L609:
C          -1 ... 0 (no inner sphere)
C           0 ... 302
C           N ... N point Lebedev grid (see AngQad).
C
C2IOp(56)
C     IOp(56) ... Whether to read in density matrix from input stream in L608.
C           0 ... No.
C           1 ... Yes.
C
C2IOp(57)
C     IOp(57) ... Whether to generate data over a grid using the total SCF density:
C           0 ... No.
C           1 ... Yes, read in name for output file.
C           2 ... Yes, also read in name for input file with a different grid and compare.
C           3 ... Output in the form of data statements.
C           4 ... Fit atomic density to Gaussians.
C           5 ... Fit atomic density to Gaussians, forcing positive definateness.
C
C2IOp(58)
C     IOp(58) ... Grid to use in generating tables of density and potential.  Must be an
C                 unpruned grid.
C           0 ... Default (99001).
C
C2IOp(59)
C     IOp(59) ... Approximations to Exc
C          -1 ... Test superposition of atomic densities using L608:
C           0 ... Do correct energies.
C           1 ... Do correct energies and 0th order approximation
C           2 ... Do correct energies and 0th-1st order approximations
C           3 ... Do correct energies and 0th-2nd order approximations
C
C2IOp(60-62)
C     IOp(60-62) ... Over-ride standard values of IRadAn, IRanWt, and IRanGd.
C
CIOp(63)
C     IOp(63) ... Suppress number of electrons test in XC quadrature in L608
C                 (for debugging with small grids):
C           0 ... Default (do test).
C           1 ... Suppress test.
C           2 ... Do test as usual.
C
C2IOp(64)
C     IOp(64) ... Natural Chemical Shielding Analysis:
C           0 ... No.
C           1 ... Yes, of isotropic value.
C           2 ... Yes, of diagonal tensor elements and isotropic value.
C           3 ... Yes, of all tensor components.
C
C2IOp(65)
C     IOp(65) ... Threshold for printing of NCS contributions
C          -1 ... Zero.
C           0 ... Default (1 pmm).
C           N ... N/1000 ppm
C
C2IOp(70)
C     IOp(70) ... Control of L610.
C
C2IOp(71)
C     IOp(71) ... XC functional in L610.
C
C2IOp(72)
C     IOp(72) ... Whether to read isotopes for hyperfine interractions and
C                 do hyperfine terms in L602:
C           0 ... Default (1).
C           1 ... Yes, if open-shell, NMR data is available, and other terms
C                 are being computed
C           2 ... No.
C           3 ... Yes, regardless of other terms.
C           4 ... Yes, reading isotopes
C
C2IOp(73)
C     IOp(73) ... Whether to save orbitals from NBO:
C           0 ... Default (No).
C           1 ... Save NBOs in place of regular MOs.
C           2 ... Save NLMOs in place of regular MOs.
C           3 ... Save NLMO occupieds and NBO virtuals.
C          10 ... Suppress re-orthogonalization.
C
C2IOp(74)
C     IOp(74) ... Whether to use Gaussian connectivity in choosing Lewis
C                 structure for NBO.
C           0 ... Default (use if present and choose is selected in NBO input).
C           1 ... Use.
C           2 ... Don't use.
C
C2IOp(75)
C     IOp(75) ... model for CM2 charges.
C
C2IOp(76)
C     IOp(76) ... Threshold for linear dependence in L607.
C           0 ... Default (1.D-6).
C           N ... 10**(-N).
C
C2IOp(77)
C     IOp(77) ... Restraint in charge fitting in L602:
C           0 ... None.
C          -1 ... 2.d-4
C           N ... N * 10^-5.
C
C2IOp(78)
C     IOp(78) ... Use MOs instead of density in AtmTab.
C           0 ... Default (2).
C           1 ... Use density.
C           2 ... Use MOs.
C
C2IOp(79)
C     IOp(79) ... Whether to calculate Hirshfeld charges.
C           0 ... Default (No).
C           1 ... Yes.
C           2 ... No.
C           3 ... Yes, do atom-atom electrostatic interactions as well.
C
C2IOp(80)
C     IOp(80) ... Whether to calculate Lowdin charges and Mayer bond orders
C           0 ... Default (No).
C           1 ... Yes.
C           2 ... No.
C
C2IOp(81)
C     IOp(81) ... Print kinetic energy of orbitals?
C           0 ... Default (yes, if doing other orbital results).
C           1 ... Yes, for the top 5 occupieds and lowest 5 virtuals.
C           2 ... No.
C           3 ... Yes, for all orbitals.
C
C2IOp(82)
C     IOp(82) ... Tensors for hyperfine spectra.
C           0 ... Default, compute if there are 100 or fewer atoms
C           1 ... Compute QEq tensors and for open-shell systems compute
C                 isotropic and anisotropic splitting tensors.
C           2 ... Do not compute tensors.
C
C2IOp(83)
C     IOp(83) ... Orbital angular momentum analysis.
C           0 ... Default (No).
C           1 ... Yes, do total angular momentum contribution to each MO.
C          10 ... Report the largest atomic d and f contributions to
C                 orbitals specified by IOp(84).
C          20 ... Report the largest transition metal atomic d and f
C                 contributions to orbitals specified by IOp(84)
C          30 ... Read a list of atoms whose d and f contributions will
C                 be analyzed.
C          90 ... Do not do atomic d and f contributions.
C         100 ... Report the population of each angular momentum on each atom.
C
C2IOp(84)
C     IOp(84) ... Orbitals to analyze for d and f contributions:
C          -1 ... All orbitals
C           0 ... Just occupied orbitals
C           N ... Occupieds plus lowest N virtuals.
C
C2IOp(85)
C     IOp(85) ... Whether to do non-equilibrium state-specific perturbation
C                 of PCM excited states (this should match IOp(73) in L914):
C           0 ... Default (Yes, if the PCM SS perturbation flag is set).
C           1 ... Yes.
C           2 ... No, use equilibrium.
C
C2IOp(86)
C     IOp(86) ... Computation of multipole moments.
C           0 ... Default (1, except for PBC and old semi-empirical)
C           1 ... Calculate with DipInt.
C           2 ... Use stored moment operators.
C?
#include "commonmol.inc"
#include "commonb.inc"
#include "commonlab.inc"
      Parameter (IRwS=514,IRwFCI=517,IRwDip=518,IRwEig=522,IRwCA=524,
     $  IRwCB=526,IRwPA=528,IRwPB=530,IRwPT=532,IRwPS=534,IODens=603,
     $  IOChrg=605,IOEDrv=619,IRwNO=637,IRwMag=695,IRwPBC=733,
     $  IOMuWt=762,MaxDen=9,MaxAtPr=50,LMax=DEFLMAX,LenInR=(2*LMax+2),
     $  MxAtHy=100,MOmega=DEFMOMEGA,NOmega=DEFNOMEGA)
      Integer MaxIOp
      Parameter (MaxIOp=DEFMAXIOP)
      Common/IOp/IOp(MaxIOp)
      Common /PhyCon/ PhyCon(30)
      Common/IO/In,IOut,IPunch
      Common/MUnit/IUnit(20)
      Integer FilNum, ChToWP
      Logical GotDen, DoPrtM, DidMOs, SCFDen, DoMul, DoBMul, SavDip,
     $  AddNuc, DoMVol, DoEKT, DoMMul, DipOK, DoHirC, DFT, DoLowd,
     $  DoHyp, AllowP(50), FMM, DoOrbL, DoAtL, PCM, CBfn, Spinor,
     $  DoBior, AsyMul, ROpen, SEHam
      Integer NSymOp,NReps,LblRep,IPrmut,FMFlag,FMFlg1
      Real*8 ChrTbl,SymOps,CelVec(3,3)
      Common/RepCom/NSymOp,NReps,LblRep(32),ChrTbl(10,16),SymOps(9,10),
     $  IPrmut(MaxAtm,10)
      Dimension V(*), ScaDFX(4), IndR(LenInR), Omega(MOmega,NOmega)
      Equivalence (ToAng,PhyCon(1)), (ToKG,PhyCon(2)), (ToE,PhyCon(3))
      Save ThrMV, RJunk, Pt5
      Data ThrMV/0.0d0/, RJunk/0.0d0/, Pt5/0.5d0/
 1000 Format(/,1X,70(1h*),/)
 1010 Format(12X,'Population analysis using the SCF density.')
 1020 Format(12X,'Population analysis using the MP first order',
     $           ' density.')
 1030 Format(12X,'Population analysis using the MP2 density.')
 1040 Format(12X,'Population analysis using the MP3 density.')
 1050 Format(12X,'Population analysis using the MP4 density.')
 1060 Format(12X,'Population analysis using the CI one-particle',
     $           ' density.')
 1070 Format(12X,'Population analysis using the CI density.')
 1080 Format(12X,'Population analysis using the QCI/CC density.')
 1090 Format(5X,'Population analysis using the density correct to',
     $           ' second order.')
 1095 Format(12X,'Population analysis using the ONIOM density.')
 1100 Format(5X,'Population analysis using transition density between',
     $          ' states',I3,' and',I3,'.')
 1110 Format('   Populations using transition density between',
     $          ' ground and state',I3,'.')
 1120 Format(' Population analysis using the one-particle density of',
     $       ' excited state',I3,'.')
 2003 Format(' Require',I8,' words for approx polar but only',I8,
     $       ' available.')
 2005 Format(' Invalid value',I10,' for IOp(22).')
 2006 Format(' Density number',I2,' requested but not found.',
     $       '  Density Mask=',I10,' IMeth=',I2,'.')
 2007 Format(' Copying SCF densities to generalized density rwf,',
     $  ' IOpCl=',I2,' IROHF=',I1,'.')
 2008 Format(' Using densities on the checkpoint file.')
 2009 Format(' Unrecognized density number',I2,'.')
 2010 Format(/,' Bonding Mulliken population analysis:')
 2011 Format(' MulPop does not work for GHF yet.')
C
      Call Drum(V,NGot)
      Call ILSW(2,40,ISpCut)
      If(NBasis.eq.0.or.ISpCut.ne.0) goto 999
      Call ILSW(2,1,IOpCl)
      Call ILSW(2,22,IROHF)
      Call DecoSC(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlX,Spinor)
      ROpen = IOpCl.eq.0.and.(IROHF.eq.1.or.
     $  (IROHF.ge.2.and.Multip.gt.1.and.ITqry(IRwPS).gt.0))
      If(ROpen) NSpBlk = 2
      Call GetPBC(NPDir,CelVec)
      Call GetNB6(NBas6D)
      IPrint = IOp(33)
      Call SetPFl(IOut,IPrint,0,IPFlag,AllowP,FMM,FMFlag,FMFlg1,NFxFlg,
     $  IHMeth,Omega,0,Junk1,Junk2,V,NGot)
      Call SetInR(1,0,0,MinL,MaxL,LenRot,IndR)
      IndRI = 1
      jStrt = IndRI + LenRot
      jIPBC = jStrt
      Call RdIPBC(1,NPDir,IRwPBC,jIPBC,jStrt,V,NGot)
      IV = jStrt
      MDV1 = NGot - IV + 1
      NMtPBC = 1
      Call ObtCel(2,V(jIPBC),2,NMtPBC)
      IfMO = Mod(IOp(7),10)
      DoBior = (Mod(IOp(7),100)/10).eq.1.and.NSpBlk.eq.2
      If(IfMO.eq.0) then
        If(NPDir.gt.0) then
          IfMO = 2
        else
          IfMO = 1
          endIf
        endIf
      IfDens = IOp(8)
      If(IfDens.eq.0) then
        If(NPDir.gt.0) then
          IfDens = 2
        else
          IfDens = 1
          endIf
        endIf
      IfPop = IOp(9)
      If(IfPop.eq.0) IfPop = 1
      IfGros = IOp(10)
      If(IfGros.eq.0) IfGros = 1
      IfCond = IOp(12)
      If(IfCond.eq.0) then
        If(IPSave(0).le.1.and.NAtoms.le.MaxAtPr) then
          IfCond = 3
        else
          IfCond = 1
          endIf
        endIf
      IUpDip = IOp(18)
      IMethR = -1
      If(IOp(28).ne.0) IMethR = 0
      Call GetPCM(IOut,IPrint,.False.,0,1,ISCRF,RJunk)
      If(PCM(ISCRF)) Call PCMIOp('load',' ',IJunk,RJunk)
      Call RWLabl(2,0)
      IDX = IOp(22)
      ICtDip = IOp(86)
C
C     Copy the HF densities from the usual HF rwfs.
C
      NBD = NBasis*NDimBl
      NTTD = NRI*(NBD*(NBD+1))/2
      If(IAbs(IDX).ge.10) then
        Write(IOut,2008)
        IUnitD = IUnit(9)
        IODen1 = FilNum(IODens,IUnitD)
      else
        IUnitD = 1
        If(IPSave(0).eq.0) Write(IOut,2007) IOpCl, IROHF
        LenD = NTTD*NMtPBC
        LenDen = NSpBlk*LenD
        IV2 = IV + LenDen
        MDV2 = NGot - IV2 + 1
        Call TstCor(IV2,NGot,'MulDrv-DenSto')
        Call FileIO(2,-IRwPA,LenD,V(IV),0)
        If(NSpBlk.eq.2) Call FileIO(2,-IRwPB,LenD,V(IV+LenD),0)
        IODen1 = IODens
        Call DenSto(IOut,IODen1,0,IMethR,LenDen,V(IV),V(IV2),MDV2)
        endIf
        write(iout,*) 'End DenSto'
      If(IOpCl.ge.4) then
        Call PrtMOS(IOut,IPrint,IOpCl,IfMO,IfDens,NBasis,NBsUse,NAE,NBE,
     $    IRwEig,IRwCA,0,IRwCB,0,IRwPA,0,IRwPB,0,V(IV),MDV1)
        SavDip = IUpDip.eq.0
        Call DQG(IOut,IPrint,ICtDip,IOpCl,SavDip,IMethR,IODens,NAtoms,
     $    NBasis,NBas6D,NMtPBC,IAn,IAtTyp,AtmChg,C,IPFlag,AllowP,FMM,
     $    FMFlag,FMFlg1,NFxFlg,IHMeth,IRwPA,IRwPB,IOEDrv,ToE,ToAng,
     $    V(jIPBC),CelVec,ISCRF,V(IV),MDV1)
c        Write(IOut,2011)
C ================================================================
CJEP
C
C        Calculate Lowdin spin-charges 
C
C      Allocate
C

       NTT   = (NBasis*(NBasis+1))/2
       NBGsq = 2*(2*NBasis)**2
       NBssq = NBasis**2
       NTTBg = 8*NTT
       jstart  = IV
       jPtmp   = jstart
       jP      = jPtmp + Max(NBGsq,NTTBg)
       jXL     = jP    + Max(NBGsq,NTTBg)
       jSSq    = jXL   + Max(NBGsq,NTTBg)
       jPS     = jSSq  + NBssq
       jTmp1   = jPS   + NBssq
       jTmp2   = jTmp1 + NBssq
       jDLwdX  = jTmp2 + NBssq
       jDLwdY  = jDLwdX+ NBssq
       jDLwdZ  = jDLwdY+ NBssq
       jMap    = jDLwdZ+ NBssq
       jChgwX = jMap  + IntoWP(NBasis)
       jChgwY = jChgwX + NAtoms
       jChgwZ = jChgwY + NAtoms
       jEnd    = jChgwZ + NAtoms
       jlenx   = jend - jstart
C       Call TstCor(jEnd,MDV,'Before GLowdin')
       Call TstCor(jEnd,MDV1,'Before GLowdin')
       Write(IOut,*)'***************************************'
       Write(IOut,*)'           GHF Population'
       Write(IOut,*)'***************************************'
       Call GLowdi(IOut,V(jPtmp),V(jP),V(jXL),V(jSSq),V(jPS),V(jTmp1),
     $   V(jTmp2),NTT,IrWPA,
     $  IOPCl,V(jDLwdX),V(jDLwdY),V(jDLwdZ),V(jMap),V(jChgwX),
     $  V(jChgwY),V(jChgwZ),V(jend),jlenx)
       Write(IOut,*)'***************************************'
       Call GMuli(IOut,V(jPtmp),V(jP),V(jXL),V(jSSq),V(jPS),V(jTmp1),
     $   V(jTmp2),NTT,IrWPA,
     $  IOPCl,V(jDLwdX),V(jDLwdY),V(jDLwdZ),V(jMap),V(jChgwX),
     $  V(jChgwY),V(jChgwZ),V(jend),jlenx)
       Write(IOut,*)'***************************************'
C------------------------------------------------------------
c        Write(IOut,2011)
C ================================================================
C        Goto 998
        Goto 999
        endIf

      write(IOut, *) 'End Endif'
      Call ConvCF(.False.,.True.,.True.,.False.,.False.,IOpCl,NAE,NBE,
     $  NBasis,NBsUse,V(IV),MDV1)
      IRotat = IOp(19)
      ISelfE = IOp(25)
      DoMul = Mod(IOp(27),10).eq.0.or.Mod(IOp(27),10).eq.2
      DoBMul = (Mod(IOp(27),100)/10).eq.2
      DoMMul = (Mod(IOp(27),1000)/100).eq.1
      AsyMul = (Mod(IOp(27),10000)/1000).eq.1
      IState = IOp(29)
      JState = IOp(30)
      IfNO = IOp(31)
      IDump = IOp(34)
      DoMVol = IOp(44).eq.1
      DoMVol = .False.
      NEKT = IOp(51)
      If(NEKT.eq.0) NEKT = -2
      DoEKT = NEKT.ne.-2
      ITabGd = IOp(57)
      IRANTab = IOp(58)
      ICM2Md = IOp(75)
      IMO2Rh = IOp(78)
      IHirsh = IOp(79)
      DoHirC = Mod(IHirsh,2).eq.1
      DoLowd = IOp(80).eq.1
      IfOKE = IOp(81)
      If(IfOKE.eq.0) IfOKE = IfMO
      IDoHyp = IOp(82)
      If(IDoHyp.eq.0.and.NAtoms.le.MxAtHy.and..not.SEHam(IHMeth))
     $  IDoHyp = 1
      DoHyp = IDoHyp.eq.1.and.NMtPBC.le.1
      Call RetDFT(IOp(60),IOp(61),IOp(62),DFT,IExCor,IXCFnc,IExCoW,
     $  IRadAn,ScaHFX,ScaDFX,IRanWt,IRanGd,ICorTp)
      If(IRadAn.eq.0) IRadAn = 4
      IOrbL = IOp(83)
      DoAtL = (Mod(IOrbL,1000)/100).eq.1
      DoOrbL = Mod(IOp(83),10).eq.1
      IDoDFO = Mod(IOp(83),100)/10
      If(IDoDFO.eq.9) IDoDFO = 0
      NDFVir = IOp(84)
      IOpE = IOp(85)
C
C     Comment applies to the PBC stuff above (by Konstantin Kudin)
C     Currently l601 is able to compute PBC mulliken populations
C     total charge, and the dipole moment
C     Many other things can be extended to work for PBC as well,
C     but require bunch of trivial changes in mulpop.F (square vs. LTT, etc)
C     Therefore at this stage such changes are omitted
C
      NTT   = (NBasis*(NBasis+1))/2
      NTT6D = (NBas6D*(NBas6D+1))/2
      NBSq  = NBasis*NBasis
      ISCF = IOpCl
      If(ISCF.eq.6) ISCF = 3
      ISCF1 = ISCF
      If(ROpen) ISCF1 = 1
      If(ISCF1.eq.0) then
        LenDen = NTT*NMtPBC
      else if(ISCF1.eq.1.or.ISCF1.eq.2) then
        LenDen = 2*NTT*NMtPBC
      else
        LenDen = 4*NTT*NMtPBC
        endIf
      LRwDip = ITqry(IRwDip)
      DipOK = NTT.gt.0.and.LRwDip.gt.0
      If(DipOK) DipOK = Mod(LRwDip,NTT).eq.0
C
C     Allocate memory.
C
C     First small scratch arrays.
C
      NAtomA = NAtAct(1,NAtoms,IAtTyp,IAn)
      NDimA = Max(NAtomA,NBas6D)
      NDim = Max(NAtoms,NBas6D)
      NDimASq = Max(NDimA**2,NDim)
      IMapBf = jStrt
      IMapL = IMapBf + NBasis
      IMapAt = IMapL + NBasis
      IMapA1 = IMapAt + InToWP(NAtomA)
      IAnMap = IMapA1 + InToWP(NAtoms)
      IATMap = IAnMap + InToWP(NAtomA)
      IACMap = IATMap + InToWP(NAtomA)
      ICMap  = IACMap + NAtomA
      IMaper = ICMap  + 3*NAtomA
      IMpDgn = IMaper + InToWP(NAtoms+1)
      IGOC   = IMpDgn + InToWP(2*NBasis)
      ILSymA = IGOC   + Max(8,NSymOp)*NDim
      ILSymB = ILSymA + ChToWP(10*NBasis)
      IScr1  = ILSymB + ChToWP(10*NBasis)
      IScr2  = IScr1  + 2*NDim
      IEnd  = IScr2  + 3*NDim - 1
      Call TstCor(IEnd,NGot,'MulDrv-1')
      Call MapAct(1,.False.,.True.,.True.,NAtoms,NAtomA,IAtTyp,IAn,
     $  AtmChg,C,NAtQM,V(IMapAt),V(IMapA1),V(IAnMap),V(IATMap),
     $  V(IACMap),V(ICMap))
      Call MakMpC(.False.,.True.,1,NBasis,I5DB1,I7FB1,NShell,ShellT,
     $  ShellC,JAn,MaxAtB,MaxLB,V(IMapBf),V(IMapL))
      If(MaxLB.lt.2) IDoDFO = 0
      IPopL = IEnd + 1
      ILabPL = IPopL + NAtomA*(MaxLB+1)*4
      ILabAt = ILabPL + InToWP(3*(MaxLB+1)) - 1
      IAsMuW = ILabAt + InToWP(3*Max(NAtomA,NBsUse))
   
      write(IOut, *) 'End MakMpc'
C
C     Possibly load the array for asymmetric Mulliken analysis.
C
      NAtTT = (NAtoms*(NAtoms+1))/2
      If(AsyMul) then
        Call TstCor(IAsMuW+NAtTT,NGot,'MulDrv-AMulLd')
C        Call AMulLd(IOut,IPrint,IOMuWt,NAtoms,IAn,V(IAsMuW))
        LAsMuW = NAtTT
      else
        LAsMuW = 0
        endIf
      write(IOut, *) 'End AMulLd'
C
C     Density matrices.
C
      IPA = IAsMuW + LAsMuW
      IPB = IPA
      If(Mod(ISCF1,2).eq.1) IPB = IPA + NTT*NMtPBC
      IPAI = IPB
      If(ISCF1.ge.2) IPAI = IPAI + NTT*NMtPBC
      IPBI = IPAI
      If(ISCF1.eq.3) IPBI = IPBI + NTT*NMtPBC
C
C     Two N**2 scratch arrays for MulPop and DQ.  These should
C     immediately follow the density matrices so that GetDen can
C     use them for scratch for excited state densities.
C
      IA = IPBI + NTT6D*NMtPBC
      IB = IA + Max(NDimASq,NTT6D*NMtPBC)
      If(IfNO.eq.0) then
        LenB = 1
      else
        LenB = 2
        endIf
      IV = IB + LenB*Max(NDimASq,NTT6D*NMtPBC)
      Call TstCor(IV,NGot,'MulDrv-2')
      MDV = NGot - IV + 1
C
C     Determine which densities to analyze.
C
      If(IDX.ge.10) then
        IDOpt = IDX - 10
      else if(IDX.le.-10) then
        IDOpt = IDX + 10
      else
        IDOpt = IDX
        If(IDOpt.eq.0.and.JState.ne.0) then
          IDOpt = -4
        else if(IDOpt.eq.0.and.IState.ne.0) then
          IDOpt = -3
          endIf
        endIf
C
C     Figure out which densities to process.
C
      Call FileIO(2,-IODen1,1,XMask,0)
      Call BNUPak(XMask,IMask,IMeth)
      If(IDOpt.le.-5) then
        Write(IOut,2005) IDX
        Call Lnk1E(0)
      else if(IDOpt.eq.-4) then
        IDSt = -2
        IDEnd = -2
      else if(IDOpt.eq.-3) then
        IDSt = -1
        IDEnd = -1
      else if(IDOpt.eq.-2) then
        IDSt = 1
        IDEnd = MaxDen + 1
      else if(IDOpt.eq.-1) then
        IDSt = IMeth + 1
        If(IMeth.eq.-1) IDSt = 1
        IDEnd = IDSt
      else if(IDOpt.le.MaxDen) then
        IDSt = IDOpt + 1
        IDEnd = IDSt
      else
        Write(IOut,2005) IDX
        Call Lnk1E(0)
        endIf
      AddNuc = IDOpt.ge.-3
C
C     Do population analysis and compute multipole moments for all
C     requested densities.
C
      write(IOut, 6666) IDOpt, IDSt, IDEnd
6666  Format(I10,I10,I10)
      DidMOs = .False.
      Do 100 IDens1 = IDSt, IDEnd
        DoPrtM = (IfNO.ne.0).or..not.DidMOs
        IDens = IDens1 - 1
        Call GetDen(IOut,IUnitD,ISCF1,IDens,IState,JState,NBasis,NMtPBC,
     $    GotDen,V(IPA))
        If(GotDen) then
          DidMOs = .True.
          Write(IOut,1000)
          If(IDens.eq.-3.and.JState.ne.0) then
            Write(IOut,1100) IState, JState
          else if(IDens.eq.-3) then
            Write(IOut,1110) IState
          else if(IDens.eq.-2) then
            Write(IOut,1120) IState
          else if(IDens.eq.0) then
            Write(IOut,1010)
          else if(IDens.eq.1) then
            Write(IOut,1020)
          else if(IDens.eq.2) then
            Write(IOut,1030)
          else if(IDens.eq.3) then
            Write(IOut,1040)
          else if(IDens.eq.4) then
            Write(IOut,1050)
          else if(IDens.eq.5) then
            Write(IOut,1060)
          else if(IDens.eq.6) then
            Write(IOut,1070)
          else if(IDens.eq.7) then
            Write(IOut,1080)
          else if(IDens.eq.8) then
            Write(IOut,1090)
          else if(IDens.eq.9) then
            Write(IOut,1095)
          else
            Write(IOut,2009) IDens
            Call Lnk1E(0)
            endIf
          Write(IOut,1000)
          If(ISCF1.ge.2) then
            Call AMove(LenDen,V(IPA),V(IA))
            Call VC2RI(NTT,V(IA),V(IPA),V(IPAI))
            If(ISCF1.eq.3) Call VC2RI(NTT,V(IA+2*NTT),V(IPB),V(IPBI))
            endIf
          If(DoMul) then
            If(DoPrtM) Call MulMOs(IOut,IPrint,IfMO,IfNO,DoBior,IROHF,
     $        ROpen,IRwS,IRwEig,IRwCA,IRwCB,IRwNO,NAtoms,Multip,NAE,NBE,
     $        NBasis,NBas6D,NBsUse,NDim,NDimA,V(IPA),V(IPB),V(IA),
     $        MDV-IA+1,V(IScr1),V(IB),V(IScr2),V(IMaper),V(IMpDgn),
     $        V(IGOC),V(ILSymA),V(ILSymB),NMtPBC,MinL,MaxL,IndR,LenRot,
     $        V(IndRI))
            Call MulPop(IOut,' ',IfDens,IfPop,IfGros,IfCond,DoAtL,IROHF,
     $        IRwS,IRwPS,IOChrg,NAtoms,NAtomA,Multip,NBasis,NDim,NDimA,
     $        IAn,C,AddNuc,V(IPA),V(IPB),V(IPAI),V(IPBI),V(IA),V(IB),
     $        V(IGOC),V(IMapBf),NMtPBC,V(IMapAt),V(IMapA1),V(IAnMap),
     $        V(IACMap),AtmChg,MaxLB,V(IMapL),V(IPopL),V(ILabPL),
     $        V(ILabAt),AsyMul,V(IAsMuW))
            endIf
          If(DoBMul) then
            Write(IOut,2010)
            Call MakMpA(NBasis,V(IScr1))
            Call BlankP(NBasis,V(IScr1),V(IPA))
            If(Mod(ISCF1,2).eq.1)
     $        Call BlankP(NBasis,V(IScr1),V(IPB))
            If(ISCF1.ge.2) Call BlankP(NBasis,V(IScr1),V(IPAI))
            If(ISCF1.ge.3) Call BlankP(NBasis,V(IScr1),V(IPBI))
            Call MulPop(IOut,' Bonding ',IfDens,IfPop,IfGros,IfCond,
     $        DoAtL,IROHF,IRwS,IRwPS,IOChrg,NAtoms,NAtomA,Multip,NBasis,
     $        NDim,NDimA,IAn,C,.False.,V(IPA),V(IPB),V(IPAI),V(IPBI),
     $        V(IA),V(IB),V(IGOC),V(IMapBf),NMtPBC,V(IMapAt),V(IMapA1),
     $        V(IAnMap),V(IACMap),AtmChg,MaxLB,V(IMapL),V(IPopL),
     $        V(ILabPL),V(ILabAt),AsyMul,V(IAsMuW))
            Call GetDen(IOut,IUnitD,ISCF1,IDens,IState,JState,NBasis,
     $        NMtPBC,GotDen,V(IPA))
            If(.not.GotDen) Call Lnk1E(0)
            endIf
C
C         Possibly compute APT charges.
C
          jQ = IV
          jDDip = jQ + NAtoms
          jEnd  = jDDip + 9*NAtoms - 1
          Call TstCor(jEnd,NGot,'MulDrv-GAPT')
          Call GAPT(IOut,IOEDrv,IUnitD,NAtoms,IAn,MicOpt,C,V(jQ),
     $      V(jDDip))
C
C         Possibly compute Lowdin and class IV charges.
C
          If(ICM2Md.gt.0.or.DoLowd) then
            jMap   = IV
            jChgL  = jMap   + InToWP(NBasis)
            jChgC4 = jChgL  + NAtoms
            jDLwd  = jChgC4 + NAtoms
            jS     = jDLwd  + NTT
            jASq   = jS     + NBSq
            jBSq   = jASq   + NBSq
            jA     = jBSq   + NBSq
            jB     = jA     + NBasis
            jV     = jB     + NBasis
            Call TstCor(jV,NGot,'MulDrv-Class4')
            MDV1 = NGot - jV
            Call Class4(IOut,NAtoms,NBasis,ICM2Md,ISCF1,IAn,AtmChg,C,
     $        V(IPA),V(IPB),V(jMap),V(jChgL),V(jChgC4),V(jDLwd),V(jS),
     $        V(jASq),V(jBSq),V(jA),V(jB),V(jV),MDV1)
            endIf
C
C         Compute multipole moments.
C
          If(Mod(ISCF1,2).eq.1)
     $      Call AAdd(NTT*NMtPBC,V(IPA),V(IPB),V(IPA))
          If(ITqry(IRwDip).ge.1) then
            Len2 = MDV / (NTT*NMtPBC)
            SavDip = IMeth.ge.0.and.IDens.eq.IMeth.and.IUpDip.eq.0
            Call AMove(NTT*NMtPBC,V(IPA),V(IA))
            Call DQ(IOut,IPrint,ICtDip,IPFlag,AllowP,FMM,FMFlag,FMFlg1,
     $        NFxFlg,IHMeth,SavDip,AddNuc,IRwDip,IOEDrv,NAtoms,NBasis,
     $        NTT,Len2,IAn,IAtTyp,AtmChg,C,V(IA),V(IV),ToE,ToAng,
     $        V(jIPBC),CelVec,ISCRF,V(IV),MDV)
            endIf
C
C         Compute electrostatic energies.
C
          SCFDen = IDens.eq.0
          Call ElEner(IOut,ISelfE,SCFDen,ISCF,IROHF,ROpen,NAE,NBE,
     $      NBasis,NBsUse,V(IPA),V(IV),MDV)
          If(Mod(ISCF1,2).eq.1)
     $      Call ASub(NTT*NMtPBC,V(IPA),V(IPB),V(IPA))
          If(DoHirC) then
            Call AMove(NTT,V(IPA),V(IA))
            If(Mod(ISCF1,2).eq.0) then
              Call AScale(NTT,Pt5,V(IPA),V(IA))
              Call AMove(NTT,V(IA),V(IA+NTT6D))
            else
              Call AMove(NTT,V(IPB),V(IA+NTT6D))
              endIf
            Call HirChg(IOut,IPrint,ISCF1,IHirsh,IRadAn,IRanWt,IRanGd,
     $        NAtoms,NBasis,IAn,IAtTyp,AtmChg,C,V(IA),V(IA+NTT6D),V(IV),
     $        MDV)
            endIf
        else if(IDens1.gt.0.and.IDens1.eq.IDSt.and.IDens1.eq.IDEnd) then
          Write(IOut,2006) IDens, IMask, IMeth
          Call Lnk1E(0)
          endIf
        If(DoEKT.and.IDens.ne.0) Call EKT(IOut,IPrint,IPFlag,AllowP,
     $    IUnitD,IDens,NEKT,NAE,NBE,NBasis,NBas6D,V(IV),MDV)
  100   Continue
C
C     Do Minimal population analysis.
C
      If(DoMMul) Call MulMin(IOut,IPrint,ISCF,IfDens,IfPop,IfGros,
     $  IfCond,DoAtL,IRwPS,IOChrg,NMtPBC,NAtoms,NAtomA,Multip,NBasis,
     $  NAE,NBE,NDim,NDimA,IRwCA,IRwCB,IAn,IAtTyp,AtmChg,C,V(IPA),
     $  V(IPB),V(IA),V(IMapBf),V(IGOC),V(IMapAt),V(IMapA1),V(IAnMap),
     $  V(IACMap),MaxLB,V(IMapL),V(IPopL),V(ILabPL),V(ILabAt),AsyMul,
     $  V(IAsMuW),V(IV),MDV)
C
C     Do orbital angular momentum contribution to each orbital.
C
      IV = IPA
      MDV1 = MDV - IV + 1
      If(DoOrbL.or.IDoDFO.gt.0) Call MulOrL(IOut,IPrint,IDump,ISCF,
     $  DoOrbL,IDoDFO,NDFVir,NAtoms,NAtomA,NBasis,NBsUse,NAE,NBE,MaxLB,
     $  IRwCA,IRwCB,IRwS,IAn,V(IAnMap),V(IATMap),V(IMapBf),V(IMapAt),
     $  V(IMapA1),V(IMapL),V(ILabAt),V(ILabPL),V(IV),MDV1)
C
C     Compute orbital kinetic energies.
C
      Call OKE(IOut,IPrint,ISCF,NPDir,IfOKE,NBasis,NBsUse,NAE,NBE,V,MDV)
C
C     Compare approximate and exact polarizabilities, if exact ones
C     are available.
C
      Call FileIO(2,-IOEDrv,1,XMask,0)
      Call BNUPak(XMask,I1,I2)
      If(IROHF.eq.0.and.Mod(I1,8).ge.4.and.DipOK) then
        IExPol = 1
        ICMO = IExPol + 6
        IT = ICMO + NBSq
        IEV = IT + NBSq
        IDip = IEV + NBSq
        IEnd = IDip + 3*NTT*NMtPBC - 1
        If(IEnd.le.NGot) then
          Ind = 1 + 3 + 9*NAtoms
          Call FileIO(2,-IOEDrv,6,V(IExPol),Ind)
          Call PrtPol(IOut,IOpCl,IROHF,IRotat,IRwDip,NAE,NBE,NBasis,
     $      NBsUse,NTT,V(IExPol),V(ICMO),V(IT),V(IEV),V(IDip))
        else
          Write(IOut,2003) IEnd, NGot
          endIf
        endIf
      If(DoMVol) Call CalcMV(IOut,IPrint,IDump,ISCF,0,NAtoms,NBasis,IAn,
     $  IAtTyp,AtmChg,C,ThrMV,V,NGot)
      If(ITabGd.ge.4) then
        Call AtmFit(IOut,IPrint,IPFlag,AllowP,FMM,FMFlag,FMFlg1,NFxFlg,
     $    IHMeth,ISCF,ITabGd,NAtoms,NAE,NBE,NBasis,NBas6D,NBsUse,IAn,
     $    IAtTyp,AtmChg,C,V,NGot)
      else if(ITabGd.ne.0) then
        Call AtmTab(In,IOut,IPrint,IPFlag,AllowP,ISCF,ITabGd,IRANTab,
     $    IMO2Rh,NAtoms,NAE,NBE,NBasis,IAn,IAtTyp,AtmChg,C,V,NGot)
        endIf
C
C     Generate various tensors related to hyperfine spectra.
C
      If(DoHyp) Call HypTen(IOut,IPrint,IPFlag,AllowP,ISCF,ISCF1,IROHF,
     $  IUnitD,IDens,IRwPS,IRwPT,IRwFCI,IRwMag,NAtoms,NBasis,NBas6D,
     $  Multip,IAn,IAtTyp,AtmChg,C,IAtSpn,IAtWgt,AtGFac,AtmWgt,V,MDV)
C
C     Do the state specific perturbation for PCM excited states.
C
      If(PCM(ISCRF)) then
        If(IOpE.eq.0) then
          IEInf = 1
        else if(IOpE.eq.2) then
          IEInf = 0
        else
          IEInf = IOpE
          endIf
        Call PCMSSP(IOut,IPrint,IEInf,IOpCl,IPFlag,AllowP,FMM,FMFlag,
     $    FMFlg1,NFxFlg,IHMeth,NBasis,NBas6D,NAtoms,IAn,IAtTyp,AtmChg,C,
     $    V,MDV)
        endIf
      Call ConvCF(.True.,.True.,.True.,.False.,.False.,IOpCl,NAE,NBE,
     $  NBasis,NBsUse,V,NGot)
  998 If(PCM(ISCRF)) Call PCMIOp('save',' ',IJunk,RJunk)
  999 Call ChainX(0)
      Return
      End
