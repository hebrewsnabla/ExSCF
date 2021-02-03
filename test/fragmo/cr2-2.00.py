#import phf
from pyscf import gto, scf, lib
#import scipy
#import numpy as np
import sys

from pyphf import util, frag

lib.num_threads(4)

xyz = 'Cr 0.0 0.0 0.0; Cr  0.0 0.0 2.00' #sys.argv[1]
#fch = 'n2.fchk' #sys.argv[2]
bas = 'def2-svp'

mol = gto.Mole()
mol.atom = xyz
mol.basis = bas
mol.build()

dm, mo, occ = frag.guess_frag(mol, [[0],[1]], [0,0], [6,-6])

#mf = util.guess_from_fchk(xyz, bas, fch)

mf = scf.UHF(mol)
mf.verbose = 6
#mf.conv_tol = 1e-2
mf.max_cycle = 2
mf.kernel(dm0 = dm)

mf2 = util.SUHF(mf)
mf2.cut_no = False
mf2.verbose = 4
mf2.diis_on = True
mf2.diis_start_cyc = 5
mf2.level_shift = 0.5
mf2.max_cycle = 50
mf2.kernel()
