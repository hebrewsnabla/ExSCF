#import phf
from pyscf import gto, scf, lib
#import scipy
#import numpy as np
import sys

from pyphf import util, frag

lib.num_threads(1)

xyz = 'N 0.0 0.0 0.0; N  0.0 0.0 0.9' #sys.argv[1]
#fch = 'n2.fchk' #sys.argv[2]
bas = 'cc-pvdz'

mol = gto.Mole()
mol.atom = xyz
mol.basis = bas
mol.build()

dm, mo, occ = frag.guess_frag(mol, [[0],[1]], [0,0], [3,-3])

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
mf2.max_cycle = 100
mf2.kernel()
