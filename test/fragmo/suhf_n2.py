#import phf
from pyscf import gto, scf, lib
#import scipy
#import numpy as np
import sys

from pyphf import util

lib.num_threads(1)

xyz = 'N 0.0 0.0 0.0; N  0.0 0.0 0.9' #sys.argv[1]
fch = 'n2.fchk' #sys.argv[2]
bas = 'cc-pvdz'


mf = util.guess_from_fchk(xyz, bas, fch)

mf2 = util.SUHF(mf)
mf2.cut_no = False
mf2.verbose = 4
mf2.diis_on = True
mf2.diis_start_cyc = 10
mf2.max_cycle = 100
mf2.kernel()
