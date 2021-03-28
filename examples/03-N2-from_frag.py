#import phf
from pyscf import lib
#import scipy
#import numpy as np
#import sys

from pyphf import util, guess

lib.num_threads(4)

xyz = 'N 0.0 0.0 0.0; N  0.0 0.0 0.9' #sys.argv[1]
#fch = 'n2.fchk' #sys.argv[2]
bas = 'cc-pvdz'

mf = guess.from_frag(xyz, bas, [[0],[1]], [0,0], [3,-3])

mf2 = util.SUHF(mf)
#mf2.cut_no = False
mf2.verbose = 4
mf2.diis_on = True
#mf2.diis_start_cyc = 5
mf2.max_cycle = 100
mf2.kernel()
