#import phf
from pyscf import gto, scf
#import scipy
#import numpy as np
import sys

from pyphf import util

xyz = sys.argv[1]
fch = sys.argv[2]
bas = xyz[:-3] + 'bas'


mf = util.guess_from_fchk(xyz, bas, fch)

mf2 = util.SUHF(mf)
mf2.cut_no = False
mf2.debug = False
mf2.max_cycle = 15
mf2.kernel()
