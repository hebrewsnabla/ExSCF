#import phf
from pyscf import gto, scf, lib
#import scipy
#import numpy as np
import sys

from pyphf import util

lib.num_threads(4)

xyz = 'eth.xyz' #sys.argv[1]
fch = 'eth_uhf.fch' #sys.argv[2]
bas = 'def2svp'


mf = util.guess_from_fchk(xyz, bas, fch)

mf2 = util.SUHF(mf)
#mf2.cut_no = False
#mf2.debug = False
mf2.diis_on = True
mf2.diis_start_cyc = 20
mf2.max_cycle = 100
mf2.tofch = True
mf2.oldfch = 'eth_uhf.fch'
mf2.kernel()
