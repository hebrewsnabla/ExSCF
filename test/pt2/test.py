#import phf
from pyscf import lib
#import scipy
#import numpy as np
#import sys

from pyphf import suscf, guess, pt2

lib.num_threads(4)

xyz = '''H 0.0 0.0 0.0; H 0.0 0.0 0.7''' #sys.argv[1]
#fch = '../test_carlos/test_uhf.fchk' #sys.argv[2]
bas = '3-21g'


#mf = util.guess_from_fchk(xyz, bas, fch)
mf = guess.mix(xyz, bas, cycle=2)

mf2 = suscf.SUHF(mf)
mf2.cut_no = False
mf2.verbose = 6
mf2.diis_on = True
mf2.max_cycle = 30
mf2.kernel()

mf4 = pt2.EMP2(mf2)
mf4.kernel()