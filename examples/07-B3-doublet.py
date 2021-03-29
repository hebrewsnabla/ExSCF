#import phf
from fch2py import fch2py
from ortho import check_orthonormal
from pyscf import lib
#import scipy
#import numpy as np
#import sys

from pyphf import util, guess

lib.num_threads(4)

xyz = '''
B1             -0.17964100        2.61877500        0.00000000
B2             -0.94741200        1.28916200        0.00000000
B3              0.58813000        1.28916200        0.00000000
'''
mf = guess.gen(xyz, 'cc-pvdz', 0, 1)

mf2 = util.SUHF(mf)
#mf2.cut_no = False
#mf2.debug = False
mf2.diis_on = True
#mf2.diis_start_cyc = 20
mf2.max_cycle = 100
#mf2.tofch = True
#mf2.oldfch = 'B3/B3_cc-pVDZ_5D7F_uhf.fch'
mf2.kernel()

# E = -73.80340789
