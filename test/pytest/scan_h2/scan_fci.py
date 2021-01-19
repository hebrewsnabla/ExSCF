#import phf
from pyscf import lib, fci
#import scipy
import numpy as np
import sys, os, contextlib

from pyphf import util, sudft

#xyz = 'h2.xyz' #sys.argv[1]
fch = 'test_uhf.fchk' #sys.argv[2]
#bas = '../h2.bas'
bas = 'cc-pvtz'

lib.num_threads(16)

for r in np.arange(0.5, 2.55, 0.1):
    output = 'fci/fci_%.1f.out' % r
    os.system("echo '\n' > %s" % output)
    with open(output, 'a', encoding='utf-8') as f:
        with contextlib.redirect_stdout(f):
            xyz = '''H 0.0 0.0 0.0; H 0.0 0.0 %f'''%r
            mf = util.guess_from_fchk(xyz, bas, fch)
            #mf = util.guess(xyz, bas)
            mf.max_cycle = 50
            mf.kernel()
            ci = fci.FCI(mf)
            print('E(FCI) = %15.8f' % ci.kernel()[0])