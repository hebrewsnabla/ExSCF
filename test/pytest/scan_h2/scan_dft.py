#import phf
from pyscf import gto, scf
#import scipy
#import numpy as np
import sys, os, contextlib

from pyphf import util, sudft

#xyz = 'h2.xyz' #sys.argv[1]
fch = '../../test_carlos/test_uhf.fchk' #sys.argv[2]
bas = '../h2.bas'

for r in [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]:
    output = 'sudft_%.1f.out' % r
    os.system("echo '\n' > %s" % output)
    with open(output, 'a', encoding='utf-8') as f:
        with contextlib.redirect_stdout(f):
            xyz = '''H 0.0 0.0 0.0; H 0.0 0.0 %f'''%r
            mf = util.guess_from_fchk(xyz, bas, fch)
            #f.flush()
            mf2 = util.SUHF(mf)
            mf2.cut_no = False
            mf2.debug = False
            mf2.max_cycle = 50
            #mf2.output = 'sudft_%.1f.out'%r
            #mf2.kernel()
            
            mf3 = sudft.SUDFT(mf2)
            mf3.kernel()
                    