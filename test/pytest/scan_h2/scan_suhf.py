#import phf
from pyscf import lib
#import scipy
import numpy as np
import sys, os, contextlib

from pyphf import util, sudft

#xyz = 'h2.xyz' #sys.argv[1]
fch = 'test_uhf.fchk' #sys.argv[2]
#bas = '../h2.bas'
bas = 'cc-pvtz'

lib.num_threads(8)

for r in np.arange(0.5, 2.55, 0.1):
    output = 'suhf/suhf_%.1f_diis.out' % r
    os.system("echo '\n' > %s" % output)
    with open(output, 'w', encoding='utf-8') as f:
        with contextlib.redirect_stdout(f):
            xyz = '''H 0.0 0.0 0.0; H 0.0 0.0 %f'''%r
            mf = util.guess_from_fchk(xyz, bas, fch)
            #mf = util.guess(xyz, bas)
            #f.flush()
            mf2 = util.SUHF(mf)
            #mf2.cut_no = False
            mf2.debug = False
            #mf2.max_cycle = 50
            mf2.diis_on = True
            #mf2.output = 'sudft_%.1f.out'%r
            mf2.kernel()
            
            #mf3 = sudft.SUDFT(mf2)
            #mf3.suxc = 'TPSS'
            #mf3.kernel()
                    
