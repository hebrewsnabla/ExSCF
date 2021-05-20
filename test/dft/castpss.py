#import phf
from pyscf import lib
#import scipy
import numpy as np
import sys, os, contextlib

from pyphf import util, casdft, guess
from pyAutoMR import autocas

#xyz = 'h2.xyz' #sys.argv[1]
#fch = 'test_uhf.fchk' #sys.argv[2]
#bas = '../h2.bas'
#bas = 'cc-pvtz'

lib.num_threads(4)

#for r in np.arange(0.8, 2.55, 0.1):
r = 2.0
output = 'castpss_%.1f.out' % r
os.system("echo '\n' > %s" % output)
#print('scan %.1f' % r)
with open(output, 'a', encoding='utf-8') as f:
    with contextlib.redirect_stdout(f):
        xyz = '''N 0.0 0.0 0.0; N 0.0 0.0 %f'''%r
        #mf = util.guess_from_fchk(xyz, bas, fch)
        bas = 'cc-pvtz'
        mf = guess.from_frag(xyz, bas, [[0],[1]], [0,0], [3,-3], cycle=50)
        guess.check_stab(mf)
        #f.flush()
        #mf2 = util.SUHF(mf)
        #mf2.cut_no = False
        #mf2.verbose = 4
        #mf2.max_cycle = 50
        #mf2.diis_on = True
        #mf2.output = 'sudft_%.1f.out'%r
        #mf2.kernel()
        mc = autocas.cas(mf, crazywfn=True)
        
        mf3 = casdft.CASDFT(mc)
        #mf3.dens = 'relaxed'
        mf3.suxc = 'TPSS'
        mf3.trunc = 'fc'
        mf3.kernel()
                    
