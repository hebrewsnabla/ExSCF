
from pyphf import suscf, guess, supdft
from pyphf import pt2_new as pt2
from pyscf import lib, mp


''' EMP2(0) for 1Delta state of O2'''

lib.num_threads(8)
xyz = '''O 0.0 0.0 0.0; O 0.0 0.0 1.208'''
bas = 'cc-pvdz'
mf = guess.mix(xyz, bas, cycle=2)

mf2 = suscf.SUHF(mf)
#mf2.verbose = 8
#mf2.noiter = True
mf2.kernel()

mf5 = pt2.EMP2(mf2)
#mf5.use_det = True
#mf5.vap = True
mf5.do_sc = False
#mf5.do_biort = True
#mf5.do_15 = True
mf5.kernel()


