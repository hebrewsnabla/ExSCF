from pyphf import suscf, guess
from pyscf import lib

lib.num_threads(4)
xyz = 'H 0.0 0.0 0.0; H 0.0 0.0 2.0'''
bas = '3-21g'
mf = guess.mix(xyz, bas, conv='tight')

mf2 = suscf.SUHF(mf)
#mf2.verbose = 8
mf2.kernel()
mo = mf2.mo_reg

mf3 = suscf.SUHF(mf)
mf3.setmom = [[],[]], [[0],[1]]
mf3.mom_reforb = mo
mf3.printmo = True
#mf3.verbose = 8
mf3.diis_on = False
mf3.kernel()


