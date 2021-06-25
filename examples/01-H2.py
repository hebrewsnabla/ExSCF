from pyphf import suscf, guess
from pyscf import lib
lib.num_threads(4)

xyz = 'H 0.0 0.0 0.0; H 0.0 0.0 2.0'''
bas = '3-21g'
mf = guess.mix(xyz, bas, conv='tight')

mf2 = suscf.SUHF(mf)
mf2.verbose = 8
mf2.kernel()

