from pyphf import util, guess

xyz = 'H 0.0 0.0 0.0; H 0.0 0.0 2.0'''
bas = '3-21g'
mf = guess.mix(xyz, bas, conv='tight')

mf2 = util.SUHF(mf)
mf2.verbose = 8
mf2.kernel()

