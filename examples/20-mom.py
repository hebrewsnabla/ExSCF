from pyphf import suscf, guess
from pyscf import lib

lib.num_threads(4)
xyz = ''' O 0.0 0.0 0.0; H 0.0 -0.757 0.587; H 0.0 0.757 0.587'''
bas = '6-31gs'
mf = guess.gen(xyz, bas, 0, 0)

mf2 = suscf.SUHF(mf)
#mf2.verbose = 8
#mf2.printmo = True
mf2.kernel()
mo = mf2.mo_reg

mf3 = suscf.SUHF(mf)
mf3.setmom =  [[4],[5]], [[],[]]
mf3.mom_reforb = mo
#mf3.printmo = True
#mf3.verbose = 8
#mf3.diis_on = False
mf3.diis_start_cyc = 10
mf3.kernel()


