from pyphf import suscf, guess
from pyscf import lib

lib.num_threads(4)
xyz = ''' O 0.0 0.0 0.0; H 0.0 -0.757 0.587; H 0.0 0.757 0.587'''
bas = '6-31gs'
mf = guess.gen(xyz, bas,0,0)

mf2 = suscf.SUHF(mf)
#mf2.verbose = 8
#mf2.printmo = True
mf2.output = 'mom'
mf2.kernel()
#mo = mf2.mo_reg


mo = lib.chkfile.load('mom_su.pchk', 'scf/mo_coeff')

mf3 = suscf.SUHF()
mf3.chkfile0 = 'mom_ges.pchk'
mf3.setmom =  [[4],[5]], [[],[]]
mf3.mom_reforb = mo
mf3.mom_start_cyc = 1
#mf3.printmo = True
#mf3.verbose = 8
#mf3.diis_on = False
mf3.output = 'mom1'
mf3.kernel()


