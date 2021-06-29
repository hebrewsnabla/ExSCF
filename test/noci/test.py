#import phf
from pyscf import lib
from pyphf import suscf, guess, noci

lib.num_threads(4)

xyz = '''
N -1.535079   -1.050802   -0.265895
H -1.424942   -1.420987    0.683068
H -1.172743   -1.798291   -0.863000
C -0.674001    0.116995   -0.396238
H -0.667420    0.432032   -1.454565
C -1.220146    1.281062    0.444439
C  0.781650   -0.165583   -0.015767
H -1.219679    1.017906    1.515968
H -0.611700    2.187271    0.312390
H -2.257523    1.488947    0.145218
O  1.176250   -1.168621    0.535350
O  1.608113    0.851853   -0.361143
H  2.499646    0.588022   -0.066072
''' #sys.argv[1]
#fch = '../test_carlos/test_uhf.fchk' #sys.argv[2]
bas = '3-21g'


#mf = util.guess_from_fchk(xyz, bas, fch)
mf = guess.gen(xyz, bas, 0, 0)

mf2 = suscf.SUHF(mf)
#mf2.cut_no = False
#mf2.verbose = 6
#mf2.diis_on = True
mf2.max_cycle = 30
mf2.kernel()

#mf4 = pt2.EMP2(mf2)
#mf4.kernel()
noci.test(mf2)