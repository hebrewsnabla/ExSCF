#import phf
from pyscf import lib
from pyphf import suscf, guess

lib.num_threads(4)

xyz = '''
 C                 -0.33039647   -0.83700439    0.00000000
 H                  0.18655089   -1.81238157    0.00000000
 H                 -1.43340448   -0.88137274    0.00000000
 O                  0.27750985    0.22075262    0.00000000
''' #sys.argv[1]
bas = 'cc-pvdz'


mf = guess.from_frag(xyz, bas, [[0,1,2],[3]], [0,0], [-2,2])

mf2 = suscf.SUHF(mf)
#mf2.cut_no = False
#mf2.verbose = 6
#mf2.diis_on = True
mf2.max_cycle = 50
mf2.kernel()
