from pyscf import lib, scf
from pyphf import guess, suscf
import autocas

lib.num_threads(8)

xyz = '''C                 -2.94294278    0.39039038    0.00000000
 C                 -1.54778278    0.39039038    0.00000000
 C                 -0.85024478    1.59814138    0.00000000
 C                 -1.54789878    2.80665038   -0.00119900
 C                 -2.94272378    2.80657238   -0.00167800
 C                 -3.64032478    1.59836638   -0.00068200
 H                 -3.49270178   -0.56192662    0.00045000
 H                 -0.99827478   -0.56212262    0.00131500
 H                  0.24943522    1.59822138    0.00063400
 H                 -0.99769878    3.75879338   -0.00125800
 H                 -3.49284578    3.75885338   -0.00263100
 H                 -4.73992878    1.59854938   -0.00086200
'''
bas = 'def2-svp'

mf = guess.mix(xyz, bas, conv='tight')

mf2 = suscf.SUHF(mf)
#mf2.level_shift = 0.3
#mf2.max_cycle=20
mf2.kernel()

