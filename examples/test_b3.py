#import phf
from fch2py import fch2py
from ortho import check_orthonormal
from pyscf import gto, scf, lib
#import scipy
#import numpy as np
import sys

from pyphf import util

lib.num_threads(4)

#xyz = 'eth.xyz' #sys.argv[1]
#fch = 'eth_uhf.fch' #sys.argv[2]
#bas = 'def2svp'


#mf = util.guess_from_fchk(xyz, bas, fch)


mol = gto.M()
#   3 atom(s)
mol.atom = '''
B1             -0.17964100        2.61877500        0.00000000
B2             -0.94741200        1.28916200        0.00000000
B3              0.58813000        1.28916200        0.00000000
'''

mol.basis = {
'B1': gto.basis.parse('''
B     S 
    0.457000000E+04    0.695554418E-03
    0.685900000E+03    0.534957299E-02
    0.156500000E+03    0.271166287E-01
    0.444700000E+02    0.101315096E+00
    0.144800000E+02    0.271880829E+00
    0.513100000E+01    0.448115931E+00
    0.189800000E+01    0.289937262E+00
    0.332900000E+00    0.143128310E-01
B     S 
    0.457000000E+04   -0.281280224E-03
    0.685900000E+03   -0.221988781E-02
    0.156500000E+03   -0.110164715E-01
    0.444700000E+02   -0.443491899E-01
    0.144800000E+02   -0.120912048E+00
    0.513100000E+01   -0.280737900E+00
    0.189800000E+01   -0.266066809E+00
    0.332900000E+00    0.109178413E+01
B     S 
    0.104300000E+00    0.100000000E+01
B     P 
    0.600100000E+01    0.541658633E-01
    0.124100000E+01    0.302379890E+00
    0.336400000E+00    0.771292217E+00
B     P 
    0.953800000E-01    0.100000000E+01
B     D 
    0.343000000E+00    0.100000000E+01
'''),
'B2': gto.basis.parse('''
B     S 
    0.457000000E+04    0.695554418E-03
    0.685900000E+03    0.534957299E-02
    0.156500000E+03    0.271166287E-01
    0.444700000E+02    0.101315096E+00
    0.144800000E+02    0.271880829E+00
    0.513100000E+01    0.448115931E+00
    0.189800000E+01    0.289937262E+00
    0.332900000E+00    0.143128310E-01
B     S 
    0.457000000E+04   -0.281280224E-03
    0.685900000E+03   -0.221988781E-02
    0.156500000E+03   -0.110164715E-01
    0.444700000E+02   -0.443491899E-01
    0.144800000E+02   -0.120912048E+00
    0.513100000E+01   -0.280737900E+00
    0.189800000E+01   -0.266066809E+00
    0.332900000E+00    0.109178413E+01
B     S 
    0.104300000E+00    0.100000000E+01
B     P 
    0.600100000E+01    0.541658633E-01
    0.124100000E+01    0.302379890E+00
    0.336400000E+00    0.771292217E+00
B     P 
    0.953800000E-01    0.100000000E+01
B     D 
    0.343000000E+00    0.100000000E+01
'''),
'B3': gto.basis.parse('''
B     S 
    0.457000000E+04    0.695554418E-03
    0.685900000E+03    0.534957299E-02
    0.156500000E+03    0.271166287E-01
    0.444700000E+02    0.101315096E+00
    0.144800000E+02    0.271880829E+00
    0.513100000E+01    0.448115931E+00
    0.189800000E+01    0.289937262E+00
    0.332900000E+00    0.143128310E-01
B     S 
    0.457000000E+04   -0.281280224E-03
    0.685900000E+03   -0.221988781E-02
    0.156500000E+03   -0.110164715E-01
    0.444700000E+02   -0.443491899E-01
    0.144800000E+02   -0.120912048E+00
    0.513100000E+01   -0.280737900E+00
    0.189800000E+01   -0.266066809E+00
    0.332900000E+00    0.109178413E+01
B     S 
    0.104300000E+00    0.100000000E+01
B     P 
    0.600100000E+01    0.541658633E-01
    0.124100000E+01    0.302379890E+00
    0.336400000E+00    0.771292217E+00
B     P 
    0.953800000E-01    0.100000000E+01
B     D 
    0.343000000E+00    0.100000000E+01
''')}

# Remember to check the charge and spin
mol.charge = 0
mol.spin = 1
mol.verbose = 4
mol.build()

mf = scf.UHF(mol)
mf.max_cycle = 1
mf.kernel()

# read MOs from .fch(k) file
nbf = mf.mo_coeff[0].shape[0]
nif = mf.mo_coeff[0].shape[1]
S = mol.intor_symmetric('int1e_ovlp')
Sdiag = S.diagonal()
alpha_coeff = fch2py('B3/B3_cc-pVDZ_5D7F_uhf.fch', nbf, nif, Sdiag, 'a')
beta_coeff = fch2py('B3/B3_cc-pVDZ_5D7F_uhf.fch', nbf, nif, Sdiag, 'b')
mf.mo_coeff = (alpha_coeff, beta_coeff)
# read done

# check if input MOs are orthonormal
check_orthonormal(nbf, nif, mf.mo_coeff[0], S)
check_orthonormal(nbf, nif, mf.mo_coeff[1], S)

dm = mf.make_rdm1()
mf.max_cycle = 10
mf.kernel(dm)


mf2 = util.SUHF(mf)
#mf2.cut_no = False
#mf2.debug = False
mf2.diis_on = True
mf2.diis_start_cyc = 20
mf2.max_cycle = 100
mf2.tofch = True
mf2.oldfch = 'B3/B3_cc-pVDZ_5D7F_uhf.fch'
mf2.kernel()
