from pyscf import gto, scf
from fch2py import fch2py
from ortho import check_orthonormal

mol = gto.M()
#   2 atom(s)
mol.atom = '''
N1              0.00000000        0.00000000        0.65000000
N2              0.00000000        0.00000000       -0.65000000
'''

mol.basis = {
'N1': gto.basis.parse('''
N     S 
    0.904600000E+04    0.701708743E-03
    0.135700000E+04    0.540299880E-02
    0.309300000E+03    0.274729510E-01
    0.877300000E+02    0.103514580E+00
    0.285600000E+02    0.279586579E+00
    0.102100000E+02    0.451317241E+00
    0.383800000E+01    0.280626875E+00
N     S 
    0.904600000E+04    0.777446797E-05
    0.309300000E+03    0.300742072E-03
    0.877300000E+02   -0.280016549E-02
    0.285600000E+02   -0.989708505E-02
    0.102100000E+02   -0.114331114E+00
    0.383800000E+01   -0.118162383E+00
    0.746600000E+00    0.109786885E+01
N     S 
    0.224800000E+00    0.100000000E+01
N     P 
    0.135500000E+02    0.589056768E-01
    0.291700000E+01    0.320461107E+00
    0.797300000E+00    0.753042062E+00
N     P 
    0.218500000E+00    0.100000000E+01
N     D 
    0.817000000E+00    0.100000000E+01
'''),
'N2': gto.basis.parse('''
N     S 
    0.904600000E+04    0.701708743E-03
    0.135700000E+04    0.540299880E-02
    0.309300000E+03    0.274729510E-01
    0.877300000E+02    0.103514580E+00
    0.285600000E+02    0.279586579E+00
    0.102100000E+02    0.451317241E+00
    0.383800000E+01    0.280626875E+00
N     S 
    0.904600000E+04    0.777446797E-05
    0.309300000E+03    0.300742072E-03
    0.877300000E+02   -0.280016549E-02
    0.285600000E+02   -0.989708505E-02
    0.102100000E+02   -0.114331114E+00
    0.383800000E+01   -0.118162383E+00
    0.746600000E+00    0.109786885E+01
N     S 
    0.224800000E+00    0.100000000E+01
N     P 
    0.135500000E+02    0.589056768E-01
    0.291700000E+01    0.320461107E+00
    0.797300000E+00    0.753042062E+00
N     P 
    0.218500000E+00    0.100000000E+01
N     D 
    0.817000000E+00    0.100000000E+01
''')}

# Remember to check the charge and spin
mol.charge = 0
mol.spin = 0
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
alpha_coeff = fch2py('uhf-0.90.fch', nbf, nif, Sdiag, 'a')
beta_coeff = fch2py('uhf-0.90.fch', nbf, nif, Sdiag, 'b')
mf.mo_coeff = (alpha_coeff, beta_coeff)
# read done

# check if input MOs are orthonormal
check_orthonormal(nbf, nif, mf.mo_coeff[0], S)
check_orthonormal(nbf, nif, mf.mo_coeff[1], S)

#dm = mf.make_rdm1()
#mf.max_cycle = 10
#mf.kernel(dm)

