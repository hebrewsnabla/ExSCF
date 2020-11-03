from pyscf import gto, scf
from fch2py import fch2py
from ortho import check_orthonormal

mol = gto.M()
#   2 atom(s)
mol.atom = '''
H1              0.00000000        0.00000000        0.00000000
H2              0.00000000        0.00000000        2.00000001
'''

mol.basis = {
'H1': gto.basis.parse('''
H     S 
    0.544717800E+01    0.156284979E+00
    0.824547240E+00    0.904690877E+00
H     S 
    0.183191580E+00    0.100000000E+01
'''),
'H2': gto.basis.parse('''
H     S 
    0.544717800E+01    0.156284979E+00
    0.824547240E+00    0.904690877E+00
H     S 
    0.183191580E+00    0.100000000E+01
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
alpha_coeff = fch2py('test_uhf.fch', nbf, nif, Sdiag, 'a')
beta_coeff = fch2py('test_uhf.fch', nbf, nif, Sdiag, 'b')
mf.mo_coeff = (alpha_coeff, beta_coeff)
# read done

# check if input MOs are orthonormal
check_orthonormal(nbf, nif, mf.mo_coeff[0], S)
check_orthonormal(nbf, nif, mf.mo_coeff[1], S)

#dm = mf.make_rdm1()
#mf.max_cycle = 10
#mf.kernel(dm)

