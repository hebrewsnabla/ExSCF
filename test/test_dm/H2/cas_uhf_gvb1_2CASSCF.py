from pyscf import gto, scf
from fch2py import fch2py
from ortho import check_orthonormal
from pyscf import mcscf, lib
from py2fch import py2fch
from shutil import copyfile
import numpy as np

lib.num_threads(1)

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

mf = scf.RHF(mol)
mf.max_cycle = 1
mf.max_memory = 4000 # MB
mf.kernel()

# read MOs from .fch(k) file
nbf = mf.mo_coeff.shape[0]
nif = mf.mo_coeff.shape[1]
S = mol.intor_symmetric('int1e_ovlp')
Sdiag = S.diagonal()
mf.mo_coeff = fch2py('cas_uhf_uno_asrot2gvb1_s.fch', nbf, nif, Sdiag, 'a')
# read done

# check if input MOs are orthonormal
check_orthonormal(nbf, nif, mf.mo_coeff, S)

#dm = mf.make_rdm1()
#mf.max_cycle = 10
#mf.kernel(dm)

mc = mcscf.CASSCF(mf,2,(1,1))
mc.fcisolver.max_memory = 2000 # MB
mc.max_memory = 2000 # MB
mc.max_cycle = 200
mc.fcisolver.spin = 0
mc.fix_spin_(ss=0)
mc.natorb = True
mc.verbose = 5
mc.kernel()

print(mc.make_rdm1() / 2.0)
print(mc.mo_coeff)
print(mc.mo_occ / 2.0)
# save NOs into .fch file
copyfile('cas_uhf_uno_asrot2gvb1_s.fch', 'cas_uhf_gvb1_2CASSCF_NO.fch')
noon = np.zeros(nif)
py2fch('cas_uhf_gvb1_2CASSCF_NO.fch',nbf,nif,mc.mo_coeff,Sdiag,'a',noon)
