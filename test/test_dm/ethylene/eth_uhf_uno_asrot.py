from pyscf import gto, scf
from fch2py import fch2py
from ortho import check_orthonormal
from pyscf import lib
import numpy as np
import os
from py2fch import py2fch
from uno import uno
from construct_vir import construct_vir
from lo import pm
from assoc_rot import assoc_rot
from shutil import copyfile

lib.num_threads(1)

mol = gto.M()
#   6 atom(s)
mol.atom = '''
C1              0.73236300       -0.00001200       -0.00002000
H1              1.30178800       -0.65477800        0.65483500
H2              1.30180201        0.65484800       -0.65476200
C2             -0.73236300       -0.00001200        0.00002000
H3             -1.30178800       -0.65477800       -0.65483600
H4             -1.30180201        0.65484800        0.65476200
'''

mol.basis = {
'C1': gto.basis.parse('''
C     S
    0.123840169E+04    0.551736505E-02
    0.186290050E+03    0.410888286E-01
    0.422511763E+02    0.182253821E+00
    0.116765579E+02    0.468284594E+00
    0.359305065E+01    0.445758173E+00
C     S
    0.402451474E+00    0.100000000E+01
C     S
    0.130901827E+00    0.100000000E+01
C     P
    0.946809706E+01    0.568883320E-01
    0.201035451E+01    0.312940593E+00
    0.547710047E+00    0.760650165E+00
C     P
    0.152686138E+00    0.100000000E+01
C     D
    0.800000000E+00    0.100000000E+01
'''),
'H1': gto.basis.parse('''
H     S
    0.130107010E+02    0.334854848E-01
    0.196225720E+01    0.234721871E+00
    0.444537960E+00    0.813770285E+00
H     S
    0.121949620E+00    0.100000000E+01
H     P
    0.800000000E+00    0.100000000E+01
'''),
'H2': gto.basis.parse('''
H     S
    0.130107010E+02    0.334854848E-01
    0.196225720E+01    0.234721871E+00
    0.444537960E+00    0.813770285E+00
H     S
    0.121949620E+00    0.100000000E+01
H     P
    0.800000000E+00    0.100000000E+01
'''),
'C2': gto.basis.parse('''
C     S
    0.123840169E+04    0.551736505E-02
    0.186290050E+03    0.410888286E-01
    0.422511763E+02    0.182253821E+00
    0.116765579E+02    0.468284594E+00
    0.359305065E+01    0.445758173E+00
C     S
    0.402451474E+00    0.100000000E+01
C     S
    0.130901827E+00    0.100000000E+01
C     P
    0.946809706E+01    0.568883320E-01
    0.201035451E+01    0.312940593E+00
    0.547710047E+00    0.760650165E+00
C     P
    0.152686138E+00    0.100000000E+01
C     D
    0.800000000E+00    0.100000000E+01
'''),
'H3': gto.basis.parse('''
H     S
    0.130107010E+02    0.334854848E-01
    0.196225720E+01    0.234721871E+00
    0.444537960E+00    0.813770285E+00
H     S
    0.121949620E+00    0.100000000E+01
H     P
    0.800000000E+00    0.100000000E+01
'''),
'H4': gto.basis.parse('''
H     S
    0.130107010E+02    0.334854848E-01
    0.196225720E+01    0.234721871E+00
    0.444537960E+00    0.813770285E+00
H     S
    0.121949620E+00    0.100000000E+01
H     P
    0.800000000E+00    0.100000000E+01
''')}

# Remember to check the charge and spin
mol.charge = 0
mol.spin = 0
mol.verbose = 4
mol.build()

mf = scf.UHF(mol)
mf.max_cycle = 1
mf.max_memory = 1000 # MB
mf.kernel()

# read MOs from .fch(k) file
nbf = mf.mo_coeff[0].shape[0]
nif = mf.mo_coeff[0].shape[1]
S = mol.intor_symmetric('int1e_ovlp')
Sdiag = S.diagonal()
alpha_coeff = fch2py('eth_uhf.fch', nbf, nif, Sdiag, 'a')
beta_coeff = fch2py('eth_uhf.fch', nbf, nif, Sdiag, 'b')
mf.mo_coeff = (alpha_coeff, beta_coeff)
# read done

# check if input MOs are orthonormal
check_orthonormal(nbf, nif, mf.mo_coeff[0], S)
check_orthonormal(nbf, nif, mf.mo_coeff[1], S)

dm = mf.make_rdm1()
mf.max_cycle = 10
mf.kernel(dm)

# transform UHF canonical orbitals to UNO
na = np.sum(mf.mo_occ[0]==1)
nb = np.sum(mf.mo_occ[1]==1)
idx, noon, alpha_coeff = uno(nbf, nif, na, nb, mf.mo_coeff[0], mf.mo_coeff[1], S)
alpha_coeff = construct_vir(nbf, nif, idx[1], alpha_coeff, S)
mf.mo_coeff = (alpha_coeff, beta_coeff)
# done transform

# save the UNO into .fch file
os.system('fch_u2r eth_uhf.fch')
os.rename('eth_uhf_r.fch', 'eth_uhf_uno.fch')
py2fch('eth_uhf_uno.fch', nbf, nif, mf.mo_coeff[0], Sdiag, 'a', noon)
# save done

# associated rotation
npair = np.int64((idx[1]-idx[0]-idx[2])/2)
idx2 = idx[0] + npair - 1
idx3 = idx2 + idx[2]
idx1 = idx2 - npair
idx4 = idx3 + npair
occ_idx = range(idx1,idx2)
vir_idx = range(idx3,idx4)
print(idx1, idx2, idx3, idx4)
occ_loc_orb = pm(mol.nbas, mol._bas[:,0], mol._bas[:,1],  mol._bas[:,3], mol.cart, nbf, npair, mf.mo_coeff[0][:,occ_idx], S, 'mulliken')
vir_loc_orb = assoc_rot(nbf, npair, mf.mo_coeff[0][:,occ_idx], occ_loc_orb, mf.mo_coeff[0][:,vir_idx])
mf.mo_coeff[0][:,occ_idx] = occ_loc_orb.copy()
mf.mo_coeff[0][:,vir_idx] = vir_loc_orb.copy()
# localization done

# save associated rotation MOs into .fch(k) file
copyfile('eth_uhf_uno.fch', 'eth_uhf_uno_asrot.fch')
noon = np.zeros(nif)
py2fch('eth_uhf_uno_asrot.fch', nbf, nif, mf.mo_coeff[0], Sdiag, 'a', noon)
