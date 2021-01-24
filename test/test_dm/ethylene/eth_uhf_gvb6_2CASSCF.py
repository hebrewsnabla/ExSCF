from pyscf import gto, scf
from fch2py import fch2py
from ortho import check_orthonormal
from pyscf import mcscf, lib
from py2fch import py2fch
from shutil import copyfile
import numpy as np

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

mf = scf.RHF(mol)
mf.max_cycle = 1
mf.max_memory = 1000 # MB
mf.kernel()

# read MOs from .fch(k) file
nbf = mf.mo_coeff.shape[0]
nif = mf.mo_coeff.shape[1]
S = mol.intor_symmetric('int1e_ovlp')
Sdiag = S.diagonal()
mf.mo_coeff = fch2py('eth_uhf_uno_asrot2gvb6_s.fch', nbf, nif, Sdiag, 'a')
# read done

# check if input MOs are orthonormal
check_orthonormal(nbf, nif, mf.mo_coeff, S)

#dm = mf.make_rdm1()
#mf.max_cycle = 10
#mf.kernel(dm)

mc = mcscf.CASSCF(mf,2,(1,1))
mc.fcisolver.max_memory = 500 # MB
mc.max_memory = 500 # MB
mc.max_cycle = 200
mc.fcisolver.spin = 0
mc.fix_spin_(ss=0)
mc.natorb = True
mc.verbose = 5
mc.kernel()

dm = mc.make_rdm1s()
print(dm[0]-dm[1])
# save NOs into .fch file
copyfile('eth_uhf_uno_asrot2gvb6_s.fch', 'eth_uhf_gvb6_2CASSCF_NO.fch')
noon = np.zeros(nif)
py2fch('eth_uhf_gvb6_2CASSCF_NO.fch', nbf, nif, mc.mo_coeff, Sdiag, 'a', noon)
