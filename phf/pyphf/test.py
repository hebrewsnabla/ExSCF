#import phf
from pyscf import gto, scf
import scipy
import numpy as np
from fch2py import fch2py
import sys

import util

xyz = sys.argv[1]
fch = sys.argv[2]
bas = xyz[:-3] + 'bas'

debug = False
cut_no = False

np.set_printoptions(precision=8, linewidth=160, suppress=True)
if debug:
    np.set_printoptions(precision=15, linewidth=200, suppress=False)

mol = gto.Mole()
#mol.atom = '''H 0. 0. 0.; H 0. 0. 2.'''
with open(xyz, 'r') as f:
    mol.atom = f.read()
print(mol.atom)
mol.basis = bas
mol.output = 'test.pylog'
mol.verbose = 4
mol.build()

mf = scf.UHF(mol)
#mf.init_guess = '1e'
mf.init_guess_breaksym = True
mf.max_cycle = 1
mf.kernel()

# read MOs from .fch(k) file
nbf = mf.mo_coeff[0].shape[0]
nif = mf.mo_coeff[0].shape[1]
S = mol.intor_symmetric('int1e_ovlp')
Sdiag = S.diagonal()
alpha_coeff = fch2py(fch, nbf, nif, Sdiag, 'a')
beta_coeff  = fch2py(fch, nbf, nif, Sdiag, 'b')
mf.mo_coeff = (alpha_coeff, beta_coeff)
# read done

#mf.mo_coeff = (
#    np.array([[0.367,-0.0839,1.25,-0.204],
#             [0.708,-0.209,-1.11,0.336],
#             [0.0212,0.171,0.165,1.30],
#             [0.0467,0.903,0.0630,-1.04]]),
#    np.array([[0.0212,0.171,0.165,1.30],
#             [0.0467,0.903,0.0630,-1.04],
#             [0.367,-0.0839,1.25,-0.204],
#             [0.708,-0.209,-1.11,0.336]])
#             )
dm = mf.make_rdm1()
mf.max_cycle = 0
mf.kernel(dm)

#Ca, Cb = mf.mo_coeff
#print(Ca)
#print(Cb)

#S = mf.get_ovlp()

mf2 = util.SUHF(mf)
mf2.cut_no = cut_no

X = mf2.X
na, nb = mf2.nelec

max_cycle = 15
cyc = 0
conv = False
while(not conv):
    print('**** Cycle %d ****' % (cyc+1))
    old_suhf = mf2.E_suhf
    hcore = mf.get_hcore()
    hcore_ortho = np.einsum('ji,jk,kl->il', X, hcore, X)
    #print(hcore_ortho)
    
    if cyc==0:
        veff = mf.get_veff(dm = dm)
    else:
        dm_reg = np.einsum('ij,tjk,lk->til', X, mf2.dm_ortho, X)
        veff = mf.get_veff(dm = dm_reg)
    veff_ortho = np.einsum('ji,tjk,kl->til', X, veff, X)
    print('dm (ortho)')
    print(mf2.dm_ortho)
    #print(veff)
    #Fa, Fb = hcore + veff
    #Fa_ortho = np.einsum('ji,jk,kl->il', X, Fa, X)
    #Fb_ortho = np.einsum('ji,jk,kl->il', X, Fb, X)
    Fa_ortho, Fb_ortho = hcore_ortho + veff_ortho
    print('Fock (ortho)',Fa_ortho, Fb_ortho)
    F_ortho = Fa_ortho, Fb_ortho

    e_uhf, e_uhf_coul = scf.uhf.energy_elec(mf, mf2.dm_ortho, hcore_ortho, veff_ortho)
    print('E(UHF) = %12.6f' % e_uhf)

    dm_no, dm_expanded, no = util.find_NO(mf2, mf2.dm_ortho, na, nb)
    Dg, Ng, Pg = util.get_Ng(mf2.grids, no, dm_no, na+nb)
    Gg, Pg_ortho = util.get_Gg(mf2.mol, Pg, no, X)
    xg, yg, ciS = util.get_xg(mf2, no, na, nb, Ng)
    mf2.xg, mf2.ciS = xg, ciS
    #yg, ciS = util.get_yg(mf2, xg)
    trHg, ciH = util.get_H(mf2, hcore_ortho, no, Pg, Gg, xg)
    S2 = util.get_S2(mf2, Pg_ortho)
    Xg, Xg_int, Yg = util.get_Yg(mf2, Dg, Ng, dm_no, na+nb)
    Feff_ortho, H_suhf, F_mod_ortho = util.get_Feff(mf2, trHg, Gg, Ng, Pg, dm_no, Dg, na+nb, Yg, Xg, no, F_ortho, mf2.dm_ortho)
    E_suhf = mf.energy_nuc() + H_suhf
    mf2.E_suhf = E_suhf
    print('E(SUHF) = %15.8f' % E_suhf)
    mo_e, mf2.dm_ortho = util.Diag_Feff(F_mod_ortho, na, nb)

    if old_suhf is not None:
        if abs(E_suhf - old_suhf)<1e-8:
            conv = True
            print('SCf converged at cycle %d' %cyc)
        
    
    cyc += 1
    if cyc >= max_cycle:
        break

