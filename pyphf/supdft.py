from pyphf import suscf, jk, sudft
#from pyscf import dft
import pyscf.dft.numint as numint

import numpy as np
from functools import partial
import time

print = partial(print, flush=True)
einsum = partial(np.einsum, optimize=True)

class PDFT():
    def __init__(self, suhf):
        self.suhf = suhf
        self.xc = 'pbe'
    def kernel(self):
        return kernel(self, self.suhf)

def kernel(pdft, suhf):
    print('\n******** %s ********' % pdft.__class__)
    mol = suhf.mol
    dm1 = suhf.suhf_dm
    dmdefm = suhf.dm_reg
    #dm2 = get_2CDM_from_2RDM(suhf.suhf_dm2, )
    print('energy decomposition')
    if suhf.debug:
        old_decomp(suhf, dm1)
    new_decomp(suhf, dm1)
    grids = sudft.set_grids(mol)
    ni = numint.NumInt()
    n, exc, vxc = ni.nr_uks(mol, grids, pdft.xc, dmdefm)
    print('E_xcdft %.6f' % exc)

def new_decomp(suhf, dm1):
    enuc = suhf.energy_nuc
    print('E_nuc %.6f' % enuc)
    dm1t = dm1[0] + dm1[1]
    Ecore = np.trace(np.dot(suhf.hcore_reg, dm1t))
    print('E_core %.6f' % Ecore)
    vj, vk = jk.get_jk(suhf.mol, dm1)
    veffj = vj[0] + vj[1] 
    veffk = -vk
    Ej = np.trace(np.dot(veffj, dm1[0]) + np.dot(veffj, dm1[1])) * 0.5
    Ek = np.trace(np.dot(veffk[0], dm1[0]) + np.dot(veffk[1], dm1[1])) * 0.5
    print('E_j %.6f' % Ej)
    print('E_k %.6f' % Ek)
    Ejk = Ej + Ek
    if suhf.debug:
        print('E_jk %.6f' % Ejk)
    Ec = suhf.E_suhf - enuc - Ecore - Ejk
    print('E_c %.6f' % Ec)

def old_decomp(suhf, dm1):
    dm1t = dm1[0] + dm1[1]
    enuc = suhf.energy_nuc
    print('E_nuc %.6f' % enuc)
    #omega, alpha, hyb = ot._numint.rsh_and_hybrid_coeff(ot.otxc, spin=spin)
    Jg, Kg = suhf.get_JKg()
    #hfx = suhf.get_EX()
    E0, Ejk, E1j, E1k = get_H(suhf, suhf.hcore_ortho, suhf.no, suhf.Pg, suhf.Gg, Jg, Kg, suhf.xg)
    E_suhf = enuc + E0 + Ejk
    print('E_suhf %.6f' % E_suhf)

def get_H(suhf, hcore_ortho, no, Pg, Gg, Jg, Kg, xg):
    #print(hcore_ortho)
    hcore_ortho = np.vstack((
        np.hstack((hcore_ortho, np.zeros(hcore_ortho.shape))),
        np.hstack((np.zeros(hcore_ortho.shape), hcore_ortho))
    ))
    hcore_no = einsum('ji,jk,kl->il', no, hcore_ortho, no)
    suhf.hcore_no = hcore_no
    if suhf.debug2:
        print(hcore_no)
    trHg0 = np.zeros(len(Pg))
    trHg1 = np.zeros(len(Pg))
    trHg1j = np.zeros(len(Pg))
    trHg1k = np.zeros(len(Pg))

    for i, pg in enumerate(Pg):
        H0 = np.trace(np.dot(hcore_no, pg)) 
        H1 = 0.5 * np.trace(np.dot(Gg[i], pg))
        H1j = 0.5 * np.trace(np.dot(Jg[i], pg))
        H1k = 0.5 * np.trace(np.dot(Kg[i], pg))
        #H = H * xg[i]
        trHg0[i] = H0 
        trHg1[i] = H1
        trHg1j[i] = H1j
        trHg1k[i] = H1k
        #print(i, H*xg[i])
    H0 = suhf.integr_beta(trHg0, fac='xg')
    H1 = suhf.integr_beta(trHg1, fac='xg')
    H1j = suhf.integr_beta(trHg1j, fac='xg')
    H1k = suhf.integr_beta(trHg1k, fac='xg')
    #H0 = suhf.integr_beta(trHg, fac='xg')
    #print('H ', H0, H1, H1j, H1k)
    print('E_core %.6f' % H0)
    print('E_jk %.6f' % H1)
    print('E_j %.6f' % H1j)
    print('E_k %.6f' % H1k)
    #print('ciH', ciH0, ciH1, ciH1j, ciH1k)
    #suhf.trHg = trHg
    return H0, H1, H1j, H1k