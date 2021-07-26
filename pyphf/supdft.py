from pyphf import suscf
from pyscf import dft
import pyscf.dft.numint as numint

import numpy as np
from functools import partial
import time

print = partial(print, flush=True)
einsum = partial(np.einsum, optimize=True)

class PDFT():
    def __init__(self, suhf):
        self.suhf = suhf
        self.ot = 'tpbe'
    def kernel(self):
        return kernel(self, self.suhf, self.ot)

def kernel(pdft, suhf, ot):
    dm1 = suhf.suhf_dm
    #dm2 = get_2CDM_from_2RDM(suhf.suhf_dm2, )

    #omega, alpha, hyb = ot._numint.rsh_and_hybrid_coeff(ot.otxc, spin=spin)
    Jg, Kg = suhf.get_JKg()
    #hfx = suhf.get_EX()
    get_H(suhf, suhf.hcore_ortho, suhf.no, suhf.Pg, suhf.Gg, Jg, Kg, suhf.xg)


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
    ciH0 = suhf.integr_beta(trHg0*xg)
    ciH1 = suhf.integr_beta(trHg1*xg)
    ciH1j = suhf.integr_beta(trHg1j*xg)
    ciH1k = suhf.integr_beta(trHg1k*xg)
    print('ciH', ciH0, ciH1, ciH1j, ciH1k)
    #suhf.trHg = trHg
    return ciH0, ciH1