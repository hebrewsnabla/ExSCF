import numpy as np
#import sympy as sym
#import scipy
from pyscf import gto, scf, mp, lib
from pyphf import suscf, util2, noci
from pyphf.suscf import count0, eig
#import os, sys
from functools import partial, reduce
import time
import copy

print = partial(print, flush=True)
einsum = partial(np.einsum, optimize=True)


class EMP2():
    def __init__(self, suhf):
        self.suhf = suhf
        self.mol = suhf.mol

        #self.cut_no = False
        self.verbose = 4
        #self.debug = False
        self.output = None

        self.norb = suhf.norb
        self.nelec = suhf.nelec
        na, nb = self.nelec
        #self.occ = na + nb
        #self.vir = self.norb - occ
    def kernel(self):
        suhf = self.suhf
        ump2 = mp.UMP2(suhf.guesshf)
        ump2.kernel()
        energy_00 = suhf.ciH
        norm_00 = suhf.ciS
        energy_01, norm_01 = get_e01(suhf, ump2)

def get_e01(suhf, ump2):
    mo = suhf.mo_reg
    grids = suhf.grids
    na, nb = suhf.nelec
    mo_expd = suscf.expd(mo, suhf.mo_occ)
    Dg, Mg, xg, Pg = noci.get_DxP(suhf, suhf.norb, mo_expd, mo_expd, na+nb)
    mo_v = mo_expd[:,na+nb:]
    mo_o = mo_expd[:,:na+nb]
    for x in xg:
        u, s, vt = np.linalg.svd(x)
        u_oo = u
        oo_diag = s
        v_oo = vt.T
        ovlp_vo = reduce(np.dot, (mo_v.T, S, mo_o, v_oo))
        ovlp_ov = reduce(np.dot, (u_oo.T, mo_o.T, S, mo_v))
        ovlp_oo = np.diag(oo_diag)
        ovlp_vv = reduce(np.dot, (mo_v.T, S, mo_v))
        co_ovlp = util2.stack22(ovlp_oo, ovlp_ov, ovlp_vo, ovlp_vv)