import numpy as np
#import sympy as sym
import scipy
from pyscf import gto, scf, mp, lib, ao2mo
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
        S = suhf.ovlp
        self.ao_ovlp = scipy.linalg.block_diag(S, S)

    def kernel(self):
        suhf = self.suhf
        #ump2 = mp.UMP2(suhf.guesshf)
        #ump2.kernel()
        ghf = copy.copy(suhf.guesshf).to_ghf()
        self.ghf = ghf
        gmp2 = mp.GMP2(ghf)
        gmp2.kernel()
        self.gmp2 = gmp2
        energy_00 = suhf.ciH
        norm_00 = suhf.ciS
        energy_01, norm_01 = get_e01(self, suhf, gmp2)

def get_e01(spt2, suhf, gmp2):
    mo = suhf.mo_reg
    grids = suhf.grids
    na, nb = suhf.nelec
    nocc = na+nb
    nvir = 2*suhf.norb - nocc
    mo_expd = suscf.expd(mo, suhf.mo_occ)
    Dg, Mg, xg, Pg = noci.get_DxP(suhf, suhf.norb, mo_expd, mo_expd, nocc)
    #mo_v = mo_expd[:,na+nb:]
    #mo_o = mo_expd[:,:na+nb]
    ghf = spt2.ghf
    mo_v = ghf.mo_coeff[:,nocc:]
    mo_o = ghf.mo_coeff[:,:nocc]
    eri_ao = suhf.mol.intor('int2e', aosym='s8')
    S = spt2.ao_ovlp
    for i, dg in enumerate(Dg):
        mg = einsum('ji,jk,kl->il', mo_o, dg, mo_o)
        u, s, vt = np.linalg.svd(mg)
        u_oo = u
        oo_diag = s
        v_oo = vt.T
        ovlp_vo = reduce(np.dot, (mo_v.T, S, dg, mo_o, v_oo))
        ovlp_ov = reduce(np.dot, (u_oo.T, mo_o.T, S, dg, mo_v))
        ovlp_oo = np.diag(oo_diag)
        ovlp_vv = reduce(np.dot, (mo_v.T, S, dg, mo_v))
        co_ovlp = util2.stack22(ovlp_oo, ovlp_ov, ovlp_vo, ovlp_vv)
        print(co_ovlp)
        u = util2.stack22( u_oo, np.zeros((nocc, nvir)), np.zeros((nvir, nocc)), np.eye(nvir))
        co_t2 = get_co_t2(gmp2.t2, v_oo)
        eri_mo = get_mo_eri(spt2.ghf.mo_coeff, eri_ao,)
        co_eri = get_co_eri(eri_mo, u_oo, nocc)

def get_co_t2(t2_mo, v_oo):
    return einsum('pi,qj, pqab -> ijab', v_oo, v_oo, t2_mo)

def get_co_eri(eri_mo, u_oo, nocc):
    #mo_eri_oovv = eri_mo.ovov.transpose(0,2,1,3)
    mo_eri_ovov = eri_mo[:nocc, nocc:, :nocc, nocc:]
    co_eri_ovov = einsum('pi,qj,paqb->iajb', u_oo, u_oo, mo_eri_ovov)
    #co_eri = co_eri_oovv.transpose(0,2,1,3)
    return co_eri_ovov

def get_mo_eri(mo_coeff, ao_eri):
    nao = mo_coeff.shape[0]//2
    nmo = mo_coeff.shape[1]
    moa = mo_coeff[:nao]
    mob = mo_coeff[nao:]
    mo_eri = ao2mo.kernel(ao_eri, (moa,moa,moa,moa), compact=False)
    mo_eri += ao2mo.kernel(ao_eri, (mob,mob,mob,mob), compact=False)
    mo_eri += ao2mo.kernel(ao_eri, (moa,moa,mob,mob), compact=False)
    mo_eri += ao2mo.kernel(ao_eri, (mob,mob,moa,moa), compact=False)
    mo_eri = mo_eri.reshape(nmo, nmo, nmo, nmo)

    return mo_eri
