import numpy as np
#import sympy as sym
import scipy
from pyscf import gto, scf, mp, lib, ao2mo
from pyphf import suscf, util2, noci, pt2
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
        self.vap = False

    def kernel(self):
        suhf = self.suhf
        #ump2 = mp.UMP2(suhf.guesshf)
        #ump2.kernel()
        self.e_elec_hf = suhf.guesshf.energy_elec()[0]
        ghf = copy.copy(suhf.guesshf).to_ghf()
        self.ghf = ghf
        ghf.converged=True
        gmp2 = mp.GMP2(ghf)
        gmp2.kernel()
        self.gmp2 = gmp2
        energy_00 = suhf.ciH
        norm_00 = suhf.ciS
        if self.vap:
            #F0 = pt2.defm_Fock(suhf.mol, suhf.hcore_ortho, suhf.dm_ortho, suhf.X )
            F0 = suhf.guesshf.get_fock(dm=suhf.dm_reg)
            F0mo = einsum('tji,tjk,tkl->til', suhf.mo_reg, F0, suhf.mo_reg)
            print(F0mo)
            self.F0mo = F0mo
        energy_01, norm_01 = get_e01(self, suhf, gmp2)
        print('e01 %.6f, S00 %.6f, S01 %.6f' % (energy_01, norm_00, norm_01))
        ecorr = energy_01 / (norm_00 + norm_01)
        self.e_corr = ecorr
        print('SUMP2 e_corr %.6f' % ecorr)

def u2g_2d(umat):
    shape = umat[0].shape
    gmat = np.zeros((shape[0]*2, shape[1]*2))
    for i in range(shape[0]):
        for j in range(shape[1]):
            gmat[2*i, 2*j] = umat[0,i,j]
            gmat[2*i+1, 2*j+1] = umat[1,i,j]
    return gmat

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
    Fov = u2g_2d(spt2.F0mo)[:nocc, nocc:]
    E01g = []
    S01g = []
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
        #print(co_ovlp)
        u = util2.stack22( u_oo, np.zeros((nocc, nvir)), np.zeros((nvir, nocc)), np.eye(nvir))
        co_t2 = get_co_t2(gmp2.t2, v_oo)
        eri_mo = get_mo_eri(spt2.ghf.mo_coeff, eri_ao)
        co_eri = get_co_eri(eri_mo, u_oo, nocc)

        term1 = get_term1(oo_diag, co_t2, ovlp_ov, spt2.e_elec_hf)
        if spt2.vap:
            co_Fov = einsum('pi, pa -> ia', v_oo, Fov)
            term15 = get_term15(co_eri, co_ovlp, oo_diag, co_t2, ovlp_vo, ovlp_ov, co_Fov)
        else:
            term15 = 0.0
        term2 = get_term2(co_eri, co_ovlp, oo_diag, co_t2, ovlp_vo, ovlp_ov)
        #print('term1 %.6f, term2 %.6f' %(term1, term2))
        print('term1 %.6f, term15 %.6f, term2 %.6f' %(term1, term15, term2))
        e_01 = term1 + term15 + term2
        norm_01 = term1 / spt2.e_elec_hf
        E01g.append(e_01)
        S01g.append(norm_01)
    print('E01g', E01g)
    print('S01g', S01g)
    E01 = suhf.integr_beta(np.array(E01g))
    S01 = suhf.integr_beta(np.array(S01g))
    return E01, S01

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

def get_term1(ovlp_oo_diag, co_t2, ovlp_ov, e_elec_hf):
    ovlp_oo_inv = 1. / ovlp_oo_diag
    term1 = 0.5 * e_elec_hf * np.prod(ovlp_oo_diag) * einsum('ijab,ia,jb,i,j->', co_t2, ovlp_ov, ovlp_ov, ovlp_oo_inv, ovlp_oo_inv)
    return term1

def get_term15(co_eri, co_ovlp, ovlp_oo_diag, co_t2, ovlp_vo, ovlp_ov, Fov):
    #nocc, nvir = co_eri.shape[:2]
    #co_eri_oovv = co_eri.transpose(0,2,1,3) - co_eri.transpose(0,2,3,1)
    ovlp_oo_inv = 1. / ovlp_oo_diag
    term15 = 0.5 *  np.prod(ovlp_oo_diag) * einsum('ijab,ia,jb,i,j,kc,k, kc->', co_t2, ovlp_ov, ovlp_ov, ovlp_oo_inv, ovlp_oo_inv, ovlp_ov, ovlp_oo_inv, Fov)
    return term15

def get_term2(co_eri, co_ovlp, ovlp_oo_diag, co_t2, ovlp_vo, ovlp_ov):
    nocc, nvir = co_eri.shape[:2]
    co_eri_oovv = co_eri.transpose(0,2,1,3) - co_eri.transpose(0,2,3,1)
    ovlp_oo_inv = 1. / ovlp_oo_diag
    
    # Case 1: i = k, j = l
    S_ijpr_matrix = get_S_ijpr_matrix(co_ovlp, ovlp_oo_diag, nocc, nvir)
    sumc_ijad = einsum('ijcd,ijca->ijad', co_eri_oovv, S_ijpr_matrix)
    sumb_ijad = einsum('ijab,ijdb->ijad', co_t2, S_ijpr_matrix)
    sumadij = einsum('ijad,ijad,i,j->', sumc_ijad, sumb_ijad, ovlp_oo_inv, ovlp_oo_inv)
    term21 = sumadij * 0.25 * np.prod(ovlp_oo_diag)
    # Case 2: i = k, j != l
    S_ipr_matrix = get_S_ipr_matrix(co_ovlp, ovlp_oo_diag, nocc, nvir)
    # 1st term.
    sumld_ic = einsum('ilcd,dl,l->ic', co_eri_oovv, ovlp_vo, ovlp_oo_inv)
    sumjb_ia = einsum('ijab,jb,j->ia', co_t2, ovlp_ov, ovlp_oo_inv)
    case2term1 = einsum('ic,ia,i,ica->', sumld_ic, sumjb_ia, ovlp_oo_inv, S_ipr_matrix)
    # 2nd term.
    sumd_ijc = einsum('ijcd,dj->ijc', co_eri_oovv, ovlp_vo)
    sumb_ija = einsum('ijab,jb->ija', co_t2, ovlp_ov)
    case2term2 = -einsum('ijc,ija,j,i,ica->', sumd_ijc, sumb_ija, 
                              np.power(ovlp_oo_inv, 2), ovlp_oo_inv, S_ipr_matrix)
    term22 = (case2term1 + case2term2) * np.prod(ovlp_oo_diag)
    #Case 3: i != k, j != l        
    # 1st term.
    sumklcd = einsum('klcd,ck,dl,k,l->', co_eri_oovv, ovlp_vo, ovlp_vo, 
                           ovlp_oo_inv, ovlp_oo_inv)
    sumijab = einsum('ijab,ia,jb,i,j->', co_t2, ovlp_ov, ovlp_ov, ovlp_oo_inv, ovlp_oo_inv)
    case3term1 = 0.25 * sumklcd * sumijab
    # 2nd term.
    sumlcd_i = einsum('ilcd,ci,dl,l->i', co_eri_oovv, ovlp_vo, ovlp_vo, ovlp_oo_inv)
    sumjab_i = einsum('ijab,ia,jb,j->i', co_t2, ovlp_ov, ovlp_ov, ovlp_oo_inv)
    case3term2 = -einsum('i,i,i->', sumlcd_i, sumjab_i, np.power(ovlp_oo_inv, 2))
    # 3rd term.
    sumcd_ij = einsum('ijcd,ci,dj->ij', co_eri_oovv, ovlp_vo, ovlp_vo)
    sumab_ij = einsum('ijab,ia,jb->ij', co_t2, ovlp_ov, ovlp_ov)
    case3term3 = 0.5 * einsum('ij,ij,i,j->', sumcd_ij, sumab_ij, 
                               np.power(ovlp_oo_inv, 2), np.power(ovlp_oo_inv, 2))
    
    term23 = (case3term1 + case3term2 + case3term3) * np.prod(ovlp_oo_diag)

    e_term2 = term21 + term22 + term23
    return e_term2

def get_S_ijpr_matrix(co_ovlp, ovlp_oo_diag, nocc, nvir):
    S_ijpr_mat = np.empty((nocc, nocc, nvir, nvir))
        
    for i in range(nocc):
        for j in range(nocc):
            inds = (i, j)
            _nocc = nocc - len(np.unique(inds))
            co_ovlp_ij = np.delete(co_ovlp, inds, axis=0) # Delete rows i, j
            co_ovlp_ij = np.delete(co_ovlp_ij, inds, axis=1) # Delete cols i, j
            co_ovlp_ov_ij = co_ovlp_ij[:_nocc, _nocc:]
            co_ovlp_vo_ij = co_ovlp_ij[_nocc:, :_nocc]
            co_ovlp_oo_diag_ij = np.delete(ovlp_oo_diag, inds, axis=0) # Delete elements i, j
            co_ovlp_oo_inv_ij = 1. / co_ovlp_oo_diag_ij
            
            # [ co_ovlp_ij.T ]_{pk} = [ co_ovlp_ij ]_{kp}
            S_ijpr_mat[i, j] = (co_ovlp[nocc:, nocc:] - 
                                einsum('pk,k,kr->pr', co_ovlp_vo_ij, co_ovlp_oo_inv_ij, co_ovlp_ov_ij))
    
    return S_ijpr_mat

def get_S_ipr_matrix(co_ovlp, ovlp_oo_diag, nocc, nvir):
    _nocc = nocc - 1
    S_ipr_mat = np.empty((nocc, nvir, nvir))
    
    for i in range(nocc):
        co_ovlp_i = np.delete(co_ovlp, i, axis=0) # Delete row i
        co_ovlp_i = np.delete(co_ovlp_i, i, axis=1) # Delete col i
        co_ovlp_ov_i = co_ovlp_i[:_nocc, _nocc:]
        co_ovlp_vo_i = co_ovlp_i[_nocc:, :_nocc]
        ovlp_oo_diag_i = np.delete(ovlp_oo_diag, i, axis=0) # Delete element i
        ovlp_oo_inv_i = 1. / ovlp_oo_diag_i
        
        # [ co_ovlp_i.T ]_{pk} = [ co_ovlp_i ]_{kp}
        S_ipr_mat[i] = (co_ovlp[nocc:, nocc:] - 
                       einsum('pk,k,kr->pr', co_ovlp_vo_i, ovlp_oo_inv_i, co_ovlp_ov_i))
    
    return S_ipr_mat