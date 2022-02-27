import numpy as np
#import sympy as sym
#import scipy
from scipy import linalg
from pyscf import gto, scf, mp, lib, ao2mo
from pyphf import suscf, util2, noci, pt2
from pyphf.suscf import count0, eig
from pyphf.elem import make_s02, make_s22
from pyphf.timing import timing
#import os, sys
from functools import partial, reduce
#import time
import copy

print = partial(print, flush=True)
einsum = lib.einsum
#einsum = partial(np.einsum, optimize=True)
#svd = partial(linalg.svd, lapack_driver='gesvd')

def svd(m):
    with lib.with_omp_threads(1):
        return linalg.svd(m, lapack_driver='gesvd')

class EMP2():
    def __init__(self, suhf):
        self.suhf = suhf
        self.mol = suhf.mol

        #self.cut_no = False
        self.verbose = suhf.verbose
        self.debug = suhf.debug
        self.debug2 = suhf.debug2
        self.output = None

        self.norb = suhf.norb
        self.nelec = suhf.nelec
        na, nb = self.nelec
        #self.occ = na + nb
        #self.vir = self.norb - occ
        S = suhf.ovlp
        self.ao_ovlp = linalg.block_diag(S, S)
        self.vap = True
        self.use_det = False
        self.do_sc = True
        self.do_biort = True
        self.do_15 = False
        #self.use_tblis = False
        self.frozen = 0

    def dump_flags(self):
        print('\n******** %s ********' % self.__class__)
        print('Do semi-canonicalization: %r' % self.do_sc)
        print('Do biorthogonalization: %r' % self.do_biort)
        print('Found TBLIS: %r' % lib.numpy_helper.FOUND_TBLIS)


    def kernel(self):

        self.dump_flags()
        np.set_printoptions(precision=6, linewidth=160, suppress=True)
        #if self.debug:
        #    np.set_printoptions(precision=10, linewidth=200, suppress=False)
        suhf = self.suhf
        #ump2 = mp.UMP2(suhf.guesshf)
        #ump2.kernel()
        #if isinstance(suhf.guesshf, scf.uhf.UHF):
        #    self.e_elec_hf = suhf.guesshf.energy_elec()[0]
        #else:
        #    suhf.guesshf.xc = 'hf'
        #    self.e_elec_hf = suhf.guesshf.energy_elec()[0]

        _hf = copy.copy(suhf.guesshf)
        _hf.mo_coeff = suhf.mo_reg
        ghf = _hf.to_ghf()
        ghf.converged=True
        self.ghf = ghf
        orbspin = ghf.mo_coeff.orbspin
        print('mo orbspin\n', orbspin)
        gmp2 = mp.GMP2(ghf)
        gmp2.frozen = self.frozen
        if not self.vap:
            gmp2.kernel()
        #e0gmp = gmp2.kernel()[0]
        #print('GMP2 before SC', e0gmp)
        self.gmp2 = gmp2
        energy_00 = suhf.ciH
        norm_00 = suhf.ciS
        do_sc = self.do_sc
        na, nb = suhf.nelec
        if self.vap:
            #F0 = pt2.defm_Fock(suhf.mol, suhf.hcore_ortho, suhf.dm_ortho, suhf.X )
            F0 = _hf.get_fock(dm=suhf.dm_reg)
            F0mo = einsum('tji,tjk,tkl->til', suhf.mo_reg, F0, suhf.mo_reg)
            dump_Fock(F0mo[0], suhf.core, suhf.act, na)
            dump_Fock(F0mo[1], suhf.core, suhf.act, nb)
            
            if self.debug:
                print('F0mo',F0mo)
            if do_sc:
                F0mo_sc, voo_sc, vvv_sc = semi_cano(F0mo, na, nb, orbspin)
                self.F0mo = F0mo_sc
                self.vsc = (voo_sc, vvv_sc)
                dump_Fock(self.F0mo, 2*suhf.core, 2*suhf.act, na+nb)
            else:
                self.F0mo = u2g_2d(F0mo, orbspin)
            if self.debug:
                print('F0mo',self.F0mo)
        if self.use_det:
            energy_01, norm_01 = get_e01_det(self, suhf, gmp2)
        else:
            energy_01, norm_01, e01terms = get_e01(self, suhf, gmp2)
        print('e01 %.6f, S00 %.6f, S01 %.6f' % (energy_01, norm_00, norm_01))
        #print(e01terms)
        print('e01.02 %.6f, e01.12 %.6f, e01.22 %.6f' % (e01terms[0], e01terms[1], e01terms[2]))
        #ecorr = energy_01 / (norm_00 + norm_01)
        E0 = suhf.E_suhf - suhf.energy_nuc
        ecorr1 = - E0 * norm_01/norm_00 
        ecorr2 = energy_01 / norm_00
        print('E0 %.6f' % E0)
        print('-E0<0|P|1>/<0|P|0> = %.6f, <0|HP|1>/<0|P|0> = %.6f' % (ecorr1, ecorr2))
        ecorr = ecorr1 + ecorr2
        self.e_corr = ecorr
        #ecorr_approx = - self.e_elec_hf * norm_01/norm_00 + (e01terms[0] + gmp2.e_corr) / norm_00
        print('SUMP2 e_corr %.6f' % ecorr)
        #print('approx. SUMP2 e_corr %.6f' % ecorr_approx)

def dump_Fock(F, core, act, occ):
    offdiago = abs(np.triu(F[:occ,:occ],1)).max()
    offdiagv = abs(np.triu(F[occ:,occ:],1)).max()
    print('max offdiag ', offdiago, offdiagv)
    print(F[core:core+act, core:core+act])

def semi_cano(F, na, nb, orbspin):
    orbspin_o = orbspin[:na+nb]
    orbspin_v = orbspin[na+nb:]
    Fasc, va_oo, va_vv = _semi_cano(F[0], na)
    Fbsc, vb_oo, vb_vv = _semi_cano(F[1], nb)
    #voo_sc = scipy.linalg.block_diag(va_oo, vb_oo)
    #vvv_sc = scipy.linalg.block_diag(va_vv, vb_vv)
    voo_sc = u2g_2d( (va_oo, vb_oo), orbspin_o )
    vvv_sc = u2g_2d( (va_vv, vb_vv), orbspin_v )
    print('voo_sc,vvv_sc\n', voo_sc, '\n', vvv_sc)
    FG = u2g_2d((Fasc,Fbsc), orbspin)
    return FG, voo_sc, vvv_sc

def _semi_cano(F, na):
    Foo = F[:na,:na]
    Fov = F[:na,na:]
    Fvo = F[na:,:na]
    Fvv = F[na:,na:]
    print('Foo,Fvv', '\n', Foo, '\n', Fvv)
    u, s, vt = svd(Foo)
    if Foo[0,0]*s[0] < 0:
        s = -s
        u_oo = u
        v_oo = -vt.T
    else:
        s = np.flip(s)
        u_oo = np.flip(u, axis=1)
        v_oo = np.flip(vt.T, axis=1)
    u2, s2, vt2 = svd(Fvv)
    if Fvv[0,0]*s2[0] < 0:
        s2 = -s2
        u_vv = u2
        v_vv = -vt2.T
    else:
        s2 = np.flip(s2) 
        u_vv = np.flip(u2, axis=1)
        v_vv = np.flip(vt2.T, axis=1)
    Fscvv = np.diag(s2)
    Fscvo = reduce(np.dot, (u_vv.T, Fvo, v_oo))
    Fscov = reduce(np.dot, (u_oo.T, Fov, v_vv))
    Fscoo = np.diag(s)
    #ovlp_vv = reduce(np.dot, (mo_v.T, S, dg, mo_v))
    Fsc = util2.stack22(Fscoo, Fscov, Fscvo, Fscvv)
    #print(Fsc)
    #print(v_oo, '\n', v_vv)
    return Fsc, v_oo, v_vv

def u2g_2d(umat, orbspin):
    shape = umat[0].shape
    gmat = np.zeros((shape[0]*2, shape[1]*2))
    gmat[np.ix_(orbspin==0, orbspin==0)] = umat[0]
    gmat[np.ix_(orbspin==1, orbspin==1)] = umat[1]
    return gmat

def get_e01_new(spt2, suhf, gmp2):
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
    print('mo\n', ghf.mo_coeff)
    eri_ao = suhf.mol.intor('int2e', aosym='s8')
    S = spt2.ao_ovlp
    if spt2.vap:
        Fov = spt2.F0mo[:nocc, nocc:]
        mo_o = np.dot(mo_o, spt2.vsc[0])
        mo_v = np.dot(mo_v, spt2.vsc[1])
        mo_energy= spt2.F0mo.diagonal() 
        print('mo_e', spt2.F0mo.diagonal())
        print('mo_sc\n', np.hstack((mo_o, mo_v)))
    E01g = []
    S01g = []
    for i, dg in enumerate(Dg):
        mg = einsum('ji,jk,kl->il', mo_o, dg, mo_o)
        print(mg)
        u, s, vt = svd(mg)
        u_oo = u
        oo_diag = s
        v_oo = vt.T
        mo_o = np.dot(mo_o, v_oo)
        mo_o2 = np.dot(mo_o, u_oo)
        ovlp_vo = reduce(np.dot, (mo_v.T, S, dg, mo_o, v_oo))
        ovlp_ov = reduce(np.dot, (u_oo.T, mo_o.T, S, dg, mo_v))
        ovlp_oo = np.diag(oo_diag)
        ovlp_vv = reduce(np.dot, (mo_v.T, S, dg, mo_v))
        u_vv = v_vv = None
        if spt2.vap:
            u2, s2, vt2 = svd(ovlp_vv)
            u_vv = u2
            v_vv = vt2.T
            ovlp_vv = np.diag(s2)
            ovlp_vo = np.dot(u_vv.T, ovlp_vo)
            ovlp_ov = np.dot(ovlp_ov, v_vv) 
            #print(ovlp_vv)
        co_ovlp = util2.stack22(ovlp_oo, ovlp_ov, ovlp_vo, ovlp_vv)
        #if spt2.debug:
        print(co_ovlp)
        print(v_oo, '\n', u_oo)
        print(v_vv, '\n', u_vv)
        u = util2.stack22( u_oo, np.zeros((nocc, nvir)), np.zeros((nvir, nocc)), np.eye(nvir))
        co_t2 = get_co_t2(gmp2.t2, v_oo, v_vv)
        eri_mo = get_mo_eri(np.hstack((mo_o, mo_v)), eri_ao)
        co_eri = get_co_eri(eri_mo, v_oo, nocc, v_vv)

        term1 = get_term1(oo_diag, co_t2, ovlp_ov, spt2.e_elec_hf)
        if spt2.vap:
            co_Fov = einsum('pi, pa, aq -> iq', u_oo, Fov, v_vv)
            term15 = get_term15(co_eri, co_ovlp, oo_diag, co_t2, ovlp_vo, ovlp_ov, co_Fov, ovlp_vv)
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

def get_e01_det(spt2, suhf, gmp2):
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
    print('mo\n', ghf.mo_coeff)
    #eri_ao = suhf.mol.intor('int2e', aosym='s8')
    S = spt2.ao_ovlp
    do_sc =False
    if spt2.vap and do_sc:
        Fov = spt2.F0mo[:nocc, nocc:]
        print('Fov before', Fov)
        mo_o = np.dot(mo_o, spt2.vsc[0])
        mo_v = np.dot(mo_v, spt2.vsc[1])
        mo_coeff = np.hstack((mo_o, mo_v))
        gmp2.kernel(mo_energy= spt2.F0mo.diagonal() , mo_coeff=np.hstack((mo_o, mo_v)))
        print('GMP2 after SC', gmp2.e_corr)
        print('mo_e', spt2.F0mo.diagonal())
        print('mo_sc\n', np.hstack((mo_o, mo_v)))
    elif spt2.vap:
        Fov = spt2.F0mo[:nocc, nocc:]
        mo_coeff = ghf.mo_coeff
        gmp2.kernel(mo_energy= spt2.F0mo.diagonal() , mo_coeff=mo_coeff)
    else:
        mo_coeff = ghf.mo_coeff
    #exit()
    E01g = []
    S01g = []
    diagv = True
    mo_occ = ghf.mo_occ
    #for i in range(nocc):
    #    for j in range(nocc):
    #        for a in range(nvir):
    #            for b in range(a+1,nvir):
    #                t = gmp2.t2[i,j,a,b]
    #                if abs(t) > 1e-6:
    #                    print([i,j,a+nocc,b+nocc], t)
    for i, dg in enumerate(Dg):
        #mg = reduce(np.dot,( mo_o.T ,S, dg, mo_o))
        #print(mg)
        
        s02 = make_s02(dg, mo_coeff, nocc, nvir, mo_occ, S)
        #co_t2 = get_co_t2(gmp2.t2, v_oo, v_vv)
        eri_mo = gmp2.ao2mo(mo_coeff).oovv
        #co_eri = get_co_eri(eri_mo, u_oo, nocc, u_vv)

        #term1 = get_term1(oo_diag, co_t2, ovlp_ov, spt2.e_elec_hf)
        print(s02.shape, gmp2.t2.shape)
        term1 = spt2.e_elec_hf * einsum('ijab,ijab->', s02, gmp2.t2)
        term15 = 0.0
        if spt2.vap:
            term15 = get_term15_det( gmp2.t2, dg, mo_coeff, nocc, nvir, mo_occ, S, Fov)
        term2 = get_term2_det(eri_mo, gmp2.t2, dg, mo_coeff, nocc, nvir, mo_occ, S) 
        #print('term1,2 %.6f %.6f' %(term1, term2))
        #continue
        print('term1 %.6f, term15 %.6f, term2 %.6f' %(term1, term15, term2))
        e_01 = term1 + term15 + term2
        norm_01 = term1 / spt2.e_elec_hf
        E01g.append(e_01)
        S01g.append(norm_01)
    #exit()
    print('E01g', E01g)
    print('S01g', S01g)
    E01 = suhf.integr_beta(np.array(E01g))
    S01 = suhf.integr_beta(np.array(S01g))
    return E01, S01

def get_term15_det( t2, dg, mo_coeff, nocc, nvir, mo_occ, S, Fov):
    # St[kc] = S[kc, ijab] * t[ijab]
    St = np.zeros((nocc, nvir))
    for k in range(nocc):
        for c in range(nvir):
            s12 = make_s22(dg, mo_coeff, nocc, nvir, mo_occ, S, [k,None,c,None])
            St[k,c] = einsum('ijab,ijab->', s12, t2)
    term15 =  einsum('kl,kl', Fov, St)
    return term15
                    
def get_term2_det(eri_mo, t2, dg, mo_coeff, nocc, nvir, mo_occ, S):
    # St[klcd] = S[klcd, ijab] * t[ijab]
    St = np.zeros((nocc, nocc, nvir, nvir))
    for k in range(nocc):
        #for l in list(range(k)) + list(range(k+1,nocc)):
        for l in range(k):
            for c in range(nvir):
                #for d in list(range(c)) + list(range(c+1,nvir)):
                for d in range(c):
                    if abs(eri_mo[k,l,c,d]) < 1e-6: 
                        continue
                    Sklcd = make_s22(dg, mo_coeff, nocc, nvir, mo_occ, S, [k,l,c,d])
                    #St[k,l,c,d] = einsum('ab,ab->', Sklcd, t2[k,l,:,:])
                    St[k,l,c,d] = einsum('ijab,ijab->', Sklcd, t2)
    term2 =  einsum('klcd,klcd', eri_mo, St)
    return term2
                    

@timing
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
    if spt2.debug:
        print('mo\n', ghf.mo_coeff)
    #eri_ao = suhf.mol.intor('int2e', aosym='s8')
    S = spt2.ao_ovlp
    if spt2.vap:
        Fov = spt2.F0mo[:nocc, nocc:]
        if spt2.debug:
            print('Fov before', Fov)
        if spt2.do_sc:
            mo_o = np.dot(mo_o, spt2.vsc[0])
            mo_v = np.dot(mo_v, spt2.vsc[1])
            mo_coeff = np.hstack((mo_o, mo_v))
            mo_energy = spt2.F0mo.diagonal()
        else:
            mo_coeff = ghf.mo_coeff
            mo_energy = spt2.F0mo.diagonal()
        e_elec_hf = ghf.energy_elec(ghf.make_rdm1(mo_coeff=mo_coeff))[0]
        spt2.e_elec_hf = e_elec_hf
        print('e_elec_hf', e_elec_hf)
        do_gmp2(gmp2, mo_coeff, mo_energy, spt2.debug)
    else:
        mo_coeff = ghf.mo_coeff
    #exit()
    if spt2.frozen > 0:
        frozen = spt2.frozen
        mo_coeff = gmp2.mo_coeff
        mo_o = mo_o[:, frozen:]
        nocc = nocc - frozen
    E01g = []
    S01g = []
    E01gterm = []
    diagv = spt2.do_biort
    for i, dg in enumerate(Dg):
        ovlp0 = reduce(np.dot,( mo_coeff.T ,S, dg, mo_coeff))
        if spt2.debug:
            print('ovlp before', ovlp0)
        mg = reduce(np.dot,( mo_o.T ,S, dg, mo_o))
        if spt2.debug:
            print(mg)
        u, s, vt = svd(mg)
        u_oo = u
        oo_diag = s
        v_oo = vt.T
        ovlp_vo = reduce(np.dot, (mo_v.T, S, dg, mo_o, v_oo))
        ovlp_ov = reduce(np.dot, (u_oo.T, mo_o.T, S, dg, mo_v))
        ovlp_oo = np.diag(oo_diag)
        ovlp_vv = reduce(np.dot, (mo_v.T, S, dg, mo_v))
        co_or = reduce(np.dot, (dg, mo_o, v_oo))
        co_ol = np.dot(mo_o, u_oo)
        if spt2.debug:
            print('co\n', co_or, '\n', co_ol)
        u_vv = v_vv = None
        if spt2.vap and diagv:
            u2, s2, vt2 = svd(ovlp_vv)
            u_vv = u2
            v_vv = vt2.T
            ovlp_vv = np.diag(s2)
            ovlp_vo = np.dot(u_vv.T, ovlp_vo)
            ovlp_ov = np.dot(ovlp_ov, v_vv) 
            #print(ovlp_vv)
            co_vr = reduce(np.dot, (dg, mo_v, v_vv))
            co_vl = np.dot(mo_v, u_vv)
            co_r = np.hstack((co_or, co_vr))
        co_ovlp = util2.stack22(ovlp_oo, ovlp_ov, ovlp_vo, ovlp_vv)
        if spt2.debug:
            print('co_ovlp\n', co_ovlp)
            print(v_oo, '\n', u_oo)
            print(v_vv, '\n', u_vv)
        u = util2.stack22( u_oo, np.zeros((nocc, nvir)), np.zeros((nvir, nocc)), np.eye(nvir))
        co_t2 = get_co_t2(gmp2.t2, v_oo, v_vv)
        if False:
            dmg = ghf.make_rdm1(mo_coeff = co_r)
            fockg = ghf.get_fock(dm=dmg)
            print('fockg', fockg)
            co_t2 = gmp2.kernel(mo_energy= fockg.diagonal() , mo_coeff=co_r)[1]
        #eri_mo = get_mo_eri(np.hstack((mo_o, mo_v)), eri_ao)
        eri_mo = get_mo_eri(gmp2, mo_coeff)
        co_eri = get_co_eri(eri_mo, u_oo, nocc, u_vv)

        term1 = get_term1(oo_diag, co_t2, ovlp_ov, spt2.e_elec_hf, spt2.debug2)
        if spt2.vap and spt2.do_15:
            co_Fov = einsum('pi, pa -> ia', u_oo, Fov)
            if diagv: co_Fov = np.dot(co_Fov, v_vv)
            if spt2.debug:
                print('Fov\n', co_Fov)
            term15 = get_term15(co_eri, co_ovlp, oo_diag, co_t2, ovlp_vo, ovlp_ov, co_Fov, ovlp_vv, spt2.debug2)
        else:
            term15 = 0.0
        term2 = get_term2(co_eri, co_ovlp, oo_diag, co_t2, ovlp_vo, ovlp_ov)
        #print('term1 %.6f, term2 %.6f' %(term1, term2))
        print('term1 %.6f, term15 %.6f, term2 %.6f' %(term1, term15, term2))
        e_01 = term1 + term15 + term2
        norm_01 = term1 / spt2.e_elec_hf
        E01g.append(e_01)
        S01g.append(norm_01)
        E01gterm.append(np.array([term1, term15, term2]))
    print('E01g', E01g)
    print('S01g', S01g)
    E01 = suhf.integr_beta(np.array(E01g))
    S01 = suhf.integr_beta(np.array(S01g))
    E01term = suhf.integr_beta(np.array(E01gterm))
    #print(E01term)
    return E01, S01, E01term

@timing
def do_gmp2(gmp2, mo_coeff, mo_energy, dbg):
    gmp2._scf.mo_coeff = mo_coeff
    gmp2._scf.mo_energy = mo_energy
    gmp2.kernel()
    #print('GMP2 after SC', gmp2.e_corr)
    if dbg:
        print('mo_e', mo_energy)
        print('mo_sc\n', mo_coeff)


@timing
def get_co_t2(t2_mo, v_oo, v_vv=None):
    co_t2 = einsum('pi,qj, pqab -> ijab', v_oo, v_oo, t2_mo)
    if v_vv is not None:
        co_t2 = einsum('ac, bd, ijab  -> ijcd', v_vv, v_vv, co_t2)
    return co_t2

@timing
def get_co_eri(eri_mo, u_oo, nocc, u_vv=None):
    mo_eri_oovv = eri_mo.oovv
    
    #mo_eri_ovov = eri_mo[:nocc, nocc:, :nocc, nocc:]
    co_eri_oovv = einsum('pi,qj,pqab->ijab', u_oo, u_oo, mo_eri_oovv)
    #co_eri = co_eri_oovv.transpose(0,2,1,3)
    if u_vv is not None:
        co_eri_oovv = einsum('ac,bd,ijab -> ijcd', u_vv, u_vv, co_eri_oovv)
    return co_eri_oovv

@timing
def get_mo_eri(gmp2, mo_coeff):
#    nao = mo_coeff.shape[0]//2
#    nmo = mo_coeff.shape[1]
#    moa = mo_coeff[:nao]
#    mob = mo_coeff[nao:]
#    mo_eri = ao2mo.kernel(ao_eri, (moa,moa,moa,moa), compact=False)
#   mo_eri += ao2mo.kernel(ao_eri, (mob,mob,mob,mob), compact=False)
#   mo_eri += ao2mo.kernel(ao_eri, (moa,moa,mob,mob), compact=False)
#   mo_eri += ao2mo.kernel(ao_eri, (mob,mob,moa,moa), compact=False)
#    mo_eri = mo_eri.reshape(nmo, nmo, nmo, nmo)
    eri_mo = gmp2.ao2mo(mo_coeff)
    return eri_mo

@timing
def get_term1(ovlp_oo_diag, co_t2, ovlp_ov, e_elec_hf, dbg):
    ovlp_oo_inv = 1. / ovlp_oo_diag
    term1 = 0.5 * np.prod(ovlp_oo_diag) * einsum('ijab,ia,jb,i,j->', co_t2, ovlp_ov, ovlp_ov, ovlp_oo_inv, ovlp_oo_inv)
    #print(term1)
    #np.set_printoptions(precision=10, linewidth=200, suppress=True)
    if dbg:
        occ = co_t2.shape[0]
        vir = co_t2.shape[2]
        print(ovlp_ov)
        for i in range(occ):
            for j in range(occ):
                for a in range(vir):
                    for b in range(vir):
                        t2 = co_t2[i,j,a,b]
                        n = ovlp_ov[i,a]*ovlp_ov[j,b]*ovlp_oo_inv[i]*ovlp_oo_inv[j]
                        if abs(n) > 1e-4:
                            print(i,j,a,b,' %2.6f %2.6f  %2.6f'%( t2, n, t2*n) )

    return term1*e_elec_hf

@timing
def get_term15(co_eri, co_ovlp, ovlp_oo_diag, co_t2, ovlp_vo, ovlp_ov, Fov, ovlp_vv, dbg):
    #nocc, nvir = co_eri.shape[:2]
    nocc = co_eri.shape[0]
    nvir = co_eri.shape[2]
    #co_eri_oovv = co_eri.transpose(0,2,1,3) - co_eri.transpose(0,2,3,1)
    ovlp_oo_inv = 1. / ovlp_oo_diag
    term15 = 0.0
    # Case 1: i = k
    #c1 = -0.5 * np.prod(ovlp_oo_diag) * einsum('ijab,ca,jb,i,j, ic->', 
    #           co_t2, ovlp_vv, ovlp_ov, ovlp_oo_inv, ovlp_oo_inv, Fov)
    #c1b = -0.5 * np.prod(ovlp_oo_diag) * einsum('ijab,ja,cb,i,j, ic->', 
    #           co_t2, ovlp_ov, ovlp_vv, ovlp_oo_inv, ovlp_oo_inv, Fov)
    #temp1 = einsum('kc,kb,k->cab', ovlp_ov, ovlp_ov, ovlp_oo_inv) \
    #   - einsum('jc,jb,j->cab', ovlp_ov, ovlp_ov, ovlp_ov, ovlp_oo_inv)
    #temp2 = temp1.transpose(0,1,3,2)*(-1)
    #print(temp1)
    temp1 = get_S_ijpr_matrix(co_ovlp, ovlp_oo_diag, nocc, nvir)
    c1 =  0.5 * np.prod(ovlp_oo_diag) * einsum('ijab, ijbc, ja, i,j, ic', co_t2,
        temp1, ovlp_ov, ovlp_oo_inv, ovlp_oo_inv, Fov)
    
    if dbg:
        for i in range(nocc):
            for j in range(nocc):
                for a in range(nvir):
                    for b in range(nvir):
                        for c in range(nvir):
                            t2 = co_t2[i,j,a,b]
                            n1 = temp1[i,j,b,c]
                            n2 = ovlp_ov[j,a]
                            n = n1*n2*ovlp_oo_inv[i]*ovlp_oo_inv[j]*Fov[i,c]
                            n *= np.prod(ovlp_oo_diag)
                            if abs(t2*n) > 1e-5:
                                print(i,j,a+nocc,b+nocc,c+nocc, 
                                      ' %2.6f %2.6f %2.6f  %2.6f %2.6f'%( t2, n, t2*n, n1, n2) )
    
    # Case 2: i != k
    temp3 = einsum('ia,jb,i,j,kc,k, kc->ijab', 
      ovlp_ov, ovlp_ov, ovlp_oo_inv, ovlp_oo_inv, ovlp_ov, ovlp_oo_inv, Fov)
    c2 = 0.5 *  np.prod(ovlp_oo_diag) * einsum('ijab,ijab->', 
      co_t2, temp3)

    term15 = c1 + c2
    print('c1, c2 %.6f %.6f'% (c1, c2))
    return term15

@timing
def get_term2(co_eri, co_ovlp, ovlp_oo_diag, co_t2, ovlp_vo, ovlp_ov):
    nocc = co_eri.shape[0]
    nvir = co_eri.shape[2]
    co_eri_oovv = co_eri
    #co_eri_oovv = co_eri.transpose(0,2,1,3) - co_eri.transpose(0,2,3,1)
    ovlp_oo_inv = 1. / ovlp_oo_diag
    
    # Case 1: i = k, j = l
    S_ijpr_matrix = get_S_ijpr_matrix(co_ovlp, ovlp_oo_diag, nocc, nvir)
    sumc_ijad = einsum('ijcd,ijca->ijad', co_eri_oovv, S_ijpr_matrix)
    sumb_ijad = einsum('ijab,ijdb->ijad', co_t2, S_ijpr_matrix)
    sumadij = einsum('ijad,ijad,i,j->', sumc_ijad, sumb_ijad, ovlp_oo_inv, ovlp_oo_inv)
    term21 = sumadij * 0.25 * np.prod(ovlp_oo_diag)
    if False:
        for i in range(nocc):
            for j in range(nocc):
                for a in range(nvir):
                    #for b in range(nvir):
                    #    for c in range(nvir):
                    for d in range(nvir):
                        #eri = co_eri_oovv[i,j,c,d]
                        #t2 = co_t2[i,j,a,b]
                        n1 = sumc_ijad[i,j,a,d]
                        n2 = sumb_ijad[i,j,a,d]
                        n = n1*n2*ovlp_oo_inv[i]*ovlp_oo_inv[j]
                        n *= np.prod(ovlp_oo_diag)
                        if abs(n) > 1e-5:
                            print(i,j,a+nocc,#b+nocc,c+nocc, 
                            d+nocc,
                              ' %2.6f %2.6f %2.6f  '%( n, n1, n2) )
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
    print('term2.x %.6f, %.6f, %.6f' % (term21, term22, term23))
    e_term2 = term21 + term22 + term23
    return e_term2

def get_S_ijbc(co_ovlp, ovlp_oo_diag, nocc, nvir):
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
            S_ijpr_mat[i, j] = einsum('pk,k,kr->pr', co_ovlp_vo_ij, co_ovlp_oo_inv_ij, co_ovlp_ov_ij)
    
    return S_ijpr_mat

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
