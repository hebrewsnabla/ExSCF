'''
EMP2 and SUPT2
Tsuchimochi, T.; Van Voorhis, T. J Chem Phys 2014, 141, 164117
Tsuchimochi, T.; Ten-no, S. L. J Chem Theory Comput 2019, 15, 6688–6702
'''
import numpy as np
#import sympy as sym
#import scipy
from pyscf import gto, scf, mp, lib
from pyphf import suscf, util2, noci
from pyphf.suscf import count0, eig
#import os, sys
from functools import partial
import time
import copy

print = partial(print, flush=True)
einsum = partial(np.einsum, optimize=True)

def defm_Fock(mol, hcore_ortho, dm, X):
    gaa, gbb = defm_G(mol, dm, X)
    #F0 = hcore_ortho 
    #norb = gaa.shape[0]
    #Faa = g[:norb, :norb]
    #Fbb = g[norb:, norb:]
    F0 = gaa + hcore_ortho, gbb + hcore_ortho
    #F0_no = hcore_no + g_no
    return F0

def defm_G(mol, dm, X):
    #G = 
    #p = dm_no

    #p_ortho = einsum('ij,jk,lk->il', no, p, no)
    #print(p_ortho)
    p_ortho = dm

    #norb = int(p_ortho.shape[0]/2)

    paa = p_ortho[0] # ortho ao
    #print(pgaa)
    #pab = p_ortho[1]
    #pba = p_ortho[norb:, :norb]
    pbb = p_ortho[1]
    # X . P(g) . X^H
    paa_ao = einsum('ij,jk,lk->il', X, paa, X) # regular ao
    #print(pgaa_ao)
    #pab_ao = einsum('ij,jk,lk->il', X, pab, X)
    #pba_ao = einsum('ij,jk,lk->il', X, pba, X)
    pbb_ao = einsum('ij,jk,lk->il', X, pbb, X)
    #print(pbb)
    #print(pbb_ao)
    #Pgaabb_ao = Pgaa_ao + Pgbb_ao
    #print(Pgaabb_ao.shape)
    #nao = Pgaabb_ao.shape[-1]
    #ndm = len(Pgab_ao)
    vj,vk = scf.hf.get_jk(mol, [paa_ao, pbb_ao], hermi=0)
    #print(vj.shape)
    #print(vj, vk)
    Gaa_ao = vj[0] + vj[1] - vk[0]
    Gbb_ao = vj[0] + vj[1] - vk[1]
    #Ggbb_ao = scf.hf.get_jk(mol, Pgbb_ao, hermi=0)
    #print(ggaa_ao)
    #Gab_ao = scf.hf.get_jk(mol, pab_ao, hermi=0)[1] *(-1)
    #Gba_ao = scf.hf.get_jk(mol, pba_ao, hermi=0)[1] *(-1)
    #ggbb_ao = scf.uhf.get_veff(mol, [pgaa_ao, pgbb_ao], hermi=0)[1]
    # X^H . G(g) . X

    gaa = einsum('ji,jk,kl->il', X, Gaa_ao, X)  # ortho ao
    #print(ggaa)
    #gab = einsum('ji,jk,kl->il', X, Gab_ao, X) 
    #gba = einsum('ji,jk,kl->il', X, Gba_ao, X) 
    gbb = einsum('ji,jk,kl->il', X, Gbb_ao, X) 
    #g = util2.stack22(gaa, gab, gba, gbb)
    #g_no = einsum('ji,jk,kl->il', no, g, no)
    #print(gaa, gbb)

    return gaa, gbb

def t1(F0, epsl):
    ''' F_ia / (e_i - e_a) '''
    l = len(epsl)
    e = epsl - epsl.reshape((l,1))
    t1 = F0 / e
    return t1

def t2():
    return 0

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
        np.set_printoptions(precision=6, linewidth=160, suppress=True)
        suhf = self.suhf
        F0 = defm_Fock(self.mol, suhf.hcore_ortho, suhf.dm_ortho, suhf.X )
        #epsl = F0.diagonal()
        print('F0', F0)
        na, nb = self.nelec
        occ = na+nb
        norb = self.norb
        vir = norb*2 - occ
        nmo = self.norb
        nvira, nvirb = nmo-na, nmo-nb
        '''
        Foo = F0_no[:occ, :occ]
        Fvv = F0_no[occ:norb*2, occ:norb*2]
        eo, co = np.linalg.eigh(Foo)
        ev, cv = np.linalg.eigh(Fvv)
        #ea = F0[0].diagonal()
        #eb = F0[1].diagonal()
        #print(suhf.hcore_no)
        print('F0_no', F0_no)
        print(eo, ev)
        '''
        mo_ortho = suhf.mo_ortho
        print(mo_ortho)
        eao, eav, cao, cav = Diag_F(F0[0], na, norb)
        ebo, ebv, cbo, cbv = Diag_F(F0[1], nb, norb)
        ca = util2.stack22(cao, np.zeros((na, nvira)), np.zeros((nvira, na)), cav)
        cb = util2.stack22(cbo, np.zeros((nb, nvirb)), np.zeros((nvirb, nb)), cbv)
        print('ca,cb', ca, cb)
        #ca_sc = einsum('ij,kj->ik', mo_ortho[0], ca)
        #cb_sc = einsum('ij,kj->ik', mo_ortho[1], cb)
        ca_sc = ca.T
        cb_sc = cb.T
        #pa_sc = einsum('ij,kj->ik', ca_sc[:,:na], ca_sc[:,:na])
        #pb_sc = einsum('ij,kj->ik', cb_sc[:,:nb], cb_sc[:,:nb])
        print('ca,cb (sc)', ca_sc, cb_sc)
        ea = np.concatenate((eao, eav))
        eb = np.concatenate((ebo, ebv))
        print(ea,eb)
        Fsc = (einsum('ji,jk,kl->il', ca, F0[0], ca), 
              einsum('ji,jk,kl->il', cb, F0[1], cb))
        print(Fsc)
        #Fsc2 = defm_Fock(self.mol, suhf.hcore_ortho, [pa_sc, pb_sc], suhf.X)
        #print(Fsc2)

        #print(e0)
        #for i in range(self.occ):
        #    for a in range(self.occ, self.norb):
#        dm_no, _, no = find_NO(suhf, suhf.dm_ortho, na, nb)
#        Dg, Ng, Pg = suscf.get_Ng(suhf.grids, no, dm_no, na+nb)
        #print('D(g) (NO)\n', Dg[0])
        #print('N(g) (NO)\n', Ng[0])
        #print('P(g) (NO)\n', Pg[0])
#        C_no = suscf.get_xg(suhf, no, suhf.mo_occ, Ng)[-1]

#        xg, yg, ciS = get_xg(suhf, C_no, Dg, na+nb)
        #get_Mg2(suhf, Dg, na+nb)
        ca_sc_ortho = einsum('ij,jk->ik', mo_ortho[0], ca_sc)
        cb_sc_ortho = einsum('ij,jk->ik', mo_ortho[1], cb_sc)
        #mo = suhf.mo_ortho
        get_tau1(suhf, na, nb, nvira, nvirb, (ca_sc_ortho, cb_sc_ortho))
        ca_sc_ao = einsum('ij,jk,kl->il', suhf.X, mo_ortho[0], ca_sc)
        cb_sc_ao = einsum('ij,jk,kl->il', suhf.X, mo_ortho[1], cb_sc)
        eris = self.ao2mo([ca_sc_ao, cb_sc_ao])
        emp2 = 0.0
        
        eia_a = ea[:na,None] - ea[None,na:]
        eia_b = eb[:nb,None] - eb[None,nb:]
        for i in range(na):
            if eris.ovov.ndim == 4:
                eris_ovov = eris.ovov[i]
            else:
                eris_ovov = np.asarray(eris.ovov[i*nvira:(i+1)*nvira])
            print(eris_ovov.shape)    
            eris_ovov = eris_ovov.reshape(nvira,na,nvira).transpose(1,0,2)
            t2i = eris_ovov.conj()/lib.direct_sum('a+jb->jab', eia_a[i], eia_a)
            emp2 += einsum('jab,jab', t2i, eris_ovov) * .5
            emp2 -= einsum('jab,jba', t2i, eris_ovov) * .5

            if eris.ovOV.ndim == 4:
                eris_ovov = eris.ovOV[i]
            else:
                eris_ovov = np.asarray(eris.ovOV[i*nvira:(i+1)*nvira])
            eris_ovov = eris_ovov.reshape(nvira,nb,nvirb).transpose(1,0,2)
            t2i = eris_ovov.conj()/lib.direct_sum('a+jb->jab', eia_a[i], eia_b)
            emp2 += einsum('JaB,JaB', t2i, eris_ovov)
        for i in range(nb):
            if eris.OVOV.ndim == 4:
                eris_ovov = eris.OVOV[i]
            else:
                eris_ovov = np.asarray(eris.OVOV[i*nvirb:(i+1)*nvirb])
            eris_ovov = eris_ovov.reshape(nvirb,nb,nvirb).transpose(1,0,2)
            t2i = eris_ovov.conj()/lib.direct_sum('a+jb->jab', eia_b[i], eia_b)
            emp2 += einsum('jab,jab', t2i, eris_ovov) * .5
            emp2 -= einsum('jab,jba', t2i, eris_ovov) * .5
        print('emp2', emp2)
        
    def ao2mo(self, mo):
        fakemp = mp.UMP2(self.suhf.guesshf, mo_coeff=mo, mo_occ=self.suhf.guesshf.mo_occ)
        return fakemp.ao2mo()


def Diag_F(F, occ, norb):
    Foo = F[:occ, :occ]
    Fvv = F[occ:, occ:]
    eo, co = np.linalg.eigh(Foo)
    ev, cv = np.linalg.eigh(Fvv)
    return eo, ev, co, cv

def get_tau1(mf, na, nb, nvira, nvirb, mo):
    tau1_aa = np.zeros((na,nvira, na,nvira))
    for i in range(na):
        for a in range(na, na+nvira):
            for j in range(na):
                if j==i: continue
                for b in range(na, na+nvira):
                    if b==a: continue
                    ciS, ciH = noci.ci_cross(mf, [[i,j],[a,b]], [[],[]], mo)
                    print(i,a,j,b, ciS)
                    tau1_aa[i,a,j,b] = ciS
    tau1_bb = np.zeros((nb,nvirb, nb,nvirb))
    for i in range(nb):
        for a in range(nb, nb+nvirb):
            for j in range(nb):
                if j==i: continue
                for b in range(nb, nb+nvirb):
                    if b==a: continue
                    ciS, ciH = noci.ci_cross(mf, [[],[]], [[i,j],[a,b]], mo)
                    print(i,a,j,b, ciS)
                    tau1_bb[i,a,j,b] = ciS
    tau1_ab = np.zeros((na,nvira, nb,nvirb))
    for i in range(na):
        for a in range(na, na+nvira):
            for j in range(nb):
                #if j==i: continue
                for b in range(nb, nb+nvirb):
                    #if b==a: continue
                    ciS, ciH = noci.ci_cross(mf, [[i],[a]], [[j],[b]], mo)
                    print('cis ',i,a,j,b, ciS)
                    tau1_bb[i,a-na,j,b-nb] = ciS
                    

def find_NO(suhf, dm, na, nb, i=0, a=0):
    cut_no = suhf.cut_no
    #dm = dm*(-1)
    #np.set_printoptions(precision=16, linewidth=200, suppress=False)
    #print(dm)
    ev_a, v_a = eig(dm[0]*(-1))
    ev_b, v_b = eig(dm[1]*(-1))
    pa = count0(ev_a)
    pb = count0(ev_b)
    if suhf.debug:
        print('NO eigenvalue')
        print(ev_a, '\n', ev_b)
    v_a1 = v_a[:,:na]
    v_a2 = v_a[:,na:]
    v_b1 = v_b[:,:nb]
    v_b2 = v_b[:,nb:]
    #print(v_a1, v_a2, v_b1, v_b2)
    v_a1 = np.vstack((v_a1, np.zeros(v_a1.shape)))
    v_a2 = np.vstack((v_a2, np.zeros(v_a2.shape)))
    v_b1 = np.vstack((np.zeros(v_b1.shape), v_b1))
    v_b2 = np.vstack((np.zeros(v_b2.shape), v_b2))
    #print(v_a1, v_a2, v_b1, v_b2)

    v = np.hstack((v_a1, v_b1, v_a2, v_b2))
    v_flip = copy.deepcopy(v)
    v_flip[:,i] = v[:,a]
    v_flip[:,a] = v[:,i]
    if cut_no:
        v = np.hstack((v_a1, v_b1, v_a2, v_b2))[:,:pa+pb]
    #v = np.hstack((v, np.zeros((v.shape[0], v.shape[0]-pa-pb))))
    if suhf.debug:
        print('NO vec')
        print(v_flip)
    dm_expd = np.hstack(
        (np.vstack((dm[0], np.zeros(dm[0].shape))), 
        np.vstack((np.zeros(dm[1].shape), dm[1])))
        )
    #print(dm_expd)
    dm_no = einsum('ji,jk,kl->il', v, dm_expd, v)
    dm_no_flip = einsum('ji,jk,kl->il', v_flip, dm_expd, v_flip)
    if suhf.debug:
        print('dm(NO)')
        print(dm_no_flip)
    #np.set_printoptions(precision=6, linewidth=160, suppress=True)
    return dm_no_flip, dm_expd, v_flip

def get_xg(suhf, C_no, Dg, occ):
    Cocc = C_no[:,:occ]
    Mg = []
    xg = []
    for dg in Dg:
        mg = einsum('ji,jk,kl->il', Cocc, dg, Cocc)
        x = np.linalg.det(mg)
        Mg.append(mg)
        xg.append(x)
    print('xg', xg)

    ciS = suhf.integr_beta(np.array(xg))
    yg = xg / ciS
    return xg, yg, ciS
    

