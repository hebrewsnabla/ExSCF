import numpy as np
import sympy as sym
import scipy
from pyscf import gto, scf
from pyscf.scf import chkfile
from pyscf.lib import temporary_env
from pyphf import sudm, util2
import os, sys
from functools import partial
import time

print = partial(print, flush=True)
einsum = partial(np.einsum, optimize=True)


def get_beta(n):
    # Gauss-Legendre quadrature
    PI = np.pi
    b0,b1 = 0.0, PI
    m = (n+1)//2
    bm = 0.5*(b1+b0)
    bl = 0.5*(b1-b0)
    grid = np.zeros(n)
    weight = np.zeros(n)
    for i in range(m):
        z = np.cos(PI*(i+1-0.25) / (n + 0.5))
        while(True):
            p1, p2 = 1.0, 0.0
            for j in range(n):
                p3 = p2
                p2 = p1
                p1 = ((2*j + 1)*z*p2 - (j)*p3)/(j+1)
            pp = n*(z*p1 - p2)/(z*z-1)
            z = z - p1/pp
            if (abs(p1/pp) < 3e-14): 
                break
        grid[i] = bm - bl*z
        grid[n-1-i] = bm + bl*z
        weight[i] = 2*bl / ((1 - z*z)*pp*pp)
        weight[n-1-i] = weight[i]
        # Note: We did not perform weight *= sin(beta) here. That's left in get_xg().
    return grid, weight

def WignerSmall(j,m):
    #j = sym.symbols('j')
    #m = sym.symbols('m')
    j = sym.sympify(j)/2
    m = sym.sympify(m)/2
    s = sym.symbols('s')
    beta = sym.symbols('beta')
    if j==0:
        d_simp = 1
    else:
        d = (-1)**s * sym.binomial(j+m,s) * sym.binomial(j-m,s) * sym.cos(beta/2)**(2*j-2*s) * sym.sin(beta/2)**(2*s)
        upper = min(j-m,j+m)
        d = sym.Sum(d, (s, 0, upper)).doit()
        d_simp = sym.trigsimp(d)
    f = sym.lambdify(beta, d_simp, 'numpy')
    return d_simp, f

def count0(vals):
    non0 = 0
    for i in vals:
        if abs(i)>1e-10:
            non0 += 1
    return non0

def eig(A):
    return np.linalg.eigh(A)
    #return scipy.linalg.eigh(A)

def find_NO(suhf, dm, na, nb):
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
    if cut_no:
        v = np.hstack((v_a1, v_b1, v_a2, v_b2))[:,:pa+pb]
    #v = np.hstack((v, np.zeros((v.shape[0], v.shape[0]-pa-pb))))
    if suhf.debug:
        print('NO vec')
        print(v)
    dm_expd = np.hstack(
        (np.vstack((dm[0], np.zeros(dm[0].shape))), 
        np.vstack((np.zeros(dm[1].shape), dm[1])))
        )
    #print(dm_expd)
    dm_no = einsum('ji,jk,kl->il', v, dm_expd, v)
    if suhf.debug:
        print('dm(NO)')
        print(dm_no)
    #np.set_printoptions(precision=6, linewidth=160, suppress=True)
    return dm_no, dm_expd, v

def get_Ng(grids, no, dm, occ):
    Dg = []
    Ng = []
    Pg = []
    #print('N(g)')
    for beta in grids:
        norb = int(len(no)/2)
        dg1 = np.eye(norb) * np.cos(beta/2)
        dg2 = np.eye(norb) * (-np.sin(beta/2))
        dg = np.vstack((
            np.hstack((dg1, -dg2)),
            np.hstack((dg2, dg1))
        ))
        #print(dg.shape)
        #print(dg)
        dg_no = einsum('ji,jk,kl->il', no, dg, no)
        #det_dg = np.linalg.det(dg_no)
        #print(det_dg)
        #print(dg_no)
        Dg.append(dg_no)
        tg = einsum('ij,jk,kl->il', dm[:occ,:], dg_no, dm[:,:occ])
        ng = np.linalg.inv(tg)
        #print(ng)
        Ng.append(ng)
        pg = einsum('ij,jk,kl,lm->im', dg_no, dm[:,:occ], ng, dm[:occ,:])
        #print(pg)
        Pg.append(pg)
    return Dg, Ng, Pg

def get_JKg(mol, Pg, no, X):
    Jg = []
    Kg = []
    Pg_ortho = []
    for pg in Pg:
        pg_ortho = einsum('ij,jk,lk->il', no, pg, no)
        #print(pg_ortho)
        Pg_ortho.append(pg_ortho)
    norb = int(Pg_ortho[0].shape[0]/2)
    Pgaa_ao = []
    Pgab_ao = []
    Pgba_ao = []
    Pgbb_ao = []
    for pg in Pg_ortho:
        pgaa = pg[:norb, :norb] # ortho ao
        #print(pgaa)
        pgab = pg[:norb, norb:]
        pgba = pg[norb:, :norb]
        pgbb = pg[norb:, norb:]
        # X . P(g) . X^H
        pgaa_ao = einsum('ij,jk,lk->il', X, pgaa, X) # regular ao
        #print(pgaa_ao)
        pgab_ao = einsum('ij,jk,lk->il', X, pgab, X)
        pgba_ao = einsum('ij,jk,lk->il', X, pgba, X)
        pgbb_ao = einsum('ij,jk,lk->il', X, pgbb, X)
        Pgaa_ao.append(pgaa_ao)
        Pgbb_ao.append(pgbb_ao)
        Pgab_ao.append(pgab_ao)
        Pgba_ao.append(pgba_ao)
    Pgaabb_ao = Pgaa_ao + Pgbb_ao
    #print(Pgaabb_ao.shape)
    ndm = len(Pgab_ao)
    vj,vk = scf.hf.get_jk(mol, Pgaabb_ao, hermi=0)
    #print(vj.shape)
    Jgaa_ao = vj[:ndm] + vj[ndm:] 
    Kgaa_ao = - vk[:ndm]
    Jgbb_ao = Jgaa_ao
    Kgbb_ao = - vk[ndm:]
    Kgab_ao = scf.hf.get_jk(mol, Pgab_ao, hermi=0)[1] *(-1)
    Kgba_ao = scf.hf.get_jk(mol, Pgba_ao, hermi=0)[1] *(-1)
    for i in range(len(Jgaa_ao)):
        jgaa = einsum('ji,jk,kl->il', X, Jgaa_ao[i], X)  # ortho ao
        kgaa = einsum('ji,jk,kl->il', X, Kgaa_ao[i], X)  # ortho ao
        #print(ggaa)
        kgab = einsum('ji,jk,kl->il', X, Kgab_ao[i], X) 
        kgba = einsum('ji,jk,kl->il', X, Kgba_ao[i], X) 
        jgbb = einsum('ji,jk,kl->il', X, Jgbb_ao[i], X) 
        kgbb = einsum('ji,jk,kl->il', X, Kgbb_ao[i], X) 
        jg = util2.stack22(jgaa, np.zeros(kgab.shape), 
                           np.zeros(kgba.shape), jgbb)
        kg = util2.stack22(kgaa, kgab, kgba, kgbb)
        jg_no = einsum('ji,jk,kl->il', no, jg, no)
        kg_no = einsum('ji,jk,kl->il', no, kg, no)
        Jg.append(jg_no)
        Kg.append(kg_no)
    return Jg, Kg, Pg_ortho

def get_jk(mol, dm, hermi=1, opt=None):
    return scf.hf.get_jk(mol, dm, hermi, opt)

def get_k(mol, dm, hermi=1, opt=None):
    with temporary_env(opt, prescreen='CVHFnrs8_vk_prescreen'):
        vk = scf.hf.get_jk(mol, dm, hermi, opt, with_j=False)[1]
    return vk

def get_Gg(mol, Pg, no, X, dm_last=None, Ggao_last=None, opt=None):
    Gg = []
    Pg_ortho = []
    for pg in Pg:
        pg_ortho = einsum('ij,jk,lk->il', no, pg, no)
        #print(pg_ortho)
        Pg_ortho.append(pg_ortho)
    norb = int(Pg_ortho[0].shape[0]/2)
    Pgaa_ao = []
    Pgab_ao = []
    Pgba_ao = []
    Pgbb_ao = []
    for pg in Pg_ortho:
        pgaa = pg[:norb, :norb] # ortho ao
        #print(pgaa)
        pgab = pg[:norb, norb:]
        pgba = pg[norb:, :norb]
        pgbb = pg[norb:, norb:]
        # X . P(g) . X^H
        pgaa_ao = einsum('ij,jk,lk->il', X, pgaa, X) # regular ao
        #print(pgaa_ao)
        pgab_ao = einsum('ij,jk,lk->il', X, pgab, X)
        pgba_ao = einsum('ij,jk,lk->il', X, pgba, X)
        pgbb_ao = einsum('ij,jk,lk->il', X, pgbb, X)
        Pgaa_ao.append(pgaa_ao)
        Pgbb_ao.append(pgbb_ao)
        Pgab_ao.append(pgab_ao)
        Pgba_ao.append(pgba_ao)
    Pgaabb_ao = Pgaa_ao + Pgbb_ao
    Pgao = [Pgaabb_ao, Pgab_ao, Pgba_ao]
    if dm_last is not None:
        old_Pgaabb_ao, old_Pgab_ao, old_Pgba_ao = dm_last
        old_vj, old_vk, old_Ggab_ao, old_Ggba_ao = Ggao_last
    #print(Pgaabb_ao.shape)
    #nao = Pgaabb_ao.shape[-1]
    if dm_last is None:
        vj,vk = get_jk(mol, Pgaabb_ao, hermi=0, opt=opt)
        Ggab_ao = get_k(mol, Pgab_ao, hermi=0, opt=opt)
        Ggba_ao = get_k(mol, Pgba_ao, hermi=0, opt=opt)
    else:
        d_aabb = util2.dmlist(Pgaabb_ao, old_Pgaabb_ao)
        d_ab = util2.dmlist(Pgab_ao, old_Pgab_ao)
        d_ba = util2.dmlist(Pgba_ao, old_Pgba_ao)
        vj,vk = get_jk(mol, d_aabb, hermi=0, opt=opt) 
        #vj = util2.dmlist(vj, old_vj, 1)
        #vk = util2.dmlist(vk, old_vk, 1)
        vj += old_vj
        vk += old_vk
        Ggab_ao = get_jk(mol, d_ab, hermi=0, opt=opt)[1] + old_Ggab_ao
        #Ggab_ao = util2.dmlist(Ggab_ao, old_Ggab_ao, 1)
        Ggba_ao = get_jk(mol, d_ba, hermi=0, opt=opt)[1] + old_Ggba_ao
        #Ggba_ao = util2.dmlist(Ggba_ao, old_Ggba_ao, 1)
    ndm = len(Pgab_ao)
    #print(vj.shape)
    #print(vj)
    Ggaa_ao = vj[:ndm] + vj[ndm:] - vk[:ndm]
    Ggbb_ao = vj[:ndm] + vj[ndm:] - vk[ndm:]
    #Ggbb_ao = scf.hf.get_jk(mol, Pgbb_ao, hermi=0)
    #print(ggaa_ao)
    Ggao = [vj,vk,Ggab_ao, Ggba_ao]
    #ggbb_ao = scf.uhf.get_veff(mol, [pgaa_ao, pgbb_ao], hermi=0)[1]
        # X^H . G(g) . X
    Ggab_ao *= -1
    Ggba_ao *= -1
    for i,ggab_ao in enumerate(Ggab_ao):
        #ggab_ao = Ggab_ao[i]
        ggba_ao = Ggba_ao[i]
        ggaa_ao = Ggaa_ao[i]
        ggbb_ao = Ggbb_ao[i]
        ggaa = einsum('ji,jk,kl->il', X, ggaa_ao, X)  # ortho ao
        #print(ggaa)
        ggab = einsum('ji,jk,kl->il', X, ggab_ao, X) 
        ggba = einsum('ji,jk,kl->il', X, ggba_ao, X) 
        ggbb = einsum('ji,jk,kl->il', X, ggbb_ao, X) 
        gg = util2.stack22(ggaa, ggab, ggba, ggbb)
        gg_no = einsum('ji,jk,kl->il', no, gg, no)
        Gg.append(gg_no)
    return Gg, Pg_ortho, Pgao, Ggao


def get_xg(suhf, no, na, nb, Ng):
    C_a, C_b = suhf.mo_ortho
    C_a1 = C_a[:,:na]
    C_a2 = C_a[:,na:]
    C_b1 = C_b[:,:nb]
    C_b2 = C_b[:,nb:]
    C_org = np.hstack((
        np.vstack((C_a1, np.zeros(C_a1.shape))),
        np.vstack((np.zeros(C_b1.shape), C_b1)),
        np.vstack((C_a2, np.zeros(C_a2.shape))),
        np.vstack((np.zeros(C_b2.shape), C_b2))
    ))
    #print(C_org)
    C_no = einsum('ji,jk->ik', no, C_org)
    if suhf.debug:
        print('C(NO)')
        print(C_no)
    occ = na+nb
    C_oo = C_no[:occ, :occ]
    detC = np.linalg.det(C_oo)
    print('detC', detC)
    detNg = []
    for ng in Ng:
        detng = np.linalg.det(ng)
        #print(detng)
        detNg.append(detng)
    detNg = np.array(detNg)
    print('detNg', detNg)
    xg = 1.0 / (detC * detNg * detC)
    print('xg', xg)
    #sinbeta = np.sin(grids)
    #ciS = einsum('i,i,i,i->', weights, sinbeta, xg, d)
    ciS = suhf.integr_beta(xg)
    #print(weights, xg, d)
    print('ciS', ciS)
    yg = xg / ciS
    return xg, yg, ciS, C_no

def integr_beta(q, d, grids, weights, fac='normal', xg=None, ciS=None):
    sinbeta = np.sin(grids)
    if fac=='xg':
        weights = weights * xg / ciS
        #print(xg, ciS, xg/ciS**2)
        #weights *= (xg / ciS**2)
        #print(weights)
    elif fac=='ci':
        weights = weights / ciS

    if q.ndim==1:
        int_q = einsum('i,i,i,i->', weights, sinbeta, q, d)
    elif q.ndim==2: 
        int_q = einsum('i,i,ij,i->j', weights, sinbeta, q, d)
    elif q.ndim==3:
        int_q = einsum('i,i,ijk,i->jk', weights, sinbeta, q, d)
    else:
        raise(ValueError, 'q must have dimension 1, 2, or 3')
    return int_q


def get_H(suhf, hcore_ortho, no, Pg, Gg, xg):
    #print(hcore_ortho)
    #hcore_ortho = suhf.hcore_ortho
    hcore_ortho = np.vstack((
        np.hstack((hcore_ortho, np.zeros(hcore_ortho.shape))),
        np.hstack((np.zeros(hcore_ortho.shape), hcore_ortho))
    ))
    hcore_no = einsum('ji,jk,kl->il', no, hcore_ortho, no)
    suhf.hcore_no = hcore_no
    if suhf.debug2:
        print(hcore_no)
    trHg = np.zeros(len(Pg))
    for i, pg in enumerate(Pg):
        H = np.trace(np.dot(hcore_no, pg)) + 0.5 * np.trace(np.dot(Gg[i], pg))
        #H = H * xg[i]
        trHg[i] = H
        #print(i, H*xg[i])
    #sinbeta = np.sin(grids)
    #ciH = einsum('i,i,i,i->', weights, sinbeta, Hg, d)
    ciH = suhf.integr_beta(trHg*xg)
    print('ciH', ciH)
    return trHg, ciH

def get_EX(suhf, no, Pg, Kg, xg):
    trXg = np.zeros(len(Pg))
    for i, pg in enumerate(Pg):
        trXg[i] = 0.5 * np.trace(np.dot(Kg[i], pg))
        #H = H * xg[i]
    #ciH = suhf.integr_beta(trHg*xg)
    #print('ciH', ciH)
    X = suhf.integr_beta(trXg, fac='xg')
    return trXg, X

def get_S2(suhf, Pg_ortho):
    norb = int(Pg_ortho[0].shape[0]/2)
    S2g = np.zeros(len(Pg_ortho))
    for i,pg in enumerate(Pg_ortho):
        pgaa = pg[:norb, :norb] # ortho ao
        #print(pgaa)
        pgab = pg[:norb, norb:]
        pgba = pg[norb:, :norb]
        pgbb = pg[norb:, norb:]
        Pc = 0.5 * (pgaa + pgbb)
        Mz = 0.5 * (pgaa - pgbb)
        My = -0.5j * (pgab - pgba)
        #print(My)
        Mx = 0.5 * (pgab + pgba)
        #print(Pc)
        trPc, trMx, trMy, trMz = list(map(np.trace, [Pc, Mx, My, Mz]))
        if suhf.debug: print(trPc, trMx, trMy, trMz)
        Pc2 = np.dot(Pc, Pc)
        Mz2 = np.dot(Mz, Mz)
        My2 = np.dot(My, My).real
        Mx2 = np.dot(Mx, Mx)
        trPc2, trMx2, trMy2, trMz2 = list(map(np.trace, [Pc2, Mx2, My2, Mz2]))
        if suhf.debug: print(trPc2, trMx2, trMy2, trMz2)
        S2g[i] = trMx**2 + (trMy**2).real + trMz**2 + 0.5*(trMx2 + trMy2 + trMz2) + 1.5*(trPc - trPc2)
    S2 = suhf.integr_beta(S2g, fac='xg')
    print('S2 = %.6f'% S2)
    return S2

def get_Yg(suhf, Dg, Ng, dm_no, occ):
    norb = len(dm_no[0])
    vir = norb - occ
    Xg = []
    for i, ng in enumerate(Ng):
        Xgi1 = einsum('ij,jk,kl->il', Dg[i], dm_no[:,:occ], Ng[i]) 
        Xgi2 = einsum('ij,jk,kl->il', Ng[i], dm_no[:occ,:], Dg[i])
        #print(Xgi1)
        #print(Xgi2)
        Xgi = np.hstack((Xgi1, np.zeros((norb, vir)))) \
            + np.vstack((Xgi2, np.zeros((vir, norb)))) 
        #print(Xgi)
        Xg.append(Xgi)
    if suhf.debug:
        print('X(g)')
        print(Xg[0])
    Xg = np.array(Xg)
    Xg_int = suhf.integr_beta(Xg, fac='xg')
    if suhf.debug:
        print('X(g) int')
        print(Xg_int)
    Yg = Xg - Xg_int
    #print(Yg.shape)
    return Xg, Xg_int, Yg

def get_Feff(suhf, trHg, Gg, Ng, Pg, Dg, occ, Yg, Xg, F_ortho):
    hcore_no = suhf.hcore_no
    dm_ortho = suhf.dm_ortho
    dm_no = suhf.dm_no
    no = suhf.no
    norb = len(Pg[0])
    vir = norb - occ
    Feff_g = [] 
    Feff0 = []
    for i, pg in enumerate(Pg):
        feff0 = Yg[i] * trHg[i]
        #print(feff0)
        fg = hcore_no + Gg[i]
        feff1 = einsum('ij,jk,kl,lm,mn->in', Ng[i], dm_no[:occ,:], fg, np.eye(norb)-pg, Dg[i])
        #print(feff1)
        feff1 = np.vstack((feff1, np.zeros((vir, norb))))
        feff2 = einsum('ij,jk,kl,lm,mn->in', np.eye(norb)-pg, fg, Dg[i], dm_no[:, :occ], Ng[i])
        #print(feff2)
        feff2 = np.hstack((feff2, np.zeros((norb, vir))))
        feff = feff0 + feff1 + feff2
        Feff_g.append(feff)
        #Feff0.append(feff0)
        #print(feff)
    Feff_g = np.array(Feff_g)
    #Feff0 = suhf.integr_beta(np.array(Feff0), fac='xg')
    #print(Feff0)
    Feff = suhf.integr_beta(Feff_g, fac='xg')
    if suhf.debug2: print('Feff\n', Feff)
    Feff_ortho = einsum('ij,jk,lk->il', no, Feff, no)
    Feff_ortho = 0.5 * (Feff_ortho + Feff_ortho.T)
    if suhf.debug2: print('Feff (ortho)\n', Feff_ortho)
    H = suhf.integr_beta(trHg, fac='xg')
    print('Hsp + Hph = ', H)

    lenf = len(F_ortho[0])
    F_ortho = np.array(F_ortho)
    Feff_ortho = np.array([Feff_ortho[:lenf, :lenf], Feff_ortho[lenf:,lenf:]])
    Q = np.eye(len(dm_ortho[0])) - dm_ortho
    Foo = einsum('tij,tjk,tkl->til', dm_ortho, F_ortho, dm_ortho)
    Fvv = einsum('tij,tjk,tkl->til', Q, F_ortho, Q)
    Fov = einsum('tij,tjk,tkl->til', dm_ortho, Feff_ortho, Q)
    Fvo = einsum('tij,tjk,tkl->til', Q, Feff_ortho, dm_ortho)
    #F_mod_ortho_a = np.vstack((np.hstack((Foo[0], Fov[0])), np.hstack((Fvo[0], Fvv[0]))))
    #F_mod_ortho_b = np.vstack((np.hstack((Foo[1], Fov[1])), np.hstack((Fvo[1], Fvv[1]))))
    F_mod_ortho_a = Foo[0] + Fov[0] + Fvo[0] + Fvv[0]
    F_mod_ortho_b = Foo[1] + Fov[1] + Fvo[1] + Fvv[1]
    F_mod_ortho_a = 0.5*(F_mod_ortho_a + F_mod_ortho_a.T)
    F_mod_ortho_b = 0.5*(F_mod_ortho_b + F_mod_ortho_b.T)
    F_mod_ortho = np.vstack((
        np.hstack((F_mod_ortho_a, np.zeros(F_mod_ortho_a.shape))),
        np.hstack((np.zeros(F_mod_ortho_b.shape), F_mod_ortho_b))
    ))
    if suhf.debug:
        print('Feff (mod,ortho)')
        print(F_mod_ortho)
    return Feff_ortho, H, F_mod_ortho

def Diag_Feff(F, na, nb):
    e_a, v_a = np.linalg.eigh(F[0])
    e_b, v_b = np.linalg.eigh(F[1])
    P_a = einsum('ij,kj->ik', v_a[:,:na], v_a[:,:na])
    P_b = einsum('ij,kj->ik', v_b[:,:nb], v_b[:,:nb])
    return [e_a, e_b], [v_a, v_b], np.array([P_a, P_b])

class SUHF():
    '''
    Attributes:
        guesshf: UHF object
        X: 
            transformation matrix (regular AO -> orthonormal AO)
            It's asymmetric, as X in G09
        S, Sz:
        nbeta: number of grids for beta
        grids: integration grids for beta
        weights:
        d_expr: SymPy expr for Wigner d 
        d_func: NumPy function for Wigner d
        d     : array
            current values of Wigner d, on given grids

        E_suhf:
        natocc: SUHF natural orbital occupation number
        natorb: SUHF natural orbital (regular basis) 
    '''

    def __init__(self, guesshf=None):
        self.guesshf = guesshf

        self.cut_no = False
        self.use_no = True
        self.verbose = 4
        #self.debug = False
        self.output = None # since 0.3.1, define output is important, it decides self.chkfile
        self.restart = False
        self.conv_tol = 1e-7 # For RMSD
        self.max_cycle = 70
        self.diis_on = True
        self.diis_start_cyc = None
        self.level_shift = None
        self.makedm = True
        self.tofch = False
        self.oldfch = None

        self.built = False

        self.E_suhf = None

    def build(self):
        print('\n******** %s ********' % self.__class__)
        print('max_cycle = %d' % self.max_cycle)
        self.debug = False
        self.debug2 = False
        if self.verbose <= 4:
            print('verbose: %d, normal' % self.verbose)
        elif self.verbose <= 6:
            self.debug = True
            print('verbose: %d, debug' % self.verbose)
        else:
            self.debug = True
            self.debug2 = True
            print('verbose: %d, debug2' % self.verbose)
        if self.output is None:
            self.output = str(os.getpid())
        if self.guesshf is not None:
            hf = self.guesshf
            self.mol = hf.mol
            self.chkfile0 = self.output + '_ges.pchk'
            #chkfile.dump_scf(hf.mol, self.chkfile0, hf.e_tot, hf.mo_energy,
            #                     hf.mo_coeff, hf.mo_occ)
            print('chkfile0: %s # the file store hf for guess' % self.chkfile0)
        elif self.chkfile0 is not None:
            pass ## todo
        elif self.chkfile is not None:
            self.mol, suinfo = util2.load(self.chkfile)
        else:
            guess = ''' 
            guesshf: a UHF object
            chkfile0: a PySCF chkfile storing a UHF result
            chkfile: SUHF chkfile
            '''
            raise AttributeError('You must provide one of below as a guess:' + guess)
        self.chkfile = self.output + '_su.pchk'
        print('chkfile: %s # the file store suhf info' % self.chkfile)
        #self.chkfile2 = self.output + '_no.pchk'
        #print('chkfile2: %s # the file store suhf NO' % self.chkfile2)
        print('conv_tol: %g' % self.conv_tol)
        if self.conv_tol > 1e-5:
            print('Warning: conv_tol too large')

        if self.diis_on:
            #assert issubclass(mf.DIIS, lib.diis.DIIS)
            #DIIS = lib.diis.SCF_DIIS
            self.diis_space = 8
            if self.diis_start_cyc is None:
                self.diis_start_cyc = 10
            self.diis_file = None
            #mf_diis.rollback = mf.diis_space_rollback
            self.diis = scf.diis.CDIIS()
            #self.max_cycle = self.diis_start_cyc + 30
            print('DIIS: %s' % self.diis.__class__)
            print('diis_start_cyc = %d' % self.diis_start_cyc)
        if self.level_shift is not None:
            shift = self.level_shift
            print('level shift: %.3f a.u.' % shift)
        #if self.output is not None:
        #    os.system("echo '' > %s" % self.output)
        #    sys.stdout = open(self.output, 'a')
        S = scf.hf.get_ovlp(self.mol)
        Ca, Cb = self.guesshf.mo_coeff
        #S_sqrt = scipy.linalg.sqrtm(S)
        if self.debug:
            print('S')
            print(S)
        #Se, Svec = scipy.linalg.eigh(S)
        #Se_msq = Se.real**(-0.5)
        #print(Se_msq)
        #S_msq2 = einsum('ji,j,jk->ik',Svec,Se_msq, Svec)
        #S_msq = scipy.linalg.fractional_matrix_power(S, -0.5)
        
        Su, Ss, Sv = scipy.linalg.svd(S)
        #print('SVD')
        #print(Su,Ss,Sv)
        X = einsum('ij,j->ij', Su, Ss**(-0.5))
        #Xd = einsum('i,ij->ij', Ss**(-0.5), Sv)
        if self.debug:
            print('SVD: S^(-1/2)')
            print(X)
        self.X = X
        #print(Xd)
        XS = np.dot(X.T,S)
        if self.debug:
            print('S^(1/2)')
            print(XS)
        self.XS = XS
        
        Ca_ortho = np.dot(XS, Ca)
        Cb_ortho = np.dot(XS, Cb)
        if self.debug:
            print('C (ortho)\n', Ca_ortho, '\n', Cb_ortho)
        self.mo_ortho = Ca_ortho, Cb_ortho
        self.dm_ortho = scf.uhf.make_rdm1(self.mo_ortho, self.guesshf.mo_occ)
        if self.debug:
            print('density matrix (ortho)')
            print(self.dm_ortho)
        self.norb = len(self.dm_ortho[0])

        self.nbeta = 8
        self.grids, self.weights = get_beta(self.nbeta)
        print('grids: ', self.grids, '\nweights: ', self.weights)
        spin = self.mol.spin
        na, nb = self.guesshf.nelec
        self.nelec = na, nb
        sz = (na-nb)/2
        print('S = %.1f, Sz = %.1f' % (spin/2, sz))
        self.S, self.Sz = spin/2, sz
        Wignerd_expr, Wignerd = WignerSmall(int(spin), int(2*sz))
        print('Wigner small d: ', Wignerd_expr)
        d = np.zeros(self.nbeta)
        for g in range(self.nbeta):
            d[g] = Wignerd(self.grids[g])
        print('value :', d)
        self.d_expr, self.d_func, self.d = Wignerd_expr, Wignerd, d
        self.E_suhf = None
        self.energy_nuc = self.mol.energy_nuc()
        hcore = hf.get_hcore()
        self.hcore_ortho = einsum('ji,jk,kl->il', X, hcore, X)
        if self.debug:
            print('hcore (ortho)\n', self.hcore_ortho)
        self.vhfopt = hf.init_direct_scf()

        self.built = True

    def integr_beta(self, q, fac='normal'):
        if fac=='xg':
            return integr_beta(q, self.d, self.grids, self.weights, 'xg', self.xg, self.ciS)
        elif fac=='ci':
            return integr_beta(q, self.d, self.grids, self.weights, 'ci', None, self.ciS)
        else:
            return integr_beta(q, self.d, self.grids, self.weights)

    
    def kernel(self):
        t_start = time.time()
        if not self.built:
            self.build()
        np.set_printoptions(precision=6, linewidth=160, suppress=True)
        if self.debug:
            np.set_printoptions(precision=10, linewidth=200, suppress=False)
        X = self.X
        na, nb = self.nelec
        norb = self.norb
        #mf = self.guesshf
        
        thresh = self.conv_tol
        max_cycle = self.max_cycle
        cyc = 1
        conv = False
        
        Pgao = None
        #print(vhfopt)
        t_pre = time.time() 
        print('time for Preparation before cyc: %.3f' % (t_pre-t_start))
        while(not conv):
            print('**** Start Cycle %d ****' % cyc)
            old_suhf = self.E_suhf
            old_dm = self.dm_ortho
            #if Pgao is not None:
            #    old_Pgao = Pgao
            #    old_Ggao = Ggao
            #else:
            old_Pgao = old_Ggao = None
            #if cyc==0:
            #    veff = mf.get_veff(dm = dm)
            #else:
            t01 = time.time()
            dm_reg = einsum('ij,tjk,lk->til', X, self.dm_ortho, X)
            veff = scf.uhf.get_veff(self.mol, dm_reg, vhfopt=self.vhfopt)
            veff_ortho = einsum('ji,tjk,kl->til', X, veff, X)
            if self.debug:
                print('dm (ortho)')
                print(self.dm_ortho)
            #Fa, Fb = hcore + veff
            #Fa_ortho = einsum('ji,jk,kl->il', X, Fa, X)
            #Fb_ortho = einsum('ji,jk,kl->il', X, Fb, X)
            Fa_ortho, Fb_ortho = self.hcore_ortho + veff_ortho
            F_ortho = Fa_ortho, Fb_ortho
            if self.debug:
                print('Fock (ortho)\n', F_ortho)
                e_uhf, e_uhf_coul = scf.uhf.energy_elec(self.guesshf, self.dm_ortho, self.hcore_ortho, veff_ortho)
                print('E(UHF) = %12.6f' % e_uhf)
            if self.use_no:
                dm_no, dm_expanded, no = find_NO(self, self.dm_ortho, na, nb)
            self.dm_no = dm_no
            self.no = no
            Dg, Ng, Pg = get_Ng(self.grids, self.no, self.dm_no, na+nb)
            self.Pg = Pg
            if self.debug:
                print('D(g) (NO)\n', Dg[0])
                print('N(g) (NO)\n', Ng[0])
                print('P(g) (NO)\n', Pg[0])
            t05 = time.time()
            print('time for NO, Ng: %.3f' % (t05-t01))
            Gg, Pg_ortho, Pgao, Ggao = get_Gg(self.mol, Pg, self.no, X, dm_last=old_Pgao, Ggao_last=old_Ggao, opt=self.vhfopt)
            if self.debug:
                print('Pg_ortho\n', Pg_ortho[0])
                print('G(g) (NO)\n' , Gg[0])
            t06 = time.time()
            print('time for Gg: %.3f' % (t06-t05))
            xg, yg, ciS, C_no = get_xg(self, self.no, na, nb, Ng)
            self.xg, self.ciS = xg, ciS
            #yg, ciS = util.get_yg(self, xg)
            trHg, ciH = get_H(self, self.hcore_ortho, self.no, Pg, Gg, xg)
            S2 = get_S2(self, Pg_ortho)
            Xg, Xg_int, Yg = get_Yg(self, Dg, Ng, self.dm_no, na+nb)
            Feff_ortho, H_suhf, F_mod_ortho = get_Feff(self, trHg, Gg, Ng, Pg, Dg, na+nb, Yg, Xg, F_ortho)
            E_suhf = self.energy_nuc + H_suhf
            self.E_suhf = E_suhf
            #print('E(SUHF) = %15.8f' % E_suhf)

            Faa = F_mod_ortho[:norb, :norb]
            Fbb = F_mod_ortho[norb:, norb:]
            F_mod_ortho = np.array([Faa,Fbb])
            if self.diis_on and cyc >= self.diis_start_cyc:
                s1e = np.eye(norb)
                F_mod_ortho = self.diis.update(s1e, self.dm_ortho, F_mod_ortho)
                print('F(mod,ortho) updated with CDIIS')
                if self.debug:
                    print(F_mod_ortho)
            if self.level_shift is not None:
                shift = self.level_shift
                s1e = np.eye(norb)
                print('level shift: %.3f a.u.' % shift)
                F_mod_ortho = lev_shift(s1e, self.dm_ortho, F_mod_ortho, shift)
            mo_e, mo_ortho, dm_ortho = Diag_Feff(F_mod_ortho, na, nb)
            util2.dump_moe(mo_e, na, nb)
            if self.debug:
                print('e_a, e_b\n', mo_e[0], '\n', mo_e[1])
                print('v_a, v_b\n', mo_ortho[0], '\n', mo_ortho[1])
                print('P_a, P_b\n', dm_ortho[0],'\n', dm_ortho[1])
            self.dm_ortho = dm_ortho
            self.mo_ortho = mo_ortho
            self.mo_e = mo_e
            t10 = time.time()
            print('time for xg, H, S2, Yg, Feff: %.3f' % (t10-t06))
        
            if old_suhf is not None:
                dE = E_suhf - old_suhf
                ddm = dm_ortho - old_dm
                conv = conv_check(E_suhf, dE, ddm, thresh, cyc)
            else:
                print(' E(SUHF) = %15.8f' % E_suhf)
            self.conv = conv
            util2.dump_chk(self.mol, self.chkfile, self.mo_e, self.mo_ortho, self.dm_ortho)
            cyc += 1
            if cyc >= max_cycle:
                print('SUHF not converged')
                break

        # extra cycle to remove level shift
        old_suhf = self.E_suhf
        old_dm = self.dm_ortho
        #old_Pgao = Pgao
        #old_Ggao = Ggao
        dm_reg = einsum('ij,tjk,lk->til', X, self.dm_ortho, X) # regular ao
        mo_ortho = np.array(self.mo_ortho)
        mo_reg = einsum('ij,tjk->tik', X, mo_ortho)
        self.dm_reg = dm_reg
        self.mo_reg = mo_reg
        if self.debug:
            print('dm_reg\n', dm_reg)
            print('mo_reg\n', mo_reg[0], '\n', mo_reg[1])
        if self.level_shift is not None:
            print('**** Extra Cycle %d ****' % cyc)
            veff = scf.uhf.get_veff(self.mol, dm_reg, vhfopt=self.vhfopt)
            veff_ortho = einsum('ji,tjk,kl->til', X, veff, X)
            if self.debug:
                print('dm (ortho)')
                print(self.dm_ortho)
            Fa_ortho, Fb_ortho = self.hcore_ortho + veff_ortho
            F_ortho = Fa_ortho, Fb_ortho
            if self.debug:
                print('Fock (ortho)\n', F_ortho)
                e_uhf, e_uhf_coul = scf.uhf.energy_elec(mf, self.dm_ortho, self.hcore_ortho, veff_ortho)
                print('E(UHF) = %12.6f' % e_uhf)
            dm_no, dm_expanded, no = find_NO(self, self.dm_ortho, na, nb)
            self.dm_no = dm_no
            self.no = no
            Dg, Ng, Pg = get_Ng(self.grids, self.no, self.dm_no, na+nb)
            self.Pg = Pg
            if self.debug:
                print('D(g) (NO)\n', Dg[0])
                print('N(g) (NO)\n', Ng[0])
                print('P(g) (NO)\n', Pg[0])
            Gg, Pg_ortho, _, _ = get_Gg(self.mol, Pg, self.no, X, opt=self.vhfopt)
            if self.debug:
                print('Pg_ortho\n', Pg_ortho[0])
                print('G(g) (NO)\n' , Gg[0])
            xg, yg, ciS, C_no = get_xg(self, self.no, na, nb, Ng)
            self.xg, self.ciS = xg, ciS
            #yg, ciS = util.get_yg(self, xg)
            trHg, ciH = get_H(self, self.hcore_ortho, self.no, Pg, Gg, xg)
            S2 = get_S2(self, Pg_ortho)
            Xg, Xg_int, Yg = get_Yg(self, Dg, Ng, self.dm_no, na+nb)
            Feff_ortho, H_suhf, F_mod_ortho = get_Feff(self, trHg, Gg, Ng, Pg, Dg, na+nb, Yg, Xg, F_ortho)
            E_suhf = self.energy_nuc + H_suhf
            self.E_suhf = E_suhf
            Faa = F_mod_ortho[:norb, :norb]
            Fbb = F_mod_ortho[norb:, norb:]
            F_mod_ortho = np.array([Faa,Fbb])
            mo_e, mo_ortho, dm_ortho = Diag_Feff(F_mod_ortho, na, nb)
            util2.dump_moe(mo_e, na, nb)
            if self.debug:
                print('e_a, e_b\n', mo_e[0], '\n', mo_e[1])
                print('v_a, v_b\n', mo_ortho[0], '\n', mo_ortho[1])
                print('P_a, P_b\n', dm_ortho[0],'\n', dm_ortho[1])
            self.dm_ortho = dm_ortho
            self.mo_ortho = mo_ortho
            self.mo_e = mo_e
            dm_reg = einsum('ij,tjk,lk->til', X, self.dm_ortho, X) # regular ao
            mo_ortho = np.array(self.mo_ortho)
            mo_reg = einsum('ij,tjk->tik', X, mo_ortho)
            self.dm_reg = dm_reg
            self.mo_reg = mo_reg
            #if old_suhf is not None:
            dE = E_suhf - old_suhf
            ddm = dm_ortho - old_dm
            max_ddm = abs(ddm).max()
            norm_ddm = np.linalg.norm(ddm) / (ddm.shape[-1]*np.sqrt(2))
            print('Extra E(SUHF) = %15.8f, delta E = %10.6g, MaxD(dm) = %10.6g, RMSD(dm) = %10.6g' % (E_suhf, dE, max_ddm, norm_ddm))
        t_aftercyc = time.time()
        print('time for cyc: %.3f' % (t_aftercyc-t_pre))

        if self.makedm:
            suhf_dm = sudm.make_1pdm(self, Dg, self.dm_no, C_no)
            t_dm = time.time()
            print('time for 1pdm: %.3f' % (t_dm-t_aftercyc))
            self.suhf_dm = suhf_dm
            self.natorb, self.natocc = sudm.natorb(self, suhf_dm)
            if self.tofch:
                S = self.mol.intor_symmetric('int1e_ovlp')
                util2.tofch(self.oldfch, self.natorb[2], self.natocc[2], S)
                #util2.tofch(self.oldfch, mo_reg, self.mo_e, S, 'SUHFMO')
            t_nat = time.time()
            print('time for natorb: %.3f' % (t_nat-t_dm))

        t_end = time.time()
        print('time tot: %.3f' % (t_end-t_start))
        return E_suhf, self.conv

    def get_JKg(self):
        return get_JKg(self.mol, self.Pg, self.no, self.X)[:2]

    def get_EX(self):
        Jg, Kg = self.get_JKg()
        return get_EX(self, self.no, self.Pg, Kg, self.xg)[1]
        

def conv_check(E_suhf, dE, ddm, thresh, cyc):
    max_ddm = abs(ddm).max()
    norm_ddm = np.linalg.norm(ddm) / (ddm.shape[-1]*np.sqrt(2))
    if abs(dE)< (thresh*1e-2) and max_ddm < (thresh*1e2) and norm_ddm < thresh:
        conv = True
        print('\n***************')
        print('SUHF converged at cycle %d' %cyc)
        print('Final E(SUHF) = %15.8f, delta E = %10.6g, MaxD(dm) = %10.6g, RMSD(dm) = %10.6g' % (E_suhf, dE, max_ddm, norm_ddm))
    else:
        conv = False
        print(' E = %15.8f, delta E = %10.6g, MaxD(dm) = %10.6g, RMSD(dm) = %10.6g' % (E_suhf, dE, max_ddm, norm_ddm))
    return conv

def lev_shift(s, dm, f, shift):
    new_f = (scf.hf.level_shift(s, dm[0], f[0], shift),
             scf.hf.level_shift(s, dm[1], f[1], shift)
             )
    return np.array(new_f)
