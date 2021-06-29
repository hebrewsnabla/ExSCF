from pyphf import suscf, deltascf
import numpy as np
from functools import partial

print = partial(print, flush=True)
einsum = partial(np.einsum, optimize=True)

def test(mf):
    na, nb = mf.nelec
    C1 = mf.mo_ortho
    occ = mf.mo_occ
    occ2 = deltascf.set_occ(occ, [[0],[1]], [[],[]])
    #C1 = Ca[:,occ[0]], Cb[:,occ[1]]
    #C2 = Ca[:,occ2[0]], Cb[:,occ2[1]]
    C1_expd = suscf.expd(C1, occ)
    C2_expd = suscf.expd(C1, occ2)
    print(C1_expd)
    print(C2_expd)

    #dm_no, _, no = suscf.find_NO(mf, mf.dm_ortho, mf.mo_occ)
    #Dg, Ng, Pg = get_Ng(mf.grids, no, dm_no)
    #C_no = suscf.get_xg(suhf, no, mf.mo_occ, Ng)[-1]
    #C2_no = suscf.get_xg(suhf, no, occ2, Ng)[-1]
    Dg, Mg, xg, Pg = get_DxP(mf.grids, mf.norb, C1_expd, C1_expd, na+nb)
    Dg, Mg, xg, Pg = get_DxP(mf.grids, mf.norb, C1_expd, C2_expd, na+nb)

def expd(mo):
    mo_expd = np.vstack((
        np.hstack((mo[0], np.zeros_like(mo[0]) )),
        np.hstack(( np.zeros_like(mo[1]), mo[1]))
    ))
    return mo_expd

def get_DxP(grids, norb, C1, C2, occ):
    C1 = C1[:,:occ]
    C2 = C2[:,:occ]
    Dg = []
    Mg = []
    xg = []
    Pg = []
    #print('N(g)')
    for beta in grids:
        #norb = int(len(no)/2)
        dg1 = np.eye(norb) * np.cos(beta/2)
        dg2 = np.eye(norb) * (-np.sin(beta/2))
        dg = np.vstack((
            np.hstack((dg1, -dg2)),
            np.hstack((dg2, dg1))
        ))
        #print(dg.shape)
        #print(dg)
#        dg_no = einsum('ji,jk,kl->il', no, dg, no)
        #det_dg = np.linalg.det(dg_no)
        #print(det_dg)
        #print(dg_no)
#        Dg.append(dg_no)
        Dg.append(dg)
#        tg = einsum('ij,jk,kl->il', dm[:occ,:], dg_no, dm[:,:occ])
#        ng = np.linalg.inv(tg)
        #print(ng)
#        Ng.append(ng)
#        pg = einsum('ij,jk,kl,lm->im', dg_no, dm[:,:occ], ng, dm[:occ,:])
        mg = einsum('ji,jk,kl->il', C1, dg, C2)
        x = np.linalg.det(mg)
        Mg.append(mg)
        xg.append(x)
        pg = einsum('ij,jk,kl,ml->im', dg, C2, np.linalg.inv(mg), C1)
        #print(pg)
        Pg.append(pg)
    print('xg', xg)

    #ciS = suhf.integr_beta(np.array(xg))
    #print('ciS', ciS)
    return Dg, Mg, xg, Pg


