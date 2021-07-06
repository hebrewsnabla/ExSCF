from pyphf import suscf, deltascf, jk, util2
import numpy as np
import scipy
from functools import partial

print = partial(print, flush=True)
einsum = partial(np.einsum, optimize=True)


def ci0(mf, excis):
    ndet = len(excis) 
    S = np.zeros((ndet+1, ndet+1))
    H = np.zeros((ndet+1, ndet+1))
    S[0,0] = mf.ciS
    H[0,0] = mf.ciH
    for r in range(ndet):
        exc = excis[r]
        s, h = ci_cross(mf, exc[0], exc[1]) 
        S[0,1+r] = s
        H[0,1+r] = h
    
    for r in range(ndet):
        excr = excis[r]
        for t in range(r,ndet):
            exct = excis[t]
            s, h = ci_2cross(mf, excr[0], excr[1], exct[0], exct[1])
            S[1+r, 1+t] = s
            H[1+r, 1+t] = h
    #print('S and H\n', S, '\n', H)
    S = S + np.tril(S.T,-1)
    H = H + np.tril(H.T,-1)
    print('S and H\n', S, '\n', H)
    e, c = scipy.linalg.eigh(H, S)
    print(e, '\n', c)
    return 0

def ci_2cross(mf, aexci, bexci, aexci2, bexci2, mo=None, doW=False):
    if mo is None:
        C1 = mf.mo_ortho
    else:
        C1 = mo
    occ = mf.mo_occ
    occ1 = deltascf.set_occ(occ, aexci, bexci)
    occ2 = deltascf.set_occ(occ, aexci2, bexci2)
    #C1 = Ca[:,occ[0]], Cb[:,occ[1]]
    #C2 = Ca[:,occ2[0]], Cb[:,occ2[1]]
    C1_expd = suscf.expd(C1, occ1)
    C2_expd = suscf.expd(C1, occ2)
#    print(C1_expd)
#    print(C2_expd)

    ciS, ciH = cross(mf, C1_expd, C2_expd, doW)
    return ciS, ciH

def all_cross(mf, na, nb, nvira, nvirb, mo=None):
    #na, nb = mf.nelec
    if mo is None:
        C1 = mf.mo_ortho
    else:
        C1 = mo
    C1_expd = suscf.expd(C1, mf.mo_occ)
    Dg, Mg, xg, Pg = get_DxP(mf, mf.norb, C1_expd, C1_expd, na+nb)
    

def ci_cross(mf, aexci, bexci, mo=None, doW=False):
    if mo is None:
        C1 = mf.mo_ortho
    else:
        C1 = mo
    occ = mf.mo_occ
    occ2 = deltascf.set_occ(occ, aexci, bexci)
    #C1 = Ca[:,occ[0]], Cb[:,occ[1]]
    #C2 = Ca[:,occ2[0]], Cb[:,occ2[1]]
    C1_expd = suscf.expd(C1, occ)
    C2_expd = suscf.expd(C1, occ2)
#    print(C1_expd)
#    print(C2_expd)

    #dm_no, _, no = suscf.find_NO(mf, mf.dm_ortho, mf.mo_occ)
    #Dg, Ng, Pg = get_Ng(mf.grids, no, dm_no)
    #C_no = suscf.get_xg(suhf, no, mf.mo_occ, Ng)[-1]
    #C2_no = suscf.get_xg(suhf, no, occ2, Ng)[-1]
    #cross(mf, C1_expd, C1_expd)
    ciS, ciH = cross(mf, C1_expd, C2_expd, doW)
#    Dg, Mg, xg, Pg = get_DxP(mf.grids, mf.norb, C1_expd, C1_expd, na+nb)
#    Gg, _, _ = jk.get_Gg_ortho(mol, Pg, X)
#    trHg, ciH = get_H(mf, mf.hcore_ortho,  Pg, Gg, xg)
    return ciS, ciH

def all_cross(mf, na, nb, nvira, nvirb, mo=None):
    #na, nb = mf.nelec
    if mo is None:
        C1 = mf.mo_ortho
    else:
        C1 = mo
    C1_expd = suscf.expd(C1, mf.mo_occ)
    Dg, Mg, xg, Pg = get_DxP(mf, mf.norb, C1_expd, C1_expd, na+nb)
    
    #Cg = Diag_Mg(Mg, C1_expd, na+nb)
    #check_ortho(Cg, Dg)
    Sia = get_Sia(mf, na, nb, nvira, nvirb)
    Siajb = get_Siajb(mf, na, nb, nvira, nvirb)

    return Sia, Siajb

def get_Sia(mf, na, nb, nvira, nvirb):
    #occ = 
    Sia = np.zeros((na, nvira))
    for i in range(na):
        for a in range(na, na+nvira):
            #occ2 = deltascf.set_occ(mo_occ, [[i],[a]], [[],[]])
            #o1 = expd_occ(mo_occ, na, nvira)
            #o2 = expd_occ(occ2m na, nvira)
            Sia[i,a-na] = ci_cross(mf, [[i],[a]], [[],[]])[0]
    print('Sia\n', Sia)
    return Sia

def get_Siajb(mf, na, nb, nvira, nvirb):
    #occ = 
    Siajb_ab = np.zeros((na, nvira, nb, nvirb))
    for i in range(na):
        for a in range(nvira):
            for j in range(nb):
                for b in range(nvirb):
                    Siajb_ab[i,a,j,b] = ci_cross(mf, [[i],[a+na]], [[j],[b+nb]])[0]
    if na > 1:
        Siajb_aa = np.zeros((na, nvira, na, nvira))
        for i in range(na):
            for a in range(na, na+nvira):
                Siajb_aa[i,a-na] = ci_cross(mf, [[i],[a]], [[],[]])[0]
    print('Siajb\n', Siajb_ab)
    return Siajb_ab

def expd_occ(occ, na, nvira):
    newocc = occ[0].copy()
     

def Diag_Mg(Mg, C, occ):
    Cg = []
    for mg in Mg:
        mg_h = mg + 2*np.triu(mg[:occ,:occ]).T
        u,s,v = scipy.linalg.svd(mg_h,lapack_driver='gesvd')
        #u,s,v = scipy.linalg.svd(mg[:occ,:occ],lapack_driver='gesvd')
        #s = np.flip(s)
        #u = np.flip(u, axis=1)
        print(s)
        print(u)
        print(v)
        #newC = C.copy()
        newC = np.dot(C[:,:occ], u)
        Cg.append(newC)
        #newmg = einsum('ji,jk,kl->il', u, mg_h, u)
        newmg = einsum('ji,jk,kl->il', u, mg[:occ,:occ], u)
        print(newC)
        print(newmg)
    return Cg

def check_ortho(Cg, Dg):
    for i in range(len(Cg)):
        s = einsum('ji,jk,kl->il', Cg[i], Dg[i], Cg[i])
        print(s)

def cross(mf, C1_expd, C2_expd, doW):
    na, nb = mf.nelec
    mol = mf.mol
    X = mf.X
    if doW:
        Dg, Mg, xg, Pg = get_DxP2(mf, mf.norb, C1_expd, C2_expd, na+nb)
    else:
        Dg, Mg, xg, Pg = get_DxP(mf, mf.norb, C1_expd, C2_expd, na+nb)
    ciS = get_ciS(mf, xg)
    Gg, _, _ = jk.get_Gg_ortho(mol, Pg, X)
    trHg, ciH = get_H(mf, mf.hcore_ortho,  Pg, Gg, xg)
    return ciS, ciH

def expd(mo):
    mo_expd = np.vstack((
        np.hstack((mo[0], np.zeros_like(mo[0]) )),
        np.hstack(( np.zeros_like(mo[1]), mo[1]))
    ))
    return mo_expd

def get_DxP(mf, norb, C1, C2, occ):
    grids = mf.grids
    C1o = C1[:,:occ]
    C2o = C2[:,:occ]
    C1v = C1[:,occ:]
    C2v = C2[:,occ:]
    Dg = []
    Mg = []
    Mgov = []
    Mgvo = []
    Mgvv = []
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
        mg = einsum('ji,jk,kl->il', C1o, dg, C2o)
#        mgov = einsum('ji,jk,kl->il', C1o, dg, C2v)
#        mgvo = einsum('ji,jk,kl->il', C1v, dg, C2o)
#        mgvv = einsum('ji,jk,kl->il', C1v, dg, C2v)
        
        x = np.linalg.det(mg)
        Mg.append(mg)
#        Mgov.append(mgov)
#        Mgvo.append(mgvo)
#        Mgvv.append(mgvv)
        xg.append(x)
        pg = einsum('ij,jk,kl,ml->im', dg, C2o, np.linalg.inv(mg), C1o)
        #print(pg)
        Pg.append(pg)
    print('xg', xg)
    print('mg\n', Mg[0], '\n', Mg[-1])
    print('pg\n', Pg[0])

#    Mg_avg, _, _, _ = get_M(mf, Mg, Mgov, Mgvo, Mgvv)
    #ciS = suhf.integr_beta(np.array(xg))
    #print('ciS', ciS)
    return Dg, Mg, xg, Pg

def get_DxP2(mf, norb, C1, C2, occ):
    grids = mf.grids
    C1o = C1[:,:occ]
    C2o = C2[:,:occ]
#    C1v = C1[:,occ:]
#    C2v = C2[:,occ:]
    Dg = []
    Mg = []
    Mgoo = []
#    Mgov = []
#    Mgvo = []
#    Mgvv = []
    xg = []
    Pg = []
    Wg = []
    #print('N(g)')
    for beta in grids:
        #norb = int(len(no)/2)
        dg1 = np.eye(norb) * np.cos(beta/2)
        dg2 = np.eye(norb) * (-np.sin(beta/2))
        dg = np.vstack((
            np.hstack((dg1, -dg2)),
            np.hstack((dg2, dg1))
        ))

#        dg_no = einsum('ji,jk,kl->il', no, dg, no)

#        Dg.append(dg_no)
        Dg.append(dg)
#        tg = einsum('ij,jk,kl->il', dm[:occ,:], dg_no, dm[:,:occ])
#        ng = np.linalg.inv(tg)
        #print(ng)
#        Ng.append(ng)
#        pg = einsum('ij,jk,kl,lm->im', dg_no, dm[:,:occ], ng, dm[:occ,:])
        mg = einsum('ji,jk,kl->il', C1, dg, C2)
#        mgov = einsum('ji,jk,kl->il', C1o, dg, C2v)
#        mgvo = einsum('ji,jk,kl->il', C1v, dg, C2o)
#        mgvv = einsum('ji,jk,kl->il', C1v, dg, C2v)
        mgoo = mg[:occ,:occ] 
        mgov = mg[:occ,occ:] 
        mgvo = mg[occ:,:occ] 
        mgvv = mg[occ:,occ:] 
        Woo = np.linalg.inv(mgoo)
        Wov = np.dot(Woo, mgov)
        Wvo = np.dot(mgvo, Woo)
        Wvv = np.linalg.inv(mgvv)
        wg = util2.stack22(Woo, Wov, Wvo, Wvv)
        Wg.append(wg)
        
        x = np.linalg.det(mgoo)
        Mg.append(mg)
        Mgoo.append(mgoo)
#        Mgvo.append(mgvo)
#        Mgvv.append(mgvv)
        xg.append(x)
        pg = einsum('ij,jk,kl,ml->im', dg, C2o, np.linalg.inv(mgoo), C1o)
        #print(pg)
        Pg.append(pg)
    print('xg', xg)
    print('mg\n', Mg[0], '\n', Mg[-1])
    print('pg\n', Pg[0])
    Wg_avg = mf.integr_beta(np.array(Wg))
    #print('Wg\n', Wg_avg)
    #print(Wg)

#    Mg_avg, _, _, _ = get_M(mf, Mg, Mgov, Mgvo, Mgvv)
#    Mg_avg = mf.integr_beta(np.array(Mg))
#    Mgoo_avg = Mg_avg[:occ,:occ]

#    print('Mg', Mg_avg)
    return Dg, Mgoo, xg, Pg

def get_ciS(suhf, xg):
    ciS = suhf.integr_beta(np.array(xg))
    print('ciS', ciS)
    return ciS

def get_M(suhf, Mg, Mgov, Mgvo, Mgvv):
    Mg_avg = suhf.integr_beta(np.array(Mg))
    Mgov_avg = suhf.integr_beta(np.array(Mgov))
    Mgvo_avg = suhf.integr_beta(np.array(Mgvo))
    Mgvv_avg = suhf.integr_beta(np.array(Mgvv))
    print('Mg', Mg_avg, '\n', Mgov_avg, '\n', Mgvo_avg, '\n', Mgvv_avg)
    Woo = np.linalg.inv(Mg_avg)
    Wov = np.dot(Woo, Mgov_avg)
    Wvo = np.dot(Mgvo_avg, Woo)
    Wvv = np.linalg.inv(Mgvv_avg)
    print('W', Woo, '\n', Wov, '\n', Wvo, '\n', Wvv)
    return Mg_avg, Mgov_avg,  Mgvo_avg, Mgvv_avg

def get_H(suhf, hcore_ortho, Pg, Gg, xg):
    #print(hcore_ortho)
    hcore_ortho = np.vstack((
        np.hstack((hcore_ortho, np.zeros(hcore_ortho.shape))),
        np.hstack((np.zeros(hcore_ortho.shape), hcore_ortho))
    ))
    #hcore_no = einsum('ji,jk,kl->il', no, hcore_ortho, no)
    #suhf.hcore_no = hcore_no
    #if suhf.debug2:
    #    print(hcore_no)
    trHg = np.zeros(len(Pg))
    for i, pg in enumerate(Pg):
        H = np.trace(np.dot(hcore_ortho, pg)) + 0.5 * np.trace(np.dot(Gg[i], pg))
        #print(H)
        #H = H * xg[i]
        trHg[i] = H
        #print(i, H*xg[i])
    ciH = suhf.integr_beta(trHg*xg)
    print('ciH', ciH)
    return trHg, ciH


