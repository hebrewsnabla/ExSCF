import numpy as np
from pyphf import util2
import sympy
from sympy.physics.quantum.cg import CG


def make_1pdm(suhf, Dg, dm_no, na, nb, C_no, no):
    occ = na+nb
    C_oo = C_no[:occ,:occ]
    Ngg = get_Ngg(Dg, dm_no, occ)
    Pgg, Pgg_ortho = get_Pgg(Dg, dm_no, Ngg, occ, no)
    xgg = get_xgg(Ngg, C_oo)

    X = suhf.X
    def ortho2reg(dm):
        return reg2ortho(dm, X, False)

    norb = int(Pgg_ortho[0].shape[0]/2)
    for i,pgg_ortho in enumerate(Pgg_ortho):
        x = xgg[i]
        pggaa = pgg_ortho[:norb, :norb]
        pggab = pgg_ortho[:norb, norb:]
        pggba = pgg_ortho[norb:, :norb]
        pggbb = pgg_ortho[norb:, norb:]

        pggaa_ao, pggab_ao, pggba_ao, pggbb_ao = map(ortho2reg, [pggaa, pggab, pggba, pggbb])
        #pgg_ao = util2.stack22(pggaa_ao, pggab_ao, pggba_ao, pggbb_ao)
        pgg1_ao = (pggaa_ao + pggbb_ao)  * x / np.sqrt(2.0)
        pgg2_ao = (pggaa_ao - pggbb_ao)  * x / np.sqrt(2.0)
        pggba_ao = pggba_ao  * (-x)
        pggab_ao = pggab_ao  * x
        print(pgg1_ao)

    #wgtf0 = suhf.d
    S, Sz = suhf.S, suhf.Sz
    cgcoeff0, cgfloat0 = get_CG(S, Sz, 0, 0, S, Sz)
    cgcoeff1, cgfloat1 = get_CG(S, Sz, 1, 0, S, Sz)
    cgcoeff2, cgfloat2 = get_CG(S, Sz, 2, 0, S, Sz)

    xggint = suhf.integr_beta(np.array(xgg), fac='ci')
    print('xggint', xggint)

    #onepdm = Pgg_ao * xgg * xgg

def get_CG(j1, m1, j2, m2, j, m):
    j1, m1, j2, m2, j, m = map(sympy.Rational, [j1, m1, j2, m2, j, m])
    CGfunc = CG(j1, m1, j2, m2, j, m)
    CGcoeff =  CGfunc.doit()
    return CGcoeff, CGcoeff.evalf()

def reg2ortho(dm, X, forward=True):
    if forward:
        return np.einsum('ji,jk,kl', X, dm, X)
    else:
        return np.einsum('ij,jk,lk', X, dm, X)


def get_Ngg(Dg, dm_no, occ):
    Ngg = []
    for dg in Dg:
        ngg_inv = np.einsum('ij,jl,lm->im', dm_no[:occ,:], dg, dm_no[:,:occ])
        ngg = np.linalg.inv(ngg_inv)
        Ngg.append(ngg)
        print('Ngg')
        print(ngg)
    return Ngg

def get_Pgg(Dg, dm_no, Ngg, occ, no):
    Pgg = []
    Pgg_ortho = []
    for i, dg in enumerate(Dg):
        pgg = np.einsum('ij,jk,kl,ln->in', dg, dm_no[:,:occ], Ngg[i], dm_no[:occ,:])
        Pgg.append(pgg)
        print('pgg')
        print(pgg)
        pgg_ortho = np.einsum('ij,jk,lk->il', no, pgg, no)
        Pgg_ortho.append(pgg_ortho)
    return Pgg, Pgg_ortho

def get_xgg(Ngg, C_oo):
    detC = np.linalg.det(C_oo)
    xgg = []
    for ngg in Ngg:
        detn = np.linalg.det(ngg)
        x = 1.0 / (detC * detn * detC)
        xgg.append(x)
    print('xgg', xgg)
    return xgg