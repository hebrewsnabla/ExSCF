import numpy as np
from pyphf import util2
import sympy
from sympy.physics.quantum.cg import CG
import os


def make_1pdm(suhf, Dg, dm_no, C_no):
    print('******* suhf density *****')
    no = suhf.no
    na,nb = suhf.nelec
    occ = na+nb
    C_oo = C_no[:occ,:occ]
    Ngg = get_Ngg(Dg, dm_no, occ)
    Pgg, Pgg_ortho = get_Pgg(Dg, dm_no, Ngg, occ, no)
    xgg = get_xgg(Ngg, C_oo)
    #wgtf0 = suhf.d
    S, Sz = suhf.S, suhf.Sz
    cgcoeff0, cgfloat0 = get_CG(S, Sz, 0, 0, S, Sz)
    cgcoeff1, cgfloat1 = get_CG(S, Sz, 1, 0, S, Sz)
    #cgcoeff2, cgfloat2 = get_CG(S, Sz, 2, 0, S, Sz)
    #print(type(cgfloat1))

    xggint = suhf.integr_beta(np.array(xgg), fac='ci')
    print('xggint', xggint)

    X = suhf.X
    def ortho2reg(dm):
        return util2.reg2ortho(dm, X, False)

    norb = int(Pgg_ortho[0].shape[0]/2)
    Onepdm_a = []
    Onepdm_b = []
    for i,pgg_ortho in enumerate(Pgg_ortho):
        x = xgg[i]
        pggaa = pgg_ortho[:norb, :norb]
        pggab = pgg_ortho[:norb, norb:]
        pggba = pgg_ortho[norb:, :norb]
        pggbb = pgg_ortho[norb:, norb:]

        pggaa_ao, pggab_ao, pggba_ao, pggbb_ao = list(map(ortho2reg, [pggaa, pggab, pggba, pggbb]))
        #print(pggbb_ao.dtype)
        #pgg_ao = util2.stack22(pggaa_ao, pggab_ao, pggba_ao, pggbb_ao)
        pgg1_ao = (pggaa_ao + pggbb_ao)  * x / np.sqrt(2.0)
        pgg2_ao = (pggaa_ao - pggbb_ao)  * x / np.sqrt(2.0)
        pggba_ao = pggba_ao  * (-x)
        pggab_ao = pggab_ao  * x
        #print(pgg1_ao)

        # Here we assume S = Sz
        onepdm_a = (pgg1_ao * cgfloat0 * cgfloat0 + pgg2_ao * cgfloat1 * cgfloat1) / np.sqrt(2.0)
        onepdm_b = (pgg1_ao * cgfloat0 * cgfloat0 - pgg2_ao * cgfloat1 * cgfloat1) / np.sqrt(2.0)
        Onepdm_a.append(onepdm_a)
        Onepdm_b.append(onepdm_b)
        #print(onepdm_a.dtype)
    int1pdm_a = suhf.integr_beta(np.array(Onepdm_a), fac='ci') / xggint
    #print(type(int1pdm_a), int1pdm_a.dtype)
    int1pdm_b = suhf.integr_beta(np.array(Onepdm_b), fac='ci') / xggint
    print('SUHF DM alpha\n', int1pdm_a, '\nSUHF DM beta\n', int1pdm_b)
    return [int1pdm_a, int1pdm_b]

def natorb(suhf, dm):
    #print(type(dm[0]), dm[0].dtype)
    XS = suhf.XS
    X = suhf.X
    # P = S^0.5 * dm * S^0.5
    P0 = np.einsum('ij,jk,lk->il', XS, dm[0], XS)
    P1 = np.einsum('ij,jk,lk->il', XS, dm[1], XS)
    natocc_a, natorb_a = np.linalg.eigh(-P0)
    natocc_b, natorb_b = np.linalg.eigh(-P1)
    natorb_a = np.dot(X, natorb_a)
    natorb_b = np.dot(X, natorb_b)
    natocc_a = -1 * natocc_a
    natocc_b = -1 * natocc_b
    print('SUHF natural orbitals alpha\n', natorb_a, '\nSUHF NO occ alpha\n', natocc_a)
    print('SUHF natural orbitals beta\n', natorb_b, '\nSUHF NO occ beta\n', natocc_b)
    return [natorb_a, natorb_b], [natocc_a, natocc_b]


def get_CG(j1, m1, j2, m2, j, m):
    j1, m1, j2, m2, j, m = list(map(sympy.Rational, [j1, m1, j2, m2, j, m]))
    CGfunc = CG(j1, m1, j2, m2, j, m)
    CGcoeff =  CGfunc.doit()
    print('Clebsch-Gordan coeff: < {} {} {} {} | {} {} > = {}'.format(j1,m1,j2,m2,j,m,CGcoeff))
    return CGcoeff, float(CGcoeff.evalf())


def get_Ngg(Dg, dm_no, occ):
    Ngg = []
    for dg in Dg:
        ngg_inv = np.einsum('ij,jl,lm->im', dm_no[:occ,:], dg, dm_no[:,:occ])
        ngg = np.linalg.inv(ngg_inv)
        Ngg.append(ngg)
        #print('Ngg')
        #print(ngg)
    return Ngg

def get_Pgg(Dg, dm_no, Ngg, occ, no):
    Pgg = []
    Pgg_ortho = []
    for i, dg in enumerate(Dg):
        pgg = np.einsum('ij,jk,kl,ln->in', dg, dm_no[:,:occ], Ngg[i], dm_no[:occ,:])
        Pgg.append(pgg)
        #print('pgg')
        #print(pgg)
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
    #print('xgg', xgg)
    return xgg