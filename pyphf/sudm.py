import numpy as np
from pyphf import util2
import sympy
from sympy.physics.quantum.cg import CG
import os
from functools import partial

print = partial(print, flush=True)
einsum = partial(np.einsum, optimize=True)

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
    if suhf.debug:
        print('SUHF DM alpha\n', int1pdm_a, '\nSUHF DM beta\n', int1pdm_b)
    return [int1pdm_a, int1pdm_b]

def natorb(suhf, dm):
    #print(type(dm[0]), dm[0].dtype)
    XS = suhf.XS
    X = suhf.X
    dmab = dm[0] + dm[1]
    # P = S^0.5 * dm * S^0.5
    P0 = einsum('ij,jk,lk->il', XS, dm[0], XS)
    P1 = einsum('ij,jk,lk->il', XS, dm[1], XS)
    Pab = P0 + P1
    natocc_a, natorb_a = np.linalg.eigh(-P0)
    natocc_b, natorb_b = np.linalg.eigh(-P1)
    natocc_ab, natorb_ab = np.linalg.eigh(-Pab)
    natorb_a = np.dot(X, natorb_a)
    natorb_b = np.dot(X, natorb_b)
    natorb_ab = np.dot(X, natorb_ab)
    natocc_a = -1 * natocc_a
    natocc_b = -1 * natocc_b
    natocc_ab = -1 * natocc_ab
    if suhf.debug:
        print('SUHF natural orbitals alpha\n', natorb_a)
        print('SUHF natural orbitals beta\n', natorb_b)
        print('SUHF natural orbitals, total\n', natorb_ab)
    print('SUHF NO occ alpha: ', dump_occ(natocc_a, 1.0)[0])
    print('SUHF NO occ beta:  ', dump_occ(natocc_b, 1.0)[0])
    occ_ab, [core, act, ext] = dump_occ(natocc_ab, 2.0)
    print('SUHF NO occ total: ', occ_ab)
    print('core %d, active %d, external %d' % (core, act, ext))
    return [natorb_a, natorb_b, natorb_ab], [natocc_a, natocc_b, natocc_ab]

def dump_occ(occ, full=1.0):
    s = ''
    core = 0
    act = 0
    ext = 0
    for i in occ:
        if i>1e-6:
            s += '%.6f  '%i
            if i>(full*0.99):
                core += 1
            elif i>(full*0.01):
                act +=1
            else:
                ext +=1
        else:
            ext += 1
    return s + '...', [core, act, ext]

def get_CG(j1, m1, j2, m2, j, m):
    j1, m1, j2, m2, j, m = list(map(sympy.Rational, [j1, m1, j2, m2, j, m]))
    CGfunc = CG(j1, m1, j2, m2, j, m)
    CGcoeff =  CGfunc.doit()
    print('Clebsch-Gordan coeff: < {} {} {} {} | {} {} > = {}'.format(j1,m1,j2,m2,j,m,CGcoeff))
    return CGcoeff, float(CGcoeff.evalf())


def get_Ngg(Dg, dm_no, occ):
    Ngg = []
    for dg in Dg:
        ngg_inv = einsum('ij,jl,lm->im', dm_no[:occ,:], dg, dm_no[:,:occ])
        ngg = np.linalg.inv(ngg_inv)
        Ngg.append(ngg)
        #print('Ngg')
        #print(ngg)
    return Ngg

def get_Pgg(Dg, dm_no, Ngg, occ, no):
    Pgg = []
    Pgg_ortho = []
    for i, dg in enumerate(Dg):
        pgg = einsum('ij,jk,kl,ln->in', dg, dm_no[:,:occ], Ngg[i], dm_no[:occ,:])
        Pgg.append(pgg)
        #print('pgg')
        #print(pgg)
        pgg_ortho = einsum('ij,jk,lk->il', no, pgg, no)
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