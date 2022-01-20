import numpy as np
from pyphf import suscf
import copy
from functools import reduce

def make_s02(dg, mo, nocc, nvir, mo_occ, ovlp):
    Siajb = np.zeros((nocc, nocc, nvir, nvir))
    for i in range(nocc):
        #for j in list(range(i)) + list(range(i+1,nocc)):
        for j in range(i):
            for a in range(nvir):
                #for b in list(range(a)) + list(range(a+1,nvir)):
                for b in range(a):
                    Siajb[i,j,a,b] = ci_cross(dg, mo, [i,a+nocc], [j,b+nocc], 
                                                 mo_occ, ovlp)
    #for i in range(nocc):
    #    for j in range(i):
    #        for a in range(nvir):
    #            for b in range(a+1,nvir):
    #                Siajb[i,j,a,b] = ci_cross(dg, mo, [[i,j],[a+nocc, b+nocc]], 
    #                                             mo_occ, ovlp)
    return Siajb

def make_s22(dg, mo, nocc, nvir, mo_occ, ovlp, leftexci):
    Siajb = np.zeros((nocc, nocc, nvir, nvir))
    #Siajb = np.zeros(( nvir, nvir))
    k,l,c,d = leftexci
    kc = [k,c+nocc]
    if l is None:
        ld = None
    else:
        ld = [l,d+nocc]
    for i in range(nocc):
        #for j in list(range(i)) + list(range(i+1,nocc)):
        for j in range(i):
    #for i in [k]:
    #    for j in [l]:
            for a in range(nvir):
                #for b in list(range(a)) + list(range(a+1,nvir)):
                for b in range(a):
                    Siajb[i,j,a,b] = ci_cross(dg, mo, [i,a+nocc], [j,b+nocc],
                                                 mo_occ, ovlp, [kc,ld])
                    #Siajb[a,b] = ci_cross(dg, mo, [i,a+nocc], [j,b+nocc],
                    #                             mo_occ, ovlp, [kc,ld])
    return Siajb

def ci_cross(dg, mo, exci1, exci2, mo_occ, S, leftexci=None):
    #mo_occ2 = set_occ(mo_occ, exci1)
    #mo_occ2 = set_occ(mo_occ2, exci1)
    #C1_expd = suscf.expd(mo, mo_occ)
    #C2_expd = suscf.expd(mo, mo_occ2)
    #nmo = len(mo_occ)
    #nocc = np.count_nonzero(mo_occ)
    #nvir = 
    if leftexci is None:
        C1_expd = mo[:,mo_occ==1]
    else:
        kc, ld = leftexci
        C1_expd = permute(mo, kc)
        if ld is not None:
            C1_expd = permute(C1_expd, ld)
        C1_expd = C1_expd[:,mo_occ==1]
    C2_expd = permute(permute(mo, exci1), exci2)[:,mo_occ==1]
    mg = reduce(np.dot, (C1_expd.T, S, dg, C2_expd))
    s = np.linalg.det(mg)
    if leftexci is None:
        if abs(s) > 1e-6:
            print(exci1, exci2, s)
    elif ld is None:
        if abs(s) > 1e-6:
            print(kc, exci1, exci2, s)
    return s

def permute(mo, exci):
    i,a = exci
    assert(i<a)
    newmo = np.hstack((mo[:,:i], mo[:,a:a+1], mo[:,i+1:a], mo[:,i:i+1], mo[:,a+1:]))
    return newmo

def set_occ(occ0, aexci=[[],[]], bexci=[[],[]]):
    ''' aexci: [[3,4],[5,6]]
    '''
#    if not mf.converged:
#        mf.kernel()
    #e0 = mf.e_tot
    #mo0 = mf.mo_coeff
    #occ0 = mf.mo_occ
    occ = copy.copy(occ0)
#    if type=='U':
    for i in aexci[0]:
        occ[i] -= 1
    for j in aexci[1]:
        occ[j] += 1
    #for i in bexci[0]:
    #    occ[1][i] -= 1
    #for j in bexci[1]:
    #    occ[1][j] += 1
    #print('Excitation: Alpha ',aexci[0], ' -> ', aexci[1])
    #print('             Beta ',bexci[0], ' -> ', bexci[1])
    length = max(aexci[1]) + 1
    #print('Former occ pattern: Alpha', occ0[0][:length])
    #print('                     Beta', occ0[1][:length])
    #print('New occ pattern:    Alpha', occ[:length])
    #print('                     Beta', occ[1][:length])
    return occ
