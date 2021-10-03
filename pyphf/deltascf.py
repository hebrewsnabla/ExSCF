import copy
#from pyscf import scf
from pyscf.scf import uhf, rohf
import numpy as np
from functools import reduce


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
        occ[0][i] -= 1
    for j in aexci[1]:
        occ[0][j] += 1
    for i in bexci[0]:
        occ[1][i] -= 1
    for j in bexci[1]:
        occ[1][j] += 1
    print('Excitation: Alpha ',aexci[0], ' -> ', aexci[1])
    print('             Beta ',bexci[0], ' -> ', bexci[1])
    length = max(aexci[1]+ bexci[1]) + 1
    print('Former occ pattern: Alpha', occ0[0][:length])
    print('                     Beta', occ0[1][:length])
    print('New occ pattern:    Alpha', occ[0][:length])
    print('                     Beta', occ[1][:length])
    return occ


def mom_occ(mf, occorb, setocc):
    '''Use maximum overlap method to determine occupation number for each orbital in every
    iteration. '''
    coef_occ_a = occorb[0][:, setocc[0]>0]
    coef_occ_b = occorb[1][:, setocc[1]>0]
    mo_energy = mf.mo_e
    mo_coeff = mf.mo_reg
    mo_occ = np.zeros_like(setocc)
    nocc_a = int(np.sum(setocc[0]))
    nocc_b = int(np.sum(setocc[1]))
    s_a = reduce(np.dot, (coef_occ_a.T, mf.ovlp, mo_coeff[0]))
    s_b = reduce(np.dot, (coef_occ_b.T, mf.ovlp, mo_coeff[1]))
    if mf.debug:
        print(s_a,'\n',s_b)
    #choose a subset of mo_coeff, which maximizes <old|now>
    ss_a = np.einsum('ij,ij->j', s_a, s_a)
    ss_b = np.einsum('ij,ij->j', s_b, s_b)
    idx_a = np.argsort(ss_a)[::-1]
    idx_b = np.argsort(ss_b)[::-1]
    mo_occ[0][idx_a[:nocc_a]] = 1.
    mo_occ[1][idx_b[:nocc_b]] = 1.
    length = max(idx_a[:nocc_a]+ idx_b[:nocc_b])
    print('Overlap: ', ss_a[:length], '\n', ss_b[:length])
    print(' New alpha occ pattern: ', mo_occ[0][:length])
    print(' New beta occ pattern:  ', mo_occ[1][:length])
    print(' Current alpha mo_energy(sorted): ', mo_energy[0][:length])
    print(' Current beta mo_energy(sorted):  ', mo_energy[1][:length])
    if (int(np.sum(mo_occ[0])) != nocc_a):
        raise ValueError('mom alpha electron occupation numbers do not match: %d, %d'%
                  (nocc_a, int(np.sum(mo_occ[0])))
                  )
    if (int(np.sum(mo_occ[1])) != nocc_b):
        raise ValueError('mom beta electron occupation numbers do not match: %d, %d'%
                  (nocc_b, int(np.sum(mo_occ[1])))
                  )
    return mo_occ
