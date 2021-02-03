import numpy as np
from py2fch import py2fch
import os
from functools import partial

print = partial(print, flush=True)
einsum = partial(np.einsum, optimize=True)

def stack22(aa, ab, ba, bb):
    tot = np.vstack((
        np.hstack((aa,ab)),
        np.hstack((ba,bb))
    ))
    return tot

def reg2ortho(dm, X, forward=True):
    if forward:
        return np.einsum('ji,jk,kl', X, dm, X)
    else:
        return np.einsum('ij,jk,lk', X, dm, X)

def tofch(oldfch, natorb, natocc, S, flag):
    fch = oldfch.split('.fch')[0] + '_' + flag + '.fch'
    os.system('cp %s %s' % (oldfch, fch))
    #S = suhf.mol.intor_symmetric('int1e_ovlp')
    Sdiag = S.diagonal()
    natorb_a, natorb_b = natorb
    natocc_a, natocc_b = natocc
    nbfa = natorb_a.shape[0]
    nifa = natorb_a.shape[1]
    py2fch(fch, nbfa, nifa, natorb_a, Sdiag, 'a', natocc_a)
    nbfb = natorb_b.shape[0]
    nifb = natorb_b.shape[1]
    py2fch(fch, nbfb, nifb, natorb_b, Sdiag, 'b', natocc_b)

def dump_moe(moe, na, nb):
    ea = moe[0]
    eb = moe[1]
    avir = len(ea) - na
    bvir = len(eb) - nb
    print('Alpha occ %d vir %d; Beta occ %d vir %d' % (na, avir, nb, bvir))
    amin = max(0, na-6)
    amax = min(na+6, len(ea))
    print('Alpha energies: ', ea[amin:na], '<- HOMO')
    print('         LUMO-> ', ea[na:amax])
    bmin = max(0, nb-6)
    bmax = min(nb+6, len(eb))
    print('Beta energies:  ', eb[bmin:nb], '<- HOMO')
    print('         LUMO-> ', eb[nb:bmax])
