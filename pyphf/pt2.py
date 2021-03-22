'''
EMP2 and SUPT2
Tsuchimochi, T.; Van Voorhis, T. J Chem Phys 2014, 141, 164117
Tsuchimochi, T.; Ten-no, S. L. J Chem Theory Comput 2019, 15, 6688–6702
'''
import numpy as np
#import sympy as sym
#import scipy
from pyscf import gto, scf
from pyphf import  util2
#import os, sys
from functools import partial
import time

print = partial(print, flush=True)
einsum = partial(np.einsum, optimize=True)

def defm_Fock(mol, hcore_no, dm_no, no, X):
    g_no = defm_G(mol, dm_no, no, X)
    F0 = hcore_no + g_no
    return F0

def defm_G(mol, dm_no, no, X):
    #G = 
    p = dm_no

    p_ortho = einsum('ij,jk,lk->il', no, p, no)

    norb = int(p.shape[0]/2)

    paa = p[:norb, :norb] # ortho ao
    #print(pgaa)
    pab = p[:norb, norb:]
    pba = p[norb:, :norb]
    pbb = p[norb:, norb:]
    # X . P(g) . X^H
    paa_ao = einsum('ij,jk,lk->il', X, paa, X) # regular ao
    #print(pgaa_ao)
    pab_ao = einsum('ij,jk,lk->il', X, pab, X)
    pba_ao = einsum('ij,jk,lk->il', X, pba, X)
    pbb_ao = einsum('ij,jk,lk->il', X, pbb, X)
    #Pgaabb_ao = Pgaa_ao + Pgbb_ao
    #print(Pgaabb_ao.shape)
    #nao = Pgaabb_ao.shape[-1]
    #ndm = len(Pgab_ao)
    vj,vk = scf.hf.get_jk(mol, [paa_ao, pbb_ao], hermi=0)
    #print(vj.shape)
    Gaa_ao = vj[0] + vj[1] - vk[0]
    Gbb_ao = vj[0] + vj[1] - vk[1]
    #Ggbb_ao = scf.hf.get_jk(mol, Pgbb_ao, hermi=0)
    #print(ggaa_ao)
    Gab_ao = scf.hf.get_jk(mol, pab_ao, hermi=0)[1] *(-1)
    Gba_ao = scf.hf.get_jk(mol, pba_ao, hermi=0)[1] *(-1)
    #ggbb_ao = scf.uhf.get_veff(mol, [pgaa_ao, pgbb_ao], hermi=0)[1]
    # X^H . G(g) . X

    gaa = einsum('ji,jk,kl->il', X, Gaa_ao, X)  # ortho ao
    #print(ggaa)
    gab = einsum('ji,jk,kl->il', X, Gab_ao, X) 
    gba = einsum('ji,jk,kl->il', X, Gba_ao, X) 
    gbb = einsum('ji,jk,kl->il', X, Gbb_ao, X) 
    g = util2.stack22(gaa, gab, gba, gbb)
    g_no = einsum('ji,jk,kl->il', no, g, no)

    return g_no

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
    
    def kernel(self):
        np.set_printoptions(precision=6, linewidth=160, suppress=True)
        suhf = self.suhf
        F0 = defm_Fock(self.mol, suhf.hcore_no, suhf.dm_no, suhf.no, suhf.X )
        epsl = F0.diagonal()
        print(suhf.hcore_no)
        print(F0, '\n', epsl)
