from pyscf import scf
from pyscf.lib import temporary_env
from pyphf import util2
from functools import partial
import numpy as np

print = partial(print, flush=True)
einsum = partial(np.einsum, optimize=True)

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
    if opt is not None:
        with temporary_env(opt, prescreen='CVHFnrs8_vk_prescreen'):
            vk = scf.hf.get_jk(mol, dm, hermi, opt, with_j=False)[1]
    else:
        vk = scf.hf.get_jk(mol, dm, hermi, opt, with_j=False)[1]

    return vk

def no2ortho(Pg, no):
    Pg_ortho = []
    for pg in Pg:
        pg_ortho = einsum('ij,jk,lk->il', no, pg, no)
        #print(pg_ortho)
        Pg_ortho.append(pg_ortho)
    return Pg_ortho

def get_Gg(mol, Pg, no, X, dm_last=None, Ggao_last=None, opt=None):
    Pg_ortho = no2ortho(Pg, no)
    Gg_ortho, Pgao, Ggao = get_Gg_ortho(mol, Pg_ortho, X, dm_last, Ggao_last, opt)
    Gg = ortho2no(Gg_ortho, no)
    return Gg, Gg_ortho, Pg_ortho, Pgao, Ggao


def get_Gg_ortho(mol, Pg_ortho, X, dm_last=None, Ggao_last=None, opt=None):
    #Pg_ortho = no2ortho(Pg)
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
    Gg_ortho = []
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
        Gg_ortho.append(gg)
    return Gg_ortho, Pgao, Ggao


def ortho2no(Gg_ortho, no):    
    Gg = []
    for gg in Gg_ortho:
        gg_no = einsum('ji,jk,kl->il', no, gg, no)
        Gg.append(gg_no)
    return Gg