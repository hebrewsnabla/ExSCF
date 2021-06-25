import numpy as np
try:
    from py2fch import py2fch
except:
    print('py2fch not found. Interface with fch is disabled. Install MOKIT if you need that.')
import os
from functools import partial
from pyscf.lib.chkfile import dump, save
from pyscf.lib.chkfile import load_mol, save_mol

print = partial(print, flush=True)
einsum = partial(np.einsum, optimize=True)

def dmlist(dmas, dmbs, fac=-1):
    ddms = []
    for i in range(len(dmas)):
        ddm = dmas[i] + fac*dmbs[i]
        ddms.append(ddm)
    return ddms

def stack22(aa, ab, ba, bb):
    tot = np.vstack((
        np.hstack((aa,ab)),
        np.hstack((ba,bb))
    ))
    return tot

def reg2ortho(dm, X, forward=True):
    if forward:
        return einsum('ji,jk,kl', X, dm, X)
    else:
        return einsum('ij,jk,lk', X, dm, X)

def tofch(oldfch, natorb, natocc, S, flag='SUHFNO'):
    fch = oldfch.split('.fch')[0] + '_' + flag + '.fch'
    os.system('fch_u2r %s' % oldfch)
    oldfch_r = oldfch.split('.fch')[0] + '_r.fch'
    os.system('mv %s %s' % (oldfch_r, fch))
    #S = suhf.mol.intor_symmetric('int1e_ovlp')
    Sdiag = S.diagonal()
    #natorb_a, natorb_b = natorb
    #natocc_a, natocc_b = natocc
    #nbfa = natorb_a.shape[0]
    #nifa = natorb_a.shape[1]
    nbf = natorb.shape[0]
    nif = natorb.shape[1]
    py2fch(fch, nbf, nif, natorb, Sdiag, 'a', natocc, True)
    #nbfb = natorb_b.shape[0]
    #nifb = natorb_b.shape[1]
    #py2fch(fch, nbfb, nifb, natorb_b, Sdiag, 'b', natocc_b)

def tofchmo(oldfch, orb, occ, S, flag='SUHFMO'):
    fch = oldfch.split('.fch')[0] + '_' + flag + '.fch'
    os.system('cp %s %s' % (oldfch, fch))
    Sdiag = S.diagonal()
    orb_a, orb_b = orb
    occ_a, occ_b = occ
    nbfa = orb_a.shape[0]
    nifa = orb_a.shape[1]
    py2fch(fch, nbfa, nifa, orb_a, Sdiag, 'a', occ_a, True)
    nbfb = orb_b.shape[0]
    nifb = orb_b.shape[1]
    py2fch(fch, nbfb, nifb, orb_b, Sdiag, 'b', occ_b, True)
    

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

def load(pchk):
    return 0

def dump_chk(mol, chkfile, e_tot, mo_e, mo, mo_occ, dm):
    save_mol(mol, chkfile)

    scf_dic = {'e_tot'    : e_tot,
               'mo_energy': mo_e,
               'mo_coeff' : mo,
               'mo_occ'   : mo_occ,
               'dm'       : dm}
    save(chkfile, 'scf', scf_dic)

def repo_info(repo_path):
    '''
    pyscf-2.0 lib/misc.py
    Repo location, version, git branch and commit ID
    '''

    def git_version(orig_head, head, branch):
        git_version = []
        if orig_head:
            git_version.append('GIT ORIG_HEAD %s' % orig_head)
        if branch:
            git_version.append('GIT HEAD (branch %s) %s' % (branch, head))
        elif head:
            git_version.append('GIT HEAD      %s' % head)
        return '\n'.join(git_version)

    repo_path = os.path.abspath(repo_path)

    if os.path.isdir(os.path.join(repo_path, '.git')):
        git_str = git_version(*git_info(repo_path))

    elif os.path.isdir(os.path.abspath(os.path.join(repo_path, '..', '.git'))):
        repo_path = os.path.abspath(os.path.join(repo_path, '..'))
        git_str = git_version(*git_info(repo_path))

    else:
        git_str = None

    info = {'path': repo_path}
    if git_str:
        info['git'] = git_str
    return info

def git_info(repo_path):
    orig_head = None
    head = None
    branch = None
    try:
        with open(os.path.join(repo_path, '.git', 'ORIG_HEAD'), 'r') as f:
            orig_head = f.read().strip()
    except IOError:
        pass

    try:
        head = os.path.join(repo_path, '.git', 'HEAD')
        with open(head, 'r') as f:
            head = f.read().splitlines()[0].strip()

        if head.startswith('ref:'):
            branch = os.path.basename(head)
            with open(os.path.join(repo_path, '.git', head.split(' ')[1]), 'r') as f:
                head = f.read().strip()
    except IOError:
        pass
    return orig_head, head, branch
