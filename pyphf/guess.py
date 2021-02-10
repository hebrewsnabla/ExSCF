from pyscf import gto, scf, dft
import numpy as np
from fch2py import fch2py

def from_fchk(xyz, bas, fch, cycle=1):
    mol = gto.Mole()
    mol.atom = xyz
    #with open(xyz, 'r') as f:
    #    mol.atom = f.read()
    #print(mol.atom)
    mol.basis = bas
    #mol.output = 'test.pylog'
    mol.verbose = 4
    mol.build()
    
    mf = scf.UHF(mol)
    #mf.init_guess = '1e'
    mf.init_guess_breaksym = True
    mf.max_cycle = 1
    mf.kernel()
    
    # read MOs from .fch(k) file
    nbf = mf.mo_coeff[0].shape[0]
    nif = mf.mo_coeff[0].shape[1]
    S = mol.intor_symmetric('int1e_ovlp')
    Sdiag = S.diagonal()
    alpha_coeff = fch2py(fch, nbf, nif, Sdiag, 'a')
    beta_coeff  = fch2py(fch, nbf, nif, Sdiag, 'b')
    mf.mo_coeff = (alpha_coeff, beta_coeff)
    # read done
    
    dm = mf.make_rdm1()
    mf.max_cycle = cycle
    mf.kernel(dm)
    return mf

def mix(xyz, bas):
    mol = gto.Mole()
    mol.atom = xyz
    #with open(xyz, 'r') as f:
    #    mol.atom = f.read()
    #print(mol.atom)
    mol.basis = bas
    #mol.output = 'test.pylog'
    mol.verbose = 4
    mol.build()
    
    mf = scf.UHF(mol)
    #mf.init_guess = '1e'
    mf.init_guess_breaksym = True
    #mf.max_cycle = 1
    mf.kernel()
    
    #dm = mf.make_rdm1()
    #mf.max_cycle = 0
    #mf.kernel(dm)
    return mf

def from_frag(xyz, bas, frags, chgs, spins, cycle=2, xc=None):
    mol = gto.Mole()
    mol.atom = xyz
    mol.basis = bas
    mol.build()
    
    dm, mo, occ = guess_frag(mol, frags, chgs, spins)
    if xc is None:
        mf = scf.UHF(mol)
    else:
        mf = dft.UKS(mol)
        mf.xc = xc
    mf.verbose = 6
    #mf.conv_tol = 1e-2
    mf.max_cycle = cycle
    mf.kernel(dm0 = dm)
    return mf


def guess_frag(mol, frags, chgs, spins):
    '''
    frags: e.g. [[0], [1]] for N2
    '''
    #mol.build()
    print('generating fragment guess')
    atom = mol.format_atom(mol.atom, unit=1)
    #print(atom)
    fraga, fragb = frags
    chga, chgb = chgs
    spina, spinb = spins
    atoma = [atom[i] for i in fraga]
    atomb = [atom[i] for i in fragb]
    print('fragments:', atoma, atomb)
    ca_a, cb_a, na_a, nb_a = do_uhf(atoma, mol.basis, chga, spina)
    ca_b, cb_b, na_b, nb_b = do_uhf(atomb, mol.basis, chgb, spinb)
    print('       na   nb')
    print('atom1  %2d   %2d' % (na_a, nb_a))
    print('atom2  %2d   %2d' % (na_b, nb_b))
    #print(mo_a)
    #print(mo_b)
    nbasa = ca_a.shape[0]
    nbasb = ca_b.shape[0]
    ca = np.vstack((
                    np.hstack((ca_a[:,:na_a], np.zeros((nbasa,na_b)), ca_a[:,na_a:], np.zeros((nbasa, ca_b.shape[1]-na_b)) )),
                    np.hstack((np.zeros((nbasb, na_a)), ca_b[:,:na_b], np.zeros((nbasb, ca_a.shape[1]-na_a)), ca_b[:,na_b:]))
                  ))
    cb = np.vstack((
                    np.hstack((cb_a[:,:nb_a], np.zeros((nbasa,nb_b)), cb_a[:,nb_a:], np.zeros((nbasa, cb_b.shape[1]-nb_b)) )),
                    np.hstack((np.zeros((nbasb, nb_a)), cb_b[:,:nb_b], np.zeros((nbasb, cb_a.shape[1]-nb_a)), cb_b[:,nb_b:]))
                  ))
    mo = np.array([ca, cb])
    na = na_a + na_b
    nb = nb_a + nb_b
    #print(ca.shape, cb.shape)
    occa = np.hstack((np.ones(na), np.zeros(ca.shape[1]-na))) 
    occb = np.hstack((np.ones(nb), np.zeros(cb.shape[1]-nb)))
    occ = np.array([occa, occb]) 
    #print(occ)
    dm = scf.uhf.make_rdm1(mo, occ)
    #print(dm.shape)
    return dm, mo, occ
    
def do_uhf(atoma, basisa, chga, spina):
    mola = gto.Mole()
    mola.atom = atoma
    mola.basis = basisa
    mola.charge = chga
    mola.spin = spina
    mola.build()
    mfa = scf.UHF(mola)
    mfa.kernel()
    #print(mfa.nelec)
    ca, cb = mfa.mo_coeff
    na, nb = mfa.nelec
    return ca, cb, na, nb
