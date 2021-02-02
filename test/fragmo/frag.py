from pyscf import gto, scf
import numpy as np

def guess_frag(mol, frags, chgs, spins):
    '''
    frags: e.g. [[0], [1]] for N2
    '''
    #mol.build()
    atom = mol.format_atom(mol.atom, unit=1)
    print(atom)
    fraga, fragb = frags
    chga, chgb = chgs
    spina, spinb = spins
    atoma = [atom[i] for i in fraga]
    atomb = [atom[i] for i in fragb]
    print(atoma)
    ca_a, cb_a, na_a, nb_a = do_uhf(atoma, mol.basis, chga, spina)
    ca_b, cb_b, na_b, nb_b = do_uhf(atomb, mol.basis, chgb, spinb)
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
    print(ca.shape, cb.shape)
    occa = np.hstack((np.ones(na), np.zeros(ca.shape[1]-na))) 
    occb = np.hstack((np.ones(nb), np.zeros(cb.shape[1]-nb)))
    occ = np.array([occa, occb]) 
    print(occ)
    dm = scf.uhf.make_rdm1(mo, occ)
    print(dm.shape)
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
    print(mfa.nelec)
    ca, cb = mfa.mo_coeff
    na, nb = mfa.nelec
    return ca, cb, na, nb
