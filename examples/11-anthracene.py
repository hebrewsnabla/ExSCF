from pyphf import suscf, guess
from pyscf import lib

lib.num_threads(24)
#from gaussian import *
#mol = load_mol_from_fch('anthracene_uhf.fch')
#mf = guess.mix(mol.atom, mol.basis, conv='tight')
mf = guess.from_fch_simp('anthracene_uhf.fch')

mf2 = suscf.SUHF(mf)
mf2.kernel()

