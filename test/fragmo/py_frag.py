from pyscf import gto, scf
import frag

mol = gto.Mole()
mol.atom = 'N 0.0 0.0 0.0; N 0.0 0.0 2.0'
mol.basis = 'cc-pvdz'
#mol.verbose = 6 
mol.build()

#mf = scf.UHF(mol)
#mf.init_guess = 'atom'
#mf.kernel()

dm, mo, occ = frag.guess_frag(mol, [[0],[1]], [0,0], [3,-3])
mf = scf.UHF(mol)
mf.verbose = 6
#mf.max_cycle = 1
#mf.kernel()

#mf.mo_coeff = mo
#mf.mo_occ = occ
#mf.max_cycle = 5
mf.kernel(dm0 = dm)
