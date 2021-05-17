from pyscf import dft
import pyscf.dft.numint as numint
#from pyphf import numint
from pyphf import sudft

import numpy as np
from functools import partial
import time

print = partial(print, flush=True)
einsum = partial(np.einsum, optimize=True)

class CASDFT():
    def __init__(self, mc):
        print('\n******** %s ********' % self.__class__)
        #suhf = util.SUHF(guesshf)
        self.mc = mc
        self.mcxc = 'tpss'
        self.grids = 'default'
        self.output = None
        #self.dens = 'deformed' # or relaxed
        self.trunc = None
        self.nref = None
        self.ncore = None

    def kernel(self):
        #if self.suhf.E_suhf is None:
        #    self.suhf.kernel()
        #if self.output is not None:
        #    sys.stdout = open(self.output, 'a')
        print('truncation: %s' % self.trunc)
        print('***** Start DFT Correlation for CAS+DFT **********')
        #print('density: %s' % self.dens)
        t1 = time.time()
        mol = self.mc.mol
        E_mc = self.mc.e_tot
        dm = self.mc.make_rdm1()

        grids = sudft.set_grids(mol, self.grids)
        ni = ks._numint
        if self.trunc is None:
            if self.mcxc == 'CS':
                mcxc = 'MGGA_C_CS'
                n, exc = sudft.get_exc(ni, mol, grids, 'HF,%s'%mcxc, dm)
            else:
                n, exc, vxc = ni.nr_uks(mol, grids, 'HF,%s'%self.mcxc, dm)
        elif self.trunc == 'f' or self.trunc == 'fc':
            natorb = self.mc.mo_coeff
            natocc = self.mc.mo_occ
            #natocc = natocc[0] + natocc[1]
            print('natocc', natocc)
            if self.nref is None:
                ref = [2.0 if occ > 1e-2 else 0.0 for occ in natocc]
            else:
                nref = self.nref
                nmo = len(natocc)
                ref = [2.0 if i < nref else 0.0 for i in range(nmo)]
            print('ref', ref)
            ref = np.array(ref)
            dm_ref = einsum('ij,j,kj -> ik', natorb, ref, natorb)
            if self.mcxc == 'CS':
                mcxc = 'MGGA_C_CS'
            else:
                mcxc = self.mcxc
            n, exc, excf = sudft.get_exc(ni, mol, grids, 'HF,%s'%mcxc, dm, trunc='f', gamma=dm, dmref=dm_ref)
        if self.trunc == 'fc':
            if self.ncore is None:
                core = [2.0 if occ > 1.99 else 0.0 for occ in natocc]
            else:
                nmo = len(natocc)
                core = [2.0 if i < self.ncore else 0.0 for i in range(nmo)]
            print('core', core)
            core = np.array(core)
            dm_core = einsum('ij,j,kj -> ik', natorb, core, natorb)
            #n1, exc1, core1 = get_exc(ni, self.suhf.mol, ks.grids, 'HF,%s'%suxc, dm_core, trunc='f', dmref=dm)
            n2, exc2, dE_fc = sudft.get_exc(ni, mol, grids, 'HF,%s'%mcxc, dm_core, trunc='fc', gamma=dm, dmref=dm_ref, dmcore=dm_core)
            #print('core1 %f, core2 %f' % (core1, core2))
            #dE_fc = core1 - core2
            print('dE_fc: %f' % dE_fc)
            excfc = excf + dE_fc

        E_mcdft = E_mc + exc
        print("E(CAS) = %15.8f" % E_mc)
        print("E_c(%s) = %15.8f" % (self.mcxc.upper(), exc))
        print("E(CAS+DFT) = %15.8f" % E_mcdft)
        if self.trunc[0] == 'f':
            E_mcfdft = E_mc + excf
            print("f E_c(%s) = %15.8f" % (self.mcxc.upper(), excf))
            print("E(CAS+fDFT) = %15.8f" % E_mcfdft)
        if self.trunc == 'fc':
            E_fcdft = E_mc + excfc
            print("fc E_c(%s) = %15.8f" % (self.mcxc.upper(), excfc))
            print("E(CAS+fcDFT) = %15.8f" % E_fcdft)
        t2 = time.time()
        print('time for DFT: %.3f' % (t2-t1))
        return exc, E_mcdft