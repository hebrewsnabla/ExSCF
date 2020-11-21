from pyphf import util
from pyscf import dft

import numpy as np


class SUDFT():
    def __init__(self, suhf):
        #suhf = util.SUHF(guesshf)
        self.suhf = suhf
        self.suxc = 'tpss'
        self.output = None

    def kernel(self):
        self.suhf.kernel()
        #if self.output is not None:
        #    sys.stdout = open(self.output, 'a')
        print('***** Start DFT Correlation for SUHF+DFT **********')
        E_suhf = self.suhf.E_suhf
        dm_ortho = self.suhf.dm_ortho
        X = self.suhf.X
        dm_reg = np.einsum('ij,tjk,lk->til', X, dm_ortho, X)

        ks = dft.UKS(self.suhf.mol)
        ni = ks._numint
        n, exc, vxc = ni.nr_uks(self.suhf.mol, ks.grids, 'HF,%s'%self.suxc, dm_reg)
        E_sudft = E_suhf + exc
        print("E(SUHF) = %15.8f" % E_suhf)
        print("E_c(%s) = %15.8f" % (self.suxc.upper(), exc))
        print("E(SUHF+DFT) = %15.8f" % E_sudft)
