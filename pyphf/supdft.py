from pyphf import util
from pyscf import dft
import pyscf.dft.numint as numint

import numpy as np
from functools import partial
import time

print = partial(print, flush=True)
einsum = partial(np.einsum, optimize=True)

class PDFT():
    def __init__(self, suhf):
        self.suhf = suhf
        self.ot = 'tpbe'
    def kernel(self):
        return kernel(self, self.suhf, self.ot)

def kernel(pdft, suhf, ot):
    dm1 = suhf.suhf_dm
    dm2 = get_2CDM_from_2RDM(suhf.suhf_dm2, )

    omega, alpha, hyb = ot._numint.rsh_and_hybrid_coeff(ot.otxc, spin=spin)

    hfx = suhf.get_EX()