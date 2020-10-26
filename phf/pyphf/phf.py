from pyscf import scf

def convert_to_suhf(mf):
    


class SUHF(scf.UHF):
    def __init__(self):
        return self




def build(mf):
    # S transform
    # get dm
    #
    # grid
    #
    # wigner d
    return mf


def kernel(mf):
    #
    while(true):
        get_fock()
        # diag dm , get dm in NO
        get_OmegaMatrices() # Ng, Pg, Gg
        get_xOmega() # Smt, Hmt
        solv_CI()
        get_S2()
        get_YOmega()
        get_Feff()
        get_E()
        Diag_Feff()
        conv = ConvTest()


def get_fock():
    # get coreH in orthomormal AO
    # X! . H . X