import numpy as np

def stack22(aa, ab, ba, bb):
    tot = np.vstack((
        np.hstack((aa,ab)),
        np.hstack((ba,bb))
    ))
    return tot

def reg2ortho(dm, X, forward=True):
    if forward:
        return np.einsum('ji,jk,kl', X, dm, X)
    else:
        return np.einsum('ij,jk,lk', X, dm, X)
