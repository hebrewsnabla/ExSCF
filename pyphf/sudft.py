from pyphf import util
from pyscf import dft
import pyscf.dft.numint as numint

import numpy as np


class SUDFT():
    def __init__(self, suhf):
        #suhf = util.SUHF(guesshf)
        self.suhf = suhf
        self.suxc = 'tpss'
        self.output = None
        self.dens = 'deformed' # or relaxed

    def kernel(self):
        self.suhf.kernel()
        #if self.output is not None:
        #    sys.stdout = open(self.output, 'a')
        print('***** Start DFT Correlation for SUHF+DFT **********')
        E_suhf = self.suhf.E_suhf
        #dm_ortho = self.suhf.dm_ortho
        #X = self.suhf.X
        #dm_reg = np.einsum('ij,tjk,lk->til', X, dm_ortho, X)
        if self.dens == 'deformed':
            dm = self.suhf.dm_reg
        elif self.dens == 'relaxed':
            dm = self.suhf.suhf_dm

        ks = dft.UKS(self.suhf.mol)
        ni = ks._numint
        if self.suxc == 'CS':
            suxc = 'MGGA_C_CS'
            n, exc = get_exc(ni, self.suhf.mol, ks.grids, 'HF,%s'%suxc, dm)
        else:
            n, exc, vxc = ni.nr_uks(self.suhf.mol, ks.grids, 'HF,%s'%self.suxc, dm)
        E_sudft = E_suhf + exc
        print("E(SUHF) = %15.8f" % E_suhf)
        print("E_c(%s) = %15.8f" % (self.suxc.upper(), exc))
        print("E(SUHF+DFT) = %15.8f" % E_sudft)
        return exc, E_sudft

def get_exc(ni, mol, grids, xc_code, dms, relativity=0, hermi=0, max_memory=2000, verbose=None):
    '''
    modified from pyscf.dft.numint.nr_uks
    '''
    xctype = ni._xc_type(xc_code)
    #if xctype == 'NLC':
    #    dms_sf = dms[0] + dms[1]
    #    nelec, excsum, vmat = nr_rks(ni, mol, grids, xc_code, dms_sf, relativity, hermi,
    #                                 max_memory, verbose)
    #    return [nelec,nelec], excsum, np.asarray([vmat,vmat])

    shls_slice = (0, mol.nbas)
    ao_loc = mol.ao_loc_nr()

    dma, dmb = numint._format_uks_dm(dms)
    nao = dma.shape[-1]
    make_rhoa, nset = ni._gen_rho_evaluator(mol, dma, hermi)[:2]
    make_rhob       = ni._gen_rho_evaluator(mol, dmb, hermi)[0]

    nelec = np.zeros((2,nset))
    excsum = np.zeros(nset)
#    vmat = np.zeros((2,nset,nao,nao), dtype=np.result_type(dma, dmb))
#    aow = None
    if xctype == 'LDA':
        ao_deriv = 0
        for ao, mask, weight, coords \
                in ni.block_loop(mol, grids, nao, ao_deriv, max_memory):
#            aow = np.ndarray(ao.shape, order='F', buffer=aow)
            for idm in range(nset):
                rho_a = make_rhoa(idm, ao, mask, xctype)
                rho_b = make_rhob(idm, ao, mask, xctype)
                exc, vxc = ni.eval_xc(xc_code, (rho_a, rho_b), spin=1,
                                      relativity=relativity, deriv=1,
                                      verbose=verbose)[:2]
                vrho = vxc[0]
                den = rho_a * weight
                nelec[0,idm] += den.sum()
                excsum[idm] += np.dot(den, exc)
                den = rho_b * weight
                nelec[1,idm] += den.sum()
                excsum[idm] += np.dot(den, exc)

#                # *.5 due to +c.c. in the end
#                #:aow = np.einsum('pi,p->pi', ao, .5*weight*vrho[:,0], out=aow)
#                aow = _scale_ao(ao, .5*weight*vrho[:,0], out=aow)
#                vmat[0,idm] += _dot_ao_ao(mol, ao, aow, mask, shls_slice, ao_loc)
#                #:aow = np.einsum('pi,p->pi', ao, .5*weight*vrho[:,1], out=aow)
#                aow = _scale_ao(ao, .5*weight*vrho[:,1], out=aow)
#                vmat[1,idm] += _dot_ao_ao(mol, ao, aow, mask, shls_slice, ao_loc)
                rho_a = rho_b = exc = vxc = vrho = None
    elif xctype == 'GGA':
        ao_deriv = 1
        for ao, mask, weight, coords \
                in ni.block_loop(mol, grids, nao, ao_deriv, max_memory):
#            aow = np.ndarray(ao[0].shape, order='F', buffer=aow)
            for idm in range(nset):
                rho_a = make_rhoa(idm, ao, mask, xctype)
                rho_b = make_rhob(idm, ao, mask, xctype)
                exc, vxc = ni.eval_xc(xc_code, (rho_a, rho_b), spin=1,
                                      relativity=relativity, deriv=1,
                                      verbose=verbose)[:2]
                den = rho_a[0]*weight
                nelec[0,idm] += den.sum()
                excsum[idm] += np.dot(den, exc)
                den = rho_b[0]*weight
                nelec[1,idm] += den.sum()
                excsum[idm] += np.dot(den, exc)

#                wva, wvb = _uks_gga_wv0((rho_a,rho_b), vxc, weight)
#                #:aow = np.einsum('npi,np->pi', ao, wva, out=aow)
#                aow = _scale_ao(ao, wva, out=aow)
#                vmat[0,idm] += _dot_ao_ao(mol, ao[0], aow, mask, shls_slice, ao_loc)
#                #:aow = np.einsum('npi,np->pi', ao, wvb, out=aow)
#                aow = _scale_ao(ao, wvb, out=aow)
#                vmat[1,idm] += _dot_ao_ao(mol, ao[0], aow, mask, shls_slice, ao_loc)
                rho_a = rho_b = exc = vxc = wva = wvb = None
    elif xctype == 'MGGA':
#        if (any(x in xc_code.upper() for x in ('CC06', 'CS', 'BR89', 'MK00'))):
#            raise NotImplementedError('laplacian in meta-GGA method')
        ao_deriv = 2
        for ao, mask, weight, coords \
                in ni.block_loop(mol, grids, nao, ao_deriv, max_memory):
#            aow = np.ndarray(ao[0].shape, order='F', buffer=aow)
            for idm in range(nset):
                rho_a = make_rhoa(idm, ao, mask, xctype)
                rho_b = make_rhob(idm, ao, mask, xctype)
                exc, vxc = ni.eval_xc(xc_code, (rho_a, rho_b), spin=1,
                                      relativity=relativity, deriv=1,
                                      verbose=verbose)[:2]
                vrho, vsigma, vlapl, vtau = vxc[:4]
                den = rho_a[0]*weight
                nelec[0,idm] += den.sum()
                excsum[idm] += np.dot(den, exc)
                den = rho_b[0]*weight
                nelec[1,idm] += den.sum()
                excsum[idm] += np.dot(den, exc)

#                wva, wvb = _uks_gga_wv0((rho_a,rho_b), vxc, weight)
#                #:aow = np.einsum('npi,np->pi', ao[:4], wva, out=aow)
#                aow = _scale_ao(ao[:4], wva, out=aow)
#                vmat[0,idm] += _dot_ao_ao(mol, ao[0], aow, mask, shls_slice, ao_loc)
#                #:aow = np.einsum('npi,np->pi', ao[:4], wvb, out=aow)
#                aow = _scale_ao(ao[:4], wvb, out=aow)
#                vmat[1,idm] += _dot_ao_ao(mol, ao[0], aow, mask, shls_slice, ao_loc)

# FIXME: .5 * .5   First 0.5 for v+v.T symmetrization.
# Second 0.5 is due to the Libxc convention tau = 1/2 \nabla\phi\dot\nabla\phi
#                wv = (.25 * weight * vtau[:,0]).reshape(-1,1)
#                vmat[0,idm] += _dot_ao_ao(mol, ao[1], wv*ao[1], mask, shls_slice, ao_loc)
#                vmat[0,idm] += _dot_ao_ao(mol, ao[2], wv*ao[2], mask, shls_slice, ao_loc)
#                vmat[0,idm] += _dot_ao_ao(mol, ao[3], wv*ao[3], mask, shls_slice, ao_loc)
#                wv = (.25 * weight * vtau[:,1]).reshape(-1,1)
#                vmat[1,idm] += _dot_ao_ao(mol, ao[1], wv*ao[1], mask, shls_slice, ao_loc)
#                vmat[1,idm] += _dot_ao_ao(mol, ao[2], wv*ao[2], mask, shls_slice, ao_loc)
#                vmat[1,idm] += _dot_ao_ao(mol, ao[3], wv*ao[3], mask, shls_slice, ao_loc)
                rho_a = rho_b = exc = vxc = vrho = wva = wvb = None

#    for i in range(nset):
#        vmat[0,i] = vmat[0,i] + vmat[0,i].conj().T
#        vmat[1,i] = vmat[1,i] + vmat[1,i].conj().T
    if isinstance(dma, np.ndarray) and dma.ndim == 2:
#        vmat = vmat[:,0]
        nelec = nelec.reshape(2)
        excsum = excsum[0]
    return nelec, excsum #, vmat