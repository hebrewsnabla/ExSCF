from pyphf import util
from pyscf import dft
import pyscf.dft.numint as numint

import numpy as np
from functools import partial
import time

print = partial(print, flush=True)
einsum = partial(np.einsum, optimize=True)


class CASDFT():
    def __init__(self, suhf):
        #suhf = util.SUHF(guesshf)
        self.mc = mc
        self.mcxc = 'tpss'
        self.output = None
        #self.dens = 'deformed' # or relaxed
        self.trunc = None

    def kernel(self):
        #if self.suhf.E_suhf is None:
        #    self.suhf.kernel()
        #if self.output is not None:
        #    sys.stdout = open(self.output, 'a')
        print('***** Start DFT Correlation for CAS+DFT **********')
        #print('density: %s' % self.dens)
        print('truncation: %s' % self.trunc)
        t1 = time.time()
        E_mc = self.mc.e_tot
        dm = self.mc.make_rdm1()

        ks = dft.UKS(self.suhf.mol)
        ni = ks._numint
        if self.trunc is None:
            if self.mcxc == 'CS':
                mcxc = 'MGGA_C_CS'
                n, exc = get_exc(ni, self.suhf.mol, ks.grids, 'HF,%s'%mcxc, dm)
            else:
                n, exc, vxc = ni.nr_uks(self.suhf.mol, ks.grids, 'HF,%s'%self.mcxc, dm)
        elif self.trunc == 'f':
            natorb = self.mc.natorb
            natocc = self.mc.mo_occ
            #natocc = natocc[0] + natocc[1]
            print('natocc', natocc)
            ref = [2.0 if occ > 1e-2 else 0.0 for occ in natocc]
            print('ref', ref)
            ref = np.array(ref)
            dm_ref = einsum('ij,j,kj -> ik', natorb, ref, natorb)
            if self.mcxc == 'CS':
                mcxc = 'MGGA_C_CS'
            else:
                mcxc = self.mcxc
            n, exc, excf = get_exc(ni, self.suhf.mol, ks.grids, 'HF,%s'%mcxc, dm, trunc='f', dmref=dm_ref)

        E_mcdft = E_mc + exc
        print("E(CAS) = %15.8f" % E_mc)
        print("E_c(%s) = %15.8f" % (self.mcxc.upper(), exc))
        print("E(CAS+DFT) = %15.8f" % E_mcdft)
        if self.trunc == 'f':
            E_mcfdft = E_mc + excf
            print("f E_c(%s) = %15.8f" % (self.mcxc.upper(), excf))
            print("E(CAS+fDFT) = %15.8f" % E_mcfdft)
        t2 = time.time()
        print('time for DFT: %.3f' % (t2-t1))
        return exc, E_mcdft

class SUDFT():
    def __init__(self, suhf):
        #suhf = util.SUHF(guesshf)
        self.suhf = suhf
        self.suxc = 'tpss'
        self.output = None
        self.dens = 'deformed' # or relaxed
        self.trunc = None

    def kernel(self):
        if self.suhf.E_suhf is None:
            self.suhf.kernel()
        #if self.output is not None:
        #    sys.stdout = open(self.output, 'a')
        print('***** Start DFT Correlation for SUHF+DFT **********')
        print('density: %s' % self.dens)
        print('truncation: %s' % self.trunc)
        t1 = time.time()
        E_suhf = self.suhf.E_suhf
        #dm_ortho = self.suhf.dm_ortho
        #X = self.suhf.X
        #dm_reg = np.einsum('ij,tjk,lk->til', X, dm_ortho, X)
        if self.dens[:3] == 'def':
            dm = self.suhf.dm_reg
        elif self.dens[:3] == 'rel':
            dm = self.suhf.suhf_dm

        ks = dft.UKS(self.suhf.mol)
        ni = ks._numint
        if self.trunc is None:
            if self.suxc == 'CS':
                suxc = 'MGGA_C_CS'
                n, exc = get_exc(ni, self.suhf.mol, ks.grids, 'HF,%s'%suxc, dm)
            else:
                n, exc, vxc = ni.nr_uks(self.suhf.mol, ks.grids, 'HF,%s'%self.suxc, dm)
        elif self.trunc == 'f' or self.trunc == 'fc':
            natorb = self.suhf.natorb[2]
            natocc = self.suhf.natocc[2]
            #natocc = natocc[0] + natocc[1]
            print('natocc', natocc)
            ref = [2.0 if occ > 1e-2 else 0.0 for occ in natocc]
            print('ref', ref)
            ref = np.array(ref)
            dm_ref = einsum('ij,j,kj -> ik', natorb, ref, natorb)
            if self.suxc == 'CS':
                suxc = 'MGGA_C_CS'
            else:
                suxc = self.suxc
            n, exc, excf = get_exc(ni, self.suhf.mol, ks.grids, 'HF,%s'%suxc, dm, trunc='f', dmref=dm_ref)
        if self.trunc == 'fc':
            core = [2.0 if occ > 1.98 else 0.0 for occ in natocc]
            print('core', core)
            core = np.array(core)
            dm_core = einsum('ij,j,kj -> ik', natorb, core, natorb)
            n1, exc1, core1 = get_exc(ni, self.suhf.mol, ks.grids, 'HF,%s'%suxc, dm_core, trunc='f', dmref=dm)
            n2, exc2, core2 = get_exc(ni, self.suhf.mol, ks.grids, 'HF,%s'%suxc, dm_core, trunc='f', dmref=dm_ref)
            dE_fc = core1 - core2
            excfc = excf + dE_fc

        E_sudft = E_suhf + exc
        print("E(SUHF) = %15.8f" % E_suhf)
        print("E_c(%s) = %15.8f" % (self.suxc.upper(), exc))
        print("E(SUHF+DFT) = %15.8f" % E_sudft)
        if self.trunc == 'f' or self.trunc == 'fc':
            E_sufdft = E_suhf + excf
            print("f E_c(%s) = %15.8f" % (self.suxc.upper(), excf))
            print("E(SUHF+fDFT) = %15.8f" % E_sufdft)
        if self.trunc == 'fc':
            E_sufcdft = E_suhf + excfc
            print("fc E_c(%s) = %15.8f" % (self.suxc.upper(), excfc))
            print("E(SUHF+fcDFT) = %15.8f" % E_sufcdft)
        
        t2 = time.time()
        print('time for DFT: %.3f' % (t2-t1))
        return exc, E_sudft

def get_exc(ni, mol, grids, xc_code, dms, trunc=None, dmref=None, 
            relativity=0, hermi=0, max_memory=2000, verbose=None):
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
    make_rho_ref    = ni._gen_rho_evaluator(mol, dmref, hermi)[0]

    nelec = np.zeros((2,nset))
    excsum = np.zeros(nset)
    excfsum = np.zeros(nset)
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
                rho_ref = make_rho_ref(idm, ao, mask, xctype)
                exc, vxc = ni.eval_xc(xc_code, (rho_a, rho_b), spin=1,
                                      relativity=relativity, deriv=1,
                                      verbose=verbose)[:2]
                vrho, vsigma, vlapl, vtau = vxc[:4]
                den_a = rho_a[0]*weight
                nelec[0,idm] += den_a.sum()
                excsum[idm] += np.dot(den_a, exc)

                den_b = rho_b[0]*weight
                nelec[1,idm] += den_b.sum()
                excsum[idm] += np.dot(den_b, exc)
                rhoall = rho_a[0] + rho_b[0]
                eta = (rho_ref[0] / rhoall)**(1.0/3)
                for i in range(len(rhoall)):
                    if rhoall[i] < 1e-45:
                        rhoall[i] = 1e-45
                        eta[i] = 1.0
                #print('f', f.shape)
                f = get_f(rhoall, eta)
                excfsum[idm] += einsum('i,i,i->', den_a+den_b, exc, f)

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
    if trunc is not None:
        return nelec, excsum , excfsum
    else:
        return nelec, excsum 

def get_f(rho, eta):
    b = [[-2.207193,     6.807648,     -6.386316,     2.860522,     -0.07466076],
         [ 1.128469,    -2.535669,      2.432821,    -1.064058,      0.03843687],
         [-0.2475593,    0.4243142,    -0.2565175,    0.08294749,   -0.003184296],
         [ 0.08616560,  -0.1715714,     0.1067547,   -0.02392882,    0.002579856],
         [-0.6500077e-2, 0.1714085e-1, -0.1462187e-1, 0.4423830e-2, -0.4427570e-3],
         [-0.2491486e-2, 0.5321373e-2, -0.3704699e-2, 0.9700054e-3, -0.9518308e-4]]
    x = (1.0/3) * np.log(3.0 / (4 * np.pi * rho))
    finv = np.zeros(rho.shape[0])
    for m in range(6):
        for n in range(5):
            finv += b[m][n] * x**m * eta**n
    f = finv**(-1)
    return f