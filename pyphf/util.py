import numpy as np
import sympy as sym
import scipy
from pyscf import gto, scf
from fch2py import fch2py


def guess_from_fchk(xyz, bas, fch):
    mol = gto.Mole()
    #mol.atom = '''H 0. 0. 0.; H 0. 0. 2.'''
    with open(xyz, 'r') as f:
        mol.atom = f.read()
    print(mol.atom)
    mol.basis = bas
    mol.output = 'test.pylog'
    mol.verbose = 4
    mol.build()
    
    mf = scf.UHF(mol)
    #mf.init_guess = '1e'
    mf.init_guess_breaksym = True
    mf.max_cycle = 1
    mf.kernel()
    
    # read MOs from .fch(k) file
    nbf = mf.mo_coeff[0].shape[0]
    nif = mf.mo_coeff[0].shape[1]
    S = mol.intor_symmetric('int1e_ovlp')
    Sdiag = S.diagonal()
    alpha_coeff = fch2py(fch, nbf, nif, Sdiag, 'a')
    beta_coeff  = fch2py(fch, nbf, nif, Sdiag, 'b')
    mf.mo_coeff = (alpha_coeff, beta_coeff)
    # read done
    
    dm = mf.make_rdm1()
    mf.max_cycle = 0
    mf.kernel(dm)
    return mf

def get_beta(n):
    # Gauss-Legendre quadrature
    PI = np.pi
    b0,b1 = 0.0, PI
    m = (n+1)//2
    bm = 0.5*(b1+b0)
    bl = 0.5*(b1-b0)
    grid = np.zeros(n)
    weight = np.zeros(n)
    for i in range(m):
        z = np.cos(PI*(i+1-0.25) / (n + 0.5))
        while(True):
            p1, p2 = 1.0, 0.0
            for j in range(n):
                p3 = p2
                p2 = p1
                p1 = ((2*j + 1)*z*p2 - (j)*p3)/(j+1)
            pp = n*(z*p1 - p2)/(z*z-1)
            z = z - p1/pp
            if (abs(p1/pp) < 3e-14): 
                break
        grid[i] = bm - bl*z
        grid[n-1-i] = bm + bl*z
        weight[i] = 2*bl / ((1 - z*z)*pp*pp)
        weight[n-1-i] = weight[i]
        # Note: We did not perform weight *= sin(beta) here. That's left in get_xg().
    return grid, weight

def WignerSmall(j,m):
    #j = sym.symbols('j')
    #m = sym.symbols('m')
    j = sym.sympify(j)/2
    m = sym.sympify(m)/2
    s = sym.symbols('s')
    beta = sym.symbols('beta')
    if j==0:
        d_simp = 1
    else:
        d = (-1)**s * sym.binomial(j+m,s) * sym.binomial(j-m,s) * sym.cos(beta/2)**(2*j-2*s) * sym.sin(beta/2)**(2*s)
        upper = min(j-m,j+m)
        d = sym.Sum(d, (s, 0, upper)).doit()
        d_simp = sym.trigsimp(d)
    f = sym.lambdify(beta, d_simp, 'numpy')
    return d_simp, f

def count0(vals):
    non0 = 0
    for i in vals:
        if abs(i)>1e-10:
            non0 += 1
    return non0

def eig(A):
    return np.linalg.eigh(A)
    #return scipy.linalg.eigh(A)

def find_NO(suhf, dm, na, nb):
    cut_no = suhf.cut_no
    #dm = dm*(-1)
    #np.set_printoptions(precision=16, linewidth=200, suppress=False)
    #print(dm)
    ev_a, v_a = eig(dm[0]*(-1))
    ev_b, v_b = eig(dm[1]*(-1))
    pa = count0(ev_a)
    pb = count0(ev_b)
    print('NO eigenvalue')
    print(ev_a, ev_b)
    v_a1 = v_a[:,:na]
    v_a2 = v_a[:,na:]
    v_b1 = v_b[:,:nb]
    v_b2 = v_b[:,nb:]
    #print(v_a1, v_a2, v_b1, v_b2)
    v_a1 = np.vstack((v_a1, np.zeros(v_a1.shape)))
    v_a2 = np.vstack((v_a2, np.zeros(v_a2.shape)))
    v_b1 = np.vstack((np.zeros(v_b1.shape), v_b1))
    v_b2 = np.vstack((np.zeros(v_b2.shape), v_b2))
    #print(v_a1, v_a2, v_b1, v_b2)

    v = np.hstack((v_a1, v_b1, v_a2, v_b2))
    if cut_no:
        v = np.hstack((v_a1, v_b1, v_a2, v_b2))[:,:pa+pb]
    #v = np.hstack((v, np.zeros((v.shape[0], v.shape[0]-pa-pb))))
    print('NO vec')
    print(v)
    dm_expd = np.hstack(
        (np.vstack((dm[0], np.zeros(dm[0].shape))), 
        np.vstack((np.zeros(dm[1].shape), dm[1])))
        )
    #print(dm_expd)
    dm_no = np.einsum('ji,jk,kl->il', v, dm_expd, v)
    print('dm(NO)')
    print(dm_no)
    #np.set_printoptions(precision=6, linewidth=160, suppress=True)
    return dm_no, dm_expd, v

def get_Ng(grids, no, dm, occ):
    # for debugging:
    #grids = grids[0:1]
    Dg = []
    Ng = []
    Pg = []
    #print('N(g)')
    for beta in grids:
        norb = int(len(no)/2)
        dg1 = np.eye(norb) * np.cos(beta/2)
        dg2 = np.eye(norb) * (-np.sin(beta/2))
        dg = np.vstack((
            np.hstack((dg1, -dg2)),
            np.hstack((dg2, dg1))
        ))
        #print(dg.shape)
        #print(dg)
        dg_no = np.einsum('ji,jk,kl->il', no, dg, no)
        #print(dg_no)
        Dg.append(dg_no)
        tg = np.einsum('ij,jk,kl->il', dm[:occ,:], dg_no, dm[:,:occ])
        ng = np.linalg.inv(tg)
        #print(ng)
        Ng.append(ng)
        pg = np.einsum('ij,jk,kl,lm->im', dg_no, dm[:,:occ], ng, dm[:occ,:])
        #print(pg)
        Pg.append(pg)
    print('D(g) (NO)')
    print(Dg[0])
    print('...')
    print('N(g) (NO)')
    print(Ng[0])
    print('...')
    print('P(g) (NO)')
    #for p in Pg: print(p)
    print(Pg[0])
    print('...')
    return Dg, Ng, Pg

def get_Gg(mol, Pg, no, X):
    Gg = []
    Pg_ortho = []
    for pg in Pg:
        pg_ortho = np.einsum('ij,jk,lk->il', no, pg, no)
        #print(pg_ortho)
        Pg_ortho.append(pg_ortho)
    print('Pg_ortho')
    print(Pg_ortho[0])
    norb = int(Pg_ortho[0].shape[0]/2)
    for pg in Pg_ortho:
        pgaa = pg[:norb, :norb] # ortho ao
        #print(pgaa)
        pgab = pg[:norb, norb:]
        pgba = pg[norb:, :norb]
        pgbb = pg[norb:, norb:]
        # X . P(g) . X^H
        pgaa_ao = np.einsum('ij,jk,lk->il', X, pgaa, X) # regular ao
        #print(pgaa_ao)
        pgab_ao = np.einsum('ij,jk,lk->il', X, pgab, X)
        pgba_ao = np.einsum('ij,jk,lk->il', X, pgba, X)
        pgbb_ao = np.einsum('ij,jk,lk->il', X, pgbb, X)
        ggaa_ao, ggbb_ao = scf.uhf.get_veff(mol, [pgaa_ao, pgbb_ao], hermi=0)
        #print(ggaa_ao)
        ggab_ao = scf.hf.get_jk(mol, pgab_ao, hermi=0)[1] *(-1)
        ggba_ao = scf.hf.get_jk(mol, pgba_ao, hermi=0)[1] *(-1)
        #ggbb_ao = scf.uhf.get_veff(mol, [pgaa_ao, pgbb_ao], hermi=0)[1]
        # X^H . G(g) . X
        ggaa = np.einsum('ji,jk,kl->il', X, ggaa_ao, X)  # ortho ao
        #print(ggaa)
        ggab = np.einsum('ji,jk,kl->il', X, ggab_ao, X) 
        ggba = np.einsum('ji,jk,kl->il', X, ggba_ao, X) 
        ggbb = np.einsum('ji,jk,kl->il', X, ggbb_ao, X) 
        gg = np.vstack((
            np.hstack((ggaa, ggab)),
            np.hstack((ggba, ggbb))
        ))
        gg_no = np.einsum('ji,jk,kl->il', no, gg, no)
        Gg.append(gg_no)
    print('G(g) (NO)' )
    print(Gg[0])
    print('...')
    return Gg, Pg_ortho

def get_xg(suhf, no, na, nb, Ng):
    C_a, C_b = suhf.mo_ortho
    C_a1 = C_a[:,:na]
    C_a2 = C_a[:,na:]
    C_b1 = C_b[:,:nb]
    C_b2 = C_b[:,nb:]
    C_org = np.hstack((
        np.vstack((C_a1, np.zeros(C_a1.shape))),
        np.vstack((np.zeros(C_b1.shape), C_b1)),
        np.vstack((C_a2, np.zeros(C_a2.shape))),
        np.vstack((np.zeros(C_b2.shape), C_b2))
    ))
    #print(C_org)
    C_no = np.einsum('ji,jk->ik', no, C_org)
    print('C(NO)')
    print(C_no)
    occ = na+nb
    C_oo = C_no[:occ, :occ]
    detC = np.linalg.det(C_oo)
    print('detC', detC)
    detNg = []
    for ng in Ng:
        detng = np.linalg.det(ng)
        print(detng)
        detNg.append(detng)
    detNg = np.array(detNg)
    xg = 1.0 / (detC * detNg * detC)
    print('xg', xg)
    #sinbeta = np.sin(grids)
    #ciS = np.einsum('i,i,i,i->', weights, sinbeta, xg, d)
    ciS = suhf.integr_beta(xg)
    #print(weights, xg, d)
    print('ciS', ciS)
    yg = xg / ciS
    return xg, yg, ciS

def integr_beta(q, d, grids, weights, fac='normal', xg=None, ciS=None):
    sinbeta = np.sin(grids)
    if fac=='xg':
        #print(xg[0], d[0], 1/ciS, weights[0]*sinbeta[0])
        #print(weights)
        weights = weights * xg / ciS
        #print(xg, ciS, xg/ciS**2)
        #weights *= (xg / ciS**2)
        #print(weights)
    if q.ndim==1:
        int_q = np.einsum('i,i,i,i->', weights, sinbeta, q, d)
    elif q.ndim==2: 
        int_q = np.einsum('i,i,ij,i->j', weights, sinbeta, q, d)
    elif q.ndim==3:
        int_q = np.einsum('i,i,ijk,i->jk', weights, sinbeta, q, d)
    else:
        raise(ValueError, 'q must have dimension 1, 2, or 3')
    return int_q


def get_H(suhf, hcore_ortho, no, Pg, Gg, xg):
    #print(hcore_ortho)
    #hcore_ortho = suhf.hcore_ortho
    hcore_ortho = np.vstack((
        np.hstack((hcore_ortho, np.zeros(hcore_ortho.shape))),
        np.hstack((np.zeros(hcore_ortho.shape), hcore_ortho))
    ))
    hcore_no = np.einsum('ji,jk,kl->il', no, hcore_ortho, no)
    suhf.hcore_no = hcore_no
    print(hcore_no)
    trHg = np.zeros(len(Pg))
    for i, pg in enumerate(Pg):
        H = np.trace(np.dot(hcore_no, pg)) + 0.5 * np.trace(np.dot(Gg[i], pg))
        #H = H * xg[i]
        trHg[i] = H
        print(i, H*xg[i])
    #sinbeta = np.sin(grids)
    #ciH = np.einsum('i,i,i,i->', weights, sinbeta, Hg, d)
    ciH = suhf.integr_beta(trHg*xg)
    print('ciH', ciH)
    return trHg, ciH

def get_S2(suhf, Pg_ortho):
    norb = int(Pg_ortho[0].shape[0]/2)
    S2g = np.zeros(len(Pg_ortho))
    for i,pg in enumerate(Pg_ortho):
        pgaa = pg[:norb, :norb] # ortho ao
        #print(pgaa)
        pgab = pg[:norb, norb:]
        pgba = pg[norb:, :norb]
        pgbb = pg[norb:, norb:]
        Pc = 0.5 * (pgaa + pgbb)
        Mz = 0.5 * (pgaa - pgbb)
        My = -0.5j * (pgab - pgba)
        #print(My)
        Mx = 0.5 * (pgab + pgba)
        #print(Pc)
        trPc, trMx, trMy, trMz = list(map(np.trace, [Pc, Mx, My, Mz]))
        print(trPc, trMx, trMy, trMz)
        Pc2 = np.dot(Pc, Pc)
        Mz2 = np.dot(Mz, Mz)
        My2 = np.dot(My, My).real
        Mx2 = np.dot(Mx, Mx)
        trPc2, trMx2, trMy2, trMz2 = list(map(np.trace, [Pc2, Mx2, My2, Mz2]))
        print(trPc2, trMx2, trMy2, trMz2)
        S2g[i] = trMx**2 + (trMy**2).real + trMz**2 + 0.5*(trMx2 + trMy2 + trMz2) + 1.5*(trPc - trPc2)
    S2 = suhf.integr_beta(S2g, fac='xg')
    print('S2 = %.6f'% S2)
    return S2

def get_Yg(suhf, Dg, Ng, dm_no, occ):
    norb = len(dm_no[0])
    vir = norb - occ
    Xg = []
    for i, ng in enumerate(Ng):
        Xgi1 = np.einsum('ij,jk,kl->il', Dg[i], dm_no[:,:occ], Ng[i]) 
        Xgi2 = np.einsum('ij,jk,kl->il', Ng[i], dm_no[:occ,:], Dg[i])
        #print(Xgi1)
        #print(Xgi2)
        Xgi = np.hstack((Xgi1, np.zeros((norb, vir)))) \
            + np.vstack((Xgi2, np.zeros((vir, norb)))) 
        #print(Xgi)
        Xg.append(Xgi)
    print('X(g)')
    print(Xg[0])
    Xg = np.array(Xg)
    Xg_int = suhf.integr_beta(Xg, fac='xg')
    print('X(g) int')
    print(Xg_int)
    Yg = Xg - Xg_int
    #print(Yg.shape)
    return Xg, Xg_int, Yg

def get_Feff(suhf, trHg, Gg, Ng, Pg, dm_no, Dg, occ, Yg, Xg,  no, F_ortho, dm_ortho):
    hcore_no = suhf.hcore_no
    norb = len(Pg[0])
    vir = norb - occ
    Feff_g = [] 
    Feff0 = []
    for i, pg in enumerate(Pg):
        feff0 = Yg[i] * trHg[i]
        #print(feff0)
        fg = hcore_no + Gg[i]
        feff1 = np.einsum('ij,jk,kl,lm,mn->in', Ng[i], dm_no[:occ,:], fg, np.eye(norb)-pg, Dg[i])
        #print(feff1)
        feff1 = np.vstack((feff1, np.zeros((vir, norb))))
        feff2 = np.einsum('ij,jk,kl,lm,mn->in', np.eye(norb)-pg, fg, Dg[i], dm_no[:, :occ], Ng[i])
        #print(feff2)
        feff2 = np.hstack((feff2, np.zeros((norb, vir))))
        feff = feff0 + feff1 + feff2
        Feff_g.append(feff)
        #Feff0.append(feff0)
        #print(feff)
    Feff_g = np.array(Feff_g)
    #Feff0 = suhf.integr_beta(np.array(Feff0), fac='xg')
    #print(Feff0)
    Feff = suhf.integr_beta(Feff_g, fac='xg')
    print('Feff')
    print(Feff)
    Feff_ortho = np.einsum('ij,jk,lk->il', no, Feff, no)
    Feff_ortho = 0.5 * (Feff_ortho + Feff_ortho.T)
    print('Feff (ortho)')
    print(Feff_ortho)
    H = suhf.integr_beta(trHg, fac='xg')
    print('Hsp + Hph = ', H)
    #F_ortho = np.vstack((
    #    np.hstack((F_ortho[0], np.zeros(F_ortho[0].shape))),
    #    np.hstack((np.zeros(F_ortho[1].shape), F_ortho[1]))
    #))

    #F_no = np.einsum('ji,jk,kl->il', no, F_ortho, no)
    #Feff_no = np.einsum('ji,jk,kl->il', no, Feff_ortho, no)
    #F_mod = np.vstack((
    #    np.hstack((F_no[:occ, :occ], Feff_no[:occ, occ:])),
    #    np.hstack((Feff_no[occ:, :occ], F_no[occ:, occ:]))
    #))
    #F_mod_ortho = np.einsum('ij,jk,lk->il', no, F_mod, no)
    
    #dm_ortho = np.vstack((
    #    np.hstack((dm_ortho[0], np.zeros(dm_ortho[0].shape))),
    #    np.hstack((np.zeros(dm_ortho[1].shape), dm_ortho[1]))
    #))
    lenf = len(F_ortho[0])
    F_ortho = np.array(F_ortho)
    Feff_ortho = np.array([Feff_ortho[:lenf, :lenf], Feff_ortho[lenf:,lenf:]])
    Q = np.eye(len(dm_ortho[0])) - dm_ortho
    Foo = np.einsum('tij,tjk,tkl->til', dm_ortho, F_ortho, dm_ortho)
    Fvv = np.einsum('tij,tjk,tkl->til', Q, F_ortho, Q)
    Fov = np.einsum('tij,tjk,tkl->til', dm_ortho, Feff_ortho, Q)
    Fvo = np.einsum('tij,tjk,tkl->til', Q, Feff_ortho, dm_ortho)
    #F_mod_ortho_a = np.vstack((np.hstack((Foo[0], Fov[0])), np.hstack((Fvo[0], Fvv[0]))))
    #F_mod_ortho_b = np.vstack((np.hstack((Foo[1], Fov[1])), np.hstack((Fvo[1], Fvv[1]))))
    F_mod_ortho_a = Foo[0] + Fov[0] + Fvo[0] + Fvv[0]
    F_mod_ortho_b = Foo[1] + Fov[1] + Fvo[1] + Fvv[1]
    F_mod_ortho_a = 0.5*(F_mod_ortho_a + F_mod_ortho_a.T)
    F_mod_ortho_b = 0.5*(F_mod_ortho_b + F_mod_ortho_b.T)
    F_mod_ortho = np.vstack((
        np.hstack((F_mod_ortho_a, np.zeros(F_mod_ortho_a.shape))),
        np.hstack((np.zeros(F_mod_ortho_b.shape), F_mod_ortho_b))
    ))
    print('Feff (mod,ortho)')
    print(F_mod_ortho)
    return Feff_ortho, H, F_mod_ortho

def Diag_Feff(Feff_ortho, na, nb):
    norb = int(len(Feff_ortho)/2)
    Faa = Feff_ortho[:norb, :norb]
    Fbb = Feff_ortho[norb:, norb:]
    e_a, v_a = np.linalg.eigh(Faa)
    e_b, v_b = np.linalg.eigh(Fbb)
    P_a = np.einsum('ij,kj->ik', v_a[:,:na], v_a[:,:na])
    P_b = np.einsum('ij,kj->ik', v_b[:,:nb], v_b[:,:nb])
    print('e_a, e_b')
    print(e_a, e_b)
    print('P_a, P_b')
    print(P_a, P_b)
    return [e_a, e_b], [P_a, P_b]

class SUHF():
    '''
    Attributes:
        guesshf: UHF object
        X: 
            transformation matrix (regular AO -> orthonormal AO)
            It's asymmetric, as X in G09
        S, Sz:
        nbeta: number of grids for beta
        grids: integration grids for beta
        weights:
        d_expr: SymPy expr for Wigner d 
        d_func: NumPy function for Wigner d
        d     : array
            current values of Wigner d, on given grids
    '''

    def __init__(self, guesshf):
        self.guesshf = guesshf
        self.mol = guesshf.mol
        S = guesshf.get_ovlp()
        Ca, Cb = guesshf.mo_coeff
        #S_sqrt = scipy.linalg.sqrtm(S)
        print(S)
        #Se, Svec = scipy.linalg.eigh(S)
        #Se_msq = Se.real**(-0.5)
        #print(Se_msq)
        #S_msq2 = np.einsum('ji,j,jk->ik',Svec,Se_msq, Svec)
        #S_msq = scipy.linalg.fractional_matrix_power(S, -0.5)
        print('SVD: S^(-1/2)')
        
        Su, Ss, Sv = scipy.linalg.svd(S)
        #print('SVD')
        #print(Su,Ss,Sv)
        X = np.einsum('ij,j->ij', Su, Ss**(-0.5))
        #Xd = np.einsum('i,ij->ij', Ss**(-0.5), Sv)
        print(X)
        self.X = X
        #print(Xd)
        print('S^(1/2)')
        XS = np.dot(X.T,S)
        print(XS)
        
        print('C (ortho)')
        Ca_ortho = np.dot(XS, Ca)
        Cb_ortho = np.dot(XS, Cb)
        print(Ca_ortho)
        print(Cb_ortho)
        
        self.mo_ortho = Ca_ortho, Cb_ortho
        self.dm_ortho = scf.uhf.make_rdm1(self.mo_ortho, guesshf.mo_occ)
        print('density matrix (ortho)')
        print(self.dm_ortho)

        self.nbeta = 8
        self.grids, self.weights = get_beta(self.nbeta)
        print(self.grids, self.weights)
        spin = self.mol.spin
        na, nb = guesshf.nelec
        self.nelec = na, nb
        sz = (na-nb)/2
        print('S = %.1f, Sz = %.1f' % (spin/2, sz))
        self.S, self.Sz = S, sz
        Wignerd_expr, Wignerd = WignerSmall(int(spin), int(2*sz))
        print('Wigner small d: ', Wignerd_expr)
        d = np.zeros(self.nbeta)
        for g in range(self.nbeta):
            d[g] = Wignerd(self.grids[g])
        print('value :', d)
        self.d_expr, self.d_func, self.d = Wignerd_expr, Wignerd, d
        self.E_suhf = None

    def integr_beta(self, q, fac='normal'):
        if fac=='xg':
            return integr_beta(q, self.d, self.grids, self.weights, 'xg', self.xg, self.ciS)
        else:
            return integr_beta(q, self.d, self.grids, self.weights)

    
    def kernel(self):
        np.set_printoptions(precision=8, linewidth=160, suppress=True)
        if self.debug:
            np.set_printoptions(precision=15, linewidth=200, suppress=False)
        X = self.X
        na, nb = self.nelec
        mf = self.guesshf
        
        max_cycle = self.max_cycle
        cyc = 0
        conv = False
        
        hcore = mf.get_hcore()
        hcore_ortho = np.einsum('ji,jk,kl->il', X, hcore, X)
        while(not conv):
            print('**** Cycle %d ****' % (cyc+1))
            old_suhf = self.E_suhf
            #print(hcore_ortho)
            
            #if cyc==0:
            #    veff = mf.get_veff(dm = dm)
            #else:
            dm_reg = np.einsum('ij,tjk,lk->til', X, self.dm_ortho, X)
            veff = mf.get_veff(dm = dm_reg)
            veff_ortho = np.einsum('ji,tjk,kl->til', X, veff, X)
            print('dm (ortho)')
            print(self.dm_ortho)
            #print(veff)
            #Fa, Fb = hcore + veff
            #Fa_ortho = np.einsum('ji,jk,kl->il', X, Fa, X)
            #Fb_ortho = np.einsum('ji,jk,kl->il', X, Fb, X)
            Fa_ortho, Fb_ortho = hcore_ortho + veff_ortho
            print('Fock (ortho)',Fa_ortho, Fb_ortho)
            F_ortho = Fa_ortho, Fb_ortho
        
            e_uhf, e_uhf_coul = scf.uhf.energy_elec(mf, self.dm_ortho, hcore_ortho, veff_ortho)
            print('E(UHF) = %12.6f' % e_uhf)
        
            dm_no, dm_expanded, no = find_NO(self, self.dm_ortho, na, nb)
            Dg, Ng, Pg = get_Ng(self.grids, no, dm_no, na+nb)
            Gg, Pg_ortho = get_Gg(self.mol, Pg, no, X)
            xg, yg, ciS = get_xg(self, no, na, nb, Ng)
            self.xg, self.ciS = xg, ciS
            #yg, ciS = util.get_yg(self, xg)
            trHg, ciH = get_H(self, hcore_ortho, no, Pg, Gg, xg)
            S2 = get_S2(self, Pg_ortho)
            Xg, Xg_int, Yg = get_Yg(self, Dg, Ng, dm_no, na+nb)
            Feff_ortho, H_suhf, F_mod_ortho = get_Feff(self, trHg, Gg, Ng, Pg, dm_no, Dg, na+nb, Yg, Xg, no, F_ortho, self.dm_ortho)
            E_suhf = mf.energy_nuc() + H_suhf
            self.E_suhf = E_suhf
            print('E(SUHF) = %15.8f' % E_suhf)
            mo_e, self.dm_ortho = Diag_Feff(F_mod_ortho, na, nb)
        
            if old_suhf is not None:
                if abs(E_suhf - old_suhf)<1e-8:
                    conv = True
                    print('\n***************')
                    print('SCF converged at cycle %d' %cyc)
                    print('Final E(SUHF) = %15.8f, delta E = %10.6g' % (E_suhf,E_suhf-old_suhf))
                    
            
            cyc += 1
            if cyc >= max_cycle:
                break


        