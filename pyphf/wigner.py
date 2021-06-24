import sympy as sym
import numpy as np

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

def wigner(spin, sz, nbeta, grids):
    Wignerd_expr, Wignerd = WignerSmall(int(spin), int(2*sz))
    print('Wigner small d: ', Wignerd_expr)
    d = np.zeros(nbeta)
    for g in range(nbeta):
        d[g] = Wignerd(grids[g])
    print('value :', d)
    return Wignerd_expr, Wignerd, d


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
