
import os
import numpy as np
import mpmath as mp ; mp.dps = 15

from utils import continuum_dir#, plot_func



def rxVex_eval(MO_data, cnj, clj, cmj, czj, cAj, cCj, Vnpt, rxV_grid):

    rxVex = {}
    for key_j in MO_data:
        jmax = len(cnj[key_j])
        rxVex.setdefault(key_j, [])
        Vtmp = np.zeros((Vnpt), dtype=complex)
        for i in range(Vnpt):
            r = rxV_grid[i]
            sum_kl = 0.*1j
            for jk in range(jmax):
                Ck = cCj[key_j][jk]
                nk = cnj[key_j][jk]
                lk = clj[key_j][jk]
                mk = cmj[key_j][jk]
                zk = czj[key_j][jk]
                for jl in range(jmax):
                    Cl = cCj[key_j][jl]
                    nl = cnj[key_j][jl]
                    ll = clj[key_j][jl]
                    ml = cmj[key_j][jl]
                    zl = czj[key_j][jl]
                    if( (lk==ll) and (mk==ml) ):
                        sum_kl += mp.conj(Ck) * Cl * rxK_kl(nk, zk, nl, zl, r)
            Vtmp[i] = sum_kl
        rxVex[key_j].append(Vtmp)

    return(rxVex)

# ---

def rxK_kl(nk, zk, nl, zl, r):
    rz = (zk+zl)*r
    k1 = mp.gammainc(nk+nl+1., 0 , rz    ) / (zk+zl)**(nk+nl+1.)
    k2 = mp.gammainc(nk+nl   , rz, mp.inf) / (zk+zl)**(nk+nl   )
    kl = (k1+k2*r) 
    return(kl)

# ---

def rxVnuc(Mn, Zn, Rn, r):
    rxv = 0.
    for n in range(Mn):
        rxv -= Zn[n] * r / max(r, Rn[n])
    return(rxv)

# ---

def generate_coulomb(dr1, VR_min, VR_max, Z):

    Vnpt = int((VR_max - VR_min) / dr1)

    rxV_grid = np.array([VR_min + i * dr1 for i in range(Vnpt)])

    Vtmp = Z * np.ones(Vnpt, dtype=complex)
    Vtmp[0] = 0.0
    
    return(rxV_grid, Vtmp)

# ---

