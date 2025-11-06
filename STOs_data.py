
import os, sys

import mpmath as mp ; mp.dps = 15

from utils import bound_dir 



def get_MO():
    
    fname = bound_dir + "MO.dat"
    f     = open(fname, 'r')
    lines = f.readlines()
    f.close()

    MO_data = {}
    for line in lines:
        p = line.split()
        MO_data[str(p[0])] = [float(p[1]), float(p[2])]

    return(MO_data)

# ---

# -------------------------------------------------------------------------------------------------
def get_STOs(MO_data):

    nj = {}
    lj = {}
    mj = {}
    zj = {}
    Aj = {}
    
    for key in MO_data:

        fname = bound_dir + "orbit_" + key + ".dat" 
        f     = open(fname, 'r')
        lines = f.readlines()
        f.close()

        nj.setdefault(key, [])
        lj.setdefault(key, [])
        mj.setdefault(key, [])
        zj.setdefault(key, [])
        Aj.setdefault(key, [])

        for line in lines:
            p = line.split()
            nj[key].append( int(  p[0]) )
            lj[key].append( int(  p[1]) )
            mj[key].append( int(  p[2]) )
            zj[key].append( float(p[3]) )
            Aj[key].append( float(p[4]) )

    return(nj, lj, mj, zj, Aj)
# -------------------------------------------------------------------------------------------------

# ---

# -------------------------------------------------------------------------------------------------
def get_cSTOs(MO_data, nj, lj, mj, zj, Aj):

    cnj = {}
    clj = {}
    cmj = {}
    czj = {}
    cAj = {}
    cCj = {}
    
    for key in MO_data:

        cnj.setdefault(key, [])
        clj.setdefault(key, [])
        cmj.setdefault(key, [])
        czj.setdefault(key, [])
        cAj.setdefault(key, [])
        cCj.setdefault(key, [])

        jmax = len(nj[key])
        for j in range(jmax):

            n = nj[key][j]
            l = lj[key][j]
            m = mj[key][j]
            z = zj[key][j]
            A = Aj[key][j]

            if( m > 0 ):

                cnj[key].append(n)
                cnj[key].append(n)
                clj[key].append(l)
                clj[key].append(l)
                cmj[key].append(-m)
                cmj[key].append(+m)
                czj[key].append(z)
                czj[key].append(z)

                A1 =            A / (2.)**0.5
                A2 = (-1.)**m * A / (2.)**0.5
                cAj[key].append(A1)
                cAj[key].append(A2)

                C1 = A1 * mp.sqrt((2.*z)**(2.*n+1.)/mp.gamma(2.*n+1.))
                C2 = A2 * mp.sqrt((2.*z)**(2.*n+1.)/mp.gamma(2.*n+1.))
                cCj[key].append(C1)
                cCj[key].append(C2)

            elif( m == 0 ):

                cnj[key].append(n)
                clj[key].append(l)
                cmj[key].append(m)
                czj[key].append(z)
                cAj[key].append(A)

                C = A * mp.sqrt((2.*z)**(2.*n+1.)/mp.gamma(2.*n+1.))
                cCj[key].append(C)

            else:

                cnj[key].append(n)
                cnj[key].append(n)
                clj[key].append(l)
                clj[key].append(l)
                cmj[key].append(+m)
                cmj[key].append(-m)
                czj[key].append(z)
                czj[key].append(z)

                A1 =             1j * A / (2.)**0.5
                A2 = -(-1.)**m * 1j * A / (2.)**0.5
                cAj[key].append(A1)
                cAj[key].append(A2)

                C1 = A1 * mp.sqrt((2.*z)**(2.*n+1.)/mp.gamma(2.*n+1.))
                C2 = A2 * mp.sqrt((2.*z)**(2.*n+1.)/mp.gamma(2.*n+1.))
                cCj[key].append(C1)
                cCj[key].append(C2)

    return(cnj, clj, cmj, czj, cAj, cCj)
# -------------------------------------------------------------------------------------------------

# ---

# -------------------------------------------------------------------------------------------------
def save_cSTOs(MO_data, cnj, clj, cmj, czj, cAj, cCj):
    
    for key in MO_data:

        f_name = bound_dir + "cplx_rep/c_orbit_" + key + ".dat" 
        os.makedirs(os.path.dirname(f_name), exist_ok=True)
        with open(f_name, 'w') as f:

            jmax = len(cnj[key])
            for j in range(jmax):

                cn = cnj[key][j]
                cl = clj[key][j]
                cm = cmj[key][j]
                cz = czj[key][j]
                cA = cAj[key][j]
                f.write(' {}  {}  {:+}  {}  {}\n'.format(cn, cl, cm, cz, cA))
# -------------------------------------------------------------------------------------------------

# ---

# -------------------------------------------------------------------------------------------------
def check_norm(MO_data, cnj, clj, cmj, czj, cAj, cCj):

    for key in MO_data:

        norm_test = 0.
        jmax = len(cnj[key])
        for j1 in range(jmax):

            n1 = cnj[key][j1]
            l1 = clj[key][j1]
            m1 = cmj[key][j1]
            z1 = czj[key][j1]
            C1 = cCj[key][j1]

            for j2 in range(jmax):

                n2 = cnj[key][j2]
                l2 = clj[key][j2]
                m2 = cmj[key][j2]
                z2 = czj[key][j2]
                C2 = cCj[key][j2]

                if( (l1==l2) and (m1==m2) ):
                    norm_test += mp.conj(C1) * C2 * mp.gamma(n1+n2+1.) / (z1+z2)**(n1+n2+1.)

        print(" check normalisation for {}: {}".format(key,norm_test))
# -------------------------------------------------------------------------------------------------

# ---

