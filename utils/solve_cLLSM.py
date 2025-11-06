
import os
import numpy as np


# ---

def run_cLLSM(Rgrid, Rmax_fit, rad_func, ng, l_max, nke_max, cusp):

    Rnpt_fit = np.where(Rgrid == Rmax_fit)
    Rnpt_fit = Rnpt_fit[0][0]

    inp_cLLSM = "cLLSM.in"
    with open(inp_cLLSM, 'w') as f:
        f.write(' {} \n'.format(ng))
        f.write(' {} \n'.format(Rnpt_fit))
        f.write(' {} \n'.format(l_max))
        f.write(' {} \n'.format(nke_max))
        f.write(' {} \n'.format(cusp))

    for l in range(l_max+1):
        f_name = "radfunc_l" + str(l) + ".dat"
        with open(f_name, 'w') as f:
            for j in range(nke_max):
                for i in range(Rnpt_fit):
                    f.write(' {:e} \n'.format(rad_func[i,l,j]))

    print('optimise linear coefficients')
    os.system("./fort/cLLSM " + " < " + inp_cLLSM)

    continuum_dir = "./output/continuum"
    os.system(f"mv ./{inp_cLLSM} {continuum_dir}")

    for l in range(l_max+1):

        dir_name = continuum_dir + "/l" + str(l) + "/"
        os.system(f"mkdir -p {dir_name}")
        f0 = "./coeff_l" + str(l) + ".dat"
        f1 = "./F_l"     + str(l) + ".dat"
        f2 = "./Fit_l"   + str(l) + ".dat"
        f3 = "./Err_l"   + str(l) + ".dat"
        dointerm = "mv " + f0 + " " + f1 + " " + f2 + " " + f3 + " " + dir_name
        os.system(dointerm)

        f_name = "./radfunc_l" + str(l) + ".dat"
        dointerm = "rm " + f_name
        os.system(dointerm)

    os.system("rm ./my_grid.dat")

# ---

def read_fits(l_max, nke_max):

    dir0 = "./output/continuum/l"
    for l in range(l_max+1):
        f_name = dir0 + str(l) + "/F_l" + str(l) + ".dat"
        with open(f_name, 'r') as file:
            lines = file.readlines()
            _n_pt = len(lines)

    fit_grd = np.zeros((_n_pt))
    fit_ref = np.zeros((_n_pt, l_max+1, nke_max), dtype=complex)
    fit_res = np.zeros((_n_pt, l_max+1, nke_max), dtype=complex)

    for l in range(l_max+1):
        f_name = dir0 + str(l) + "/F_l" + str(l) + ".dat"
        with open(f_name, 'r') as file:
            lines = file.readlines()
            for i,line in enumerate(lines):
                _line = line.split()
                fit_grd[i] = _line[0]

    for l in range(l_max+1):

        f_name = dir0 + str(l) + "/F_l" + str(l) + ".dat"
        with open(f_name, 'r') as file:
            lines = file.readlines()
            for i,line in enumerate(lines):
                _line = line.split()
                for _nke in range(nke_max):
                    fit_ref[i,l,_nke] = float(_line[2 * _nke + 1]) + 1j * float(_line[2 * _nke + 2])

        f_name = dir0 + str(l) + "/Fit_l" + str(l) + ".dat"
        with open(f_name, 'r') as file:
            lines = file.readlines()
            for i,line in enumerate(lines):
                _line = line.split()
                for _nke in range(nke_max):
                    fit_res[i,l,_nke] = float(_line[2 * _nke + 1]) + 1j * float(_line[2 * _nke + 2])

    return fit_grd, fit_ref, fit_res

# ---

