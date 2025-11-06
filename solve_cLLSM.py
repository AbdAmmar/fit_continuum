
import os
import numpy as np

from utils import continuum_dir 


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

    os.system("ifort -mkl -O2 cLLSM.f90 -o cLLSM")
    print('optimise linear coefficients')
    dointerm = "./cLLSM " + " < " + inp_cLLSM
    os.system(dointerm)

    dointerm = "mv " + inp_cLLSM + " " + continuum_dir
    os.system(dointerm)

    for l in range(l_max+1):

        dir_name = continuum_dir + "/l" + str(l) + "/"
        f0 = "coeff_l" + str(l) + ".dat"
        f1 = "F_l"     + str(l) + ".dat"
        f2 = "Fit_l"   + str(l) + ".dat"
        f3 = "Err_l"   + str(l) + ".dat"
        dointerm = "mv " + f0 + " " + f1 + " " + f2 + " " + f3 + " " + dir_name
        os.system(dointerm)

        f_name = "radfunc_l" + str(l) + ".dat"
        dointerm = "rm " + f_name
        os.system(dointerm)

    os.system("rm cLLSM my_grid.dat")
# -------------------------------------------------------------------------------------------------

# ---

