
from INPUT_OPT import u_Rmin, u_Rmax, u_dr
from INPUT_OPT import ng, Rmax_fit, cusp
from INPUT_OPT import l_max, nke_max, ke_grid, Ee_grid
from INPUT_OPT import pot_type

from potent_model import generate_coulomb
from solve_radial import radial_inputfiles, execute_radial, get_radial_fcts, man_rad, plot_rad_func
from solve_cLLSM import run_cLLSM

import numpy as np
import time


if __name__=="__main__":
    t_beg = time.time()

    print(" Local current time: {}".format(time.asctime(time.localtime(time.time()))))


    t0 = time.time()
    if(pot_type == "Coul"):

        print(" --- coulomb potential --- ")
        # radial potential ( r x V(r) )
        rxV_grid, rxV = generate_coulomb(dr1=0.01, VR_min=0.0, VR_max=20.0, Z=-1.0)

    else:
        print(f" {pot_type} not implemented yet")
        quit()
    

    # RADIAL input
    u_npt = int((u_Rmax - u_Rmin) / u_dr)
    u_grid = np.array([u_Rmin + i*u_dr for i in range(u_npt)])

    eps = 1e-15

    print(" --- call RADIAL --- ")
    radial_inputfiles(nke_max, Ee_grid, l_max, eps, rxV_grid, rxV, u_grid)
    execute_radial(rxV, l_max, FC="gfortran")
    rad_func, rad_shif_Inner, rad_shif_Coul = get_radial_fcts(rxV, u_npt, l_max, nke_max)
    man_rad(rxV, l_max, nke_max)
    plot_rad_func(u_grid, rad_func, l_max, nke_max)
    print(" ")

    # optimise linear coefficients
    # if cusp = 1, add r^{l+1} in the Gaussian representation
    run_cLLSM(u_grid, Rmax_fit, rad_func, ng, l_max, nke_max, cusp)

    print(" --- done work after {} minutes --- ".format((time.time()-t_beg)/60.))



