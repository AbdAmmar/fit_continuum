
import os
import numpy as np



def radial_inputfiles(nke_max, Ee_grid, l_max, eps, rxV_grid, rxV, Rgrid):

    with open('my_grid.dat', 'w') as f:
        for i in range(Rgrid.shape[0]):
            f.write(' {:e}\n'.format(Rgrid[i]))

    pot_file = 'potential.dat'

    with open(pot_file, 'w') as f:
        for i in range(rxV_grid.shape[0]):
            f.write(' {:e}  {:e}\n'.format(rxV_grid[i], rxV[i].real))

    for l in range(l_max+1):
        inp_file = 'numpoten_l' + str(l) + ".in"

        with open(inp_file, 'w') as f:
            f.write('{}\n'.format(pot_file))
            f.write('my_grid.dat\n')
            for i in range(nke_max):
                # option for RADIAL 
                f.write('2\n')
                # energy      momentum    tolerance for RADIAL
                Ee = Ee_grid[i]
                f.write('{:e} {:d} {:e}\n'.format(Ee, l, eps))
                # name of output script of function
                f_name = 'l' + str(l)
                f.write('{}_k{}.dat\n'.format(f_name,i+1))
            f.write('-1\n')

# ---

def execute_radial(rxV, l_max):

    for l in range(l_max+1):
        print(' executing RADIAL for l = {}'.format(l))

        inp_file = 'numpoten_l' + str(l) + ".in"
        out_file = 'numpoten_l' + str(l) + ".out"
        dointerm = "./fort/numpoten < " + inp_file + " > " + out_file
        os.system(dointerm)

    continuum_dir = "./output/continuum/"
    os.makedirs(os.path.dirname(continuum_dir), exist_ok=True)
    dointerm = "mv SPLERR.dat potential.dat resn.dat pot-spline.dat " + continuum_dir
    os.system(dointerm)

# ---

def get_radial_fcts(rxV, u_npt, l_max, nke_max):

    func       = np.zeros((u_npt, l_max+1, nke_max))
    shif_Inner = np.zeros((l_max+1, nke_max))
    shif_Coul  = np.zeros((l_max+1, nke_max))
    for l in range(l_max+1):

        for i in range(nke_max):
            f_name = 'l' + str(l) + '_k' + str(i+1) + ".dat"
            lines  = []
            f      = open(f_name, 'r')
            lines  = f.readlines()
            f.close()

            for j in range(u_npt):
                y = lines[8+j].split()[1]
                func[j,l,i] = float(y)

            if(len(lines[4].split()[3]) == 6):
                #sh_in = lines[4].split()[4][0:17]
                sh_in = lines[4].split()[4]
            else:
                #sh_in = lines[4].split()[3][6:24]
                sh_in = lines[4].split()[3][6:]

            if(len(lines[5].split()[3]) == 6):
                #sh_cb = lines[5].split()[4][0:17]
                sh_cb = lines[5].split()[4]
            else:
                #sh_cb = lines[5].split()[3][6:24]
                sh_cb = lines[5].split()[3][6:]

            shif_Inner[l,i] = float(sh_in)
            shif_Coul [l,i] = float(sh_cb)

        psh_file = 'PhaseShift_l' + str(l) + ".dat"
        with open(psh_file, 'w') as f:
            for i in range(nke_max):
                f.write('{}   {}\n'.format(shif_Inner[l,i], shif_Coul[l,i]))

    return(func, shif_Inner, shif_Coul)

# ---

def man_rad(rxV, l_max, nke_max):

    continuum_dir = "./output/continuum"
    os.makedirs(os.path.dirname(continuum_dir), exist_ok=True)

    for l in range(l_max+1):
        inp_file = './numpoten_l'   + str(l) + ".in"
        out_file = './numpoten_l'   + str(l) + ".out"
        psh_file = './PhaseShift_l' + str(l) + ".dat"
        dir_name = continuum_dir + "/l" + str(l) + "/"
        os.makedirs(os.path.dirname(dir_name), exist_ok=True)
        dointerm = "mv " + inp_file + " " + out_file + " " + psh_file + " " + dir_name
        os.system(dointerm)

        for i in range(nke_max):
            f_name = './l' + str(l) + '_k' + str(i+1) + ".dat"
            dointerm = "mv " + f_name + " " + dir_name
            os.system(dointerm)

# ---

