
import os
import matplotlib.pyplot as plt


# ---

def plot_func(fgrid, f, fig_output):
    plt.xlabel('r (a.u.)')
    plt.plot(fgrid, f)
    plt.savefig(fig_output)
    plt.close()

# ---

def plot_2func(fgrid, f_ref, f, fig_output):
    plt.xlabel('r (a.u.)')
    plt.plot(fgrid, f_ref, "-", label="ref")
    plt.plot(fgrid, f, "--")
    plt.legend()
    plt.savefig(fig_output)
    plt.close()

# ---


def plot_rad_func(Rgrid, rad_func, l_max, nke_max):

    continuum_dir = "./output/continuum"

    for l in range(l_max+1):

        dir_name = continuum_dir + "/l" + str(l) + "/"
        os.makedirs(os.path.dirname(dir_name), exist_ok=True)

        for i in range(nke_max):
            f_name     = 'l' + str(l) + '_k' + str(i+1) + ".pdf"
            fig_output = dir_name + f_name
            plot_func(Rgrid, rad_func[:,l,i].real, fig_output)

# ---

def plot_rad_func(Rgrid, rad_func, l_max, nke_max):

    continuum_dir = "./output/continuum"

    for l in range(l_max+1):

        dir_name = continuum_dir + "/l" + str(l) + "/"
        os.makedirs(os.path.dirname(dir_name), exist_ok=True)

        for i in range(nke_max):
            f_name     = 'l' + str(l) + '_k' + str(i+1) + ".pdf"
            fig_output = dir_name + f_name
            plot_func(Rgrid, rad_func[:,l,i].real, fig_output)

# ---

def plot_fit(fit_grid, fit_ref, fit_res, l_max, nke_max):

    continuum_dir = "./output/continuum"

    for l in range(l_max+1):

        dir_name = continuum_dir + "/l" + str(l) + "/"
        os.makedirs(os.path.dirname(dir_name), exist_ok=True)

        for i in range(nke_max):
            f_name     = 'Fit_l' + str(l) + '_k' + str(i+1) + ".pdf"
            fig_output = dir_name + f_name
            plot_2func(fit_grid, fit_ref[:,l,i].real, fit_res[:,l,i].real, fig_output)

# ---




