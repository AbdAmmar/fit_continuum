
import numpy as np
from utils.CST import eV_to_au


# ---

# set 1
l_max = 5
# set 2
#l_max = 4

# ---

# nb of Gaussians
ng = 30
# if cusp = 1, add r^{l+1} in the Gaussian representation
cusp = 1
# range of fit
Rmax_fit = 25.0

# ---

# Coul or Dist
pot_type = "Coul"
#pot_type = "Dist"

# ---

# range of radial fcts
u_Rmin, u_Rmax, u_dr = 0., 50., 0.01
u_npt = int((u_Rmax - u_Rmin) / u_dr)
u_grid = np.array([u_Rmin + i*u_dr for i in range(u_npt)])

# ---

# energy grid
nke_max = 3

Ee_grid = [
    Ee * eV_to_au          # in a.u.
    for Ee in [12, 37, 74] # in eV
]

ke_grid = [(2.*Ee)**0.5  for Ee in Ee_grid] # in a.u.

# ---



