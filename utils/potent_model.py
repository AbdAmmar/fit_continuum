
import numpy as np



# ---

def generate_coulomb(dr1, VR_min, VR_max, Z):

    Vnpt = int((VR_max - VR_min) / dr1)

    rxV_grid = np.array([VR_min + i * dr1 for i in range(Vnpt)])

    Vtmp = Z * np.ones(Vnpt, dtype=complex)
    Vtmp[0] = 0.0
    
    return(rxV_grid, Vtmp)

# ---

