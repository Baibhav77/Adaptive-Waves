import numpy as np
from nload import nload

def set_stack(mat, thick, lam, angle, *varargin):
    """
    FUNCTION set_stack
    creates dictionary for the multilayered stack

    INPUT:
    mat - list of strings with names for files containing n&k for materials
    thick - list of thicknesses of layers, notice first and last layer are infinite
    lam - list of wavelengths
    angle - list of angles of incidence
    varargin - optional arguments:
               1. incoh - threshold thickness for incoherence in layers

    OUTPUT:
    stack - dictionary of stack ready for further analysis
    """
    if isinstance(lam, np.ndarray) and lam.ndim == 1:
        lam = lam[:, np.newaxis]

    # Interpolating N&K values
    mat_unique, ic = np.unique(mat, return_inverse=True)
    n_unique = np.zeros((len(lam), len(mat_unique)), dtype=complex)

    for i, material in enumerate(mat_unique):
        if material == "air":
            n_unique[:, i] = np.ones(len(lam), dtype=complex)  # Handle air directly
        else:
            n_unique[:, i] = nload(material, lam.flatten())

    n = n_unique[:, ic]

    # Creating stack dictionary
    if len(varargin) == 1:
        nincoh = thick < varargin[0]
    else:
        nincoh = np.ones_like(thick, dtype=bool)

    stack = {
        "nk": n,
        "thick": thick,
        "nincoh": nincoh,
        "wavelength": lam.flatten(),
        "angle": angle
    }

    return stack
