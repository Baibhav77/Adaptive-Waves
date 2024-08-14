import numpy as np
from nload import nload

def set_stack(mat, thick, lam, angle, *varargin):
    """
    Creates dictionary for the multilayered stack.

    Parameters:
    mat (list of str): Names for files containing n&k for materials.
    thick (list of float): Thicknesses of layers, notice first and last layer are infinite.
    lam (np.ndarray): Wavelengths.
    angle (float): Angle of incidence.
    varargin (optional): Additional arguments, e.g., threshold thickness for incoherence in layers.

    Returns:
    dict: Dictionary of stack ready for further analysis.
    """

    if lam.ndim == 1:
        lam = lam[:, np.newaxis]

    # Interpolating N&K values
    mat_unique, indices = np.unique(mat, return_inverse=True)
    n_unique = np.zeros((len(lam), len(mat_unique)), dtype=complex)

    for imat, material in enumerate(mat_unique):
        n_unique[:, imat] = nload(material, lam.flatten())

    n = n_unique[:, indices]

    # Creating stack dictionary
    if len(varargin) == 1:
        nincoh = np.array(thick) < varargin[0]
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

# Example usage (you should replace the following with your actual data)
if __name__ == "__main__":
    mat = ["SiO2", "TiO2", "SiO2", "TiO2"]
    thick = [112.25, 23.9375, 27, 43.875]
    lam = np.linspace(380, 1080, 851)
    angle = 0
    stack = set_stack(mat, thick, lam, angle, 50)
    print(stack)
