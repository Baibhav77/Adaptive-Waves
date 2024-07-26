import numpy as np
from scipy.interpolate import interp1d
import os

def nload(file, lam):
    """
    FUNCTION NLOAD
    interpolates n&k data from file to a specific wavelength range
    INPUT:
    file - 1. name of file with n&k data, three columns: wavelength (nm) | n | k
           2. constant real number if starts with "const="
           3. n=1 if set to "air"
    lam - wavelength (nm) for interpolation

    OUTPUT:
    ref - interpolated complex refractive index: ref = n + ik
    """
    if file == "air":
        ref = np.ones(len(lam), dtype=complex)
    elif file.startswith("const="):
        if "j" not in file:
            ref = float(file.split("=")[1]) * np.ones(len(lam), dtype=complex)
        else:
            realn = float(file.split("=")[1].split("+")[0])
            imagn = float(file.split("+")[1].replace("j", ""))
            ref = complex(realn, imagn) * np.ones(len(lam), dtype=complex)
    else:
        # Ensure the file has a .csv extension
        if not file.endswith(".csv"):
            file += ".csv"
        # Construct the full path to the CSV file
        file_path = os.path.join(os.path.dirname(__file__), file)
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File '{file}' not found at {file_path}.")
        
        data = np.loadtxt(file_path, delimiter=',')
        wavelength = data[:, 0]
        n_values = data[:, 1]
        k_values = data[:, 2]
        
        interp_n = interp1d(wavelength, n_values, kind='cubic', fill_value='extrapolate')
        interp_k = interp1d(wavelength, k_values, kind='cubic', fill_value='extrapolate')
        
        n_interp = interp_n(lam)
        k_interp = interp_k(lam)
        
        ref = n_interp + 1j * k_interp

    return ref
