import numpy as np
import scipy.io
from scipy.interpolate import interp1d

def nload(file, lam):
    """
    Interpolates n&k data from file to a specific wavelength range.

    Parameters:
    file (str): Name of file with n&k data, constant real number if starts with "const=", or "air".
    lam (np.ndarray): Wavelength (nm) for interpolation.

    Returns:
    np.ndarray: Interpolated complex refractive index: ref = n + ik.
    """
    if file == "air":
        ref = np.ones(len(lam), dtype=complex)
    elif file.startswith("const="):
        if "j" not in file:
            ref = float(file.split("=")[1]) * np.ones(len(lam), dtype=complex)
        else:
            parts = file.split("=")[1].split("+")
            realn = float(parts[0])
            imagn = float(parts[1][:-1])
            ref = complex(realn, imagn) * np.ones(len(lam), dtype=complex)
    else:
        data = scipy.io.loadmat(file)
        print(f"Loaded data from {file}: {data}")

        # Access the required data directly
        material_key = file.split('/')[-1].split('.')[0]  # Extract key based on filename
        ntable = data[material_key]
        print(f"Data under key '{material_key}': {ntable}")

        if isinstance(ntable, np.ndarray):
            print(f"Regular ndarray data: {ntable}")
        else:
            raise TypeError("Loaded data is not an ndarray.")

        # Ensure the data is in float format
        ntable = ntable.astype(float)
        n_interp = interp1d(ntable[:, 0], ntable[:, 1], kind='cubic', fill_value="extrapolate")
        k_interp = interp1d(ntable[:, 0], ntable[:, 2], kind='cubic', fill_value="extrapolate")
        ref = n_interp(lam) + 1j * k_interp(lam)

    return ref

# Example usage (you should replace the following with your actual data)
if __name__ == "__main__":
    lam = np.linspace(380, 1080, 851)
    file = "SiO2.mat"  # Replace with your actual file path
    ref_index = nload(file, lam)
    print(ref_index)
