import numpy as np
from scipy.interpolate import interp1d

def heat(spectrum, emis):
    """
    Calculates heat-related characteristics of low-E glasses
    INPUT:
    spectrum - numpy array with columns [lambda, T, R] (wavelength, transmittance, reflectance)
    emis - emissivity, a portion of absorbance converted to heat

    OUTPUT:
    Tvis - visible light transmittance
    SHGC - solar heat gain coefficient
    LSG - light-to-solar gain ratio
    """
    # Load data
    AM15 = np.loadtxt('AM15.txt')  # solar illumination
    Veye = np.loadtxt('Veye.txt')  # human response eye function

    # Tvis calculation
    interp_T = interp1d(spectrum[:, 0], spectrum[:, 1], kind='cubic', fill_value='extrapolate')
    interp_AM15_vis = interp1d(AM15[:, 0], AM15[:, 1], kind='cubic', fill_value='extrapolate')
    Tvis_spectrum = interp_T(Veye[:, 0])
    AM15_vis = interp_AM15_vis(Veye[:, 0])
    Tvis = np.sum(Tvis_spectrum * Veye[:, 1] * AM15_vis) / np.sum(Veye[:, 1] * AM15_vis)

    # SHGC calculation
    lam_IR = np.linspace(300, 2000, 1700)
    AT = spectrum[:, 1] + emis * (1 - spectrum[:, 1] - spectrum[:, 2])
    interp_AM15_IR = interp1d(AM15[:, 0], AM15[:, 1], kind='cubic', fill_value='extrapolate')
    interp_AT_IR = interp1d(spectrum[:, 0], AT, kind='cubic', fill_value='extrapolate')
    AM15_IR = interp_AM15_IR(lam_IR)
    AT_IR = interp_AT_IR(lam_IR)
    SHGC = np.sum(AT_IR * AM15_IR) / np.sum(AM15_IR)

    # LSG calculation
    LSG = Tvis / SHGC

    return Tvis, SHGC, LSG

# Example usage
# spectrum = np.array(...)  # replace with your actual spectrum data
# emis = ...  # replace with your actual emissivity value
# Tvis, SHGC, LSG = heat(spectrum, emis)
