import numpy as np
from scipy.interpolate import interp1d

def cielab(lam, spectrum):
    """
    cielab function calculates L*a*b parameters for a given spectrum
    INPUT
    lam - wavelength
    spectrum - spectrum of interest

    OUTPUT
    - L* = coordinate of brightness, visible gamut range [0, 100] (L = 100 indicates diffuse white; specular white might be higher);
    - a* = coordinate of chromaticity for red-green, gamut range [-128, 127];
    - b* = coordinate of chromaticity for yellow-blue, gamut range [-128, 127];

    REFERENCES
    - http://www.brucelindbloom.com/index.html?Equations.html
    - https://scipython.com/blog/converting-a-spectrum-to-a-colour/
    """

    # Load data
    colMatching = np.loadtxt('colMatching.txt')
    D65 = np.loadtxt('D65.txt')
    Veye = np.loadtxt('Veye.txt')

    # Rearrange input
    mask = (lam >= 380) & (lam <= 780)
    lam = lam[mask]
    spectrum = spectrum[mask]

    lamtab = np.arange(380, 781)
    colMatching = interp1d(lamtab, colMatching, axis=0)(lam)
    D65 = interp1d(lamtab, D65)(lam)

    # Calculate L*a*b
    XYZn = np.array([0.94811, 1, 1.07304])  # 10 deg observer

    Par2 = D65 * colMatching[:, 1]
    ContrX = (spectrum * D65) * colMatching[:, 0]
    ContrY = (spectrum * D65) * colMatching[:, 1]
    ContrZ = (spectrum * D65) * colMatching[:, 2]

    X = np.sum(ContrX) / np.sum(Par2)
    Y = np.sum(ContrY) / np.sum(Par2)
    Z = np.sum(ContrZ) / np.sum(Par2)

    Ref_X = X / XYZn[0]
    Ref_Y = Y / XYZn[1]
    Ref_Z = Z / XYZn[2]

    if Ref_X > 0.008856:
        f_X = Ref_X**(1/3)
    else:
        f_X = 7.787 * Ref_X + 16 / 116

    if Ref_Y > 0.008856:
        f_Y = Ref_Y**(1/3)
        L = 116 * f_Y - 16
    else:
        L = 903.3 * Ref_Y
        f_Y = 7.787 * Ref_Y + 16 / 116

    if Ref_Z > 0.008856:
        f_Z = Ref_Z**(1/3)
    else:
        f_Z = 7.787 * Ref_Z + 16 / 116

    a = 500 * (f_X - f_Y)
    b = 200 * (f_Y - f_Z)

    return L, a, b
