import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import trapz

def get_electric(lam, T):
    """
    FUNCTION ATR1D 
    calculates absorptance, transmittance and reflectance for 1D multilayered stack
    INPUT
    T - Transmittance
    OUTPUT
    IV - Dictionary containing jsc, Voc, V, I, P, and Pmax
    """
    # Constants
    nm = 1e-9
    ec = 1.60217663e-19
    h = 6.62607015e-34
    c = 3e8
    k = 1.380649e-23
    T_room = 300
    n = 1
    I0 = 55e-12
    cjsc = ec / (h * c)
    cvoc = n * k * T_room / ec

    # Ensure lam is a column vector
    lam = np.array(lam).flatten()

    # Load data
    AM15 = np.loadtxt('AM15.txt', dtype={'names': ('wavelength', 'direct'), 'formats': (np.float, np.float)})
    IQE = np.loadtxt('IQE.txt')

    AM15_interp = interp1d(AM15['wavelength'], AM15['direct'], kind='linear', fill_value='extrapolate')
    IQE_interp = interp1d(IQE[:, 0], IQE[:, 1], kind='linear', fill_value='extrapolate')

    AM15_lam = AM15_interp(lam)
    IQE_lam = IQE_interp(lam)

    # Calculations
    jsc = trapz(lam * nm, cjsc * AM15_lam * IQE_lam * T * lam / 10)
    Voc = cvoc * np.log(jsc / I0 + 1)
    V = np.linspace(0, Voc, 1000)
    I = jsc - I0 * (np.exp(V / cvoc) - 1)
    P = I * V
    Pmax = np.max(P)

    # Output dictionary
    IV = {
        'jsc': jsc,
        'Voc': Voc,
        'V': V,
        'I': I,
        'P': P,
        'Pmax': Pmax
    }

    return IV


