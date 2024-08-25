import numpy as np
import matplotlib.pyplot as plt
from set_stack import set_stack
from ATR1D import ATR1D
from get_electric import get_electric

# GENERAL INPUT DATA
lam = np.linspace(380, 1080, 851)  # Wavelength range
dgls = 3200000  # Thickness of glass layer
dAZO = [5]  # Thickness of AZO layer
dEVA = [2000000]  # Thickness of EVA layer
dARC = [61.8]  # Thickness of ARC layer

Bragg1 = ["MgF2"]  # First Bragg stack material
Bragg2 = []  # Second Bragg stack material
matPSC = Bragg1 + ["GLS_NEW"] + Bragg2  # Combining Bragg1 and Bragg2 with glass in between

dBragg1 = [110.18]  # Thicknesses for Bragg1 (SiO2/TiO2 x2)
dBragg2 = []  # Thicknesses for Bragg2 (TiO2/SiO2 x2)
dPSC = dBragg1 + [dgls] + dBragg2
incoh = 1e3  # Incoherent layer is 1000nm thickness
theta = 0  # Angle of incidence

# PSC with AZO and EVA on PV

#This is for Full Stack

addLayer1 = "SiO2"
addLayer2 = "EVA"
addLayer3 = "const=2.4254"
addLayer4 = "mSi"

materials = ["air"] + matPSC + [addLayer1, addLayer2, addLayer3, addLayer4]
materials = ["air"] + matPSC + ["air", "air", "air", "air"]
d = dPSC + dAZO + dEVA + dARC  # Combine all thicknesses

# Creating the stack for the structure
stack = set_stack(materials, d, lam, theta, incoh)

# Calculating reflectance, transmittance, and absorptance
R_PSC_AZO_EVA_PV, T_PSC_AZO_EVA_PV, _ = ATR1D(stack)
A_PSC_AZO_EVA_PV = 1 - T_PSC_AZO_EVA_PV["sp"] - R_PSC_AZO_EVA_PV["sp"]  # Absorptance
# PLOTTING TRANSMITTANCE (T)
plt.figure()
plt.plot(lam, T_PSC_AZO_EVA_PV["sp"], 'b', linewidth=1.5, label='Transmittance (T)')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Transmittance')
plt.legend()
plt.title('Transmittance Plot')
plt.autoscale()
plt.show()

# PLOTTING REFLECTANCE (R)
plt.figure()
plt.plot(lam, R_PSC_AZO_EVA_PV["sp"], 'r', linewidth=1.5, label='Reflectance (R)')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Reflectance')
plt.legend()
plt.title('Reflectance Plot')
plt.autoscale()
plt.show()

# PLOTTING ABSORPTANCE (A)
plt.figure()
plt.plot(lam, A_PSC_AZO_EVA_PV, 'g', linewidth=1.5, label='Absorptance (A)')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Absorptance')
plt.legend()
plt.title('Absorptance Plot')
plt.autoscale()
plt.show()
