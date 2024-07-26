import os
import numpy as np
import matplotlib.pyplot as plt
from ATR1D import ATR1D
from set_stack import set_stack
from nload import nload

# Check if the file exists
file_path = "plain_1.csv"
if not os.path.exists(file_path):
    print(f"Error: File '{file_path}' not found.")
else:
    print(f"File '{file_path}' found.")

    # INPUT DATA
    # -------------------------------------------------------------------------
    # wavelengths of interest
    wavelengths = np.linspace(200, 1000, 801)
    theta = 0
    materials = ["air","plain_1","air"]  # Ensure this matches the file you are loading
    thicknesses = [100]  # Ensure this is correctly set up for your layers

    # Create stack
    stack = set_stack(materials, thicknesses, wavelengths, theta)

    # GETTING ABSORBANCE, TRANSMITTANCE, AND REFLECTANCE
    # -------------------------------------------------------------------------
    try:
        A, T, R = ATR1D(stack)
    except Exception as e:
        print(f"Error in ATR1D: {e}")
    else:
        # Ensure A, T, R are numpy arrays and plot them
        if isinstance(A, np.ndarray) and isinstance(T, np.ndarray) and isinstance(R, np.ndarray):
            plt.figure()
            plt.plot(wavelengths, A, label='Absorbance', linewidth=2)
            plt.plot(wavelengths, T, label='Transmittance', linewidth=2)
            plt.plot(wavelengths, R, label='Reflectance', linewidth=2)
            plt.legend()
            plt.xlim([200, 1000])
            plt.ylim([0, 1])
            plt.xlabel('Wavelength, nm', fontsize=16)
            plt.ylabel('Values', fontsize=16)
            plt.title('Absorbance, Transmittance, Reflectance')
            plt.grid(True)
            plt.show()
        else:
            print("Error: A, T, and R should be numpy arrays.")
