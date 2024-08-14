import numpy as np

def ATR1D(stack):
    """
    Calculates absorptance, transmittance and reflectance for 1D multilayered stack.

    Parameters:
    stack (dict): Stack dictionary containing nk, thick, nincoh, wavelength, and angle.

    Returns:
    tuple: Absorptance, transmittance, reflectance.
    """
    lam = stack["wavelength"]
    d = np.array(stack["thick"])
    theta0 = stack["angle"]
    n = stack["nk"]

    num_lam = len(lam)
    num_layers = len(d) + 1

    rs = np.zeros((num_lam, num_layers))
    rp = np.zeros((num_lam, num_layers))
    ts = np.zeros((num_lam, num_layers))
    tp = np.zeros((num_lam, num_layers))
    t_s = [[None] * num_layers for _ in range(num_lam)]
    t_p = [[None] * num_layers for _ in range(num_lam)]
    ph = [[None] * len(d) for _ in range(num_lam)]

    R_s = np.zeros(num_lam)
    R_p = np.zeros(num_lam)
    T_s = np.zeros(num_lam)
    T_p = np.zeros(num_lam)
    I_s = [None] * num_lam
    I_p = [None] * num_lam

    theta = np.zeros_like(n)
    theta[:, 0] = theta0

    for i in range(1, n.shape[1]):
        theta[:, i] = np.arcsin(n[:, i - 1] / n[:, i] * np.sin(theta[:, i - 1]))

    # Ensure d is broadcast correctly
    d_expanded = d[np.newaxis, :]
    cos_theta = np.cos(theta[:, 1:-1])
    d_cos_theta = 1
    lam_expanded = lam[:, np.newaxis]
    delt = 2 * np.pi * n[:, 1:-1] * d_cos_theta / lam_expanded
    edp = np.exp(1j * delt)
    edm = 1 / edp

    for i in range(num_layers):
        rs[:, i] = (n[:, i] * np.cos(theta[:, i]) - n[:, i + 1] * np.cos(theta[:, i + 1])) / (
                    n[:, i] * np.cos(theta[:, i]) + n[:, i + 1] * np.cos(theta[:, i + 1]))
        rp[:, i] = (n[:, i] * np.cos(theta[:, i + 1]) - n[:, i + 1] * np.cos(theta[:, i])) / (
                    n[:, i] * np.cos(theta[:, i + 1]) + n[:, i + 1] * np.cos(theta[:, i]))

        ts[:, i] = 2 * n[:, i] * np.cos(theta[:, i]) / (
                    n[:, i] * np.cos(theta[:, i]) + n[:, i + 1] * np.cos(theta[:, i + 1]))
        tp[:, i] = 2 * n[:, i] * np.cos(theta[:, i]) / (
                    n[:, i] * np.cos(theta[:, i + 1]) + n[:, i + 1] * np.cos(theta[:, i]))

    for i in range(num_lam):
        t_s[i][0] = np.array([[1, rs[i, 0]], [rs[i, 0], 1]]) / ts[i, 0]
        t_p[i][0] = np.array([[1, rp[i, 0]], [rp[i, 0], 1]]) / tp[i, 0]
        for j in range(1, num_layers):
            ph[i][j - 1] = np.array([[edm[i, j - 1], 0], [0, edp[i, j - 1]]])
            t_s[i][j] = np.array([[1, rs[i, j]], [rs[i, j], 1]]) / ts[i, j]
            t_p[i][j] = np.array([[1, rp[i, j]], [rp[i, j], 1]]) / tp[i, j]

    cfTs = np.real(n[:, -1] * np.cos(theta[:, -1]) / (n[:, 0] * np.cos(theta0)))
    cfTp = np.real(np.conj(n[:, -1]) * np.cos(theta[:, -1]) / (np.conj(n[:, 0]) * np.cos(theta0)))

    nincoh = stack["nincoh"]
    for i in range(num_lam):
        I_s[i] = np.eye(2)
        I_p[i] = np.eye(2)
        j = 0
        while j < len(d):
            tmp_s = t_s[i][j]
            tmp_p = t_p[i][j]
            while nincoh[j]:
                tmp_s = tmp_s @ ph[i][j] @ t_s[i][j + 1]
                tmp_p = tmp_p @ ph[i][j] @ t_p[i][j + 1]
                j += 1
                if j >= len(d):
                    I_s[i] = I_s[i] @ np.abs(tmp_s) ** 2
                    I_p[i] = I_p[i] @ np.abs(tmp_p) ** 2
                    break
            if j < len(d) and not nincoh[j]:
                I_s[i] = I_s[i] @ np.abs(tmp_s) ** 2 @ np.abs(ph[i][j]) ** 2
                I_p[i] = I_p[i] @ np.abs(tmp_p) ** 2 @ np.abs(ph[i][j]) ** 2
                j += 1
                if j >= len(d):
                    I_s[i] = I_s[i] @ np.abs(t_s[i][j]) ** 2
                    I_p[i] = I_p[i] @ np.abs(t_p[i][j]) ** 2
                    break
                elif not nincoh[j]:
                    I_s[i] = I_s[i] @ np.abs(t_s[i][j]) ** 2 @ np.abs(ph[i][j]) ** 2
                    I_p[i] = I_p[i] @ np.abs(t_p[i][j]) ** 2 @ np.abs(ph[i][j]) ** 2
                    j += 1
                    if j >= len(d):
                        I_s[i] = I_s[i] @ np.abs(t_s[i][j]) ** 2
                        I_p[i] = I_p[i] @ np.abs(t_p[i][j]) ** 2
                        break

        R_s[i] = I_s[i][1, 0] / I_s[i][0, 0]
        R_p[i] = I_p[i][1, 0] / I_p[i][0, 0]
        T_s[i] = cfTs[i] / I_s[i][0, 0]
        T_p[i] = cfTp[i] / I_p[i][0, 0]

    T_sp = (T_s + T_p) / 2
    R_sp = (R_s + R_p) / 2
    A_s = 1 - T_s - R_s
    A_p = 1 - T_p - R_p
    A_sp = 1 - T_sp - R_sp

    return {"A": A_sp, "T": T_sp, "R": R_sp}

# Example usage (you should replace the following with your actual data)
