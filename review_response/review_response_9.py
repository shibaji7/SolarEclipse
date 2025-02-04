import numpy as np

def electron_density_to_plasma_frequency(n_e, output_in_hertz=True):
    """
    Converts electron density to plasma frequency.

    Parameters:
        n_e (float): Electron density in electrons per cubic meter (m^-3).
        output_in_hertz (bool, optional): If True, output plasma frequency in Hertz (Hz),
            otherwise in radians per second (rad/s). Defaults to True.

    Returns:
        float: Plasma frequency in Hz or rad/s
    """
    e = 1.602e-19  # Elementary charge (Coulombs)
    epsilon_0 = 8.854e-12  # Permittivity of free space (F/m)
    m_e = 9.109e-31  # Electron rest mass (kg)
    pi = np.pi

    if n_e < 0:
      raise ValueError(f"Electron density n_e must be non-negative: {n_e}")

    omega_p = np.sqrt((n_e * e**2) / (epsilon_0 * m_e))

    if output_in_hertz:
        return omega_p / (2 * pi)
    else:
        return omega_p

ed = 76.5 * 1e6
print(electron_density_to_plasma_frequency(ed)/1e3)