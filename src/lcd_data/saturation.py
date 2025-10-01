'''
This module provides functionality to calculate quantities related to water vapor pressure.
'''


def rh(temperature: float, temperature_dew_point: float) -> float:
    '''
    Compute relative humidity (RH) from air temperature and dew point temperature.

    This function calculates RH as the ratio of the saturation vapor pressure
    at the dew point temperature to that at the air temperature, expressed as a
    percentage.

    Parameters
    ----------
    temperature : float
        Air temperature in degrees Celsius.
    temperature_dew_point : float
        Dew point temperature in degrees Celsius.

    Returns
    -------
    float
        Relative humidity in percent (%).

    Notes
    -----
    - RH is computed as:
          RH = 100 × esatw(Td) / esatw(T)
      where `esatw` is the saturation vapor pressure (hPa).
    - Input temperatures should be in degrees Celsius.
    - Results outside the validity range of `esatw` (approximately –85 °C to +70 °C)
      may be inaccurate.

    '''

    return 100 * esatw(temperature_dew_point) / esatw(temperature)


def esatw(temperature: float) -> float:
    '''
    Compute the saturation water vapor pressure over a flat liquid water surface.

    This function uses an empirical 8th-order polynomial approximation (based on
    Flatau et al. 1992) to estimate the equilibrium vapor pressure of water above
    liquid water. The input is temperature in degrees Celsius, and the output is
    saturation vapor pressure in hectopascals (hPa).

    Parameters
    ----------
    temperature : float
        Air temperature in degrees Celsius.

    Returns
    -------
    float
        Saturation water vapor pressure (over liquid water) in hectopascals (hPa).

    Notes
    -----
    - The polynomial coefficients are derived from fits valid approximately over
      the range –85 °C to +70 °C (for water vapor over liquid) (Flatau et al. 1992).
    - Outside that temperature range the approximation error may grow significantly.
    - The formula assumes a flat liquid water surface (i.e., no curvature or
      capillarity effects) and standard atmospheric conditions.

    References
    ----------
    Flatau, P. J., Walko, R. L., & Cotton, W. R. (1992). Polynomial fits to saturation
    vapor pressure. *Journal of Applied Meteorology.*
    (See also documentation in NASA POWER methodology: valid –85 °C to +70 °C)
    '''

    a0, a1, a2, a3, a4, a5, a6, a7, a8 = (
        6.11239921,
        0.443987641,
        0.142986287e-1,
        0.264847430e-3,
        0.302950461e-5,
        0.206739458e-7,
        0.640689451e-10,
        -0.952447341e-13,
        -0.976195544e-15,
    )

    # Temperature relative to the triple point of water
    dt = temperature + 273.15 - 273.16

    satpw = a0 + dt * (a1 + dt * (a2 + dt * (a3 + dt * (a4 + dt * (a5 + dt * (a6 + dt * (a7 + a8 * dt)))))))

    return satpw
