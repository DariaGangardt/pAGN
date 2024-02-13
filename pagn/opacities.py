"""
opacities
Date: 05/06/23

Authors: Alessandro Trani, Daria Gangardt, Clément Bonnerot

Script containing methods that return an opacity value for given values of density and temperature, to be then used by
the AGN models in Sirko.py and Thompson.py.

"""
import numpy as np
from scipy.interpolate import RectBivariateSpline
import os

def electron_scattering_opacity(X=0.7):
    """
    Electron scattering opacity that depends on the hydrogen abundance fraction (see 17.1 of Kippenhahn)
    Args:
        X: hydrogen abundance fraction

    Returns: opacity in SI units
    """
    k = 0.20 * (1+X) * 0.1  # m^2 / kg
    return k

def get_combined_opacity():
    """
    Returns a function that interpolates the Semenov 2003 and Badnel 2005 opacities.
    Note: the function doest not use rho and T, as arguments, but T and R, where
    R = rho/T6^3, and T6 = 1e-6 T.
    The function takes *base 10 logarithm* of T and R, and returns *base 10 logarithm* of kappa
    units used: SI

    Returns
    -------
    kappa: obj
        Interpolated spline function of opacity
    """
    this_dir, this_filename = os.path.split(__file__)

    RTCmat = np.loadtxt(os.path.join(this_dir, "opacity_tables/logRT_combined_table.dat"))
    RTClogT = np.loadtxt(os.path.join(this_dir, "opacity_tables/logRT_combined_Taxis.dat"))
    RTClogR = np.loadtxt(os.path.join(this_dir, "opacity_tables/logRT_combined_Raxis.dat"))


    # Convert Kappa from cm^2 / g to m^2/kg
    RTCmat = RTCmat - 1
    # Convert density from g/cm^3 / to kg/m^3
    RTClogR = RTClogR + 3
    kappa = RectBivariateSpline(RTClogT, RTClogR, RTCmat, kx=2, ky=2, s=0.0) #need scipy ver 1.10.0
    return kappa

def get_semenov2003_opacity():
    """
    Returns a function that interpolates the Semenov 2003 opacities.
    Note: The function takes *base 10 logarithm* of T and R, and returns *base 10 logarithm* of kappa
    Units used: SI

    Returns
    -------
    kappa_bar: obj
        Interpolated spline function of opacity
    """
    this_dir, this_filename = os.path.split(__file__)
    table = np.loadtxt(os.path.join(this_dir, "opacity_tables/opacity_table.dat"))
    table = table*1e-1 #opacities given in cgs units
    temp = np.logspace(1, np.log10(9999), 1001)
    rho = np.logspace(-18, -7, 12)
    kappa = RectBivariateSpline(rho, temp, table.T, kx=2, ky=2, s=0.0)
    def kappa_bar(logT, logR, grid = False):
        T = 10**logT
        rho_SI = (10**logR)*((1e-6)*T)**3
        rho_cgs = rho_SI*1e-3
        kap = kappa(rho_cgs, T, grid = grid) #opacity grid obtained using cgs units
        kappa_b = np.log10(kap)
        return kappa_b
    return kappa_bar

def get_custom_opacity(table, rho, temp):
    """
    Returns a function that interpolates the given opacity table.
    Note: The returne function takes *base 10 logarithm* of T and R, and returns *base 10 logarithm* of kappa
    Units used: SI

    Parameters
    ----------
    table: 2D array
        Table of opacity values in SI units, given over a density x temperature grid.
    rho: 1D array
        Array of density values in SI units.
    temp: 1D array
        Array of temperature values in SI units.

    Returns
    -------
    kappa_bar: obj
        Interpolated spline function of opacity.

    """

    table = table.T
    kappa = RectBivariateSpline(rho, temp, table.T, kx=2, ky=2, s=0.0)
    def kappa_bar(logT, logR, grid = False):
        T = 10**logT
        rho_SI = (10**logR)*((1e-6)*T)**3
        kap = kappa(rho_SI, T, grid = grid)
        kappa_b = np.log10(kap)
        return kappa_b
    return kappa_bar


