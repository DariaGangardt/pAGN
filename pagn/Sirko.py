"""
Sirko
Date: 01/06/23

Authors: Alessandro A. Trani, Daria Gangardt

Module that uses the equations from Sirko & Goodman 2003 to produce a self-gravitating Active Galactic Nuclei disc

Contains SirkoAGN class, which creates an AGN disc object, documented below.

"""

import matplotlib.pyplot as plt
import numpy as np
from . import constants as ct
from .opacities import electron_scattering_opacity, get_semenov2003_opacity, get_combined_opacity, get_custom_opacity
from scipy.optimize import root
import matplotlib.colors as mpl_col

class SirkoAGN:
    def __init__(self, Mbh=1e8*ct.MSun, le=0.5, alpha=0.01, b=0, Mdot=None,
                 eps=0.1, opacity = "combined", X=0.7, debug=False, xtol=1e-10,
                 rootm="lm"):
        """ Class that creates an AGN disc object using the equations from Sirko and Goodman 2003.

        Parameters
        ----------
        Mbh: float, optional (default: 1e8*ct.MSun)
            Mass of the Super Massive Black Hole (SMBH) fuelling the AGN in kg.
        le: float, optional (default: 0.5)
            Luminosity ratio L0/LE, where L0 is the non self gravitating luminosity and LE is the Eddington luminosity of the SMBH.
        alpha: float, optional (default: 0.01)
            Viscosity parameter of the inner region of the AGN.
        b: float, optional (default: 0.)
            Power index for viscosity-gas pressure relation, can only be 0 or 1.
        Mdot: float, optional (default: None)
            Mass accretion rate, taken as constant throughout the disc, in units of kg/s. Can be given instead of le.
        eps: float, optional (default: 0.1)
            SMBH radiative efficiency parameter.
        opacity: string/tuple, optional (default: 'combined')
            Which opacity table to use. Options are "semenov", which extrapolates values from Semenov et al. 2003 only,
            "combined", which combines the Semenov et al. 2003 values with the Badnell et al. 2005 values, or a custom
            table of values entered as a tuple made up of three arrays:
             (opacity 2D array, density 1D array, temperature 1D array). For the custom table, inputs must be in SI units.
        X: float, optional (default: 0.7)
            Hydrogen abundance of gas in AGN disk.
        debug: bool, optional (default: False)
            Enter debug regime. Code will check root finder operations by printing solutions at every value of r.
        xtol: float, optional (default: 1e-10)
            Tolerance in the root finding method.
        rootm: string, optional (default: 'lm')
            Root finding method to be used for both star and no star formation regimes.

        """
        # Disk parameters
        self.Mbh = Mbh
        self.alpha = alpha
        self.b = b

        # b can only be 0 or 1
        if b not in [0, 1]:
            raise ValueError(f"b can only be 0 or 1 (found {b})")

        # Auxiliary parameters
        self.eps = eps
        self.X = X
        self.debug = debug
        self.xtol = xtol
        self.rootm = rootm

        Ledd = self.L_Edd()
        if Mdot is None and le is None:
            raise ValueError("Please provide either an accretion rate or an Eddington ratio.")
        elif Mdot is None:
            Mdot_edd = Ledd / (self.eps * ct.c ** 2)
            self.le = le
            self.Mdot = self.le * Mdot_edd
        else:
            self.Mdot = Mdot
            self.le = self.Mdot * ct.c ** 2 * self.eps / Ledd

        self.Rs = 2 * self.Mbh * ct.G / (ct.c ** 2)
        self.Rmin = self.Rs / (4 * eps)

        self.Rmax = 1e7 * self.Rs

        # Opacity table variables
        if type(opacity) != str:
            # custom opacity values
            if len(opacity) == 3:
                self.opacity = "custom"
                (opac, rho, temp) = opacity
                if opac.shape == (len(rho), len(temp)):
                    self.kappaLogTLogRopac = get_custom_opacity(opac, rho, temp,)
                else:
                    raise ValueError("Please enter a tuple of 3 arrays: (a numpy 2D array of kappa values in SI units,"
                                     " a numpy 1D array of density values in SI units, a numpy 1D array of temperature values in SI units), "
                                     "making sure that the kappa values were obtained over a density x temperature grid.")
            else:
                raise ValueError("Please enter a tuple of 3 arrays: (a numpy 2D array of kappa values in SI units,"
                                 " a numpy 1D array of density values in SI units, a numpy 1D array of temperature values in SI units), "
                                 "making sure that the kappa values were obtained over a density x temperature grid.")
        elif type(opacity) == str:
            self.opacity = opacity.lower()
            if self.opacity == "combined":
                self.kappaLogTLogRopac = get_combined_opacity()
            elif self.opacity == "semenov":
                self.kappaLogTLogRopac = get_semenov2003_opacity()
        else:
            raise ValueError("Available opacities are combined, semenov or custom.")

        self.i = 0  # current radial index
        self.isf = -1  # radial index where star formation starts
        self.Mdisk = 0.0  # We will integrate this

        print("### Sirko & Goodman 2003 parameters ###")
        print("Mbh = {:e} MSun".format(self.Mbh / ct.MSun))
        print("Mdot = {:e} MSun/yr".format(self.Mdot * (ct.yr / ct.MSun)))
        print("le =", self.le)
        print("Rs = {:e} pc".format(self.Rs / ct.pc))
        print("Rmin = {:e} Rs".format(self.Rmin / self.Rs))
        print("Rmax = {:e} Rs, {:e} pc".format(self.Rmax / self.Rs, self.Rmax / ct.pc))
        print("alpha =", self.alpha)
        print("b =", self.b)
        print("eps =", self.eps)
        print("X =", self.X)
        print("Opacity =", self.opacity)

        print("\ndebug =", self.debug)
        print("xtol =", self.xtol)
        print("root method =", self.rootm)

    def solve_disk(self, N=1e4):
        """ Method to evolve the AGN disc, from outer boundary inwards, using Sirko and Goodman 2003 equations.

        Parameters
        ----------
        N: float, optional (default: 1e4)
            Disc radial resolution.

        """
        N = int(N)
        Rmin = 6 * self.Rs  # Cannot go too close to Rmin.
        self.Redge = np.logspace(np.log10(Rmin), np.log10(self.Rmax), N + 1)
        self.R = (self.Redge[:-1] + self.Redge[1:]) * 0.5
        self.deltaR = self.Redge[1:] - self.Redge[:-1]

        # Solutions
        self.Omega = np.sqrt(ct.G * self.Mbh / (self.R*self.R*self.R))
        # Only valid for no star formation
        Teff4 = 3 * self.Mdot * (1 - np.sqrt(self.Rmin / self.R)) * self.Omega * self.Omega / (8 * np.pi)
        self.Teff4 = Teff4 / ct.sigmaSB
        self.rho = np.zeros_like(self.R)
        self.T = np.zeros_like(self.R)
        self.h = np.zeros_like(self.R)
        self.tauV = np.zeros_like(self.R)
        self.kappa = np.zeros_like(self.R)
        self.cs = np.zeros_like(self.R)
        self.Q = np.zeros_like(self.R)

        root_options = {'col_deriv': True, 'factor': 0.1, 'xtol': self.xtol}
        if self.rootm == "hybr": root_options['maxfev'] = 10000
        x_guess = np.log10(np.array([self.Omega[self.i] * self.R[self.i] * 1e-2, self.Teff4[self.i] ** 0.25]))

        if self.debug: print("x_guess, log10(T) and log10(Teff) at innermost R:", x_guess)

        # Inside out, starting with no star formation
        for self.i in range(N):
            sol = root(self.no_starformation, x_guess, method=self.rootm,
                       options=root_options,)

            if sol.success:
                # T, cs, tauV, rho and kappa are assigned in the function itself
                # Calculate everything else
                self.Q[self.i] = self.Omega[self.i] ** 2 / (2 * np.pi * ct.G * self.rho[self.i])
                self.h[self.i] = self.cs[self.i] / self.Omega[self.i]
                self.Mdisk += 2 * np.pi * 2 * self.h[self.i] * self.deltaR[self.i] * self.R[self.i] * self.rho[self.i]

                # next guess
                x_guess = np.log10((self.cs[self.i], self.T[self.i]))

                if self.debug:
                    print("Looking for roots in pressure equilibrium, thermal equilibrium equations; no star formation")
                    sol_check_zeros = self.no_starformation(sol.x)
                    print("Zero residuals at i =", self.i, " :", sol_check_zeros)

                #Check if still in star formation regime
                if self.Q[self.i] < 1:
                    print("Q<1 at i={:d} (R={:1.2e} Rs)".format(self.i, self.R[self.i] / self.Rs))
                    self.isf = self.i
                    break
            else:
                sol_check_zeros = self.no_starformation(sol.x)
                sol_check_zeros = sol_check_zeros
                print("No no star formation regime solution found at index", self.i, "\nZero residuals:", sol_check_zeros)
                print("Last solution:", sol.x)
                print("Showing root finding solutions for no star formation regime")
                self.plot_roots_no_star(guess=sol.x, zoomguess=True)
                break

            print("{:d}/{:d}".format(self.i, N), end="\r")

        # Inside out, extending the disc with star formation
        if self.isf > 0:
            x_guess = np.log10((self.Teff4[self.isf] ** 0.25, self.T[self.isf]))
            print("Beginning star formation at index", self.isf)
            for self.i in range(self.isf, N):
                #Star formation maintains disc from vertical gravitational collapse, ensuring Q=1
                self.Q[self.i] = 1
                self.rho[self.i] = self.Omega[self.i] ** 2 / (2 * np.pi * ct.G)
                if not self.b:
                    #if b = 0, we can find cs and h before entering the root finder
                    Mdotprime = self.Mdot * (1 - (self.Rmin / self.R[self.i]) ** 0.5)
                    self.cs[self.i] = (Mdotprime * self.Omega[self.i] ** 2 / (
                            6 * np.pi * self.alpha * self.rho[self.i])) ** (1 / 3)
                    self.h[self.i] = self.cs[self.i] / self.Omega[self.i]

                sol = root(self.yes_starformation, x_guess, method=self.rootm,
                           options=root_options,)

                if sol.success:
                    self.Mdisk += 2 * np.pi * 2 * self.h[self.i] * self.deltaR[self.i] * self.R[self.i] * self.rho[self.i]
                    # next guess
                    x_guess = np.log10((self.Teff4[self.i]**0.25, self.T[self.i]))

                    if self.debug:
                        print("Looking for roots in pressure equilibrium, thermal equilibrium equations; yes star formation")
                        sol_check_zeros = self.yes_starformation(sol.x)
                        print("Zero residuals at i =", self.i, " :", sol_check_zeros)

                else:
                    sol_check_zeros = self.yes_starformation(sol.x)
                    sol_check_zeros = sol_check_zeros
                    print("No star formation regime solution found at index", self.i, "\nZero residuals:", sol_check_zeros)
                    print("Last solution:", sol.x)
                    print("Showing root finding solutions for star formation regime")
                    self.plot_roots_yes_star(guess=sol.x, zoomguess=True)
                    break

                print("{:d}/{:d}".format(self.i, N), end="\r")
        self.R_AGN = self.R[self.isf]
        print("Mdisk =", self.Mdisk / ct.MSun, "Msun")
        print("Mdisk/Mbh =", self.Mdisk / self.Mbh)

    def L_Edd(self,):
        """ Method to calculate the Eddington luminosity given a mass and hydrogen abundance

        Returns
        -------
        L: float
            Eddington luminosity of the Super Massive Black Hole (SMBH) fuelling the AGN in watts.
        """
        kes = electron_scattering_opacity(X=self.X)
        L = 4 * np.pi * ct.G * ct.c * self.Mbh / kes
        return L

    def plot(self, params="all"):
        """ Simple disk plotting method, can be used to check that disk has been evolved correctly.

        Parameters
        ----------
        params: str or list of str, optional (default: "all")
            List of parameters to plot, can be any parameter combination from ['h', 'rho', 'tau', 'T']

        """
        all_params = ['h', 'rho', 'tau', 'T']
        param_dic = {'h': [r"$\log_{10}{h/r}$", np.log10(self.h / self.R)],
                     'rho': [r"$\log_{10}{\rho \, [\mathrm{g cm^{-3}}]}$", np.log10(self.rho * ct.SI_to_gcm3)],
                     'tau': [r"$\log_{10}{\tau}$", np.log10(self.tauV)],
                     'T': [r"$\log_{10}{T \, [\mathrm{K}]}$", np.log10(self.T)], }
        if params == "all" or set(all_params).issubset(params):
            f, ax = plt.subplots(4, 1, figsize=(10, 15), sharex=True, gridspec_kw=dict(hspace=0), tight_layout=True)
            ax[0].plot(np.log10(self.R / self.Rs), param_dic['h'][1])
            ax[0].set_ylabel(param_dic['h'][0])
            ax[1].plot(np.log10(self.R / self.Rs), param_dic['rho'][1])
            ax[1].set_ylabel(param_dic['rho'][0])
            ax[2].plot(np.log10(self.R / self.Rs), param_dic['tau'][1])
            ax[2].set_ylabel(param_dic['tau'][0])
            ax[3].plot(np.log10(self.R / self.Rs), param_dic['T'][1])
            ax[3].set_ylabel(param_dic['T'][0])
            ax[3].set_xlabel(r"$\log_{10}{R \, [R_S]}$")
            for a in ax:
                a.axvline(np.log10(self.R_AGN / self.Rs), -100, 100)

        elif set(params).issubset(all_params):
            f, ax = plt.subplots(len(params), 1, figsize=(10, 10 + int(len(params)) * (2 / 3)), sharex=True,
                                 gridspec_kw=dict(hspace=0), tight_layout=True)
            ax[-1].set_xlabel(r"$\log_{10}{R \, [R_S]}$")
            for i_param, param in enumerate(params):
                ax[i_param].plot(np.log10(self.R / self.Rs), param_dic[param][1])
                ax[i_param].set_ylabel(param_dic[param][0])
                ax[i_param].axvline(np.log10(self.R_AGN / self.Rs), -100, 100)

        plt.show()

    def plot_roots_no_star(self, guess=None, zoomguess=False):
        """ Plots current solution space for root finder in the no star formation regime.

        Parameters
        ----------
        guess: array, optional (default: None)
            Array of log10(cs), log10(T) guesses to be plotted on solution space.
        zoomguess: bool, optional (default: False)
            Flag to zoom in on given guess in plot.
        """
        if guess is None: zoomguess = False

        print("### Plotting no star formation regime solutions ###")
        print("i = {:} ".format(self.i))
        print("R = {:e} Rs = {:e} pc".format(self.R[self.i] / self.Rs, self.R[self.i] / ct.pc, ))
        if guess is not None:
            print("guess = (log cs, log T) = ({:} , {:})".format(guess[0], guess[1]))
        else:
            print("guess = {}".format(guess))
        print("zoom guess = {}".format(zoomguess))
        print("###")

        if zoomguess:
            cs_arr = np.linspace(guess[0]-1, guess[0]+1, 100)
            T_arr = np.linspace(guess[1]-1, guess[1]+1, 100)
        else:
            cs_arr = np.linspace(3, 10, 200)
            T_arr = np.linspace(1.5, 9, 200)

        cs_axx, T_axx = np.meshgrid(cs_arr, T_arr)
        Zp = np.zeros(shape=(len(cs_arr), len(T_arr)))
        Zr = np.zeros(shape=(len(cs_arr), len(T_arr)))
        for i in range(len(cs_arr)):
            for j in range(len(T_arr)):
                f = self.no_starformation(np.array([cs_axx[i,j], T_axx[i,j]]), set=False)
                Zp[i,j] = f[0] #np.log10(f[0]) if f[0] > 0 else -np.log10(-f[0])
                Zr[i, j] = f[1] #np.log10(f[1]) if f[1] > 0 else -np.log10(-f[1])

        f, ax = plt.subplots(1, 3, figsize=(16, 7), sharey=True, tight_layout=True) #gridspec_kw=dict(hspace=0),
        pcol = ax[0].pcolormesh(cs_axx, T_axx, Zp, cmap="coolwarm", norm=mpl_col.CenteredNorm(halfrange=1))
        f.colorbar(pcol, ax=ax[0], label = r'$\mathrm{Balancing \, Pressures}$')
        ax[0].set_xlabel(r'$\log{c_s}$')
        ax[0].set_ylabel(r'$\log{T}$')

        pcol = ax[1].pcolormesh(cs_axx, T_axx, Zr, cmap="coolwarm", norm=mpl_col.CenteredNorm(halfrange=1))
        f.colorbar(pcol, ax=ax[1], label = r'$\mathrm{Balancing \, Temperatures}$')
        ax[1].set_xlabel(r'$\log{c_s}$')
        ax[1].set_ylabel(r'$\log{T}$')

        CS = ax[2].contour(cs_axx, T_axx, Zr, [0.0])
        CS = ax[2].contour(cs_axx, T_axx, Zp, [0.0], linestyles = '--')
        ax[2].set_xlabel(r'$\log{c_s}$')
        ax[2].set_ylabel(r'$\log{T}$')

        if guess is not None:
            for axx in ax:
                axx.axvline(guess[0], c="green", alpha=0.5)
                axx.axhline(guess[1], c="green", alpha=0.5)

        plt.show()

    def plot_roots_yes_star(self, guess=None, zoomguess=False):
        """ Plots current solution space for root finder in the star formation regime.

        Parameters
        ----------
        guess: array, optional (default: None)
           Array of log10(T), log10(eta) guesses to be plotted on solution space.
        zoomguess: bool, optional (default: False)
           Flag to zoom in on given guess in plot.
        """
        print("### Plotting yes star formation regime solutions ###")
        print("i = {:} ".format(self.i))
        print("R = {:e} Rs = {:e} pc".format(self.R[self.i] / self.Rs, self.R[self.i] / ct.pc, ))
        if guess is not None:
            print("guess = (log Teff, log T) = ({:} , {:})".format(guess[0], guess[1]))
        else:
            print("guess = {}".format(guess))
        print("zoom guess = {}".format(zoomguess))
        print("###")
        if zoomguess:
            Teff_arr = np.linspace(guess[0]-1, guess[0]+1, 100)
            T_arr = np.linspace(guess[1]-1, guess[1]+1, 100)
        else:
            Teff_arr = np.linspace(1, 5, 200)
            T_arr = np.linspace(1.5, 6, 200)

        Teff_axx, T_axx = np.meshgrid(Teff_arr, T_arr)
        Zp = np.zeros(shape=(len(Teff_arr), len(T_arr)))
        Zr = np.zeros(shape=(len(Teff_arr), len(T_arr)))
        for i in range(len(Teff_arr)):
            for j in range(len(T_arr)):
                f = self.yes_starformation(np.array([Teff_axx[i,j], T_axx[i,j]]), set=False)
                Zp[i,j] = f[0] #np.log10(f[0]) if f[0] > 0 else -np.log10(-f[0])
                Zr[i, j] = f[1] #np.log10(f[1]) if f[1] > 0 else -np.log10(-f[1])

        f, ax = plt.subplots(1, 3, figsize=(16, 7), sharey=True, tight_layout=True) #gridspec_kw=dict(hspace=0),
        pcol = ax[0].pcolormesh(Teff_axx, T_axx, Zp, cmap="coolwarm", norm=mpl_col.CenteredNorm(halfrange=1))
        f.colorbar(pcol, ax=ax[0])
        ax[0].set_xlabel(r'$\log{T_{\rm eff}}$')
        ax[0].set_ylabel(r'$\log{T}$')

        pcol = ax[1].pcolormesh(Teff_axx, T_axx, Zr, cmap="coolwarm", norm=mpl_col.CenteredNorm(halfrange=1))
        f.colorbar(pcol, ax=ax[1])
        ax[1].set_xlabel(r'$\log{T_{\rm eff}}$')
        ax[1].set_ylabel(r'$\log{T}$')

        CS = ax[2].contour(Teff_axx, T_axx, Zr, [0.0])
        CS = ax[2].contour(Teff_axx, T_axx, Zp, [0.0], linestyles='--')
        ax[2].set_xlabel(r'$\log{T_{\rm eff}}$')
        ax[2].set_ylabel(r'$\log{T}$')

        if guess is not None:
            for axx in ax:
                axx.axvline(guess[0], c="green", alpha=0.5)
                axx.axhline(guess[1], c="green", alpha=0.5)

        plt.show()

    def no_starformation(self, x, set=True):
        """The system of equations to find the [log10(cs), log10(T)] values, assuming no starformation,
        pulling known quantities from member variables, using self.i as grid index

        Parameters
        ----------
        x: array
            Array of [log10(T), log10(rho)] guesses.
        set: bool, optional (default: True)
            Flag on whether to set member values with current solutions. Set to False when debugging.

        Returns
        -------
        sols: array
            2d array of equilibrium solutions that needs to be zero
        """
        Logcs, LogT = x
        cs, T = 10 ** x

        R = self.R[self.i]
        Omega = self.Omega[self.i]
        Teff4 = self.Teff4[self.i]
        Mdotprime = self.Mdot * (1 - (self.Rmin / R) ** 0.5)

        # Sound speed has a different expression depending on b=0,1
        if not self.b:
            rho = Mdotprime * Omega * Omega / (6 * np.pi * self.alpha * cs * cs * cs)
        else:
            rho = Mdotprime * Omega * Omega * ct.massU / (6 * np.pi * self.alpha * cs * T * ct.Kb)
        cs2 = cs * cs
        Logrho = np.log10(rho)

        LogRopac = Logrho - LogT * 3 + 18
        kappa = 10 ** self.kappaLogTLogRopac(LogT, LogRopac, grid=False)

        tauV = kappa * rho * cs / Omega

        if set:
            self.kappa[self.i] = kappa
            self.tauV[self.i] = tauV
            self.rho[self.i] = rho
            self.cs[self.i] = cs
            self.T[self.i] = T

        Ppgas = (ct.Kb / ct.massU) * T / cs2
        PPrad = tauV * ct.sigmaSB / (2 * ct.c) * Teff4 / (rho * cs2)

        # Rewritten as dimensionless, dividing by cs2
        zero_pressure = Ppgas + PPrad - 1
        # Same, dividing by Teff4
        opacfac = (0.375 * tauV + 0.5 + 0.25 / tauV)
        T4 = T * T * T * T
        T4ratio = Teff4 / T4
        Rr = T4ratio * opacfac
        zero_radiation = Rr - 1

        return np.array([zero_pressure, zero_radiation])

    def yes_starformation(self, x, set=True):
        """ The system of equations to find the [log10(Teff), log10(T)] values, assuming star formation,
        pulling known quantities from member variables, using self.i as grid index

        Parameters
        ----------
        x: array
            Array of [log10(T), log10(eta)] guesses.
        set: bool, optional (default: True)
            Flag on whether to set member values with current solutions. Set to False when debugging.

        Returns
        -------
        sols: array
            2d array of equilibrium solutions that needs to be zero

        """
        LogTeff, LogT = x
        Teff, T = 10 ** x
        Teff4 = Teff ** 4

        Omega = self.Omega[self.i]

        if not self.b:
            rho = self.rho[self.i]
            cs = self.cs[self.i]
        else:
            rho = self.rho[self.i]
            Mdotprime = self.Mdot * (1 - (self.Rmin / self.R[self.i]) ** 0.5)
            self.cs[self.i] = cs = Mdotprime * Omega * Omega * ct.massU / (6 * np.pi * self.alpha * rho * ct.Kb * T)
            self.h[self.i] = cs/Omega

        Logrho = np.log10(rho)
        LogRopac = Logrho - LogT * 3 + 18
        kappa = 10 ** self.kappaLogTLogRopac(LogT, LogRopac, grid=False)
        tauV = kappa * rho * cs / Omega

        cs2 = cs * cs
        Ppgas = (ct.Kb / ct.massU) * T / cs2
        PPrad = tauV * ct.sigmaSB / (2 * ct.c) * Teff4 / (rho * cs2)

        # Rewritten as dimensionless, dividing by cs2
        zero_pressure = Ppgas + PPrad - 1
        # Same, dividing by Teff4
        opacfac = (0.375 * tauV + 0.5 + 0.25 / tauV)
        T4 = T * T * T * T
        T4ratio = Teff4 / T4
        Rr = T4ratio * opacfac
        zero_radiation = Rr - 1

        if set:
            self.Teff4[self.i] = Teff4
            self.T[self.i] = T
            self.kappa[self.i] = kappa
            self.tauV[self.i] = tauV

        return np.array([zero_pressure, zero_radiation])
