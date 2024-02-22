"""
Thompson
Date: 03/06/23

Authors: Alessandro A. Trani, Daria Gangardt

Module that uses the equations from Thompson 2005 to produce a self-gravitating Active Galactic Nuclei disc

Contains ThompsonAGN class, which creates an AGN disc object, documented below.

"""
import numpy as np
from . import constants as ct
from .opacities import electron_scattering_opacity, get_semenov2003_opacity, get_combined_opacity, get_custom_opacity
from scipy.optimize import root
import matplotlib.colors as mpl_col
import matplotlib.pyplot as plt


class ThompsonAGN:
    def __init__(self, Mbh=1e8*ct.MSun, sigma=None, epsilon=1e-3, m=0.2, xi=1., X = 0.7,
                 Rout=200*ct.pc, Mdot_out=320*ct.MSun/ct.yr, Rin=None, opacity="combined", debug=False, xtol=1e-10,
                 rootm="lm"):
        """ Class that creates an AGN disc object using the equations from Thompson et al. 2005.

        Parameters
        ----------
        Mbh: float, optional (default: 1e8*ct.MSun)
            Mass of the Super Massive Black Hole (SMBH) fuelling the AGN in kg.
        sigma: float, optional (default: None)
            Stellar velocity dispersion. If None, calculated using M-Sigma relation. Input in m/s.
        epsilon: float, optional (default: 0.1)
            Star formation radiative efficiency parameter.
        m: float, optional (default: 0.2)
            Fraction of local sound speed, conveys angular momentum transport efficiency in AGN.
        xi: float, optional (default: 1.)
            Fraction of pressure from supernovae in total radiated pressure.
        X: float, optional (default: 0.7)
            Hydrogen abundance of gas in AGN disk.
        Rout: float, optional (default: 200*ct.pc)
            Radius at which accretion rate is no longer constant in outer region of AGN disc.
            If None, set to 1e7 Schwarzchild radii. Input in m.
        Mdot_out: float, optional (default: 320*ct.MSun/ct.yr)
            Accretion rate at Rout. If None, scaled to critical accretion rate. Input is in kg/s.
        Rin: float, optional (default: None)
            Radius of inner boundary of AGN disc. If None, set as 6 gravitational radii of SMBH. Else, input is in m.
        opacity: string/tuple, optional (default: 'combined')
            Which opacity table to use. Options are "semenov", which extrapolates values from Semenov et al. 2003 only,
            "combined", which combines the Semenov et al. 2003 values with the Badnell et al. 2005 values, or a custom
            table of values entered as a tuple made up of three arrays:
             (opacity 2D array, density 1D array, temperature 1D array). For the custom table, inputs must be in SI units.
        debug: bool, optional (default: False)
            Enter debug regime. Code will check root finder operations by printing solutions at every value of r.
        xtol: float, optional (default: 1e-10)
            Tolerance in the root finding method.
        rootm: string, optional (default: 'lm')
            Root finding method to be used for both star and no star formation regimes.
        """
        # Disk parameters
        if Mbh is None and sigma is None:
            raise ValueError("Please provide either a BH mass or a stellar dispersion velocity")
        if Mbh is None:
            print("M from sigma using M-sigma relation")
            self.Mbh = (1.3e8) * ((sigma / 200) ** 4.24) * ct.MSun # using Gultekin 2009 numbers
        else:
            self.Mbh = Mbh
        if sigma is None:
            print("sigma from M using M-sigma relation")
            self.sigma = (200 * 1e3) * (Mbh / (1.3e8*ct.MSun)) ** (1 / 4.24)  # using Gultekin 2009 numbers
        else:
            self.sigma = sigma

        self.epsilon = epsilon
        self.m = m
        self.X = X
        self.xi = xi

        self.Rs = 2 * self.Mbh * ct.G / ct.c ** 2

        if Rout is None:
            self.Rout = 1e7*self.Rs
        else:
            self.Rout = Rout

        if Mdot_out is None:
            self.Mdot_out = 320*(ct.MSun/ct.yr)*(self.Rout/(95*ct.pc)) * (self.sigma/(188e3))**2 #suggested scaling
        else:
            self.Mdot_out = Mdot_out

        if Rin is None:
            self.Rin = 3 * self.Rs  # 3 Schwarzchild radii of SMBH
        else:
            self.Rin = Rin

        # Auxiliary parameters
        self.debug = debug

        self.Rmax = self.Rout #Can be set to higher value than Rout

        # Opacity table variables
        if type(opacity) != str:
            #custom opacity values
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

        # disc shell
        self.i = -1  # current radial index
        self.iswitch = 0  # radial index where star formation ends
        self.rootm = rootm
        self.xtol = xtol

        print("### Thompson et al. 2005 parameters ###")
        print("Mbh = {:e} MSun".format(self.Mbh / ct.MSun))
        print("Mdot_out = {:e} MSun/yr".format(self.Mdot_out * (ct.yr / ct.MSun)))
        print("Rs = {:e} pc".format(self.Rs / ct.pc))
        print("Rin = {:e} Rs".format(self.Rin / self.Rs))
        print("Rout = {:e} Rs = {:e} pc".format(self.Rout / self.Rs, self.Rout/ct.pc))
        print("sigma = {:e} km/s".format(self.sigma / 1e3))
        print("epsilon = ", self.epsilon)
        print("m = ", self.m)
        print("xi = ", self.xi)
        print("Opacity =", self.opacity)

        print("\ndebug =", self.debug)
        print("xtol =", self.xtol)
        print("root method =", self.rootm)

    def solve_disk(self, N=1e4):
        """ Method to evolve the AGN disc, from outer boundary inwards, using Thompson et al. 2005 equations.

        Parameters
        ----------
        N: float, optional (default: 1e4)
            Disc radial resolution.
        """
        N = int(N)

        Rmin = self.Rin
        self.Redge = np.logspace(np.log10(Rmin), np.log10(self.Rout), N + 1)


        self.R = (self.Redge[:-1] + self.Redge[1:]) * 0.5
        self.deltaR = self.Redge[1:] - self.Redge[:-1]

        # Solutions
        self.Omega = np.sqrt(
            ct.G * self.Mbh / (self.R * self.R * self.R) + 2 * (self.sigma * self.sigma) / (self.R * self.R))

        self.rho = np.zeros_like(self.R)
        self.T = np.zeros_like(self.R)
        self.h = np.zeros_like(self.R)
        self.tauV = np.zeros_like(self.R)
        self.kappa = np.zeros_like(self.R)
        self.cs = np.zeros_like(self.R)
        self.Q = np.zeros_like(self.R)
        self.Mdot = np.zeros_like(self.R)
        self.Teff4 = np.zeros_like(self.R)
        self.eta = np.zeros_like(self.R)

        self.Mdot[self.i] = self.Mdot_out  # Set Mdot to constant rate at outer boundary

        # Guessing star formation
        fg0 = ((2 ** 1.5) * self.Mdot_out * ct.G
               / (self.m * (self.Omega[self.i] * self.Omega[self.i] * self.Omega[self.i] * self.R[self.i] * self.R[
                    self.i] * self.R[self.i]))) ** 0.5  # gas fraction at boundary
        T0 = np.sqrt(fg0 * self.sigma * self.sigma / self.R[self.i]) * (
                (3 * ct.c) ** (0.25) * ((2 ** 3.5) * np.pi * ct.G * ct.sigmaSB) ** (
            -0.25))  # optically thick limit for temperature
        rho0 = self.Omega[self.i] ** 2 / (np.sqrt(2) * np.pi * ct.G)

        x_guess = np.log10([T0, rho0]) #initial guess

        h0 = self.Rout*fg0/(2**1.5)
        self.Mdotmax = 8 * np.pi * rho0 * h0 * self.sigma * self.sigma * self.Rout / (self.epsilon * ct.c)
        if self.Mdot_out > self.Mdotmax:
            print("Maximum accretion rate: ", self.Mdotmax * ct.yr / ct.MSun, "MSun / yr")
            print("Outer boundary accretion rate: ", self.Mdotmax * ct.yr / ct.MSun, "MSun / yr")
            raise ValueError("Accretion rate at outer boundary higher than maximum, no disk would form. "
                             "Please enter a smaller Mdot_out.")
            exit()

        if self.debug: print("x_guess, log10(T) and log10(rho) at outermost R:", x_guess)

        print(" ### Beginning integration at Rmax ###")

        #root solver options
        root_options = {'col_deriv': True, 'factor': 0.1, 'xtol': self.xtol, 'ftol': 0.01}
        if self.rootm == "hybr": root_options['maxfev'] = 10000

        self.nostar = False
        self.R_AGN = 0.
        self.Mdisk = 0.

        # Outside in, no star formation
        for self.i in range(len(self.R) - 1, -1, -1):
            # Firstly, no star formation guess
            Teff4 = 3 * self.Mdot[self.i] * (1 - np.sqrt(self.Rin / self.R[self.i])) * self.Omega[self.i] * self.Omega[
                self.i] / (8 * np.pi)
            Teff4 = Teff4 / ct.sigmaSB
            self.Teff4[self.i] = Teff4
            sol = root(self.no_starformation, x_guess, method=self.rootm,
                       options=root_options,)

            if sol.success:
                # Calculate Q
                self.Q[self.i] = self.Omega[self.i] ** 2 / (np.sqrt(2) * np.pi * ct.G * self.rho[self.i])
                if self.Q[self.i] <= 1:
                    #enter star formation regime
                    self.Q[self.i] = 1
                    if self.i == len(self.R) - 1:
                        # Here we find the guess if we are at the boundary using optically thick solutions
                        Logrho0 = np.log10(rho0)
                        LogT0 = np.log10(T0)
                        LogRopac0 = Logrho0 - LogT0 * 3 + 18

                        kappa0 = 10 ** self.kappaLogTLogRopac(LogT0, LogRopac0, grid=False)
                        tauV0 = kappa0 * self.sigma * self.sigma * fg0 / (2 * np.pi * ct.G * self.R[self.i])

                        Teff4 = ((4 / 3) * (T0 ** 4)) / (tauV0 + (2 / (3 * tauV0)) + (4 / 3))
                        Sigmastar0 = 2 * ct.sigmaSB * Teff4 / (self.epsilon * ct.c * ct.c)

                        h0 = np.sqrt(self.Mdot_out / (4 * np.pi * self.R[self.i] * self.Omega[self.i] * self.m * rho0))
                        eta0 = Sigmastar0 / (2 * rho0 * h0 * self.Omega[self.i])

                        x_guess = np.array([np.log10(T0), np.log10(eta0)])
                    elif self.eta[self.i + 1] < 1e-8:
                        #ensure no log10(0) errors
                        x_guess = np.array([np.log10(self.T[self.i + 1]), -8])
                    else:
                        x_guess = np.array([np.log10(self.T[self.i + 1]), np.log10(self.eta[self.i + 1])])

                    sol = root(self.yes_starformation, x_guess, method=self.rootm,
                               options=root_options,)
                    if sol.success:
                        #set values for next r solutions
                        if self.debug:
                            print(
                                "Looking for roots in pressure equilibrium, thermal equilibrium equations; yes star formation")
                            sol_check_zeros = self.yes_starformation(sol.x)
                            print("Zero residuals at i =", self.i, " :", sol_check_zeros)
                        # Calculate the next Mdot
                        Macc = 4 * np.pi * self.R[self.i + 1:] * self.rho[self.i + 1:] * self.h[self.i + 1:] \
                               * self.Omega[self.i + 1:] * self.eta[self.i + 1:] * self.deltaR[self.i + 1:]
                        Macc = Macc.sum()
                        self.Mdot[self.i - 1] = self.Mdot_out - Macc
                        self.Mdisk += 2 * np.pi * 2 * self.h[self.i] * self.deltaR[self.i] * self.R[self.i] * self.rho[self.i]

                        # Next no_starformation guess is the current starformation temperature
                        x_guess = np.array(np.log10([self.T[self.i], self.rho[self.i]]))
                    else:
                        print(f"No yes_starformation solution found at i = {self.i}, log10(R/Rs) = {np.log10(self.R[self.i]/self.Rs)}. Plotting roots.")
                        self.plot_roots_yes_star(guess = sol.x)
                        self.plot()
                        exit()
                else:
                    #definitely in no star formation regime, set values for next r value
                    if self.debug:
                        print(
                            "Looking for roots in pressure equilibrium, thermal equilibrium equations; no star formation")
                        sol_check_zeros = self.no_starformation(sol.x)
                        print("Zero residuals at i =", self.i, " :", sol_check_zeros)
                    self.nostar = True
                    if not self.R_AGN:
                        #set radius at which we enter non star forming regime
                        self.R_AGN = self.R[self.i]
                        self.iswitch = self.i
                        print("### Switching to no star formation regime at i = {}, R = {}Rs ###".format(self.i, self.R_AGN/self.Rs))
                    # Here if Q>1, we keep using the previous solution
                    if self.i != 0:
                        self.Mdot[self.i - 1] = self.Mdot[self.i]  # No starformation, the next Mdot remains the same
                    self.Mdisk += 2 * np.pi * 2 * self.h[self.i] * self.deltaR[self.i] * self.R[self.i] * self.rho[
                        self.i]
                    x_guess = np.log10([self.T[self.i], self.rho[self.i]])
            else:
                print(f"No no_starformation solution found at i = {self.i}, "
                      f"log10(R/Rs) = {np.log10(self.R[self.i] / self.Rs)}. Plotting roots.")
                self.plot_roots_no_star(guess = sol.x)
                self.plot()
                exit()
            print("{:d}/{:d}".format(self.i, N), end="\r")

        #equations can be solved outside of Rout using the optically thick, constant accretion rate approximation
        if self.Rmax > self.Rout:
            R_out = np.logspace(np.log10(self.Rout), np.log10(self.Rmax), int(N/2))
            Omega_out = np.sqrt(ct.G * self.Mbh / (R_out* R_out * R_out) + 2 * (self.sigma * self.sigma) / (R_out * R_out))
            Q_out = np.ones_like(R_out)
            rho_out = (Omega_out**2)/(np.sqrt(2)*np.pi*ct.G*Q_out)
            fg_out = ((2 ** 1.5) * self.Mdot_out * ct.G
                   / (Q_out*self.m * (Omega_out * Omega_out* Omega_out * R_out * R_out * R_out))) ** 0.5  # gas fraction at boundary
            h_out = R_out*fg_out*Q_out/(2**1.5)
            cs_out = self.sigma*fg_out*Q_out/2
            Mdot_out = self.Mdot_out*np.ones_like(R_out)

            T_out = np.sqrt(fg_out * self.sigma * self.sigma / R_out)* (
                (3 * ct.c* Q_out)**0.25) * (((2 ** 3.5) * np.pi * ct.G * ct.sigmaSB) ** -0.25)  # optically thick limit for temperature

            Logrho0 = np.log10(rho_out)
            LogT0 = np.log10(T_out)
            LogRopac0 = Logrho0 - LogT0 * 3 + 18
            kappa_out = 10 ** self.kappaLogTLogRopac(LogT0, LogRopac0, grid=False)
            tauV_out = kappa_out * self.sigma * self.sigma * fg_out / (2 * np.pi * ct.G * R_out)

            Teff4_out = ((4 / 3) * (T_out ** 4)) / (tauV_out + (2 / (3 * tauV_out)) + (4 / 3))
            Sigmastar_out = 2 * ct.sigmaSB * Teff4_out / (self.epsilon * ct.c * ct.c)
            eta_out = Sigmastar_out / (2 * rho_out * h_out * Omega_out)

            self.Omega = np.concatenate([self.Omega, Omega_out])
            self.R  = np.concatenate([self.R, R_out])
            self.rho = np.concatenate([self.rho, rho_out])
            self.T = np.concatenate([self.T, T_out])
            self.h = np.concatenate([self.h, h_out])
            self.tauV = np.concatenate([self.tauV, tauV_out])
            self.kappa = np.concatenate([self.kappa, kappa_out])
            self.cs = np.concatenate([self.cs, cs_out])
            self.Q = np.concatenate([self.Q, Q_out])
            self.Mdot =np.concatenate([self.Mdot, Mdot_out])
            self.Teff4 = np.concatenate([self.Teff4, Teff4_out])
            self.eta = np.concatenate([self.eta, eta_out])

        self.Mc = 4 * np.pi * self.R * self.R * self.h * self.rho * self.eta * self.Omega
        self.Mdotmax = 8 * np.pi * self.rho * self.h * self.sigma * self.sigma * self.R / (self.epsilon * ct.c)
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

    def plot_mdot(self):
        """ Simple accretion rate plotting method, can be used to check that disk has been evolved correctly and that
        the accretion rate is high enough to form a luminous AGN disk
        """

        f, ax = plt.subplots(1, 1, figsize=(10, 5), tight_layout=True)
        Mstar = 2 * np.pi * self.R * self.R * self.rho * self.h * self.Omega * self.eta
        Mc = self.Mc
        Medd = self.L_Edd()/(0.1*ct.c*ct.c)

        print(" ### Checking Accretion Rates ###")
        print("Mdot_Edd = {:e} Msun per year".format(Medd * ct.yr / ct.MSun))
        print("Mdot_c (r = Rout) = {:e} Msun per year \nMdot_out = {:e} Msun per year".format(
            Mc[-1]*ct.yr/ct.MSun, self.Mdot_out*ct.yr/ct.MSun))
        print("Mdot (r = Rin) = {:e} Msun per year = {:e} Mdot_Edd".format(self.Mdot[0] * ct.yr / ct.MSun, self.Mdot[0]/Medd))

        ax.plot(np.log10(self.R / self.Rs), np.log10(self.Mdot/Medd), label = r'$\dot{M}_{\rm total}$ ')
        ax.plot(np.log10(self.R[Mstar>0] / self.Rs), np.log10(Mstar[Mstar >0]/Medd), label =  r'$\dot{M}_{\rm star}$')
        ax.plot(np.log10(self.R[Mstar>0]/ self.Rs), np.log10(Mc[Mstar>0]/Medd),
                label=r'$\dot{M}_{\rm C}$')
        Mdotmax = 8*np.pi*self.rho*self.h*self.sigma*self.sigma*self.R/(self.epsilon*ct.c)
        ax.plot(np.log10(self.R[self.R > self.R_AGN]/ self.Rs), np.log10(Mdotmax[self.R > self.R_AGN]/ Medd), label=r'$\dot{M}_{\rm max}$')
        ax.axvline( np.log10(self.R_AGN/self.Rs), -100, 100, c = 'k')

        ax.set_ylabel(r'$\log_{10}{\dot{M}/\dot{M}_{\rm Edd}}$')
        ax.set_xlabel(r'$\log_{10}{R/R_{\rm S}}$')
        ax.legend()
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
            print("guess = (log rho, log T) = ({:} , {:})".format(guess[1], guess[0]))
        else:
            print("guess = {}".format(guess))
        print("zoom guess = {}".format(zoomguess))
        print("###")


        #zoom in on guess
        if zoomguess:
            rho_arr = np.linspace(guess[0] - 1, guess[0] + 1, 100)
            T_arr = np.linspace(guess[1] - 1, guess[1] + 1, 100)
        else:
            rho_arr = np.linspace(-30, 10., 200)
            T_arr = np.linspace(0, 10., 200)

        T_axx, rho_axx = np.meshgrid(T_arr, rho_arr)
        Zp = np.zeros(shape=(len(T_arr), len(rho_arr)))
        Zr = np.zeros(shape=(len(T_arr), len(rho_arr)))
        for i in range(len(rho_arr)):
            for j in range(len(T_arr)):
                f = self.no_starformation(np.array([T_axx[i,j], rho_axx[i,j]]), set=False)
                Zp[i,j] = f[0]
                Zr[i, j] = f[1]


        f, ax = plt.subplots(1, 3, figsize=(20, 8), sharey=True, sharex =True, tight_layout=True) #gridspec_kw=dict(hspace=0),
        pcol = ax[0].pcolormesh(T_axx, rho_axx, Zp, cmap="coolwarm", norm=mpl_col.CenteredNorm(halfrange=1))
        f.colorbar(pcol, ax=ax[0], label = r'$\mathrm{Balancing \, Pressure}$')
        ax[0].set_xlabel(r'$\log{T}$')
        ax[0].set_ylabel(r'$\log{\rho}$')

        pcol = ax[1].pcolormesh(T_axx, rho_axx,  Zr, cmap="coolwarm", norm=mpl_col.CenteredNorm(halfrange=1))
        f.colorbar(pcol, ax=ax[1], label = r'$\mathrm{Balancing \, Temperatures}$')
        ax[1].set_xlabel(r'$\log{T}$')
        ax[1].set_ylabel(r'$\log{\rho}$')

        CS = ax[2].contour(T_axx, rho_axx,  Zr, [0.0])
        CS = ax[2].contour(T_axx, rho_axx, Zp, [0.0], linestyles = '--')
        ax[2].set_xlabel(r'$\log{T}$')
        ax[2].set_ylabel(r'$\log{\rho}$')

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
        if guess is None: zoomguess = False

        print("### Plotting yes star formation regime solutions ###")
        print("i = {:} ".format(self.i))
        print("R = {:e} Rs = {:e} pc".format(self.R[self.i] / self.Rs, self.R[self.i] / ct.pc, ))
        if guess is not None:
            print("guess = (log eta, log T) = ({:} , {:})".format(guess[1], guess[0]))
        else:
            print("guess = {}".format(guess))
        print("zoom guess = {}".format(zoomguess))
        print("###")

        if zoomguess:
            T_arr = np.linspace(guess[0] - 1, guess[0] + 1, 100)
            eta_arr = np.linspace(guess[1] - 1, guess[1] + 1, 100)
        else:
            T_arr = np.linspace(1., 6., 200)
            eta_arr = np.linspace(1., -5., 200)

        T_axx, eta_axx = np.meshgrid(T_arr, eta_arr)
        Zp = np.zeros(shape=(len(T_arr), len(eta_arr)))
        Zr = np.zeros(shape=(len(T_arr), len(eta_arr)))
        for i in range(len(T_arr)):
            for j in range(len(eta_arr)):
                f = self.yes_starformation(np.array([T_axx[i, j], eta_axx[i, j], ]), set=False)
                Zp[i, j] = f[0]  # np.log10(f[0]) if f[0] > 0 else -np.log10(-f[0])
                Zr[i, j] = f[1]  # np.log10(f[1]) if f[1] > 0 else -np.log10(-f[1])

        f, ax = plt.subplots(1, 3, figsize=(16, 7), sharey=True, tight_layout=True)  # gridspec_kw=dict(hspace=0),
        pcol = ax[0].pcolormesh(T_axx, eta_axx,  Zp, cmap="coolwarm", norm=mpl_col.CenteredNorm(halfrange=1))
        f.colorbar(pcol, ax=ax[0], label = r'$\mathrm{Balancing \, Pressures}$')
        ax[0].set_xlabel(r'$\log{T}$')
        ax[0].set_ylabel(r'$\log{eta}$')

        pcol = ax[1].pcolormesh(T_axx, eta_axx, Zr, cmap="coolwarm", norm=mpl_col.CenteredNorm(halfrange=1))
        f.colorbar(pcol, ax=ax[1], label = r'$\mathrm{Balancing \, Temperatures}$')
        ax[1].set_xlabel(r'$\log{T}$')
        ax[1].set_ylabel(r'$\log{eta}$')

        CS = ax[2].contour(T_axx, eta_axx, Zr, [0.0])
        CS = ax[2].contour(T_axx, eta_axx, Zp, [0.0], linestyles='--')
        ax[2].scatter(T_axx[amin_res], eta_axx[amin_res])
        ax[2].set_xlabel(r'$\log{T}$')
        ax[2].set_ylabel(r'$\log{eta}$')

        if guess is not None:
            for axx in ax:
                axx.axvline(guess[0], c="green", alpha=0.5)
                axx.axhline(guess[1], c="green", alpha=0.5)

        plt.show()

    def no_starformation(self, x, set=True):
        """The system of equations to find the [log10(T), log10(rho)] values, assuming no starformation,
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
        LogT, Logrho = x
        T, rho = 10 ** x
        Omega = self.Omega[self.i]
        Teff4 = self.Teff4[self.i]
        Mdot = self.Mdot[self.i]
        r = self.R[self.i]

        LogRopac = Logrho - LogT * 3 + 18
        kappa = 10 ** self.kappaLogTLogRopac(LogT, LogRopac, grid=False)
        h = np.sqrt(Mdot / (4 * np.pi * r * Omega * self.m * rho))
        tauV = kappa * rho * h

        opacfac = 0.75 * (tauV + 4 / 3 + 2 / (3 * tauV))
        T4 = T * T * T * T
        T4ratio = Teff4 / T4

        Rr = T4ratio * opacfac
        zero_radiation = Rr - 1
        press_factor = rho * h * h * Omega * Omega
        Pa = ct.Kb * T / (ct.massU * h * h * Omega * Omega)
        Pb = ct.sigmaSB * Teff4 * tauV / (ct.c * press_factor)
        zero_pressure = Pa + Pb - 1
        sols = np.array([zero_pressure, zero_radiation])

        if set:
            self.kappa[self.i] = kappa
            self.tauV[self.i] = tauV
            self.rho[self.i] = rho
            self.h[self.i] = h
            self.T[self.i] = T
            self.cs[self.i] = h * Omega

        return sols

    def yes_starformation(self, x, set=True):
        """ The system of equations to find the [log10(T), log10(eta)] values, assuming star formation,
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
        LogT, Logeta = x
        T, eta = 10 ** x

        Omega = self.Omega[self.i]
        Mdot = self.Mdot[self.i]
        rho = (Omega ** 2) / (np.sqrt(2) * np.pi * ct.G)
        h = np.sqrt(Mdot / (4 * np.pi * self.R[self.i] * Omega * self.m * rho))
        Mdotprime = Mdot * (1 - np.sqrt(self.Rin / self.R[self.i]))

        Logrho = np.log10(rho)
        LogRopac = Logrho - LogT * 3 + 18

        kappa = 10 ** self.kappaLogTLogRopac(LogT, LogRopac, grid=False)
        tauV = kappa * rho * h

        opacfac = 0.75 * (tauV + 4 / 3 + 2 / (3 * tauV))
        T4 = T * T * T * T
        Teff4 = T4 / opacfac

        factorT = T * ct.Kb / ct.massU
        zero_pressure = 1 + 2 * h * eta * Omega * self.epsilon * ct.c * (
                tauV / 2 + self.xi) / factorT - h * h * Omega * Omega / factorT

        factorTeff = ct.sigmaSB * Teff4
        zero_radiation = 1 - rho * h * Omega * eta * self.epsilon * (ct.c ** 2) / factorTeff - \
                         3 * Omega * Omega * Mdotprime / (8 * np.pi * factorTeff)
        if set:
            self.T[self.i], self.eta[self.i] = T, eta
            self.Teff4[self.i] = Teff4
            self.kappa[self.i] = kappa
            self.tauV[self.i] = tauV
            self.rho[self.i] = rho
            self.h[self.i] = h
            self.cs[self.i] = h * Omega
        return np.array([zero_pressure, zero_radiation])
