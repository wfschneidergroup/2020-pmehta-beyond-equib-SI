import numpy as np
from ase import units
from scipy.optimize import fsolve
from mpmath import mp, findroot
from scipy.integrate import odeint

J2eV = 6.24150974E18            # eV/J
h = 6.626068E-34 * J2eV
wave2eV = 0.00012398

R = 82.057338  # cm3-atm / K / mol
Na = 6.022140e23  # 1 / mol


class simpleMkm:
    """
    Class to run a simple microkinetic model for,
    A2 + B <-> 2AB
    Parametrized for ammonia synthesis with
    A = N, and B = 3/2 H2
    Elementary rate equations described in main text    
    """

    def __init__(self,
                 T,
                 EA,
                 E_rxn=-0.9514,
                 S=[0.001986,
                    0.001354,
                    0.001998],
                 stoich=[-1, -3, 2],
                 A2_frac=1/4.,
                 total_site_density=1.e-3,  # mols / cm3
                 dps=25,
                 plasma_on=False,
                 model='site',
                 concentration_based=False,
                 nu0=10./60.,  # volumetric flow (assuming forward rxn) cm3/ s
                 V=1.0,  # reactor volume cm3
                 k_eimpact_A2=1e-8,  # cm3 / s. For A2 activation
                 k_eimpact_AB=0.0,  # cm3 / s. For AB activation
                 n_e=1e12,  # electron number density, 1/cm3
                 rxn2_barrier=False,  # Add a barrier for A* + B <--> AB?
                 excite_type='diss', # diss / vib (no / reduced barrier for excited state)
                 Evib_A2=0.3, # Barrier reduced by
                 Evib_AB=0.3,
                 p0=1.0,
                 reverse=False,
                 include_plasma_rxns=False,
                 **kwargs):

        # Extra kwargs passed to enable key sharing between
        # simpleMkm and MkmRunner

        if dps is not None:
            mp.dps = dps

        self.T = T
        self.EA = EA
        self.p0 = p0

        self.A2_frac = A2_frac
        self.C_0 = 1. / R / T
        self.total_site_density = total_site_density
        self.V = V
        self.stoich = stoich

        self._nu0 = nu0
        self._pA2_0 = A2_frac * p0

        self.reverse = reverse

        if not reverse:
            self.pA2_0 = A2_frac * p0
            self.pB_0 = (1 - A2_frac) * p0
            self.pAB_0 = 0.0
            self.nu0 = self.get_vol_flowrate(0.0)

        else:
            self.pA2_0 = 0.0
            self.pB_0 = 0.0
            self.pAB_0 = 1.0
            self.nu0 = self.get_vol_flowrate(1.0)

        self.n_e = n_e
        self.C_e = n_e / Na
        self.k_eimpact_A2 = k_eimpact_A2
        self.k_eimpact_AB = k_eimpact_AB

        self.model = model
        self.include_plasma_rxns = include_plasma_rxns

        self.concentration_based = concentration_based
        self.rxn2_barrier = rxn2_barrier
        self.excite_type = excite_type
        if excite_type == 'vib':
            self.Evib_A2 = Evib_A2
            self.Evib_AB = Evib_AB

        self.E_rxn = E_rxn
        self.S = S
        self.dS_rxn = np.dot(np.array(S), np.array(stoich))
        self.G_rxn = E_rxn - T * self.dS_rxn
        self.Keq = np.exp(-self.G_rxn / units.kB / T)
        self.plasma_on = plasma_on
        self.dE = self.get_rxn_energies()
        self.dS = self.get_rxn_entropies()
        self.Ea = self.get_Eacts()
        self.dSTS = self.get_TS_entropies()

        kf, kr = self.get_rate_constants()
        self.kf = kf
        self.kr = kr

        if concentration_based:
            kf_conc, kr_conc = self.convert_rate_constants()
            self.kf_conc = kf_conc
            self.kr_conc = kr_conc

    def get_pressures(self, X):
        """
        Given a conversion and an initial
        pressure, report the new pressures.

        Reported converseions are based on the
        forward reaction, with no initial concentration
        of product
        """
        if isinstance(X, tuple):
            X = X[0]

        p0 = self.p0
        stoich = self.stoich

        A2_frac = self.A2_frac * p0
        B_frac = (1 - A2_frac) * p0

        nA2 = A2_frac + stoich[0] * X * A2_frac
        nB = B_frac + stoich[1] * X * A2_frac
        nAB = stoich[2] * X * A2_frac

        ntot = nA2 + nB + nAB

        pA2 = p0 * nA2 / ntot
        pB = p0 * nB / ntot
        pAB = p0 * nAB / ntot

        return [pA2, pB, pAB]

    def get_conversion(self, pA2, pA2_0=1./4., pA2_prime=0.0,
                       pAB_prime=0.0, X0=[0.1, 0.01, 0.01]):
        """
        Given the final and inital pressures
        report the conversion
        """
        if not self.plasma_on:
            def obj(X):
                return pA2 - 4. * pA2_0 * (1. - X) / (4. - 2. * X)
            if not self.reverse:
                X = findroot(obj, 0.005)
            else:
                X = findroot(obj, 0.995)

        else:
            def obj(X):
                x, y, z = X
                f1 = pA2 - 4 * pA2_0 * (1. - x - y - z) / (4. - 2 * x - 2 * z)
                f2 = pA2_prime - 4 * pA2_0 * y / (4. - 2 * x - 2 * z)
                f3 = pAB_prime - 4 * pA2_0 * 2 * z / (4. - 2 * x - 2 * z)

                return [f1, f2, f3]

            X = fsolve(obj, X0)

        return X

    def get_vol_flowrate(self, X):
        """
        Given a conversion,
        get the new volumetric flow rate
        """
        return self._nu0 * (4. - 2. * X) / 4.

    def get_concentrations(self, Ps):
        return list(np.array(Ps) / R / self.T)

    @np.vectorize
    def get_molar_flowrate(self, P, nu):
        """
        Given a pressure, and volumetric flow rate, get the molar flowrate
        """
        F = P * nu / R / self.T
        return F

    @staticmethod
    def _eqb_objective_func(X, self):

        pA2, pB, pAB = self.get_pressures(X)
        Keq = self.Keq

        return Keq - pAB ** 2 / pA2 / pB ** 3

    def get_eqb_conversion(self, X0=0.1):

        def get_findroot_eqns(*args):
            return simpleMkm._eqb_objective_func(args, self)

        X = fsolve(get_findroot_eqns, X0)
        return X

    def get_rxn_energies(self):
        EA = self.EA
        dE = []

        dE = np.zeros(2)
        dE[0] = 2 * EA
        dE[1] = (self.E_rxn - dE[0]) / 2.

        return dE

    def get_rxn_entropies(self):

        S = self.S
        dS = []

        dS = np.zeros(2)
        dS[0] = 0 - S[0]
        dS[1] = S[2] - 3. / 2. * S[1]

        return dS

    def get_Eacts(self):

        Ea = []
        dE = self.dE

        Ea = np.zeros(2)

        Ea[0] = 1.57 / 2. * dE[0] + 1.56

        # From Grabow
        if self.rxn2_barrier:
            Ea[1] = -0.195 * dE[0] + 1.24
        return Ea

    def get_TS_entropies(self):

        dSTS = []
        dS = self.dS

        dSTS = np.zeros(2)
        # Assume all entropy is lost at TS
        dSTS[0] = dS[0]
        # Assume TS is initial state like
        dSTS[1] = 0

        return dSTS

    def get_rate_constants(self):
        T = self.T
        kbT = units.kB * T
        dE = self.dE
        dS = self.dS
        Ea = self.Ea
        dSTS = self.dSTS

        K = np.zeros(len(dE))           # equilibrium constants
        kf = np.zeros(len(dE))                # forward rate constants
        kr = np.zeros(len(dE))             # reverse rate constants

        for i in range(len(dE)):
            dG = dE[i] - T * dS[i]
            K[i] = mp.exp(-dG / kbT)

            Ea[i] = max([0, dE[i], Ea[i]])

            kf[i] = kbT / h * mp.exp(dSTS[i] / units.kB) * mp.exp(-Ea[i] / kbT)
            kr[i] = kf[i] / K[i]

        return kf, kr


    def convert_rate_constants(self):
        """
        The goal here is to convert rate constants
        from having units of 1 /s to rate constants compatible
        with concentration based expressions
        """

        C_0 = self.C_0
        total_site_density = self.total_site_density

        kf = self.kf
        kr = self.kr

        kf_conc = kf.copy()
        kr_conc = kr.copy()

        kf_conc[0] = kf_conc[0] / C_0 / total_site_density
        kr_conc[0] = kr_conc[0] / total_site_density

        kf_conc[1] = kf_conc[1] / C_0 ** 1.5
        kr_conc[1] = kr_conc[1] / C_0

        if self.plasma_on:

            # A2 + e -> A2' + e
            # Note that these units are compatible with electron number density
            # and not with electron molar concentration
            kf2 = self.k_eimpact_A2
            kr2 = 0.0

            kf4 = self.k_eimpact_AB
            kr4 = 0.0

            kbT = units.kB * self.T

            if self.excite_type == 'diss':
                # A2' + 2* -> 2A*
                kf3 = kbT / h * mp.exp(self.dSTS[0] / units.kB) \
                      / C_0 / total_site_density
                kr3 = 0.0

                # AB' + * -> A* + B
                kf5 = kbT / h / C_0 * mp.exp(-self.dS[1] / units.kB)
                kr5 = 0.0

            elif self.excite_type == 'vib':
                dE = self.dE
                Ea = self.Ea
                dS = self.dS
                dSTS = self.dSTS

                dE3 = dE[0] - self.Evib_A2
                Ea3 = Ea[0] - self.Evib_A2

                dG3 = dE3 - self.T * dS[0]
                K3 = mp.exp(- dG3 / kbT)

                Ea3 = max([0, dE3, Ea3])
                kf3 = kbT / h * mp.exp(dSTS[0] / units.kB) \
                      * mp.exp(-Ea3 / kbT)

                kr3 = kf3 / K3

                kf3 = kf3 / C_0 / total_site_density
                kr3 = kr3 / total_site_density

                kf4 = 0.0
                kr4 = 0.0

                kf5 = 0.0
                kr5 = 0.0

            kf_conc = np.append(kf_conc, [kf2, kf3, kf4, kf5])

            # Both these steps are assumed to be irreversible
            kr_conc = np.append(kr_conc, [kr2, kr3, kr4, kr5])

        return kf_conc, kr_conc

    def get_rates(self, theta, Ps):
        """
        Uses rate constants in 1 / s
        """
        if not self.plasma_on:
            # if not self.concentration_based:
                kf = self.kf
                kr = self.kr

                if isinstance(theta, tuple):
                    tA = theta[0]
                else:
                    tA = theta

                pA2, pB, pAB = Ps
                tstar = 1 - tA
                rate = np.zeros(2)

                rate[0] = kf[0] * pA2 * tstar ** 2 - kr[0] * tA ** 2
                rate[1] = kf[1] * tA * pB ** 1.5 - kr[1] * pAB * tstar
        return rate

    def get_rates_conc(self, Cs, rate_type='net'):
        """
        Uses rate constants in concentration based units
        """
        
        kf = self.kf_conc
        kr = self.kr_conc

        if not self.plasma_on:
            CA2, CB, CAB, CAstar = Cs
            Cstar = self.total_site_density - CAstar
            rate = np.zeros(2)

            if rate_type == 'net':
                rate[0] = kf[0] * CA2 * Cstar ** 2 - kr[0] * CAstar ** 2
                rate[1] = kf[1] * CAstar * CB ** 1.5 - kr[1] * CAB * Cstar

            if rate_type == 'forward':
                rate[0] = kf[0] * CA2 * Cstar ** 2
                rate[1] = kf[1] * CAstar * CB ** 1.5

            if rate_type == 'reverse':
                rate[0] = kr[0] * CAstar ** 2
                rate[1] = kr[1] * CAB * Cstar

        elif self.plasma_on:
            rate = np.zeros(6)
            CA2, CA2_prime, CB, CAB, CAB_prime, CAstar = Cs
            Cstar = self.total_site_density - CAstar

            if rate_type == 'net':
                rate[0] = kf[0] * CA2 * Cstar ** 2 - kr[0] * CAstar ** 2
                rate[1] = kf[1] * CAstar * CB ** 1.5 - kr[1] * CAB * Cstar
                rate[2] = kf[2] * self.n_e * CA2 - kr[2] * CA2_prime * self.n_e
                rate[3] = kf[3] * CA2_prime * Cstar ** 2 - kr[3] * CAstar ** 2
                rate[4] = kf[4] * self.n_e * CAB - kr[4] * CAB_prime * self.n_e
                rate[5] = kf[5] * CAB_prime * Cstar \
                          - kr[5] * CAstar * CB ** 1.5

            elif rate_type == 'forward':
                rate[0] = kf[0] * CA2 * Cstar ** 2
                rate[1] = kf[1] * CAstar * CB ** 1.5
                rate[2] = kf[2] * self.n_e * CA2
                rate[3] = kf[3] * CA2_prime * Cstar ** 2
                rate[4] = kf[4] * self.n_e * CAB
                rate[5] = kf[5] * CAB_prime * Cstar

            elif rate_type == 'reverse':
                rate[0] = kr[0] * CAstar ** 2
                rate[1] = kr[1] * CAB * Cstar
                rate[2] = kr[2] * CA2_prime * self.n_e
                rate[3] = kr[3] * CAstar ** 2
                rate[4] = kr[4] * CAB_prime * self.n_e
                rate[5] = kr[5] * CAstar * CB ** 1.5

        return rate

    @staticmethod
    def get_odes(theta, t, X, self):
        """
        ODEs for site-based model
        """
        if not self.plasma_on:
            if self.model == 'site':

                Ps = self.get_pressures(X)
                rate = self.get_rates(theta, Ps)

                dthetaAdt = 2 * rate[0] - rate[1]
                return dthetaAdt

    @staticmethod
    def get_odes_CSTR(Concs, t, self):
        """
        ODEs for CSTR model
        """
        if self.model == 'reactor':
            T = self.T
            V = self.V
            nu0 = self.nu0
            sites = self.total_site_density

            if self.reverse:
                X = 1.0
            else:
                X = 0.0

            pA2_0 = self.pA2_0
            pB_0 = self.pB_0
            pAB_0 = self.pAB_0

            CA2_0 = pA2_0 / R / T
            CB_0 = pB_0 / R / T
            CAB_0 = pAB_0 / R / T

            if not self.concentration_based:

                CA2, CB, CAB, thetaA = Concs

                pA2 = CA2 * R * T
                pB = CB * R * T
                pAB = CAB * R * T

                rate = self.get_rates(thetaA, [pA2, pB, pAB])

                X = self.get_conversion(pA2, self._pA2_0)
                # print X
                nu = self.get_vol_flowrate(X)

                dthetaAdt = 2 * rate[0] - rate[1]
                dCA2dt = (CA2_0 * nu0 - CA2 * nu - rate[0] * sites * V) / V
                dCBdt = (CB_0 * nu0 - CB * nu - 1.5 * rate[1] * sites * V) / V
                dCABdt = (CAB_0 * nu0 - CAB * nu + rate[1] * sites * V) / V

                return [dCA2dt, dCBdt, dCABdt, dthetaAdt]

            if self.concentration_based:
                if not self.plasma_on:
                    CA2, CB, CAB, CAstar = Concs

                    pA2 = CA2 * R * T
                    pB = CB * R * T
                    pAB = CAB * R * T

                    rate = self.get_rates_conc(Concs)

                    # Conversion is always measured in foward direction
                    X = self.get_conversion(pA2, self._pA2_0)
                    # print X, pA2
                    nu = self.get_vol_flowrate(X)

                    dCAstardt = 2 * rate[0] - rate[1]
                    dCA2dt = (CA2_0 * nu0 - CA2 * nu - rate[0] * V) / V
                    dCBdt = (CB_0 * nu0 - CB * nu - 1.5 * rate[1] * V) / V
                    dCABdt = (CAB_0 * nu0 - CAB * nu + rate[1] * V) / V

                    return [dCA2dt, dCBdt, dCABdt, dCAstardt]

                elif self.plasma_on:
                    CA2_prime_0 = 0.
                    CAB_prime_0 = 0.

                    CA2, CA2_prime, CB, CAB, CAB_prime, CAstar = Concs

                    pA2 = CA2 * R * T
                    pB = CB * R * T
                    pAB = CAB * R * T
                    pA2_prime = CA2_prime * R * T
                    pAB_prime = CAB_prime * R * T

                    rate = self.get_rates_conc(Concs)

                    X = self.get_conversion(pA2, self._pA2_0,
                                            pA2_prime=pA2_prime,
                                            pAB_prime=pAB_prime)

                    nu = self.get_vol_flowrate(X[0])

                    dCAstardt = 2 * (rate[0] + rate[3]) - rate[1] + rate[5]

                    dCA2dt = (CA2_0 * nu0 - CA2 * nu) / V \
                             - rate[0] - rate[2]

                    dCA2_primedt = (CA2_prime_0 * nu0 - CA2_prime * nu) / V \
                                    + rate[2] - rate[3]

                    dCBdt = (CB_0 * nu0 - CB * nu) / V \
                            - (rate[1] + rate[5]) * 1.5

                    dCABdt = (CAB_0 * nu0 - CAB * nu) / V \
                             + rate[1] - rate[4]

                    dCAB_primedt = (CAB_prime_0 * nu0 - CAB_prime * nu) / V \
                                   + rate[4] - rate[5]

                    return [dCA2dt, dCA2_primedt, dCBdt,
                            dCABdt, dCAB_primedt, dCAstardt]

    def find_steady_state_roots(self,
                                theta0=[0.],
                                Concs0=np.zeros(4),
                                X=0.01,
                                tol=1e-22,
                                solver='secant'):
        mp.dps = 25
        mp.pretty = True

        def get_findroot_eqns(*args):
            if self.model=='site':
                return self.get_odes(args, 0, X, self)
            else:
                # reactor based model
                return self.get_odes_CSTR(args, 0, self)

        if not self.model == 'reactor':
            args = tuple(theta0)
            md = False
        else:
            args = tuple(Concs0)
            md = True

        soln = findroot(get_findroot_eqns,
                        args,
                        solver=solver,
                        tol=tol,
                        multidimensional=md)

        self.steady_state_soln = soln
        return soln

    def integrate_odes(self,
                       theta0=0.01,
                       Concs0=np.zeros(4),
                       X=0.01,
                       timespan=[0, 1e8],
                       h0=1e-20,
                       mxstep=200000,          # maximum number of steps
                       rtol=1E-12,            # relative tolerance
                       atol=1E-18,           # Absolute tolerance
                       full_output=1):

        # if not self.plasma_on:
        if self.model == 'site':
            self.theta0 = theta0
            theta, out = odeint(simpleMkm.get_odes,
                                theta0,
                                timespan,
                                args=(X, self,),
                                h0=h0,
                                mxstep=mxstep,  # maximum number of steps
                                rtol=rtol,      # relative tolerance
                                atol=atol,      # Absolute tolerance
                                full_output=full_output)
            self.integrated_thetas = theta
            self.ode_output = out

            return theta[-1, :]

        elif self.model == 'reactor':
            self.Concs0 = Concs0
            args = self,
            Concs, out = odeint(simpleMkm.get_odes_CSTR,
                                Concs0,
                                timespan,
                                args=args,
                                h0=h0,
                                mxstep=mxstep,  # maximum number of steps
                                rtol=rtol,      # relative tolerance
                                atol=atol,      # Absolute tolerance
                                full_output=full_output)
            self.integrated_Concs = Concs
            self.ode_output = out

            return Concs[-1, :]
