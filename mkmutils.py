import numpy as np
from mpmath import mp
import json

mp.dps = 50
J2eV = 6.24150974E18            # eV/J
h = 6.626068E-34 * J2eV

R = 82.057338  # cm3-atm / K / mol
Na = 6.022140e23  # 1 / mol

# Spectroscopic constants of nitrogen
# we, wexe, weye, weze cite:heijkers-2015-co2con-microw
# ws = np.array([2372.45, 18.1017, 1.27552e-2, -7.95949])  # in cm-1
# cm1eV = 0.00012398426
# E_Nvibs = cm1eV * ws


def make_PES(Hs,
             Eas,
             ntimes,
             names,
             col='C0',
             width=1,
             fontsize=22,
             label=None,
             axis_labels=True,
             Eref=0.,
             ls='-',
             IS_start=0.,
             legend=False):
    import matplotlib.pyplot as plt
    # Initialize
    x = width
    E = Eref
    line, = plt.plot([IS_start, x], [E, E], c=col, ls=ls)
    # col = line.get_color()
    for i, nEaHt in enumerate(zip(names, Eas, Hs, ntimes)):
        name, Ea, H, times = nEaHt

        for time in range(1, times+1):

            if Ea != 0:
                # barrier
                xs = [x, x + 0.5 * width, x + width]
                es = [E, E + Ea, E+H]
                # Fitting polynomial
                p = np.polyfit(xs, es, 2)

                xfit = np.linspace(x, x + width, 10000)
                efit = np.polyval(p, xfit)

                plt.plot(xfit, efit, c=col, ls=ls)
                plt.text(x + 0.8 * width,
                         E + Ea - 0.1, name, fontsize=fontsize - 4)

            else:
                # No barrier
                plt.plot([x, x + width], [E, E+H], c=col, ls=ls)
                plt.text(x, E + 0.1, name, fontsize=fontsize - 4)
            x += width
            E += H
            plt.plot([x, x + width], [E, E], c=col, ls=ls)
            x += width

    if axis_labels:
        plt.ylabel('Potential Energy (eV)', fontsize=fontsize + 2)
        plt.xlabel('Reaction Coordinate', fontsize=fontsize + 2)

    if label:
        line, = plt.plot([], [], ls=ls, label=label)
        if legend:
            plt.legend(fontsize=fontsize)
        return line
    plt.xticks([])
    plt.yticks(fontsize=fontsize)


def plot_conversion(Ts, Xs,
                    Xeqs, label='',
                    legend=True, prefix='',
                    upper_ylim=0.2, Agg=True):

    if Agg:
        import matplotlib
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.style.use('seaborn-talk')
    plt.figure()

    plt.plot(Ts, Xs, 'o--', label=label)
    plt.plot(Ts, Xeqs, 'k--')

    plt.xlabel('Temperature (K)')
    plt.ylabel('A$_{2}$ Conversion (%)')
    plt.title(prefix)
    plt.ylim(0., upper_ylim)
    if legend:
        plt.legend()
    plt.savefig('{0}-X-v-T.png'.format(prefix), dpi=200)


def load_variables(filename):
    with open(filename) as f:
        d = json.load(f)
    return d


qscript_template = '''#!/bin/bash
#$ -N {1}
#$ -pe {2} {3}
#$ -q {4}
#$ -cwd
source ~/.bash_profile
python {0}.py 1> {0}.out 2> {0}.err
'''

py_template = '''import numpy as np
from mkmutils import mkmRunner, load_variables

R = 82.057338 # cm3-atm / K / mol

prefix = '{0}'

inp = load_variables('{0}.mkminp'.format(prefix))

EA = inp['EA']

runner = mkmRunner(**inp)
runner.run()
'''


class mkmRunner:
    """
    Utility class for running microkinetic model calculations
    using simpleMkm
    """
    def __init__(self,
                 EA,
                 **kwargs):

        # defaults
        self.allowed_keys = dict(Tstart=300.,
                                 Tstop=1000.,
                                 npts=30,
                                 E_rxn=-0.9514,
                                 S=[0.001986,
                                    0.001354,
                                    0.001998],
                                 stoich=[-1, -3, 2],                                 
                                 total_site_density=1.e-3,  # mols / cm3
                                 plasma_on=False,
                                 nu0=10./60,  # Initial volumetric flow cm3/ s
                                 V=1.0,  # reactor volume cm3
                                 k_eimpact_A2=1e-8,  # cm3 / s
                                 k_eimpact_AB=0.0,  # cm3 / s
                                 n_e=1e12,  # electron number density, 1/cm3
                                 rxn2_barrier=False,
                                 concentration_based=True,
                                 model='reactor',
                                 excite_type='diss',
                                 Evib_A2=0.3,
                                 Evib_AB=0.3,
                                 p0=1.0,
                                 A2_frac=1/4.,
                                 reverse=False,
                                 report_conversion=True,
                                 include_plasma_rxns=False)

        input_params = {}
        input_params['EA'] = EA
        for key, val in self.allowed_keys.iteritems():
            if key in kwargs:
                input_params[key] = kwargs[key]
            else:
                input_params[key] = val

        self.input_params = input_params
        prefix = self.get_prefix()

        self.prefix = prefix
        self.input_file = '{0}.mkminp'.format(prefix)
        self.output_file = '{0}.mkmout'.format(prefix)
        self.py_file = '{0}.py'.format(prefix)
        self.qscript = '{0}.sh'.format(prefix)

    def get_prefix(self):

        inp = self.input_params
        prefix = ''
        keys2include = ['EA',
                        'plasma_on',
                        'rxn2_barrier',
                        'total_site_density',
                        'nu0',
                        'V',
                        'reverse']

        if inp['plasma_on']:
            keys2include.extend(['excite_type',
                                 'n_e',
                                 'k_eimpact_A2',
                                 'k_eimpact_AB'])

            if inp['excite_type'] == 'vib':
                keys2include.extend(['Evib_A2',
                                    'Evib_AB'])

        for key in keys2include:
            val = inp[key]
            if isinstance(val, float):
                valstr = '{0:1.2e}'.format(val)
            else:
                valstr = str(val)
            prefix += '{0}={1}--'.format(key, valstr)

        return prefix.strip('--')

    def write_json(self, file_type='inp'):

        if file_type == 'inp':
            fname = self.input_file
            params = self.input_params

        elif file_type == 'out':
            fname = self.output_file
            params = self.output

        with open(fname, 'w') as f:
            json.dump(params, f)

    def write_input(self,
                    jobname='MkmRunner',
                    pe='smp',
                    ncores=2,
                    q='long'):

        self.write_json()
        py = py_template.format(self.prefix)

        with open(self.py_file, 'w') as f:
            f.write(py)

        with open(self.qscript, 'w') as f:
            f.write(qscript_template.format(self.prefix,
                                            jobname,
                                            pe,
                                            ncores,
                                            q))

    def run_job(self):
        """
        Submit job to queue
        """
        from subprocess import Popen, PIPE
        qscript = self.qscript
        p = Popen(['qsub', qscript], stdin=PIPE, stdout=PIPE, stderr=PIPE)

        out, err = p.communicate()
        jobid = None
        if out == '' or err != '':
            raise Exception('something went wrong in qsub:\n\n{0}'.format(err))

        jobid = out.split()[2]

        with open('{0}.jobid'.format(self.prefix), 'w') as f:
            f.write(jobid)

        return jobid

    def run(self):
        """
        Run temperature sweep model
        """
        
        from simplemkm import simpleMkm as mkm

        inp = self.input_params
        EA = inp['EA']
        Ts = np.linspace(inp['Tstop'], inp['Tstart'], inp['npts'])

        plasma_on = inp['plasma_on']

        allConcs = []
        allpressures = []
        eqpressures = []
        Xeqs = []
        Xs = []
        all_net_rates = []
        all_forward_rates = []
        all_reverse_rates = []        

        if inp['reverse']:
            Xeq = 0.8
        else:
            Xeq = 0.1

        for T in Ts:

            mod = mkm(T,
                      **inp)

            Xeq = mod.get_eqb_conversion(X0=Xeq)
            Xeqs.append(Xeq[0])
            eqp = mod.get_pressures(Xeq[0])
            eqpressures.append([mp.nstr(p, 25)
                                for p in eqp])

            if T == Ts[0]:
                if not inp['reverse']:
                    CA2 = inp['A2_frac'] * inp['p0'] / R / T
                    CB = (1. - inp['A2_frac']) * inp['p0'] / R / T
                    CAB = 0.0

                else:
                    CAB = inp['p0'] / R / T
                    CA2 = 0.0
                    CB = 0.0

                if plasma_on:
                    Concs = np.array([CA2, 0.0, CB, CAB, 0.0, 0.0])
                else:
                    Concs = np.array([CA2, CB, CAB, 0.0])

                Concs = mod.integrate_odes(Concs0=Concs)
            else:
                Concs[0:-1] = pressures / R / T

            try:
                Concs = mod.find_steady_state_roots(Concs0=Concs)
            except Exception as E:
                print "{0} K. Root finder failed once: {1}".format(T,
                                                                   E)
                Concs = mod.integrate_odes(Concs0=Concs)
                try:
                    Concs = mod.find_steady_state_roots(Concs0=Concs)
                except Exception as E:
                    print "{0} K. Root finder failed again: {1}".format(T,
                                                                        E)

            allConcs.append([mp.nstr(C, 25) for C in Concs])
            pressures = Concs[0:-1] * R * T
            allpressures.append([mp.nstr(p, 25) for p in pressures])

            net_rates = mod.get_rates_conc(Concs)
            forward_rates = mod.get_rates_conc(Concs,
                                               rate_type='forward')
            reverse_rates = mod.get_rates_conc(Concs,
                                               rate_type='reverse')

            all_net_rates.append([mp.nstr(r) for r in net_rates])
            all_forward_rates.append([mp.nstr(r) for r in forward_rates])
            all_reverse_rates.append([mp.nstr(r) for r in reverse_rates])

            if inp['report_conversion']:
                X = [0.1, 0.1, 0.1]

                if not plasma_on:
                    X = mod.get_conversion(pressures[0], mod._pA2_0)
                    Xs.append(mp.nstr(X, 25))

                else:
                    X = mod.get_conversion(pressures[0],
                                           mod._pA2_0,
                                           pressures[1],
                                           pressures[4],
                                           X0=X)

                    Xs.append([mp.nstr(x, 25) for x in X])

        self.output = {}
        self.output['EA'] = EA
        self.output['T'] = Ts.tolist()
        self.output['Concs'] = allConcs
        self.output['pressures'] = allpressures
        self.output['Xeq'] = Xeqs
        self.output['eq_pressures'] = eqpressures
        self.output['X'] = Xs
        self.output['forward_rates'] = all_forward_rates
        self.output['reverse_rates'] = all_reverse_rates
        self.output['net_rates'] = all_net_rates

        self.write_json('out')
