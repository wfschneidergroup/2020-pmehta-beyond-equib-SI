import numpy as np
from mkmutils import mkmRunner, load_variables

R = 82.057338 # cm3-atm / K / mol

prefix = 'EA=-1.20e+00--plasma_on=True--rxn2_barrier=True--total_site_density=1.00e-03--nu0=1.67e-01--V=1.00e+00--reverse=False--excite_type=vib--n_e=1.00e+12--k_eimpact_A2=1.00e-14--k_eimpact_AB=0.00e+00--Evib_A2=1.00e+00--Evib_AB=3.00e-01'

inp = load_variables('EA=-1.20e+00--plasma_on=True--rxn2_barrier=True--total_site_density=1.00e-03--nu0=1.67e-01--V=1.00e+00--reverse=False--excite_type=vib--n_e=1.00e+12--k_eimpact_A2=1.00e-14--k_eimpact_AB=0.00e+00--Evib_A2=1.00e+00--Evib_AB=3.00e-01.mkminp'.format(prefix))

EA = inp['EA']

runner = mkmRunner(**inp)
runner.run()
