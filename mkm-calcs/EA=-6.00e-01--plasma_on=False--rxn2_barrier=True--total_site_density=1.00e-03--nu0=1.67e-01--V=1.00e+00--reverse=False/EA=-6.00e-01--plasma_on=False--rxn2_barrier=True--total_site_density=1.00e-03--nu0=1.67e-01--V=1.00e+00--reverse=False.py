import numpy as np
from mkmutils import mkmRunner, load_variables

R = 82.057338 # cm3-atm / K / mol

prefix = 'EA=-6.00e-01--plasma_on=False--rxn2_barrier=True--total_site_density=1.00e-03--nu0=1.67e-01--V=1.00e+00--reverse=False'

inp = load_variables('EA=-6.00e-01--plasma_on=False--rxn2_barrier=True--total_site_density=1.00e-03--nu0=1.67e-01--V=1.00e+00--reverse=False.mkminp'.format(prefix))

EA = inp['EA']

runner = mkmRunner(**inp)
runner.run()
