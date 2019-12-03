#!/bin/bash
#$ -N MkmRunner
#$ -pe smp 2
#$ -q long
#$ -cwd
source ~/.bash_profile
python EA=-6.00e-01--plasma_on=True--rxn2_barrier=True--total_site_density=1.00e-03--nu0=1.67e-01--V=1.00e+00--reverse=False--excite_type=vib--n_e=1.00e+12--k_eimpact_A2=1.00e-14--k_eimpact_AB=0.00e+00--Evib_A2=1.00e+00--Evib_AB=3.00e-01.py 1> EA=-6.00e-01--plasma_on=True--rxn2_barrier=True--total_site_density=1.00e-03--nu0=1.67e-01--V=1.00e+00--reverse=False--excite_type=vib--n_e=1.00e+12--k_eimpact_A2=1.00e-14--k_eimpact_AB=0.00e+00--Evib_A2=1.00e+00--Evib_AB=3.00e-01.out 2> EA=-6.00e-01--plasma_on=True--rxn2_barrier=True--total_site_density=1.00e-03--nu0=1.67e-01--V=1.00e+00--reverse=False--excite_type=vib--n_e=1.00e+12--k_eimpact_A2=1.00e-14--k_eimpact_AB=0.00e+00--Evib_A2=1.00e+00--Evib_AB=3.00e-01.err
