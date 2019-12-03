#!/bin/bash
#$ -N MkmRunner
#$ -pe smp 2
#$ -q long
#$ -cwd
source ~/.bash_profile
python EA=-6.00e-01--plasma_on=False--rxn2_barrier=True--total_site_density=1.00e-03--nu0=1.67e-01--V=1.00e+00--reverse=False.py 1> EA=-6.00e-01--plasma_on=False--rxn2_barrier=True--total_site_density=1.00e-03--nu0=1.67e-01--V=1.00e+00--reverse=False.out 2> EA=-6.00e-01--plasma_on=False--rxn2_barrier=True--total_site_density=1.00e-03--nu0=1.67e-01--V=1.00e+00--reverse=False.err
