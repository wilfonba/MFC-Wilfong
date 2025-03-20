This case file can be ran from the MFC/ directory via:

./mfc.sh run examples/3D_IGR_perf_test/case.py -N <N> -n <N> --case-optimization -j 32 --gpu <other options> -- -dim <dim> -nt <nt> -ns <ns>

where:

<dim> - inlet direction (1, 2, or 3, default 1)
<nt> - number of time steps (default 100)
<ns> - number of saves (default 1)

This case has a hardcoded size of 400 x 400 x 400 = 64 million cells. This can
be changed by updating line 33 in the case file
