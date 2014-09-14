#### python-pbs
A python module that helps with submitting jobs to different supercomputers. The detailed wrapper is cluster-dependent, so different modules for different clusters are collected here. The main two commands are `submit_serial(command, job name, account to use, ...) ` and `submit_mpi(command, nprocs, job name, ...)`.

See `pbs_example.py` for a usage example.

- `pbs_flux.py` - flux wrapper
