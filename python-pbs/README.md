#### python-pbs
A python module that helps with submitting jobs to different supercomputers. The detailed wrapper is cluster-dependent, so different modules for different clusters are collected here. The main two commands are `submit_serial(command, job name, account to use, ...) ` and `submit_mpi(command, nprocs, job name, ...)`.

See `pbs_example.py` for a usage example.

- `pbs_flux.py` - umich Flux wrapper (PBS).
- `pbs_stampede.py` - TACC Stampede wrapper (SLURM).

##### Usage 
The best way to use it is to create a directory somewhere in home location, say `${HOME}/pythonmodules` and add the line 
`export PYTHONPATH=${HOME}/pythonmodules:$PYTHONPATH` to `.bash_profile`. Afterwards any python script will find this module automatically.
