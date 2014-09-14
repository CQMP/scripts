# import pbs tools
from pbs_flux import *

print "Total number of jobs of current user =",njobs()
print "njobs named test1 =",njobs("test1")
# submit serial job
# dry_run means create a pbs file (default = "task.pbs") but not run qsub
submit_serial(args = ["date\ntouch 1.dat\ndate > 1.dat"], account = "lsa_flux", name="test1", cpu_time="00:01:00", use_scratch=True, dry_run = True)
# submit mpi job
submit_mpi(args = ["date\ntouch 1.dat\ndate > 1.dat"], add_mpirun = True, nprocs = 24, account = "lsa_flux", name="test1", cpu_time="00:01:00", use_scratch=True, pbs_file = "taskmpi.pbs", dry_run = True)
