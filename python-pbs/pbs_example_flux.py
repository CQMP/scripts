# import pbs tools
from pbs_flux import *

print "Total number of jobs of current user =",njobs()
print "njobs named test1 =",njobs("test1")
# submit serial (1-cpu) job
# dry_run means create a pbs file (default = "task.pbs") but not run qsub
submit_serial(commands = "data", queue = "flux", account = "egull_flux", name="test1", cpu_time="00:01:00", use_scratch=True, dry_run = True)
submit_serial(commands = ["data"], queue = "flux", account = "egull_flux", name="test1", cpu_time="00:01:00", use_scratch=True, dry_run = True)
submit_serial(commands = ["touch","1.dat"], queue = "flux", account = "egull_flux", name="test1", cpu_time="00:01:00", use_scratch=True, dry_run = True)
submit_serial(commands = [["date"],["touch","1.dat"],["date > 1.dat"]], queue = "flux", account = "egull_flux", name="test1", cpu_time="00:01:00", use_scratch=True, dry_run = True)

# this will raise an error:
# submit_serial(commands = [[["data"]]], queue = "flux", account = "egull_flux", name="test1", cpu_time="00:01:00", use_scratch=True, dry_run = True)

# submit mpi jobs
submit_mpi(commands = [["date"],["touch","1.dat"],["date > 1.dat"]], add_mpirun = True, nprocs = 24, queue = "flux", account = "egull_flux", name="test1", cpu_time="00:01:00", use_scratch=True, pbs_file = "taskmpi.pbs", dry_run = True)
submit_mpi(commands = [["date"],["touch","1.dat"],["date > 1.dat"]], add_mpirun = True, add_mpirun_each = True, nprocs = 24, priority = 1, queue = "flux", account = "egull_flux", name="test1", cpu_time="00:01:00", use_scratch=True, pbs_file = "taskmpi2.pbs", dry_run = True)
