# import pbs tools
from pbs_osg import *
account= "TG-DMR130036"

print "Total number of jobs of current user =",njobs()
print "njobs named test1 =",njobs("test1")
# submit serial (1-cpu) job
# dry_run means create a pbs file (default = "task.pbs") but not run qsub
submit_serial(commands = [["/bin/date"], ["/bin/ls"]], name="test1", cpu_time="00:01:00", dry_run = True)
#submit_serial(commands = ["data"], queue = "serial", name="test1", cpu_time="00:01:00", use_scratch=True, dry_run = True)
#submit_serial(commands = ["touch","1.dat"], queue = "serial", name="test1", cpu_time="00:01:00", use_scratch=True, dry_run = True)
#submit_serial(commands = [["date"],["touch","1.dat"],["date > 1.dat"]], queue = "serial", name="test1", cpu_time="00:01:00", use_scratch=True, dry_run = True)

