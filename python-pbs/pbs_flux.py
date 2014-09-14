import subprocess
import os
import random
import string
import time

""" Python wrapper around pbs """

__all__ = ["which", "njobs", "submit_serial", "submit_mpi"]

def which(file):
    """ Equivalent of bash which command """
    for path in os.environ["PATH"].split(":"):
        if os.path.exists(path + "/" + file):
                return path + "/" + file

    return None

def njobs(name = ""):
    """ Wrapper around qstat - gets the total number of tasks with a given name or account """ 
    user = os.environ.get('USER')
    ntasks = int(subprocess.check_output("qstat | grep "+user+" | grep '"+name+"' | wc -l", shell=True))
    return ntasks

# Submit mpi job
def submit_mpi(   commands,  # commands to submit
                  name, # name of the job
                  nprocs, # number of cpus
                  account, # which account to use
                  queue = "flux", # queue name
                  prefix = "", # what will be executed before the command. Smth like "export LD_LIBRARY_PATH="..."
                  postfix = "", # what will be executed after the run
                  add_mpirun = False, # add mpirun in front
                  mpi_exec = "mpirun", # which mpirun command
                  use_scratch = False, # work in scratch
                  scratch_location = os.environ["SCRATCH"], # location of scratch fs 
                  job_input = "*.ini *.dat *.h5", # what to copy to scratch as input
                  output_stream_file = "stdout", # where to gather output 
                  error_stream_file = "stderr", # where to gather errors
                  cpu_time = "164:00:00", # how long to run
                  ram_size = "4096mb", # memory allocation
                  file_size="4096mb",
                  pbs_file = "task.pbs", # name of the pbs file 
                  dry_run = False ):
    """ Submit mpi job """
    print "Submitting",' '.join(commands),"to account",account
    start_dir = os.getcwd()
    scratch_dir = scratch_location+os.path.sep
    scratch_dir = scratch_dir + time.strftime("%Y-%m-%d:%H:%M:%S", time.gmtime())
    scratch_dir = scratch_dir +'_' +''.join(random.choice(string.lowercase) for x in range(4))
    print "scratch dir: ", scratch_dir
    bash_path = which("bash")
    if os.path.exists(output_stream_file):
        print "Removing",output_stream_file
        os.remove(output_stream_file)
    if os.path.exists(error_stream_file):
        print "Removing",error_stream_file
        os.remove(error_stream_file)
    f = open(pbs_file,"w")
    f.write("#!"+bash_path+os.linesep)
    f.write("#PBS -S "+bash_path+os.linesep) # load bash
    f.write("#PBS -N "+name+os.linesep) # job name
    f.write("#PBS -A "+account+os.linesep) # account
    f.write("#PBS -q "+queue+os.linesep) # queue
    f.write("#PBS -l qos=flux"+os.linesep)
    f.write("#PBS -l procs="+str(nprocs)+",walltime="+cpu_time+os.linesep)
    f.write("#PBS -l pmem="+ram_size+os.linesep)
    f.write("#PBS -o "+os.path.join(start_dir, output_stream_file)+os.linesep) #output
    f.write("#PBS -e "+os.path.join(start_dir, error_stream_file)+os.linesep) #errors
    f.write("#PBS -V "+os.linesep) # Takes your current environment (like PATH and other variables) and sends them along with your job to the compute node.
    f.write(os.linesep)
    f.write("cd "+start_dir+os.linesep)
    if use_scratch:
        f.write("mkdir -p "+scratch_dir+os.linesep)
        f.write("cp -v "+job_input+" "+scratch_dir+os.path.sep+os.linesep)
        f.write("cd "+scratch_dir+os.linesep)

    f.write("cat $PBS_NODEFILE" + os.linesep)
    f.write(prefix + os.linesep)

    mpi_exec = which(mpi_exec)
    call_str = [mpi_exec + " --np " + str(nprocs)+" " if add_mpirun else ""][0] 
    call_str = call_str + ' '.join(commands)

    f.write(call_str+os.linesep+os.linesep)
    f.write(postfix + os.linesep)
    if use_scratch:
        f.write("rsync -a "+scratch_dir+os.path.sep+"* "+start_dir+os.path.sep+os.linesep)
        f.write("rm -rfv "+scratch_dir+os.linesep)
    f.close()
    
    qsub = which("qsub")
    full_command=qsub+" "+pbs_file

    print full_command
    if not dry_run:
        print subprocess.check_output(full_command, shell=True)

# Submit serial job
def submit_serial(commands,  # commands to submit
                  name, # name of the job
                  account, # which account to use
                  queue = "flux", # queue name
                  prefix = "", # what will be executed before the command. Smth like "export LD_LIBRARY_PATH="..."
                  postfix = "", # what will be executed after the run
                  use_scratch = False, # work in scratch
                  scratch_location = os.environ["SCRATCH"], # location of scratch fs 
                  job_input = "*.ini *.dat *.h5", # what to copy to scratch as input
                  output_stream_file = "stdout", # where to gather output 
                  error_stream_file = "stderr", # where to gather errors
                  cpu_time = "164:00:00", # how long to run
                  ram_size = "4096mb", # memory allocation
                  file_size="4096mb",
                  pbs_file = "task.pbs", # name of the pbs file 
                  dry_run = False ):
 
    """ Submit serial job """
    submit_mpi(commands = commands, 
               name = name, 
               nprocs = 1, 
               account = account, 
               queue = queue, 
               add_mpirun = False,
               prefix = prefix, 
               postfix = postfix, 
               output_stream_file = output_stream_file, 
               error_stream_file = error_stream_file, 
               cpu_time = cpu_time, 
               ram_size = ram_size, 
               file_size = file_size, 
               use_scratch = use_scratch, 
               scratch_location = scratch_location,
               pbs_file = pbs_file,
               dry_run = dry_run)



