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
    ntasks = int(subprocess.check_output("condor_q -submitter aantipov | wc -l", shell=True))
    return ntasks

# Submit mpi job
def submit_mpi(   commands,  # commands to submit
                  name, # name of the job
                  nprocs, # number of cpus
                  account = "TG-DMR130036",
                  reqs = "(CVMFS_oasis_opensciencegrid_org_REVISION >= 3620)",  
                  prefix = "", # what will be executed before the command. Smth like "export LD_LIBRARY_PATH="..."
                  postfix = "", # what will be executed after the run
                  job_input = "", # what to copy to scratch as input
                  job_output = "", # what to copy to scratch as input
                  output_stream_file = "stdout", # where to gather output 
                  error_stream_file = "stderr", # where to gather errors
                  log_file="log", 
                  cpu_time = "164:00:00", # how long to run
                  ram_size = 1024, # memory allocation MB
                  file_size=512,
                  pbs_file = "task.job", # name of the pbs file 
                  dry_run = False ):
    linesep = '\n' # see https://docs.python.org/2/library/os.html?highlight=linesep#os.linesep
    scratch_dir = os.getcwd() 
    """ Submit mpi job """
    # First analyze the input commands - convert to a 2d array and raise an error if there are more
    commands = [commands] if not isinstance(commands, list) else commands
    commands = [commands] if not isinstance(commands[0], list) else commands
    if isinstance(commands[0][0], list):
        raise NameError('Commands should be a list with depth <=2') 
    # By now we have a 2d array of commands. Convert them to '\n' separated commands
    print "Submitting :", linesep.join([' '.join(x) for x in commands]), "to account",account
    start_dir = os.getcwd()
    bash_path = which("bash")
    if os.path.exists(output_stream_file):
        print "Removing",output_stream_file
        os.remove(output_stream_file)
    if os.path.exists(error_stream_file):
        print "Removing",error_stream_file
        os.remove(error_stream_file)
    f = open(pbs_file,"w")
    f.write("universe = vanilla" + linesep)
    f.write('requirements = (FileSystemDomain != "") && (Arch == "X86_64") && ' + reqs + linesep)
    f.write('+ProjectName = ' + account + linesep)
    f.write("output = " + output_stream_file + linesep) #output
    f.write("error = " +  error_stream_file  + linesep) #errors
    f.write("log = "+os.path.join(start_dir, log_file)+linesep) #errors
    f.write("notification = NEVER"+linesep)

    f.write("should_transfer_files = YES" + linesep)
    #f.write("when_to_transfer_output = ON_EXIT" + linesep)
    f.write("when_to_transfer_output = ON_EXIT_OR_EVICT" + linesep)
    f.write("transfer_input_files = " + job_input + linesep)
    f.write("transfer_output_files = " + job_output + linesep)

    f.write("request_memory = " + str(ram_size).rstrip("mb") + linesep)
    f.write("request_disk = " + str(file_size).rstrip("mb") + linesep)
    f.write("request_cpus = " + str(nprocs) + linesep)

    #f.write("#SBATCH -t "+cpu_time+linesep)
    f.write(linesep)
    #f.write(prefix + linesep)

    for x in commands:
        cmd = x[0]
        args = " ".join(x[1:len(x)-1])
        f.write("executable = " + cmd + linesep)
        f.write("arguments = " + args + linesep)
        f.write("queue" + linesep + linesep)
    #call_str = linesep.join([' '.join(x) for x in commands])
    #f.write(call_str+linesep+linesep)
    #f.write(postfix + linesep)
    f.close()
    
    sbatch = which("condor_submit")
    full_command=sbatch+" "+pbs_file

    print full_command
    if not dry_run:
        print subprocess.check_output(full_command, shell=True)

# Submit serial job
def submit_serial(commands,  # commands to submit
                  name, # name of the job
                  account = "TG-DMR130036",
                  reqs = "(CVMFS_oasis_opensciencegrid_org_REVISION >= 3620)",  
                  prefix = "", # what will be executed before the command. Smth like "export LD_LIBRARY_PATH="..."
                  postfix = "", # what will be executed after the run
                  job_input = "", # what to copy to scratch as input
                  job_output = "", # what to copy to scratch as input
                  output_stream_file = "stdout", # where to gather output 
                  error_stream_file = "stderr", # where to gather errors
                  log_file="log", 
                  cpu_time = "164:00:00", # how long to run
                  ram_size = 1024, # memory allocation MB
                  file_size=512,
                  pbs_file = "task.job", # name of the pbs file 
                  dry_run = False ):
    """ Submit serial job """
    submit_mpi(commands = commands, 
               name = name, 
               nprocs = 1, 
               account = account, 
               reqs = reqs,
               prefix = prefix, 
               postfix = postfix, 
               job_input = job_input,
               job_output = job_output,
               output_stream_file = output_stream_file, 
               error_stream_file = error_stream_file, 
               log_file = log_file,
               cpu_time = cpu_time, 
               ram_size = ram_size, 
               file_size = file_size, 
               pbs_file = pbs_file,
               dry_run = dry_run)



