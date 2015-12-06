import subprocess
import shutil
import os

__all__ = ["which", "get_ntasks", "get_total_ntasks", "submit_serial", "submit_mpi"]

def which(file):
    ''' Equivalent of bash which command '''
    for path in os.environ["PATH"].split(":"):
        if os.path.exists(path + "/" + file):
                return path + "/" + file

    return None

def get_ntasks(name):
    ''' Wrapper around qstat - get's total number of tasks with a given name or queue '''
    user = os.environ.get('USER')
    ntasks = int(subprocess.check_output("qstat | grep "+user+" | grep "+name+" | wc -l", shell=True))
    return ntasks

def get_total_ntasks():
    ''' Wrapper around qstat - get's total number of tasks of current user '''
    user = os.environ.get('USER')
    ntasks = int(subprocess.check_output("qstat | grep "+user+" | wc -l", shell=True))
    return ntasks


def submit_serial(args, queue, name, 
                  output_stream_file = "stdout", 
                  error_stream_file = "stderr", 
                  cpu_time = "164:00:00", 
                  ram_size = "4096M",
                  file_size="4096M",
                  dry_run = False ):
    ''' Submits a serial job to the cluster '''
    print "Submitting",' '.join(args),"to queue",queue
    start_dir = os.getcwd()
    bash_path = which("bash")
    if os.path.exists(output_stream_file):
        print "Removing",output_stream_file
        os.remove(output_stream_file)
    if os.path.exists(error_stream_file):
        print "Removing",error_stream_file
        os.remove(error_stream_file)
    f = open("task.pbs","w")
    f.write("#PBS -S "+bash_path+os.linesep)
    f.write("#PBS -M Andrey.E.Antipov@gmail.com"+os.linesep)
    f.write("#PBS -o "+os.path.join(start_dir, output_stream_file)+os.linesep)
    f.write("#PBS -e "+os.path.join(start_dir, error_stream_file)+os.linesep)
    f.write("#PBS -l nodes=1"+os.linesep)
    f.write("#PBS -N "+name+os.linesep)
    f.write(os.linesep)
    f.write("export LD_LIBRARY_PATH=${HOME}/_local_"+queue+"/lib:${LD_LIBRARY_PATH}"+os.linesep) # user-specific
    f.write("cd "+start_dir+os.linesep)

    call_str = ' '.join(args)
    f.write(call_str+os.linesep)
    f.close()
    
    qsub = which("qsub")
    command=qsub+" task.pbs"
    print command
    if not dry_run:
        print subprocess.check_output(command, shell=True)

def submit_mpi(nprocs, args, name, 
                  queue = "regular",
                  output_stream_file = "stdout", 
                  error_stream_file = "stderr", 
                  add_mpirun = True,
                  preamb = "",
                  #mpi_exec = "mpirun",
                  cpu_time = "164:00:00", 
                  ram_size = "4096M",
                  file_size="4096M",
                  dry_run = False ):
    ''' Submits a serial job to the cluster '''
    print "Submitting",' '.join(args),"to queue",queue
    start_dir = os.getcwd()
    bash_path = which("bash")
    if os.path.exists(output_stream_file):
        print "Removing",output_stream_file
        os.remove(output_stream_file)
    if os.path.exists(error_stream_file):
        print "Removing",error_stream_file
        os.remove(error_stream_file)
    f = open("task.pbs","w")
    f.write("#PBS -N "+name+os.linesep)
    f.write("#PBS -q "+queue+os.linesep)
    f.write("#PBS -V "+os.linesep)
    f.write("#PBS -o "+os.path.join(start_dir, output_stream_file)+os.linesep)
    f.write("#PBS -e "+os.path.join(start_dir, error_stream_file)+os.linesep)
    f.write("#PBS -l mppwidth="+str(nprocs)+os.linesep)
    f.write(os.linesep)
    #f.write("export LD_LIBRARY_PATH=${HOME}/_local_"+queue+"/lib:${LD_LIBRARY_PATH}"+os.linesep) # user-specific
    f.write(preamb+os.linesep)
    f.write("cd "+start_dir+os.linesep)

    call_str = ["aprun -n " + str(nprocs)+" " if add_mpirun else ""][0] 
    call_str = call_str + ' '.join(args)
    print call_str
    f.write(call_str+os.linesep)
    f.close()
    
    qsub = which("qsub")
    command=qsub+" task.pbs"
    print command
    if not dry_run:
        print subprocess.check_output(command, shell=True)
