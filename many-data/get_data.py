import os
import numpy as np
import subprocess
import time
import sys
import shutil
import itertools
#from cluster_tools import *
import argparse
import collections
from cluster_tools import *

data_dir_array = ["pd"]
L_array = [ 16,8,12 ] # array of system sizes
#U_array =[2.0] # array of U values
U_array = [ 0.25, 0.5, 1.0, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 16.0, 18.0, 20.0 ] 
#U_array = [ 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0 ] 
#U_array = [10.0]

dry_run = False # if true - will not make actual calculation
queue = "mid10" # queue for execution
python_path = "/opt/local/bin/python2.7"
mklroot = "/opt/intel/composer_xe_2013.3.171/mkl"
nprocs = 16

print "Output_dirs     :",data_dir_array
print "System sizes :",L_array
print "U values     :",U_array

def main():
    import mc_scan 
    counter = 0
    origdir = os.getcwd()
    for data_dir in data_dir_array:
      for L in L_array:
        for U in U_array:
            print "U=",U,"L=",L
            counter = counter + 1
            name = "U"+str(U)+"L"+str(L)

            call_args = [ python_path, 
                 os.getcwd()+os.path.sep+"mc_scan.py",
                 "--U",str(U),
                 "--L",str(L),
                 "--nprocs",str(nprocs),
                 "--data_dir", str(data_dir)
               ]
            #mc_scan.run(U = U, L = L, nprocs = nprocs, data_dir = data_dir, dry_run = False)
            
            submit_mpi(
                nprocs = nprocs,
                args = call_args,
                add_mpirun = False,
                queue=queue,
                name=name,
                dry_run=dry_run,
                output_stream_file = "stdout_"+name,
                error_stream_file = "stderr_"+name,
                cpu_time = "144:00:00",
                file_size = "40G",
                ram_size = "12G",
                #preamb="source "+mklroot+"/bin/mklvars.sh intel64"
                preamb = "source "+os.environ["HOME"]+"/.bashrc"
                )     
            
    print counter, "calcs submitted"

if __name__ == "__main__":
    main()
