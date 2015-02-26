import os
import numpy as np
import subprocess
import time
import sys
import shutil
import itertools
import argparse
import collections
from pbs_stampede import *

data_dir_array = ["data"]

Gamma_array = [0.3]
# Kondo regime
U_array = [2.0]
ed_array = [-0.1]

#T_array = [0.1, 0.13, 0.15, 0.2, 0.5, 1.0]
#T_array = [0.1, 0.13, 0.15, 0.2, 0.5, 1.0, 0.05]
T_array = [0.1, 0.07, 0.15, 0.2, 0.3, 0.5, 1.0]
T_array = [0.1, 0.07, 0.2]
#T_array = [0.3, 0.5, 1.0]
#T_array = [0.1]
#T_array=[0.05]
bw_array = [5.0]
tmax_array=[5.0]#, 0.1, 0.2, 0.3, 0.4, 0.5] #np.arange(0,0.5,0.1)
voltage_array = np.arange(0,0.21, 0.05)
print "voltage array : ",voltage_array 

run_script="bold_hyb_run.py"

max_nca_order = 100

run_qmc = False
run_nca = True
clean = False
dry_run = False # True # False # True # False

cycle_len = 64

nprocs_min = 8
nprocs_max = 8

nca_mixing = 1.0

def main():
    counter = 0
    origdir = os.getcwd()
    for data_dir in data_dir_array:
      for Gamma in Gamma_array:
        Gamma = float(Gamma)
        for T in T_array:
            T=float(T)
            beta = 1.0/T
            for bw in bw_array:
                bw = float(bw)
                for U in U_array:
                    U=float(U)
                    for ed in ed_array:
                        ed=float(ed)

                        tk = Tk(U, ed, Gamma)
                        print "Kondo temperature =", tk
                        tk_round = round(tk, 2)
                        print tk_round

                        #h_step = round(tk_round / 2., 3)
                        #h_array = np.arange(0,3*tk_round,h_step)
                        for voltage in voltage_array:
                            for tmax in tmax_array:
                                tmax = float(tmax)
                                print "U=",U,"tmax=",tmax,"beta = ",beta, "T=",1/beta, "bw = ",bw, "Gamma =", Gamma

                                levels = imp_levels(U, ed, 0)
                                print "levels =", levels 

                                print "Kondo temperature =", tk

                                dirname = origdir+"/"+data_dir+"/G"+str(Gamma)+"/T"+str(T)+"/bw"+str(bw)+"/U"+str(U)+"/ed"+str(ed)+"/v"+str(voltage)+"/tmax"+str(tmax)
                                print dirname
                                os.makedirs(dirname) if not os.path.exists(dirname) else None
                                os.chdir(dirname)

                                np.savetxt("levels.dat", np.array(levels).transpose())
                                np.savetxt("Tk.dat", np.array([tk]))

                                nprocs = nprocs_min
                                ncycles = ncycles_f(tmax, bw, beta = beta, Gamma = Gamma, cycle_len = cycle_len, nprocs = nprocs) 
                                while ncycles>1e9 and nprocs<nprocs_max:
                                    nprocs*=2
                                    ncycles = ncycles_f(tmax, bw, beta = beta, Gamma = Gamma, cycle_len = cycle_len, nprocs = nprocs)
                                print "Doing",ncycles," cycles on",nprocs,"cpus" 

                                name = "b"+str(round(beta,2))+"t"+str(round(tmax,2))+"v"+str(round(voltage,2))

                                npts_real = max(round(tmax/0.02), 250)
                                npts_imag = max(round(beta/0.02), 250)

                                args = [which("python"), origdir+os.path.sep+run_script,
                                        "--T", str(float(T)),
                                        "--gamma", str(float(Gamma)), 
                                        "--tmax", str(float(tmax)),
                                        "--run_qmc", str(int(run_qmc)),
                                        "--run_nca", str(int(run_nca)),
                                        "--npts_real", str(int(npts_real)),
                                        "--npts_imag", str(int(npts_imag)),
                                        "--max_nca_order", str(int(max_nca_order)), 
                                        "--levels",  " ".join([str(x) for x in levels]),
                                        "--half_bandwidth"  , str(float(bw)),
                                        "--mode"  , "det",
                                        "--mixing"  , str(float(nca_mixing)),
                                        "--ncycles", str(ncycles),
                                        "--voltage", str(voltage), 
                                        "--cycle_len", str(cycle_len),
                                        "--clean", str(int(clean))
                                        ]
                                print ' '.join(args)
                                if 1 == 1: # not os.path.exists("output.h5"):
                                    while njobs()>48:
                                        print "too many jobs in queue - pending"
                                        time.sleep(300)
                                    submit_mpi(
                                        commands = args, 
                                        prefix="module load bold_hyb",
                                        add_mpirun = False, 
                                        nprocs = nprocs, 
                                        queue = "normal", 
                                        name=name, 
                                        #cpu_time="24:00:00", 
                                        cpu_time="00:30:00", 
                                        ram_size = "31000mb",
                                        file_size = "10000mb",
                                        use_scratch=False, 
                                        pbs_file = "taskmpi.pbs", 
                                        dry_run = dry_run
                                    )
                                counter = counter + 1

                                os.chdir(origdir)
                
    print counter, "calcs submitted"
    if not dry_run:
        send_final_message("aantipov@gmail.com", txt = "Submitted " + str(counter) + " calcs in " + origdir)

def imp_levels(U,ed,h):
    return np.array([0, ed - h/2., ed + h/2, 2*ed + U])
def Tk(U,ed,Gamma):
    tl = np.sqrt(Gamma * U / 2) * np.exp( -np.pi * np.abs(ed) * np.abs(ed + U) / (2.*Gamma*U ))
    tk = 0.4 * tl
    return tk
 
def nsteps_f(t, bw, beta, Gamma):
    t0 = 0.4
    beta0 = 1.0
    bw0 = 10 
    bw = bw * Gamma
    ncycles_t0 = 2**25 * 16 * 16 # empirical, the calc was for ncycles=2**27 x nprocs=16 x cycle_len=16
    return int(round(max( (beta / beta0)**2 * ncycles_t0*np.exp(2.2*(bw*t-bw0*t0)), ncycles_t0))) 
def ncycles_f(t, bw, beta, Gamma, cycle_len, nprocs):
    return int(round(nsteps_f(t, bw, beta, Gamma) / cycle_len / nprocs))
def send_final_message(email, txt = "I've finished"):
    import smtplib
    from email.MIMEMultipart import MIMEMultipart
    from email.MIMEText import MIMEText
    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.ehlo()
    server.starttls()
    server.ehlo()
    server.login("andrey.e.antipov@gmail.com", "groyxrvltxgzcyvs")
    msg = MIMEMultipart()
    fromaddr = "andrey.e.antipov@gmail.com"
    toaddr = "aantipov@umich.edu"
    msg['From'] = fromaddr 
    msg['To'] = toaddr
    msg['Subject'] = "Stampede Herald"
    msg.attach(MIMEText(txt, 'plain'))
    server.sendmail(fromaddr, toaddr, msg.as_string())

if __name__ == "__main__":
    main()
