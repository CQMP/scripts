import os
import numpy as np
import subprocess


def run(U,L,data_dir,dry_run,nprocs):

    ncycles = 2**13 # number of measurements
    nwarmup = 1 # number of warmup steps
    cycle_len = 2**6 # how many mc steps in a measurement
    random_seed = True # set completely random seed
    save_plaintext = False # save data to plaintext files
    overwrite_dir = False #True # overwrite old files (useful for recalculation)
    chebyshev = False
    cheb_prefactor = 2.0
    calc_history = False # True # this gets errorbars of dos 

    precision = 4
    step = 10**(-precision)*2

    origdir = os.getcwd()

    T_approx = 0.1
    print "Tc approx =", T_approx
    # dirty scan
    #T_array = np.arange(max(round(T_approx*0.7,precision), 0.005), round(T_approx*1.8,precision), 0.005)

    # prec scan
    T_array = np.arange(T_approx - 0.005, T_approx + 0.005, step)

    T_array = T_array[::-1]

    print "U=",U
    print "T=",T_array
    print "nprocs = ",nprocs

    #if dry_run:
    #    exit()

    for T in T_array:
        beta=1.0/T
        print
        print "============="
        print "Making calc in",
        dirname=os.path.abspath(origdir+os.path.sep+"{0}/L{1}/U{2}/T{3}".format(data_dir,L,U,T))
        print dirname
        dir_exists = os.path.exists(dirname)
        make_calc = True
        if not dir_exists:
            os.makedirs(dirname)
            make_calc = True
        else:
            make_calc = not os.path.exists(dirname+os.path.sep+"output.h5") or overwrite_dir
            
        print "Making calc:",make_calc

        if make_calc:
            os.chdir(dirname)
            args = [
            "--U", str(float(U)),
            "--T", str(float(T)),
            "--mu", str(float(U)/2.0),
            "--L", str(L),
            "--seed" if random_seed else "",
            "-p" if save_plaintext else "",
            "--chebyshev" if chebyshev else "",
            "--calc_history", str(int(calc_history)),
            "--ncycles", str(ncycles),
            "--cyclelen", str(cycle_len),
            "--cheb_prefactor", str(cheb_prefactor), 
            "--nwarmup", str(nwarmup),
            "--dos_width", str(max(6.0, 2.5*U))
            ]

            print os.getenv("QPREFIX")
            call_args = [
            "mpirun", "--np", str(nprocs), os.getenv("QPREFIX")+"/bin/fk_mc_cubic2d"] + args
            print call_args
            print ' '.join(call_args)
            if not dry_run:
                subprocess.call(' '.join(call_args),shell=True)
            else:
                print "skipping run"
            os.chdir(origdir)

        #    run_fkmc.submit(args)
        print "\n"
        #    counter = counter + 1
 

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='FK U loop')
    parser.add_argument('--U', '-U', help='value of U', type=float, default = 1.0)
    parser.add_argument('--data_dir', help='data directory', default = "output")
    parser.add_argument('--L', help='L',type=int, default = 12)
    parser.add_argument('--nprocs', help='n procs',type=int, default = 8)
    parser.add_argument('--dry_run', help='dry run', action="store_true", dest="dry_run", default=False)
    args = parser.parse_args()

    run(U = args.U, L = args.L, data_dir = args.data_dir, dry_run = args.dry_run, nprocs = args.nprocs)
