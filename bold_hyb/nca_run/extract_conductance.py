import os
import numpy as np
import scipy
import scipy.optimize
import string
import itertools
import h5py

def main():
    data_dir="nca_stats"
    os.chdir(data_dir)

    datafiles = filter(
        lambda x: x.find("current_")!=-1 and x.find(".dat")!=-1, 
        os.listdir(os.getcwd()))
    print datafiles
    tmax_array = np.unique(sorted(filter(lambda x: x>0, map(lambda x: float(x.lstrip("current_").rstrip('.dat').split("t")[1]), datafiles))))
    print "tmax :", tmax_array

    for tmax in tmax_array:
        print "tmax =", tmax
        dfiles2 = filter(lambda x : x.find("t"+str(tmax))!=-1, datafiles)
        print dfiles2
        T_array = np.unique(sorted(filter(lambda x: x>0, map(lambda x: float(x.lstrip("current_").split('t'+str(tmax)+'.dat')[0].split("T")[1]), dfiles2))))
        print "T :", T_array
        streams = dict()
        streams["cond"] = open("cond_"+"t"+str(tmax)+".dat","w")

        for T in T_array:
            dataf = filter(lambda x : x.find("T"+str(T)+"t")!=-1, dfiles2)[0]
            data = np.loadtxt(dataf, unpack=True)
            conductance_vals = np.ediff1d(data[1])/np.ediff1d(data[0])
            cond = conductance_vals[0]
            print "T =",T,"conductance = ", cond
            streams["cond"].write(str(T)+"   "+str(cond)+ "\n")
            



if __name__ == "__main__":
    main()
