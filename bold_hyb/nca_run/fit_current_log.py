import os
import numpy as np
import scipy
import scipy.optimize
from scipy.optimize import curve_fit
import string
import itertools
import h5py
import collections

def main():

    obs_binned = ["dm0", "dm1", "dm2", "dm3", "n0", "n1", "current"] 
    obs_mcdata = ["order", "sign"]
    obs_composite = ["m", "n"]

    data_dir = "data"
    top_string_array=["nca"]

    origdir=os.getcwd()
    h5filename = "output.h5"

    current_threshold = 5e-4

    for top_string in top_string_array:
        out_dir = top_string + "_stats"
        out_dir = os.path.abspath(out_dir)
        print "=================\n",data_dir, "data_dir\n","=================\n"
        print "Saving output to", out_dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        os.chdir(data_dir)

        def tree():
            return collections.defaultdict(tree)

        fit_vals = tree()
        currents = tree()

        Gamma_array = sorted(map(lambda x: eval(x.strip("G")), filter(lambda x: x.find("G")!=-1,os.listdir(os.getcwd()))),key = lambda x: float(x))
        for Gamma in Gamma_array:
          os.chdir('G'+str(Gamma))
          T_array = sorted(map(lambda x: eval(x.strip("T")), filter(lambda x: x.find("T")!=-1,os.listdir(os.getcwd()))),key = lambda x: float(x))
          for T in T_array:
            os.chdir('T'+str(T))
            beta = 1.0/T
            bw_array = sorted(map(lambda x: eval(x.strip("bw")), filter(lambda x: x.find("bw")!=-1,os.listdir(os.getcwd()))),key = lambda x: float(x))
            for bw in bw_array:
                os.chdir('bw'+str(bw))
                U_array = sorted(map(lambda x: eval(x.strip("U")), filter(lambda x: x.find("U")!=-1,os.listdir(os.getcwd()))),key = lambda x: float(x))
                for U in U_array:
                    os.chdir('U'+str(U))
                    ed_array = sorted(map(lambda x: eval(x.strip("ed")), filter(lambda x: x.find("ed")!=-1,os.listdir(os.getcwd()))),key = lambda x: float(x))
                    for ed in ed_array:
                        os.chdir('ed'+str(ed))
                        v_array = sorted(map(lambda x: eval(x.strip("v")), filter(lambda x: x.find("v")!=-1,os.listdir(os.getcwd()))),key = lambda x: float(x))
                        v_array = filter(lambda x : x > 0, v_array)
                        for v in v_array:
                            os.chdir('v'+str(v))
                            tmax_array = [5.0]
                            for tmax in tmax_array:
                                os.chdir('tmax'+str(tmax))

                                print "============="
                                print "Statistics for",
                                print os.path.relpath(os.getcwd(),origdir)


                                if os.path.exists(h5filename):
                                    h5f = h5py.File(h5filename, "r")
                                    top = h5f["nca/stats"]
                                    if top["current_t"]:
                                        current_vals = top["current_t"]["data"][()]
                                        t_vals = top["current_t"]["grids/0/values"][()]
                                        steady_c = current_vals[len(current_vals) - 1]
                                        print "steady state current =", steady_c
                                        # rescale current - divide over steady state
                                        rescaled_current = np.abs(current_vals-steady_c)/np.abs(steady_c)
                                        filter_indices = rescaled_current >= current_threshold
                                        rescaled_current1 = rescaled_current[filter_indices]
                                        t_vals1 = t_vals[filter_indices]
                                        # take a log
                                        log_current = np.log(rescaled_current1)
                                        #print rescaled_current

                                        fit1 = lambda x, a, b : a - b*x
                                        fit_data = curve_fit(fit1, t_vals1, log_current)
                                        (a,b) = fit_data[0]
                                        (a_err, b_err) = (np.sqrt(fit_data[1][0][0]), np.sqrt(fit_data[1][1][1]))
                                        fit_vals[Gamma][bw][ed][tmax][v][T] = (a,b, a_err, b_err)
                                        currents[Gamma][bw][ed][tmax][v][T] = (t_vals, rescaled_current)
                                        print b, "exp(-", a, "x)"
                                        
                                    h5f.close()
                                os.chdir("..")  #exit tmax
                            os.chdir("..")  #exit beta
                        os.chdir("..")  #exit bw
                    os.chdir("..")  #exit h
                os.chdir("..")  #exit ed
            os.chdir("..")  #exit U
          os.chdir("..")  #exit U

        streams=dict()
        # finished fits
        # now save fitted data
        for G in fit_vals.keys():
          for bw in fit_vals[G].keys():
            for ed in fit_vals[G][bw].keys():
                for tmax in fit_vals[G][bw][ed].keys():
                    for v in fit_vals[G][bw][ed][tmax].keys():
                        streams["current"] = open(out_dir+os.path.sep+"cur_fit"+"_U"+str(U)+"ed"+str(ed)+"bw"+str(bw)+"v"+str(v)+".dat","w")
                        for T in sorted(fit_vals[G][bw][ed][tmax][v].keys()):
                            print "U"+str(U)+"ed"+str(ed)+"bw"+str(bw)+"T"+str(T)+"v"+str(v)
                            print fit_vals[G][bw][ed][tmax][v][T]
                            (a,b,a_err,b_err) = fit_vals[G][bw][ed][tmax][v][T]
                            streams["current"].write(str(T)+"   "+str(a)+" "+str(a_err)+"  "+str(b)+" "+str(b_err)+"\n")
                            dname = out_dir+os.path.sep+"current_t"+"_U"+str(U)+"ed"+str(ed)+"bw"+str(bw)+"T"+str(T)+"v"+str(v)+".dat"
                            t_vals = currents[G][bw][ed][tmax][v][T][0]
                            c_vals = currents[G][bw][ed][tmax][v][T][1]
                            f_v = np.exp(a - b * t_vals)
                            print "->", dname
                            np.savetxt(dname, np.vstack([t_vals, c_vals, f_v]).transpose())

        os.chdir("..")  #exit data_dir

def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
def find_nearest(array,value):
    return array[find_nearest_index(array,value)]

if __name__ == "__main__":
    main()
