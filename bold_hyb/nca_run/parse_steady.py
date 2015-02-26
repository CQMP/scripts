import os
import numpy as np
import scipy
import scipy.optimize
import string
import itertools
import h5py

def main():

    obs_binned = ["dm0", "dm1", "dm2", "dm3", "n0", "n1", "current"] 
    obs_mcdata = ["order", "sign"]
    obs_composite = ["m", "n"]

    data_dir = "data"
    top_string_array=["qmc", "nca"]

    origdir=os.getcwd()
    h5filename = "output.h5"

    for top_string in top_string_array:
        out_dir = top_string + "_stats"
        out_dir = os.path.abspath(out_dir)
        print "=================\n",data_dir, "data_dir\n","=================\n"
        print "Saving output to", out_dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        os.chdir(data_dir)
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
                        tmax_array = [5.0]
                        v_array = sorted(map(lambda x: eval(x.strip("v")), filter(lambda x: x.find("v")!=-1,os.listdir(os.getcwd()))),key = lambda x: float(x))
                        for tmax in tmax_array:
                            streams = dict()
                            for obs in obs_binned + obs_mcdata + obs_composite:
                                streams[obs] = open(out_dir+os.path.sep+obs+"_U"+str(U)+"ed"+str(ed)+"bw"+str(bw)+"T"+str(T)+"t"+str(tmax)+".dat","w")
                            for v in v_array:
                                os.chdir('v'+str(v))
                                os.chdir('tmax'+str(tmax))

                                print "============="
                                print "Statistics for",
                                print os.path.relpath(os.getcwd(),origdir)


                                if os.path.exists(h5filename):
                                  h5f = h5py.File(h5filename, "r")
                                  if top_string in h5f:
                                    top = h5f[top_string]
                                    h5data = top["stats"]
                                    # parse observables from stats
                                    for obs in obs_binned:
                                        try:
                                            obs_data = h5data[obs]
                                            (value, error) = obs_data
                                            value = abs(value) if abs(value) < 1e-3 else value
                                            print obs,":", value,"+/-",error
                                            streams[obs].write(str(v)+"   "+str(value)+" "+str(error)+"\n")
                                        except:
                                            print "Couldn't get", obs
                                    m_data = [(h5data["n0"][0] - h5data["n1"][0]), h5data["n0"][1]]
                                    streams["m"].write(str(v)+"   "+str(m_data[0])+" "+str(m_data[1]) + "\n")
                                    print "m:", m_data
                                    n_data = [(h5data["n0"][0] + h5data["n1"][0]), h5data["n1"][1]]
                                    streams["n"].write(str(v)+"   "+str(n_data[0])+" "+str(n_data[1]) + "\n")
                                    print "n:", n_data

                                    # parse observables from mcdata
                                    if top_string == "qmc":
                                        h5data = top["mcdata"]
                                        for obs in obs_mcdata:
                                            try:
                                                obs_data = h5data[obs]
                                                value = obs_data["mean"]["value"][()]
                                                error = obs_data["mean"]["error"][()]
                                                print obs,":", value,"+/-",error
                                                streams[obs].write(str(v)+"   "+str(value)+" "+str(error)+"\n")
                                            except:
                                                print "Couldn't get", obs
                                    
                                  h5f.close()
                                os.chdir("..")  #exit tmax
                                os.chdir("..")  #exit beta
                        os.chdir("..")  #exit bw
                    os.chdir("..")  #exit h
                os.chdir("..")  #exit ed
            os.chdir("..")  #exit U
          os.chdir("..")  #exit U
        os.chdir("..")  #exit data_dir

def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
def find_nearest(array,value):
    return array[find_nearest_index(array,value)]

if __name__ == "__main__":
    main()
