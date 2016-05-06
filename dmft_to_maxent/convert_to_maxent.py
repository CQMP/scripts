# a script to convert DMFT output for input to maxent
# also works with plain input files

from __future__ import print_function,division
import numpy as np
import os
import h5py
from math import pi

#inspired from a script by Andrey

def main(params):
    print("Parsing data with the following parameters:\n","\n ".join([str(key)+" : "+str(value) for key,value in params.items()]))

    use_hdf5 = not params["no_hdf5"]
    ph_sym = params["ph"]
    const_sigma = params["const_sigma"]

    in_file = params["filename"]
    
    # filename is script run directory (working directory) 
    # + basename that is passed in
    out_file = os.getcwd()
    out_file += '/'+os.path.basename(in_file) 
    print(out_file)

    if not use_hdf5:
        #grab from textfile
        (gw_im_grids, gw_im_data) = read_txt(params["gw_imag"])
        #hard code in the first column of the file.
        gw_im_data = gw_im_data[0]
        gw_im_grids = gw_im_grids[0]
        print("using first column of data")

        #make real vectors just in case
        gw_re_data = []
        gw_re_grids = []
        if not ph_sym:
            (gw_re_grids, gw_re_data) = read_txt(params["gw_real"]) if os.path.exists(params["gw_real"]) else (gw_im_grids,gw_im_data*0)
            gw_re_data = gw_re_data[0]
        if const_sigma >0:
            print("generating constant sigma values")
            sigma_im = len(gw_im_data)*[const_sigma]
            sigma_re = []
            if not ph_sym: sigma_re = len(gw_re_data)*[const_sigma]
        else:
            print("Error! Please supply constant error value")
            exit(1)

        write_to_file(out_file,gw_im_grids,gw_re_data,gw_im_data,sigma_re,sigma_im,ph_sym)

     
    else:
        #grab from hdf5
        f = h5py.File(in_file,"r")
        beta = float(f["/parameters/BETA"].value)
        print("found beta =",beta)
        print("using spin up data")
        gw_re_data = f["/simulation/results/G_omega_up_re0/mean/value"].value
        gw_im_data = f["/simulation/results/G_omega_up_im0/mean/value"].value
        print("found",len(gw_re_data), "data points")

        if const_sigma >0:
            print("generating constant sigma values")
            sigma_re = len(gw_re_data)*[const_sigma]
            sigma_im = len(gw_im_data)*[const_sigma]
        else:
            print("using h5 errors")
            sigma_re = f["/simulation/results/G_omega_up_re0/mean/error"].value
            sigma_im = f["/simulation/results/G_omega_up_im0/mean/error"].value

        gw_grid = iwn_grid(beta,len(gw_im_data)) 

        write_to_file(out_file,gw_grid,gw_re_data,gw_im_data,sigma_re,sigma_im,ph_sym)
        create_param_file(out_file,beta,len(gw_im_data),ph_sym)

def write_to_file(filename,gw_grid,gw_re_data,gw_im_data,sigma_re,sigma_im,ph_sym):
    '''print out proper data format'''

    #constants for writing a file
    s = ' '
    nl = '\n'

    fname, f_extension = os.path.splitext(filename)
    data_file = open(fname+"_maxent_in","w+")
    print("creating data file: "+fname+"_maxent_in")

    for i in range(len(gw_im_data)):
       to_write = str(gw_grid[i])+s
       if not ph_sym:
           to_write += str(gw_re_data[i])+s+str(sigma_re[i])+s
       to_write += str(gw_im_data[i])+s+str(sigma_im[i])+nl
       data_file.write(to_write)

    data_file.close()

def create_param_file(filename,beta,ndat,ph_sym):
    '''create minimum param file for maxent'''
    fname, f_extension = os.path.splitext(filename)
    param_file = open(fname+"_maxent.param","w+")
    data_str = fname+"_maxent_in"
    print("creating param file: "+fname+"_maxent.param")
    #write values
    nl = '\n'
    param_file.write("#Generated param file"+nl)
    param_file.write("DATASPACE=frequency"+nl)
    param_file.write("BETA="+str(beta)+nl)
    param_file.write("NDAT="+str(ndat)+nl)
    param_file.write("DATA="+data_str+nl)
    if ph_sym:
        param_file.write("PARTICLE_HOLE_SYMMETRY=1"+nl)

def convert_multidimensional(data):
    ''' convert the flat grid+data structure into (grids,data) tuple. '''
    ncols = data.shape[0]
    nrows = data.shape[1]

    shape_out=np.array([])
    grids = []

    for i in range(ncols):
        grid = np.unique(data[i])
        current_dim = len(grid)
        #print "dim[",i,"] = ", current_dim
        shape_out = np.append(shape_out,[current_dim])
        grids.insert(ncols-i, grid)
        nrows = nrows / current_dim
        if (nrows == 1): 
            break
    dims = len(shape_out)
    data_dims = ncols - dims
    shape_out = np.insert(shape_out,0,data_dims) if data_dims > 1 else None
    complex_data = False # (ncols - dims)==2
    print(str(dims)+"d data,",shape_out, "dimensions,", ("complex" if complex_data else "real"), "data")
    data_out = np.vstack([data[dims+x] for x in range(data_dims)])
    array_out = np.reshape(data_out,shape_out)
    return (grids,array_out) 

def read_txt(fname):
    ''' read data from txt file. '''
    if not os.path.exists(fname):
        print ("No such file:", fname + ".","Exiting.")
        exit(1)

    print("--> Loading", fname)
    firstline=open(fname,'rb').readline()
    has_header = not np.all([isfloat(x) for x in open(fname,'rb').readline().split()])
    print("file has",("a" if has_header else "no"), "header;",)
    
    all_data = None
    # Slow genfromtxt from numpy
    all_data = np.genfromtxt(fname, dtype = np.float, skip_header = int(has_header), autostrip=True, unpack = True)
    print("shape :", all_data.shape)

    return convert_multidimensional(all_data)

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def iwn_grid(beta,n):
    '''generate grid of matsubara frequencies of length n'''
    grid = []
    for i in range(n):
        grid.append((2*i+1)*pi/beta)
    return grid

if __name__ == "__main__":
    import argparse
    desc = '''DMFT to Maxent'''
    ep= "pass an HDF5 file from DMFT or the seperate real and imaginary text files. For HDF5 input a param file will be created"
    parser = argparse.ArgumentParser(description=desc,epilog=ep,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--filename',default="sim.h5",help="HDF5 input file")
    parser.add_argument('--no_hdf5', help='The input is not an HDF5 file',default = False,action='store_const', const=True)
    parser.add_argument("--gw_real", default = "G_omegareal.dat", help = "file with the real part of gw"
                                                                         +"\n not needed for HDF5 input"
                                                                         +"\n (default: G_omegareal.dat)")
    parser.add_argument("--gw_imag", default = "G_omega.dat", help = "file with the imag part of gw"
                                                                     +"\n not needed for HDF5 input"
                                                                         +"\n (default: G_omega.dat)")
    parser.add_argument("--ph",default=False,action='store_const', const=True, help="if data is Particle-Hole Symmetric")
    parser.add_argument("--const_sigma",type=float,default=-1,help="add a constant error bar")
    args = parser.parse_args()
    main(vars(args))
