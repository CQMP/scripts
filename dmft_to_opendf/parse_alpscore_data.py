# A converter script between DMFT output and opendf input
import sys
sys.path.append('/home/jpfleblanc/.local/lib/python2.7/site-packages')


import numpy as np
print( np.version.version)
#from itertools import izip
import os
from numpy import pi as PI
import h5py

# Try to use fast table IO routines from pandas (does not exist on every machine -> optional)
use_pandas = False
try:
    import pandas
    use_pandas = True
except:
    print("No pandas found, results in slow loading")
    print( "do 'pip install pandas'")
    use_pandas = False

def main(params):
    # read input data from txt
    print("Parsing data with the following parameters\n","\n".join([str(key)+" : "+str(value) for key,value in params.items()]))
    (gw_grids, gw_data_in) = read_txt(params["gw"])
    (sigma_grids, sigma_data_in)   = read_txt(params["sigma"])
    (vertex_grids, vertex_data_in) = read_txt(params["vertex"])

    nflavors = gw_data_in.shape[0]/2
    print("There are", nflavors, "flavors")
    if ((sigma_data_in.shape[0]) / 2 != nflavors):
        print("Flavor mismatch between sigma and g :", nflavors, "!=", (sigma_data_in.shape[0]) / 2)
    if (not np.all(gw_grids[0] == sigma_grids[0])):
        print("Grid mismatch between sigma and g")

    fgrid_in = gw_grids[0]
    beta = 2*PI/(fgrid_in[1] - fgrid_in[0]) 
    print("beta =", beta)

    # process the data - combine real and imaginary parts
    mu = params["mu"]
    iw = fgrid_in*1j
    gw_data = np.zeros((gw_data_in.shape[0]/2, gw_data_in.shape[1]), dtype=np.complex)
    for x in range(nflavors):
        gw_data[x] = gw_data_in[2*x] + gw_data_in[2*x+1]*1j

    sigma_data = np.zeros(gw_data.shape, dtype=np.complex)
    vertex_data = np.zeros(np.append([nflavors], vertex_data_in[0].shape), dtype=np.complex)
    delta_data = sigma_data*0
    for x in range(nflavors):
        sigma_data[x] = sigma_data_in[2*x] + sigma_data_in[2*x+1]*1j
        vertex_data[x] = vertex_data_in[2*x] + vertex_data_in[2*x+1]*1j
        delta_data[x] = iw + mu - sigma_data[x] - 1./gw_data[x]

    # Postprocess data 
    F00 = -vertex_data[0]  # minus sign comes from difference between DGA and DF notations
    F01 = -vertex_data[1]
    bvertex_index_grid_in, fvertex_index_grid_in, _ = vertex_grids 

    # Create grids for df data
    wbmax = min(int(params["nbosonic"]), int(np.amax(bvertex_index_grid_in)) + 1);
    wfmax = min(int(params["nfermionic"]), int(np.amax(fvertex_index_grid_in)) + 1);
    b_index_grid = np.arange(-wbmax+1, wbmax, 1);
    f_index_grid = np.arange(-wfmax, wfmax, 1);
    bgrid = b_index_grid * 2 * PI / beta * 1j
    fgrid = (f_index_grid * 2 + 1) * PI / beta * 1j

    # Extract vertices for required number of fermionic/bosonic freqs
    fvertex_inds = np.array([np.searchsorted(fvertex_index_grid_in, x) for x in f_index_grid])
    bvertex_inds = np.array([np.searchsorted(bvertex_index_grid_in, x) for x in b_index_grid])
    vertex_ind_mesh = np.ix_(bvertex_inds, fvertex_inds, fvertex_inds)
    F00_filtered = F00[vertex_ind_mesh]
    F01_filtered = F01[vertex_ind_mesh]
    vertex_grids_out = (bgrid, fgrid, fgrid)

    # Symmetrize gf
    f_index_grid_gw = np.round((fgrid_in * beta / np.pi - 1 ) / 2).astype(int)
    sigma = np.zeros([nflavors, len(fgrid)], dtype=np.complex)
    delta = sigma*0
    gw = sigma*0
    for s in range(nflavors):
        for (obj,obj_orig) in zip((delta, sigma, gw),(delta_data,sigma_data,gw_data)):
            obj[s, wfmax : 2*wfmax] = obj_orig[s][0:wfmax]
            obj[s, 0:wfmax] = np.conjugate(obj_orig[s][0:wfmax][::-1])
    
    # output to hdf5
    data = h5py.File("qmc_output.h5", "w")
    top = data.create_group("dmft")
    save_grid_object(F00_filtered, vertex_grids_out , "F00", top)
    save_grid_object(F01_filtered, vertex_grids_out , "F01", top)
    for s in range(nflavors):
        save_grid_object(delta[s], [fgrid], "delta"+str(s), top) 
        save_grid_object(gw[s], [fgrid], "gw"+str(s), top) 
        save_grid_object(sigma[s], [fgrid], "sigma"+str(s), top) 
        
    data.close()
    
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def convert_multidimensional(data):
    ''' convert the flat grid+data structure into (grids,data) tuple. '''
    ncols = data.shape[0]
    nrows = data.shape[1]
    print(data.shape)

    shape_out=np.array([])
    #grids=np.array([[]])
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
    print(str(dims)+"d data,",shape_out, "dimensions,", ("complex" if complex_data else "real"), "data" )
    data_out = np.vstack([data[dims+x] for x in range(data_dims)])
    shape_out=map(int,shape_out)	
    array_out = np.reshape(data_out,shape_out)
    return (grids,array_out) 


def read_txt(fname):
    ''' read data from txt file. '''
    if not os.path.exists(fname):
        print("No such file:", fname + ".","Exiting.")
        exit(1)

    print("--> Loading", fname)
    f = open(fname,'rb')
    firstline=f.readline()
    f.close()
    has_header = not np.all([isfloat(x) for x in open(fname,'rb').readline().split()])
    print("file has",("a" if has_header else "no"), "header;",)
    
    all_data = None
    if use_pandas: 
        # Fast pandas read
        data1 = pandas.read_table(fname, parse_dates = False, header=None, comment = "#", delim_whitespace = True) #dayfirst = False, keep_date_col = True)
        all_data = data1.values.transpose() 
    else:
        # Slow genfromtxt from numpy
        all_data = np.genfromtxt(fname, dtype = np.float, skip_header = int(has_header), autostrip=True, unpack = True)
    print("shape :", all_data.shape)

    return convert_multidimensional(all_data)

def save_grid_object(data,grids,name,h5group):
    ''' dump (grids,data) structure to a gftools/ALPSCore compatible input '''
    print("--> saving", name)
    out = h5group.create_group(name)
    data_is_complex = (data.dtype == np.complex)
    data = data.astype(np.complex).view(np.float).reshape(np.append(data.shape,[2])) if data_is_complex else data 
    data_d = out.create_dataset("data", data = data)
    grids_d = out.create_group("grids")
    data_d.attrs["__complex__"] = int(data_is_complex)
    for (i,grid) in zip(range(len(grids)),grids):
        #print i, grid.shape
        grid_is_complex = grids[i].dtype == np.complex
        grid_g = grids_d.create_group(str(i))
        grid_data = grids[i] if not grid_is_complex else grids[i].astype(np.complex).view(np.float).reshape([len(grids[i]),2])
        grid_dataset = grid_g.create_dataset("values",data=grid_data)
        grid_dataset.attrs["__complex__"] = int(grid_is_complex)



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='bold_hyb')
    parser.add_argument('--mu', help='chemical potential', type=float, default = 0.000)
    parser.add_argument("--plaintext", default = 0, help = "save additionally to plaintext files")
    parser.add_argument("--nbosonic", default = 1024, help = "max number of non-negative bosonic Matsubara frequencies")
    parser.add_argument("--nfermionic", default = 1024, help = "max number of positive fermionic Matsubara frequencies")
    parser.add_argument("--vertex", default = "vertexF.dat", help = "vertex file") 
    parser.add_argument("--gw", default = "gw.dat", help = "gw file")
    parser.add_argument("--sigma", default = "sigma.dat", help = "self-energy file")
    args = parser.parse_args()
    main(vars(args))
    

