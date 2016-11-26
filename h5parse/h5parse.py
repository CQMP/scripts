#!/usr/local/bin/python3

import h5py
import argparse
import numpy as np
import sys

def main():
    parser = argparse.ArgumentParser(description="dump hdf5 to txt", usage="h5parse [-h] file datasets[:col,col,col] [datasets ...] [-o OUTPUT]")
    parser.add_argument('file',     type=str,             help="hdf5 archive")
    parser.add_argument('datasets', type=str, nargs="+",  help="datasets from hdf5 to dump to txt")
    parser.add_argument('-o', "--output", type=str,       help="output file")
    args = parser.parse_args()

    try:
        f = h5py.File(args.file, "r")
    except OSError:
        print("error: file %s does not exist" % args.file)
        return 1

    n = len(args.datasets)
    data = [0]*n

    for i in range(n):
        arg_in = args.datasets[i].split(":")
        datapath = arg_in[0]
        try:
            load_data = f[datapath].value
        except KeyError:
            print("error: dataset %s does not exist" % datapath)
            return 1

        if len(arg_in) == 1:
            data[i] = load_data
        elif len(arg_in) == 2:
            cols = [int(i) for i in arg_in[1].split(",")]
            if all(c < load_data.shape[1] for c in cols):
                data[i] = load_data[:, cols]
            else:
                print("error: column out of bounds in dataset %s" % datapath)
                return 1
        

    nrows = data[0].shape[0]
    ncols = 0
    for i in range(n):
        if len(data[i].shape) > 1:
            ncols += data[i].shape[1]
        else:
            ncols += 1

        if data[i].shape[0] != nrows:
            print("error datasets have incompatible shapes")
            return 0

    out = np.zeros((nrows, ncols))
    col = 0
    for i in range(n):
        if len(data[i].shape) > 1:
            for j in range(data[i].shape[1]):
                out[:, col] = data[i][:,j]
                col += 1
        else:
            out[:, col] = data[i]
            col += 1


    if (args.output):
        np.savetxt(args.output, out)
    else:
        np.savetxt(sys.stdout.buffer, out)

    return 0

if __name__ == "__main__":
    main()
