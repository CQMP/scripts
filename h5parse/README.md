h5parse
======
python script for dumping hdf5 files to text.
includes bash completion script for datasets within hdf5 file.
completion requires bash version > 4.0

### Usage
$ h5parse [-h] file.h5 dataset[:col,col,col] [dataset ...] [-o OUTPUT]

### Examples
h5parse output.h5 /qmc/gtau0/grids/0/values /qmc/gtau0/data:1 /qmc/gtau0_err/data:0 -o gtau.dat
