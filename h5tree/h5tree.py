# -*- coding: utf-8 -*-
from __future__ import print_function
#*************************************
#*    h5tree - a python script       *
#*  to view an h5 file as a tree     *
#************************************
try:
    import h5py
except:
    print("This script requires h5py to be installed")
    print("visit: www.h5py.org")
    exit(1)

def print_tree(name):
    indent_level = name.count('/')
    basename = name.split("/")[-1]

    indents = ''
    #skip a level so it's easier to read
    if indent_level > 1:
        indents = '│     '
        indent_level-=2
    
    indents += '│  '*(indent_level) +'├─ '

    if name.count('/') > 0:
        print("{}{}".format(indents,basename))
    else:
        #add bold to top level names
        print("{}\033[1m{}\033[0m".format(indents,basename))
        


def main(fname):
    f = h5py.File(fname,'r')
    f.visit(print_tree)
    f.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='h5tree - view h5 file in a tree structure')
    parser.add_argument("filename",help="the h5 file")
    args = parser.parse_args()
    main(args.filename)
