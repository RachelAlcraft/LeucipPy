'''
This file interfaces with biopython to create biopython structures which are used in the geopemtryt calculations
'''

from os import listdir
from os.path import isfile, join
import os
import Bio.PDB as bio
import pandas as pd


def loadPdbStructures(pdbs,directory,extension='pdb',prefix='',log=0,count=0):
    # First load all the pdb files in the given folder as biopython structures
    if directory == '':
        directory = str(os.getcwd())
    directory += '/'
    if log > 0:
        print('LeucipPy(1) Loading BioPython in', directory, extension, prefix)
    parser = bio.PDBParser()
    strucs = []
    if pdbs == []: #Then we load all matching files in the directory
        onlyfiles = [f for f in listdir(directory) if isfile(join(directory, f))]
        length = len(onlyfiles)
        if count > 0:
            length = min(count,length)
        for c in range(0, length):
            fn = onlyfiles[c]
            if ('.' + extension) in fn:
                fnp = directory + fn
                fndot = fn.split('.')
                fnfile = fndot[0].split('/')
                fnnam = fnfile[len(fnfile) - 1:][0]
                try:
                    if log > 1:
                        print('LeucipPy(2)',fnp)
                    struc = parser.get_structure(fnnam,fnp)
                    strucs.append(struc)
                except:
                    print('!LeuciPy BioPythonMaker:structure failed',fnp)
    else:
        for pdb in pdbs:
            filename = directory + prefix + pdb + '.' + extension
            try:
                if log > 1:
                    print('LeucipPy(2)', filename)
                struc = parser.get_structure(pdb, filename)
                strucs.append(struc)
            except:
                print('!LeuciPy BioPythonMaker:structure failed', filename)

    return strucs