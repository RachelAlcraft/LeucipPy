'''
This file interfaces with biopython to create biopython structures which are used in the geopemtryt calculations
'''

from os import listdir
from os.path import isfile, join
import os
import Bio.PDB as bio
import pandas as pd
from joblib import Parallel, delayed
from Bio import BiopythonWarning
import warnings
warnings.simplefilter("ignore",BiopythonWarning)

try:
    from LeucipPy import GeometryMaker as gm  
except:
    import GeometryMaker as gm
    
def loadPdbStructures_loop(onlyfiles,directory,log,length,j,chunk,geos,exc_hetatm):    
    #https://stackoverflow.com/questions/9786102/how-do-i-parallelize-a-simple-python-loop
    parser = bio.PDBParser()
    rets = []
    for i in range(j,j+chunk):
        if i < len(onlyfiles):
            fn = onlyfiles[i]    
            fnp = directory + fn
            fndot = fn.split('.')
            fnfile = fndot[0].split('/')
            fnnam = fnfile[len(fnfile) - 1:][0]
            #try:
            
            if True:
                if log > 1:
                    print('LeucipPy(2) BIO',str(i+1) + '/' + str(length),'-',fnp)
                struc = parser.get_structure(fnnam,fnp)        
                if geos != []:
                    geo = gm.GeometryMaker([struc], log=0,exc_hetatm=exc_hetatm)
                    data = geo.calculateGeometry(geos, log=0)
                    rets.append(data)
                else:
                    rets.append(struc)
            #except:
            #    print('!LeuciPy BioPythonMaker:structure failed',fnp)
    return rets
    
    
def loadPdbStructures(pdbs,directory,extension='pdb',prefix='',log=0,count=0,start=0,jobs=1,geos=[],exc_hetatm=True):
    # First load all the pdb files in the given folder as biopython structures
    if directory == '':
        directory = str(os.getcwd())
    if directory[-1:] != '/':
        directory += '/'
    if log > 0:
        print('LeucipPy(1) Loading BioPython in', directory, extension, prefix)    
    strucs = []
    onlyfiles = []
    if pdbs == []: #Then we load all matching files in the directory
        onlyfilesa = [f for f in listdir(directory) if isfile(join(directory, f))]
        for c in range(start, len(onlyfilesa)):
            fn = onlyfilesa[c]
            if ('.' + extension) in fn:
                onlyfiles.append(fn)
    else:        
        for pdb in pdbs:            
            onlyfiles.append(prefix + pdb + '.' + extension)                    
    length = len(onlyfiles)
    if count > 0:
        length = min(count+start,length)
    #if jobs > 1:        
    chunk = 10
    strucs = Parallel(n_jobs=jobs)(delayed(loadPdbStructures_loop)(onlyfiles,directory,log,length,i,chunk,geos,exc_hetatm) for i in range(start,length,chunk))
    structs = []
    for st in strucs:
        for s in st:
            structs.append(s)

    if geos == []:   
        return structs
    else: # then we are actually returning a dataframe to save on memory
        df = pd.concat(structs)
        return df
