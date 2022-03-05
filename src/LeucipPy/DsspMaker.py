# !/usr/bin/python3
"""
GeoDataFrame (LeucipPy.GeoDataFrame))
------------------------------------------------------
This class manipulates given biopython structures and creates dataframes oft he geometry
"""
import os
import Bio.PDB as bio
from os import listdir
from os.path import isfile, join

import pandas as pd
try:
    from LeucipPy import GeoPdb as pdb
except:
    import GeoPdb as pdb


class DsspMaker:
    '''
    '''
    def __init__(self,pdbs, directory, extension='pdb', prefix='', log=0, count=0,start=0):
        self.dssp_dic = {}
        struc_pdb_files = []
        if directory == '':
            directory = str(os.getcwd())
        if directory[-1:] != '/':
            directory += '/'
        if log > 0:
            print('LeucipPy(1) Loading dssp in', directory, extension, prefix)
        parser = bio.PDBParser()
        if pdbs == []:  # Then we load all matching files in the directory
            onlyfiles = [f for f in listdir(directory) if isfile(join(directory, f))]
            length = len(onlyfiles)
            if count > 0:
                length = min(count+start, length)
            for c in range(start, length):
                fn = onlyfiles[c]
                if ('.' + extension) in fn:
                    fnp = directory + fn
                    fndot = fn.split('.')
                    fnfile = fndot[0].split('/')
                    fnnam = fnfile[len(fnfile) - 1:][0]
                    try:
                        if log > 1:
                            print('LeucipPy(2) DSSP', str(c+1) + '/' + str(length),'-', fnp)
                        struc = parser.get_structure(fnnam, fnp)
                        struc_pdb_files.append([struc,fnnam,fnp])
                    except:
                        print('!LeuciPy BioPythonMaker:structure failed', fnp)
        else:
            count = 0
            for pdb in pdbs:
                count += 1
                filename = directory + prefix + pdb + '.' + extension
                try:
                    if log > 1:
                        print('LeucipPy(2) DSSP', str(count)+'/'+str(len(pdbs)),'-', filename)
                    struc = parser.get_structure(pdb, filename)
                    struc_pdb_files.append([struc,pdb,filename])
                except:
                    print('!LeuciPy BioPythonMaker:structure failed', filename)

        from Bio.PDB.DSSP import DSSP
        for structure,pdb,filename in struc_pdb_files:
            model = structure[0]
            dssp = DSSP(model, filename)
            for akey in list(dssp.keys()):
                chain = akey[0]
                res_no = akey[1][1]
                row = dssp[akey]
                ss = row[2]
                self.dssp_dic[pdb+chain+str(res_no)] = ss

    def addDsspColumn(self,data,pdb_col = 'pdb_code',chain_col = 'chain',res_col = 'rid',col_name = 'dssp'):
        dssp_col = []
        for i in range(0,len(data.index)):
            pdb = data[pdb_col].values[i]
            chain = data[chain_col].values[i]
            res = str(data[res_col].values[i])
            if pdb+chain+res in self.dssp_dic:
                dssp_col.append(self.dssp_dic[pdb+chain+res])
            else:
                dssp_col.append('-')
        data[col_name] = dssp_col
        return data





