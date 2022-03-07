'''
RSA: 15th Feb 2022
This is a creation utility for training sets for a secondary structure classifier
'''

import pandas as pd
import BioPythonMaker as bpm
import GeometryMaker as dfm
import DsspMaker as dm

### INPUTS ################################
dsspPrintPath = 'Lists/'
pdbListPath = 'Lists/Pdb100.csv'
#geos = ['N:CA','CA:C','C:N+1','C:O','N:CA:C:N+1','C-1:N:CA:C','CA-1:C-1:N:CA','CA:C:N+1:CA+1','CA-1:CA:CA+1','CA-1:CA+1','CA:CA+1','CA:CA-1']
geos = ['N:O+2','N:O+3','N:O+4','O:N+2','O:N+3','O:N+4','O:{N}+2','N:{O}+2','N:CA:C','N:CA:C:N+1','C-1:N:CA:C']
tag = 'hb'
### System paths ########################## (Rachel's linux laptop)
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
#pdbDataPathLx =  '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
pdbDataPathLx =  'PdbData/'
outList = dsspPrintPath + 'SSTrainingSet_' + tag + '.csv'
### Get pdb list #################################
pdbdata = pd.read_csv(pdbListPath)
print(pdbdata)
pdblist = pdbdata['PDB'].tolist()[0:]
pdblist.sort()
pdblist = pdblist[:10]

### Load structures ######################
#strucs = bpm.loadPdbStructures(pdblist,pdbDataPathLx,extension='ent',prefix='pdb',log=1)
strucs = bpm.loadPdbStructures([],'PdbData/',extension='ent',prefix='pdb')
geo = dfm.GeometryMaker(strucs)
data = geo.calculateGeometry(geos)
print(data)

### add DSSP ######################
dssp = dm.DsspMaker(pdblist,pdbDataPathLx,extension='ent',prefix='pdb',log=1)
data = dssp.addDsspColumn(data)
print(data)

### Clean it into something for the training set ######################
cols = ['dssp','aa']
cols = cols + geos
data = data.query('dssp != "-"')
data = data[cols]

### save ######################
data.to_csv(outList,index=False)
print(data)










