# 1 ######## TEST THE HELPER FUNCTIONS ##############################
import LeucipPy as leu


# 2 ######## TEST FE search ##############################
import GeoDataFrame as gdf
import Bio.PDB as bio
import os
from urllib.request import urlretrieve
parser = bio.PDBParser()
strucs = []
#for pdb_code in ['4u9h','5jsk','6rk0']:
for pdb_code in ['5jsk']:
  pdb_file = 'pdb' + pdb_code + '.ent'
  pdb_loc = ('https://www.ebi.ac.uk/pdbe/entry-files/download/pdb' + pdb_code + '.ent')
  if not os.path.exists(pdb_file):
      urlretrieve(pdb_loc, pdb_file)
  struc = parser.get_structure(pdb_code,pdb_file)
  strucs.append(struc)

geo = gdf.GeoDataFrame(strucs)
#df = geo.calculateGeometry(['FE:{O}','FE:{O@2}'])
#print(df)

#df = geo.calculateGeometry(['FE:{O,N,NE2}','FE:{O,N@2}'])
#print(df)

#df = geo.calculateGeometry(['FE:{O,N}+1','FE:{O,N,NE2@2}+1'])
#print(df)

#df = geo.calculateGeometry(['FE:{O,N,FE}+1','FE:{O,N,FE@2}+1'],log=0)
#print(df)
'''
df = geo.calculateGeometry(['FE:{NE2}'])
print(df)
df = geo.calculateGeometry(['FE:{O}'])
print(df)
df = geo.calculateGeometry(['FE:{NE2,O}'])
print(df)
df = geo.calculateGeometry(['FE:{O,NE2}+1'],log=0)
print(df)
df = geo.calculateGeometry(['FE:{O,NE2@2}+1'],log=0)
print(df)
df = geo.calculateGeometry(['FE:{O,NE2,FE}+1'],log=0)
print(df)
df = geo.calculateGeometry(['FE:{O,NE,FE@2}+1'],log=0)
print(df)
df = geo.calculateGeometry(['FE:(S)'],log=0)
print(df)
'''
df = geo.calculateGeometry(['(FE):(N,O)+1'],log=0)
print(df)
df = geo.calculateGeometry(['(FE):(N,O@2)+1'],log=0)
print(df)

#df = geo.calculateGeometry(['FE:{O,N,NE2,FE@2}+1'],log=4)
#print(df)