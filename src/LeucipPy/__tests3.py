
import GeoDataFrame as gdf
import Bio.PDB as bio
import os
from urllib.request import urlretrieve
parser = bio.PDBParser()
strucs = []

for pdb_code in ['6zx4','3u7z','6vb2']:
  pdb_file = 'pdb' + pdb_code + '.ent'
  pdb_loc = ('https://www.ebi.ac.uk/pdbe/entry-files/download/pdb' + pdb_code + '.ent')
  if not os.path.exists(pdb_file):
      urlretrieve(pdb_loc, pdb_file)
  struc = parser.get_structure(pdb_code,pdb_file)
  strucs.append(struc)

geo = gdf.GeoDataFrame(strucs)
geoA = 'SG:{SG}+1'
geoB = 'SG:{N}+2'
geoC = 'SG:{SG@1}'
geoD = 'SG:{SG@2}'
geoE = 'SG:SG'
df = geo.calculateGeometry([geoA,geoB],log=1)
print(df)


import GeoHTML as ghm
import StatsThings as st
rep = ghm.GeoHTML("Testing LeucipPy Distances","test_dist.html")
rep.addPlot2d(df,'seaborn',geo_x=geoA,geo_y=geoB,hue='pdb_code',title='Cross Link NOS')
rep.addPlot2d(df,'scatter',geo_x=geoA,geo_y=geoB,hue='bfactor',title='Cross Link NOS',crange=[2,3],alpha=0.5)
rep.addPlot2d(df,'hist2d',geo_x=geoA,geo_y=geoB,hue='count',title='Cross Link NOS')
rep.addPlot2d(df,'probability',geo_x=geoA,geo_y=geoB,hue='count',title='Cross Link NOS')
a,tst = st.normalityShapiroWilk(df[geoA].values,0.05,htmlise=True)
rep.addBoxComment(tst)
rep.printReport()

print(df)
#dfa = geo.filferDataFrame(df,inclusions={'aa':['LYS','CYS']})
dfb = geo.filterDataFrame(df,inclusions={'aa':['CYS']})
#dfc = geo.filferDataFrame(df,exclusions={'aa':['CYS','GLY']})
#dfd = geo.filferDataFrame(df,exclusions={'aa':['PRO']})
#print(dfa)
print(dfb)
#print(dfc)
#print(dfd)