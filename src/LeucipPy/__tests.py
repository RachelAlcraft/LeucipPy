# 1 ######## TEST THE HELPER FUNCTIONS ##############################
import LeucipPy as leu
print(leu.getPdbLink('1ejg'))
print(leu.getEbiLink('1ejg'))
print(leu.getEnghHuberStatistics())


# 2 ######## TEST THE DATAFRAME ##############################
import GeoDataFrame as gdf
import Bio.PDB as bio
import os
from urllib.request import urlretrieve
parser = bio.PDBParser()
strucs = []
for pdb_code in ['3nir','1ejg']:
  pdb_file = 'pdb' + pdb_code + '.ent'
  pdb_loc = ('https://www.ebi.ac.uk/pdbe/entry-files/download/pdb' + pdb_code + '.ent')
  if not os.path.exists(pdb_file):
      urlretrieve(pdb_loc, pdb_file)
  struc = parser.get_structure(pdb_code,pdb_file)
  strucs.append(struc)

geo = gdf.GeoDataFrame(strucs)
df = geo.calculateGeometry(['C:O','C:N+1'])
#print(df)

# 2 ######## TEST THE DATAFRAME DATA ##############################
df2 = geo.calculateData()
#print(df2)

# 3 ####### Write an html report
import GeoHTML as ghm
rep = ghm.GeoHTML("Testing LeucipPy","test_leu.html")
rep.addPlot2d(df,'scatter','C:O','C:N+1','bfactor',title='test')
rep.addPlot2d(df,'seaborn','C:N+1','C:O','pdb_code')
rep.addLineComment('Something to be said')
rep.changeColNumber(4)
rep.addPlot2d(df,'seaborn','C:O','C:N+1','aa')
rep.addBoxComment('Comment in the box')
rep.addPlot1d(df,'histogram','C:O',hue='rid',title='Hist')
rep.addPlot1d(df,'histogram','C:O',hue='aa',title='Hist')
rep.addPlot1d(df,'histogram','C:O',hue='bfactor',title='Hist')
#rep.addSeries(df['C:O'].describe(),'C:O')
rep.addSeries(df['C:O'].describe(),'C:O', transpose=True)
rep.changeColNumber(1)
rep.addDataFrame(leu.getEnghHuberStatistics(),'E&H')



rep.addBoxComment('Finished')
rep.printReport()

