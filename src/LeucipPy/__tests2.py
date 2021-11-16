# 1 ######## TEST THE HELPER FUNCTIONS ##############################
import LeucipPy as leu
print(leu.getPdbLink('3nir'))



# 2 ######## TEST THE DATAFRAME ##############################
import GeoDataFrame as gdf
import Bio.PDB as bio
import os
from urllib.request import urlretrieve
parser = bio.PDBParser()
strucs = []
for pdb_code in ['3nir']:
  pdb_file = 'pdb' + pdb_code + '.ent'
  struc = parser.get_structure(pdb_code,pdb_file)
  strucs.append(struc)
geo = gdf.GeoDataFrame(strucs)
df = geo.calculateGeometry(['C:O'])

# 3 ####### Write an html report
import GeoHTML as ghm
rep = ghm.GeoHTML("Testing LeucipPy","test_leu2.html",remove_strip=False)
rep.addPlot1d(df,'histogram','C:O',hue='rid',title='Hist')
rep.addPlot1d(df,'histogram','C:O',hue='aa',title='Hist')
rep.addPlot1d(df,'histogram','C:O',hue='bfactor',title='Hist')
rep.changeColNumber(2)
rep.addSeries(df['C:O'].describe(),'C:O', transpose=True)
df = df.sort_values(by='C:O',ascending=True)
rep.addDataFrame(df,'C:O')
rep.addPlotPi(df,'C:O','aa','me')


# 4 ############# Test CifFile
filename = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/mmcif_data/1ejg-sf.cif'
import CifFile as cif
cf = cif.CifFile(filename)
mtx = cf.getMatrix('_refln','index_l','index_k','F_meas_au',zz=[0,'index_h'])
mtx2 = cf.getMatrix('_refln','index_l','index_k','F_meas_au',zz=[50,'index_h'])
cub = cf.getCube('_refln','index_h','index_l','index_k','F_meas_au')
rep.changeColNumber(3)
rep.addSurface(mtx2,alpha=0.9,overlay=False)
rep.addSurface(mtx,alpha=0.9,overlay=False)
rep.addSurface(mtx,alpha=0.9,overlay=True)
rep.addSurface(mtx2,alpha=0.5,overlay=False,palette='bone')

rep.printReport()