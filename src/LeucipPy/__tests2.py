# 1 ######## TEST THE HELPER FUNCTIONS ##############################
import LeucipPy as leu
print(leu.getPdbLink('3nir'))



# 2 ######## TEST THE DATAFRAME ##############################
import DataFrameMaker as gdf
import Bio.PDB as bio
import os
from urllib.request import urlretrieve
parser = bio.PDBParser()
strucs = []
for pdb_code in ['3nir']:
  pdb_file = 'pdb' + pdb_code + '.ent'
  struc = parser.get_structure(pdb_code,pdb_file)
  strucs.append(struc)
geo = gdf.DataFrameMaker(strucs)
df = geo.calculateGeometry(['C:O'])

# 3 ####### Write an html report
import HtmlReportMaker as ghm
rep = ghm.HtmlReportMaker("Testing LeucipPy","test_leu2.html",remove_strip=False)
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
print("Cell entry", cf.getValue("_cell","entry_id"))
mtx = cf.getMatrix('_refln','index_l','index_k','F_meas_au',zz=[0,'index_h'])
mtx2 = cf.getMatrix('_refln','index_l','index_k','F_meas_au',zz=[50,'index_h'])
cub = cf.getCube('_refln','index_h','index_l','index_k','F_meas_au')
rep.changeColNumber(3)
rep.addSurface(mtx2,alpha=0.9,overlay=False,cmin=15,cmax=25)
rep.addSurface(mtx2,alpha=0.9,overlay=False,cmin=15)
rep.addSurface(mtx2,alpha=0.9,overlay=False,cmax=25,colourbar=True)
rep.addContours(mtx2,title='C',style='both',contourlabel=True,alpha=0.5,cmin=15,cmax=25)
rep.addContours(mtx2,title='C',style='both',contourlabel=True,alpha=0.5,cmin=15)
rep.addContours(mtx2,title='C',style='both',contourlabel=True,alpha=0.5,cmax=25)
rep.addSurface(mtx,alpha=0.9,overlay=False)
rep.addSurface(mtx,alpha=0.9,overlay=False,cap=0.5,title='capped')
rep.addSurface(mtx,alpha=0.9,overlay=True)
rep.addSurface(mtx2,alpha=0.5,overlay=False,palette='bone')

rep.addContours(mtx2,title='A',style='filled',contourlabel=True,alpha=0.5,overlay=True,palette='bone')
rep.addContours(mtx2,title='A',style='lines',contourlabel=True,alpha=0.9,colourbar=False,palette='rainbow')
rep.addContours(mtx2,title='B',style='filled',contourlabel=False,alpha=0.9,levels=[12,20,28])
rep.addContours(mtx2,title='B',style='filled',contourlabel=False,alpha=0.9,levels=[12,20,28],cmin=15,cmax=25)
rep.addContours(mtx2,title='C',style='both',contourlabel=True,alpha=0.5)
rep.addContours(mtx2,title='A',style='lines',contourlabel=False,alpha=0.9)
rep.addContours(mtx2,title='B',style='filled',contourlabel=False,alpha=0.9)
rep.addContours(mtx2,title='C',style='both',contourlabel=False,alpha=0.5,centred=True)

rep.printReport()