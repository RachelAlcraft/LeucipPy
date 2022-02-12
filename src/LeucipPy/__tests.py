# 1 ######## TEST THE HELPER FUNCTIONS ##############################
import LeucipPy as leu
print(leu.getPdbLink('1ejg'))
print(leu.getEbiLink('1ejg'))
print(leu.getEnghHuberStatistics())


# 2 ######## TEST THE DATAFRAME ##############################
import DataFrameMaker as gdf
import Bio.PDB as bio
import os
from urllib.request import urlretrieve
parser = bio.PDBParser()
strucs = []
for pdb_code in ['3nir','1ejg','4u9h','5jsk']:
  pdb_file = 'Data/pdb' + pdb_code + '.ent'
  pdb_loc = ('https://www.ebi.ac.uk/pdbe/entry-files/download/pdb' + pdb_code + '.ent')
  if not os.path.exists(pdb_file):
      urlretrieve(pdb_loc, pdb_file)
  struc = parser.get_structure(pdb_code,pdb_file)
  strucs.append(struc)

geo = gdf.DataFrameMaker(strucs)
df = geo.calculateGeometry(['C:O','C:N+1'])
df3 = geo.calculateGeometry(['FE:{O}','FE:{O@2}'])
print(df3)
#print(df)

# 2 ######## TEST THE DATAFRAME DATA ##############################
df2 = geo.calculateData()
#print(df2)

# 3 ####### Write an html report
import HtmlReportMaker as ghm
rep = ghm.HtmlReportMaker("Testing LeucipPy","test_leu.html")

rep.addPlot1d(df,'histogram','C:O',hue='rid',title='Hist',overlay=False)
rep.addPlot1d(df3,'histogram','FE:{O}',hue='rid',title='Hist',overlay=False)
rep.addPlot2d(df3,'scatter','FE:{O}','FE:{O@2}','bfactor',title='test',overlay=False,alpha=1)

rep.addPlot2d(df,'scatter','C:O','C:N+1','bfactor',title='test',overlay=False,alpha=0.1)
rep.addPlot2d(df,'scatter','C:O','C:N+1','bfactor',title='test',overlay=False,alpha=0.8)
rep.addPlot2d(df,'seaborn','C:N+1','C:O','pdb_code')
rep.addLineComment('Something to be said')
rep.changeColNumber(4)
rep.addPlot2d(df,'seaborn','C:O','C:N+1','aa')
rep.addBoxComment('Comment in the box')
rep.addPlot1d(df,'histogram','C:O',hue='rid',title='Hist',overlay=True)
rep.addPlot1d(df,'histogram','C:O',hue='aa',title='Hist',palette='seagreen',alpha=0.5)
rep.addPlot1d(df,'barstacked',['C:O','C:N+1'],hue='bfactor',title='Hist',xrange=[1.2,1.6],palette=['red','blue'],alpha=0.5)
#rep.addSeries(df['C:O'].describe(),'C:O')
rep.addSeries(df['C:O'].describe(),'C:O', transpose=True)
rep.changeColNumber(1)
rep.addDataFrame(leu.getEnghHuberStatistics(),'E&H')

import StatsThings as sts
same, res = sts.compareDistributionsMannWhitney(df['C:O'].values,df['C:N+1'].values,0.05, htmlise=True)


rep.addBoxComment(res)
rep.addBoxComment('Finished')
rep.printReport()

