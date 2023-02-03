# 1 ######## TEST THE HELPER FUNCTIONS ##############################
import LeucipPy as leu
print(leu.getPdbLink('3nir'))



# 2 ######## TEST THE DATAFRAME ##############################
import GeometryMaker as gdf
import BioPythonMaker as bpm
import os
from urllib.request import urlretrieve

code_file_location_hry = os.path.dirname(os.path.realpath(__file__)).split("/")
data_location_path = "/".join(code_file_location_hry) + "/PdbData/"
print("Data location",data_location_path)

strucs = bpm.loadPdbStructures(["2ayw"],data_location_path,extension='ent',prefix='pdb')  
geo_mak = gdf.GeometryMaker(strucs,log=0,exc_hetatm=True)
df = geo_mak.calculateGeometry(['C:N+1'])
df = df.sort_values(by='C:N+1',ascending=True)
print(df[['C:N+1','rid','aa','aa+1']])

# 3 ####### Write an html report
import HtmlReportMaker as ghm
rep = ghm.HtmlReportMaker("Testing LeucipPy","test_leu10.html",remove_strip=False)
rep.addPlot1d(df,'histogram','C:N+1',hue='rid',title='Hist')
rep.addPlot1d(df,'histogram','C:N+1',hue='aa',title='Hist')
rep.addPlot1d(df,'histogram','C:N+1',hue='bfactor',title='Hist')
rep.changeColNumber(2)
rep.addSeries(df['C:N+1'].describe(),'C:N+1', transpose=True)
df = df.sort_values(by='C:N+1',ascending=True)
rep.addDataFrame(df,'C:N+1')
rep.addPlotPi(df,'C:N+1','aa','me')

rep.printReport()