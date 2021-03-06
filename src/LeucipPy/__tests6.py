import WilliamsDivergenceMaker as wdm
import BioPythonMaker as bpm
import GeometryMaker as dfm
import HtmlReportMaker as hrm
import DsspMaker as dm

strucs = bpm.loadPdbStructures([],'Data/',extension='ent',prefix='pdb',log=2)
geo = dfm.GeometryMaker(strucs,log=2)
data = geo.calculateGeometry(['N:CA'])
dssp = dm.DsspMaker([],'Data/',extension='ent',prefix='pdb',log=2)
data = dssp.addDsspColumn(data)
print(data)
df = geo.calculateGeometry(['FE:{O}','FE:{O@2}'],log=2)
print(df)

#df = geo.calculateGeometry(['FE:{O,N,NE2}','FE:{O,N@2}'])
#print(df)

#df = geo.calculateGeometry(['FE:{O,N}+1','FE:{O,N,NE2@2}+1'])
#print(df)

#df = geo.calculateGeometry(['FE:{O,N,FE}+1','FE:{O,N,FE@2}+1'],log=0)
#print(df)

