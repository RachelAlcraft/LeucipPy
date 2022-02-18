
'''
TESTING creating geoemtry with atoms file
'''
import pandas as pd
import BioPythonMaker as bpm
import GeometryMaker as dfm
import HtmlReportMaker as hrm
import DsspMaker as dm

dic_a = {'pdb_code':[],'rid':[],'aa':[],'atom':[],'x':[],'y':[],'z':[]}
dic_a['pdb_code'].append('pdb0'),dic_a['rid'].append(1),dic_a['aa'].append('GLY'),dic_a['atom'].append('N'),dic_a['x'].append(7.224000),dic_a['y'].append(2.499000),dic_a['z'].append(8.390)
dic_a['pdb_code'].append('pdb0'),dic_a['rid'].append(1),dic_a['aa'].append('GLY'),dic_a['atom'].append('CA'),dic_a['x'].append(6.524000),dic_a['y'].append(1.546000),dic_a['z'].append(7.525)
dic_a['pdb_code'].append('pdb0'),dic_a['rid'].append(1),dic_a['aa'].append('GLY'),dic_a['atom'].append('C'),dic_a['x'].append(7.322801),dic_a['y'].append(2.562225),dic_a['z'].append(6.717)
dic_a['pdb_code'].append('pdb0'),dic_a['rid'].append(1),dic_a['aa'].append('GLY'),dic_a['atom'].append('O'),dic_a['x'].append(7.145820),dic_a['y'].append(2.797427),dic_a['z'].append(5.505)
dic_a['pdb_code'].append('pdb0'),dic_a['rid'].append(2),dic_a['aa'].append('GLY'),dic_a['atom'].append('N'),dic_a['x'].append(8.373636),dic_a['y'].append(3.023054), dic_a['z'].append(7.407)
dic_a['pdb_code'].append('pdb0'),dic_a['rid'].append(2),dic_a['aa'].append('GLY'),dic_a['atom'].append('CA'),dic_a['x'].append(9.373636),dic_a['y'].append(3.023054),dic_a['z'].append(7.407)
dic_a['pdb_code'].append('pdb0'),dic_a['rid'].append(2),dic_a['aa'].append('GLY'),dic_a['atom'].append('O'),dic_a['x'].append(11.359937),dic_a['y'].append(3.656131),dic_a['z'].append(4.628)
dic_a['pdb_code'].append('pdb0'),dic_a['rid'].append(3),dic_a['aa'].append('GLY'),dic_a['atom'].append('N'),dic_a['x'].append(10.306552),dic_a['y'].append(2.225277),dic_a['z'].append(6.024)
dic_a['pdb_code'].append('pdb0'),dic_a['rid'].append(3),dic_a['aa'].append('GLY'),dic_a['atom'].append('CA'),dic_a['x'].append(10.929683),dic_a['y'].append(1.276727),dic_a['z'].append(5.086)
dic_a['pdb_code'].append('pdb0'),dic_a['rid'].append(3),dic_a['aa'].append('GLY'),dic_a['atom'].append('C'),dic_a['x'].append(11.690921),dic_a['y'].append(1.540267),dic_a['z'].append(3.779)
dic_a['pdb_code'].append('pdb0'),dic_a['rid'].append(3),dic_a['aa'].append('GLY'),dic_a['atom'].append('O'),dic_a['x'].append(11.439531),dic_a['y'].append(1.047436),dic_a['z'].append(2.638)

adf = pd.DataFrame.from_dict(dic_a)
geo = dfm.GeometryMaker([],init_biopython=False,atoms_data=adf)
data = geo.calculateGeometry(['N:CA'],log=3)
print(data)

