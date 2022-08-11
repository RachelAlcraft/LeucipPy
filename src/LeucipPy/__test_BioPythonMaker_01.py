# 1 ######## TEST THE CLASS BioPythonMaker ##############################

print("")
print("")
print("### UNIT TESTING BIOPYTHONMAKER ###")
import os
import BioPythonMaker as bpm
from urllib.request import urlretrieve
from Bio import BiopythonWarning
import warnings
warnings.simplefilter("ignore",BiopythonWarning)
import time

pdb_codes = ['3nir','1ejg','4u9h','5jsk','3nir','1ejg','4u9h','5jsk','3nir','1ejg','4u9h','5jsk']#,'3nir','1ejg','4u9h','5jsk','3nir','1ejg','4u9h','5jsk','3nir','1ejg','4u9h','5jsk','3nir','1ejg','4u9h','5jsk','3nir','1ejg','4u9h','5jsk','3nir','1ejg','4u9h','5jsk','3nir','1ejg','4u9h','5jsk','3nir','1ejg','4u9h','5jsk','3nir','1ejg','4u9h','5jsk','3nir','1ejg','4u9h','5jsk','3nir','1ejg','4u9h','5jsk','3nir','1ejg','4u9h','5jsk']

prefix = "pdb"
suffix = "ent"


code_file_location_hry = os.path.dirname(os.path.realpath(__file__)).split("/")
print(code_file_location_hry)
data_location_path = "/".join(code_file_location_hry) + "/PdbData/"

for pdb_code in pdb_codes:
  pdb_file = data_location_path + prefix + pdb_code + '.' + suffix
  pdb_loc = ('https://www.ebi.ac.uk/pdbe/entry-files/download/pdb' + pdb_code + '.' + suffix)
  if not os.path.exists(pdb_file):
    urlretrieve(pdb_loc, pdb_file)
      
jobs = [1,2,4,8]
for job in jobs:
  print('### Load structures from BioPython USING jobs=',job,"####")
  t = time.time()
  strucs = bpm.loadPdbStructures(pdb_codes,data_location_path,prefix=prefix,extension=suffix,log=1,jobs=job)
  print(time.time()-t)
  print("### Loaded",len(strucs),"structures from BioPython #############")
  #for struc in strucs:
  #  print(struc)
