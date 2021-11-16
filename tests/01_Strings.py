import os
from urllib.request import urlretrieve

import Bio.PDB as bio

parser = bio.PDBParser()
strucs = []

for pdb_code in pdb_codes:
    pdb_file = 'pdb' + pdb_code + '.ent'
    pdb_loc = ('https://www.ebi.ac.uk/pdbe/entry-files/download/pdb' + pdb_code + '.ent')
    if not os.path.exists(pdb_file):
        urlretrieve(pdb_loc, pdb_file)
    struc = parser.get_structure(pdb_code, pdb_file)
    strucs.append(struc)
