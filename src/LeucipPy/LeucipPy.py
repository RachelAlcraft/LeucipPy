# !/usr/bin/python3
"""
LeucipPy is a protein geometry library
------------------------------------------------------
These top level functions are helper functions to retrieve data from the various websites
"""

import pandas as pd

def getEbiLink(pdb_code):
    """Returns the html path to the pdb file on the ebi server
    """

    pdb_loc = 'https://www.ebi.ac.uk/pdbe/entry/pdb/' + pdb_code
    return pdb_loc

def getPdbLink(pdb_code):
    """Returns the html path to the pdb file on the ebi server
    """

    file_name = 'pdb' + pdb_code + '.ent'
    pdb_loc = 'https://www.ebi.ac.uk/pdbe/entry-files/download/' + file_name
    return file_name, pdb_loc

def getElectronDensityLink(pdb_code, diff=False):
    """Returns the html path to the 2FoFc electron density file on the ebi server
    :param diff: Optionaly, defaults to False, if true gets the Fo-Fc, if Flase gets the 2Fo-Fc
    """

    file_name = pdb_code + '.ccp4'
    if diff:
        file_name = pdb_code + '_diff.ccp4'
    pdb_loc = 'https://www.ebi.ac.uk/pdbe/coordinates/files/' + file_name
    return file_name, pdb_loc

def getStructureFactorLink(pdb_code):
    """Returns the html path to the structure factors file on the ebi server
    """

    file_name = 'r' + pdb_code + 'sf.ent'
    pdb_loc = 'https://www.ebi.ac.uk/pdbe/entry-files/download/' + file_name
    return file_name, pdb_loc

def getMtzLink(pdb_code):
    """Returns the html path to the mtz file on the pdb server
    """

    file_name = 'r' + pdb_code + '_phases.mtz'
    pdb_loc = 'https://edmaps.rcsb.org/coefficients/' + file_name
    return file_name, pdb_loc

def getPdbCifLink(pdb_code):
    """Returns the html path to the cif file on the pdb server
    """

    file_name = pdb_code + '.cif'
    pdb_loc = 'https://files.rcsb.org/download/' + file_name
    return file_name, pdb_loc

def getStructureFactorCifLink(pdb_code):
    """Returns the html path to the cif file on the pdb server
    """

    file_name = pdb_code + '-sf.cif'
    pdb_loc = 'https://files.rcsb.org/download/' + file_name
    return file_name, pdb_loc

def getEnghHuberStatistics(type='BackboneBondLengths'):
    """Returns a dataframe of the E&H summary statistics from 1991

    These results are taken from the Engh and Huber International Tables for Ctsyallography (2006). Vol F, Chapter 18.3, pp 382-192 (RA Engh and R Huber)

    :param: type currently only have backone bond lengths
    """

    cols = ['statistic','year']
    cols.extend(['N:CA','CA:C', 'C:O', 'C:CB', 'C:N+1'])
    cols.extend(['N:CA:PRO','CA:C:PRO','C:O:PRO','C:CB:PRO','C:N+1:PRO'])
    cols.extend(['N:CA:GLY', 'CA:C:GLY', 'C:O:GLY', 'C:CB:GLY', 'C:N+1:GLY'])
    vals = []
    vals.append(['mean','2001',    1.459, 1.525, 1.229, 1.532, 1.336,  1.468, 1.524, 1.228, 1.531, 1.338,  1.456, 1.514, 1.232, 1.532, 1.326])
    vals.append(['mean-1', '2001', 1.439, 1.499, 1.210, 1.499, 1.313,  1.451, 1.504, 1.208, 1.511, 1.319,  1.441, 1.498, 1.216, 1.499, 1.308])
    vals.append(['mean+1', '2001', 1.479, 1.551, 1.248, 1.563, 1.359,  1.485, 1.544, 1.248, 1.551, 1.357,  1.471, 1.530, 1.248, 1.563, 1.344])
    vals.append(['mean', '1999',   1.458, 1.525, 1.231, 1.530, 1.329,  1.466, 1.525, 1.231, 1.530, 1.341,  1.451, 1.516, 1.231, 1.530, 1.329])
    vals.append(['mean-1', '1999', 1.439, 1.504, 1.211, 1.510, 1.305,  1.451, 1.504, 1.211, 1.510, 1.325,  1.435, 1.498, 1.211, 1.510, 1.305])
    vals.append(['mean+1', '1999', 1.477, 1.546, 1.251, 1.550, 1.343,  1.481, 1.546, 1.251, 1.550, 1.357,  1.467, 1.534, 1.251, 1.550, 1.343])
    df = pd.DataFrame(vals, columns=cols)
    return df

def loadSimpleMatrix(filepath):
    import numpy as np
    with open(filepath,'r') as f:
        ed_data = f.read().splitlines()
    rows = len(ed_data)
    ed_slice = np.zeros((rows,rows))
    for i in range(0,rows):
        row = ed_data[i].split(',')
        for j in range(0, rows):
            val = float(row[j])
            ed_slice[i,j] = val
    return ed_slice