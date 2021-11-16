# !/usr/bin/python3
"""
GeoResidue (LeucipPy.GeoResidue))
------------------------------------------------------
Internal class to manage a single residue
"""


class GeoResidue:
    """Class for a single atom

    :data amino acid:
    :data atoms:
    """

    def __init__(self, amino_acid,rid,ridx):
        """Initialises a GeoPdb with a biopython structure

        :param biopython_structure: A list of structures that has been created from biopython
        """

        self.atoms = {}
        self.amino_acid = amino_acid
        self.rid = rid
        self.ridx = ridx






