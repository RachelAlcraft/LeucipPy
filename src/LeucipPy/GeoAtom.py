# !/usr/bin/python3
"""
GeoAtom (LeucipPy.GeoAtom))
------------------------------------------------------
Internal class to manage a single atom
"""

class GeoAtom:
    """Class for a single atom

    :data atom_type:
    :data atom_name:
    :data x:
    :data y:
    :data z:
    :data disordered:
    :data occupancy:
    :data bfactor:
    """

    def __init__(self, atom_type,atom_name,atom_no,disordered,occupancy,bfactor,x,y,z):
        """Initialises a GeoPdb with a biopython structure

        :param biopython_structure: A list of structures that has been created from biopython
        """

        self.info = {}
        self.atom_type = atom_type
        self.atom_name = atom_name
        self.atom_no = atom_no
        self.disordered = disordered
        self.occupancy = occupancy
        self.bfactor = bfactor
        self.x = x
        self.y = y
        self.z = z




