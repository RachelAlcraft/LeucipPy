# !/usr/bin/python3
"""
GeoPdb (LeucipPy.GeoPdb))
------------------------------------------------------
Internal class to manage a single pdb
"""

import pandas as pd
try:
    from LeucipPy import GeoAtom as atm
    from LeucipPy import GeoResidue as res
    from LeucipPy import GeoCalculate as calc
except:
    import GeoAtom as atm
    import GeoResidue as res
    import GeoCalculate as calc

class GeoPdb:
    """CLass for a single pdb

    :data pdb_code: The pdb code
    :data resolution: The resolution of the structure
    :data chains: a dictionary of the chains which is in turn a list of residues

    """
    def __init__(self,biopython_structure):
        """Initialises a GeoPdb with a biopython structure

        :param biopython_structure: A list of structures that has been created from biopython
        """

        self.bio_struc = biopython_structure
        self.chains = {}
        self.createAtomsList()


    def calculateGeometry(self,chain, resno, geo_atoms,log=0):
        """Creates the geoemtry from the structure in the class for 1 geo

        :param chain: The chain id
        :param rid: The residue number
        :param geo_atoms: A list of tuples that is the atom, the displacement

        :returns: [bool,float,bfactor, occupancy, atoms, GeoAtom] a bool for if it could be calculated, and the value, and the reference atom
        """
        total_rid = 0
        total_ridx = 0
        other = ''
        is_nearest = False

        rids = self.chains[chain]
        all_there = True
        atoms = []
        for atom,offset in geo_atoms:
            this_rid = resno + offset
            if this_rid not in rids:
                return False,0,0,0,0,0,0,None,''
            else:
                rid = rids[this_rid]
                if '{' in atom and '}' in atom: #then we need to do a distance search first and iut must be the second atom
                    #print('dbg1', atom)
                    is_nearest = True
                    refatm = atoms[0]
                    atomlist = atom[1:len(atom) - 1]
                    nearest = 1
                    if "@" in atomlist:
                        #print('dbg1', atomlist)
                        ats = atomlist.split('@')
                        #print('dbg1',ats)
                        atomlist = ats[0]
                        nearest =ats[1]
                    else:
                        atomlist = atom[1:len(atom)-1]
                    atomlist = atomlist.split(',')
                    debugatomlist = atom[1:len(atom)-1].split(',')
                    last_atm = None
                    last_dis = 10000
                    last_other = ""
                    for atmtype in atomlist:
                        dis,this_atm,other2 = self._getNearestAtom(refatm, resno,chain,atmtype,offset,nearest,log)

                        if dis < last_dis:
                            last_atm = this_atm
                            last_dis = dis
                            last_other = other2

                    other += "_" + last_other
                    if last_other != "":
                        atoms.append(last_atm)
                    else:
                        return False, 0, 0, 0, 0, 0, 0, None,''
                elif atom not in rid.atoms:
                    return False,0,0,0,0,0,0,None,''
                else:
                    total_rid += rid.rid;
                    total_ridx += rid.ridx;
                    atm = rid.atoms[atom]
                    other += "_" + str(rid.amino_acid) + str(rid.rid) + str(chain) + "|" + atm.atom_name
                    atoms.append(atm)

        val = 0
        if len(atoms) == 2:
            if not is_nearest and atoms[0].atom_name == atoms[1].atom_name and atoms[0].atom_no == atoms[1].atom_no :#then we actually just want the distance magnitude of the single atom
                val = calc.getMagnitude(atoms[0].x,atoms[0].y,atoms[0].z)
            else:
                val = calc.getDistance(atoms[0].x, atoms[0].y, atoms[0].z,
                                       atoms[1].x, atoms[1].y, atoms[1].z)
        elif len(atoms) == 3:
            val = calc.getAngle(atoms[0].x, atoms[0].y, atoms[0].z,
                                atoms[1].x, atoms[1].y, atoms[1].z,
                                atoms[2].x, atoms[2].y, atoms[2].z)
        elif len(atoms) == 4:
            val = calc.getDihedral(atoms[0].x, atoms[0].y, atoms[0].z,
                                   atoms[1].x, atoms[1].y, atoms[1].z,
                                   atoms[2].x, atoms[2].y, atoms[2].z,
                                   atoms[3].x, atoms[3].y, atoms[3].z)

        total_bfactor = 0
        total_occupancy = 0
        for atm in atoms:
            total_bfactor += atm.bfactor
            total_occupancy  += atm.occupancy



        return True, val, total_bfactor, total_occupancy, total_rid, total_ridx, len(atoms), atoms[0],other[1:]


    def _getNearestAtom(self,refatm, ref_rid,ref_chain,atmtype,resmax,nearest,log=0):
        if log > 1:
            print('LeucipPy(2) nearest:',atmtype,"within res num=",resmax,"nearest=",nearest)
        dic_res = {}
        last_dis = 10000
        last_atom = None
        other = ''
        for chain,resdic in self.chains.items():
            for no,res in resdic.items():
                for attype,atm in res.atoms.items():
                    if attype == atmtype:
                        distance = calc.getDistance(refatm.x, refatm.y,refatm.z,atm.x,atm.y,atm.z)
                        if distance < last_dis and (abs(no - ref_rid) >= resmax or ref_chain != chain):
                            other = str(res.amino_acid) + str(no) + str(chain) + "|" + atm.atom_name
                            dic_res[distance] = [atm,other]

                        #    last_dis = distance
                        #    last_atom = atm


        count = 1
        sorted_dic = dict(sorted(dic_res.items()))
        found = False
        for dis,st in sorted_dic.items():
            if log > 2:
                print("LeucipPy(3) nearest", count,nearest,round(dis,4),st[0].atom_name)
            if int(count) == int(nearest):
                last_dis = dis
                last_atom = st[0]
                other = st[1]
                found = True
                #print("...LeucipPy found")
                break
            count +=1

        if found:
            #print("...LeucipPy return",last_dis,other)
            return last_dis,last_atom,other
        else:
            return 0, "", ""

    def createAtomsList(self):
        self.resolution = self.bio_struc.header['resolution']
        self.pdb_code = self.bio_struc.id
        atomNo = 0
        ridx = 0
        for model in self.bio_struc:
            for chain in model:
                self.chains[chain.id] = {} ####  a chain is a dictionary of residue number to GeoResidue ############################################################
                for residue in chain:
                    aa = residue.get_resname()
                    rid = residue.get_full_id()[3][1]
                    resd = res.GeoResidue(aa,rid,ridx)
                    chain = residue.get_full_id()[2]
                    hetatm = residue.get_full_id()[3][0]
                    ridx = ridx + 1

                    for atom in residue:
                        disordered = 'N'
                        if atom.is_disordered():
                            disordered = 'Y'
                            if atom.disordered_has_id("A"):
                                atom.disordered_select("A")
                        if atom.get_occupancy() < 1:
                            disordered = 'Y'
                        atomNo += 1
                        atom_name = atom.get_name()
                        occupant = atom.get_full_id()[4][1]
                        if occupant == ' ':
                            occupant = 'A'
                        x = atom.get_vector()[0]
                        y = atom.get_vector()[1]
                        z = atom.get_vector()[2]
                        bfactor = atom.get_bfactor()
                        occupancy = atom.get_occupancy()
                        one_atom = atm.GeoAtom(atom_name[0],atom_name,atomNo,disordered,occupancy,bfactor,x,y,z)
                        resd.atoms[atom_name] = one_atom
                    self.chains[chain][rid] = resd




