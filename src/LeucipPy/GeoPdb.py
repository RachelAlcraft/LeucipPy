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
    def __init__(self,biopython_structure,init_bioython=True,pdb_code='',residues=[]):
        """Initialises a GeoPdb with a biopython structure

        :param biopython_structure: A list of structures that has been created from biopython
        """
        self.bio_struc = biopython_structure
        self.chains = {}
        if init_bioython:
            self.createAtomsList()
        else:
            self.createAtomsList2(pdb_code,residues)

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

        empty_return = [False, 0, 0, 0, 0, 0, 0, 0, None,'']

        rids = self.chains[chain]
        all_there = True
        atoms = []
        first = True
        for atom,offset in geo_atoms:
            #print(atom,offset)
            #if '(' in atom and ')' in atom:
            #    atom = atom[1:len(atom) - 1]
            this_rid = resno + offset
            if False:#int(this_rid) not in rids:
                print('False')
                #print(this_rid, resno)
                #return False,0,0,0,0,0,0,None,''
            else:
                #rid = rids[this_rid]
                brackets = False
                elements = False
                if '{' in atom and '}' in atom:
                    #print(atom)
                    brackets = True
                elif '(' in atom and ')' in atom:
                    #print(atom)
                    brackets = True
                    elements = True

                if brackets and not first: #then we need to do a distance search first and iut must be the second atom
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
                    #for atmtype in atomlist:
                    dis,this_atm,other2 = self._getNearestAtom(refatm, resno,chain,atomlist,offset,nearest,log,elements)
                    if dis < last_dis:
                        last_atm = this_atm
                        last_dis = dis
                        last_other = other2

                    other += "_" + last_other
                    if last_other != "":
                        atoms.append(last_atm)
                    else:
                        if log > 1:
                            print('...LeucipPy(2) Not found',resno,chain,atomlist,offset,nearest)
                        return empty_return
                elif int(this_rid) not in rids:
                    #print(this_rid, resno)
                    return empty_return

                else:
                    rid = rids[this_rid]
                    if brackets:
                        atom = atom[1:len(atom) - 1]
                        isit, atom = self.elementInList(atom, rid.atoms)
                        if isit:
                            atm = rid.atoms[atom]
                            total_rid += rid.rid;
                            total_ridx += rid.ridx;
                            other += "_" + str(rid.amino_acid) + "|" + str(rid.rid) + str(chain) + "|" + atm.atom_name
                            atoms.append(atm)
                        else:
                            return empty_return
                    else:
                        if atom not in rid.atoms:
                            return empty_return
                        total_rid += rid.rid;
                        total_ridx += rid.ridx;
                        atm = rid.atoms[atom]
                        other += "_" + str(rid.amino_acid)  + "|" + str(rid.rid) + str(chain) + "|" + atm.atom_name
                        atoms.append(atm)

            first = False
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
        bfactor = 0
        total_occupancy = 0
        for atm in atoms:
            if bfactor == 0:
                bfactor = atm.bfactor
            total_bfactor += atm.bfactor
            total_occupancy  += atm.occupancy



        return True, val, bfactor,total_bfactor, total_occupancy, total_rid, total_ridx, len(atoms), atoms[0],other[1:]


    def _getNearestAtom(self,refatm, ref_rid,ref_chain,atmtypes,resmax,nearest,log=0,elements=False):
        if log > 1:
            print('LeucipPy(2) nearest:',atmtypes,"within res num=",resmax,"nearest=",nearest)
        dic_res = {}
        last_dis = 10000
        last_atom = None
        other = ''
        for chain,resdic in self.chains.items():
            for no,res in resdic.items():
                for attype,atm in res.atoms.items():
                    inlist = False
                    if elements:
                        if attype[:1] in atmtypes:
                            inlist=True
                        elif len(attype) > 1:
                            if attype[:2] in atmtypes:
                                inlist = True
                    else:
                        if attype in atmtypes:
                            inlist=True
                    if inlist:
                        distance = calc.getDistance(refatm.x, refatm.y,refatm.z,atm.x,atm.y,atm.z)
                        if (abs(int(no) - int(ref_rid)) >= int(resmax) or ref_chain != chain):
                            other = str(res.amino_acid) + "|" + str(no) + str(chain) + "|" + atm.atom_name
                            dic_res[distance] = [atm,other]
                            if log > 2:
                                print("LeucipPy(3) to dictionary", other)

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
                if log > 3:
                    print("...LeucipPy(4) found")
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


    def createAtomsList2(self,pdb_code,residues):
        self.resolution = 0
        self.pdb_code = pdb_code
        self.chains['A'] = {}
        for res in residues:
            self.chains['A'][res.rid] = res

    def elementInList(self,element, atomlist):
        for atm in atomlist:
            if element in atm:
                return True, atm
        return False, None



