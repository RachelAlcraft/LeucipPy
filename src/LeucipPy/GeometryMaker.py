# !/usr/bin/python3
"""
GeoDataFrame (LeucipPy.GeoDataFrame))
------------------------------------------------------
This class manipulates given biopython structures and creates dataframes oft he geometry
"""

import pandas as pd
try:
    from LeucipPy import GeoPdb as pdb
    from LeucipPy import GeoAtom as atm
    from LeucipPy import GeoResidue as res
except:
    import GeoPdb as pdb
    import GeoAtom as atm
    import GeoResidue as res


class GeometryMaker:
    def __init__(self,biopython_structures,log=0,init_biopython=True,atoms_data=pd.DataFrame()):
        """Initialises a GeoDataFrame with the list of biopython structures

        :param biopython_structures: A list of structures that has been created from biopython
        """

        self.bio_strucs = []
        if init_biopython:
            count = 0
            for bio_struc in biopython_structures:
                count += 1
                if log > 0:
                    print('LeucipPy(1): load pdb atoms',str(count)+'/'+str(len(biopython_structures)),'-', bio_struc.id)
                geopdb = pdb.GeoPdb(bio_struc)
                self.bio_strucs.append(geopdb)
        else:
            print('Then we init from a dataframe of atoms instead to match up with our synthetic geometry format')
            '''
            cols = pdb_code, rid, aa, atom, x, y, z
            '''
            pdb_old = ''
            res_old = -1
            resd = None
            residues = []

            for i in range(0,len(atoms_data.index)):
                pdb_code = atoms_data['pdb_code'].values[i]
                rid = atoms_data['rid'].values[i]
                aa = atoms_data['aa'].values[i]
                atom = atoms_data['atom'].values[i]
                x = atoms_data['x'].values[i]
                y = atoms_data['y'].values[i]
                z = atoms_data['z'].values[i]
                if pdb_old != '' and pdb_code != pdb_old:
                    # then we had a pdb code finish in the last step
                    pro = pdb.GeoPdb(None,False,pdb_code,residues)
                    self.bio_strucs.append(pro)
                    residues = []
                if rid != res_old:
                    resd = res.GeoResidue(aa,rid,i)
                    residues.append(resd)
                res_old = rid
                pdb_old = pdb_code
                one_atom = atm.GeoAtom(atom,atom,i,False,1,0, x, y, z)
                resd.atoms[atom] = one_atom

            pro = pdb.GeoPdb(None, False, pdb_old, residues)
            self.bio_strucs.append(pro)

    def calculateGeometry(self,geos,hues=['pdb_code','resolution','chain','aa-1','aa','aa+1','id','rid','ridx','avgrid','avgridx','bfactor','avgbfactor','occupancy','info'],log=0):
        """Creates the geoemtry from the structures in the class

        :param geos: A list of geometric measures to calculate in the format 2,3 or 4 atoms for distance, angle or dihedral, e.g. 'N:CA', 'N:CA:C', or 'N:CA:C:N+1'
        :param hues:A list of hues hat will associate with the geoemtric values, can be bfactor, amino acid (aa), residue number (rid) etc see docs
        :returns: the pandas dataframe with a r per geoemtric calculation per residue wh columns of geoemtric measures and hues
        """

        vals = []
        used_hues = []
        for geo in geos:
            hues.append('info' + geo)

        count = 1
        for geopdb in self.bio_strucs:
            hue_pdb = geopdb.pdb_code
            if log > 0:
                print('LeucipPy(1) df calc for ' + hue_pdb, count, '/', len(self.bio_strucs))
                count += 1

            hue_res = geopdb.resolution
            ridx = 1
            for chain, res in geopdb.chains.items():
                for rid, resd in res.items():
                    avg_bfactor = 0
                    bfactor = 0
                    avg_occupancy = 0
                    avg_rid = 0
                    avg_ridx = 0
                    num_atoms = 0
                    hue_aa = resd.amino_acid
                    tuplerow = []
                    aam1 = ''
                    aap1 = ''
                    refatom = None
                    if len(resd.atoms) > 0 and 'CA' in resd.atoms:
                        refatom = resd.atoms['CA']

                    all_geos_ok = True
                    all_hues = {}
                    for geo in geos:
                        geo_as_atoms = self.geoToAtoms(geo)
                        ok,val,bfac,avgbfac,occ,rno, rnox,num,refatom,other = geopdb.calculateGeometry(chain,rid,geo_as_atoms,log)
                        avg_bfactor += avgbfac
                        avg_occupancy += occ
                        avg_rid += rno
                        avg_ridx += rnox
                        num_atoms += num
                        if bfactor == 0:
                            bfactor = bfac
                        if ok:
                            tuplerow.append(val)
                            all_hues['info' + geo] = other

                            geosaa = geo.split(':')
                            infosaa = other.split('_')
                            for i in range(0,len(geosaa)):
                                geoaa = geosaa[i]
                                if 'CA+1' == geoaa or 'N+1' == geoaa or 'C+1' == geoaa or 'O+1' == geoaa:
                                    aap1 = infosaa[i][:3]
                                if 'CA-1:' == geoaa or 'N-1' == geoaa or 'C-1' == geoaa or 'O-1' == geoaa:
                                    aam1 = infosaa[i][:3]

                        else:
                            all_geos_ok = False
                    #Append hues
                    if all_geos_ok:#we are only adding complete rows
                        all_hues['pdb_code'] =hue_pdb
                        all_hues['resolution'] =hue_res
                        all_hues['chain'] = chain
                        all_hues['aa'] =hue_aa
                        all_hues['aa+1'] = aap1
                        all_hues['aa-1'] = aam1
                        all_hues['avgrid'] = avg_rid/ num_atoms
                        all_hues['avgridx'] =avg_ridx/ num_atoms
                        all_hues['rid'] = rid
                        all_hues['ridx'] = resd.ridx
                        all_hues['avgbfactor'] =avg_bfactor / num_atoms
                        all_hues['bfactor'] = bfactor
                        all_hues['occupancy'] =avg_occupancy / num_atoms
                        used_hues = []
                        for hue in hues:
                            if hue in all_hues and hue not in used_hues:
                                used_hues.append(hue)
                                tuplerow.append(all_hues[hue])
                        # TODO special optional hue of type aa-1:aa:aa+1
                        vals.append(tuplerow)
                    ridx += 1

        geos2 = []
        for geo in geos:
            geos2.append(geo)

        geos2.extend(used_hues)
        df = pd.DataFrame(vals,columns=geos2)
        return df

    def geoToAtoms(self, geo):
        """Converts a single geo to the residue and atoms reuired
        """

        geo_atoms = []
        geo_split = geo.split(':')
        for g in geo_split:
            ref=0
            if '+' in g:
                g_split = g.split('+')
                g_atom = g_split[0]
                ref = int(g_split[1])
            elif '-' in g:
                g_split = g.split('-')
                g_atom = g_split[0]
                ref = int(g_split[1]) * -1
            else:
                g_atom = g
            geo_atoms.append([g_atom,ref])
        return geo_atoms

    def calculateData(self,hues=['pdb_code','resolution','chain','aa','rid','ridx','atom_no','atom_name','element','bfactor','occupancy','x','y','z'],log=0):
        """Creates the geoemtry from the structures in the class

        :param geos: A list of geometric measures to calculate in the format 2,3 or 4 atoms for distance, angle or dihedral, e.g. 'N:CA', 'N:CA:C', or 'N:CA:C:N+1'
        :param hues:A list of hues hat will associate with the geoemtric values, can be bfactor, amino acid (aa), residue number (rid) etc see docs
        :returns: the pandas dataframe with a r per geoemtric calculation per residue wh columns of geoemtric measures and hues
        """

        vals = []
        count = 1
        used_hues = []
        for geopdb in self.bio_strucs:
            pdb = geopdb.pdb_code
            if log > 0:
                print('LeucipPy(1) df calc for ' + pdb, count,'/',len(self.bio_strucs))
                count += 1
            reso = geopdb.resolution

            for chain, res in geopdb.chains.items():
                for rid, resd in res.items():
                    aa = resd.amino_acid
                    for nm,atm in resd.atoms.items():
                        tuplerow = []
                        all_hues = {}
                        all_hues['pdb_code'] =pdb
                        all_hues['resolution'] =reso
                        all_hues['chain'] = chain
                        all_hues['aa'] =aa
                        all_hues['rid'] = resd.rid
                        all_hues['ridx'] =resd.ridx
                        all_hues['atom_no'] = atm.atom_no
                        all_hues['atom_name'] = atm.atom_name
                        all_hues['element'] = atm.atom_type
                        all_hues['bfactor'] =atm.bfactor
                        all_hues['occupancy'] =atm.occupancy
                        all_hues['x'] = atm.x
                        all_hues['y'] = atm.y
                        all_hues['z'] = atm.z
                        used_hues = []
                        for hue in hues:
                            if hue in all_hues:
                                used_hues.append(hue)
                                tuplerow.append(all_hues[hue])
                        vals.append(tuplerow)
        df = pd.DataFrame(vals,columns=used_hues)
        return df

    def filterDataFrame(self,data, inclusions={},exclusions={}):
        df = data
        for ky,vls in inclusions.items():
            df = df[df[ky].isin(vls)]

        for ky,vls in exclusions.items():
            df = df[~df[ky].isin(vls)]
        return df
