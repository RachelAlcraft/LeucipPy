# !/usr/bin/python3
"""
GeoAtom (LeucipPy.GeoAtom))
------------------------------------------------------
Internal class to manage a single atom
"""

import numpy as np

class CifFile:
    """Class for a cif file of any type

    """

    def __init__(self, filename):
        self.filename = filename
        self.LoopedItems = {}
        self.NonLoopedItems = {}
        with open(filename,'r') as f:
            lines = f.read().splitlines()

        on_loop = False
        on_item = False
        item_id = ''
        headers = []
        for line in lines:
            if line[0] == '#' and on_item: #Closing an item
                on_item = False
                on_loop = False
                item_id = ''
            elif line[0] == 'l' and not on_item:#We are entering a loop set
                on_item = True
                on_loop = True
                headers = []
            elif line[0] == '_' and not on_loop: # we are starting or continuing a non loop set
                on_item = True
                on_loop = False
                id = self.getIdFromLine(line)
                nam = self.getNameFromLine(line)
                val = self.gettValueFromLine(line)
                self.NonLoopedItems[id] = {nam:val}
            elif line[0] == '_' and on_item and on_loop:#adding headers to a loop set
                id = self.getIdFromLine(line)
                nam = self.getNameFromLine(line).strip()
                item_id = id
                headers.append(nam)
                if id not in self.LoopedItems:
                    self.LoopedItems[id] = {nam:[]}
                else:
                    self.LoopedItems[id][nam] = []
            elif line[0] != '#' and line[0] != '_' and on_item and on_loop:#then we are adding a row of values to a loop
                vals = self.getValuesFromLine(line)
                for i in range(len(headers)):
                    header = headers[i]
                    val = vals[i]
                    self.LoopedItems[item_id][header].append(val)

    def getMinMax(self, id, x):
        col_x = [int(xx) for xx in self.LoopedItems[id][x]]
        return min(col_x),max(col_x)

    def getMatrix(self,id,x,y,v,zz=[-1,'']):
        col_x = [int(xx) for xx in self.LoopedItems[id][x]]
        col_y = [int(yy) for yy in self.LoopedItems[id][y]]
        if zz[0] != -1:
            col_z = [int(zzz) for zzz in self.LoopedItems[id][zz[1]]]
        col_v = [float(vv) for vv in self.LoopedItems[id][v]]
        # x and y should be indices that match, but we will make a matrix based on these.ab002  nDhjRELh
        max_idx = max(max(col_x),max(col_y)) + 1
        min_idx = min(min(col_x), min(col_y))
        ax = int(max_idx - min_idx)
        mtx = np.zeros((ax,ax))
        for i in range(len(col_x)):
            if zz[0] == -1:
                x = col_x[i] - min_idx
                y = col_y[i]- min_idx
                v= col_v[i]
                mtx[x,y] = v
            else:
                zi = col_z[i]
                x = col_x[i] - min_idx
                y = col_y[i] - min_idx
                v = col_v[i]
                if zi == zz[0]:
                    mtx[x, y] = v
        return mtx

    def getCube(self,id,x,y,z,v):
        col_x = [int(xx) for xx in self.LoopedItems[id][x]]
        col_y = [int(yy) for yy in self.LoopedItems[id][y]]
        col_z = [int(zz) for zz in self.LoopedItems[id][z]]
        col_v = [float(vv) for vv in self.LoopedItems[id][v]]
        # x and y should be indices that match, but we will make a matrix based on these.ab002  nDhjRELh
        max_idx = max(max(col_x),max(col_y),max(col_z)) + 1
        min_idx = min(min(col_x),min(col_y),min(col_z))
        ax = int(max_idx - min_idx)
        mtx = np.zeros((ax,ax,ax))
        for i in range(len(col_x)):
            x = col_x[i] - min_idx
            y = col_y[i]- min_idx
            z = col_z[i] - min_idx
            v= col_v[i]
            mtx[x,y,z] = v
        return mtx

    ### Helper functions ############################################################################################################
    def getIdFromLine(self,line):
        vals = line.split('.')
        return vals[0]

    def getNameFromLine(self,line):
        vals = line.split('.')
        rhs = vals[1]
        vals = rhs.split(' ')
        return vals[0]

    def gettValueFromLine(self,line):
        vals = line.split('.')
        rhs = vals[1]
        vals = rhs.split(' ')
        rhs = vals[0]
        quotes = [pos for pos,char in enumerate(rhs) if char == "'"]
        if len(quotes)<2:
            rhs = rhs.strip()
        else:
            rhs = rhs[quotes[0]:quotes[1]]
        return rhs

    def getValuesFromLine(self,line):
        bits = []
        quotes = [pos for pos, char in enumerate(line) if char == "'"]
        #if there are no quotes we put it all into bits otherwise we put the quotes into bits as tuple
        while len(quotes) >=2:
            a = quotes[0]
            b = quotes[0+1]
            pre = line[0:a]
            bits.append([False,pre])
            mid = line[a:b]
            bits.append([True, mid])
            line = line[b:]
            quotes = [pos for pos, char in enumerate(line) if char == "'"]
        bits.append([False,line])

        final_bits = []
        for done, line in bits:
            if done:
                final_bits.append(line)
            else:
                splits = line.split(' ')
                for split in splits:
                    if len(split) > 0:
                        final_bits.append(split)
        return final_bits


