# !/usr/bin/python3
"""
GeoAtom (LeucipPy.LeuciFile))
------------------------------------------------------
LeuciFile is the format from LeuciPlus
BEING_DATAITEM
csv
END_DATAITEM
it allows multiple csvs to be output from a batch processs
"""

import numpy as np
import pandas as pd
import math

class LeuciFile:
    """Class for a leuci file of any type

    """

    def __init__(self, filename):
        self.filename = filename
        with open(filename,'r') as f:
            self.text = f.read()

    def getDataFrame(self,item):
        dataResults = {}
        datalst = []
        ID_start = 'BEGIN_' + item
        startPos = self.text.find(ID_start) + len(ID_start)
        ID_end = 'END_' + item
        endPos = self.text.find(ID_end)
        data = (self.text[startPos:endPos]).strip()
        datalsta = data.split('\n')
        firstRow = True
        for row in datalsta:
            if not firstRow:
                lst = row.split(',')
                datalst.append(lst)
            else:
                headers = row.split(',')
            firstRow = False
        df = pd.DataFrame(datalst, columns=headers)
        return df

    def saveText(self,item, filename):
        ID_start = 'BEGIN_' + item
        startPos = self.text.find(ID_start) + len(ID_start)
        ID_end = 'END_' + item
        endPos = self.text.find(ID_end)
        data = (self.text[startPos:endPos]).strip()
        f = open(filename,"w")
        f.write(data)
        f.close()

    def getDataFrameAsMatrix(self,data, hue):
        vmin = 10000
        vmax = -10000
        real_len = len(data[hue].values)
        sq_len = int(math.sqrt(real_len))
        mtx = data[hue].values.reshape(int(sq_len), int(sq_len))
        numtx = np.zeros((sq_len, sq_len))
        for a in range(sq_len):
            for b in range(sq_len):
                numtx[a, b] = float(mtx[a, b])
                vmin = min(vmin, float(mtx[a, b]))
                vmax = max(vmax, float(mtx[a, b]))

        return numtx, vmin, vmax

    def placeDataFrameOnMatrix(self,data, hue,mtx):
        idxs = data.index
        for idx in idxs:
            i = int(data['i'][idx])
            j = int(data['j'][idx])
            val = float(data[hue][idx])
            mtx[i, j] += val
        return mtx

