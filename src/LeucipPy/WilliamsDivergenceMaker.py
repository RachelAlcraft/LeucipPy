# !/usr/bin/python3
"""
WilliamsCoeffcientsMaker (LeucipPy.WilliamsCoeffcientsMaker)
------------------------------------------------------
This class takes multiple measures and a dtaframe and calculates relative information coeffcient by convolution (Mark Williams' method)
"""
import math
import random
import numpy as np
import pandas as pd
import scipy


class Divergence:
    """Class contains data results of divergence relationship
    """
    def __init__(self,geoA,geoB,stat,histAB,diffAB,convAB):
        self.geoA = geoA
        self.geoB = geoB
        self.stat = stat
        self.histA = []#histA
        self.histB = []#histB
        self.histAB = histAB
        self.diffAB = diffAB
        self.convAB = convAB
        self.p_value = float('NaN')
        self.p_hist = {'divergence':[]}
        self.p_mean = float('NaN')
        self.p_std = float('NaN')
        self.p_zvalue = float('NaN')
        self.stat2 = float('NaN')



class WilliamsDivergenceMaker:
    """Class to manage the data for cross correlations

    """

    def __init__(self,data,geos,cut_off=25,density=5,log=0,norm=False,pval_iters=100,delay_load=False):
        """Initialises Williams Coefficients with dataframe, geos and dimensions
        """
        self.data = data
        self.geos = geos
        self.log = log
        self.density = density
        self.samples = len(data.index)
        self.bins = int(math.sqrt(self.samples/self.density))
        self.norm = norm
        self.pval_calcs = pval_iters
        self.cut_off = cut_off
        self.correlations2d = {}
        self.correlations1d = {}
        #self.pvalues = {} #mean, sd, histogram
        # init with all the cross correlation data - currently only for 2d
        if not delay_load:
            for geoA in geos:
                self.calculateCoefficientData1D(geoA,self.cut_off,norm)
                for geoB in geos:
                    if geoA != geoB:
                        self.calculateCoefficientData2D(geoA,geoB,self.cut_off,norm)

    def getCoefficientsDataFrame(self,as_pivot=False,fraction=1,asc=False,filter=0): #rab0
        df = {}
        df['geoA'] = []
        df['geoB'] = []
        df['stat'] = []
        for geogeo,div in self.correlations2d.items():
            df['geoA'].append(div.geoA)
            df['geoB'].append(div.geoB)
            df['stat'].append(div.stat)

        dfdf = pd.DataFrame.from_dict(df)
        if filter > 0:
            query = 'stat >= ' + str(filter)
            dfdf = dfdf.query(query)

        if fraction < 1:
            dfdf = dfdf.sort_values(by='stat',ascending=asc)
            dfdf = dfdf.head(int(len(dfdf.index)*fraction))
        if as_pivot:
            #heat = pd.pivot_table(dfdf, index='geoA', columns='geoB', values=['stat'])
            #cols = heat.index
            #heat = heat.reindex(heat['stat'].sort_values(by=[cols[0]], ascending=False).index)
            #heat.columns = cols
            heat = dfdf.pivot('geoA', 'geoB', 'stat')
            return heat
        else:
            return dfdf

    def getRelativePlotCoefficientsDataFrame(self,geo,as_pivot=False,fraction=1,asc=False,filter=0):
        corrs = {}
        for geogeo, div in self.correlations2d.items():
            if geo == div.geoA:
                corrs[div.geoB] = div
            if geo == div.geoB:
                corrs[div.geoA] = div

        #now I have a dictionary of all the correlations with geo, I want to compare them

        #corrs = self.correlations2d[geo]
        df = {}
        df['geoA'] = []
        df['geoB'] = []
        df['stat'] = []
        # the tuples contain: geoB,stat,histAB,diffAB,convAB
        for geoA,divA in corrs.items():
            for geoB,divB in corrs.items():
                if geoA != geoB:
                    diffA = divA.diffAB
                    diffB = divB.diffAB

                    stat = 0
                    for x in range(0, diffA.shape[0]):
                        for y in range(0, diffA.shape[1]):
                            stat += abs(diffA[x, y] -diffB[x, y])

                    df['geoA'].append(geoA)
                    df['geoB'].append(geoB)
                    df['stat'].append(stat/2)

        dfdf = pd.DataFrame.from_dict(df)
        if filter > 0:
            query = 'stat >= ' + str(filter)
            dfdf = dfdf.query(query)
        if fraction < 1:
            dfdf = dfdf.sort_values(by='stat', ascending=asc)
            dfdf = dfdf.head(int(len(dfdf.index) * fraction))
        if as_pivot:
            #heat = pd.pivot_table(dfdf, index='geoA', columns='geoB', values=['stat'])
            heat = dfdf.pivot('geoA', 'geoB', 'stat')
            return heat
        else:
            return dfdf

    def calculateCoefficientData2D(self,geoA,geoB,trim,norm):
        if geoA+'_'+geoB not in self.correlations2d:# and geoB+'_'+geoA not in self.correlations2d: ineffiecient to do it both ways
            if self.log > 0:
                print('LeucipPy(1): WilliamsDivergenceMaker 2d init for geoA,geoB=', geoA, geoB, 'bins=', self.bins, 'samples=', len(self.data.index), 'pval iters=', self.pval_calcs)
            temp_df = self.data[[geoA,geoB]]
            temp_df = temp_df.sort_values(by=geoA,ascending=True)
            temp_df = temp_df.iloc[trim:,:]
            temp_df = temp_df.sort_values(by=geoA, ascending=False)
            temp_df = temp_df.iloc[trim:, :]
            temp_df = temp_df.sort_values(by=geoB, ascending=True)
            temp_df = temp_df.iloc[trim:, :]
            temp_df = temp_df.sort_values(by=geoB, ascending=False)
            temp_df = temp_df.iloc[trim:, :]
            stat,histAB,diffAB,convAB = self._calculateCorrelation2D(temp_df,geoA,geoB,norm)
            div = Divergence(geoA,geoB,stat,histAB,diffAB,convAB)
            #p_value = self.getPValue(geoA,geoB,stat)
            #self.getPValueInfo(geoA,geoB)
            #div.p
            self.correlations2d[geoA+'_'+geoB] = div
            #now we can calculate the paired values which will go into the same div
            if self.pval_calcs > 0:
                self.calculatePairedPValues(geoA, geoB, self.cut_off, norm)
                div.p_value = self.getPValue(geoA, geoB, stat)
                div.p_zvalue = self.getZValue(geoA,geoB,stat)

#                print('###### testing ########')
                #fifty = self.getCriticalValue(geoA,geoB,0.5)
                #print('mean=',div.p_mean,'val=',fifty)
                #nintyfive = self.getCriticalValue(geoA, geoB, 0.95)
                #prob = self.getPValue(geoA, geoB, nintyfive)
                #print(nintyfive,prob)


    def _calculateCorrelation2D(self,temp_df,geoA,geoB,norm):
        histA, binsA = np.histogram(temp_df[geoA], bins=self.bins, density=False)
        histB, binsB = np.histogram(temp_df[geoB], bins=self.bins, density=False)
        histABtmp,xe,ye = np.histogram2d(temp_df[geoA],temp_df[geoB],bins=self.bins,density=False)

        data_len = len(histA)
        histAB = np.zeros((data_len,data_len))
        convAB = np.zeros((data_len, data_len))
        diffAB = np.zeros((data_len, data_len))

        sum_hist = 0
        sum_conv = 0
        for x in range(0,data_len):
            for y in range(0, data_len):
                histAB[x,y] = histABtmp[y,x] #we just need to reverse this
                convAB[x,y] = histA[y] * histB[x]
                sum_hist += histAB[x,y]
                sum_conv += convAB[x,y]

        #normalise and calc stat
        stat = 0
        for x in range(0,data_len):
            for y in range(0, data_len):
                histAB[x,y] = histAB[x,y]/sum_hist
                convAB[x,y] = convAB[x,y]/sum_conv
                diffAB[x,y] = histAB[x,y]-convAB[x,y]
                stat += abs(diffAB[x,y])

        stat = stat/2 #make it between 0 and 1
        if norm:
            stat = stat * 100/self.bins
        return [stat,histAB,diffAB,convAB]

    def calculatePairedPValues(self,geoA,geoB,trim,norm):
        #if geoA+'_'+geoB not in self.correlations2d and geoB+'_'+geoA not in self.correlations2d:
        temp_df = self.data[[geoA,geoB]]
        temp_df = temp_df.sort_values(by=geoA,ascending=True)
        temp_df = temp_df.iloc[trim:,:]
        temp_df = temp_df.sort_values(by=geoA, ascending=False)
        temp_df = temp_df.iloc[trim:, :]
        temp_df = temp_df.sort_values(by=geoB, ascending=True)
        temp_df = temp_df.iloc[trim:, :]
        temp_df = temp_df.sort_values(by=geoB, ascending=False)
        temp_df = temp_df.iloc[trim:, :]
        hist = []
        for i in range(0,self.pval_calcs):
            rand_df = self.randomiseData(temp_df,[geoA,geoB])
            stat,histAB, diffAB, convAB = self._calculateCorrelation2D(rand_df, geoA, geoB, norm)
            hist.append(stat)
        mean = np.mean(hist)
        sd = np.std(hist)
        hist_dic = {}
        hist_dic['divergence'] = hist
        self.correlations2d[geoA+'_'+geoB].p_mean=mean
        self.correlations2d[geoA+'_'+geoB].p_std=sd
        self.correlations2d[geoA+'_'+geoB].p_hist=hist_dic

    def calculateCoefficientData1D(self,geoA,trim,norm):
        if geoA not in self.correlations1d:
            if self.log > 0:
                print('LeucipPy(1): WilliamsDivergenceMaker 1d init for geo=', geoA, 'bins=', self.bins, 'samples=', len(self.data.index))
            temp_df = self.data[[geoA]]
            temp_df = temp_df.sort_values(by=geoA,ascending=True)
            temp_df = temp_df.iloc[trim:,:]
            temp_df = temp_df.sort_values(by=geoA, ascending=False)
            temp_df = temp_df.iloc[trim:, :]
            histA, binsA = np.histogram(temp_df[geoA], bins=self.bins, density=False)
            data_len = len(histA)

            sum_hist = 0.0
            for x in range(0,data_len):
                sum_hist += float(histA[x])

            #normalise and calc stat
            histAA = np.zeros_like(histA, dtype=float)
            histB = np.zeros_like(histA,dtype=float)
            histDiff = np.zeros_like(histA,dtype=float)
            stat = 0.0
            for x in range(0,data_len):
                histAA[x] = float(histA[x])/float(sum_hist)
                histB[x] = float(1/float(data_len))
                histDiff[x] = float(abs(histAA[x]-histB[x]))
                stat += histDiff[x]

            stat = stat/2 #make it between 0 and 1
            if norm:
                stat = stat * 100/self.bins
            self.correlations1d[geoA] = stat

    def randomiseData(self,data,geos):
        dic_cut= {}
        for geo in geos:
            cut_data = list(data[geo].values)
            random.shuffle(cut_data)
            dic_cut[geo] = cut_data
        df_cut = pd.DataFrame.from_dict(dic_cut)
        return df_cut

    def getCorrelation(self,geos):
        if len(geos) == 1:
            self.calculateCoefficientData1D(geos[0],self.cut_off,self.norm)
            return self.correlations1d[geos[0]]
        elif len(geos) == 2:
            geo1 = geos[0]+'_'+geos[1]
            self.calculateCoefficientData2D(geos[0],geos[1],self.cut_off,self.norm)
            if geo1 in self.correlations2d:
                return self.correlations2d[geo1]

        return 0

    def getZValue(self,geoA,geoB,stat):
        div = self.getCorrelation([geoA, geoB])
        zvalue = (stat-div.p_mean)/div.p_std
        return zvalue

    def getPValue(self,geoA,geoB,stat):
        zvalue = self.getZValue(geoA,geoB,stat)
        #pval = NormalDist().cdf(zvalue) 3.8
        pval = scipy.stats.norm.cdf(abs(zvalue)) #one sided
        #pval = scipy.stats.norm.pdf(abs(zvalue))*2  # two sided
        if pval > 0.5:
            pval = 1-pval
        return pval

    def getCriticalValue(self,geoA,geoB,prob):
        div = self.getCorrelation([geoA, geoB])
        zval = scipy.stats.norm.ppf(abs(prob))
        pval = zval*div.p_std + div.p_mean
        return pval




