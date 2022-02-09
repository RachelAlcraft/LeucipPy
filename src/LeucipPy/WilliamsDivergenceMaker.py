# !/usr/bin/python3
"""
WilliamsCoeffcientsMaker (LeucipPy.WilliamsCoeffcientsMaker)
------------------------------------------------------
This class takes multiple measures and a dtaframe and calculates relative information coeffcient by convolution (Mark Williams' method)
"""
import random
import numpy as np
import pandas as pd
import scipy


class WilliamsDivergenceMaker:
    """Class to manage the data for cross correlations

    """

    def __init__(self,data,geos,cut_off=25,bins=10,log=0,norm=False,pval_iters=100):
        """Initialises Williams Coefficients with dataframe, geos and dimensions
        """
        self.data = data
        self.geos = geos
        self.log = log
        self.bins = bins
        self.pval_calcs = pval_iters
        self.cut_off = cut_off
        self.correlations2d = {}
        self.correlations1d = {}
        self.pvalues = {} #mean, sd, histogram
        # init with all the cross correlation data - currently only for 2d
        for geoA in geos:
            if log > 0:
                print('LeucipPy(1): WilliamsDivergenceMaker 1d init for geo=',geoA,'bins=',bins,'samples=',len(data.index))
            self.calculateCoefficientData1D(geoA,self.cut_off,norm)
            for geoB in geos:
                if geoA != geoB:
                    if log > 0:
                        print('LeucipPy(1): WilliamsDivergenceMaker 2d init for geoA,geoB=', geoA,geoB, 'bins=', bins, 'samples=', len(data.index), 'pval iters=', pval_iters)
                    self.calculatePairedPValues(geoA, geoB, self.cut_off, norm)
                    self.calculateCoefficientData2D(geoA,geoB,self.cut_off,norm)

    def getCoefficientsDataFrame(self,as_pivot=False,fraction=1,asc=False,filter=0): #rab0
        df = {}
        df['geoA'] = []
        df['geoB'] = []
        df['stat'] = []
        for geogeo,tpl in self.correlations2d.items():
            geos = geogeo.split('_')
            geoA,geoB = geos[0],geos[1]
            stat = tpl[0]
            df['geoA'].append(geoA)
            df['geoB'].append(geoB)
            df['stat'].append(stat)

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
        for geogeo, tpl in self.correlations2d.items():
            geos = geogeo.split('_')
            print(geogeo,geos)
            if geo in geos:
                geoB = geos[0]
                if geoB == geo:
                    geoB = geos[1]
                corrs[geoB] = tpl

        #now I have a dictionary of all the correlations with geo, I want to compare them

        #corrs = self.correlations2d[geo]
        df = {}
        df['geoA'] = []
        df['geoB'] = []
        df['stat'] = []
        # the tuples contain: geoB,stat,histAB,diffAB,convAB
        for geoA,tplA in corrs.items():
            for geoB,tplB in corrs.items():
                if geoA != geoB:
                    diffA = tplA[2]
                    diffB = tplB[2]

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
        if geoA + geoB not in self.correlations2d and geoB + geoA not in self.correlations2d:
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
            p_value = self.getPValue(geoA,geoB,stat)
            self.correlations2d[geoA+'_'+geoB] = [stat,p_value,histAB,diffAB,convAB]

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
        if geoA+geoB not in self.pvalues and geoB + geoA not in self.pvalues:
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
            hist_dic['p_value'] = hist
            self.pvalues[geoA+'_'+geoB] = [mean,sd,hist_dic]

    def calculateCoefficientData1D(self,geoA,trim,norm):
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
            return self.correlations1d[geos[0]]
        elif len(geos) == 2:
            geo1 = geos[0]+'_'+geos[1]
            geo2 = geos[1]+'_'+geos[0]
            if geo1 in self.correlations2d:
                return self.correlations2d[geo1]
            if geo2 in self.correlations2d:
                return self.correlations2d[geo2]
        return 0

    def getPValueInfo(self,geoA,geoB):
        geo1 = geoA+'_'+geoB
        geo2 = geoB+'_'+geoA
        if geo1 in self.pvalues:
            return self.pvalues[geo1]
        if geo2 in self.pvalues:
            return self.pvalues[geo2]
        return 0,1,[]

    def getPValue(self,geoA,geoB,stat):
        mean,sd,hist = self.getPValueInfo(geoA,geoB)
        zvalue = (stat-mean)/sd
        #pval = NormalDist().cdf(zvalue) 3.8
        pval = scipy.stats.norm.pdf(abs(zvalue)) #one sided
        pval = scipy.stats.norm.pdf(abs(zvalue))*2  # two sided#
        if pval > 0.5:
            pval = 1-pval
        return pval




