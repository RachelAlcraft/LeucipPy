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
        #self.p_mean = float('NaN')
        #self.p_std = float('NaN')
        #self.p_zvalue = float('NaN')
        #self.stat2 = float('NaN')



class WilliamsDivergenceMaker:
    """Class to manage the data for cross correlations

    """

    def __init__(self,data,geos,cut_off=25,bins=5,density=0,log=0,norm=True,pval_iters=100,delay_load=False,p_resample=True):
        """Initialises Williams Coefficients with dataframe, geos and dimensions
        """
        norm=True#overwrite for now to stop tesing confusion
        p_resample = True

        self.data = data
        self.geos = geos
        self.log = log
        self.density = density
        self.p_resample = p_resample
        self.samples = len(data.index)
        if density == 0:
            self.bins = bins
        else:
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
                self.calculateCoefficientData1D(geoA)
                for geoB in geos:
                    if geoA != geoB:
                        self.calculateCoefficientData2D(geoA,geoB)

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

    def calculateCoefficientData2D(self,geoA,geoB):
        if geoA+'_'+geoB not in self.correlations2d:# and geoB+'_'+geoA not in self.correlations2d: ineffiecient to do it both ways
            if self.log > 0:
                print('LeucipPy(1): WilliamsDivergenceMaker 2d init for geoA,geoB=', geoA, geoB, 'bins=', self.bins, 'samples=', len(self.data.index), 'pval iters=', self.pval_calcs)
            #temp_df = self.data[[geoA,geoB]]
            #temp_df = temp_df.sort_values(by=geoA,ascending=True)
            #temp_df = temp_df.iloc[trim:,:]
            #temp_df = temp_df.sort_values(by=geoA, ascending=False)
            #temp_df = temp_df.iloc[trim:, :]
            #temp_df = temp_df.sort_values(by=geoB, ascending=True)
            #temp_df = temp_df.iloc[trim:, :]
            #temp_df = temp_df.sort_values(by=geoB, ascending=False)
            #temp_df = temp_df.iloc[trim:, :]
            stat,histAB,diffAB,convAB = self.compareToConvolved(self.data,geoA,geoB)
            div = Divergence(geoA,geoB,stat,histAB,diffAB,convAB)
            #p_value = self.getPValue(geoA,geoB,stat)
            #self.getPValueInfo(geoA,geoB)
            #div.p
            self.correlations2d[geoA+'_'+geoB] = div
            #now we can calculate the paired values which will go into the same div
            if self.pval_calcs > 0:
                self.calculatePairedPValues(geoA, geoB)
                #div.p_value = self.getPValue(geoA, geoB, stat)
                #div.p_zvalue = self.getZValue(geoA,geoB,stat)

#                print('###### testing ########')
                #fifty = self.getCriticalValue(geoA,geoB,0.5)
                #print('mean=',div.p_mean,'val=',fifty)
                #nintyfive = self.getCriticalValue(geoA, geoB, 0.95)
                #prob = self.getPValue(geoA, geoB, nintyfive)
                #print(nintyfive,prob)


    def compareToConvolved(self,temp_df,geoA,geoB):
        histA, binsA = np.histogram(temp_df[geoA], bins=self.bins, density=False)
        histB, binsB = np.histogram(temp_df[geoB], bins=self.bins, density=False)
        histAB,xe,ye = np.histogram2d(temp_df[geoA],temp_df[geoB],bins=self.bins,density=False)

        data_len = len(histA)
        convAB = np.zeros((data_len, data_len))

        for x in range(0,data_len):
            for y in range(0, data_len):
                convAB[x,y] = histA[x] * histB[y]

        return self.compareTwoHistograms(histAB,convAB)

    def compareToRandom(self,df_compare, geoA, geoB):
        df_random = self.randomiseData(df_compare[[geoA,geoB]])
        return self.compareTwoDistributions(df_compare,df_random,geoA,geoB)

    def compareToObserved(self,df_compare,geoA,geoB):
        return self.compareTwoDistributions(self.data,df_compare,geoA,geoB)
        in_num = len(df_compare.index)

    def compareTwoDistributions(self,df_compareA, df_compareB,geoA, geoB):

        amin = min(df_compareA[geoA].min(), df_compareB[geoA].min())
        amax = max(df_compareA[geoA].max(), df_compareB[geoA].max())
        bmin = min(df_compareA[geoB].min(), df_compareB[geoB].min())
        bmax = max(df_compareA[geoB].max(), df_compareB[geoB].max())

        hist2Ax, xe, ye = np.histogram2d(df_compareA[geoA], df_compareA[geoB], bins=self.bins, density=False, range=[[amin, amax], [bmin, bmax]])
        hist2Bx, xxe, xye = np.histogram2d(df_compareB[geoA], df_compareB[geoB], bins=self.bins, density=False, range=[[amin, amax], [bmin, bmax]])
        return self.compareTwoHistograms(hist2Ax,hist2Bx)

    def compareTwoHistograms(self,histAx,histBx):
        data_len = len(histAx)
        hist2A = np.zeros((data_len, data_len))
        hist2B = np.zeros((data_len, data_len))
        diffAB = np.zeros((data_len, data_len))

        sum_histA = 0
        sum_histB = 0
        for x in range(0, data_len):
            for y in range(0, data_len):
                hist2A[x, y] = histAx[y, x]  # we just need to reverse this
                hist2B[x, y] = histBx[y, x]  # we just need to reverse this
                sum_histA += hist2A[x, y]
                sum_histB += hist2B[x, y]

        # normalise and calc stat
        stat = 0
        for x in range(0, data_len):
            for y in range(0, data_len):
                hist2A[x, y] = hist2A[x, y] / sum_histA
                hist2B[x, y] = hist2B[x, y] / sum_histB
                diffAB[x, y] = hist2A[x, y] - hist2B[x, y]
                stat += abs(diffAB[x, y])

        stat = stat / 2  # make it between 0 and 1
        if self.norm:
            stat = stat / (1 - (1 / self.bins))
        return [stat, hist2A, diffAB, hist2B]


    def calculatePairedPValues(self,geoA,geoB):
        #if geoA+'_'+geoB not in self.correlations2d and geoB+'_'+geoA not in self.correlations2d:
        #temp_df = self.data[[geoA,geoB]]
        #temp_df = temp_df.sort_values(by=geoA,ascending=True)
        #temp_df = temp_df.iloc[trim:,:]
        #temp_df = temp_df.sort_values(by=geoA, ascending=False)
        #temp_df = temp_df.iloc[trim:, :]
        #temp_df = temp_df.sort_values(by=geoB, ascending=True)
        #temp_df = temp_df.iloc[trim:, :]
        #temp_df = temp_df.sort_values(by=geoB, ascending=False)
        #temp_df = temp_df.iloc[trim:, :]
        hist_resamp = []
        hist_shuffle = []
        for i in range(0,self.pval_calcs):
            #if self.p_resample:
            samp_df, shuffled_df = self.resampleData(self.data[[geoA,geoB]])
            #else:
            #    rand_df = self.randomiseData(temp_df,[geoA, geoB])
            stat_shuff,histAB1, diffAB1, convAB1 = self.compareToConvolved(shuffled_df, geoA, geoB)
            stat_samp,histAB2, diffAB2, convAB2 = self.compareToConvolved(samp_df, geoA, geoB)
            hist_shuffle.append(stat_shuff)
            hist_resamp.append(stat_samp)
        mean_resamp = np.mean(hist_resamp)
        mean_shuffle = np.mean(hist_shuffle)
        sd_resamp = np.std(hist_resamp)
        sd_shuffle = np.std(hist_shuffle)
        hist_dic = {}
        hist_dic['divergence_shuffled'] = hist_shuffle
        hist_dic['divergence_resampled'] = hist_resamp

        v_shuf = self.getCriticalValue(hist_shuffle, 0.99)
        v_samp = self.getCriticalValue(hist_resamp, 0.01)
        inter = (v_shuf + v_samp)/2
        p_shuf = 1-self.getPValue(hist_shuffle,inter)
        p_samp = 1-self.getPValue(hist_resamp, inter)
        #print(inter,v_shuf,v_samp,p_shuf,p_samp)

        p_value = (p_shuf + p_samp)/(2-(p_shuf+p_samp))
        self.correlations2d[geoA+'_'+geoB].p_value = p_value

        #self.correlations2d[geoA+'_'+geoB].p_mean=mean_shuffle
        #self.correlations2d[geoA+'_'+geoB].p_std=sd_shuffle
        self.correlations2d[geoA+'_'+geoB].p_hist=hist_dic

    def calculateCoefficientData1D(self,geoA):
        if geoA not in self.correlations1d:
            if self.log > 0:
                print('LeucipPy(1): WilliamsDivergenceMaker 1d init for geo=', geoA, 'bins=', self.bins, 'samples=', len(self.data.index))
            temp_df = self.data[[geoA]]
            #temp_df = temp_df.sort_values(by=geoA,ascending=True)
            #temp_df = temp_df.iloc[trim:,:]
            #temp_df = temp_df.sort_values(by=geoA, ascending=False)
            #temp_df = temp_df.iloc[trim:, :]
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
            if self.norm:
                stat = stat * 100/self.bins
            self.correlations1d[geoA] = stat

    def randomiseData(self,data):
        dic_cut= {}
        for geo in data.columns:
            cut_data = list(data[geo].values)
            random.shuffle(cut_data)
            dic_cut[geo] = cut_data
        df_cut = pd.DataFrame.from_dict(dic_cut)
        return df_cut

    def resampleData(self,data):
        dataresampled = data.sample(frac=1,replace=True)
        datashuffled = self.randomiseData(dataresampled)
        return dataresampled,datashuffled

    def getCorrelation(self,geos):
        if len(geos) == 1:
            self.calculateCoefficientData1D(geos[0],self.cu)
            return self.correlations1d[geos[0]]
        elif len(geos) == 2:
            geo1 = geos[0]+'_'+geos[1]
            self.calculateCoefficientData2D(geos[0],geos[1])
            if geo1 in self.correlations2d:
                return self.correlations2d[geo1]

        return 0

    def getZValue(self,mean,sd,stat):
        zvalue = (stat-mean)/sd
        return zvalue

    def getPValue(self,hist,stat):
        mean = np.mean(hist)
        sd = np.std(hist)
        zvalue = self.getZValue(mean,sd,stat)
        pval = scipy.stats.norm.cdf(abs(zvalue)) #one sided
        return pval

    def getCriticalValue(self,hist,prob):
        mean = np.mean(hist)
        sd = np.std(hist)
        zval = scipy.stats.norm.ppf(abs(prob))
        pval = zval*sd + mean
        return pval




