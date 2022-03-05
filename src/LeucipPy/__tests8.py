import math
import random

import softwareproperties.ppa

import WilliamsDivergenceMaker as wdm
import BioPythonMaker as bpm
import GeometryMaker as dfm
import HtmlReportMaker as hrm

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

rep = hrm.HtmlReportMaker("WCCC","williams_new.html",cols=7)

############### REAL data #######################
rep.addLineComment('Real pdb data')
strucs = bpm.loadPdbStructures([],'Data/',extension='ent',prefix='pdb')
geo_mak = dfm.GeometryMaker(strucs,log=0)
geos = ['N:CA','CA:C','C:O','C:N+1','C-1:N','N:N+1','N:CA:C:N+1']
data = geo_mak.calculateGeometry(geos)
cm = wdm.WilliamsDivergenceMaker(data,geos,bins=10,log=1,norm=True,p_resample=True,pval_iters=1000)

rep.addLineComment('Scatters')
rep.changeColNumber(6)
print('Creating scatters')
for i in range(0,len(geos)):
    geoA = geos[i]
    for j in range(i+1,len(geos)):
        geoB = geos[j]
        if geoA != geoB:
            print('Scatter',geoA,geoB)

            samp_df, shuffled_df = cm.resampleData(data[[geoA,geoB]])

            stat1, hist2A1, diffAB1, hist2B1 = cm.compareToConvolved(samp_df, geoA, geoB)
            stat2, hist2A2, diffAB2, hist2B2 = cm.compareToConvolved(shuffled_df, geoA, geoB)
            stat3, hist2A3, diffAB3, hist2B3 = cm.compareToConvolved(data, geoA, geoB)

            div = cm.getCorrelation([geoA, geoB])
            stat, p_value, A, D, B,pphist = div.stat,div.p_value,div.histAB,div.diffAB,div.convAB,div.p_hist
            hist = div.p_hist

            maxV = max(np.max(A),np.max(B))
            rep.addLineComment(geoA + ' ' + geoB + ' stat=' + str(round(stat,3)) + ' p_value=' + str(round(p_value,3)))
            rep.addPlot2d(data, 'scatter', geo_x=geoA, geo_y=geoB, hue=geoA,title='observed ' + str(round(stat3,3)))
            rep.addPlot2d(shuffled_df, 'scatter', geo_x=geoA, geo_y=geoB, hue=geoA,title='shuffled ' + str(round(stat2,3)))
            #rep.addPlot2d(samp_df, 'scatter', geo_x=geoA, geo_y=geoB, hue=geoA,title='resampled ' + str(round(stat1,3)))
            if len(hist['divergence_shuffled']) > 0:
                rep.addPlot1d(hist,'histogram',geo_x='divergence_shuffled',title='',overlay=True,alpha=0.5)
                rep.addPlot1d(hist, 'histogram', geo_x='divergence_resampled', title='',alpha=0.5,palette='mediumseagreen')
            rep.addSurface(A,'Original Data',cmin=0,cmax=maxV,palette='Blues')
            rep.addSurface(D, 'Difference Data, '+'stat='+str(round(stat,3)), cmin=-1*maxV, cmax=maxV, palette='RdBu')
            rep.addSurface(B, 'Convolved Data', cmin=0, cmax=maxV, palette='Reds')




rep.printReport()




