import math
import random
import WilliamsDivergenceMaker as wdm
import BioPythonMaker as bpm
import GeometryMaker as dfm
import HtmlReportMaker as hrm
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
############### REAL data #######################
rep = hrm.HtmlReportMaker("WCCC","Results/kullback_leibler.html",cols=6)
rep.addLineComment('Real pdb data')
strucs = bpm.loadPdbStructures([],'PdbData/',extension='ent',prefix='pdb')
geo_mak = dfm.GeometryMaker(strucs,log=0)
geos = ['N:CA','CA:C','C:O','C:N+1','C-1:N','N:N+1','N:CA:C:N+1']
data = geo_mak.calculateGeometry(geos)
cm = wdm.WilliamsDivergenceMaker(data,geos,bins=10,log=1,norm=False,p_resample=True,pval_iters=1000)

rep.addLineComment('Scatters')
rep.changeColNumber(6)
print('Creating scatters')
for i in range(0,len(geos)):
    geoA = geos[i]
    for j in range(i+1,len(geos)):
        geoB = geos[j]
        if geoA != geoB:
            print('Scatter',geoA,geoB)
            df_rand = cm.randomiseData(data[[geoA,geoB]])

            div = cm.getCorrelation([geoA, geoB])
            tpl = cm.compareToObserved(data, geoA, geoB)

            stat, p_value, A, D, B,phist = div.stat,div.p_value,div.histAB,div.diffAB,div.convAB,div.p_hist
            maxV = max(np.max(A),np.max(B))
            rep.addLineComment(geoA + ' ' + geoB + ' stat=' + str(round(stat,3)) + ' p_value=' + str(round(p_value,3)))
            rep.addPlot2d(data, 'scatter', geo_x=geoA, geo_y=geoB, hue=geoA)
            rep.addPlot2d(df_rand, 'scatter', geo_x=geoA, geo_y=geoB, hue=geoA)
            if len(phist['divergence_shuffled']) > 0:
                A_B, B_A = cm.calculateKullbackLeibler(phist['divergence_shuffled'],phist['divergence_resampled'])
                #print(A_B,B_A)
                print(str(round(A_B, 4)),str(round(B_A, 4)))
                title = 'Kullback-Leibler Divergence=[' + str(round(A_B,4)) + ',' +str(round(B_A,4)) + ']'
                rep.addPlot1d(phist,'histogram',geo_x='divergence_shuffled',title='',overlay=True,alpha=0.5)
                rep.addPlot1d(phist, 'histogram', geo_x='divergence_resampled', title=title,alpha=0.5,palette='mediumseagreen')
            rep.addSurface(A,'Original Data',cmin=0,cmax=maxV,palette='Blues')
            rep.addSurface(D, 'Difference Data', cmin=-1*maxV, cmax=maxV, palette='RdBu')
            rep.addSurface(B, 'Convolved Data', cmin=0, cmax=maxV, palette='Reds')




rep.printReport()




