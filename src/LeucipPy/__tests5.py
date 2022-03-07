import math
import random

#import softwareproperties.ppa

import WilliamsDivergenceMaker as wdm
import BioPythonMaker as bpm
import GeometryMaker as dfm
import HtmlReportMaker as hrm

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

rep = hrm.HtmlReportMaker("WCCC","Results/williams.html",cols=6)
################ FALSE DATA #########################
rep.addLineComment('Artificial pdb data')
dic_fake = {}
dic_fake['A'] = []
dic_fake['B'] = []
dic_fake['C'] = []
dic_fake['D'] = []
dic_fake['E'] = []
dic_fake['F'] = []
for i in range(0,500):
    dic_fake['A'].append(i)
    dic_fake['B'].append(5)
    dic_fake['C'].append(random.randint(0,10))
    dic_fake['D'].append(math.sin(i/500))
    dic_fake['E'].append(math.cos(i/500))
    dic_fake['F'].append(np.random.normal(0,1))

fake_geos =['A','B','C','D','E','F']
df_fake = pd.DataFrame.from_dict(dic_fake)
cm_fake = wdm.WilliamsDivergenceMaker(df_fake,fake_geos,bins=10,log=1,delay_load=True,pval_iters=1000)
rep.addLineComment('Histograms')
print('Creating histograms')
for geo in fake_geos:
    stat = cm_fake.getCorrelation([geo])
    rep.addPlot1d(df_fake,'histogram',geo,hue='',bins=50,title=str(round(stat,3)))


rep.addLineComment('Scatters')
rep.changeColNumber(7)
print('Creating scatters')
for i in range(0,len(fake_geos)):
    geoA = fake_geos[i]
    for j in range(i+1,len(fake_geos)):
        geoB = fake_geos[j]
        if geoA != geoB:
            print('Scatter',geoA,geoB)
            df_rand = cm_fake.randomiseData(df_fake[[geoA,geoB]])
            div = cm_fake.getCorrelation([geoA, geoB])
            stat, p_value, A, D, B = div.stat,div.p_value,div.histAB,div.diffAB,div.convAB
            hist = div.p_hist
            maxV = max(np.max(A),np.max(B))
            rep.addLineComment(geoA + ' ' + geoB + ' stat=' + str(round(stat,3)) + ' p-value=' + str(round(p_value,3)))
            rep.addPlot2d(df_fake, 'scatter', title=str(round(stat,3)),geo_x=geoA, geo_y=geoB, hue=geoA)
            rep.addPlot2d(df_rand, 'scatter', geo_x=geoA, geo_y=geoB, hue=geoA)
            if len(hist['divergence_shuffled']) > 0:
                rep.addPlot1d(hist,'histogram',geo_x='divergence_shuffled',title='' ,overlay=True,alpha=0.5)
                rep.addPlot1d(hist, 'histogram', geo_x='divergence_resampled', title='',alpha=0.5,palette='mediumseagreen')
            rep.addSurface(A,'Original Data',cmin=0,cmax=maxV,palette='Blues')
            rep.addSurface(D, 'Difference Data', cmin=-1*maxV, cmax=maxV, palette='RdBu')
            rep.addSurface(B, 'Convolved Data', cmin=0, cmax=maxV, palette='Reds')

############### REAL data #######################
rep.addLineComment('Real pdb data')
strucs = bpm.loadPdbStructures([],'PdbData/',extension='ent',prefix='pdb')
geo_mak = dfm.GeometryMaker(strucs,log=0)
geos = ['N:CA','CA:C','C:O','C:N+1','C-1:N','N:N+1','N:CA:C:N+1']
data = geo_mak.calculateGeometry(geos)
cm = wdm.WilliamsDivergenceMaker(data,geos,bins=10,log=1,norm=False,p_resample=True,pval_iters=1000)
cm.getRelativePlotCoefficientsDataFrame('N:CA')


singlemap = cm.getRelativePlotCoefficientsDataFrame('N:CA')
single = cm.getRelativePlotCoefficientsDataFrame('N:CA',as_pivot=True)
print(singlemap)
print(single)
fig,ax = plt.subplots()
splt = sns.heatmap(single,annot=False,fmt='d',linewidth=.5,ax=ax,cmap='inferno_r')
plt.title('N:CA')
rep.addPlotOnly(fig,ax)


heatmap = cm.getCoefficientsDataFrame(filter=0.1)
heat = cm.getCoefficientsDataFrame(as_pivot=True,fraction=0.5,asc=False)
fig,ax = plt.subplots()
splt = sns.heatmap(heat,annot=False,fmt='d',linewidth=.5,ax=ax,cmap='inferno_r')
plt.title('All')
rep.addPlotOnly(fig,ax)

rep.addLineComment('Histograms')
print('Creating histograms')
for geo in geos:
    stat = cm.correlations1d[geo]
    rep.addPlot1d(data,'histogram',geo,hue='pdb_code',bins=50,title=str(round(stat,3)))

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
            if len(hist['divergence_shuffled']) > 0:
                A_B, B_A = cm.calculateKullbackLeibler(phist['divergence_shuffled'],phist['divergence_resampled'])
                title = 'Kullback-Leibler Divergence=[' + str(round(A_B,4)) + ',' +str(round(B_A,4)) + ']'
                rep.addPlot1d(phist,'histogram',geo_x='divergence_shuffled',title='',overlay=True,alpha=0.5)
                rep.addPlot1d(phist, 'histogram', geo_x='divergence_resampled', title=title,alpha=0.5,palette='mediumseagreen')
            rep.addSurface(A,'Original Data',cmin=0,cmax=maxV,palette='Blues')
            rep.addSurface(D, 'Difference Data', cmin=-1*maxV, cmax=maxV, palette='RdBu')
            rep.addSurface(B, 'Convolved Data', cmin=0, cmax=maxV, palette='Reds')




rep.printReport()




