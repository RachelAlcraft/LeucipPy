# !/usr/bin/python3
"""
GeoHTML (LeucipPy.GeoHTML)
------------------------------------------------------
This class manipulates given dataframes into preformatted html reports using matplotlib.
The class assists in high volume output for visual analysis
"""
import base64
import gc
import io

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import gaussian_kde


class GeoHTML:
    """Class for a single atom

    :data html_string:
    :data outpath:
    :data title:
    """

    def __init__(self, title,outpath, cols=3, remove_strip=False):
        """Initialises a GeoPdb with a biopython structure

        :param biopython_structure: A list of structures that has been created from biopython
        """

        #self.df = dataframe
        self.outpath = outpath
        self.title=title
        self.cols = cols
        self.num_col = 0
        self.html_string = self.getHeaderString(self.title, remove_strip)
        self.table_open = False
        self.overlay_open = False
        self.fig = None
        self.ax = None


    def printReport(self):
        if self.table_open:
            self.html_string += '</table>\n'
            self.table_open = False
        self.html_string += self.getFooterString()
        f = open(self.outpath, "w+")
        f.write(self.html_string)
        f.close()

    def changeColNumber(self,cols):
        self.cols = cols
        if self.table_open:
            self.html_string += '</table>\n'
            self.table_open = False
            self.num_col = 0

    def addLineComment(self,comment):
        if self.table_open:
            self.html_string += '</tr></table>\n'
            self.table_open = False
            self.num_col=0
        self.html_string += '<p>'+comment+'</p>\n'

    def addBoxComment(self,comment):
        self.HTMLIncrement()
        self.html_string += '<td width=' + str(int(100 / self.cols)) + '%>'+comment+'</td>\n'

    def addPlot2d(self,data,plottype,geo_x,geo_y,hue,title='',palette='viridis',overlay=False,alpha=1,xrange=[None,None],yrange=[None,None],crange=[None,None]):
        self.incrementOverlay(overlay)
        self.ax.grid(b=True, which='major', color='Gainsboro', linestyle='-')
        self.ax.set_axisbelow(True)

        count = len(data[geo_x])
        minX,maxX = 0,0
        minY, maxY = 0, 0
        cmin, cmax = 0, 0
        try:
            minX, maxX = min(data[geo_x]), max(data[geo_x])
            minY, maxY = min(data[geo_y]), max(data[geo_y])
            cmin, cmax = min(data[hue]), max(data[hue])
        except:
            pass
        if xrange[0] != None:
            minx=xrange[0]
        if xrange[1] != None:
            maxx=xrange[1]
        if yrange[0] != None:
            miny=yrange[0]
        if yrange[1] != None:
            maxy=yrange[1]
        if crange[0] != None:
            cmin=crange[0]
        if crange[1] != None:
            cmax=crange[1]

        if plottype == 'scatter':
            g = self.ax.scatter(data[geo_x], data[geo_y], c=data[hue], cmap=palette, edgecolor='silver', alpha=alpha, linewidth=0.5, s=20,vmin=cmin,vmax=cmax)
            cb = plt.colorbar(g)
            cb.set_label(hue)
        elif plottype == 'seaborn':

            im = sns.scatterplot(x=geo_x, y=geo_y, hue=hue, data=data, alpha=alpha, legend='brief',palette=palette, edgecolor='silver', linewidth=0.5)
            # https://stackoverflow.com/questions/53437462/how-do-i-remove-an-attribute-from-the-legend-of-a-scatter-plot
            # EXTRACT CURRENT HANDLES AND LABELS
            h, l = self.ax.get_legend_handles_labels()
            # COLOR LEGEND (FIRST guess at size ITEMS) we don;t want to plot the distanceinc
            huelen = len(data.sort_values(by=hue, ascending=True)[hue].unique()) + 1
            col_lgd = plt.legend(h[:huelen], l[:huelen], loc='upper left', bbox_to_anchor=(1.05, 1), fancybox=True, shadow=True, ncol=1)
            plt.gca().add_artist(col_lgd)
            self.ax.set_xlabel('')
            self.ax.set_ylabel('')
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)  # Put the legend out of the figure
        elif plottype == 'hist2d':
            x = data[geo_x]
            y = data[geo_y]
            bins=100
            gridsize=50
            hb = plt.hexbin(x, y, bins=bins, cmap=palette, gridsize=gridsize,extent=[minX,maxX,minY,maxY])
        elif plottype == 'probability':
            contours = 12
            bins = 50
            minX, maxX = min(data[geo_x]),max(data[geo_x])
            minY, maxY = min(data[geo_y]),max(data[geo_y])
            plt.axis([minX, maxX, minY, maxY])
            xgrid = np.linspace(minX, maxX, bins)
            ygrid = np.linspace(minY, maxY, bins)
            xdata = data[geo_x]
            ydata = data[geo_y]
            data = np.vstack([xdata, ydata])
            xgrid = np.linspace(minX, maxX, bins)
            ygrid = np.linspace(minY, maxY, bins)
            Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
            grid_sized = np.vstack([Xgrid.ravel(), Ygrid.ravel()])
            # fit an array of size [Ndim, Nsamples]
            kde = gaussian_kde(data, bw_method=0.1)
            # evaluate on a regular grid
            Z = kde.evaluate(grid_sized)
            zgrid = Z.reshape(Xgrid.shape)
            self.ax.grid(True, which='major', axis='both', linestyle='-', color=(0.5, 0.5, 0.5), alpha=0.1)
            im = plt.pcolormesh(xgrid, ygrid, zgrid, shading='gouraud', cmap=palette, alpha=alpha)
            cs = plt.contour(xgrid, ygrid, zgrid, contours, colors='0.7', linewidths=0.4, alpha=alpha,extent=[minX,maxX,minX,maxY])
            self.ax.set_axisbelow(True)

        self.ax.set_xlabel(geo_x + "\nCount=" + str(count))
        self.ax.set_ylabel(geo_y)
        plt.title(title)
        #Having plotted we now need to get the image data from the plt
        if not overlay:
            encoded = self.getPlotImage(self.fig, self.ax)
            htmlstring = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
            self.HTMLIncrement()
            self.html_string += '<td width=' + str(int(100/self.cols)) + '%>' + htmlstring + '</td>\n'

    def addPlot1d(self, data,plottype, geo_x,hue='',title='',palette='crimson',overlay=False,alpha=1,xrange=[None,None],cumulative=False,density=False,bins=20):
        self.incrementOverlay(overlay)
        if xrange[0] == None:
            xrange[0] = min(data[geo_x])
        if xrange[1] == None:
            xrange[1] = max(data[geo_x])

        if len(geo_x[0]) > 1: #todo evident bug if the x is a 1 lne stingf
            mydata = data[geo_x]
            mydata = mydata.transpose()
        else:
            mydata = data[geo_x]

        if plottype == 'histogram':
            g = plt.hist(mydata, edgecolor='k', bins=bins, color=palette, alpha=alpha, label='geo_x',range=xrange,cumulative=cumulative,density=density)
        elif plottype == 'step':
            g = plt.hist(mydata, edgecolor='k', bins=bins, color=palette, alpha=alpha, label='geo_x',range=xrange,cumulative=cumulative,density=density,histtype='step')
        elif plottype == 'stepfilled':
            g = plt.hist(mydata, edgecolor='k', bins=bins, color=palette, alpha=alpha, label='geo_x',range=xrange,cumulative=cumulative,density=density,histtype='stepfilled')
        elif plottype == 'barstacked':
            g = plt.hist(mydata, edgecolor='k', bins=bins, color=palette, alpha=alpha, label='geo_x',range=xrange,cumulative=cumulative,density=density,histtype='barstacked')
        elif plottype == 'violin':
            g = plt.violinplot(mydata,showmeans=False,showextrema=True,showmedians=True,quantiles=[0.25,0.75])

        #cb = plt.colorbar(g)
        #hue is used to find the outliers
        if hue != '':
            geos = []
            if len(geo_x[0]) > 1:  # todo evident bug if the x is a 1 lne stingf
                geos = geo_x
            else:
                geos = [geo_x]

            for geo in geos:
                data = data.sort_values(by=geo,ascending=True)
                firstval = data.head(1)[hue].values[0]
                try:
                    firstval = str(round(firstval))
                except:
                    pass
                firstgeo = round(data.head(1)[geo].values[0],3)
                data = data.sort_values(by=geo, ascending=False)
                lastval = data.head(1)[hue].values[0]
                try:
                    lastval = str(round(lastval))
                except:
                    pass
                lastgeo = round(data.head(1)[geo].values[0],3)
                title += '\n' + geo + ' ' + hue + ': ' + str(firstval)+ '=' + str(firstgeo) + ' ' + str(lastval) + '=' + str(lastgeo)

        self.ax.set_xlabel(geo_x)
        plt.title(title)
        # Having plotted we now need to get the image data from the plt
        if not overlay:
            encoded = self.getPlotImage(self.fig, self.ax)
            htmlstring = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
            self.HTMLIncrement()
            self.html_string += '<td width=' + str(int(100 / self.cols)) + '%>' + htmlstring + '</td>\n'

    def addPlotPi(self, data,geo_x,hue,title='',colors=[],overlay=False,alpha=1):
        self.incrementOverlay(overlay)
        if colors == []:
            plt.pie(data[geo_x],labels=data[hue])
        else:
            plt.pie(data[geo_x], labels=data[hue],colors=colors,alpha=alpha)
        plt.title(title)
        # Having plotted we now need to get the image data from the plt
        if not overlay:
            encoded = self.getPlotImage(self.fig, self.ax)
            htmlstring = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
            self.HTMLIncrement()
            self.html_string += '<td width=' + str(int(100 / self.cols)) + '%>' + htmlstring + '</td>\n'

    def addSeries(self,data,title='',transpose=False):
        if transpose:
            self.addSeries_Transposed(data,title)
        else:
            self.addSeries_Not(data,title)

    def addSeries_Transposed(self, data, title=''):
        rows = len(data.index)
        idx_names = list(data.index)
        innerhtml = "<table>\n"
        innerhtml += "<tr>\n"
        for r in range(0, rows):
            innerhtml += "<td>" + str(idx_names[r]) + "</td>\n"
        innerhtml += "</tr>\n"
        innerhtml += "<tr>"
        for r in range(0, rows):
            header = idx_names[r]
            innerhtml += "<td>"
            try:#make a guess at how to round it if it is a number
                number = float(data[r])
                strnumber = str(round(number, 1))
                if abs(number) < 10:
                    strnumber = str(round(number, 3))
                elif abs(number) < 100:
                    strnumber = str(round(number, 2))
                elif abs(number) > 1000:
                    strnumber = str(int(round(number, 0)))
                # print(header, number,strnumber)
                innerhtml += strnumber

            except:
                innerhtml += str(data[r])
            innerhtml += "</td>\n"

        innerhtml += "</tr>\n"
        innerhtml += "</table>\n"

        self.HTMLIncrement()
        self.html_string += '<td width=' + str(int(100 / self.cols)) + '%>' + title + '</br>' + innerhtml + '</td>\n'

    def addSeries_Not(self, data, title=''):
        rows = len(data.index)
        idx_names = list(data.index)
        innerhtml = "<table>\n"
        for r in range(0, rows):
            innerhtml += "<tr><td>" + str(idx_names[r]) + "</td>\n"
            header = idx_names[r]
            innerhtml += "<td>"
            try:#make a guess at how to round it if it is a number
                number = float(data[r])
                strnumber = str(round(number, 1))
                if abs(number) < 10:
                    strnumber = str(round(number, 3))
                elif abs(number) < 100:
                    strnumber = str(round(number, 2))
                elif abs(number) > 1000:
                    strnumber = str(int(round(number, 0)))
                # print(header, number,strnumber)
                innerhtml += strnumber
            except:
                innerhtml += str(data[r])
            innerhtml += "</td>\n"
            innerhtml += "</tr>\n"
        innerhtml += "</table>\n"

        self.HTMLIncrement()
        self.html_string += '<td width=' + str(int(100 / self.cols)) + '%>' + title + '</br>' + innerhtml + '</td>\n'

    def mtxCap(self,mtx,cap):
        if cap < 0:
            return self.mtxCapMin(mtx,cap)
        else:
            mtxret = np.zeros(mtx.shape)
            mmin,mmax = np.amin(mtx),np.max(mtx)
            capped = mmin + cap*(mmax-mmin)
            a,b = mtx.shape
            for i in range(a):
                for j in range(b):
                    if mtx[i,j] > capped:
                        mtxret[i, j] = capped
                    else:
                        mtxret[i, j] = mtx[i, j]
            return mtxret

    def mtxCapMin(self,mtx,cap):
        mtxret = np.zeros(mtx.shape)
        mmin,mmax = np.amin(mtx),np.max(mtx)
        capped = mmax + cap*(mmax-mmin)
        a,b = mtx.shape
        for i in range(a):
            for j in range(b):
                if mtx[i,j] < capped:
                    mtxret[i, j] = capped
                else:
                    mtxret[i, j] = mtx[i, j]
        return mtxret


    def incrementOverlay(self,overlay):
        if overlay:
            if not self.overlay_open:
                self.fig, self.ax = plt.subplots()
                self.overlay_open = True
        else:
            if not self.overlay_open:
                self.fig, self.ax = plt.subplots()
            self.overlay_open = False

    def addSurface(self, mtx, title='',palette='inferno',alpha=0.9,overlay=False,cmin=None,cmax=None,cap=0,colourbar=False,centred=False):
        self.incrementOverlay(overlay)

        if cap != 0:
            mtx = self.mtxCap(mtx,cap)

        vmin, vmax = np.amin(mtx), np.max(mtx)
        if cmin != None:
            vmin = cmin
        if cmax != None:
            vmax = cmax

        image = plt.imshow(mtx, cmap=palette, interpolation='nearest', origin='lower', aspect='equal', alpha=alpha,vmin=vmin, vmax=vmax)

        plt.axis('off')
        plt.title(title)

        if colourbar:
            self.fig.colorbar(image, ax=self.ax)
        # Having plotted we now need to get the image data from the plt
        if not overlay:
            encoded = self.getPlotImage(self.fig, self.ax)
            htmlstring = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
            self.HTMLIncrement()
            self.html_string += '<td width=' + str(int(100 / self.cols)) + '%>' + htmlstring + '</td>\n'

    def addContours(self, mtx, title='',style='filled', contourlabel=False,levels=10,centred=False,palette='inferno',alpha=0.9,colourbar=True,overlay=False,cmin=None,cmax=None,cap=0,linewidth=None):
        '''
        :param mtx:
        :param title:
        :param style: Can be filled, lines or both
        :param label: label the contours
        :param levels: can be number or cohosen list
        :param centres: sets 0 as the middle in the colours (use divergent colour map)
        :param palette:
        :param alpha:
        :param overlay:
        :return:
        '''
        self.incrementOverlay(overlay)

        # Need to turn the matrix into something countourable
        if cap != 0:
            mtx = self.mtxCap(mtx,cap)

        vmin,vmax = np.amin(mtx),np.max(mtx)
        if cmin != None:
            vmin = cmin
        if cmax != None:
            vmax = cmax
        elif centred:
            vmax = max(abs(vmin),abs(vmax))
            vmin = -1 * vmax

        #Need the X and Z
        a,b = mtx.shape
        x = np.arange(0,a,1)
        y = np.arange(0, b, 1)
        X,Y = np.meshgrid(x,y)
        Z = mtx

        self.ax.set_aspect('equal')
        self.ax.axes.xaxis.set_ticklabels([])
        self.ax.axes.yaxis.set_ticklabels([])

        if style=='filled':
            cf = self.ax.contourf(X,Y,Z,cmap=palette,levels=levels,alpha=alpha,vmin=vmin,vmax=vmax,extend='both')
        elif style == 'lines':
            cf = self.ax.contour(X, Y, Z, cmap=palette, levels=levels, alpha=alpha,vmin=vmin,vmax=vmax,extend='both',linewidths=linewidth)
        else:
            cf = self.ax.contourf(X, Y, Z, cmap=palette, alpha=alpha,vmin=vmin,vmax=vmax,extend='both')
            cf2 = self.ax.contour(X, Y, Z, colors='black', vmin=vmin,vmax=vmax,extend='both',linewidths=linewidth)

        if contourlabel:
            if style=='both':
                self.ax.clabel(cf2,inline=True,fontsize=8)
            elif style=='lines':
                self.ax.clabel(cf, inline=True, fontsize=8)
            # can't put contours names without lines

        if colourbar:
            self.fig.colorbar(cf,ax=self.ax)

        plt.axis('off')
        plt.title(title)
        # Having plotted we now need to get the image data from the plt
        if not overlay:
            encoded = self.getPlotImage(self.fig, self.ax)
            htmlstring = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
            self.HTMLIncrement()
            self.html_string += '<td width=' + str(int(100 / self.cols)) + '%>' + htmlstring + '</td>\n'

    def addDataFrame(self, data, title=''):
        rows = len(data.index)
        cols = data.columns
        innerhtml = "<table>\n"
        innerhtml += "<tr>\n"
        for col in cols:
            innerhtml += "<td>" + str(col) + "</td>\n"
        innerhtml += "</tr>\n"

        for r in range(0, rows):
            innerhtml += "<tr>"
            for col in cols:
                innerhtml += "<td>"
                collist = data[col].tolist()
                strval = str(collist[r])
                innerhtml += strval
                innerhtml += "</td>\n"

        innerhtml += "</tr>\n"
        innerhtml += "</table>\n"

        self.HTMLIncrement()
        self.html_string += '<td width=' + str(int(100 / self.cols)) + '%>' + title + '</br>' + innerhtml + '</td>\n'

    def HTMLIncrement(self):
        if not self.table_open:
            self.html_string += '<table><tr>\n'
            self.table_open = True
        if self.num_col == self.cols:
            self.html_string += '</tr><tr>\n'
            self.num_col = 0
        self.num_col += 1

    def getPlotImage(self,fig, ax):
        img = io.BytesIO()
        fig.savefig(img, format='png', bbox_inches='tight')
        img.seek(0)
        encoded = base64.b64encode(img.getvalue())
        plt.close('all')
        gc.collect()
        return encoded

    def getHeaderString(self,title, remove_strip):
        html = '<!DOCTYPE html><html lang="en"><head><title>LeucipPy Report</title>\n'
        # html += '<style> body {background-color:SeaShell;} table {table-layout:fixed;display:table;margin:0 auto;}td {border:1px solid RosyBrown;background-color:SeaShell;}</style>'
        # html += '<style> body {background-color:HoneyDew;} table {background-color:HoneyDew;} .innertable td {border:1px solid MistyRose;background-color:MintCream;}</style>'
        html += '<style> body {text-align:center;background-color:seashell ;} img {width:95% }'
        html += 'table {font-size:0.8vw;width:95%;table-layout:fixed;display:table;margin:0 auto;background-color:lightgrey ;}'
        html += ' td {border:1px solid MistyRose;background-color:lavenderblush;}</style>'
        html += '</head>\n'
        html += '<body>'
        if not remove_strip:
            html += '<hr/>'
            html += '<div style="background-color:thistle;padding:10px"> LeucipPy: Protein Geometry Correlations </div>'
        html += '<hr/>'
        html += '<div style="background-color:lightgrey;padding:5px">'
        html += '<h1>' + title + '</h1>\n'
        html += '</div>'
        html += '<hr/>\n'

        return html

    def getFooterString(self):
        html = '<hr/><div style = "background-color:thistle;padding:10px" >'
        html += '<a href = "https://rachelalcraft.github.io/leucippy.html" title = "LeucipPy" target = "_self">LeucipPy</a>'
        html += ' by <a href = "mailto:rachelalcraft@gmail.com">Rachel Alcraft</a>'
        html += ' ~ supervisor <a href = "http://people.cryst.bbk.ac.uk/~ubcg66a/">Mark A Williams</a>'
        html += ' ~ Birkbeck, University of London (2021)'
        html += ' ~ Follow <a href = "https://rachelalcraft.github.io/citing.html"> citation guidance </a> </br></div><hr/>'
        return html


