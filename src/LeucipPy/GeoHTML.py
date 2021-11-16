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

    def addPlot2d(self,data,plottype,geo_x,geo_y,hue,title='',palette='viridis'):
        fig,ax = plt.subplots()
        ax.grid(b=True, which='major', color='Gainsboro', linestyle='-')
        ax.set_axisbelow(True)
        if plottype == 'scatter':
            g = ax.scatter(data[geo_x], data[geo_y], c=data[hue], cmap=palette, edgecolor='silver', alpha=0.75, linewidth=0.5, s=20)
            cb = plt.colorbar(g)
            cb.set_label(hue)
        elif plottype == 'seaborn':
            alpha = 0.65
            im = sns.scatterplot(x=geo_x, y=geo_y, hue=hue, data=data, alpha=0.75, legend='brief',palette=palette, edgecolor='silver', linewidth=0.5)
            # https://stackoverflow.com/questions/53437462/how-do-i-remove-an-attribute-from-the-legend-of-a-scatter-plot
            # EXTRACT CURRENT HANDLES AND LABELS
            h, l = ax.get_legend_handles_labels()
            # COLOR LEGEND (FIRST guess at size ITEMS) we don;t want to plot the distanceinc
            huelen = len(data.sort_values(by=hue, ascending=True)[hue].unique()) + 1
            col_lgd = plt.legend(h[:huelen], l[:huelen], loc='upper left', bbox_to_anchor=(1.05, 1), fancybox=True, shadow=True, ncol=1)
            plt.gca().add_artist(col_lgd)
            ax.set_xlabel('')
            ax.set_ylabel('')
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)  # Put the legend out of the figure

        ax.set_xlabel(geo_x)
        ax.set_ylabel(geo_y)
        plt.title(title)
        #Having plotted we now need to get the image data from the plt
        encoded = self.getPlotImage(fig, ax)
        htmlstring = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
        self.HTMLIncrement()
        self.html_string += '<td width=' + str(int(100/self.cols)) + '%>' + htmlstring + '</td>\n'

    def addPlot1d(self, data,plottype, geo_x,hue='',title='',palette='crimson'):
        fig, ax = plt.subplots()
        if plottype == 'histogram':
            plt.hist(data[geo_x], EdgeColor='k', bins=20, color=palette, alpha=1, label='geo_x')

        #hue is used to find the outliers
        if hue != '':
            data = data.sort_values(by=geo_x,ascending=True)
            firstval = data.head(1)[hue].values[0]
            try:
                firstval = str(round(firstval))
            except:
                pass
            firstgeo = round(data.head(1)[geo_x].values[0],3)
            data = data.sort_values(by=geo_x, ascending=False)
            lastval = data.head(1)[hue].values[0]
            try:
                lastval = str(round(lastval))
            except:
                pass
            lastgeo = round(data.head(1)[geo_x].values[0],3)
            title += '\n' + hue + ': ' + str(firstval)+ '=' + str(firstgeo) + ' ' + str(lastval) + '=' + str(lastgeo)

        ax.set_xlabel(geo_x)
        plt.title(title)
        # Having plotted we now need to get the image data from the plt
        encoded = self.getPlotImage(fig, ax)
        htmlstring = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
        self.HTMLIncrement()
        self.html_string += '<td width=' + str(int(100 / self.cols)) + '%>' + htmlstring + '</td>\n'

    def addPlotPi(self, data,geo_x,hue,title='',colors=[]):
        fig, ax = plt.subplots()
        if colors == []:
            plt.pie(data[geo_x],labels=data[hue])
        else:
            plt.pie(data[geo_x], labels=data[hue],colors=colors)
        plt.title(title)
        # Having plotted we now need to get the image data from the plt
        encoded = self.getPlotImage(fig, ax)
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

    def addSurface(self, mtx, title='',palette='inferno',alpha=0.9,overlay=False):
        if overlay:
            if not self.overlay_open:
                self.fig, self.ax = plt.subplots()
                self.overlay_open = True
        else:
            if not self.overlay_open:
                self.fig, self.ax = plt.subplots()
            self.overlay_open = False

        image = plt.imshow(mtx,cmap=palette,interpolation='nearest',origin='lower',aspect='equal',alpha=alpha)
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


