#useful functions for plotting stats

import sys, numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
   
#colors=['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
colors=[(0,0,0),(230/255.0,159/255.0,0), (86/255.0,180/255.0,233/255.0), (240/255.0,228/255.0,66/255.0), (0,114/255.0,178/255.0), (0,58/255.0,115/255.0), (214/255.0,94/255.0,0), (204/255.0,121/255.0,167/255.0)]

def scatterPlot(x_values, x_label, y_values, y_label, title, outputPrefix, outputFormat):
    fig=plt.figure()
    plot=fig.add_subplot(111)
    plot.loglog(x_values, y_values, 'o')
    plot.set_xlabel(x_label)
    plot.set_ylabel(y_label)
    plot.set_title(title)
    outName=outputPrefix+"."+outputFormat
    plt.savefig(outName, dpi=200, facecolor='w', edgecolor='w', format=outputFormat)
    plt.clf()

def stackedBar(valueLists, y_label, title, x_labels, legend, outputPrefix, outputFormat):
    ind = numpy.arange(len(valueLists[0]))    # the x locations for the groups
    width = 0.5       # the width of the bars: can also be len(x) sequence
    
    colorIndex=0
    plotLayers=[]
    bottoms=[0]*len(valueLists[0])
    bottoms=numpy.array(bottoms)
    for i in range(len(valueLists)):
        #print len(valueLists[i])
        plotLayers.append(plt.bar(ind, valueLists[i], width, bottom=bottoms, color=colors[colorIndex%len(colors)]))
        bottoms=bottoms+numpy.array(valueLists[i])
        colorIndex+=1

    plt.ylabel(y_label)
    plt.title(title)
    plt.xticks(ind+width/2., x_labels, rotation=85 )

    plotZero=[]
    for plotl in plotLayers:
        plotZero.append(plotl[0])
    
    plt.legend(plotZero, legend, loc = 'right', bbox_to_anchor = (1.3, 0.7))
    plt.ylim(0, 1) 
    outName=outputPrefix+"."+outputFormat
    plt.subplots_adjust(bottom=0.38, right=0.8)
    plt.savefig(outName, dpi=200, facecolor='w', edgecolor='w', format=outputFormat)
    plt.clf()

