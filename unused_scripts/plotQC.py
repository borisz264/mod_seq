#Created 6/4/2010
#Boris Zinshteyn
#Takes as input any number of folders containing the ouput of ribosomal footprinting pipelines, finds the  QC files based on endings, and plots statistics.
#History:
#   1.0.1  08/28/2013   PV  Modified this to work with pipeline v60 upwards. Changed how the script looks for QC files

import sys
import os

from unused_scripts import plottingFuncs as plot


def plotAnnotations(annotationFiles, directories, outputPrefix):
    plotDict={}
    
    for fileName in annotationFiles:
        fileDict={}
        f=open(fileName, 'r')
        for line in f:
            ll=line.split('\t')
            fileDict[ll[0]]=float(ll[1])
        totalCounts=fileDict['totCounts']
        del fileDict['totCounts']   
        sumCounts=sum(fileDict.values())
        other=(totalCounts)-sumCounts
        print fileName + " : " + str(totalCounts/1000000) + "M"
        for type in fileDict.keys():
            if type not in plotDict:
                plotDict[type]=[]
            plotDict[type].append(fileDict[type]/totalCounts)
        if 'other' not in plotDict:
            plotDict['other']=[]
        plotDict["other"].append(other/totalCounts)
        f.close()
    #print plotDict
    legend=plotDict.keys()
    valueLists=[]
    for label in legend:
        valueLists.append(plotDict[label])
        #print sum(plotDict[label])
    y_label='Fraction of Reads'
    x_labels=[]
    for dir in directories:
        #print "Dire: " + dir
        x_labels.append(dir.split('/')[0])
        
    outputFormat='pdf'
    plot.stackedBar(valueLists, y_label, "Annotation of Mapped Reads for Total RNA Libraries", x_labels, legend, outputPrefix, outputFormat)
    

def plotMappings(mappingFiles, directories, outputPrefix):
    plotDict={}
    for fileName in annotationFiles:
        fileDict={}
        f=open(fileName, 'r')
        for line in f:
            ll=line.split('\t')
            fileDict[ll[0]]=float(ll[1])
        totalCounts=fileDict['totCounts']
        del fileDict['totCounts']   
        sumCounts=sum(fileDict.values())
        other=(totalCounts)-sumCounts
        
        for type in fileDict.keys():
            if type not in plotDict:
                plotDict[type]=[]
            plotDict[type].append(fileDict[type]/totalCounts)
        if 'other' not in plotDict:
            plotDict['other']=[]
        plotDict["other"].append(other/totalCounts)
        f.close()
    #print plotDict
    legend=plotDict.keys()
    valueLists=[]
    for label in legend:
        valueLists.append(plotDict[label])
        #print sum(plotDict[label])
    y_label='Fraction of Reads'
    x_labels=[]
    for dir in directories:
        x_labels.append(dir.split('/')[-1])
    outputFormat='png'
    plot.stackedBar(valueLists, y_label, "Annotation of Mapped Reads", x_labels, legend, outputPrefix, outputFormat)


def main():
    if len(sys.argv) < 3 :
        print "There aren't enough input arguments to run the script!"
        print "usage:"
        print "\t python plotQC.py <<output_prefix>> <<Annotation Counts Files containing Directiories 1-N>>"
        sys.exit()
    
    outputPrefix = sys.argv[1]
    directories = sys.argv[2:]
        
    annotationFiles = []
    mappingFiles = []
    biasFiles = []
    dirs = []
    
    #annReg = re.compile('\*.annotations')
    #biasReg = re.compile('\w*.3pbias')
    #mapReg = re.compile('\w*.mappingcounts')
    for directory in directories:
        if not directory.endswith("/"):
            directory += "/"
        qcdir = directory + "qc"
        if os.path.isdir(qcdir):
            dirList = os.listdir(qcdir)
            for fileName in dirList:
                if fileName.endswith('annotation_counts'):
                    annotationFiles.append(qcdir+'/'+fileName)
                #if fileName.endswith('.mappingcounts'):
                #    mappingFiles.append(directory+'/'+fileName)
                #if fileName.endswith('.3pbias'):
                #    biasFiles.append(directory+'/'+fileName)
            dirs.append(qcdir)
    plotAnnotations(annotationFiles, dirs, outputPrefix)

main()
