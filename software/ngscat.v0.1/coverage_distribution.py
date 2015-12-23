#!/usr/bin/python


import optparse
import sys
import os

import matplotlib
matplotlib.use('Agg')
import xlwt
import numpy

sys.path.append('/home/javi/MGP/utils')

import bam_file


def draw_histogram(distributions, legend, dirout):
    """************************************************************************************************************************************************************
    Task: draws a histogram for a set of distributions.
    Inputs:
        distributions: list of lists. Each sublist contains values of a distribution.
        legend: list of strings, each containing the label to put under each tick in the x axis.
        dirout: string containing the full path to the output directory.
    Output: <dirout>/Coverage_histo.png. It will contain a histogram for the set of distributions.
    ************************************************************************************************************************************************************"""
    
    fig=matplotlib.pyplot.figure(figsize=(13,6))
    ax = fig.add_subplot(111)
    rects = []
    
    # Draw a (overlapped) histogram for each distribution
    bins = numpy.arange(0, 160, 4)
    colors = ["#d7d7d7", "#46a246", "#c65309"]
    for i,dist in enumerate(distributions):    
        rects.append(ax.hist(dist, bins=bins, facecolor=colors[i%3], alpha=0.7)[2])

    #Quitamos los ejes de la derecha y arriba
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    #########################
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Coverage')    
    
    # Shink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Add graphic legend
    ax.legend( tuple([rect[0] for rect in rects]), tuple(legend), loc="upper left", bbox_to_anchor=(1,1) )
    
    fig.savefig(dirout+'/Coverage_histo.png')
    matplotlib.pyplot.close(fig)
    
    
    
    
    
def draw_boxplot(distributions, legend, dirout):
    """************************************************************************************************************************************************************
    Task: draws a box plot for a set of distributions.
    Inputs:
        distributions: list of lists. Each sublist contains values of a distribution.
        legend: list of strings, each containing the label to put under each tick in the x axis.
        dirout: string containing the full path to the output directory.
    Output: <dirout>/Coverage_boxp.png. It will contain a box plot for the set of distributions.
    ************************************************************************************************************************************************************"""
    
    fig2=matplotlib.pyplot.figure(figsize=(13,6))
    ax2 = fig2.add_subplot(111)
    bp = matplotlib.pyplot.boxplot(distributions, notch=0, sym='+', vert=1, whis=1.5)
    matplotlib.pyplot.setp(bp['boxes'], color='black')
    matplotlib.pyplot.setp(bp['whiskers'], color='black')
    matplotlib.pyplot.setp(bp['fliers'], color='red', marker='+')
    ax2.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',alpha=0.5)
    ax2.set_axisbelow(True)
    
    #Quitamos los ejes de la derecha y arriba
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.get_xaxis().tick_bottom()
    ax2.get_yaxis().tick_left()
    #########################
    ax2.set_xticklabels(legend)

    fig2.savefig(dirout+'/Coverage_boxp.png')
    matplotlib.pyplot.close(fig2)




    
def coverage_distribution(bams,beds,dirout,labels,normalize):
    """************************************************************************************************************************************************************
    Task: calculates coverage distribution for a set of bam and bed (capture) files.
    Inputs:
        bams: list of strings with the paths to the bam files.
        beds: list of strings with the paths to the bed files.
        dirout: string containing the full path to the directory where results will be stored.
        labels: list of strings with the labels to name each sample (bam) in the graph.
        normalize: {True,False} to indicate whether normalization should be applied.                
    Output: <dirout>/Coverage_histo.png with the coverage histogram, <dirout>/Coverage_boxp.png with the boxplots and <dirout>/Coverage_stats.xls with quartiles,
        mean, maximum and minimum values.
    ************************************************************************************************************************************************************"""
    
    # Chek output directory exists. In case not, create it
    if(not os.path.isdir(dirout)):
        print 'WARNING: directory '+dirout+' does not exist. Creating...'
        os.mkdir(dirout)
       
    distributions = []
    bamlist = [] 
    # Indexes each bam file and creates the corresponding bam_file objects
    for i,bamfilename in enumerate(bams):
        # Check indexing of the bam file, needed for pysam use
        if(not os.path.isfile(bamfilename+'.bai') and not os.path.isfile(bamfilename.replace('.bam','.bai'))):
            print 'WARNING: index not found for '+bamfilename+'. Indexing...'
            pysam.index(bamfilename)
            print '    Done.'

        bam = bam_file.bam_file(bamfilename, 'rb')
        bamlist.append(bam)
        
    sizes = numpy.array([bam.nreads() for bam in bamlist])
    minsize = sizes.min()
    
    print 'The smaller bam is '+bamlist[sizes.argmin()].filename+' and contains '+str(minsize)+' reads.'
        
    # Process each file and store counting results
    print 'Counting covered bases...'
    for i,bam in enumerate(bamlist):
        print '    '+bam.filename
        
        # Check whether normalization should be applied
        if(normalize): normalizedbam = bam.normalize(minsize)
        else: normalizedbam = bam
            
        distributions.append(normalizedbam.get_coverage_distribution(beds[i]))
                
    draw_histogram(distributions, labels, dirout)
    draw_boxplot(distributions, labels, dirout)
        
    # Initialize the workbook and sheet
    wb = xlwt.Workbook()
    ws = wb.add_sheet('Coverage distribution')

    # Create header font
    header_style = xlwt.easyxf('font: bold on')    

    ws.write(0,0,'Sample',header_style);ws.write(0,1,'Q1',header_style);ws.write(0,2,'Q2',header_style);ws.write(0,3,'Q3',header_style);
    ws.write(0,4,'Max. coverage',header_style); ws.write(0,5,'Min. coverage',header_style);ws.write(0,6,'Mean coverage',header_style);
    
    # Calculate distribution stats for each of the bams
    for i,dist in enumerate(distributions):

        #Sacamos estadisticas
        ndist=numpy.array(dist)
        p25=numpy.percentile(ndist, 25)
        p50=numpy.percentile(ndist, 50)
        p75=numpy.percentile(ndist, 75)

        maximum=numpy.max(ndist)
        minimum=numpy.min(ndist)
        mean=numpy.average(ndist)
    
        ws.write(i+1,0,labels[i]);ws.write(i+1,1,p25);ws.write(i+1,2,p50);ws.write(i+1,3,p75);ws.write(i+1,4,maximum);ws.write(i+1,5,minimum);ws.write(i+1,6,mean);

    wb.save(dirout+'/Coverage_stats.xls')
    
    
def main():   
    ################################################

    #### Options and arguments #####################

    ################################################
    usage="""    
    ************************************************************************************************************************************************************
    Task: calculates coverage distribution for a set of bam and bed (capture) files.
    Output: <dirout>/Coverage_histo.png with the coverage histogram, <dirout>/Coverage_boxp.png with the boxplots and <dirout>/Coverage_stats.xls with quartiles,
        mean, maximum and minimum values.
    ************************************************************************************************************************************************************


    
    Usage: %prog --bams <filenames> --beds <filenames> --dirout <dirname> --labels <listlabels>
    """    
    
    parser = optparse.OptionParser(usage)
    parser.add_option("--bams", dest="bams", help="""String containing a comma separated list of paths to the bam files.""")
    parser.add_option("--beds", dest="beds", help="""String containing a comma separated list of paths to the bed files.""")
    parser.add_option("--dirout", dest="dirout", help="""String containing the full path to the directory where results will be stored.""")
    parser.add_option("--labels", dest="labels", help="""String containing a comma separated list of labels to name each sample (bam) in the graph.""")
    parser.add_option("--normalize", dest="normalize", help="""Optional. {y,n} to indicate whether bam files should be normalized or not.""", default='y')
    

    (options, args) = parser.parse_args()

    # Check number of arguments    
    if len(sys.argv) < 9:
        parser.print_help()
        sys.exit(1)
   
    # call core function
    coverage_distribution(options.bams.split(','), options.beds.split(','), options.dirout, options.labels.split(','), options.normalize=='y')


    print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Finished <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'


if __name__=='__main__':
    main()
