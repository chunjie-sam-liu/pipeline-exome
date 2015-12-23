#!/usr/bin/python


import string
import os
import sys
import numpy
import optparse

try:
    import progressbar
except ImportError:
    print 'WARNING: module progressbar was not loaded.'

import math
import xlwt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

sys.path.append('/home/javi/MGP/utils/')

import bam_file


def count_lines(filename):
    print 'Calculating file size...'
    tmp = os.popen('wc -l '+filename)
    nlines = string.atof(tmp.readline().split(' ')[0])
    
    if(tmp.close()<>None):
        print 'Error: some error occurred while running '
        print '    wc -l '+self.filename
        print 'at tsv_file.py'
        print 'Exiting'
        sys.exit(1)        
    print '    Done.'
    
    return nlines





def load_coverage_file(filename):
    """*********************************************************************************************************************************************************
    Task: loads coverage counts from a .coverage file. Only exons of size > 3 are considered.
    Inputs:
        filename: string containing the name of the files to be loaded. Must contain the coverage
            per position per exon. For an example see:
            /home/fjavier/MGP/capture_methods/data/coverage/GU_20120719_FC1_6L1S_AL_01_3376_BC1_AA_F3.filtered.singleHits.realigned.recalibrated.bam.coverage
    Outputs: 
        exon: numpy array with the list of exon identifier as they appear in the .coverage file.
        coverage: numpy array with the list of coverage values for each position in each exon. Contains as many items as the exon array.
        length: list containing the length of each exon.
    *********************************************************************************************************************************************************"""
    
    exon = []
    coverage = []
    length = []
    curr_exon = None
    next_exon_idx = 0
    
    # A progress bar is initialized
    widgets = ['Loading exon coverage: ', progressbar.Percentage(), ' ', 
                progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
    pbar = progressbar.ProgressBar(widgets=widgets, maxval=count_lines(filename)).start() 
    
    # One line for each entry
    fd = file(filename)
    for i,line in enumerate(fd):
        parts = line.split('\t') # Depending on the target file, the number of columns are different they share in common the following columns: 1st) chromosome, 2nd) Exon start position, 3rd) Exon end position and last) coverage in a single base
        exon_length = string.atoi(parts[2])-string.atoi(parts[1]) # parts[2] contains the exon end+1
        # Only exons longer than 3 nts
        if(exon_length > 3):                                          
            if(curr_exon<>(parts[0]+"-"+parts[1]+"-"+parts[2])): # New exon found
                curr_exon = parts[0]+"-"+parts[1]+"-"+parts[2]                  
                next_exon_idx += 1		
                
            exon.append(next_exon_idx) # Exon ID (integer number)
            coverage.append(string.atof(parts[-1])) #Coverage (last column)
            length.append(exon_length)
                
        pbar.update(i+1)        
        
    pbar.finish()
    fd.close()
    
    exon = numpy.array(exon)
    coverage = numpy.array(coverage)
    
    return [exon, coverage, length]





def load_coverage_file_per_exon(filename):
    """*********************************************************************************************************************************************************
    JPFLORIDO
    Task: loads normalized standard deviation of coverage per exon .coverage file. Only exons of size > 3 are considered.
    Inputs:
        filename: string containing the name of the files to be loaded. Must contain the coverage
            per position per exon. For an example see:
            /home/fjavier/MGP/capture_methods/data/coverage/GU_20120719_FC1_6L1S_AL_01_3376_BC1_AA_F3.filtered.singleHits.realigned.recalibrated.bam.coverage
    Outputs: 
        exon: numpy array with the list of exon identifier as they appear in the .coverage file.
        coverage: numpy array with the list of coverage values for each position in each exon. Contains as many items as the exon array.
        length: list containing the length of each exon.
    *********************************************************************************************************************************************************"""
    

    coverage = []
    curr_exon = None
    next_exon_idx = 0
    
    # A progress bar is initialized
    print 'Loading exon coverage...'
#    widgets = ['Loading exon coverage: ', progressbar.Percentage(), ' ', 
#                progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
#    pbar = progressbar.ProgressBar(widgets=widgets, maxval=count_lines(filename)).start() 
    
    # One line for each entry
    dictCoverage={} # Dictionary: key (exon id), value (std value of coverage for that exon)
    
    fd = file(filename)
    for i,line in enumerate(fd):
        parts = line.split('\t') # Depending on the target file, the number of columns are different they share in common the following columns: 1st) chromosome, 2nd) Exon start position, 3rd) Exon end position and last) coverage in a single base
        exon_length = string.atoi(parts[2])-string.atoi(parts[1]) # parts[2] contains the exon end+1
        # Only exons longer than 3 nts
        if(exon_length > 3):                                                         
            if(curr_exon<>(parts[0]+"-"+parts[1]+"-"+parts[2])): # New exon found                                  
                curr_exon = parts[0]+"-"+parts[1]+"-"+parts[2]                  
                next_exon_idx += 1   
                if(next_exon_idx>=2): #Measure normalized standard deviation values
                    coverage=numpy.array(coverage)
                    if(len((coverage>0).nonzero()[0])): # Avoid exons with zero coverage values in all bases of exons
                        dictCoverage[next_exon_idx-1]=coverage.std()/float(coverage.mean())                                                                                          
                    coverage=[]                   
            coverage.append(string.atof(parts[-1])) #Coverage (last column)                
#        pbar.update(i+1)                
#    pbar.finish()
    fd.close()
    
    # For the last exon
    coverage=numpy.array(coverage)
    dictCoverage[next_exon_idx]=numpy.array(coverage).std()/float(numpy.array(coverage).mean())        
    
    return dictCoverage




def exon_coverage_std_lite(groups, fileoutprefix, legend=None, executiongranted=None, status=None, coveragestd=None, warnthreshold=0.3):
    """************************************************************************************************************************************************************
        JPFLORIDO
        Task: generates the distribution of coverage standard deviation across exons. Improved version -> coverage file is read exon by exon to improve memory usage
        Inputs:
            groups: list of sublists. Each sublist contains coverage filenames of samples related somehow, e.g. samples sequenced in the same run.    
            fileoutprefix: String containing the fileout prefix.
            legend: list of strings containing descriptions describing each of the groups that will be processed. These descriptions will form the legend of the bar plot.    
            target: target file used   
        Output: two .png figures are generated. One containing the distributions of coverage standard deviation across exons
            and a box plot of such distributions.
        ************************************************************************************************************************************************************"""

    if(executiongranted<>None):
        executiongranted.acquire()
        
    fig = pyplot.figure(figsize=(13,6))
    ax = fig.add_subplot(111)
    boxplot = pyplot.figure()
    axb = boxplot.add_subplot(111)
    
    rects = []
    colours = ['#46a246', '#ff0000', '#00ff00', '#0000ff', '#cc0011', '#007722', '#110066']
    global_stdsampling = []

    # Process each group and draw the corresponding histogram in the graph    
    for colouridx,filelist in enumerate(groups):
        # Samples std for each exon in current file
        for filename in filelist:        
            # Loads exons and coverage per position       
            filename=filename.replace("\n","")
                                                                    
            dictCoverage=load_coverage_file_per_exon(filename)
                   
            std_sampling = []                          
            # Goes through each exon position
            for exon in sorted(dictCoverage.keys()):                            
                std_sampling.append(dictCoverage[exon])               


#        print '# exons < 0.028 - '+legend[colouridx]+': '+str(len((numpy.array(std_sampling)<0.028).nonzero()[0]))
        bins = numpy.arange(0, 1, 0.007)
        rects.append(ax.hist(std_sampling, bins, alpha=0.5, facecolor=colours[colouridx])[2])
        std_sampling = numpy.array(std_sampling)
        global_stdsampling.append(std_sampling)

                                                                    
    # add some
#    fig.suptitle('Distribution of coverage standard deviations (normalized) across exons', fontsize=14, fontweight='bold')
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Normalized standard deviation')
    ax.set_xlim(0, 1)
    
#    boxplot.suptitle('Distribution of coverage standard deviations (normalized) across exons', fontsize=14, fontweight='bold')
    axb.boxplot(global_stdsampling)    

    if(legend<>None and len(legend)>1): 
        if(len(legend[0])>25):
            axb.set_xticklabels([tag[:25]+'...' for tag in legend])
        else:
            axb.set_xticklabels(legend)

        # Shink current axis by 20%       
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width, box.height*0.9])

        # Add graphic legend
        ax.legend( tuple([rect[0] for rect in rects]), tuple(legend), loc="lower left", bbox_to_anchor=(0,1.03) )
        
    
#            ax.legend( tuple([rect[0] for rect in rects]), tuple(legend), loc="upper left", bbox_to_anchor=(1,1) )
        
    fig.savefig(fileoutprefix+'/std_distribution.png')
    matplotlib.pyplot.close(fig)
    boxplot.savefig(fileoutprefix+'/std_boxplot.png')
    matplotlib.pyplot.close(boxplot)
    
    #Imprimimos las estadisticas
    # Initialize the workbook and sheet
    wb = xlwt.Workbook()
    ws = wb.add_sheet('Coverage variation within exons')
    # Create header font
    header_style = xlwt.easyxf('font: bold on')
    ws.write(0,0,'Sample',header_style);
    ws.write(0,1,'Q1',header_style);ws.write(0,2,'Q2',header_style);ws.write(0,3,'Q3',header_style);ws.write(0,4,'Maximum',header_style);
    ws.write(0,5,'Minimum',header_style);ws.write(0,6,'Mean',header_style)

    status.value = True
    for i,std_sampling in enumerate(global_stdsampling):
        std_sampling = numpy.array(std_sampling)
        p25=numpy.percentile(std_sampling, 25)
        p50=numpy.percentile(std_sampling, 50)
        p75=numpy.percentile(std_sampling, 75)
        maximum=numpy.max(std_sampling)
        minimum=numpy.min(std_sampling)
        mean=numpy.average(std_sampling)              

        if(coveragestd<>None):
            coveragestd[i] = mean
                    
        if(status<>None):
            status.value = status.value and (mean<=warnthreshold)
             
        if(legend<>None):
            ws.write(i+1,0,legend[i])
        ws.write(i+1,1,p25);ws.write(i+1,2,p50);ws.write(i+1,3,p75);ws.write(i+1,4,maximum);ws.write(i+1,5,minimum);ws.write(i+1,6,mean);
    wb.save(fileoutprefix+'/std_wexons.xls')
    
    if(executiongranted<>None):
        executiongranted.release()
        
    
 
        

#exon_coverage_std_lite([['/tmp//2338.coverage']], '/tmp/')


def exon_coverage_std(groups, fileoutprefix, bedfilename, legend=None, normalize=True):
    """************************************************************************************************************************************************************
    Task: generates the distribution of coverage standard deviation across exons.
    Inputs:
        groups: list of sublists. Each sublist contains bam filenames of samples related somehow, e.g. samples sequenced in the same run.    
        fileoutprefix: String containing the fileout prefix.
        bedfilename: string containing the name of the bed with the regions to analyze.
        legend: list of strings containing descriptions describing each of the groups that will be processed. These descriptions will form the legend of the bar plot.    
        normalize: {True, False} to indicate whether bam files should be normalized
    Output: two .png figures are generated. One containing the distributions of coverage standard deviation across exons
        and a box plot of such distributions.
    ************************************************************************************************************************************************************"""
       
    minsize = 1000000000000000
    minbamfilename = None
    bamgroups = []
    # Process each group and draw the corresponding histogram in the graph    
    for colouridx,filelist in enumerate(groups):
        bamlist = [] 
        # Samples std for each exon in current file
        for filename in filelist:                
            # Check indexing of the bam file, needed for pysam use
            if(not os.path.isfile(filename+'.bai') and not os.path.isfile(filename.replace('.bam','.bai'))):
                print 'WARNING: index not found for '+filename+'. Indexing...'
                pysam.index(filename)
                print '    Done.'
    
            bam = bam_file.bam_file(filename, 'rb')
            bamlist.append(bam)
            
            # Find the bam with the minimum number of reads
            if(bam.nreads() < minsize):
                minsize = bam.nreads()
                minbamfilename = bam.filename
                
        bamgroups.append(bamlist)                    
    
    print 'The smaller bam is '+minbamfilename+' and contains '+str(minsize)+' reads.'
        
    fig = pyplot.figure(figsize=(13,6))
    ax = fig.add_subplot(111)
    boxplot = pyplot.figure()
    axb = boxplot.add_subplot(111)
    
    rects = []
    colours = ['#ff0000', '#00ff00', '#0000ff', '#cc0011', '#007722', '#110066']
    global_stdsampling = []

    # Process each group and draw the corresponding histogram in the graph    
    for colouridx,filelist in enumerate(bamgroups):
        # Samples std for each exon in current file
        for bam in filelist:
            print '    '+bam.filename
            
            # Check whether normalization should be applied
            if(normalize): normalizedbam = bam.normalize(minsize)
            else: normalizedbam = bam
                    
            std_sampling = normalizedbam.region_coverage_std(bedfilename)

#        print '# exons < 0.028 - '+legend[colouridx]+': '+str(len((numpy.array(std_sampling)<0.028).nonzero()[0]))
        bins = numpy.arange(0, 1, 0.007)
        rects.append(ax.hist(std_sampling, bins, alpha=0.5, facecolor=colours[colouridx])[2])
        std_sampling = numpy.array(std_sampling)
        global_stdsampling.append(list(numpy.log10(std_sampling[(std_sampling>0)])))

                                                                    
    # add some
    fig.suptitle('Distribution of coverage standard deviations (normalized) across exons', fontsize=14, fontweight='bold')
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Normalized standard deviation')
    ax.set_xlim(0, 1)
    
    boxplot.suptitle('Distribution of coverage standard deviations (normalized) across exons', fontsize=14, fontweight='bold')
    axb.boxplot(global_stdsampling)    

    # Check whether graph legend should be included or not
    if(legend<>None): 
        axb.set_xticklabels(legend)                

        # Shink current axis by 20%       
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # Add graphic legend
        ax.legend( tuple([rect[0] for rect in rects]), tuple(legend), loc="upper left", bbox_to_anchor=(1,1) )
        
    fig.savefig(fileoutprefix+'/std_distribution.png')
    matplotlib.pyplot.close(fig)
    boxplot.savefig(fileoutprefix+'/std_boxplot.png')
    matplotlib.pyplot.close(boxplot)
    print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Finished <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
    
    
    
    
def main():   
    ################################################

    #### Options and arguments #####################

    ################################################
    usage="""
    ************************************************************************************************************************************************************
    Task: generates the distribution of coverage standard deviation across exons.
    Output: two .png figures are generated. One containing the distributions of coverage standard deviation across exons
        and a box plot of such distributions.
    ************************************************************************************************************************************************************
    
    
    Usage: %prog --groups <coverage_files_groups> --fileoutprefix <prefix> --bed <filename> --graphlegend <legend>    
    """         

    parser = optparse.OptionParser(usage)
    parser.add_option("--groups", dest="groups", help="""String containing a comma separated list of files, each containing the list of bam files for a given group.""")    
    parser.add_option("--fileoutprefix", dest="fileoutprefix", help="""String containing the fileout prefix.""")
    parser.add_option("--bed", dest="bed", help="""String containing the file name of the bed containing the regions to analyze.""")       
    parser.add_option("--graphlegend", dest="graphlegend", help="""String containing a comma separated list of descriptions describing each of the files that will be processed. These descriptions will form the legend of the bar plot.""")    
    parser.add_option("--normalize", dest="normalize", help="""Optional. {y,n} to indicate whether bam files should be normalized or not.""", default='y')
         
    (options, args) = parser.parse_args()

    # Check number of arguments    
    if len(sys.argv) < 9:
        parser.print_help()
        sys.exit(1)
   
      
    # call core function
    exon_coverage_std([group.split(',') for group in options.groups.split('-')], options.fileoutprefix, options.bed, options.graphlegend.split(','), options.normalize=='y')


if __name__=='__main__':
    main()
