#!/usr/bin/python

import pysam
import sets
import sys
import os
import optparse
import string
import numpy

try:
    import progressbar
except ImportError:
    print 'WARNING: module progressbar was not loaded.'

import xlwt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot 

import bam_file




def draw_graph(fileout, values, xticklabels, xlabel, ylabel, legend=None):
    fig = pyplot.figure(figsize=(13,6))
    ax = fig.add_subplot(111)
    ind = numpy.arange(2*len(values[0]), step=2)
    
    colours = ['#46a246', '#ff0000', '#00ff00', '#0000ff', '#cc0011', '#007722', '#110066', '#c1c1c1', '#544db1', '#aa5198', '#bbd1e9', '#f1c4ab', '#24687a']
    
    width = 0.5       # the width of the bars
    rects = []
    tick_pos = 0
    for i,value_list in enumerate(values):
#        rects.append(ax.bar(ind+tick_pos*width, value_list, width, color='r'))
        rects.append(ax.bar(ind+tick_pos*width, value_list, width, color=colours[i%12]))
        tick_pos += 1
    
    # add some
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
#    ax.set_title('Scores by group and gender')
    ax.set_xticks(ind+(len(values)*width)/2.0)
    ax.set_xticklabels( xticklabels )
    
    ax.set_ylim((0,100))

    if(legend<>None):    
        # Shink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width, box.height*0.9])
    
#            ax.legend( tuple([rect[0] for rect in rects]), tuple(legend), loc="upper left", bbox_to_anchor=(1,1) )
        ax.legend( tuple([rect[0] for rect in rects]), tuple(legend), loc="lower left", bbox_to_anchor=(0,1.03) )
    
#        ax.legend( tuple([rect[0] for rect in rects]), tuple(legend), loc="upper left", bbox_to_anchor=(1,1) )

    fig.savefig(fileout)
    matplotlib.pyplot.close(fig)





def draw_graph_wreplicates(fileout, values_wreplicates, xticklabels, xlabel, ylabel, legend_wreplicates):
    """************************************************************************************************************************************************************
    Task: draws a bar graph with the percentage of covered positions for each sample/replicate.
    Inputs:
        fileout: string containing the full path to the output file.
        values_wreplicates: list of sublists. Each sublist contains the percentage of covered positions at each coverage for a given sample.
        xticklabels: list of strings, each containing the label to put under each tick in the x axis.
        xlabel, ylabel: strings, each containing the label for the x and the y axis respectively.
    Output: fileout will be created. It will contain a bar graph: two bars for each sample, one indicating the number of on-target reads and another one indicating
        the number of off-target reads.
    ************************************************************************************************************************************************************"""
    
    # Merges replicated samples
    values = []
    values_std = []
    legend = []
    for replicate_type in sets.Set(legend_wreplicates):

        # Creates a new list of sublists to store the values of current sample type for each coverage threshold. In other words, each sublist will contain the set of counts
        # for a given coverage threshold
        curr_replicate = [[] for i in range(len(values_wreplicates[0]))]
        nreplicates = 0
        # Passes through each sample looking for those of the same type than "replicate_type"
        for i,another_replicate_type in enumerate(legend_wreplicates):            
            # If sample i is of type replicate_type, add its counts to curr_replicate for posterior aggregation
            if(another_replicate_type==replicate_type):
                # Each sublist contains the counts of all the replicates of this type. There is one sublist for each coverage threshold.
                for j in range(len(curr_replicate)): 
                    curr_replicate[j].append(values_wreplicates[i][j])

        # Calculates mean and std of the counts at each threshold
        new_value_list = map(numpy.mean, curr_replicate)
        new_value_std = map(numpy.std, curr_replicate)        
        values.append(new_value_list)
        values_std.append(new_value_std)
        legend.append(replicate_type)
                                       
   
    fig = pyplot.figure(figsize=(13,6))
    ax = fig.add_subplot(111)
    ind = numpy.arange(2*len(values[0]), step=2)
    
    colours = ['#c65309', '#92b46a', '#b11919', '#166bc5', '#b767a6', '#110066', '#c1c1c1', '#544db1', '#aa5198', '#bbd1e9', '#f1c4ab', '#24687a']
    
    width = 0.15       # the width of the bars
    rects = []
    tick_pos = 0
    for i,value_list in enumerate(values):
#        rects.append(ax.bar(ind+tick_pos*width, value_list, width, color='r'))
        rects.append(ax.bar(ind+tick_pos*width, value_list, width, color=colours[i%12], yerr=values_std[i]))
        tick_pos += 1
    
    # add some
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
#    ax.set_title('Scores by group and gender')
    ax.set_xticks(ind+(len(values)*width)/2.0)
    ax.set_xticklabels( xticklabels )

    # Add legend box
    if(legend<>None):    
        # Shink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
        ax.legend( tuple([rect[0] for rect in rects]), tuple(legend), loc="upper left", bbox_to_anchor=(1,1) )

    fig.savefig(fileout)
    matplotlib.pyplot.close(fig)





def count_lines(filename):
    print 'Calculating file size...'
#    tmp = os.popen('wc -l '+filename)
#    nlines = string.atof(tmp.readline().split(' ')[0])
#    
#    if(tmp.close()<>None):
#        print 'Error: some error occurred while running '
#        print '    wc -l '+filename
#        print 'at bam_file.py'
#        print 'Exiting'
#        sys.exit(1)        
#    print '    Done.'
    
    return len(file(filename).readlines())





def precalculated_target_coverage(coveragefile, coveragelist):
    """************************************************************************************************************************************************************
    Task: draws statistics about the percentage of covered exons and transcripts at different coverage levels. A transcript is considered to be covered when
        at least the 90% of its positions present a coverage greater than the threshold.
    Inputs:
        filelist: list of strings indicating those files to be processed. For a file format example see
            /home/javi/MGP/capture_methods/data/coverage/GU_20120719_FC1_6L1S_AL_01_3376_BC1_AA_F3.filtered.singleHits.realigned.recalibrated.bam.coverage
        coveragelist: list of values with coverage thresholds to use.
        graph_legend: list of descriptions describing each of the files that will be processed. These descriptions will form the legend of the bar plot. 
        dirout: string containing the full path to the directory where data will be saved.
    Output: a summary .xls file and two bar plots depicting coverage vs. %covered-positions and coverage vs. #covered transcripts. Figures will be saved as
        <dirout>/coverage_summary.xls, <dirout>/covered_positions.png and <dirout>/covered_transcripts.png
    ************************************************************************************************************************************************************"""

    ntotal_positions = count_lines(coveragefile)

    # covered_positions_per_depth: list of integers. There will be a position for each coverage threshold. Each value will be the count of positions
    #     covered for the corresponding threshold.
    # ccds_counts: dictionary. Keys are transcript ids. values are lists of two elements. The first element of this list will contain the length of the
    #     transcript. The second element will be a list of integers with as many positions as coverage thresholds, being each value the count of positions
    #     covered for the corresponding threshold. 
    covered_positions_per_depth = [0 for i in range(len(coveragelist))]

    # A progress bar is initialized
    print 'Counting covered positions...'
#    widgets = ['Counting: ', progressbar.Percentage(), ' ', 
#                progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
#    pbar = progressbar.ProgressBar(widgets=widgets, maxval=ntotal_positions).start() 

    # Each line contains the coverage for a given position
    fd = file(coveragefile)
    for k,line in enumerate(fd):
        parts = line.split('\t')            
                    
        # Check whether the coverage is over each threshold
        for i,cov in enumerate(coveragelist):
            current_coverage = string.atof(parts[-1])
            # In case coverage is over the threshold, add 1 to the global number of covered positions and to the counts of the current transcript 
            if(current_coverage>=cov): 
                covered_positions_per_depth[i] += 1                   
        
#        pbar.update(k+1)
        
#    pbar.finish()
    fd.close()
                   
    return [ntotal_positions, covered_positions_per_depth]





def target_coverage_lite(filenames, coveragelist, outprefix, graph_legend=None, xticklabels=None, executiongranted=None, status=None, coveredbases=None, 
                         warnthreshold=90):
    """************************************************************************************************************************************************************
    Task: draws statistics about the percentage of covered exons and transcripts at different coverage levels. A transcript is considered to be covered when
        at least the 90% of its positions present a coverage greater than the threshold.
    Inputs:
        coveragelist: list of values with coverage thresholds to use.
        graph_legend: list of descriptions describing each of the files that will be processed. These descriptions will form the legend of the bar plot.
        dirout: string containing the full path to the directory where data will be saved.
    Output: 
    ************************************************************************************************************************************************************"""

    if(executiongranted<>None):
        executiongranted.acquire()
        
    ntotal_positions_list = []
    covered_positions_per_depth_list = []
    bamlist = []
    
    for filename in filenames:
        ntotal_positions,covered_positions_per_depth = precalculated_target_coverage(filename, coveragelist)
        ntotal_positions_list.append(ntotal_positions)
        covered_positions_per_depth_list.append(covered_positions_per_depth)

    # Initialize the workbook and sheet
    wb = xlwt.Workbook()
    ws = wb.add_sheet('Bases')
    
    # Create header font
    header_style = xlwt.easyxf('font: bold on')
    
    # Write table header
    for i,cov in enumerate(coveragelist):    
        ws.write(0,i*2+1, 'Positions with coverage >='+str(cov)+'x', style=header_style)
        ws.write(0,i*2+1+1, '% of positions with coverage >='+str(cov)+'x', style=header_style)

    # Write one line for each file
    for k,filename in enumerate(filenames):
        ws.write(k+1, 0, filename) #CONTINUAR ARREGLANDO POR AQUI
                                      
        # Write count of covered positions in each file for each coverage threshold    
        # Write counts for current file
        for j,value in enumerate(covered_positions_per_depth_list[k]):
            ws.write(k+1,j*2+1,value)
            ws.write(k+1,j*2+1+1,value*100.0/ntotal_positions_list[k])               
            
        # Calculate percentage of covered positions. Pass through the results of each file.
        # Divide each count by the total number of positions
        for j,value in enumerate(covered_positions_per_depth_list[k]):
            covered_positions_per_depth_list[k][j] = value*100.0/ntotal_positions_list[k]        
      
    # Check whether the output directory is already created
    if(not os.path.isdir(os.path.dirname(outprefix))):
        print 'WARNING: directory '+os.path.dirname(outprefix)+' not found. Creating new directory.'
        os.mkdir(os.path.dirname(outprefix))
          
    if(xticklabels==None): xticklabels = ['>='+str(cov)+'x' for cov in coveragelist]
    
    # Save .xls file and generate the two bar plots.
    wb.save(outprefix+'/coverage_summary.xls')
    draw_graph(outprefix+'covered_positions.png', covered_positions_per_depth_list, xticklabels, 'Coverage threshold', 
               '% covered positions', graph_legend)
    
    print 'Graph saved at '+outprefix+'covered_positions.png'

    if(coveredbases<>None):
        for i,covered_positions_per_depth in enumerate(covered_positions_per_depth_list):
            coveredbases[i] = covered_positions_per_depth[0]
        
    if(executiongranted<>None):
        executiongranted.release()
    
    if(status<>None):
        i=0
        while((i<len(filenames)) and (covered_positions_per_depth_list[i][0] >= warnthreshold)):
            i+=1
        status.value = (i==len(filenames))
        
#target_coverage_lite(['/tmp/test1.coverage', '/tmp/test2.coverage'], [10, 20, 30], '/tmp/')


def target_coverage(filelist, targetfiles, coveragelist, graph_legend, outprefix, xticklabels=None, normalize=False):
    """************************************************************************************************************************************************************
    Task: draws statistics about the percentage of covered exons and transcripts at different coverage levels. A transcript is considered to be covered when
        at least the 90% of its positions present a coverage greater than the threshold.
    Inputs:
        filelist: list of strings indicating those files to be processed. For a file format example see
            /home/javi/MGP/capture_methods/data/coverage/GU_20120719_FC1_6L1S_AL_01_3376_BC1_AA_F3.filtered.singleHits.realigned.recalibrated.bam.coverage
        coveragelist: list of values with coverage thresholds to use.
        graph_legend: list of descriptions describing each of the files that will be processed. These descriptions will form the legend of the bar plot. These 
            labels will also be used to identify sample replicates. Replicates will be merged in one bar in the bar plot.
        outprefix: string containing the full path to the directory where data will be saved.
        xticklabels: list of strings with labels for the ticks in the x axis.
        normalize: boolean to indicate whether bam files should be normalized or not.        
    Output: a summary .xls file and a bar plot depicting coverage vs. %covered-positions. Figures will be saved as
        <dirout>/coverage_summary.xls, <dirout>/covered_positions.png
    ************************************************************************************************************************************************************"""

    numpy.random.seed(1)
    covered_positions = []    
    ntotal_positions = []
    bamlist = []
    
    # Process each file and store counting results
    for filename in filelist:
        # Check whether index already exists for the bam file, needed for pysam use
        if(not os.path.isfile(filename+'.bai')):
            print 'Creating index for '+filename
            pysam.index(filename)
            print '    Done.'
                        
        bamlist.append(bam_file.bam_file(filename))
    sizes = numpy.array([bam.nreads() for bam in bamlist])
    minsize = sizes.min()
    
    print 'The smaller bam is '+filelist[sizes.argmin()]+' and contains '+str(minsize)+' reads.'
        
    # Process each file and store counting results
    print 'Counting covered bases...'
    for i,bam in enumerate(bamlist):
        print '    '+bam.filename
        
        # Check whether normalization should be run
        if(normalize): normalizedbam = bam.normalize(minsize)
        else: normalizedbam = bam
        
        ntotal_positions_tmp,covered_positions_per_depth = normalizedbam.target_coverage(coveragelist, targetfiles[i])
        covered_positions.append(covered_positions_per_depth)
        ntotal_positions.append(ntotal_positions_tmp)

    # Initialize the workbook and sheet
    wb = xlwt.Workbook()
    ws = wb.add_sheet('Bases')
    
    # Create header font
    header_style = xlwt.easyxf('font: bold on')
    
    for i,cov in enumerate(coveragelist):
        ws.write(0,i*2+1, 'Coverage >='+str(cov)+'x', style=header_style)
        ws.write(0,i*2+1+1, '%', style=header_style)
    
    # Write count of covered positions in each file for each coverage threshold    
    for i,value_list in enumerate(covered_positions):
        # Use graph legend elements for row identifiers
        if(graph_legend<>None):
            ws.write(i+1,0, graph_legend[i], style=header_style)
        else:
            ws.write(i+1,0, os.path.basename(filelist[i]), style=header_style)
            
        # Write counts for current file
        for j,value in enumerate(value_list):
            ws.write(i+1,j*2+1,value)
            ws.write(i+1,j*2+1+1,value*100.0/ntotal_positions[i])               
        
    # Calculate percentage of covered positions. Pass through the results of each file.
    for i in range(len(covered_positions)):
        # Divide each count by the total number of positions
        for j,value in enumerate(covered_positions[i]):
            covered_positions[i][j] = value*100.0/ntotal_positions[i]        
      
    # Check whether the output directory is already created
    if(not os.path.isdir(os.path.dirname(outprefix))):
        print 'WARNING: directory '+os.path.dirname(outprefix)+' not found. Creating new directory.'
        os.mkdir(os.path.dirname(outprefix))
          
    # If x labels are not provided, generate ad hoc labels
    if(xticklabels==None): xticklabels = ['>='+str(cov)+'x' for cov in coveragelist]
    
    # Save .xls file and generate the two bar plots.
    wb.save(outprefix+'coverage_summary.xls')
    draw_graph_wreplicates(outprefix+'covered_positions.png', covered_positions, xticklabels, 'Coverage threshold', 
               '% covered positions', graph_legend)
    













def main():   
    ################################################

    #### Options and arguments #####################

    ################################################
    usage="""    
    ************************************************************************************************************************************************************
    Task: draws statistics about the percentage of covered positions.
    Output: a summary .xls file and a bar plot depicting coverage vs. %covered-positions. Figures will be saved as
        <dirout>/coverage_summary.xls, <dirout>/covered_positions.png
    ************************************************************************************************************************************************************


    
    usage: %prog --filelist bamlist --target bedlist --coveragelist coveragelist --graphlegend labellistt --dirout dirname --xticklabels xlabels --normalize {y,n}"""
    
    parser = optparse.OptionParser(usage)
    parser.add_option("--filelist", dest="filelist", help="""String containing a comma separated list indicating those bam files to be processed.""")
    parser.add_option("--target", dest="target", help="""String containing a comma separated list indicating the full path to the bed files containing the targets. Each bam will be compared with the corresponding bed file.""")
    parser.add_option("--coveragelist", dest="coveragelist", help="""String containing a comma separated list of values with coverage thresholds to use.""")
    parser.add_option("--graphlegend", dest="graphlegend", help="""Optional. String containing a comma separated list of descriptions describing each of the files that will be processed. These descriptions will form the legend of the bar plot. These labels will also be used to identify sample replicates. Replicates will be merged in one bar in the bar plot. Default = None.""", default=None)
    parser.add_option("--dirout", dest="dirout", help="""String containing the full path to the directory where data will be saved.""")
    parser.add_option("--xticklabels", dest="xticklabels", help="""Optional. Comma separated strings of labels for the ticks in the x axis. Default: coverage list.""", default=None)
    parser.add_option("--normalize", dest="normalize", help="""Optional. {y,n} to indicate whether bam files should be normalized or not. Default = y.""", default='y')
    (options, args) = parser.parse_args()

    # Check number of arguments    
    if len(sys.argv) < 11:
        parser.print_help()
        sys.exit(1)

    # Check whether graphlegend should be included
    if(options.graphlegend<>None): 
        legend = options.graphlegend.split(',')
    else:
        legend = None
           
    # call core function with/without customized x labels
    if(options.xticklabels<>None):
        target_coverage(options.filelist.split(','), options.target, map(string.atof,options.coveragelist.split(',')), options.graphlegend.split(','), 
                        options.dirout, options.xticklabels.split(','), normalize=(options.normalize=='y'))
    else:
        target_coverage(options.filelist.split(','), options.target.split(','), map(string.atof,options.coveragelist.split(',')), legend, 
                        options.dirout, normalize=(options.normalize=='y'))

    print '>>>>>>>>>>>>>>>>>>>>>> Finished <<<<<<<<<<<<<<<<<<<<<<<<<<'



if __name__=='__main__':
    main()

