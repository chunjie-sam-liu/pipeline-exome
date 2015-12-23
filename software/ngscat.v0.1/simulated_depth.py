#!/usr/bin/python


import optparse
import string
import os
import sys
import random
import glob

try:
    import progressbar
except ImportError:
    print 'WARNING: module progressbar was not loaded.'

sys.path.append('/home/javi/MGP/utils')

import bam_file

#TMP = '/data/analisis/MGP/capture_methods_tmp/'
TMP = '/tmp/'


def run(command):
    """************************************************************************************************************************************************************
    Task: launches a system call
    Inputs:
        command: string containing the system call.
    ************************************************************************************************************************************************************"""

    print 'CMD: '+command
    # Checks whether an error occurred during the execution of the system call    
    fd = os.popen(command)
    if(fd.close()<>None):
        print 'Some error occurred while executing: '
        print '    '+command





def count_lines(filename):
    """************************************************************************************************************************************************************
    Task: counts the number of lines in a text file. 
    Inputs:
        filename: string containing the name of the file.
    Output:
        nlines: integer containing the number of lines
    ************************************************************************************************************************************************************"""
    
    print 'Calculating file size...'
    tmp = os.popen('wc -l '+filename)
    nlines = string.atof(tmp.readline().split(' ')[0])
    
    # Check whether any error occurred in the system call
    if(tmp.close()<>None):
        print 'Error: some error occurred while running '
        print '    wc -l '+filename
        print 'at tsv_file.py'
        print 'Exiting'
        sys.exit(1)        
    print '    Done.'
    
    return nlines





def simulated_depth(bam, target, depth, coveragethreshold, fileout, executiongranted=None, tmpdir=None):
    """************************************************************************************************************************************************************
    Task: randomly selects a number of reads from a given bam and calculates target coverage. 
    Inputs:
        pipelinehome: String containing the home dir where pipeline output is stored. E.g.: /data/pipeline_outputs/solid/parana/11847_2012-09-14_bfast_190408/
        target: String containing the full path to the bed file.
        depth: Integer containing the run depth in number of reads (millions).
        fileout: String containing the name of the file where results will be stored.
    Output: generates a text file (fileout) with a tab separated line: <dept>\t<ncovered positions>\t<%covered positions>                     
    ************************************************************************************************************************************************************"""
    
    global TMP
    
    if(tmpdir<>None):
        TMP = tmpdir
        
    if(executiongranted<>None):
        executiongranted.acquire()
        
    pid = str(os.getpid())
        
    bam = bam_file.bam_file(bam, 'rb')
    [positions,coverage,chromosomes,processedbed] = bam.myCoverageBed(target,depth*1000000, tmpdir=TMP)
    
#    totalregions = sum([len(processedbed.chrs[chr]) for chr in processedbed.chrs])
    
    # A progress bar is initialized
    print 'Loading coverage...'
#    widgets = ['Progress: ', progressbar.Percentage(), ' ', 
#                progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
#    pbar = progressbar.ProgressBar(widgets=widgets, maxval=totalregions).start() 

    nregions = 0
    npositions = 0
    ncovered_positions = 0
    for chr in processedbed.chrs:
        positionsidx = chromosomes[chr][0]        
        for i,region in enumerate(processedbed.chrs[chr]):
            npositions += (region[1]-region[0]+1)
            while((positionsidx+1)<=chromosomes[chr][1] and positions[positionsidx+1]<=region[1]):
                if(coverage[positionsidx] >= coveragethreshold): 
                    ncovered_positions += (positions[positionsidx+1]-positions[positionsidx])
                positionsidx += 1
                
            if(coverage[positionsidx] >= coveragethreshold): 
                ncovered_positions += (region[1]-positions[positionsidx]+1)
                
            positionsidx += 1
            
            nregions += 1
#            pbar.update(nregions)
            
#    pbar.finish()
                 
    
    print 'Writing results at '+fileout+' ...'
    fd = file(fileout, 'w')
    fd.write(os.path.basename(bam.filename)+'\n')
    fd.write(str(min(bam.nreads(),depth*1000000))+'\t'+str(ncovered_positions)+'\t'+str(ncovered_positions*100.0/npositions))
    fd.close()
    print '    Done.'

    if(executiongranted<>None):
        executiongranted.release()
    









def simulated_depth_old(tmpsam, target, fastq, depth, coverage, fileout):
    """************************************************************************************************************************************************************
    Task: randomly selects a number of reads from a given bam and calculates target coverage. 
    Inputs:
        pipelinehome: String containing the home dir where pipeline output is stored. E.g.: /data/pipeline_outputs/solid/parana/11847_2012-09-14_bfast_190408/
        target: String containing the full path to the bed file.
        depth: Integer containing the run depth in number of reads (millions).
        fileout: String containing the name of the file where results will be stored.
    Output: generates a text file (fileout) with a tab separated line: <dept>\t<ncovered positions>\t<%covered positions>                     
    ************************************************************************************************************************************************************"""
    
    pid = str(os.getpid())
        
    nreads = count_lines(fastq)/4.0
#    nreads = 2000000
    
    newsam = TMP+'/'+pid+'.sam'
    
#    tmpsam = '/home/javi/tmp/test.sam'
    
    nheader_lines = 0
    fdr = file(tmpsam)
    fdw = file(newsam, 'w')
    nextline = fdr.readline()
    while(nextline[0]=='@'):
        fdw.write(nextline)
        nheader_lines += 1
        nextline = fdr.readline()
        
    nmapped_reads = count_lines(tmpsam)-nheader_lines
    depth = min(depth*1000000, nreads)/1000000.0
    
    probmapped = nmapped_reads*1.0/nreads
    probaccept = depth*1000000.0/nreads
    
    # Includes one read (unmapped reads are not written in the new sam file) until the required depth is reached 
    print 'Generating simulated sam...'
    written_reads = 0
    while(written_reads<(depth*1000000) and nextline<>''):
#    while(written_reads<(depth*500000) and nextline<>''):
        # Select mapped/unmapped read
        if(random.uniform(0,1)<=probmapped):
            # Select whether to include current read or not
            if(random.uniform(0,1)<=probaccept):                
                fdw.write(nextline)
                written_reads += 1
                nextline = fdr.readline() 
            else: 
                nextline = fdr.readline()
        else:
            # Select whether to "include" an unmapped read or not. I is not actually included in the bam file, but the counter is increased indicating that a new unmapped
            # read has been considered
            if(random.uniform(0,1)<=probaccept): 
                written_reads += 1    
    fdw.close()
    fdr.close()    
    print '    Done.'
       
    print 'Generating bam...'
    run('samtools view -bS '+newsam+' -o '+TMP+'/'+pid+'.bam')
    print '    Done.'
    os.remove(newsam)

    print 'Calculating count per position...'
    run("""coverageBed -d -abam """+TMP+'/'+pid+""".bam -b """+target+' > "'+TMP+'/'+pid+'.coverage"')
    print '    Done.'
    
    print 'Counting covered positions...'
    fd = file(TMP+'/'+pid+'.coverage')
    ncovered_positions = 0
    npositions = 0
    # Pass through the coverage file counting those positions with coverage greater or equal than 10x
    for i,line in enumerate(fd):
        npositions += 1
        if(string.atof(line.split('\t')[-1]) >= coverage): ncovered_positions += 1        
#        if(i==1000000): break;
    fd.close()        
    print '    Done.'
    
    os.remove(TMP+'/'+pid+'.bam')
    os.remove(TMP+'/'+pid+'.coverage')
    
    print 'Writing results...'
    fd = file(fileout, 'w')
    fd.write(str(written_reads)+'\t'+str(ncovered_positions)+'\t'+str(ncovered_positions*100.0/npositions))
    fd.close()
    print '    Done.'





def main():   
    ################################################

    #### Options and arguments #####################

    ################################################
    usage="""
    ************************************************************************************************************************************************************
    Task: randomly selects a number of reads from a given bam and calculates target coverage.       
    Output: generates a text file (-o option) with a tab separated line: <dept>\t<ncovered positions>\t<%covered positions>                     
    ************************************************************************************************************************************************************
    
    
    Usage: %prog -b <bamfilename> -t <bedfilename> -d <depth> -c <coverage> -o <fileout>    
    """         

    parser = optparse.OptionParser(usage)
    parser.add_option("-b", dest="bam", help="""String containing the full path to the bam file.""")
    parser.add_option("-t", dest="target", help="""String containing the full path to the bed file.""")  
    parser.add_option("-d", dest="depth", help="""Integer containing the run depth in number of reads (millions).""")
    parser.add_option("-c", dest="coverage", help="""Integer containing the coverage threshold per position.""")
    parser.add_option("-o", dest="out", help="""String containing the name of the file where results will be stored.""")
    (options, args) = parser.parse_args()

    # Check number of arguments    
    if len(sys.argv) < 11:
        parser.print_help()
        sys.exit(1)
   
    # call core function
#    simulated_depth(options.sam, options.target, options.fastq, string.atoi(options.depth), string.atoi(options.coverage), options.out)
    simulated_depth(options.bam, options.target, string.atoi(options.depth), string.atoi(options.coverage), options.out)





if __name__=='__main__':
    main()
