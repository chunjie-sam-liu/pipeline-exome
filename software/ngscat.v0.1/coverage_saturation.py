#!/usr/bin/python


import optparse
import string
import os
import sys
import glob

import multiprocessing
import pysam

import simulated_depth
import draw_saturation_curve

sys.path.append('/home/javi/MGP/utils/')
#import job_queue

import bam_file

#TMP = '/data/analisis/MGP/capture_methods_tmp/'
TMP = '/tmp/'




def run(command, launch=True):
    """************************************************************************************************************************************************************
    Task: launches a system call
    Inputs:
        command: string containing the system call.
    ************************************************************************************************************************************************************"""

    print 'CMD: '+command
    
    if(launch):
        # Checks whether an error occurred during the execution of the system call    
        fd = os.popen(command)
        if(fd.close()<>None):
            print 'Some error occurred while executing: '
            print '    '+command





def coverage_saturation_local(bamlist, targets, depthlist, coverage, legend, fileout, executiongranted=None, status=None, slopes=None, tmpdir=None, warnthreshold=1e-5):
    """************************************************************************************************************************************************************
    Task: calculates and draws coverage saturation plots for a list of samples. Just the same as the one below but in multithreading mode.
    Inputs:
        bamlist: list of strings with the names of the bams to process.
        targets: list of strings with the names of the beds containing the targets for each run.
        depthlist: list of integers containing the run depths to test (millions of reads).
        legend: list of descriptions describing each of the files that will be processed. These descriptions will form the legend of the bar plot.
        fileout: String containing the name of the file where the plot will be saved.
    Outputs:       
    ************************************************************************************************************************************************************"""

    # Check whether a temporary directory is provided as an argument
    if(tmpdir<>None):
        TMP = tmpdir
        
    pid = str(os.getpid())
    simulated_depth_processes = []
#    executiongranted = multiprocessing.Semaphore(2)
    
    # Launches one thread for each sample and depth for calculating % of covered positions
    result_files = []
    for i,bam in enumerate(bamlist):
        
        # Check whether there is an index for current bam
        if(not os.path.isfile(bam+'.bai') and not os.path.isfile(bam.replace('.bam','.bai'))):
            print 'WARNING: index not found for '+bam+'. Indexing...'
            pysam.index(bam)
            
        # Threads are launched for each bam and depth point. If provided depth values are greater than the number of reads in the bam file, the maximum depth
        # value to be used will be the number of reads in the bam and no more threads will be launched.
        nreads_bam=bam_file.bam_file(bam).nreads()
        sorteddepths = depthlist
        sorteddepths.sort()
        if(nreads_bam>=(sorteddepths[1]*1000000)):
            endreached=False
            j=0
            while(j<len(depthlist) and not endreached):
                depth = depthlist[j]
                # If a legend is provided, use it to differentiate job ids
                if(legend<>None):
                    jobid = 'coverage_'+pid+'_'+str(depth)+'_'+legend[i].lower()
                else:
                    jobid = 'coverage_'+pid+'_'+str(depth)+'_'+os.path.basename(bamlist[i])
            
                print "Submitting depth "+str(depth)+", file "+bam   
                
                # Activate the flag to indicate that following depth values are greater than the number of reads in the bam
                if((depth*1000000)>=nreads_bam):
                    endreached = True
                
    #            queue.wait()            
                newprocess = multiprocessing.Process(target=simulated_depth.simulated_depth, args=(bam,targets[i],depth,coverage,TMP+'/'+jobid,executiongranted,
                                                                                                   TMP,))
                simulated_depth_processes.append(newprocess)
                newprocess.start()
    #            queue.push(newprocess)            
                
                result_files.append(TMP+'/'+jobid)
                j += 1           
        else:
            print 'WARNING: the number of reads in '+str(bam)+' is '+str(nreads_bam)
            print '    The set of depths provided for coverage saturation calculus is 10e6*'+str(depthlist)
            print '    At least two depths equal or lower than the number of mapped reads are required.' 
            
    if(len(simulated_depth_processes)>0):
        # Wait for all the processess to finish
        for process in simulated_depth_processes:
            process.join()
            process.terminate()
            
        print 'Submitting draw saturation curve...'
        slope_status,tmpslopes = draw_saturation_curve.draw_saturation_curve(result_files,'% covered positions',fileout,legend,warnthreshold=warnthreshold)
        
        if(slopes<>None):
            for i,slope in enumerate(tmpslopes):
                slopes[i] = slope
        
        # Calculate status flag as an OR among the flags for each bam file
        if(status<>None):
            status.value = (sum(slope_status)==len(bamlist))
        
    
        # Remove temporary files
        for afile in result_files: os.remove(afile)
    else:
        status.value = False
    

#coverage_saturation_local(['/tmp/GU_20130213_FC1_2L6S_4L4S_JA_MAM1_F3.filtered.singleHits.realigned.recalibrated.bam',
#                           '/tmp/GU_20130213_FC1_2L6S_4L4S_JA_MAM2_F3.filtered.singleHits.realigned.recalibrated.bam'],
#                           ['/tmp/120718_HG19_DeafReces_MAMP_EZ.g1k.target.bed','/tmp/120718_HG19_DeafReces_MAMP_EZ.g1k.target.bed'], 
#                           depthlist=[0.1,0.2,0.5,1,2,3,4,5,6,7],coverage=10,legend=None,fileout='/tmp/test.png')


def coverage_saturation(bamlist, targets, depthlist, coverage, legend, fileout, submit=True):
    """************************************************************************************************************************************************************
    Task: calculates and draws coverage saturation plots for a list of samples.
    Inputs:
        dirlist: list of strings with the names of the pipeline output directories.
        targets: list of strings with the names of the beds containing the targets for each run.
        depthlist: list of integers containing the run depths to test (millions of reads).
        legend: list of descriptions describing each of the files that will be processed. These descriptions will form the legend of the bar plot.
        fileout: String containing the name of the file where the plot will be saved.
    Outputs:       
    ************************************************************************************************************************************************************"""

    pid = str(os.getpid())
    
    # Launches one job for each sample and depth for calculating % of covered positions
    result_files = ''
    hold_jid_counts=''
    for i,bam in enumerate(bamlist):        
        for depth in depthlist:
            jobid = 'coverage_'+pid+'_'+str(depth)+'_'+legend[i].lower()
            print "Submitting "+jobid                  
            run("""qsub -N """+jobid+""" -o """+TMP+"""/logs/ -e """+TMP+"""/logs/ -l h_vmem=3.0G -M fjavier@bioinfomgp.org -m beas -R y \
                    do_simulated_depth.sh """+bam+""" """+targets[i]+""" """+str(depth)+""" """+str(coverage)+""" """+TMP+ """/""" +jobid, submit)
            hold_jid_counts += ' -hold_jid '+jobid
            result_files += (TMP+'/'+jobid+',')            
            
    jobid='draw_saturation_'+pid
    print 'Submitting '+jobid    
    run("""qsub -N """+jobid+""" -o """+TMP+"""/logs/ -e """+TMP+'/logs/'+hold_jid_counts+""" -l h_vmem=1.0G -M fjavier@bioinfomgp.org -m beas -R y \
            do_draw_saturation_curve.sh """+result_files[:-1]+' "% covered positions" '+fileout, submit)
    
   
            
    
    
    
    
def main():   
    ################################################

    #### Options and arguments #####################

    ################################################
    usage="""
    ************************************************************************************************************************************************************
    Task: calculates and draws coverage saturation plots for a list of samples.
    Output: a png figure saved at --fileout containing saturation plots for all of the samples indicated in --dirlist.                        
    ************************************************************************************************************************************************************   
    
    Usage: %prog --dirlist <dirlist> --depthlist <depthlist> --legend <legend> --out <fileout> 
    """         

    parser = optparse.OptionParser(usage)
    parser.add_option("--bamlist", dest="bamlist", help="""String containing comma separated names of the bam files.""")
    parser.add_option("--targets", dest="targets", help="""String containing comma separated names of the beds containing the targets for each run.""")
    parser.add_option("--depthlist", dest="depthlist", help="""Comma separated integers containing the run depths to test (millions of reads).""")
    parser.add_option("--coverage", dest="coverage", help="""Integer containing the coverage threshold per position to test.""")
    parser.add_option("--legend", dest="legend", help="""String containing a comma separated list of descriptions describing each of the files that will be processed. These descriptions will form the legend of the bar plot.""")
    parser.add_option("--fileout", dest="fileout", help="""String containing the name of the file where the graph will be saved.""")
    parser.add_option("--multithread", dest="multithread", help="""{y,n} to indicate whether to run as jobs in SGE queue or as threads in the local machine.""")
    parser.add_option("--submit", dest="submit", help="""Optional. {y,n} to indicate whether jobs should be submitted or not. Default = y""", default='y')       

    (options, args) = parser.parse_args()

    # Check number of arguments    
    if len(sys.argv) < 15:
        parser.print_help()
        sys.exit(1)
   

    if(options.multithread=='y'):
        # call core function
        coverage_saturation_local(options.bamlist.split(','), options.targets.split(','), map(string.atoi, options.depthlist.split(',')), string.atoi(options.coverage), 
                            options.legend.split(','), options.fileout)    
    else:
        # call core function
        coverage_saturation(options.bamlist.split(','), options.targets.split(','), map(string.atoi, options.depthlist.split(',')), string.atoi(options.coverage), options.legend.split(','), options.fileout, 
                            options.submit=='y')


    print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Finished <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'



if __name__=='__main__':
    main()
