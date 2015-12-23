#!/usr/bin/python

import sys
import re
import os
import sys
import optparse
import string
import numpy

try:
	import progressbar
except ImportError:
	print 'WARNING: module progressbar was not loaded.'

import glob
import pysam

try:
	import pybedtools
except ImportError:
	print 'WARNING: module pybedtools was not imported'
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
 
### NEW PACKAGES ADDED for lite version of gcbias #####
from pylab import plot, figure, imshow, xlabel, ylabel, cm, show
from scipy import stats, mgrid, c_, reshape, random, rot90
#########################

sys.path.append('/home/javi/MGP/utils')
import bam_file
import bed_file


CCDS_BED = '/home/javi/MGP/capture_methods/lib/CCDS_UCSC_3-09-2012.bed'
REF = '/usr/local/reference_genomes/human/human_g1k_v37.fasta'
TMP = '/tmp/'
BEDTOOLSPATH = '/usr/local/bedtools/bin/'
HOME = '/home/javi/'




def run(command):
	"""************************************************************************************************************************************************************
	Task: launches a system call
	Inputs:
		command: string containing the system call.
	************************************************************************************************************************************************************"""

	# Checks whether an error occurred during the execution of the system call	
	fd = os.popen(command)
	if(fd.close()<>None):
		print 'Some error occurred while executing: '
		print '	'+command





def count_lines(filename):
	print 'Calculating file size...'
	tmp = os.popen('wc -l '+filename)
	nlines = string.atof(tmp.readline().split(' ')[0])
	
	if(tmp.close()<>None):
		print 'Error: some error occurred while running '
		print '	wc -l '+filename
		print 'at bam_file.py'
		print 'Exiting'
		sys.exit(1)		
	print '	Done.'
	
	return nlines





def region_coverage(coveragefile):
	
	
	# A progress bar is initialized
#	print 'Calculating mean coverage per region...'

#	widgets = ['Calculating mean coverage per region: ', progressbar.Percentage(), ' ', 
#				progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
#	pbar = progressbar.ProgressBar(widgets=widgets, maxval=count_lines(coveragefile)).start() 
	
	i = 1
	coverage = {}
	curr_region = None
	fd = file(coveragefile)
	for line in fd:
		parts = line.split('\t')
		newregion = (parts[0], string.atoi(parts[1]), string.atoi(parts[2]))

		if(curr_region<>newregion):
			if(curr_region<>None):
				coverage[curr_region] = numpy.mean(sampling)

			sampling = []
			curr_region = newregion 

		sampling.append(string.atoi(parts[-1]))
#		pbar.update(i)
		i+=1

	coverage[curr_region] = numpy.mean(sampling)
				   
#	pbar.finish() 
	fd.close()
	
	return coverage
	
	
	
	
	
def gcbias(filelist, fileoutlist, bedfilelist):
	"""************************************************************************************************************************************************************
	Task: draws coverage as a function of gc content
	Input:
		filelist: list of strings, each containing the full path of the bam file to analyze.
		fileoutlist: list of strings, each containing the full path of the png file where the corresponding figure will be saved.
		bedfilelist: 
	Output: a bmp file will be created named "fileout" where a graph that compares gc content and mean coverage will be saved.	
	************************************************************************************************************************************************************"""
	
	pid = str(os.getpid())
	
	numpy.random.seed(1)
	ntotal_positions = []
	bamlist = []
	
	# Process each file and store counting results
	for filename in filelist:
		# Check whether index already exists for the bam file, needed for pysam use
		if(not os.path.isfile(filename+'.bai')):
			print 'Creating index for '+filename
			pysam.index(filename)
			print '	Done.'
						
		bamlist.append(bam_file.bam_file(filename))
	sizes = numpy.array([bam.nreads() for bam in bamlist])
	minsize = sizes.min()
	
	print 'The smaller bam is '+filelist[sizes.argmin()]+' and contains '+str(minsize)+' reads.'
		
	# Process each file and store counting results
	for i,bamfile in enumerate(bamlist):
	
		print 'Processing '+bamfile.filename
		print 'Results will be written at '+fileoutlist[i]
		
		# Check whether normalization should be run
		if(normalize): normalizedbam = bamfile.normalize(minsize)
		else: normalizedbam = bamfile
		
		coveragefile = TMP+'/'+pid+'.coverage'
		print 'Calculating coverage per position...'
		run(BEDTOOLSPATH+'coverageBed -d -abam '+normalizedbam.filename+' -b '+bedfilelist[i]+' > '+coveragefile)   
	
		coverage = region_coverage(coveragefile)
	
		print 'Calculating nt content...'
		bedfd = pybedtools.BedTool(bedfilelist[i])
		pybedtools._bedtools_installed = True
		pybedtools.set_bedtools_path(BEDTOOLSPATH)	
		ntcontent = bedfd.nucleotide_content(REF)
		
		# Each entry in ntcontent is parsed to extract the gc content of each exon
		gccontent = {}
		for entry in ntcontent:
			gccontent[(entry.fields[0], string.atoi(entry.fields[1]), string.atoi(entry.fields[2]))] = string.atof(entry.fields[-8])*100
		print '	Done.'		
			
		fig = pyplot.figure(figsize=(13,6))
		ax = fig.add_subplot(111)
		
		region_ids = coverage.keys()
		coveragearray = numpy.array([coverage[id] for id in region_ids])
		gccontentarray = numpy.array([gccontent[id] for id in region_ids]) # Values in [0,1]
	
		xmin = gccontentarray.min()
		xmax = gccontentarray.max() # Due to the imshow sentence, we need to rescale gccontent from [0,1] to [0,100]
		ymin = coveragearray.min()
		ymax = coveragearray.max()
		 
		# Perform a kernel density estimator on the results
		X, Y = mgrid[xmin:xmax:100j, ymin:ymax:100j]
		positions = c_[X.ravel(), Y.ravel()]
		values = c_[gccontentarray, coveragearray]
		kernel = stats.kde.gaussian_kde(values.T)
		Z = reshape(kernel(positions.T).T, X.T.shape)
		
		
		fig = pyplot.figure(figsize=(6,6))
		ax = fig.add_subplot(111)
		sc=ax.imshow(rot90(Z),cmap=cm.gist_earth_r,extent=[xmin, 100, ymin, ymax], aspect="auto") # Due to the imshow sentence, we need to rescale gccontent from [0,1] to [0,100]
		cbar=fig.colorbar(sc,ticks=[numpy.min(Z),numpy.max(Z)])
		cbar.ax.set_yticklabels(['Low','High'])
		cbar.set_label('Density')
		ax.set_xlabel('GC content (%)')
		ax.set_ylabel('Mean coverage')
		fig.savefig(fileoutlist[i])
		matplotlib.pyplot.close(fig)
	
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Finished <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'





def measureGCbias(wholeChromosome,currentChromosome,bedFile):
    """************************************************************************************************************************************************************
    Task: measures GC content in a given sequence according to the target file
    Input:
        wholeChromosome: sequence corresponding to a given chromosome
        currentChromosome: chromosome to which the sequence belongs
        bedFile: target file object        
    Output: a dictionary with the percentage of GC for each region of the target file    
    **"""
    
    gccontent = {}
    

    # For each exon, count the gc content according to the reference loaded for currentChromosome (wholeChromosome)
    for currentExon in bedFile.chrs[str(currentChromosome)]:
        exon_init=int(currentExon[0])
        exon_end=int(currentExon[1])
        
        currentSequence=wholeChromosome[(exon_init-1):exon_end] # The first component is at 0 index and the last component is at n-1 index
        countsG=currentSequence.count('G')
        countsC=currentSequence.count('C')
        
    
        gccontent[currentChromosome,exon_init,exon_end]=((countsG+countsC)/float(len(currentSequence)))*100.0 # Percentage of GC content
                
    return gccontent





def gcbias_lite(coveragefile, bedfilename, reference, fileout, graphtitle=None, executiongranted=None, status=None, bedTools=False):
	"""************************************************************************************************************************************************************
	Task: draws coverage as a function of gc content. IMPROVED VERSION of gcbias that avoids the use of bedtools (pybedtools)
	Input:
		coveragefile: string containing the full path of the bam.coverage file to analyze. This file has been built according to 1-base format
		bedfilename: target file -> assumes original-standard bed file
		reference: fasta file with reference genome
		fileout: string containing the full path of the bmp file where the restulting figure will be saved.
		bedTools: whether pybedtools are used instead of the own method
	Output: a png file will be created named "fileout" where a graph that compares gc content and mean coverage will be saved.	
	************************************************************************************************************************************************************"""
	   
	if(executiongranted<>None):
		executiongranted.acquire()
	
	pid = str(os.getpid())
 
#	print 'Processing '+coveragefile
#	print 'Results will be written at '+fileout
	coverage = region_coverage(coveragefile) # Calculate mean coverage per region
	
##	fdw=file('regionCoverage.txt','w')	
##	for element in sorted(coverage.keys()):
##		fdw.write(str(element)+'\n')		
##	fdw.close()

	if(len(coverage)>1):	
		
		if not bedTools:   # Own method
#			print 'Own method'
			chromosomes={}	 
			allKeys=coverage.keys()
			
			for currentKey in allKeys:
				chromosomes[currentKey[0]]=1 # Stores all chromosomes to be examined (the ones contained in the target file)
						
			# Load BED file -> since coverage information is in 1-base format, BED format must be transformed to 1-base
			bed=bed_file.bed_file(bedfilename)
			sortedBed=bed.my_sort_bed() # Sort bed avoiding bedtools
			nonOverlappingBed=sortedBed.non_overlapping_exons(1) # Base 1!!! # This generates a BED file in base 1 (Non-standard BED)
			finalBed=nonOverlappingBed.my_sort_bed() # BED file in base 1 (Non-standard BED)
			finalBed.load_custom(-1) # Load chromosome and positions in base 1....(finalBed is in base 1 -> Non-standard BED)	
	
						
			#Load FASTA file		
			fastaFile=file(reference,'r')
			
			storeSequence=False
			wholeChromosome=''
			currentChromosome=''
			gccontent={}		
	
		
			for line in fastaFile: # Read each line of the fasta file
				if line.startswith('>'): # New chromosome starts -> reading a new line until another '>' is found
#					print 'Processing ' +line+'\n' 
					if storeSequence: # a chromosome has been read run gc bias				
						currentGCcontent=measureGCbias(wholeChromosome,currentChromosome,finalBed)
						gccontent.update(currentGCcontent) # Update dictionary
						storeSequence=False
					currentChromosome=re.split(' +',line)[0] # Format: >1 dna:chromosome chromosome:GRCh37:1:1:249250621:1
					currentChromosome=currentChromosome.split('>')[1].strip() # Chromosome string
					if(currentChromosome in chromosomes): # If current chromosome read in the FASTA file is in the list of chromosomes in the BED file
						storeSequence=True
					wholeChromosome='' # To store whole sequence for the current chromosome
				elif (not re.search('>',line) and storeSequence):
					wholeChromosome=wholeChromosome+line.rstrip() # Remove '\n' from current line and concatenates to wholeChromosome
					
	
			if(storeSequence): # For the last chromosome
					currentGCcontent=measureGCbias(wholeChromosome,currentChromosome,finalBed)
					gccontent.update(currentGCcontent)  # Update dictionary
					
			fastaFile.close()  
			region_ids=[]					
			region_ids = coverage.keys()
			
			if(len(gccontent)==0):
				print 'ERROR: G+C content values can not be calculated. Probably the provided reference file '+reference+' does not match with '
				print '	the target file '+bedfilename+'. That is, sequences of regions in the target file are probably not included within the'
				print '	reference file.'
				sys.exit(1)
			   
		else:			
			print 'Calculating nt content by means of pybedtools...'
			bed=bed_file.bed_file(bedfilename)
			sortedBed=bed.my_sort_bed() # Sort bed avoiding bedtools
			nonOverlappingBed=sortedBed.non_overlapping_exons(1) # base one!!! 
			finalBed=nonOverlappingBed.my_sort_bed() # BED file in base 1
			bedfd = pybedtools.BedTool(finalBed.filename)
			bedfd=bedfd.remove_invalid() # Remove negative coordinates or features with length=0, which do not work with bedtools
			pybedtools._bedtools_installed = True
			pybedtools.set_bedtools_path(BEDTOOLSPATH)	
			ntcontent = bedfd.nucleotide_content(reference)
				
			# Each entry in ntcontent is parsed to extract the gc content of each exon
			gccontent = {}
			for entry in ntcontent:
				gccontent[(entry.fields[0], string.atoi(entry.fields[1]), string.atoi(entry.fields[2]))] = string.atof(entry.fields[-8])*100
			print '	Done.'						
			# gccontent keys in dictionary: chromosome, exon init, exon end   
			
			region_ids=[]
			for currentKey in coverage.keys(): # Pybedtools does not work with regions with zero length -> remove them (there are a few of them)
				if currentKey[1]!=currentKey[2]:
					region_ids.append(currentKey)
						
		
##		
##		fdw=file('gcContent.txt','w')	
##		for element in sorted(gccontent.keys()):
##			fdw.write(str(element)+'\n')		
##		fdw.close()
##			
		#region_ids = gccontent.keys()
		coveragearray = numpy.array([coverage[id] for id in region_ids])
		gccontentarray = numpy.array([gccontent[id] for id in region_ids]) # Values in [0,1]	
				
#		fig = pyplot.figure(figsize=(6,6))
#		ax = fig.add_subplot(111)
#		
#		ax.hist(gccontentarray,bins=100)
#		fig.suptitle('Dsitribution of GC content regardless of coverage value')	
#		ax.set_ylabel('Frequency')
#		ax.set_xlabel('GC content')
#		ax.set_xlim(0, 100)
#		fig.savefig('distribution.png')										
					
		xmin = gccontentarray.min()
		xmax = gccontentarray.max() # Due to the imshow sentence, we need to rescale gccontent from [0,1] to [0,100]
		ymin = coveragearray.min()
		ymax = coveragearray.max()
		 
		# Perform a kernel density estimator on the results
		X, Y = mgrid[xmin:xmax:100j, ymin:ymax:100j]
		positions = c_[X.ravel(), Y.ravel()]
		values = c_[gccontentarray, coveragearray]
		kernel = stats.kde.gaussian_kde(values.T)
		Z = reshape(kernel(positions.T).T, X.T.shape)
		
		
		fig = pyplot.figure(figsize=(6,6))
		ax = fig.add_subplot(111)
		sc=ax.imshow(rot90(Z),cmap=cm.gist_earth_r,extent=[xmin, 100, ymin, ymax], aspect="auto") # Due to the imshow sentence, we need to rescale gccontent from [0,1] to [0,100]
		cbar=fig.colorbar(sc,ticks=[numpy.min(Z),numpy.max(Z)])
		cbar.ax.set_yticklabels(['Low','High'])
		cbar.set_label('Density')
		ax.set_xlabel('GC content (%)')
		ax.set_ylabel('Mean coverage')
		
		if(len(graphtitle)>25):
			ax.set_title(graphtitle[:25]+'...')
		else:
			ax.set_title(graphtitle)
			
		fig.savefig(fileout)
		matplotlib.pyplot.close(fig)
		
		if(status<>None):
			meanvalue = gccontentarray.mean()
			status.value = (meanvalue>=45 and meanvalue<=55)
		

	else:
		print 'WARNING: only one region found in the bed file. Skipping GC bias calculation.'
		
	if(executiongranted<>None):
		executiongranted.release()

	
#gcbias_lite('/tmp//GU_20130110_FC1_6L4S_JA_C148_F3.filtered.realigned.recalibrated.24796.coverage', '/tmp/test.bed', 
#		    '/home/javi/MGP/mitochondria/lib/human_g1k_v37.MT.fasta', '/tmp/test.png')

#gcbias_lite('/tmp/3314.coverage', '/home/javi/MGP/capture_methods/lib/SeqCap_EZ_Exome_v3_primary.g1k.bed', '/tmp/test.png')
#gcbias_lite('/tmp/13483.coverage', '/tmp/ref_mt.bed', '/tmp/test.png')

def main():   
	################################################

	#### Options and arguments #####################

	################################################
	usage="""	
	************************************************************************************************************************************************************
	Task: draws coverage as a function of gc content		
	************************************************************************************************************************************************************


	
	usage: %prog -i <bamfile> -o <fileout>"""
	
	parser = optparse.OptionParser(usage)
	parser.add_option("-i", dest="bamfilelist", help="""String containing a comma-separated list with the names of the bams to analyze.""")
	parser.add_option("-b", dest="bedfilelist", help="""String containing a comma-separated list with the names of the bed files that contain the regions to analyze.""")
	parser.add_option("-o", dest="fileoutlist", help="""String containing the name of the (png) file where the figure will be saved.""")
	(options, args) = parser.parse_args()

	# Check number of arguments	
	if len(sys.argv) < 7:
		parser.print_help()
		sys.exit(1)
   
	# call core function
	#gcbias(options.bamfile, options.fileout)
	gcbias(options.bamfilelist.split(','), options.fileoutlist.split(','), options.bedfilelist.split(',')) # Improved version





if __name__=='__main__':
	main()

