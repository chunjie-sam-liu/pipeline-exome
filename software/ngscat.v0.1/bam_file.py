try:
	import numpy
except ImportError:
	print 'WARNING: module numpy was not loaded.'	

import gc
import time
import random
import pysam
import copy
import os
import sys
import string
import sets

try:
	import progressbar
except ImportError:
	print 'WARNING: module progressbar was not loaded.'
	
import xlwt

import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot


import bed_file
import multiprocessing

#TMP = '/home/fjavier/tmp/'
TMP = '/tmp/'
#CHR_LENGTHS = '/data/reference_genomes/human/human_g1k_v37.genome'
CHR_LENGTHS = '/home/javi/MGP/data/reference_genomes/human/human_g1k_v37.genome'




class bam_file(pysam.Samfile):
	
	
	
	
	
	def __init__(self, _filename=None, mode='rb', header=None):
#		self.filename = _filename
		self._nreads = None
		
		pysam.Samfile.__init__(self, _filename, mode, header=header)
			
	
	
	
	  
		
		
		
		
	def issorted(self):
		print 'Checking sorting of '+self.filename
		
		visitedcontigs = sets.Set()
		previouscontig = None
		previousstart = None
		previousend = None
		
		try:
			read = self.next()
			readsavailable=True
			sorted=True
			currend=read.pos+sum([nbases[1] for nbases in read.cigar if nbases[0]<>1])-1
		except StopIteration:
			readsavailable=False		
			
		rc = 1
#		while(readsavailable and ((previouscontig<>read.rname and (read.rname not in visitedcontigs)) or 
#								 (previousstart<read.pos or (previousstart==read.pos and previousend<=currend)))):
		while(readsavailable and ((previouscontig<>read.rname and (read.rname not in visitedcontigs)) or 
								 (previousstart<=read.pos))):

			if(not rc%1000000):
				print str(rc)+' reads checked'
				
			visitedcontigs.add(read.rname)
			previouscontig=read.rname
			previousstart=read.pos 
			previousend=currend
			
			try:
				read = self.next()
				readsavailable=True
				currend=read.pos+sum([nbases[1] for nbases in read.cigar if nbases[0]<>1])-1				
			except StopIteration:
				readsavailable=False
			except TypeError:
				if(read.is_unmapped):
					print 'ERROR: unmapped reads found at '+self.filename
				else:
					print 'ERROR: incorrect bam format'
				print '	Read position: '+str(rc)
				print '	Alignment entry: '+str(read)
				print '	Exiting.'
				sys.exit(1)
				
			rc += 1

		print '	Done.'		
		
		if(readsavailable):
			print 'WARNING: bam not sorted.'
			print '	Check read '+ str(read)
			print '	and the read before.'
			
			return False
		else:
			return True
		
		
		
	def run(self, command):
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
			print '	'+command





	def nreads(self):
		return self.mapped+self.unmapped
#		print 'Calculating number or reads...'
#
#		if(self._nreads<>None): return self._nreads
#		
#		tmp = os.popen('samtools view '+self.filename+' | wc -l')
#		nlines = string.atof(tmp.readline().split(' ')[0])
#		
#		if(tmp.close()<>None):
#			print 'Error: some error occurred while running '
#			print 'samtools view '+self.filename+' | wc -l'
#			print 'at bam_file.py'
#			print 'Exiting'
#			sys.exit(1)		
#		print '	Done.'
#		
#		return nlines
		



		

	def mappingsize(self):
		totalsize=0
		
		for contig in self.header['SQ']:
			totalsize += contig['LN']
			
		return totalsize
		
		
		
		
		
	def skip_header(self, fd):
		line = fd.readline()
		while(line.startswith('@')):
			line = fd.readline()
			
		return fd
	
	
	
	
	
	def get_mapped_reads(self):
	
		print 'Loading mapped reads...'
		nunmapped = 0
		mapped_reads = {}
		for read in self:
			if(read.is_unmapped):
				nunmapped += 1
			else:
				mapped_reads[read.qname] = read

		print str(len(mapped_reads))+' mapped reads loaded'
		print str(nunmapped)+' unmapped reads found'
				
		return mapped_reads





	def get_mapped_reads_ids(self):
	
		print 'Loading mapped reads ids...'
		nunmapped = 0
		mapped_reads = {}
		for read in self:
			if(read.is_unmapped):
				nunmapped += 1
			else:
				mapped_reads[read.qname] = None

		print str(len(mapped_reads))+' mapped reads loaded'
		print str(nunmapped)+' unmapped reads found'
				
		return mapped_reads
		
		
		
		
	def count_lines(self, filename):
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





	def target_coverage(self, coveragelist, targetfile):
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
	
		pid = str(os.getpid())
		coveragefile = TMP+'/'+pid+'.coverage'
		
		print 'Calculating coverage per target position...'
		self.run('coverageBed -d -abam '+self.filename+' -b '+targetfile+' > '+coveragefile)
		
		ntotal_positions = self.count_lines(coveragefile)

		# covered_positions_per_depth: list of integers. There will be a position for each coverage threshold. Each value will be the count of positions
		#	 covered for the corresponding threshold.
		# ccds_counts: dictionary. Keys are transcript ids. values are lists of two elements. The first element of this list will contain the length of the
		#	 transcript. The second element will be a list of integers with as many positions as coverage thresholds, being each value the count of positions
		#	 covered for the corresponding threshold. 
		covered_positions_per_depth = [0 for i in range(len(coveragelist))]

		# A progress bar is initialized
		widgets = ['Counting: ', progressbar.Percentage(), ' ', 
					progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
		pbar = progressbar.ProgressBar(widgets=widgets, maxval=ntotal_positions).start() 

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
			
			pbar.update(k+1)
			
		pbar.finish()
		fd.close()
					   
		return [ntotal_positions, covered_positions_per_depth]
	
	
	


	def generate_gff(self, chrs, sizes, binsize):
		pid = str(os.getpid())
		gfffile = TMP+'/'+pid+'.gff'
		fd = file(gfffile, 'w')
		
		for chr in chrs:
			for base in range(0,sizes[chr],binsize):
				fd.write(chr+'\ttmp\texon\t'+str(base+1)+'\t'+str(base+binsize)+'\t0\t+\t.\tcoord "'+chr+':'+str(base+1)+'";\n')
				
		fd.close()
		
		
		
	
	def coverage_distribution(self, fileout, binsize=10000000):
		pid = str(os.getpid())
		chrs = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
		base0 = {}
		lengths = {}
		
		tmpbed = bed_file.bed_file(TMP+'/'+pid+'.bed')
		coveragefile = TMP+'/'+pid+'.coverage'
			   
		print 'Loading chr lengths...'
		fd = file(CHR_LENGTHS)
		for line in fd:
			parts = line.split('\t') 
			lengths[parts[0]] = string.atoi(parts[1])
		print '	Done.'
		fd.close()
			
		print 'Calculating base 0 coordinates for each chr...'
		base0['1']=0
		for i in range(1,len(chrs)):
			base0[chrs[i]]=base0[chrs[i-1]]+lengths[chrs[i-1]]
		print '	Done.'
			
#		totallength = base0[chrs[-1]]+lengths[chrs[i-1]]
		print 'Calculating histogram...'
		points = []
		gfffile = TMP+'/'+pid+'.gff'
		countsfile = TMP+'/'+pid+'.counts'

		self.generate_gff(chrs, lengths, binsize)
		self.run("""samtools view """+self.filename+""" | htseq-count -q -i coord - """+gfffile+""" > """+countsfile)
		
		fd = file(countsfile)
		line = fd.readline()
		while('no_feature' not in line):
			chrcoord,counts = line.split()
			chr,coord = chrcoord.split(':')
			points.append((base0[chr]+string.atoi(coord)+binsize/2, string.atoi(counts)))
			line = fd.readline()
					   
		points.sort()	
		fig = pyplot.figure(figsize=(13,6))
		ax = fig.add_subplot(111)
		ax.plot([point[0] for point in points], [point[1] for point in points])
		
		ax.set_xlim(right=lengths['Y']+base0['Y'])
		ax.set_ylabel('# reads')
		ax.set_xlabel('Chr position')
		
				
		for chr in base0: ax.axvline(base0[chr], color='#ff0000', linestyle='--', linewidth=0.5, alpha=0.5)
		
		fig.savefig(fileout)
		matplotlib.pyplot.close(fig)
			
			
		
		
		
	def enrichment(self, ontarget, offtarget, dirout):
		"""*******************************************************************************************************************************************
		Task: calculates the enrichment in mapped reads of targeted regions as (#reads/kbase_ontarget)/(#reads/kbase_offtarget) 
		Inputs:
			ontarget, offtarget: strings containing the full path to the bed files containing on-target and off-target regions. We implemented
				this method for testing the enrichment of target regions compared to whole genome gencode regions, that is why off-target reads
				are not simply reads not mapped on-target.
		Outputs: 
			nreads_on, notonbam.nreads(), nreads_off: number of reads on target, number of reads that do not map on target, and number of reads
				that do not map on target but do map on the regions indicated at "offtarget".
			Reads that do no map on-target but do map on "offtarget" per Kbase
			Reads that map on-target per Kbase
			Enrichment, calculated as (#reads/kbase_ontarget)/(#reads/kbase_offtarget)
		*******************************************************************************************************************************************"""
		
		pid = str(os.getpid())
		
		on = bed_file.bed_file(ontarget)
		off = bed_file.bed_file(offtarget)	   

		onbam = bam_file(TMP+'/'+pid+'.on.bam')
		notonbam = bam_file(TMP+'/'+pid+'.noton.bam')
		offbam = bam_file(TMP+'/'+pid+'.off.bam')
										
		self.run(""" intersectBed -abam """+self.filename+""" -b """+ontarget+""" > """+onbam.filename)
		
		# The set of reads that do not map on-taret is obtained. This set is then intersected with the "offtarget" bed.
		self.run(""" intersectBed -v -abam """+self.filename+""" -b """+ontarget+""" > """+notonbam.filename)
		self.run(""" intersectBed -abam """+notonbam.filename+""" -b """+offtarget+""" > """+offbam.filename)
		
		nreads_off = offbam.nreads()
		offbam.coverage_distribution(dirout+'/'+os.path.basename(self.filename).replace('.bam','.png'))
		nreads_on = onbam.nreads()
		nreads_noton = notonbam.nreads()
		ontarget_size = on.size()
		offtarget_size = off.size()

		
		os.remove(onbam.filename)
		os.remove(notonbam.filename)
		os.remove(offbam.filename)
		
		return [nreads_on, nreads_noton, nreads_off, nreads_off*1000.0/offtarget_size, nreads_on*1000.0/ontarget_size, 
				(nreads_on*1.0/ontarget_size)/(nreads_off*1.0/offtarget_size)]
		




	def draw_coverage_distribution(self, on_points, noton_points, off_points, base0,lengths, fileout):
		fig = pyplot.figure(figsize=(13,6))
		ax = fig.add_subplot(111)
		
		
		rects = []		
		rects.append(ax.plot([point[0] for point in on_points], [point[1] for point in on_points], color="#ff0000"))
		rects.append(ax.plot([point[0] for point in noton_points], [point[1] for point in noton_points], color="#00ff00"))
		rects.append(ax.plot([point[0] for point in off_points], [point[1] for point in off_points], color="#0000ff"))
	
		legend = ['On-target', 'Off-target', 'Off-target in gencode']
			
		ax.set_xlim(right=lengths['Y']+base0['Y'])
		ax.set_ylabel('# reads')
		ax.set_xlabel('Chr position')
		
		# Shink current axis by 20%
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	
		ax.legend( tuple([rect[0] for rect in rects]), tuple(legend), loc="upper left", bbox_to_anchor=(1,1) )
				
		for chr in base0: ax.axvline(base0[chr], color='#555555', linestyle='--', linewidth=0.5, alpha=0.5)
		
		fig.savefig(fileout)
		matplotlib.pyplot.close(fig)
		
		
		
		
		
		
	def coverage_distribution2(self, binsize=10000000):
		pid = str(os.getpid())
		chrs = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
		base0 = {}
		lengths = {}
		
		tmpbed = bed_file.bed_file(TMP+'/'+pid+'.bed')
		coveragefile = TMP+'/'+pid+'.coverage'
			   
		print 'Loading chr lengths...'
		fd = file(CHR_LENGTHS)
		for line in fd:
			parts = line.split('\t') 
			lengths[parts[0]] = string.atoi(parts[1])
		print '	Done.'
		fd.close()
			
		print 'Calculating base 0 coordinates for each chr...'
		base0['1']=0
		for i in range(1,len(chrs)):
			base0[chrs[i]]=base0[chrs[i-1]]+lengths[chrs[i-1]]
		print '	Done.'
			
#		totallength = base0[chrs[-1]]+lengths[chrs[i-1]]
		print 'Calculating histogram...'
		points = []
		gfffile = TMP+'/'+pid+'.gff'
		countsfile = TMP+'/'+pid+'.counts'

		self.generate_gff(chrs, lengths, binsize)
		self.run("""samtools view """+self.filename+""" | htseq-count -q -i coord - """+gfffile+""" > """+countsfile)
		
		fd = file(countsfile)
		line = fd.readline()
		while('no_feature' not in line):
			chrcoord,counts = line.split()
			chr,coord = chrcoord.split(':')
			points.append((base0[chr]+string.atoi(coord)+binsize/2, string.atoi(counts)))
			line = fd.readline()
					   
		points.sort()	

		return [base0,lengths, points]




	def enrichment2(self, ontarget, offtarget, dirout):
		"""*******************************************************************************************************************************************
		Task: calculates the enrichment in mapped reads of targeted regions as (#reads/kbase_ontarget)/(#reads/kbase_offtarget) 
		Inputs:
			ontarget, offtarget: strings containing the full path to the bed files containing on-target and off-target regions. We implemented
				this method for testing the enrichment of target regions compared to whole genome gencode regions, that is why off-target reads
				are not simply reads not mapped on-target.
		Outputs: 
			nreads_on, notonbam.nreads(), nreads_off: number of reads on target, number of reads that do not map on target, and number of reads
				that do not map on target but do map on the regions indicated at "offtarget".
			Reads that do no map on-target but do map on "offtarget" per Kbase
			Reads that map on-target per Kbase
			Enrichment, calculated as (#reads/kbase_ontarget)/(#reads/kbase_offtarget)
		*******************************************************************************************************************************************"""
		
		# pid: integer containing current process id
		pid = str(os.getpid())
		
		# on, off: bed_file objects representing on and off targets
		on = bed_file.bed_file(ontarget)
		off = bed_file.bed_file(offtarget)	   

		# notonbamfilename: string containing the name of the .bam where reads that do not overlap with the target will be stored.
		# offbamfilename: string containing the name of the .bam where reads that do not overlap with the target BUT that do overlap with the offtarget bed will be stored.
		# onbamfilename: string containing the name of the .bam containing those reads that overlap with ontarget.
		notonbamfilename = TMP+'/'+pid+'.noton.bam'
		offbamfilename = TMP+'/'+pid+'.off.bam'
		onbamfilename = TMP+'/'+pid+'.on.bam'
										
		# Extract those reads overlapping with the target
		self.run(""" intersectBed -abam """+self.filename+""" -b """+ontarget+""" > """+onbamfilename)	   
		
		# The set of reads that do not map on-tarfet is obtained. This set is then intersected with the "offtarget" bed.
		self.run(""" intersectBed -v -abam """+self.filename+""" -b """+ontarget+""" > """+notonbamfilename)
		self.run(""" intersectBed -abam """+notonbamfilename+""" -b """+offtarget+""" > """+offbamfilename)

		# Index the three new bam. Neccessary to avoid pysam to fail.
		pysam.index(onbamfilename)
		pysam.index(notonbamfilename)
		pysam.index(offbamfilename)
		
		# Generate bam_file objects for the three bam
		onbam = bam_file(onbamfilename, 'rb')
		notonbam = bam_file(notonbamfilename, 'rb')
		offbam = bam_file(offbamfilename, 'rb')
		
		# nreads_off: integer containing the number of reads that do not overlap with ontarget BUT that do overlap with offtarget
		# base0: dictionary where keys are chromsomes and values are integers, each containing the coordinate that represents the base 0 of the chromome, considering that all the chromosomes will be drawn one before another.
		# lengths: dictionary where keys are chromosomes and values the length of each chromosome.
		# onpoints: tuples of integers (x,y) representing the read count (y) at genomic position x
		# nreads_on: integer containing the number of reads that do overlap with the target
		# nreads_noton: integer containing the number of reads that do not overlap with the target
		# ontarget_size: integer containing the number of bases covered by the ontarget. Bases falling in overlapped regions are counted just once.
		# offtarget_size: integer containing the number of bases covered by the offtarget. Bases falling in overlapped regions are counted just once.
		nreads_off = offbam.mapped
		base0,lengths,offpoints = offbam.coverage_distribution2()
		nreads_on = onbam.mapped
		base0,lengths,onpoints = onbam.coverage_distribution2()
		nreads_noton = notonbam.mapped
		base0,lengths,notonpoints = notonbam.coverage_distribution2()
		self.draw_coverage_distribution(onpoints, notonpoints, offpoints, base0,lengths, dirout+'/'+os.path.basename(self.filename).replace('.bam','.png'))
		ontarget_size = on.size()
		offtarget_size = off.size()

		# Remove temporary bams
		os.remove(onbam.filename)
		os.remove(notonbam.filename)
		os.remove(offbam.filename)
		
		return [nreads_on, nreads_noton, nreads_off, nreads_off*1000.0/offtarget_size, nreads_on*1000.0/ontarget_size, 
				(nreads_on*1.0/ontarget_size)/(nreads_off*1.0/offtarget_size)]
	
	
	
	
	
	def select_reads(self, n):
		numpy.random.seed(1)
		
		selected = numpy.random.uniform(size=self.nreads())<=(n*1.0/self.nreads())
		nselected = len(selected.nonzero()[0])

		i=0
		while(i<len(selected) and nselected<n):
			if(not selected[i]): 
				selected[i]=True
				nselected += 1
			i+=1
				
		while(i<len(selected) and nselected>n):
			if(selected[i]):
				selected[i]=False
				nselected -= 1
			i+=1
				
		return selected.nonzero()[0]
		
		
		
		

	def select_reads_boolean(self, n):
		"""*******************************************************************************************************************************************
		JPFLORIDO
		Task:  Randomly selected a given set of reads from a BAM file
		Inputs: BAM file and number of reads to be selected
		Outputs: Boolean array indicating whether i-th read has been selected or not
		*******************************************************************************************************************************************"""
	
		numpy.random.seed(1)
		
		selected = numpy.random.uniform(size=self.nreads())<=(n*1.0/self.nreads())
		nselected = len(selected.nonzero()[0])

		i=0
		while(i<len(selected) and nselected<n):
			if(not selected[i]): 
				selected[i]=True
				nselected += 1
			i+=1
				
		while(i<len(selected) and nselected>n):
			if(selected[i]):
				selected[i]=False
				nselected -= 1
			i+=1
				
		return selected




		
	def sort_bam(self):
		"""*******************************************************************************************************************************************
		JPFLORIDO. 
		Task:  Sort BAM file by position and creat the corresponding index
		Inputs: BAM file to be sorted
		Outputs: Sorted file (object) with the corresponding index created
		*******************************************************************************************************************************************"""
		pid = str(os.getpid())
		print 'Sorting BAM according to position...'
		sortedBAMfilename = TMP+'/'+pid+ os.path.basename(self.filename)+".sorted"
#		self.run('samtools sort '+self.filename+' '+sortedBAMfilename)
		pysam.sort(self.filename, sortedBAMfilename)
		
		# Index sorted BAM	
#		command='samtools index '+sortedBAMfilename+'.bam'
#		fd = os.popen(command)
#		if(fd.close()<>None):
#			print 'Some error occurred while executing: '
#			print '	'+command

		pysam.index(sortedBAMfilename+'.bam')

		return bam_file(sortedBAMfilename+'.bam','rb')





	def normalize(self, n):
		if(n>=self.nreads()):
			print 'File: '+self.filename+'. No actual normalization was done, current number of reads is '+str(self.nreads())+'. Required number of reads is '+str(n)+'.' 
			return bam_file(self.filename, 'rb')
		
		pid = str(os.getpid())
		selectionprob = n*1.0/self.nreads()
		normalizedbam = bam_file(TMP+'/'+pid+'.bam', 'wb', header=self.header)
		
		i = self.fetch()
		currread = 0
		nselected = 0
		niterated = 0

		remove = self.select_reads(self.nreads()-n)
		
		# A progress bar is initialized
		widgets = ['Writing selected reads: ', progressbar.Percentage(), ' ', 
					progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
		pbar = progressbar.ProgressBar(widgets=widgets, maxval=self.nreads()).start() 

		try:
			previous = 0
			for next in remove:
				for j in range(previous,next):
					normalizedbam.write(i.next())
					nselected += 1
					pbar.update(j)
				i.next()
				previous = next+1									   
				
			for j in range(previous, self.nreads()): 
				normalizedbam.write(i.next())
				nselected += 1
			
		except StopIteration, stop:
			print 'ERROR: method normalize at bam_file.py.'
			print '	n = '+str(n)   
			print '	self.nreads()-n = '+str(self.nreads()-n)		 
			print '	j = '+str(j)
			print '	nselected = '+str(nselected)
			print '	previous = '+str(previous)
			print '	next = '+str(next)
			print 'Iteration should have stopped before end of bam.'
			sys.exit(1)
						
		pbar.finish()
		
		print str(nselected)+' reads written'
		
		normalizedbam.close()
		pysam.index(TMP+'/'+pid+'.bam')
		
		return bam_file(TMP+'/'+pid+'.bam', 'rb')
		
		 
	
	
	
	def ncovered_reads(self, target):
		pid = str(os.getpid())
		newbamfilename = TMP+'/'+pid
		self.run(""" intersectBed -abam """+self.filename+""" -b """+target+""" > """+newbamfilename)
		pysam.index(newbamfilename)
		bam = bam_file(newbamfilename, 'rb')
		
		return bam.nreads()
	
	
	

	def myCoverageBed(self,target,numberReads=None,writeToFile=None,executiongranted=None,tmpdir=None,bedGraphFile=None):
		"""*******************************************************************************************************************************************
		JPFLORIDO
		Task:  Custom method equivalent to bedtool's coverage bed
		Inputs: BAM file, target file (BED file), number of desired reads to be taken into account from the BAM file (default all) and option to write results to file the same way as bedtool does
		Outputs: 
			positionArray_intersect: a vector of bp coordinates in which coverage changes according to the intersection of BAM file and target file
			coverageArray_intersect: coverage (number of reads) for a given bp position 
			chromosomeCoordinates: a dictionary in which each entry is related to a chromosome and the content indicates the starting and ending position in positionArray_intersect and coverageArray_intersect arrays
			finalBed: target file with overlapping exons removed and sorted (object)		
		Other issues: it is intented to return also positions off target. However, coverage around exons changes a lot, so there is an important increasing
		# in the memory used 
		*******************************************************************************************************************************************"""

		global TMP
		
		if(executiongranted<>None):
			executiongranted.acquire()
			
		if(tmpdir<>None):
			TMP = tmpdir
			
		pid = str(os.getpid())
		
		# Check whether BAM file is sorted
		command='samtools view -H '+ self.filename+' | grep SO:'
		fd=os.popen(command)
		outputCommand=fd.read()
		
		
		
		if ((fd.close()<>None) or('coordinate' not in outputCommand.split('SO:')[-1])): # BAM not sorted -> Sort and indexing BAM
#			print 'ERROR: bam file must be sorted.'
#			print '	Exiting.'
#			sys.exit(1)
			sortedBam=self.sort_bam()
		else:
			sortedBam=self
			
			
		if(writeToFile!=None): # Results written to output file
			fdw = file(writeToFile, 'w')
			
		
		positionArray_intersect=[]
		coverageArray_intersect=[]
		chromosomeCoordinates={} # Dictionary that controls, for a given chromosome, its starting position and end position in the previous intersection arrays
		#positionArray_offTarget=[] # Stores positions that are off target -> Finally not used 
		#coverageArray_offTarget=[] # Stores coverage that are off target -> finally not used		
		#chromosomeCoordinatesOffTarget={}# Dictionary that controls, for a given chromosome, its starting position and end position in the previous off target arrays -> finally not used

		# Load target file, remove overlapping regions, sort it and load it
		bed=bed_file.bed_file(target)
		sortedBed=bed.my_sort_bed(tmpdir=TMP)
		nonOverlappingBed=sortedBed.non_overlapping_exons(1,tmpdir=TMP) # Base 1!!! # This generates a BED file in base 1 (Non-standard BED)
		finalBed=nonOverlappingBed.my_sort_bed(tmpdir=TMP) # BED file in base 1 (Non-standard BED)
		finalBed.load_custom(-1) # Load chromosome and positions in base 1....(finalBed is in base 1 -> Non-standard BED)					
		

		
		
		# Move along all chromosomes
		for currentChromosome in finalBed.chrs.keys():
			# Get all reads of current chromosome it there are reads in such chromosome
			
			positionArray=[] # Structure that stores bp positions along the current chromosome
			coverageArray=[] # Structure that stores coverage values for a related position in the chromosome (positionArray)


			if(currentChromosome in sortedBam.references):					
				if(sortedBam.count(str(currentChromosome))>0): # There might me information about the chromosome in the BAM header but no reads
					allReads=sortedBam.fetch(str(currentChromosome))				
				
					initPositions=[] # Structure that stores initial positions of each read
					endPositions=[] # Structure that stores end positions of each read
					for currentRead in allReads:
						if not currentRead.cigarstring:
							continue
						initPositions.append(int(currentRead.pos)+1) # Fetch is 0-base indexing!!! We are working on 1-base indexing 
						endPositions.append(int(currentRead.aend))
					
					# Select a given number of reads according to numReadsDesired
					if(numberReads==None):
						numReadsDesired=sortedBam.nreads()
					else:
						numReadsDesired=numberReads
						
					numpy.random.seed(1)
					selected = numpy.random.uniform(size=len(initPositions))<=(numReadsDesired*1.0/sortedBam.nreads())
				
					# Convert to numpy arrays
					initPositions=numpy.array(initPositions,dtype=numpy.uint32)
					endPositions=numpy.array(endPositions,dtype=numpy.uint32)
				
					selectedPositions_init=initPositions[selected]
					selectedPositions_end=endPositions[selected]
					selectedPositions_end+=1 # At the end of the position, the coverage counts. Sum 1 to say that at this position the coverage decreases
				
					# Sort each vector independtly	
					selectedPositions_init.sort()
					selectedPositions_end.sort()
					
					totalPositions_init=len(selectedPositions_init)
					totalPositions_end=len(selectedPositions_end)
					
					#widgets = ['Examination of reads in chromosome '+ str(currentChromosome)+' ', progressbar.Percentage(), ' ', 
					#		progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
					#pbar = progressbar.ProgressBar(widgets=widgets, maxval=totalPositions_init+totalPositions_end).start() 
					pbarIter=0
#					print 'Examination of reads in chromosome '+ str(currentChromosome)
					
					
					# First iteration is done manually	
					positionArray.append(selectedPositions_init[0]) # Init position of the first read
					coverageArray.append(1)# A single read (coverage=1)
					
					indexInit=1 # Controls index along selectedPositions_init (fist position in this array has been read)
					indexEnd=0 # Controls index along selectedPositions_end
		
					while indexEnd<totalPositions_end:	# While there are reads to be visited	
						if(indexInit< totalPositions_init): # There are still reads that have to be visited
							if(selectedPositions_init[indexInit]<selectedPositions_end[indexEnd]): # If current position in init is smaller than current position in end
								position=selectedPositions_init[indexInit]
								indexInit+=1
								partialSum=1
							elif(selectedPositions_init[indexInit]>selectedPositions_end[indexEnd]): # If current position in end is greater than current position in init
								position=selectedPositions_end[indexEnd]
								indexEnd+=1
								partialSum=-1
							else: # If current position in init is equal to the position in end
								position=selectedPositions_end[indexEnd]
								indexInit+=1
								indexEnd+=1
								partialSum=0									
						else: # All starting positions for all reads have been visited
								position=selectedPositions_end[indexEnd]
								indexEnd+=1
								partialSum=-1				
												
						# Check whether position is already in the vector of positions
						if(position==positionArray[-1]): # More than a read starts or ends at the same time
							coverageArray[-1]+=partialSum					
						elif(partialSum!=0): # If partialSum==0, then a read ends and a new read start -> do not update information
							positionArray.append(position)
							coverageArray.append(coverageArray[-1]+partialSum)
						
						pbarIter=indexInit+indexEnd	
						#pbar.update(pbarIter)
					#pbar.finish()
					
					
					# Transform positionArray and coverageArray to numpy arrays to save memory
					positionArray=numpy.array(positionArray,dtype=numpy.uint32)
					coverageArray=numpy.array(coverageArray,dtype=numpy.uint16)
													
					numPositions=len(positionArray)
					
					# Create a bedgraph with coverage per position for all reads contained in the current chromosome
					if(bedGraphFile!=None):												
						extension='.'
						onlyName=os.path.basename(bedGraphFile)
						components=onlyName.split(extension)
						prefixFile=extension.join(components[:-1])
						newFileName=prefixFile+'.'+str(currentChromosome)+'.'+components[-1]
						if(len(os.path.dirname(bedGraphFile))>0):
							newFileNameFULL=os.path.dirname(bedGraphFile)+'/'+newFileName
						else:
							newFileNameFULL=newFileName
						
						fdw_bedGraph=file(newFileNameFULL,'w')
						fdw_bedGraph.write('track type=bedGraph name="coverage_chr_'+str(currentChromosome)+self.filename+'" description="coverage per position for '+self.filename+' chromosome '+str(currentChromosome)+'"\n')
															
						positionArray_bedGraph=positionArray-1 # BED GRAPH displays in base 1, although data must use 0-base indexing (http://genome.ucsc.edu/FAQ/FAQtracks#tracks1). See format of bedgraph http://genome.ucsc.edu/goldenPath/help/bedgraph.html
						coverageArray_bedGraph=coverageArray					
						for index in range(0,len(positionArray_bedGraph)-1):
							if(coverageArray_bedGraph[index]!=0):
								fdw_bedGraph.write(str(currentChromosome)+' '+str(positionArray_bedGraph[index])+' '+str(positionArray_bedGraph[index+1])+' '+str(coverageArray_bedGraph[index])+'\n')							
						
						fdw_bedGraph.close()
				else:
					numPositions=0
			else: # There are no reads for the current chromosome				
				numPositions=0

			
			# Load exons and get the positions in positionArray that intersects with the regions (exons)	
			#widgets = ['Intersection of reads with BED in chromosome '+ str(currentChromosome)+' ', progressbar.Percentage(), ' ', 
			#			progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
			#pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(finalBed.chrs[str(currentChromosome)])).start() 
#			print 'Intersection of reads with BED in chromosome '+ str(currentChromosome) 
			pbarIter=0
			
						
			indexCoverage=0
			nextPositionChange=0
			previousIndexCoverage=indexCoverage
			firstPosition=len(positionArray_intersect) # First position in the intersection vectors (positionArray_intersect and coverageArray_intersect)
			#firstPositionOffTarget=len(coverageArray_offTarget) # Finally not used (off target reads)


				
			for currentExon in finalBed.chrs[str(currentChromosome)]:
				exon_init=int(currentExon[0])
				exon_end=int(currentExon[1])
				
					
				# Search for the first position that matches with the beginning of currentExon
				thereAreReads=False
				#firstIndexCoverage=indexCoverage		
				readsAfterExonBeg=False # Whether there are reads after the starting position of the exon ##NEW
				
				# Check if change in coverage from the previous iteration is before exon_init. If not, move one position back so that previous coverage value is recovered
				if(indexCoverage >0 and indexCoverage<numPositions and positionArray[nextPositionChange]>exon_init):
					indexCoverage=indexCoverage-1
				
				if(indexCoverage<numPositions and positionArray[indexCoverage]>exon_init and positionArray[indexCoverage]<=exon_end): # There are reads in the exon, but the first read starts after the beginning of the exon 
					startPosition=indexCoverage
					readsAfterExonBeg=True 
				else: # Reads start before the beginning of the exon
					while(indexCoverage<numPositions and positionArray[indexCoverage]<=exon_init):
						indexCoverage+=1			
					startPosition=indexCoverage-1
				
				#if(firstIndexCoverage!=startPosition): # Store coverage of read mapped out of exons (finally not used -> off target reads)
				#	positionArray_offTarget.extend(positionArray[firstIndexCoverage:(startPosition)]) # Interested in positions below startPosition
				#	coverageArray_offTarget.extend(coverageArray[firstIndexCoverage:(startPosition)])
				#	# Except for the last exon, reads that map after curentExon will be visited in the next exon (positions before the exon)
				
				if(indexCoverage!=previousIndexCoverage):
					thereAreReads=True
				# Search for the first position that matches with the end of currentExon
				while(indexCoverage<numPositions and positionArray[indexCoverage]<=exon_end):
					indexCoverage+=1
				endPosition=indexCoverage-1
				
				nextPositionChange=indexCoverage # positionArray[indexCoverage] contains next position after exon_end in which coverage changes
				
				
				if(endPosition!=startPosition or thereAreReads):		# Puede darse el caso de que un exon este en una sola coordenada de positionArray -> tenerlo en cuenta!!!!	
					numElements=endPosition-startPosition+1 # Num of coverage/bases values to be interseted		
					if(writeToFile==None):		
						if(readsAfterExonBeg):
							positionArray_intersect.extend([exon_init])
							coverageArray_intersect.extend([0])	
																	
						positionArray_intersect.extend(positionArray[startPosition:(endPosition+1)])
						coverageArray_intersect.extend(coverageArray[startPosition:(endPosition+1)])
					
						# It may happen that positionArray[startPosition] and positionArray[endPosition+1] are not equal to the beginning and end of the exon respectively
						# Thus, it is forced to modify positionArray_intersect[indexIntersect] and positionArray_intersect[indexIntersect_end-1] to the beggining and end of the exon respectively
						if(not readsAfterExonBeg):
							positionArray_intersect[-numElements]=exon_init
						#positionArray_intersect[-1]=exon_end # No es valido si la posicion final del exon tiene el mismo coverage que otra posicion anterior									
					else:# Write to file base per base current exon
						keys=range(exon_init,exon_end+1)
						dicExon=dict(zip(keys,[-1]*len(keys)))
						
						currentPositionArray=positionArray[startPosition:(endPosition+1)]
						if(not readsAfterExonBeg): 
							currentPositionArray[0]=exon_init
						else: # If reads after after exon_init, place a zero coverage to the beginning of the exon
							dicExon[exon_init]=0							
						
						currentCoverageArray=coverageArray[startPosition:(endPosition+1)]
						dicExon.update(zip(currentPositionArray,currentCoverageArray)) # Fill positions where coverage changes
						
						previousKey=sorted(dicExon.keys())[0]
						
						for currentKey in sorted(dicExon.keys()):
							if(dicExon[currentKey]==-1): # Take into account the coverage of the last position (key) that has coverage
								#dicExon[currentKey]=dicExon[previousKey]
								fdw.write(str(currentChromosome)+'\t'+str(exon_init)+'\t'+str(exon_end)+'\t'+str(dicExon[previousKey])+'\n')
							else:
								fdw.write(str(currentChromosome)+'\t'+str(exon_init)+'\t'+str(exon_end)+'\t'+str(dicExon[currentKey])+'\n')					
								previousKey=currentKey
																																			
					previousIndexCoverage=indexCoverage
				else: # A exon is not covered -> startPosition moves until the end but endPosition could not move more. However, we want to point that, although there are no reads for
					# the current exon, the initial position of the exon and its coverage (zero) are stored
					if(writeToFile==None):
						positionArray_intersect.append(exon_init) # Initial position of exon
						coverageArray_intersect.append(0) # Zero coverage
					else:
						keys=range(exon_init,exon_end+1)
						dicExon=dict(zip(keys,[0]*len(keys)))
						
						for currentKey in sorted(dicExon.keys()):
							fdw.write(str(currentChromosome)+'\t'+str(exon_init)+'\t'+str(exon_end)+'\t'+str(dicExon[currentKey])+'\n')
						
					
					indexCoverage=previousIndexCoverage # Move backwards in positionArray to search reads for the next exon (all exons have not visited yet)
			
				
				#pbarIter+=1	
				#pbar.update(pbarIter)
				
			# Check if there are reads after the end of the last exon (finally not used -> off target reads
			#if(endPosition!=(numPositions-1)):
			#	positionArray_offTarget.extend(positionArray[endPosition+1:(numPositions)]) 
			#	coverageArray_offTarget.extend(coverageArray[endPosition+1:(numPositions)])
	
			#pbar.finish()	
	
			# Get the positions in the intersection vector for the current chromosome
			chromosomeCoordinates[currentChromosome]=(firstPosition,len(positionArray_intersect)-1)
			#chromosomeCoordinatesOffTarget[currentChromosome]=(firstPositionOffTarget,len(positionArray_offTarget)-1)# Finally not used (off target reads))
			
			#chromosomeCoordinatesOffTarget[currentChromosome]=str(firstPositionOffTarget)+'-'+str(len(positionArray_offTarget)-1) 
		
			del positionArray
			gc.collect()
	
			del coverageArray
			gc.collect()
		
		positionArray_intersect=numpy.array(positionArray_intersect,dtype=numpy.uint32)
		coverageArray_intersect=numpy.array(coverageArray_intersect,dtype=numpy.uint16)
		
		if(writeToFile!=None):
			fdw.close()
			
		if(executiongranted<>None):
			executiongranted.release()
		
		return [positionArray_intersect,coverageArray_intersect,chromosomeCoordinates,finalBed]
	
	
	
	
		
	def myReadsOnTarget(self,target):
		"""*******************************************************************************************************************************************
		JPFLORIDO
		Task:  Custom method equivalent to bedtool's intersect bed. It is made to count the number of reads on/off target and the number of reads on/off target that start and end at the same position
		Inputs: BAM file (self) and target file 
		Outputs: 
			dicOnTarget: dictionary containing the number of reads on target per chromosome
			dicOffTarget: dictionary containing the number of reads off target per chromosome
			dicTotal: dictionary containing the number of total reads (on + off target) per chromosome
			duplicatesOnTarget: array containing number of reads on target that start and end at the same position at different levels of reads per start/end position; duplicatesOnTarget[0] = number of reads that start/end at a unique position;  duplicatesOnTarget[1] = number of reads that start and end at the same position (at most two reads per start/end positions); uplicatesOnTarget[n] = no.of reads that start and end at the same position (at most n+1 reads per start/end positions)
			duplicatesOffTarget: array containing number of reads off target that start and end at the same position at different levels of reads per start/end position; duplicatesOffTarget[0] = number of reads that start/end at a unique position;  duplicatesOffTarget[1] = number of reads that start and end at the same position (at most two reads per start/end positions); uplicatesOffTarget[n] = no.of reads that start and end at the same position (at most n+1 reads per start/end positions)
			
			
		Other issues: if a read targets two regions in the BED file, the read is count twice. Pysam works in 0-base indexing, so, the overlapping regions 
		# will be removed from the target file as they are -> no conversion to real zero or one-base indexing
		*******************************************************************************************************************************************"""

		global TMP
		
		
		# Load target file, remove overlapping regions, sort it and load it
		bed=bed_file.bed_file(target)
		sortedBed=bed.my_sort_bed(tmpdir=TMP)
		nonOverlappingBed=sortedBed.non_overlapping_exons(1,tmpdir=TMP) # Base-1 indexing
		finalBed=nonOverlappingBed.my_sort_bed(tmpdir=TMP) # BED file in base 1 (Non-standard BED)
		finalBed.load_custom(-1) # Load chromosome and positions as they are				
	

		dicOnTarget={}
		dicOffTarget={}
		dicTotalReads={}
		countDuplicatesOnTarget={}
		countDuplicatesOffTarget={}
		previousChromosome='0' # False chromosome
		dicOnTargetChr={}
		dicOffTargetChr={}
		readsOnTarget=0 # fjavier: Overall number of reads on target 
		readsOnTargetChr=0
		readsOffTargetChr=0
		thereisInfo=False
		duplicatesOnTarget=[] # Contains reads that appear once, two times, three times...on target
		duplicatesOffTarget=[] # Contains reads that appear once, two times, three times...off target
		

		
		for currentRead in self:
			# Get chromosome info from currentRead
			if not currentRead.cigarstring:
				continue
			try:
				currentChromosome=self.getrname(currentRead.tid)
			except ValueError:
				if(currentRead.tid<0):
					print "\nPlease, check that BAM file has only mapped reads"
					exit()
			if(currentChromosome != previousChromosome):
#				print "Examining Chromosome "+str(currentChromosome)+'... \n'											
				if(thereisInfo): # A new chromosome is started. Store information about previous chromosome
#					print "Saving information of chromosome "+ str(previousChromosome)+'\n'
					readsOnTarget += readsOnTargetChr
					dicOnTarget[str(previousChromosome)]=readsOnTargetChr # Count reads on target for the current chromosome
					dicOffTarget[str(previousChromosome)]=readsOffTargetChr # Count reads off target for the current chromosome
					dicTotalReads[str(previousChromosome)]=readsOnTargetChr+readsOffTargetChr # Count all reads for the current chromosome
					
					# Get counts of each read for the current chromosome and accumulate them in two dictionaries
					countReadsOnTarget=numpy.array(dicOnTargetChr.values())
					if(len(countReadsOnTarget)>0):
						for count in range(min(countReadsOnTarget),max(countReadsOnTarget)+1):				
							if(countDuplicatesOnTarget.has_key(count)):
								countDuplicatesOnTarget[count]+=len(numpy.where(countReadsOnTarget==count)[0])		
							else:
								countDuplicatesOnTarget[count]=len(numpy.where(countReadsOnTarget==count)[0])
							
						
					countReadsOffTarget=numpy.array(dicOffTargetChr.values())
					if(len(countReadsOffTarget)>0):
						for count in range(min(countReadsOffTarget),max(countReadsOffTarget)+1):				
							if(countDuplicatesOffTarget.has_key(count)):
								countDuplicatesOffTarget[count]+=len(numpy.where(countReadsOffTarget==count)[0])
							else:
								countDuplicatesOffTarget[count]=len(numpy.where(countReadsOffTarget==count)[0])
			
					del dicOnTargetChr
					gc.collect()

					del dicOffTargetChr
					gc.collect()
				
				if(currentChromosome in finalBed.chrs.keys()): # Check whether chromosome related to the read is in the BED file
					targets_chr=numpy.array(finalBed.chrs[str(currentChromosome)]) # all rows, 1st column -> exon_init positions; all rows, 2nd column -> exon_end positions
					numRegions=len(targets_chr[:,0])
									
				dicOnTargetChr={}
				dicOffTargetChr={}
				readsOnTargetChr=0
				readsOffTargetChr=0
				startPos=0
				previousChromosome=currentChromosome
				thereisInfo=False
			
			
			if(currentChromosome in finalBed.chrs.keys()):	# Go through the process only if currentchromosome is present in the region
				thereisInfo=True
				if not currentRead.cigarstring:
					continue
				currentReadInit=int(currentRead.pos)+1 # BED has been transformed to base 1
				currentReadEnd=int(currentRead.aend)					
				currentPos=startPos
				stop=False
	
				while((not stop) and currentPos<numRegions):
					if((currentReadInit>=targets_chr[currentPos,0] and currentReadInit<=targets_chr[currentPos,1]) or
					(currentReadEnd>=targets_chr[currentPos,0] and currentReadEnd<=targets_chr[currentPos,1])): # Read on target: read init inside the region or read end inside the region
	
						if(dicOnTargetChr.has_key((currentReadInit,currentReadEnd))): # Read on target
							dicOnTargetChr[(currentReadInit,currentReadEnd)]+=1
						else:
							dicOnTargetChr[(currentReadInit,currentReadEnd)]=1
						
						readsOnTargetChr+=1	
	
						stop=True
						startPos=currentPos								
					else: # Check to the next region, unless a read ends before the current region
						if(currentReadEnd<targets_chr[currentPos,0]): # A read ends before the current region -> stop searching. Read off target
							if(dicOffTargetChr.has_key((currentReadInit,currentReadEnd))):
								dicOffTargetChr[(currentReadInit,currentReadEnd)]+=1
							else:
								dicOffTargetChr[(currentReadInit,currentReadEnd)]=1
		
							readsOffTargetChr+=1
							stop=True
						else:
							currentPos+=1
				
				if(currentPos==numRegions): # Read off target -> we moved along all regions and none of them mathes with the read			
					if(dicOffTargetChr.has_key((currentReadInit,currentReadEnd))):
						dicOffTargetChr[(currentReadInit,currentReadEnd)]+=1
					else:
						dicOffTargetChr[(currentReadInit,currentReadEnd)]=1
	
					readsOffTargetChr+=1
			
						
		if(thereisInfo):					
			# Store information for the last chromosome		
#			print "Saving information of chromosome "+ str(currentChromosome)+'\n'
			readsOnTarget += readsOnTargetChr
			dicOnTarget[str(currentChromosome)]=readsOnTargetChr # Count reads on target for the current chromosome
			dicOffTarget[str(currentChromosome)]=readsOffTargetChr # Count reads off target for the current chromosome
			dicTotalReads[str(currentChromosome)]=readsOnTargetChr+readsOffTargetChr # Count all reads for the current chromosome
		
			# Get counts of each read for the current chromosome and accumulate them in two dictionaries
			countReadsOnTarget=numpy.array(dicOnTargetChr.values())
			if(len(countReadsOnTarget)>0):
				for count in range(min(countReadsOnTarget),max(countReadsOnTarget)+1):				
					if(countDuplicatesOnTarget.has_key(count)):
						countDuplicatesOnTarget[count]+=len(numpy.where(countReadsOnTarget==count)[0])		
					else:
						countDuplicatesOnTarget[count]=len(numpy.where(countReadsOnTarget==count)[0])
				
			
			countReadsOffTarget=numpy.array(dicOffTargetChr.values())
			if(len(countReadsOffTarget)>0):
				for count in range(min(countReadsOffTarget),max(countReadsOffTarget)+1):				
					if(countDuplicatesOffTarget.has_key(count)):
						countDuplicatesOffTarget[count]+=len(numpy.where(countReadsOffTarget==count)[0])
					else:
						countDuplicatesOffTarget[count]=len(numpy.where(countReadsOffTarget==count)[0])
			
		# Store information about duplicates in an array		
		
		if(len(countDuplicatesOnTarget.keys())>0):
			duplicatesOnTarget=numpy.zeros(max(countDuplicatesOnTarget.keys()))
			for currentKey in sorted(countDuplicatesOnTarget.keys()):
				duplicatesOnTarget[currentKey-1]=countDuplicatesOnTarget[currentKey]*currentKey
				# It is multiplied by currentKey, since countDuplicatesOnTarget[currentKey] has the number of start/end positions that have currentKey reads. We want the number of total reads.
		
		if(len(countDuplicatesOffTarget.keys())>0):
			duplicatesOffTarget=numpy.zeros(max(countDuplicatesOffTarget.keys()))
			for currentKey in sorted(countDuplicatesOffTarget.keys()):
				duplicatesOffTarget[currentKey-1]=countDuplicatesOffTarget[currentKey]*currentKey


		# Check whether there are no reads for a given chromosome contained in BED file
		
		for currentChromosome in finalBed.chrs.keys():
			if(not currentChromosome in dicOnTarget.keys()): # If chromosome is not in the dictionary is because there are no reads for that chromosome
				dicOnTarget[str(currentChromosome)]=0			
				dicTotalReads[str(currentChromosome)]=0
		#FJAVIER: dicOffTarget is not actually needed yet, that is why it is not being returned
		return [readsOnTarget,dicOnTarget,dicTotalReads,duplicatesOnTarget,duplicatesOffTarget]
			
	
	
	
	
	
	def coverageperbase(self,currentChromosome):
		"""*******************************************************************************************************************************************
		Task:  almost the same code as the first part of myCoverageBed. Gets the coverage per position for a given chromosome without focussing on
			a given target.
		Inputs:
			currentChromosome: string containing the identifier of the chromosome that will be inspected. 
		Outputs: 
			positionArray: a vector of bp coordinates in which coverage changes according to the reads of this BAM
			coverageArray: coverage (number of reads) for a given bp position 
		Other issues: it is intented to return also positions off target. However, coverage around exons changes a lot, so there is an important increasing
		# in the memory used 
		*******************************************************************************************************************************************"""

		pid = str(os.getpid())
		
		# Check whether BAM file is sorted
		command='samtools view -H '+ self.filename+' | grep SO:'
		fd=os.popen(command)
		outputCommand=fd.read()
		removetmp = False
		
		# Sort bed in case it is not sorted		
		if ((fd.close()<>None) or('coordinate' not in outputCommand.split('SO:')[-1])): # BAM not sorted -> Sort and indexing BAM
			sortedBam=self.sort_bam()
			removetmp = True
		else:
			sortedBam=self
							
		positionArray=[] # Structure that stores bp positions along the current chromosome
		coverageArray=[] # Structure that stores coverage values for a related position in the chromosome (positionArray)

		# Get all reads of current chromosome if there are reads in such chromosome
		if(currentChromosome in sortedBam.references):	
			allReads=sortedBam.fetch(str(currentChromosome))
		
			initPositions=[] # Structure that stores initial positions of each read
			endPositions=[] # Structure that stores end positions of each read
			for currentRead in allReads:
				initPositions.append(int(currentRead.pos)+1) # Fetch is 0-base indexing!!! We are working on 1-base indexing 
				endPositions.append(int(currentRead.aend))
			
			# Convert to numpy arrays
			initPositions=numpy.array(initPositions,dtype=numpy.uint32)
			endPositions=numpy.array(endPositions,dtype=numpy.uint32)				
			endPositions+=1 # At the end of the position, the coverage counts. Sum 1 to say that at this position the coverage decreases
		
			# Sort each vector independtly	
			initPositions.sort()
			endPositions.sort()
			
			totalPositions_init=len(initPositions)
			totalPositions_end=len(endPositions)

			# This happens if there are no reads in this chromosome
			if(totalPositions_init+totalPositions_end == 0):
				print 'ERROR: no reads found for chr '+currentChromosome+'. However, this chromosome id was found in the header of the bam file, indicating that there have to be reads in it.'
				print '	The problem is most probably due to a incorrect .bai file. Please, regenerate the index and try again.'
				sys.exit(1)
				
#			widgets = ['Examination of reads in chromosome '+ str(currentChromosome)+' ', progressbar.Percentage(), ' ', 
#					progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
#			pbar = progressbar.ProgressBar(widgets=widgets, maxval=totalPositions_init+totalPositions_end).start()
			print 'Examining reads at chromosome '+str(currentChromosome)
			pbarIter=0
			
			
			# First iteration is done manually
			if(initPositions[0]>1):
				positionArray.append(1) # Init position of the first read
				coverageArray.append(0)# A single read (coverage=1)					
			positionArray.append(initPositions[0]) # Init position of the first read
			coverageArray.append(1)# A single read (coverage=1)
			
			indexInit=1 # Controls index along initPositions (fist position in this array has been read)
			indexEnd=0 # Controls index along endPositions

			while indexEnd<totalPositions_end:	# While there are reads to be visited	
				if(indexInit< totalPositions_init): # There are still reads that have to be visited
					if(initPositions[indexInit]<endPositions[indexEnd]): # If current position in init is smaller than current position in end
						position=initPositions[indexInit]
						indexInit+=1
						partialSum=1
					elif(initPositions[indexInit]>endPositions[indexEnd]): # If current position in end is greater than current position in init
						position=endPositions[indexEnd]
						indexEnd+=1
						partialSum=-1
					else: # If current position in init is equal to the position in end
						position=endPositions[indexEnd]
						indexInit+=1
						indexEnd+=1
						partialSum=0									
				else: # All starting positions for all reads have been visited
						position=endPositions[indexEnd]
						indexEnd+=1
						partialSum=-1				
										
				# Check whether position is already in the vector of positions
				if(position==positionArray[-1]): # More than a read starts or ends at the same time
					coverageArray[-1]+=partialSum					
				elif(partialSum!=0): # If partialSum==0, then a read ends and a new read start -> do not update information
					positionArray.append(position)
					coverageArray.append(coverageArray[-1]+partialSum)
				
				pbarIter=indexInit+indexEnd	
#				pbar.update(pbarIter)
#			pbar.finish()
			
			
			# Transform positionArray and coverageArray to numpy arrays to save memory
			positionArray=numpy.array(positionArray,dtype=numpy.uint32)
			coverageArray=numpy.array(coverageArray,dtype=numpy.uint16)
											
			numPositions=len(positionArray)
					
		else: # There are no reads for the current chromosome				
			numPositions=0
		
		if(removetmp):
			os.remove(sortedBam.filename)
			
		return [positionArray,coverageArray]





	def coveragetrack(self, bedfilename, windowsize, overlap, fileout):
		"""************************************************************************************************************************************************************
		Task: generates a track of coverage for this file and a given target.		 
		Inputs:
			bedfilename: String containing the full path to a bed file containing the target region.
			windowsize: integer containing the size of the windows that will be used to promediate coverage.
			overlap: integer containing the size of each step when moving the window through the chromosome.
			fileout: string containing the full path to the bedgraph file where the track will be stored.
		 Ouputs: a bedgraph file named fileout will be created containing the track of coverage.
		************************************************************************************************************************************************************"""
		
		rawbed = bed_file.bed_file(bedfilename)
		sortedBed=rawbed.my_sort_bed() # Sort bed avoiding bedtools
		bed=sortedBed.non_overlapping_exons(1) # Cargar BED tal cual (standard BED)
		bed.load_custom(-1) # Load chromosome and positions as they are (standard BED)	

		counts = bed.getWindows(windowsize,overlap)
		
		# Initializes all counts to 0
		coords = bed.sortcountkeys(counts)
		previouschr = None
		for coord in counts:
			counts[coord] = 0	
		
		# Calculates mean coverage within each window		
		i=0
		while(i<len(coords)):
			currchr = coords[i][0]		
			# Calculate coverage for current chromosome
			if(currchr<>previouschr):
#				PREGUNTARLE A JAVI EN Q BASE VIENEN DADOS ESTOS DOS VECTORES
				changepositions,coverage = self.coverageperbase(currchr)
				previouschr = currchr
				coverageidx = 0				

			if(coverageidx<len(changepositions)):
				# Move the pointer until the coordinate in changepositions is within current window
				# Check if the pointer is below the lower limit of the window				
				if(changepositions[coverageidx]<=coords[i][1]):
					# Increase the pointer until it falls within the window. If coverage changes in the first base of the window, the value
					# is saved to be used later
					try:				
						while(coverageidx<len(changepositions) and changepositions[coverageidx]<=coords[i][1]):
							lastcoverage = coverage[coverageidx]
							coverageidx += 1
					except IndexError:
						print 'len(changepositions) = '+str(len(changepositions))
						print 'coverageidx = '+str(coverageidx)
						print 'i = '+str(i)
						print 'len(coords) = '+str(len(coords))
						print 'changepositions[coverageidx-1] = '+str(changepositions[coverageidx-1])
						print 'coords[i][1] = '+str(coords[i][1])
						sys.exit(1)
				# Check if the pointer is above the lower limit of the window. This happens if previous window overlaps with current window					
				elif(changepositions[coverageidx]>coords[i][1]):
					# Decreases the pointer until reaching the lower limit of the window				
					while(changepositions[coverageidx]>coords[i][1]): #ME SALTA UN ERROR AQUI AL LLAMAR DESDE BAM_UTILS.aggregatedcoveragetrack
						coverageidx -= 1
					# Save the "first" coverage of the window					
					lastcoverage = coverage[coverageidx]
					# Move pointer to the first change above the lower limit of the window				
					coverageidx+=1
				
				# Check whether we are already at the end of the changepositions vector. In other words, there are no more change positions for the remaining
				# windows in coords.
				if(coverageidx<len(changepositions)):					
					# Pass through each interval of constant coverage summing up the corresponding coverage of each base								
					lastbase = coords[i][1]
					while(coverageidx<len(changepositions) and changepositions[coverageidx]<=coords[i][2]):
						counts[coords[i]] += (changepositions[coverageidx]-lastbase)*lastcoverage
						lastbase = changepositions[coverageidx]
						lastcoverage = coverage[coverageidx]
						coverageidx += 1
	
					counts[coords[i]] += (coords[i][2]-lastbase+1)*lastcoverage
					counts[coords[i]] = counts[coords[i]]*1.0/(coords[i][2]-coords[i][1]+1)
				else:
					counts[coords[i]] = 0

					# The piece of code in this 'else' is only run when windows span beyond coverage limits in 'changepositions'. It would
					# be really improbable that the last coverage change is different from 0.
					if(lastcoverage<>0):
						print 'WARNING: please, check that all data is correct at bam_file.coveragetrack.'
						print '	Execution stopped.'
						sys.exit(1)
					
			else:
				counts[coords[i]] = 0
				
				# The piece of code in this 'else' is only run when windows span beyond coverage limits in 'changepositions'. It would
				# be really improbable that the last coverage change is different from 0.
				if(lastcoverage<>0):
					print 'WARNING: please, check that all data is correct at bam_file.coveragetrack.'
					print '	Execution stopped.'
					sys.exit(1)
					
			i += 1
						
		bed.windows2bedgraph(counts, fileout, startshift=-1)	
				
		
		
		

	def coveragetrackchr(self, chr, bedfilename, windowsize=300, overlap=100):
		"""************************************************************************************************************************************************************
		Task: generates a track for a given chromsome containing genomic positions and coverage in the form of a list of tuples.		 
		Inputs:
			chr: string containing the chromsome identifier.
			bedfilename: String containing the full path to a bed file indicating target regions.
			windowsize: integer containing the size of the windows that will be used to promediate coverage.
			overlap: integer containing the size of each step when moving the window through the chromosome.
		 Ouputs: a list of tuples of the form (position, coverage), both integers. 
		************************************************************************************************************************************************************"""
		
		rawbed = bed_file.bed_file(bedfilename)
		sortedBed=rawbed.my_sort_bed() # Sort bed avoiding bedtools
		bed=sortedBed.non_overlapping_exons(1) # Cargar BED tal cual (standard BED)
		bed.load_custom(-1) # Load chromosome and positions as they are (standard BED)	

		counts = bed.getWindows(windowsize,overlap,chrs=[chr])				
		coords = bed.sortcountkeys(counts)
		
		# Initializes each window count to 0
		for coord in counts:
			counts[coord] = 0	

		changepositions,coverage = self.coverageperbase(chr)		
		coverageidx = 0
		# Calculates mean coverage within each window
		for coord in coords:
#				PREGUNTARLE A JAVI EN Q BASE VIENEN DADOS ESTOS DOS VECTORES
			# Move the pointer until the coordinate in changepositions is within current window
			# Check if the pointer is below the lower limit of the window	
			if(changepositions[coverageidx]<=coord[1]):
				# Increase the pointer until it falls within the window. If coverage changes in the first base of the window, the value
				# is saved to be used later
				while(changepositions[coverageidx]<=coord[1]):
					lastcoverage = coverage[coverageidx]
					coverageidx += 1
			# Check if the pointer is above the lower limit of the window. This happens if previous window overlaps with current window
			elif(changepositions[coverageidx]>coord[1]):
				# Decreases the pointer until reaching the lower limit of the window
				while(changepositions[coverageidx]>coord[1]):
					coverageidx -= 1
				# Save the "first" coverage of the window
				lastcoverage = coverage[coverageidx]
				# Move pointer to the first change above the lower limit of the window
				coverageidx+=1

			# Pass through each interval of constant coverage summing up the corresponding coverage of each base				
			lastbase = coord[1]
			while(changepositions[coverageidx]<=coord[2]):
				counts[coord] += (changepositions[coverageidx]-lastbase)*lastcoverage
				lastbase = changepositions[coverageidx]
				lastcoverage = coverage[coverageidx]
				coverageidx += 1
			counts[coord] += (coord[2]-lastbase+1)*lastcoverage
			counts[coord] = counts[coord]*1.0/(coord[2]-coord[1]+1)
				
	
					
		return [((coord[1]+coord[2])/2.0,counts[coord]) for coord in coords]	





	def coveredpositions(self, chr, threshold):
		chrinfo = self.header['SQ'][0]
		i = 1
		while(i<len(self.header['SQ']) and chrinfo['SN']<>chr):
			chrinfo = self.header['SQ'][i]
			i+=1
			
		if(chrinfo['SN']<>chr):
			print 'ERROR: '+self.filename+' does not contain mapped reads on contig '+chr
			print '	Exiting.'
			sys.exit(1)
			
		chrsize = int(chrinfo['LN'])
		
		changepositions,coverage = self.coverageperbase(chr)
		
		print 'Analyzing coverage per base...'
		fullcoverage = []
		previousbase=1
		lastcoverage=0
		for i,base in enumerate(changepositions):
			if(lastcoverage>=threshold): 
				fullcoverage += [True for j in range(previousbase, base)]
			else:
				fullcoverage += [False for j in range(previousbase, base)]
			
			previousbase = base
			lastcoverage = coverage[i]
			
		if(lastcoverage>=threshold): 
			fullcoverage += [True for j in range(previousbase, chrsize+1)]
		else:
			fullcoverage += [False for j in range(previousbase, chrsize+1)]			
		print '	Done.'
			

		return numpy.array(fullcoverage)
		
	
	
	
			
	def nonofftarget(self,bed):
		"""************************************************************************************************************************************************************
		Task: returns reads on/off target
		Inputs:
			bed: bed file with capture coordinates
		************************************************************************************************************************************************************"""
	
		capture=bed_file.bed_file(bed)
		bamOT=capture.intersectbam(self, "ontarget.bam")
		tread=self.nreads()
		nread=bamOT.nreads()
		
		return [nread, tread-nread]
		
						
		
		
		
		
		
	def generate_ontarget_results(self, bamlist, nread, tread, onperchr, totalperchr, enrichment, percontarget, outdir, legend):
		"""************************************************************************************************************************************************************
		Task: this method is dependant on "reads_on_target". Generates a graph and an xls file that summarize on/off target read results.
		Inputs:
			bamlist: list of bam_file objects.
			nread: list of integers, each containing the number of on-target reads in the corresponding bam file.
			tread: list of integers, each containing the total number of reads in the corresponding bam file.			
			onperchr: list of dictionaries. There must be one dictionary for each bam file in bamlist. Each dictionary contains data about the number of on-target 
				reads per chromosome. Keys of the dictionary are contig ids (e.g. 'chr1'). Elements are the number of on-target reads mapped in that contig.
			totalperchr: list of dictionaries. There must be one dictionary for each bam file in bamlist. Each dictionary contains data about the total number of 
				reads mapped in each chromosome. Keys of the dictionary are contig ids (e.g. 'chr1'). Elements are the total number of reads mapped in that contig.
			enrichment: list or multiprocessing.Array objecto containing the enrichment value for each bam (on-target reads per Kb)/(off target reads per Kb)
			percontarget: list or multiprocessing.Array object containing the percentaje of reads on target for each bam file
			outdir: string containing the full path to the directory were result files will be saved.
			legend: list of strings containing the label to tag each bam file in the xls and png files that will be created.
		Outputs: a new bar plot will be created named outdir/reads_on_target.png. This is a bar plot which indicates for each bam, the percentaje of on-target
			reads in each chromosome. In addition, a new outdir/reads_on_target.xls file will be created containing the numbers used to build the bar plot.					
		************************************************************************************************************************************************************"""
		
		contigids = onperchr[0].keys()
		contigids.sort()

		colours = ['#46a246', '#ff0000', '#00ff00', '#0000ff', '#cc0011', '#007722', '#110066', '#c1c1c1', '#544db1', '#aa5198', '#bbd1e9', '#f1c4ab', '#24687a']		
		fig=pyplot.figure(figsize=(13,7))
		ax = fig.add_subplot(111)
		ind = numpy.arange(2*len(contigids), step=2)			

		# Create a set of bars for each bam file		
		width = 0.7	   # the width of the bars
		rects = []
		tick_pos = 0
		for i in range(len(onperchr)):
			# Calculate on-target percentajes per contig (chr)
			y = []
			for chr in contigids:
				# If there are no reads in current contig, let y be the lowest value (0) to avoid division by 0
				if(totalperchr[i][chr]>0):
					y.append(onperchr[i][chr]*100.0/totalperchr[i][chr])
				else:
					y.append(0)
			rects.append(ax.bar(ind+tick_pos*width, y, width, color=colours[i%12]))
			tick_pos += 1
		
		ax.set_xlim(0)
		ax.set_xticks(ind+(len(bamlist)*width)/2.0)
		
		ax.set_xticklabels(map(str, contigids), rotation='vertical')
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().tick_left()
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.set_ylabel("% on-target reads")
		ax.set_xlabel("Chromosomes")

		# Just create the legend if more than one bam file are provided		
		if(len(bamlist)>1):
			
			# If legend is not provided generate an automatic legend
			if(legend==None):
				legend = []
				# The label for each bam will be its basename
				for bam in bamlist:
					legend.append(os.path.basename(bam.filename))
				 	
			# Shink current axis by 20%
			box = ax.get_position()
			ax.set_position([box.x0, box.y0, box.width, box.height*0.9])	
			ax.legend( tuple([rect[0] for rect in rects]), tuple(legend), loc="lower left", bbox_to_anchor=(0,1.03) )		
		
		# Initialize the workbook and sheet
		wb = xlwt.Workbook()
		
		# A sheet is created in the xls for each bam file
		for j,bam in enumerate(bamlist):
			ws = wb.add_sheet(str(j+1)+'-'+legend[j][:10])
			
			# Create header font
			header_style = xlwt.easyxf('font: bold on')
			
			ws.write(0,0,'Enrichment: ', header_style); ws.write(0,1,enrichment[j])		
			ws.write(4,1,'Reads on target', header_style);ws.write(4,2,'Reads off target', header_style);ws.write(4,3,'% reads on target', header_style);
			ws.write(4,4,'% reads off target', header_style);
			ws.write(5,0,'Total', header_style);
			
			percontarget[j] = nread[j]*100.0/tread[j]
			ws.write(5,1,nread[j]);ws.write(5,2,tread[j]-nread[j]);ws.write(5,3,percontarget[j]);ws.write(5,4,(tread[j]-nread[j])*100.0/tread[j]);
			
			# Write 
			for i,chr in enumerate(contigids):
				ws.write((i+6),0,chr, header_style);
				
				# Leave an empty cell if the number of reads mapped in current contig is 0 (avoid division by zero)
				if(totalperchr[j][chr]>0):
					ws.write((i+6),1,onperchr[j][chr]);ws.write((i+6),2,totalperchr[j][chr]-onperchr[j][chr]);ws.write((i+6),3,onperchr[j][chr]*100.0/totalperchr[j][chr]);
					ws.write((i+6),4,(totalperchr[j][chr]-onperchr[j][chr])*100.0/totalperchr[j][chr]);
				else:
					ws.write((i+6),1,onperchr[j][chr]);ws.write((i+6),2,totalperchr[j][chr]-onperchr[j][chr]);ws.write((i+6),3,'');
					ws.write((i+6),4,'');
			
		wb.save(outdir+'/reads_on_target.xls')
		fig.savefig(outdir+'/reads_on_target.png')
		matplotlib.pyplot.close(fig)





	def generate_duplicates_results(self, bamlist, ontargetreads, offtargetreads, onduplicates, offduplicates, outdir, legend):
		"""************************************************************************************************************************************************************
		Task: this method is dependant on "reads_on_target". Generates a graph and an xls file that summarize duplicated read results.
		Inputs:
			bamlist: list of bam_file objects.
			ontargetreads: list of integers, each containing the number of on-target reads in the corresponding bam file.
			offtarget reads: list of integers, each containing the number of off-target reads in the corresponding bam file.
			onduplicates: list of sublists. Each sublist contains data about the number of duplicates in the corresponding bam among the on-target reads. Indexes in each sublist indicate
				the number of replications (index 0 - unique reads, index 1 - duplicated reads, index 2 - triplicated reads, etc.). Integers within each position
				of the sublists indicate the number of reads that appear with such a number of replications.
			offduplicates: list of sublists. Each sublist contains data about the number of duplicates in the corresponding bam among the off-target reads. Indexes in each sublist indicate
				the number of replications (index 0 - unique reads, index 1 - duplicated reads, index 2 - triplicated reads, etc.). Integers within each position
				of the sublists indicate the number of reads that appear with such a number of replications. 				 
			outdir: string containing the full path to the directory were the new files will be saved.
			legend: list of strings containing the label to tag each bam file in the xls and png files that will be created.
		Outputs: a new bar plot will be created named outdir/duplicates.png. This is a bar plot with two bars per number of replications. Green bars indicate numbers
			of replications among on-target reads. Red bars indicate numbers of replications among off-target reads. In addition, a new outdir/duplicates.xls file
			will be created containing the numbers used to build the bar plot.					
		************************************************************************************************************************************************************"""
	
		MAXREPLICATES = 3
		
		# Initialize the workbook and sheet
		wb = xlwt.Workbook()
		header_style = xlwt.easyxf('font: bold on')		
	
		# Processess each bam file
		for j,bam in enumerate(bamlist):		
			fig=pyplot.figure(figsize=(13,6))
			ax = fig.add_subplot(111)		
			
			width = 0.3	   # the width of the bars
			rects = []
			y = []
#			sys.exit(1)
			xlim = max(len(onduplicates[j]), len(offduplicates[j]))
			
			i=0
			# Check that there are on-target reads to avoid a division by zero. Then build bars for on-target reads.
			if(ontargetreads[j]>0):
				ind = []				
				overmaxreplicates = 0
				# Build y values at each number of replicates
				for i in range(int(xlim)):
					# It may happen that for example there are triplicates among on-target reads, but that there are not triplicates among off-target reads and
					# viceversa
					if(i>=len(onduplicates[j])):
						on = 0
						off = offduplicates[j][i]
					elif(i>=len(offduplicates[j])):
						on = onduplicates[j][i]
						off = 0
					else:
						on = onduplicates[j][i]
						off = offduplicates[j][i]

					if(i<MAXREPLICATES):					
						# ind will contain x-tick labels
						ind.append(i+1)
						y.append(on*100.0/ontargetreads[j])
					else:
						overmaxreplicates += on
						
				if(i>=MAXREPLICATES):
					# ind will contain x-tick labels
					ind.append(MAXREPLICATES+1)
					y.append(overmaxreplicates*100.0/ontargetreads[j])
					
				ind = numpy.array(ind)
				rects.append(ax.bar(ind, y, width, color='#46a246'))
			else:
				print 'WARNING: no ON-target reads were found at bam '+bam.filename
				print '	The most probable reason for this is an incorrect target (bed) file. Check that contig ids are exactly the same that appear in the bam'
				print '	file, e.g. check whether chromosome ids start with "chr" or not'

			# Check that there are off-target reads to avoid a division by zero. Then build bars for off-target reads.
			if(offtargetreads[j]>0):
				ind = []				
				y = []
				# Build y values at each number of replicates
				for i in range(int(xlim)):
					overmaxreplicates = 0
					# It may happen that for example there are triplicates among on-target reads, but that there are not triplicates among off-target reads and
					# viceversa
					if(i>=len(onduplicates[j])):
						on = 0
						off = offduplicates[j][i]
					elif(i>=len(offduplicates[j])):
						on = onduplicates[j][i]
						off = 0
					else:
						on = onduplicates[j][i]
						off = offduplicates[j][i]
					if(i<MAXREPLICATES):	
						ind.append(i+1)									
						y.append(off*100.0/offtargetreads[j])
					else:
						overmaxreplicates += off					
				if(i>=MAXREPLICATES):
					ind.append(MAXREPLICATES+1)
					y.append(overmaxreplicates*100.0/offtargetreads[j])
					
				ind = numpy.array(ind)														
				rects.append(ax.bar(ind+width, y, width, color='#ff0000'))
			
			ax.set_xticks(ind+width)	
			
			if(xlim>MAXREPLICATES):		
				ax.set_xticklabels(map(str, ind[:-1])+['>='+str(MAXREPLICATES+1)])
			else:
				ax.set_xticklabels(map(str, ind))
						
			ax.get_xaxis().tick_bottom()
			ax.get_yaxis().tick_left()
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.set_ylabel("% of reads")
			ax.set_xlabel("# times duplicated")
			
			# Shink current axis by 20%
			box = ax.get_position()
			ax.set_position([box.x0, box.y0, box.width*0.8, box.height])		
			ax.legend( tuple([rect[0] for rect in rects]), ('On target', 'Off target'), loc="upper left", bbox_to_anchor=(1,1) )
			
			# Include graph title to distinguish graphs from different samples
			if(len(bamlist)>1):				
				title = legend[j]		
				# Truncate the number of characters in the title to 70
				if(len(title)> 70):
					title = title[:70]+'...'
				ax.set_title(title)
		
			fig.savefig(outdir+'/duplicates'+str(j)+'.png')
			matplotlib.pyplot.close(fig)
			
					
			ws = wb.add_sheet(str(j+1)+'-'+legend[j][:10])			
			# Create header font
			ws.write(1,0,'# on'); ws.write(2,0,'# off'); ws.write(3,0,'# total'); ws.write(4,0,'% on'); ws.write(5,0,'% off')
						
			# Write the numbers used to build the bar plot in a spreadsheet
			for n in range(int(xlim)):
				overmaxreplicateson = 0
				overmaxreplicatesoff = 0

				# It may happen that for example there are triplicates among on-target reads, but that there are not triplicates among off-target reads and
				# viceversa			
				if(n>=len(onduplicates[j])):
					on = 0
					off = offduplicates[j][n]						
				elif(n>=len(offduplicates[j])):
					on = onduplicates[j][n]
					off = 0
				else:
					on = onduplicates[j][n]
					off = offduplicates[j][n]
				
				if(n<MAXREPLICATES):
					ws.write(0,n+1,n+1,header_style)
					total = on+off
					ws.write(1,n+1,on);ws.write(2,n+1,off);ws.write(3,n+1,total);
				
					# Check that there are on-target reads to avoid a division by zero
					if(ontargetreads[j]>0):
						ws.write(4,n+1,on*100/ontargetreads[j]);
					else:
						ws.write(4,n+1,0);
					
					# Check that there are off-target reads to avoid a division by zero	
					if(offtargetreads[j]>0):
						ws.write(5,n+1,off*100/offtargetreads[j]);
					else:
						ws.write(5,n+1,0)
				else:
					overmaxreplicateson += on
					overmaxreplicatesoff += off
				
			if(n>=MAXREPLICATES):
				ws.write(0,MAXREPLICATES+1,'>='+str(MAXREPLICATES+1),header_style)
				total = overmaxreplicateson+overmaxreplicatesoff
				ws.write(1,MAXREPLICATES+1,on);ws.write(2,MAXREPLICATES+1,off);ws.write(3,MAXREPLICATES+1,total);				
				# Check that there are on-target reads to avoid a division by zero
				if(ontargetreads[j]>0):
					ws.write(4,MAXREPLICATES+1,overmaxreplicateson*100/ontargetreads[j]);
				else:
					ws.write(4,MAXREPLICATES+1,0);
				
				# Check that there are off-target reads to avoid a division by zero	
				if(offtargetreads[j]>0):
					ws.write(5,MAXREPLICATES+1,overmaxreplicatesoff/offtargetreads[j]);
				else:
					ws.write(5,MAXREPLICATES+1,0)
					
				n = MAXREPLICATES
							
			ws.write(0,n+2,'Total reads', header_style)
			ws.write(1,n+2,ontargetreads[j])
			ws.write(2,n+2,offtargetreads[j])
		wb.save(outdir+'/duplicates.xls')

		
		
		
		
		
				
	def reads_on_target(self, bed, outdir, bamlist=[], legend=None, executiongranted=None, onoff_status=None, duplicates_status=None, 
						retonduplicates=None, retoffduplicates=None,enrichment=None, percontarget=None, tmpdir=None, warnthreshold=80):
		"""************************************************************************************************************************************************************
		Task: Print reads on traget and off target
		Inputs:
			bed: bed file with capture coordinates
			outdir: Output folder
			bamlist: list of bam_file objects representing bam files which are also wanted to be analyzed.
			legend: list of strings containing the label to tag each bam file in the xls and png files that will be created.
			executiongranted: multiprocessing.Semaphore object to control the use of machine resources.
			onoff_status: multiprocessing.Value object to return whether the number of on-target reads is extremely low (False) or not (True)
			duplicates_status: multiprocessing.Value object to return whether the number of duplicated reads on-target is greater than the number of duplicated
				off-target (False) or not (True).
			enrichment: multiprocessing.Array objecto to return the enrichment value for each bam (on-target reads per Kb)/(off target reads per Kb)
			percontarget: multiprocessing.Array object to return percentaje of reads on target for each bam file
			tmpdir: string containing the path to a temporary directory where temporary files will be stored.
		Outputs: two new files named reads_on_target.png and reads_on_target.xls will be generated at 'outdir'. In addition 'status' will be modified to indicate
			whether the number of on-target reads is extremely low (False) or not (True)		
		************************************************************************************************************************************************************"""
			
		global TMP
		
		# Check whether the method can be executed or not
		if(executiongranted<>None):
			executiongranted.acquire()
		
		# Check whether to use default tmpdir
		if(tmpdir<>None):
			TMP = tmpdir	
			
		# Create a list of bam_file objects which includes self. Generates an array with the total number of reads in each bam.
		bamlist = [self]+bamlist
		tread = []		
		for bam in bamlist:
			tread.append(bam.nreads())

		# Calculate number of reads and duplicated reads on/off target per chromosome 
		nread = []; onperchr = []; totalperchr = []; onduplicates = []; offduplicates = []		
		for bam in bamlist:
			nread_tmp,onperchr_tmp,totalperchr_tmp,onduplicates_tmp,offduplicates_tmp= bam.myReadsOnTarget(bed)
			nread.append(nread_tmp);onperchr.append(onperchr_tmp);totalperchr.append(totalperchr_tmp);onduplicates.append(onduplicates_tmp);
			offduplicates.append(offduplicates_tmp)

		bedobj = bed_file.bed_file(bed)
		targetsize = bedobj.size()
		
		# If enrichment was not provided as a parameter a new list of floats is created
		if(enrichment==None):
			enrichment = [0 for i in range(len(bamlist))]
			
		# Calculates enrichment for each bam		
		for i in range(len(bamlist)):
			# If all reads are on-target reads let enrichment=-1 to avoid division by zero 
			if(tread[i]==nread[i]):
				enrichment[i] = -1				
			else:
				enrichment[i] = (nread[i]*1000.0/targetsize)/((tread[i]-nread[i])*1000.0/(self.mappingsize()-targetsize)) # The same target must had been used for all bam files
		
		if(percontarget==None):
			percontarget = [None for i in range(len(bamlist))]
			
		self.generate_ontarget_results(bamlist, nread, tread, onperchr, totalperchr, enrichment, percontarget, outdir, legend)
		self.generate_duplicates_results(bamlist, nread, [tread[i]-nread[i] for i in range(len(nread))], onduplicates, offduplicates, outdir, legend)
		
						
		# Check whether semaphore should be released to allow the execution of some other process
		if(executiongranted<>None):
			executiongranted.release()
			
		# If the percentage of on-target reads is >= 80 results are assumed to be ok
		if(onoff_status<>None):			
			i=0
			while(i<len(nread) and (tread[i]>0) and ((nread[i]*100.0/tread[i])>=warnthreshold)):
				i+=1
			onoff_status.value = (i==len(nread))

		# If the percentage of duplicates on-target is greater than the percentage of duplicates off-target results are assumed to be ok
		if(duplicates_status<>None):
			i=0
			duplicates_status.value=True
			while(i<len(onduplicates) and duplicates_status.value):
				duplicates_status.value = (sum(onduplicates[i][1:])>sum(offduplicates[i][1:]))				
				i+=1				
				
		if(retonduplicates<>None):	
			for i in range(len(bamlist)):
				retonduplicates[i] = sum(onduplicates[i][1:])*100.0/tread[i]
				retoffduplicates[i] = sum(offduplicates[i][1:])*100.0/tread[i]
			
			

				
				
				
	def get_coverage_distribution(self, bedfilename):

		coveragefile = TMP+'/'+str(os.getpid())+'.coverage'

		print 'Calculating coverage per target position...'
		self.run('coverageBed -d -abam '+self.filename+' -b '+bedfilename+' > '+coveragefile)

		dist = []
		
		# A progress bar is initialized
		widgets = ['Loading coverage distribution: ', progressbar.Percentage(), ' ', 
					progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
		pbar = progressbar.ProgressBar(widgets=widgets, maxval=self.count_lines(coveragefile)).start() 
		
		i = 1
		fd = file(coveragefile)
		# Recorremos el archivo y guardamos la ultima columna en un vector
		for line in fd:
			dist.append(string.atoi(line.split("\t")[-1]))
			pbar.update(i)
			i+=1
			
		pbar.finish()
		fd.close()
		
		return dist
	
	
	

	def load_coverage_file(self, filename):
		"""*********************************************************************************************************************************************************
		Task: loads coverage counts from a .coverage file. Only regions of size > 3 are considered.
		Inputs:
			filename: string containing the name of the file to be loaded. Must contain the coverage
				per position per exon. 
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
		print 'Calculating file size...'
		pbar = progressbar.ProgressBar(widgets=widgets, maxval=self.count_lines(filename)).start() 
		
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




	
	def region_coverage_std(self, bedfilename):
		"""************************************************************************************************************************************************************
		Task: gets the sampling distribution of coverage standard deviation across regions.
		Inputs:
			bedfilename: string containing the name of the bed with the regions to analyze.
		Output: 
			std_sampling: list of real values containing the coverage std for each region in the bed file.
		************************************************************************************************************************************************************"""
		
		coveragefile = TMP+'/'+str(os.getpid())+'.coverage'

		print 'Calculating coverage per target position...'
		self.run('coverageBed -d -abam '+self.filename+' -b '+bedfilename+' > '+coveragefile)
		
		# Loads exons and coverage per position	   
		[loaded_exon, loaded_coverage, length] = self.load_coverage_file(coveragefile)
		
		# A progress bar is initialized
		widgets = ['Generating left position matrix: ', progressbar.Percentage(), ' ', 
					progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
		pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(loaded_exon)).start() 
		
		std_sampling = []
		positions_coverage = []		
		currexon = None
				
		# Goes through each region position
		for i in range(len(loaded_exon)):
			exon = loaded_exon[i]			   

			# If it is not, the first position of a new region has been found			 
			if(exon<>currexon):
				positions_coverage = numpy.array(positions_coverage)
				# Skip regions with 0 coverage in all positions
				if(len((positions_coverage>0).nonzero()[0])): std_sampling.append(positions_coverage.std()/positions_coverage.mean()) # Normalized standard deviation		 
				currexon = exon
				positions_coverage = [] 
					
			positions_coverage.append(loaded_coverage[i])
						  
			pbar.update(i+1)
		pbar.finish()
		
		return std_sampling
	
	
	
	
	
	def fastq(self):		
		pysam.bam2fq(self.filename)
		
		return fastq_file.fastq_file(self.filename.replace('.bam','.fq'))
	
	
	
	
	def minusbam(self, bam, fileout=None):
		
		pid = str(os.getpid())
		if(fileout==None): fileout=TMP+'/'+pid+'.minus.bam'
		
		tosubstract = bam.get_mapped_reads_ids()
		
		print str(len(tosubstract))+' reads found in right bam'
		
		bamout = bam_file(fileout, "wb", header=self.header)
		
		removed = 0
		written = 0
		for read in self:
			if(read.qname not in tosubstract):
				bamout.write(read)
				written += 1
			else:
				removed += 1
				
		print str(removed) + ' reads removed'
		print str(written) + ' reads remain at '+fileout
		bamout.close()
		pysam.index(fileout)
		
		return bam_file(fileout)





	def getMQ_regions(self,regions):
		"""*******************************************************************************************************************************************
		JPFLORIDO
		Task:  Get mapping quality for each region given the set of regions as input. For each region, the mean mapping quality for all reads related to the region is given
		Inputs: BAM file (self) and set of regions
		Outputs: 
			set of regions with the MQ value		
		*******************************************************************************************************************************************"""
	
		#totalRegions=len(regions.keys())
		#initRegion=1
		for currentRegion in sorted(regions.keys()):				
			# fetch works using 0-based indexing. However, when a region is give (init,end), the end position is not taken into account. For example if a read starts at position end, the read is not returned.
			# That is the reason the end value has a unit added
			allReads=self.fetch(str(currentRegion[0]),currentRegion[1],currentRegion[2]+1)
			allMQ=[]
			for currentRead in allReads: # Get mapping qualities for all reads
				allMQ.append(currentRead.mapq)
			if(len(allMQ)>0):		
				regions[currentRegion]=numpy.array(allMQ).mean()
				
		return regions



	def getBCQ_regions(self,regions,phredValue):
		"""*******************************************************************************************************************************************
		JPFLORIDO
		Task:  Get an score based on base calling quality for all reads and each region given the set of regions as input. For each region, the average base calling value is measured for each position taking into account all reads that map that position. Finally, an average value of all these measurements is done
		Inputs: BAM file (self) and set of regions and value for Phred offset (33 or 64)
		Outputs: 
			set of regions with the average Base Calling Quality value		
		*******************************************************************************************************************************************"""

		for currentRegion in sorted(regions.keys()):
			
			# Get all reads contained in the current region
			initRegion=currentRegion[1]
			endRegion=currentRegion[2]
			
			# fetch works using 0-based indexing. However, when a region is give (init,end), the end position is not taken into account. For example if a read starts at position end, the read is not returned.
			# That is the reason the end value has a unit added

			allReads=self.fetch(str(currentRegion[0]),initRegion,endRegion+1) # Fetch no devuelve la read que comienza justo en endRegion, asi que sumo 1. Las ventanas estan ya en base 1
			
			dicBCQ={} # Will contain all BCQ for each position in the region
			allBCQ=[] # Will contain average BCQ values per position in the region
			
			for currentRead in allReads:
				newQB=currentRead.qual # Vector of BCQ values for the current read
				if(len(currentRead.cigar)>1): # There are deletions and/or insertions. Build base calling vector according to these deletions/insertions										
					newQB=[]
					previousIndex=0
					for index in range(0,len(currentRead.cigar)): # Check all values of cigar informatiomn
						currentCode=currentRead.cigar[index]						
						if(currentCode[0]==0): # MATCH
							newQB.extend(currentRead.qual[previousIndex:(previousIndex+currentCode[1])])
							previousIndex+=currentCode[1]
						elif(currentCode[0]==1): # Insertion -> remove those bases involved in the insertion								
							previousIndex+=currentCode[1]							
						elif(currentCode[0]==2): # DEletion -> there is no BCQ for a deletion
							newQB.extend(['NUL']*currentCode[1])	
				
				currentReadInit=currentRead.pos # Base 0
				currentReadEnd=currentReadInit+len(newQB)-1 # Adapt the length of the read according to their insertions and or deletions (if any)
				
				if(currentReadInit<initRegion): # Read starts before the region -> avoid all values before the beginning of the region
					initPos=initRegion
					indexBCQ=initRegion-currentReadInit
				else:
					initPos=currentReadInit
					indexBCQ=0
				
				currentPos=initPos				
				while(currentPos<=currentReadEnd): # Move along all BCQ values from the beginning to the end
					if(newQB[indexBCQ]!='NUL'): # Due to the deletions
						if(currentPos in dicBCQ.keys()):
							dicBCQ[currentPos].append(ord(newQB[indexBCQ])-phredValue) # Base calling quality (based on Phred+33)
						else:
							dicBCQ[currentPos]=[ord(newQB[indexBCQ])-phredValue]
					indexBCQ+=1
					currentPos+=1
			
			# At this point, for each position between initRegion and endRegion, there is a list of BCQ	values -> get the average value per position
			for currentKey in dicBCQ.keys():
				allBCQ.append(numpy.array(dicBCQ[currentKey]).mean())
				
			if(len(allBCQ)>0):
				regions[currentRegion]=numpy.array(allBCQ).mean() # Finally, get the average value of BCQ for the whole region
		
	
			
		return regions





	def drawcoveragefeatures(self, coverage, sampling, fileout):
		"""************************************************************************************************************************************************************
		Task: generates the graph that compares coverage and features distribution.		 
		Inputs:
			coverage: list of tuples each containing a pair (position, coverage), both integers.
			sampling: list of integers containing the positions within current chromosome where features occur.
			fileout: string containing the full path to a png file where the graph will be saved.
		 Ouputs: a png file will be generated named fileout containing a graph that compares coverage and variant distribution. 
		************************************************************************************************************************************************************"""
		
		fig = pyplot.figure(figsize=(23,15))
		ax = fig.add_subplot(111)
		width = 0.5

		x = [values[0] for values in coverage]	
#		coverage_rect = ax.plot(x, numpy.log10(numpy.array([values[1] for values in coverage])), color='#ff683e', linewidth=5)
		axcoverage = ax.twinx()
		nvars_rect = ax.hist(sampling, bins=100, color='#004586')[2]
		coverage_rect = axcoverage.plot(x, numpy.array([values[1] for values in coverage]), color='#ff683e', linewidth=5)		
		
		# add some
#		ax.yaxis.set_ticks_position("right")
#		axcoverage.yaxis.set_ticks_position("left")		
		axcoverage.set_ylabel('coverage')
		ax.set_ylabel('# probes')
		ax.set_xlabel('Genomic position')
		ax.set_xlim(right=max(x))
#		axcoverage.set_ylim(top=1000)
	
		
		# Shink current axis by 20%
		box = ax.get_position()
	
		ax.legend( tuple([nvars_rect[0], coverage_rect[0]]), tuple(['# probes', 'Coverage']), loc="upper center", prop={'size':30} )
	
		# Changes font size of each text element in the graph
		for item in ([ax.xaxis.label, ax.yaxis.label, axcoverage.yaxis.label] +
				 ax.get_xticklabels() + ax.get_yticklabels() + axcoverage.get_yticklabels()):
			item.set_fontsize(30)
			item.set_weight('bold')
			
		fig.savefig(fileout)
		matplotlib.pyplot.close(fig)





	def coveragefeaturesdistr(self, chr, featuresbedfilename, targetbedfilename, dirout):
		"""************************************************************************************************************************************************************
		Task: generates a graph that compares coverage and features distribution.		 
		Inputs:
			featuresbedfilename: String containing the full path to a bed file		
			dirout: String containing the output directory where graphs will be saved (one per chromosome).
			targetbedfilename: String containing the full path to a bed file indicating target regions.
		 Ouputs: a set of png files will be generated at dirout containing graphs that compare coverage and features distribution. 
		************************************************************************************************************************************************************"""

		
		coverage = self.coveragetrackchr(chr,targetbedfilename,windowsize=480,overlap=160)
		sampling = bed_file.bed_file(featuresbedfilename).get_centers()['MT']
		self.drawcoveragefeatures(coverage, sampling, dirout+'/'+chr+'.png')
		
		
		
		
		
	def unannotated(self, gtffilename, fileout):
		pid = str(time.time())
		tmpsam = TMP+'/'+pid+'.sam'
		print """CMD: samtools view -h """+self.filename+""" | htseq-count --mode intersection-strict -t exon -i gene_id -o """+tmpsam+""" - """+gtffilename
		os.system("""samtools view -h """+self.filename+""" | htseq-count --mode intersection-strict -t exon -i gene_id -o """+tmpsam+""" - """+gtffilename)
		
		toremove = {}
		fd = file(tmpsam)
		for line in fd:
			if ('XF:Z:no_feature' not in line and 'XF:Z:ambiguous' not in line): 
				toremove[line.split('\t')[0]] = None				
		fd.close()
		
		newbam = bam_file(fileout, 'wb', header=self.header)
		
		for read in self:
			if(read.qname not in toremove):
				newbam.write(read)
			
		newbam.close()
		pysam.index(fileout)
		
		return bam_file(fileout)
		
		
		
		
	def mergealignments(self, bamfilename, out):

		alreadywritten = {}
				
		newbam = bam_file(out, 'wb', header=self.header)
		
		# A progress bar is initialized
		widgets = ['Initializing new bam: ', progressbar.Percentage(), ' ', 
					progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
		pbar = progressbar.ProgressBar(widgets=widgets, maxval=self.nreads()).start() 
		i = 0
				
		for read in self:
			alreadywritten[read.qname] = None
			newbam.write(read)
			i+=1
			pbar.update(i)
		pbar.finish()

		bam = bam_file(bamfilename, 'rb')
		
		chrids = [seqinfo['SN'] for seqinfo in self.header['SQ']]
		chridmap = [chrids.index(seqinfo['SN']) for seqinfo in bam.header['SQ']] 
						
		# A progress bar is initialized
		widgets = ['Merging alignments from '+bamfilename+': ', progressbar.Percentage(), ' ', 
					progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
		pbar = progressbar.ProgressBar(widgets=widgets, maxval=bam.nreads()).start() 
		i = 0
		mergedalignments = 0
		for read in bam:
			if(read.qname not in alreadywritten):
				read.rname = chridmap[read.rname]
				newbam.write(read)
				mergedalignments +=1				
		bam.close()
		pbar.finish()
		
		
		newbam.close()
		
		print str(mergedalignments)+' extra alignments found in '+bamfilename
		print str(self.nreads()+mergedalignments)+' total alignments included in '+out
		
		pysam.index(out)
		
		return bam_file(out)
		
		
		
		
		

		





#bam_file('/tmp/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.mapped.bam', 'rb').reads_on_target('/tmp/20130108.exome.targets.bed','/tmp',legend=['test'])
#bam_file('/tmp/test.bam').issorted()
#bam_file('/home/javi/MGP/course/data/GU_20130110_FC1_6L4S_JA_C148_F3.filtered.realigned.recalibrated.chr22.bam').reads_on_target('/home/javi/MGP/course/data/SeqCap_EZ_Exome_v3_primary.g1k.chr22.bed', '/tmp/course/')

#enrichment = multiprocessing.Array('f', 2)
#percontarget = multiprocessing.Array('f', 2)
#a = bam_file('/tmp/GU_20121002_FC2_6L4S_AL_C121_F3.filtered.singleHits.realigned.recalibrated.bam')
#b = bam_file('/tmp/shrimp-GU_20121002_FC2_6L4S_AL_C121_F3.filtered.realigned.recalibrated.bam')
#a.reads_on_target('/tmp/test1.bed', '/tmp/', bamlist=[b], enrichment=enrichment, percontarget=percontarget)
#a.reads_on_target('/tmp/test1.bed', '/tmp/', enrichment=enrichment, percontarget=percontarget)
#b = bam_file('/tmp/test2.sorted.bam')
#a.minusbam(b, '/tmp/minus_test.bam')
