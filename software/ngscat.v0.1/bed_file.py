import os
import numpy
import string
import region

try:
	import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot

except ImportError:
	print 'WARNING: module pyplot was not loaded.'
	
try:
	import progressbar
except ImportError:
	print 'WARNING: module progressbar was not loaded.'

import sets

#import bam_file

import pysam

import bam_file


CHR_LENGTHS = '/usr/local/reference_genomes/human/human_g1k_v37.1-22XYM.genome'
TMP = '/tmp/'
BEDTOOLS='/usr/local/bedtools/bin/'

class bed_file:
	
	
	
	
	def __init__(self, _filename):
		self.filename = _filename
		self.chrs = None
		self.nregions = None
		
		
		

	def checkformat(self):
		"""************************************************************************************************************************************************************
		Task: checks the format of the bed file. The only requirements checked are that each line presents at least 3 tab separated columns, the
			two on the right must present integer values indicating the start/end position respectively. Right value must be greater than the
			left value.
		Outputs:
			err: string containing the detected error. Empty string in case of a correct format.
		************************************************************************************************************************************************************"""
		
		fd = file(self.filename)

		line = fd.readline()
		fields = line.split('\t')
		lc = 1
		err = ''
		
		# Checks that the two columns on the right contain integer values
		try:
			# Parses each line and checks that there are at least 3 fields, the two on the right containing integer values and being the right one
			# greater than the left one
			while(line<>'' and len(fields)>2 and int(fields[1])<int(fields[2])):
				lc += 1
				line = fd.readline()
				fields = line.split('\t')				
		except ValueError:
			err += 'Incorrect start/end values at line '+str(lc)+'\n'
			err += 'Start/End coordinates must be indicated with integer values. The right value must be greater than the left value.\n'
			err += 'Line found: '+line
			fd.close()
		
			return err

		# If it get to this point means that either the file ended or there is a line with less than 3 fields
		if(line<>''):
			err += 'Incorrect line format at line '+str(lc)+'\n'
			err += 'At least three columns are expected in each line\n'
			err += 'The right value must be greater than the left value.\n'
			err += 'Line found: '+line
			fd.close()
		
		return err
		
		
		
		
		
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
	   
		
		
		
		
				
	def count_lines(self, filename):
#		print 'Calculating file size...'
#		tmp = os.popen('wc -l '+filename)
#		nlines = string.atof(tmp.readline().split(' ')[0])
#		
#		if(tmp.close()<>None):
#			print 'Error: some error occurred while running '
#			print '	wc -l '+filename
#			print 'at bam_file.py'
#			print 'Exiting'
#			sys.exit(1)		
#		print '	Done.'
		
		return len(file(self.filename).readlines())





	def __sub__(self, other):
		pid = str(os.getpid())
		newbed = bed_file(TMP+'/'+pid+os.path.basename(self.filename)[:-4]+'-'+os.path.basename(other.filename))
		self.run(BEDTOOLS+"bedtools subtract -a "+self.filename+" -b "+other.filename+" > "+newbed.filename)
		
		return newbed
	
	
	
	
	def intersect(self, other, fileout=None):
		pid = str(os.getpid())
		
		if(fileout==None):
			newbed = bed_file(TMP+'/'+pid+os.path.basename(self.filename)[:-4]+'_intersect_'+os.path.basename(other.filename))
		else:
			newbed = bed_file(fileout)
			
		self.run(BEDTOOLS+"bedtools intersect -a "+self.filename+" -b "+other.filename+" > "+newbed.filename)
		
		return newbed
		
		
		
		
		
	def intersectbam(self, bam):
		"""************************************************************************************************************************************************************
		Task: IntersectBam
		Inputs:
			bam: Bam_file type 
		************************************************************************************************************************************************************"""
		pid = str(os.getpid())
		newbam = TMP+'/'+pid+'.intersect.bam'
		self.run(BEDTOOLS+"intersectBed -abam "+bam.filename+" -b "+self.filename+" > "+newbam)	   
		pysam.index(newbam)
		
		return bam_file.bam_file(newbam,"rb")





	def sum_region_size(self):
		sum = 0
		fd = file(self.filename)
		for line in fd:
			parts = line.split('\t')
			sum += string.atoi(parts[2])-string.atoi(parts[1])
		fd.close()
		
		return sum





	def listids(self):
		fd = filename(self.filename)
		listids = []
		for line in fd:
			listids.append(line.split('\t')[3])
		fd.close()
			
		return listids




	def load_chr_lengths(self):
		lengths = {}
		
		print 'Loading chr lengths...'
		fd = file(CHR_LENGTHS)
		for line in fd:
			parts = line.split('\t') 
			lengths[parts[0]] = string.atoi(parts[1])
		print '	Done.'
		fd.close()		
		
		return lengths
		
		
		
	def extend(self, n, fileout=None):
		"""*******************************************************************************************************************************************
		Task: generates a new bed file in which regions of this bed are extended +-n bases.  
		Inputs:
			n: integer with the number of bases to extend.
			fileout: string containing the full path to the new bed file.
		Outputs: a new bed file will be created named fileout. In case fileout is not provided, a new file will be created named with the prefix of
			self.filename and ended in .extended.bed 
		*******************************************************************************************************************************************"""	
		
		lengths = self.load_chr_lengths()
		
		# If an output filename is not provided generates one
		if(fileout==None):
			fileout=self.filename.replace('.bed', '.extended'+str(n)+'.bed')
		
		# Each region in each line is extended +-n bases
		fd = file(self.filename)
		fdw = file(fileout, 'w')
		for line in fd:
			parts = line.split('\t')
			fdw.write(parts[0]+'\t'+str(max(0,string.atoi(parts[1])-n))+'\t'+str(min(lengths[parts[0]], string.atoi(parts[2])+n))+'\n')
		fd.close()
		fdw.close()
			
		return bed_file(self.filename.replace('.bed', '.extended'+str(n)+'.bed'))





	def extendnoref(self, n, fileout=None):
		"""*******************************************************************************************************************************************
		Task: generates a new bed file in which regions of this bed are extended +-n bases.  
		Inputs:
			n: integer with the number of bases to extend.
			fileout: string containing the full path to the new bed file.
		Outputs: a new bed file will be created named fileout. In case fileout is not provided, a new file will be created named with the prefix of
			self.filename and ended in .extended.bed 
		*******************************************************************************************************************************************"""	
				
		# If an output filename is not provided generates one
		if(fileout==None):
			fileout=self.filename.replace('.bed', '.extended'+str(n)+'.bed')
		
		# Each region in each line is extended +-n bases
		fd = file(self.filename)
		fdw = file(fileout, 'w')
		for line in fd:
			parts = line.split('\t')
			fdw.write(parts[0]+'\t'+str(max(0,string.atoi(parts[1])-n))+'\t'+str(string.atoi(parts[2])+n)+'\n')
		fd.close()
		fdw.close()
			
		return bed_file(fileout)
		
		
		
		
		
	def load_and_return(self):
		chrs = {}		
		nregions = self.count_lines(self.filename)
		
		# A progress bar is initialized
		widgets = ['Loading bed regions: ', progressbar.Percentage(), ' ', 
					progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
		pbar = progressbar.ProgressBar(widgets=widgets, maxval=nregions).start() 
		
		fd = file(self.filename)
		for i,line in enumerate(fd):
			parts = line.split('\t')
			if(parts[0] not in chrs):
				chrs[parts[0]] = [[], []]

			chrs[parts[0]][0].append(string.atoi(parts[1]))
			chrs[parts[0]][1].append(string.atoi(parts[2])-1)
			pbar.update(i+1)
		fd.close()				   
		pbar.finish()
		
		return chrs
		
		
		
		
		
	def load(self):
		"""************************************************************************************************************************************************************
		Task: loads data into self.chrs and self.nregions.
		Output: self.chrs and self.nregions are modified. 
			self.chrs: dictionary. Each key represents a chromosome. Values are lists of tuples (start,end) indicating each of the regions in the chromosome.
				>>>>>>> WARNING: ending coordinate is also transformed to base 0!!!! <<<<<<<
			self.nregions = self.count_lines(self.filename)
		************************************************************************************************************************************************************"""
		
		# self.chrs: dictionary. Each key represents a chromosome. Values are lists of tuples (start,end) indicating each of the regions in the chromosome.
		#	 >>>>>>> WARNING: ending coordinate is also transformed to base 0!!!! <<<<<<<
		# self.nregions: total number of regions in the bed file.
		self.chrs = {}		
		self.nregions = self.count_lines(self.filename)
		
		# A progress bar is initialized
#		widgets = ['Loading bed regions: ', progressbar.Percentage(), ' ', 
#					progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
#		pbar = progressbar.ProgressBar(widgets=widgets, maxval=self.nregions).start() 
		
#		print 'Loading bed regions...'
		
		# Each line is parsed to load the region it contains
		fd = file(self.filename)
		for i,line in enumerate(fd):
			parts = line.split('\t')
			# No coordinates for the chr in this line were already loaded
			if(parts[0] not in self.chrs):
				self.chrs[parts[0]] = []

			# Current coordinates are added to current chromosome. >>>> WARNING: ending coordinate is also transformed to base 0!!! <<<<<
			self.chrs[parts[0]].append((string.atoi(parts[1]),string.atoi(parts[2])-1))
#			pbar.update(i+1)
		fd.close()				   
#		pbar.finish()
		print '	Done.'
		




		
		
		
		
	def load_custom(self,base):
		"""************************************************************************************************************************************************************
		JPFLORIDO
		Task: loads data into self.chrs and self.nregions according to the base indicated as argument
		Inputs: BED file (self) and base (1 or 0). If base is not 1 or 0, bed file is loaded as it is....
		Output: self.chrs and self.nregions are modified. 
			self.chrs: dictionary. Each key represents a chromosome. Values are lists of tuples (start,end) indicating each of the regions in the chromosome.
			self.nregions = self.count_lines(self.filename)
			bed file loaded with coordinates in base "base"
		 		
		************************************************************************************************************************************************************"""
		
		self.chrs = {}		
		self.nregions = self.count_lines(self.filename)
		
		# A progress bar is initialized
#		print 'Loading bed regions...'
#		widgets = ['Loading bed regions: ', progressbar.Percentage(), ' ', 
#					progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
#		pbar = progressbar.ProgressBar(widgets=widgets, maxval=self.nregions).start() 
		
		if(base==0): # Base 0 -> end position -1
			initOffset=0
			endOffset=-1
		elif(base==1): # Base 1-> start position +1
			initOffset=1
			endOffset=0
		else: # Any other value...
			initOffset=0 # As it is...
			endOffset=0
				
		# Each line is parsed to load the region it contains
		fd = file(self.filename)
		for i,line in enumerate(fd):
			parts = line.split('\t')
			# No coordinates for the chr in this line were already loaded
			if(parts[0] not in self.chrs):
				self.chrs[parts[0]] = []

			# Current coordinates are added to current chromosome. >>>> WARNING: check coordinates accoridng to "base" argument <<<<<<<<<<<<
			self.chrs[parts[0]].append((string.atoi(parts[1])+initOffset,string.atoi(parts[2])+endOffset))
#			pbar.update(i+1)
#		fd.close()				   
#		pbar.finish()




	
	def covered_bases(self):
		"""************************************************************************************************************************************************************
		Task: calculates the number of bases covered by regions in this bed appropriately merging overlaps. It is required that self.chrs is already loaded by 
			calling self.load()
		Output:
			nbases: integer containing the number of bases covered by regions in this bed. Overlapping regions are merged to count each "overlapped" base just once.
		************************************************************************************************************************************************************"""

		nbases = 0
					   
		# Passes through the regions in each chromosome
		for chr in self.chrs:
			# Sorts regions in current chromosome by first coordinate (and second coordinate in case of same first coordinate)
#			print 'Chr '+chr
#			print '	Sorting regions...'
			self.chrs[chr].sort()
			
			# A progress bar is initialized
#			print 'Ignoring overlaps...'

#			widgets = ['Ignoring overlaps: ', progressbar.Percentage(), ' ', 
#						progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
#			pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(self.chrs[chr])).start() 
			
			# Passes through each tuple (start,end) and processes its overlaps
			curr=0
			while(curr<len(self.chrs[chr])):
				
				# "i" is the index of the next region to compare with curr
				# start: integer containing the starting position of the next set of overlapped regions
				i=curr+1
				start = self.chrs[chr][curr][0]
				entered = 0
				
				if(chr=='20'):
					a=1				
				# Check whether starting position of next region falls within current region
				while(i<len(self.chrs[chr]) and self.chrs[chr][curr][0]<=self.chrs[chr][i][0] and self.chrs[chr][i][0]<=self.chrs[chr][curr][1]):
					if(self.chrs[chr][curr][0]<>self.chrs[chr][i][0] or self.chrs[chr][curr][1]<>self.chrs[chr][i][1]):
						a=1
					# curr region is always the one with the higher ending position. If the ending point of "i" region is higher than the one of "curr" region, update
					# curr and make "i" be the next region
					if(self.chrs[chr][curr][1]<self.chrs[chr][i][1]):						
						curr = i
						i = curr+1
					# otherwise  process the next overlapping region and compare it with curr
					else:
						i+=1
				
					entered += 1
				
				if(entered>1):
					a=1
				# The number of bases is calculated as the ending position of the region with the highest ending coordinate within this set of overlapping regions, minus
				# the starting coordinate of this whole set of overlapped regions. 
				# curr is updated to "i", which is the first region that does not overlap with curr according to the "while" conditions
				nbases += (self.chrs[chr][curr][1]-start+1)
				curr = i  
#				pbar.update(curr)
		
#		pbar.finish()
		
		return nbases		





		 



	def non_overlapping_exons(self,baseCodification,outputFile=None,tmpdir=None):
		"""************************************************************************************************************************************************************
		JPFLORIDO
		Task: Get exons of a given bed file removing overlapping areas
		Inputs: 
		    self: bed file
		    baseCodification: whether exons are in "real" base 0, 1 or as it is....
		Output:
			A set of tuples for each exon: chromosome, exon begin position, exon end position
		Requirements: 	WARNING: BED FILE MUST BE SORTED BEFORE
		************************************************************************************************************************************************************"""
	   		
	   	global TMP
	   	
	   	if(tmpdir<>None):
	   		TMP=tmpdir
	   		
		chromosomes=[]
		start_positions=[]
		end_positions=[]
				
		if(self.chrs==None): 
			if(baseCodification==1):
				self.load_custom(1) # "Real" Base 1
			elif(baseCodification==0):
				self.load_custom(0) # "Real" Base 0
			else:
				self.load_custom(-1) # Base 0 as it is described for standard BED format
					   
		# Passes through the regions in each chromosome
		for chr in self.chrs:
			# Sorts regions in current chromosome by first coordinate (and second coordinate in case of same first coordinate)
#			print 'Chr '+chr
#			print '	Sorting regions...'
			self.chrs[chr].sort()
			
			# A progress bar is initialized
#			print 'Ignoring overlaps...'

#			widgets = ['Ignoring overlaps: ', progressbar.Percentage(), ' ', 
#						progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
#			pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(self.chrs[chr])).start() 
			
			# Passes through each tuple (start,end) and processes its overlaps
			curr=0
			while(curr<len(self.chrs[chr])):
				
				# "i" is the index of the next region to compare with curr
				# start: integer containing the starting position of the next set of overlapped regions
				i=curr+1
				start = self.chrs[chr][curr][0]
				entered = 0
				
				if(chr=='20'):
					a=1				
				# Check whether starting position of next region falls within current region
				while(i<len(self.chrs[chr]) and self.chrs[chr][curr][0]<=self.chrs[chr][i][0] and self.chrs[chr][i][0]<=self.chrs[chr][curr][1]):
					if(self.chrs[chr][curr][0]<>self.chrs[chr][i][0] or self.chrs[chr][curr][1]<>self.chrs[chr][i][1]):
						a=1
					# curr region is always the one with the higher ending position. If the ending point of "i" region is higher than the one of "curr" region, update
					# curr and make "i" be the next region
					if(self.chrs[chr][curr][1]<self.chrs[chr][i][1]):						
						curr = i
						i = curr+1
					# otherwise  process the next overlapping region and compare it with curr
					else:
						i+=1
				
					entered += 1
				
				if(entered>1):
					a=1
				# The number of bases is calculated as the ending position of the region with the highest ending coordinate within this set of overlapping regions, minus
				# the starting coordinate of this whole set of overlapped regions. 
				# curr is updated to "i", which is the first region that does not overlap with curr according to the "while" conditions				
				chromosomes.append(chr)
				start_positions.append(start)
				end_positions.append(self.chrs[chr][curr][1]) # The ending coordinate depends on the base codification (baseCodification argument)
				curr = i  
#				pbar.update(curr)
		
#		pbar.finish()

		# Write new BED file 
		pid = str(os.getpid())
		
		if(outputFile==None):
			outputFile = TMP+'/'+pid+os.path.basename(self.filename.replace('.bed', '_noOverlapping.bed'))
			
		fdw=file(outputFile,'w')
		for i,currentChromosome in enumerate(chromosomes):
			fdw.write(str(chromosomes[i])+'\t'+str(start_positions[i])+'\t'+str(end_positions[i])+'\n')
		fdw.close()
			
		return bed_file(outputFile)


	def sort_bed(self,newbedfilename=None):
		"""************************************************************************************************************************************************************
		JPFLORIDO
		Task: Sort a BED file by chromosome and start position
		Input: BED file (self)
		Output:
			An object with the sorted BED file
		************************************************************************************************************************************************************"""
	   	
		if(newbedfilename==None):
			pid = str(os.getpid())			
			newbedfilename = TMP+'/'+pid+os.path.basename(self.filename.replace('.bed', '_sorted.bed'))
			
		self.run(""" sortBed -i """+self.filename+""" > """+newbedfilename)
		return bed_file(newbedfilename)	 





	def my_sort_bed(self, newbedfilename=None, tmpdir=None):
		"""************************************************************************************************************************************************************
		JPFLORIDO
		Task: Sort a BED file by chromosome and start position. This function does not use bedtools
		Input: BED file (self)
		Output:
			An object with the sorted BED file
		************************************************************************************************************************************************************"""
		
		global TMP
		
		if(tmpdir<>None):
			TMP = tmpdir
			
		if(newbedfilename==None):
			pid = str(os.getpid())			
			newbedfilename = TMP+'/'+pid+os.path.basename(self.filename.replace('.bed', '_sorted.bed'))
		
		# Load bed file (chromosome, start and end information are stored -> standard bed format)
		self.load_custom(-1) # as it is....
		# Sort it
		# Passes through the regions in each chromosome
		fdw=file(newbedfilename,'w')			
		for chr in self.chrs:
			# Sorts regions in current chromosome by first coordinate (and second coordinate in case of same first coordinate)
#			print 'Chr '+chr
#			print '	Sorting regions...'
			self.chrs[chr].sort()
			
			# Store it
			curr=0
			while(curr<len(self.chrs[chr])): 
				fdw.write(str(chr)+'\t'+str(self.chrs[chr][curr][0])+'\t'+str(self.chrs[chr][curr][1]) +'\n')
				curr+=1
				
		fdw.close()
					
		# Return new bed file
		return bed_file(newbedfilename)	 





	def size(self):
		"""************************************************************************************************************************************************************
		Task: calculates the number of bases covered by regions in this bed appropriately merging overlaps.
		Output:
			nbases: integer containing the number of bases covered by regions in this bed. Overlapping regions are merged to count each "overlapped" base just once.
		************************************************************************************************************************************************************"""
		
		# Check whether self.chrs is already loaded
		if(self.chrs==None): self.load()
		nbases = 0
		
		
		
#		print 'Calculating bed size...'
		return self.covered_bases()

			
			
			
			
	def meansize(self):
		total = 0
		nregions = 0
		fd = file(self.filename)
		for line in fd:
			fields = line.split('\t')
			total+= float(fields[2])-float(fields[1])
			nregions += 1
		fd.close()
		
		return total/nregions
   




	def mediansize(self):
		distr = []
		fd = file(self.filename)
		for line in fd:
			fields = line.split('\t')
			distr.append(float(fields[2])-float(fields[1]))
		fd.close()
		
		return numpy.median(distr)





	def sizepercentile(self, n):
		distr = []
		fd = file(self.filename)
		for line in fd:
			fields = line.split('\t')
			distr.append(float(fields[2])-float(fields[1]))
		fd.close()
		
		return len((numpy.array(distr)<n).nonzero()[0])*100.0/len(distr)




	def percentilesize(self, p):
		distr = []
		fd = file(self.filename)
		for line in fd:
			fields = line.split('\t')
			distr.append(float(fields[2])-float(fields[1]))
		fd.close()
		
		return numpy.percentile(distr, p)
		
		



	def sizemax(self):
		distr = []
		fd = file(self.filename)
		for line in fd:
			fields = line.split('\t')
			distr.append(float(fields[2])-float(fields[1]))
		fd.close()
		
		return max(distr)





	def histogram(self, fileout):
		"""*********************************************************************************************************************************************************
		Task: 
		Inputs:
		*********************************************************************************************************************************************************"""
				
		distr = []
		fd = file(self.filename)
		for line in fd:
			fields = line.split('\t')
			distr.append(float(fields[2])-float(fields[1]))
		fd.close()
			
		fig = pyplot.figure()
		ax = fig.add_subplot(111)
		n, bins, patches = ax.hist(distr, bins=1000, facecolor='green')
		ax.set_xlabel('Region length')
		ax.set_ylabel('Frequency')
#		ax.set_xlim(lengths.min(), lengths.max())
		ax.set_xlim(0, 1000)
 #	   ax.set_ylim(0, 2000)
#		ax.set_title(histogram_title)
		ax.grid(True)
		
		fig.savefig(fileout)
		matplotlib.pyplot.close(fig)
		

			
			
			
			
	def load_base1(self):
		"""************************************************************************************************************************************************************
		JPFLORIDO
		Task: loads data into self.chrs and self.nregions.
		Output: self.chrs and self.nregions are modified. 
			self.chrs: dictionary. Each key represents a chromosome. Values are lists of tuples (start,end) indicating each of the regions in the chromosome.
			self.nregions = self.count_lines(self.filename)
		COORDINATES ARE IN BASE 1
		************************************************************************************************************************************************************"""
		
		# self.chrs: dictionary. Each key represents a chromosome. Values are lists of tuples (start,end) indicating each of the regions in the chromosome.
		#	 >>>>>>> WARNING: ending coordinate is in base 1!!!! <<<<<<<
		# self.nregions: total number of regions in the bed file.
		self.chrs = {}		
		self.nregions = self.count_lines(self.filename)
		
		# A progress bar is initialized
#		print 'Loading bed regions...'

#		widgets = ['Loading bed regions: ', progressbar.Percentage(), ' ', 
#					progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
#		pbar = progressbar.ProgressBar(widgets=widgets, maxval=self.nregions).start() 
		
		# Each line is parsed to load the region it contains
		fd = file(self.filename)
		for i,line in enumerate(fd):
			parts = line.split('\t')
			# No coordinates for the chr in this line were already loaded
			if(parts[0] not in self.chrs):
				self.chrs[parts[0]] = []

			# Current coordinates are added to current chromosome. >>>> WARNING: ending coordinate is in base 1 !!!!!! <<<<<<<<<<<<
			self.chrs[parts[0]].append((string.atoi(parts[1]),string.atoi(parts[2])))
#			pbar.update(i+1)
		fd.close()				   
#		pbar.finish()
		
		
		
		
		
	def getWindows(self,windowSize,offset,chrs=None):
		"""************************************************************************************************************************************************************
		JPFLORIDO
		Task: get the set of coordinates (chromosome, start, end) of each window for the current BED file
		Inputs:
			windowSize: size of the window
			offset: offset to get a new window
			baseCodification: base to be used ("real" base 1, "real" base 0 or base 0 as it is described for standard BED format)
		Output:
			nbases: integer containing the number of bases covered by regions in this bed. Overlapping regions are merged to count each "overlapped" base just once.
		************************************************************************************************************************************************************"""						
				
		windowSet={}
		
		if(chrs==None):
			chrs = self.chrs
			
		for chr in chrs: # For each chromosome, get regions invovled in that chromosome
			
			print 'Current chromosome->'+ chr
			for region in self.chrs[chr]:							
				initRegion=region[0]
				endRegion=region[1]				
				if (((endRegion-initRegion)+1)<=windowSize):
					windowSet[(chr,initRegion,endRegion)]=-1 # Current region is shorter than window size					
				else: # The region is greater than the window size
					addWindow=True					
					initWindow=initRegion
					endWindow=initRegion+windowSize-1																		
					while(addWindow):		
						windowSet[(chr,initWindow,endWindow)]=-1
																				
						if(endWindow<endRegion): # More windows can be added							
							initWindow=initWindow+offset
							endWindow=endWindow+offset
						else: # window ends goes beyond the end of the region (modify the last element)
							del windowSet[(chr,initWindow,endWindow)]
							windowSet[(chr,initWindow,endRegion)]=-1												
							addWindow=False							
	
		return windowSet
		
		
		
		
		
	def sortcountkeys(self, counts):
		"""*******************************************************************************************************************************************
		Task: obtains a sorted list with the keys in counts. 
		Inputs:
			counts: dictionary. Keys are tuples of the form ('chr', start,end), where chr is a string and start/end are integers. The elements
				of the dictionary are numbers.
		Outputs: sorted list containing the keys of 'counts'.
		*******************************************************************************************************************************************"""	
		
		# Transform chromosome identifiers to integers to allow a correct sorting
		keys = map(list,counts.keys())
		for coordinate in keys:
			# To account for 'X' and 'Y' chromosomes
			try:
				coordinate[0] = string.atoi(coordinate[0])
			except ValueError:
				continue;
		keys.sort()
		
		# Transform chromsome identifiers to strings
		for coordinate in keys:
			coordinate[0] = str(coordinate[0])			
			
		return map(tuple, keys)





	def windows2bedgraph(self,windows,fileout,startshift=0,endshift=0):
		fd = file(fileout, 'w')
		fd.write('browser pack refGene encodeRegions\n')
		fd.write('browser full altGraph\n')
		fd.write('track type=bedGraph name="BedGraph Format" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n')
		
		sortedkeys = self.sortcountkeys(windows)
		for region in sortedkeys:
			fd.write(region[0]+'\t'+str(region[1]+startshift)+'\t'+str(region[2]+endshift)+'\t'+str(windows[region])+'\n')
		fd.close()





	def windows2circosformat(self,windows,fileout,startshift=0,endshift=0):
		fd = file(fileout, 'w')
		
		sortedkeys = self.sortcountkeys(windows)
		for region in sortedkeys:
			fd.write('hs'+region[0]+' '+str(region[1]+startshift)+' '+str(region[2]+endshift)+' '+str(windows[region])+'\n')
		fd.close()





	def get_centers(self):
		centers = {}
		fd = file(self.filename)
		for line in fd:
			fields = line.split('\t')
			if(fields[0] not in centers):
				centers[fields[0]] = [((int(fields[1])+1)+(int(fields[2])))/2]
			else:
				centers[fields[0]].append(((int(fields[1])+1)+(int(fields[2])))/2)
		
		fd.close()
		
		return centers
	
	
	
	
	
	def get_region(self):
		"""************************************************************************************************************************************************************
		Task: return a list of regions with all intervals
		Output:
			regions: list of regions
		************************************************************************************************************************************************************"""
		
		fd=file(self.filename)
		regions=[]
		for line in fd:
			aline=line.replace('\n','').split('\t')
			#new region 
			r=region.region(aline[0],aline[1],aline[2])
			regions.append(r)
		return regions
		
		fd.close()

#bed_file('/tmp/test.bed').checkformat()
