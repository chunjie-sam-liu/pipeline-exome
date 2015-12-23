
try:
    import numpy
except ImportError:
    print 'WARNING: numpy module was not imported.'
try:
    from scipy import stats
except ImportError:
    print 'WARNING: scipy.stats module was not imported.'

try:
    import progressbar
except ImportError:
    print 'WARNING: progressbar was not imported'
    
import os
import sys
import bed_file
import region

from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot


class bedgraph_file:
    
    def __init__(self, _filename):
        self.filename = _filename
        
        
        
        
    
    def scores(self):
        sampling = []
        
        print 'Loading scores from '+self.filename+'...'
        fd = file(self.filename)
        fd.readline();fd.readline();fd.readline()
        for i,line in enumerate(fd):
            value = (line.split('\t')[-1])
            try:
                sampling.append(float(value))
            except ValueError:
                if(value<>'none\n' and value<>'None\n'):
                    print 'ERROR: invalid score value at line '+str(i+4)
                    sys.exit(1)
        fd.close()
        print '    Done.'

        return numpy.array(sampling)
    
    
    
    
    def windowsscores(self): #ESTA FUNCION NO LA HE COMPROBADO
        data = []
        print 'Loading data from '+self.filename+' ...'
        fd = file(self.filename)
        fd.readline();fd.readline();fd.readline()
        for i,line in enumerate(fd):
            fields = (line.split('\t'))
            try:
                data.append(fields[0],fields[1],fields[2],float(fields[-1]))
            except ValueError:
                if(fields[-1]<>'none\n' and fields[-1]<>'None\n'):
                    print 'ERROR at bedgraph_file.windowsscores: invalid score value at line '+str(i+4)
                    sys.exit(1)
                    
        return data
    
    
    
        
    def spearman(self, bedgraph):
        sampling1 = self.scores()
        sampling2 = bedgraph.scores()
        
        return stats.spearmanr(sampling1, sampling2)
        
        
        
    def pearson(self, bedgraph):
        sampling1 = self.scores()
        sampling2 = bedgraph.scores()
        
        return stats.pearsonr(sampling1, sampling2)
        
        
        
        
    def outliers(self, fileout, prob=0.01):
        """************************************************************************************************************************************************************
        Task: selects those regions with a score <= percentile '100*prob' or >= percentile '100-prob*100'.   
        Inputs:            
            fileout: string containing the full path to the bedgraph file were outlier regions will be saved.
            prob: float indicating the percentile (/100) to use as the threshold that determines a region to be an outlier.
        Ouputs: a new bedgraph file will be created named fileout.
        ************************************************************************************************************************************************************"""
        
        sampling = self.scores()
        inferior = numpy.percentile(sampling, int(prob*100))
        superior = numpy.percentile(sampling, 100-int(prob*100))
        
        print 'Lower bond: '+str(inferior)
        print 'Higher bond: '+str(superior)
        
        print 'Selecting outliers...'
        writtenlow = 0
        writtenhigh = 0       
        fdw = file(fileout, 'w')
        fdr = file(self.filename)
        fdw.write(fdr.readline());fdw.write(fdr.readline());fdw.write(fdr.readline())
        for i,line in enumerate(fdr):
            # Scores may be 'none'
            try:
                score = float(line.split('\t')[-1])
                #Check whether the score is lower/higher than the percentiles used as thresholds
                if(score<=inferior):
                    fdw.write(line)
                    writtenlow += 1
                elif(score>=superior):
                    fdw.write(line)
                    writtenhigh += 1
            except ValueError:
                continue;
            
        fdr.close()
        fdw.close()
        print '    Done.'
        
        print str(writtenlow)+' outliers below the lower bond'
        print str(writtenhigh)+ ' outliers above the upper bond'
        
        return bedgraph_file(fileout)





    def loutliers(self, fileout, prob=0.01):
        """************************************************************************************************************************************************************
        Task: selects those regions with a score <= percentile '100*prob'   
        Inputs:            
            fileout: string containing the full path to the bedgraph file were outlier regions will be saved.
            prob: float indicating the percentile (/100) to use as the threshold that determines a region to be an outlier.
        Ouputs: a new bedgraph file will be created named fileout.
        ************************************************************************************************************************************************************"""
        
        sampling = self.scores()
        inferior = numpy.percentile(sampling, int(prob*100))
        
        print 'Lower bond: '+str(inferior)
        
        print 'Selecting outliers...'
        writtenlow = 0      
        fdw = file(fileout, 'w')
        fdr = file(self.filename)
        fdw.write(fdr.readline());fdw.write(fdr.readline());fdw.write(fdr.readline())
        for i,line in enumerate(fdr):
            # Scores may be 'none'            
            try:
                score = float(line.split('\t')[-1])
                #Check whether the score is lower than the percentile used as the threshold                
                if(score<=inferior):
                    fdw.write(line)
                    writtenlow += 1
            except ValueError:
                continue;
            
        fdr.close()
        fdw.close()
        print '    Done.'
        
        print str(writtenlow)+' outliers below the lower bond'
        
        return bedgraph_file(fileout)
    
    
    

    def houtliers(self, fileout, prob=0.01):
        """************************************************************************************************************************************************************
        Task: selects those regions with a score >= percentile '100-prob*100'.   
        Inputs:            
            fileout: string containing the full path to the bedgraph file were outlier regions will be saved.
            prob: float indicating the percentile (/100) to use as the threshold that determines a region to be an outlier.
        Ouputs: a new bedgraph file will be created named fileout.
        ************************************************************************************************************************************************************"""
        
        sampling = self.scores()
        superior = numpy.percentile(sampling, 100-int(prob*100))
        
        print 'Higher bond: '+str(superior)
        
        print 'Selecting outliers...'
        writtenhigh = 0       
        fdw = file(fileout, 'w')
        fdr = file(self.filename)
        fdw.write(fdr.readline());fdw.write(fdr.readline());fdw.write(fdr.readline())
        for i,line in enumerate(fdr):
            # Scores may be 'none'                        
            try:
                score = float(line.split('\t')[-1])
                #Check whether the score is higher than the percentile used as the threshold                                
                if(score>=superior):
                    fdw.write(line)
                    writtenhigh += 1
            except ValueError:
                continue;
            
        fdr.close()
        fdw.close()
        print '    Done.'
        
        print str(writtenhigh)+ ' outliers above the upper bond'
        
        return bedgraph_file(fileout)




    
    def draw_another_bed_score(self, bedgraph, fileout):
        """************************************************************************************************************************************************************
        Task: draws the distribution of the scores of the regions in 'bedgraph' that overlap with the regions in self. 
        Inputs:
            bedgraph: bedgraph_file object representing the other bedgraph file.
            fileout: string containing the full path to the png file were the distribution will be saved.
        Ouputs: a new png file will be created named fileout.
        ************************************************************************************************************************************************************"""
        
        print 'Loading coordinates from '+self.filename
        starts1 = {}
        ends1 = {}
        fd = file(self.filename)
        fd.readline(); fd.readline(); fd.readline()
        for line in fd:
            fields = line.split('\t')
            # Check whether there is already an entry for current chromosome in starts1            
            if(fields[0] not in starts1):
                starts1[fields[0]] = [float(fields[1])]
                ends1[fields[0]] = [float(fields[2])]
            else:
                starts1[fields[0]].append(float(fields[1]))
                ends1[fields[0]].append(float(fields[2]))
        fd.close()
        print '    Done.'

        print 'Loading coordinates and scores from '+bedgraph.filename
        starts2 = {}
        ends2 = {}
        scores2 = {}        
        fd = file(bedgraph.filename)
        fd.readline(); fd.readline(); fd.readline()
        for line in fd:
            fields = line.split('\t')
            # Check whether there is already an entry for current chromosome in starts2            
            if(fields[0] not in starts2):
                starts2[fields[0]] = [float(fields[1])]
                ends2[fields[0]] = [float(fields[2])]
                scores2[fields[0]] = [float(fields[3])]
            else:
                starts2[fields[0]].append(float(fields[1]))
                ends2[fields[0]].append(float(fields[2]))
                scores2[fields[0]].append(float(fields[3]))
        fd.close()
        print '    Done.'
        
        # Transform the lists of each chromosome in a numpy.array        
        for chr in starts1:
            starts1[chr] = numpy.array(starts1[chr])
            ends1[chr] = numpy.array(ends1[chr])
            
        # Transform the lists of each chromosome in a numpy.array            
        for chr in starts2:            
            starts2[chr] = numpy.array(starts2[chr])
            ends2[chr] = numpy.array(ends2[chr])

        widgets = ['Comparing windows: ', progressbar.Percentage(), ' ', 
                   progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
        pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(starts1)).start() 
        
        exactcount = 0
        x=[]; y=[]
        # Compare regions in each chromosome independently        
        for k,chr in enumerate(starts1):
            # Check that there are regions at current chromosome in 'bedgraph'            
            if(chr in starts2):
                # For each region of current chromosome in self, look for overlaps in bedgraph                
                for i in range(len(starts1[chr])):
                    # Check whether any region in bedgraph overlaps with the 'start' limit of current region                                                                
                    overlap = ((starts1[chr][i]>=starts2[chr]) * (starts1[chr][i]<=ends2[chr])).nonzero()[0]
                    if(len(overlap)>0):
                        j = 0
                        # Check whether any of the overlapping regions is actually the same exact region                        
                        while(j<len(overlap) and (starts1[chr][i]<>starts2[chr][overlap[j]] or ends1[chr][i]<>ends2[chr][overlap[j]])):
                            j += 1
                        
                        # Check whether an exact match was found in the while above                        
                        if(j<len(overlap)):
                            exactcount += 1
                            x.append(scores2[chr][overlap[j]])
                        else:
                            x.append(scores2[chr][overlap[0]])                     
                    else:
                        # Check whether any region in bedgraph overlaps with the 'end' limit of current region                                                
                        overlap = ((ends1[chr][i]>=starts2[chr]) * (ends1[chr][i]<=ends2[chr])).nonzero()[0]
                        if(len(overlap)>0):
                            x.append(scores2[chr][overlap[0]])
            pbar.update(k+1)
        print '    Done.'
                
        pbar.finish()
        print str(len(x))+' windows overlap between both tracks'
        print str(exactcount)+' windows perfectly match between both tracks'
                               
        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        n, bins, patches = ax.hist(x, 50, facecolor='green')
        ax.set_xlabel(os.path.basename(bedgraph.filename))
        ax.set_ylabel('Frequency')
        fig.savefig(fileout)
        matplotlib.pyplot.close(fig)




    
    def compare_graph(self, bedgraph, fileout):
        """************************************************************************************************************************************************************
        Task: compares the scores of this bedgraph file with the scores in another bedgraph file. 
        Inputs:
            bedgraph: bedgraph_file object representing the other bedgraph file.
            fileout: string containing the full path to the png file were the comparative graph will be saved.
        Ouputs: a new png file will be created named fileout.
        ************************************************************************************************************************************************************"""
        
        starts1 = {}
        ends1 = {}
        scores1 = {}
        fd = file(self.filename)
        fd.readline(); fd.readline(); fd.readline()
        # Loads chr,start,end,score from each line
        for i,line in enumerate(fd):
            fields = line.split('\t')
            # Check whether there is already an entry for current chromosome in starts1
            if(fields[0] not in starts1):
                try:
                    scores1[fields[0]] = [float(fields[3])]
                    starts1[fields[0]] = [float(fields[1])]
                    ends1[fields[0]] = [float(fields[2])]                    
                except ValueError:
                    print 'WARNING: '+fields[3]+' at line '+str(i+4)+' of file '+self.filename+' is not a number.'
                    
            else:
                try:
                    scores1[fields[0]].append(float(fields[3]))
                    starts1[fields[0]].append(float(fields[1]))
                    ends1[fields[0]].append(float(fields[2]))                    
                except ValueError:
                    print 'WARNING: '+fields[3]+' at line '+str(i+4)+' of file '+self.filename+' is not a number.'


        fd.close()

        starts2 = {}
        ends2 = {}
        scores2 = {}        
        fd = file(bedgraph.filename)
        fd.readline(); fd.readline(); fd.readline()
        # Loads chr,start,end,score from each line
        for i,line in enumerate(fd):
            fields = line.split('\t')
            # Check whether there is already an entry for current chromosome in starts2
            if(fields[0] not in starts2):
                starts2[fields[0]] = [float(fields[1])]
                ends2[fields[0]] = [float(fields[2])]
                
                try:
                    scores2[fields[0]] = [float(fields[3])]
                except ValueError:
                    print 'WARNING: '+fields[3]+' at line '+str(i+4)+' of file '+bedgraph.filename+' is not a number.'
                
            else:
                starts2[fields[0]].append(float(fields[1]))
                ends2[fields[0]].append(float(fields[2]))
                
                try:
                    scores2[fields[0]].append(float(fields[3]))
                except ValueError:
                    print 'WARNING: '+fields[3]+' at line '+str(i+4)+' of file '+bedgraph.filename+' is not a number.'
                
        fd.close()

        # Transform the lists of each chromosome in a numpy.array
        for chr in starts1:
            starts1[chr] = numpy.array(starts1[chr])
            ends1[chr] = numpy.array(ends1[chr])
        
        # Transform the lists of each chromosome in a numpy.array    
        for chr in starts2:            
            starts2[chr] = numpy.array(starts2[chr])
            ends2[chr] = numpy.array(ends2[chr])

        widgets = ['Comparing windows: ', progressbar.Percentage(), ' ', 
                   progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
        pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(starts1)).start() 

        exactcount = 0
        x=[]; y=[]
        # Compare regions in each chromosome independently
        for k,chr in enumerate(starts1):
            # Check that there are regions at current chromosome in 'bedgraph'
            if(chr in starts2):
                # For each region of current chromosome in self, look for overlaps in bedgraph
                for i in range(len(starts1[chr])):
                    # Check whether any region in bedgraph overlaps with the 'start' limit of current region                                
                    overlap = ((starts1[chr][i]>=starts2[chr]) * (starts1[chr][i]<=ends2[chr])).nonzero()[0]
                    if(len(overlap)>0):
                        j = 0
                        # Check whether any of the overlapping regions is actually the same exact region
                        while(j<len(overlap) and (starts1[chr][i]<>starts2[chr][overlap[j]] or ends1[chr][i]<>ends2[chr][overlap[j]])):
                            j += 1
                        
                        # Check whether an exact match was found in the while above
                        if(j<len(overlap)):
                            exactcount += 1
                            x.append(scores1[chr][i])
                            y.append(scores2[chr][overlap[j]])
                        else:
                            x.append(scores1[chr][i])
                            y.append(scores2[chr][overlap[0]])                     
                    else:
                        # Check whether any region in bedgraph overlaps with the 'end' limit of current region                        
                        overlap = ((ends1[chr][i]>=starts2[chr]) * (ends1[chr][i]<=ends2[chr])).nonzero()[0]
                        if(len(overlap)>0):
                            x.append(scores1[chr][i])
                            y.append(scores2[chr][overlap[0]])
                pbar.update(k+1)
                
        pbar.finish()                            
                
        print str(len(x))+' windows overlap between both tracks'
        print str(exactcount)+' windows perfectly match between both tracks'
                               
        fig = pyplot.figure(figsize=(13,10))
        ax = fig.add_subplot(111)

        ax.scatter(x,y)
        ax.set_ylabel(os.path.basename(bedgraph.filename))
        ax.set_xlabel(os.path.basename(self.filename))
        
#        fig = pyplot.figure()
#        ax = fig.add_subplot(111, projection='3d')
#        hist, xedges, yedges = numpy.histogram2d(x, y, bins=4)
#        
#        elements = (len(xedges) - 1) * (len(yedges) - 1)
#        xpos, ypos = numpy.meshgrid(xedges+0.25, yedges+0.25)
#        
#        xpos = xpos.flatten()
#        ypos = ypos.flatten()
#        zpos = numpy.zeros(elements)
#        dx = 0.5 * numpy.ones_like(zpos)
#        dy = dx.copy()
#        dz = hist.flatten()
#        
#        ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')
#        
#        ax.set_ylabel(os.path.basename(bedgraph.filename))
#        ax.set_xlabel(os.path.basename(self.filename))
        
#        ax.bar3d(dx, dy, dz, color='b', zsort='average')
        
        
        fig.savefig(fileout)
        matplotlib.pyplot.close(fig)
        
        
    
    
    
    def rmscore(self, score, out):
        """************************************************************************************************************************************************************
        Task: removes all of those regions that present score 'score'. 
        Inputs:       
            out: string containing the full path to the bedgraph file were unfiltered regions will be saved.
            score: float indicating the score of the regions that should be removed.
        Ouputs: a new bedgraph file will be created containing the unfiltered regions.
        ************************************************************************************************************************************************************"""
        
        fdr = file(self.filename)
        fdw = file(out, 'w')
        
        removed = 0
        remain = 0
        print 'Filtering...'
        fdw.write(fdr.readline());fdw.write(fdr.readline());fdw.write(fdr.readline())
        for line in fdr:
            fields = line.split('\t')
            if(float(fields[-1])<>score):
                fdw.write(line)
                remain += 1
            else:
                removed += 1
        print '    Done.'
        
        fdw.close()
        fdr.close()
    
        print str(removed) + ' regions removed'
        print str(remain) + ' regions remain'
        print str(removed+remain) + ' regions in total'
        
        
        
        
        
    def get_region(self):
        """************************************************************************************************************************************************************
        Task: returns a list of regions with all intervals
        Output:
            regions: list of regions
        ************************************************************************************************************************************************************"""
        
        fd=file(self.filename)
        regions=[]
        for line in fd:
            aline=line.split('\t')
            #new region 
            r=region.region(aline[0],aline[1],aline[2],aline[3])
            regions.append(r)
        return regions
    
    
    def get_batch(self,fd,size):
        """************************************************************************************************************************************************************
        Task: returns a list of n lines 
        Input:
            fd:file handler
            size:~size of batch in bytes
        Output:
            batch: list of lines
            fd:file handler
        ************************************************************************************************************************************************************"""
        batch=[]
        batch=fd.readlines(size)
        return batch,fd
        
        
    
    def getOffTarget(self,offset,coverageThreshold,target,outfile,tmpdir=None):
        """************************************************************************************************************************************************************
        Task: selects off-tareget(+offset) regions with a coverage >  coverageThreshold
        Inputs:       
            offset: integer indicating the number of bases to extend the target.
            coverageThreshold: integer indicating the coverage threshold to select the region
            target: ROIs bed file
        Ouputs: a new bedgraph file will be created containing selected regions.
        ************************************************************************************************************************************************************"""
              
        pid = str(os.getpid())
        tmpbed = tmpdir+'/'+pid+'.extended.bed'
        
        bed=bed_file.bed_file(target)
        extendedBed=bed.extendnoref(offset,tmpbed)
        sortedBed=extendedBed.my_sort_bed()
        nonOverlappingBed=sortedBed.non_overlapping_exons(-1) # Base 0, it is a standard BED
        finalBed=nonOverlappingBed.my_sort_bed() # BED file in base 0
        finalBed.load_custom(-1) # Load chromosome and positions in base 0                 
        bed_region=finalBed.get_region()
        bed_index=0 #index to control bed_region position
        
        
        fd=file(self.filename)
        header=fd.readline()
        reading=True #boolean to control while loop
        chr_found=False
        batch_n=1
        fdw=file(outfile,'w')
        
        while reading:
            batch,fd=self.get_batch(fd, 10000000)
#            print batch_n
            batch_n=batch_n+1
            
            if batch==[]:
                reading=False
            else:
                for line in batch:
                    aline=line.replace('\n','').split(' ')
                    #new region 
                    r=region.region(aline[0],aline[1],aline[2],aline[3])
                    search_open=True
 
                    while search_open:
                        type_overlap=r.overlap_type(bed_region[bed_index])

                        
                        if type_overlap==0: #bed region comes before bedgraph region
                            search_open=True
                            
                            if bed_index+1<len(bed_region) and (chr_found==False or (chr_found==True and r.chrom==bed_region[bed_index].chrom)):
                                bed_index=bed_index+1
                            elif r.value>=coverageThreshold:
                                search_open=False 
                                for region_selected in r-bed_region[bed_index]:
                                    fdw.write(str(region_selected))
                            else:
                                search_open=False 
                                                               
                        
                        elif type_overlap==-1: #bed region comes after bedgraph region
                            search_open=False
                            chr_found=True
                            if r.value>=coverageThreshold:
                                for region_selected in r-bed_region[bed_index]:
                                    fdw.write(str(region_selected))
                                
                        else:
                            search_open=False
                            chr_found=True
                            if r.value>=coverageThreshold:
                                for region_selected in r-bed_region[bed_index]:
                                    fdw.write(str(region_selected))
        fd.close()
                
        
        
