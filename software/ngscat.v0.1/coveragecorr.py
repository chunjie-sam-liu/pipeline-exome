
import sys
import numpy

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

from pylab import plot, figure, imshow, xlabel, ylabel, cm, show, polyfit, polyval
from scipy import stats, mgrid, c_, reshape, random, rot90


def region_coverage(coveragefile):
    """************************************************************************************************************************************************************
    Task: calculates mean coverage per region.
    Inputs:
        coveragefile: string with the full path to a text file where coverage per target position is specified. Format:
            20    68318    68439    93
            20    68318    68439    95
            20    68318    68439    96
            20    68318    68439    96
            20    68318    68439    96
            20    68318    68439    94
            ...            
    Outputs: 
        coverage: dictionary. Keys are region coordinates in the form of tuples. Elements are floats containing mean coverage per region.
    ************************************************************************************************************************************************************"""   
    
    print 'Calculating mean coverage per region...'

    # Each line contains the coverage of one base    
    i = 1
    coverage = {}
    curr_region = None
    fd = file(coveragefile)
    for line in fd:
        parts = line.split('\t')
        newregion = (parts[0], int(parts[1]), int(parts[2]))

        # Check if the region of current base is the same as the one of the base before
        if(curr_region<>newregion):
            # In case a new region is found, save mean coverage for previous region
            if(curr_region<>None):
                coverage[curr_region] = numpy.mean(sampling)

            sampling = []
            curr_region = newregion 

        sampling.append(int(parts[-1]))
        i+=1

    coverage[curr_region] = numpy.mean(sampling)
                   
    fd.close()
    
    return coverage





def coveragecorr(coveragefiles,fileout,legend,executiongranted=None,status=None,corr=None,warnthreshold=0.9):
    """************************************************************************************************************************************************************
    Task: calculates coverage correlation between two samples.
    Inputs:
        coveragefiles: list of two strings, each containing the full path to a .coverage file.            
        fileout: String indicating the full path to the png file that will be created.
        legend: list of two strings, each indicating a tag to label the axis of the correlation plot. (First element in legend must correspond to
            first file in coveragefiles, second element in legend must correspond to the second element in coveragefiles).
        executiongranted: multiprocessing.Semaphore object to control the use of machine resources.
        status: multiprocessing.Value object to return whether both samples present correlated coverages or not.            
    Outputs: a new png file will be created displaying coverage correlation between both samples.
    ************************************************************************************************************************************************************"""
    
    # If executiongranted is not set proceed with the execution
    if(executiongranted<>None):
        executiongranted.acquire()

    coverage1 = region_coverage(coveragefiles[0]) # Calculate mean coverage per region
    coverage2 = region_coverage(coveragefiles[1]) # Calculate mean coverage per region
#    coverage1 = [1,2]

    # More than one target region is necessary for correlation calculation
    if(len(coverage1)>1):
        x = []
        y = []
        # Load the values that will be plotted in x and y
        for region in coverage1:
            x.append(coverage1[region])
            y.append(coverage2[region])
        x = numpy.array(x)
        y = numpy.array(y)
        
#        x = numpy.array(list(numpy.random.uniform(20,30,1000))+list(numpy.random.uniform(120,130,1000)))
#        y = numpy.array(list(numpy.random.uniform(20,30,1000))+list(numpy.random.uniform(120,130,1000)))

#        x = numpy.array(list(numpy.random.normal(30,5,3000))+list(numpy.random.normal(130,5,3000)))
#        y = numpy.array(list(numpy.random.normal(30,5,3000))+list(numpy.random.normal(130,5,3000)))
        
#        x = n
        xmin = x.min()
        xmax = x.max() # Due to the imshow sentence, we need to rescale gccontent from [0,1] to [0,100]
        ymin = y.min()
        ymax = y.max()
         
        # Perform a kernel density estimator on the results
        X, Y = mgrid[xmin:xmax:100j, ymin:ymax:100j]
        positions = c_[X.ravel(), Y.ravel()]
        values = c_[x, y]
        kernel = stats.kde.gaussian_kde(values.T)
        Z = reshape(kernel(positions.T).T, X.T.shape)
        
        fig = pyplot.figure(figsize=(6,6))
        ax = fig.add_subplot(111, axisbg="#000000")      
    
        m,b,r_value,p_value,std_err=stats.linregress(x,y)
        ax.plot(range(0,int(xmax)), m*numpy.arange(0,int(xmax))+b, linewidth=0.6) 
        
        # Fills in the 'd' vector which represents the "function" that determines the colour of each point in the graph
        xfactor = 100/(xmax-xmin)
        yfactor = 100/(ymax-ymin)
        d = []
        for i in range(len(x)):                
            # Z dimensions are (ymax-ymin+1)x(xmax-xmin+1)        
            xvalue = int(round((x[i]-xmin)*xfactor))-1
            # Rounding may yield negative values after substracting -1. Truncate to 0.
            if(xvalue<0): xvalue=0
            yvalue = int(round((y[i]-ymin)*yfactor))-1
            # Rounding may yield negative values after substracting -1. Truncate to 0.
            if(yvalue<0): yvalue=0              
            
            d.append(Z[xvalue][yvalue])
            
        sc=ax.scatter(x,y,cmap=cm.Reds,c=d,s=1,edgecolors='none')
        ax.text(3,ymax+13,'R=%.3f'%r_value, fontsize=12)
        cbar=fig.colorbar(sc,ticks=[numpy.min(Z)+0.00001,numpy.max(Z)])
        cbar.ax.set_yticklabels(['Low','High'])
        cbar.set_label('Density')
        
        
        xlabel = 'Coverage - '+legend[0]
        ylabel = 'Coverage - '+legend[1]
        
        # Truncate x label to 50 characters
        if(len(xlabel)>50):
            xlabel = xlabel[:50]+'...'
        # Truncate y label to 50 characters            
        if(len(ylabel)>50):
            ylabel = ylabel[:50]+'...'
                
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(left=0,right=xmax+10)
        ax.set_ylim(bottom=0,top=ymax+10)
        fig.savefig(fileout)
        matplotlib.pyplot.close(fig)
        
        if(corr<>None):
            corr.value = r_value
            
        # In case status variable is set, return the status. True if correlation >= 0.9
        if(status<>None):
            status.value = r_value>=warnthreshold
                
    else:
        print 'WARNING: only one region found in the bed file. Skipping coverage correlation calculation.'
        
        # In case status variable is set, return False as correlation could not be calculated
        if(status<>None):
            status.value = False

    # Release semaphore            
    if(executiongranted<>None):
        executiongranted.release()
        
def coveragecorrtest():
    
    x = numpy.array(list(numpy.random.normal(30,5,1000))+list(numpy.random.normal(130,5,1000)))
    y = numpy.array(list(numpy.random.normal(30,5,1000))+list(numpy.random.normal(130,5,1000)))
    
    xmin = x.min()
    xmax = x.max() # Due to the imshow sentence, we need to rescale gccontent from [0,1] to [0,100]
    ymin = y.min()
    ymax = y.max()
     
    # Perform a kernel density estimator on the results
    X, Y = mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = c_[X.ravel(), Y.ravel()]
    values = c_[x, y]
    kernel = stats.kde.gaussian_kde(values.T)
    Z = reshape(kernel(positions.T).T, X.T.shape)
    
    
    # Fills in the 'd' vector which represents the "function" that determines the colour of each point in the graph
    xfactor = 100/(xmax-xmin)
    yfactor = 100/(ymax-ymin)
    d = []
    for i in range(len(x)):                
        # Z dimensions are (ymax-ymin+1)x(xmax-xmin+1)        
        xvalue = int(round((x[i]-xmin)*xfactor))-1
        # Rounding may yield negative values after substracting -1. Truncate to 0.
        if(xvalue<0): xvalue=0
        yvalue = int(round((y[i]-ymin)*yfactor))-1
        # Rounding may yield negative values after substracting -1. Truncate to 0.
        if(yvalue<0): yvalue=0              
        
        d.append(Z[xvalue][yvalue])
    
    fig = pyplot.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    sc=ax.scatter(x,y,cmap=cm.Reds,c=d,s=1,edgecolors='none')
#    sc=ax.imshow(rot90(Z),cmap=cm.Reds,extent=[xmin, xmax, ymin, ymax]) # Due to the imshow sentence, we need to rescale gccontent from [0,1] to [0,100]    
#    sc=ax.imshow(rot90(Z),cmap=cm.gist_earth_r,extent=[xmin, xmax, ymin, ymax]) # Due to the imshow sentence, we need to rescale gccontent from [0,1] to [0,100]
    cbar=fig.colorbar(sc, ticks=(numpy.min(Z)+0.00001,numpy.max(Z)))
  
#    cbar.set_ticks((numpy.min(Z),numpy.max(Z)))
#    cbar=fig.colorbar(sc,ticks=None)
#    cbar=fig.colorbar(sc,ax=ax, ticks=[0,1])
#    cbar.set_ticks(())
    cbar.set_ticklabels(['Low','High'])
#    cbar.update_ticks()

    cbar.set_label('Density')
    ax.set_xlabel('GC content (%)')
    ax.set_ylabel('Mean coverage')
    
    fig.savefig('/tmp/tmp.png')
    matplotlib.pyplot.close(fig)
	
    

    
#coveragecorr(None, '/tmp/tmp.png', ['a','b'])
#coveragecorrtest()
