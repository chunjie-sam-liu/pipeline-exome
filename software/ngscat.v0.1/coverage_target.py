import sys
import string
import numpy
import matplotlib
matplotlib.use('Agg')
import pylab

sys.path.append('/home/javi/MGP/utils/')
import bed_file


def print_coverage(coverageFiles, npoints, outdir, legend=None, executiongranted=None, status=None, lowcovbases=None, warnregionsize=100, warncoveragethreshold=6):
    """*******************************************************************************************************************************************
        Task: Plot coverage per position
        Inputs:
            coverageFiles: List of strings with the paths to standard bedtools coverageBed files
            npoints: Number points in all chromosome graphs            
            outdir: Out directory
            legend: list of strings, each containing the label to tag each sample in the chromosome graphs
            executiongranted: multiprocessing.Semaphore object to control the use of machine resources.
            status: multiprocessing.Value object to return whether there are more than three points with coverage <=6x (False) or not (True)
        Outputs: it will generate one png file named outdir+chr+'_Ontarget_Coverage.png' for each chromosome. In addition, a txt file will be
            created listing target intervals with 0 coverage.
            status: may be modified to indicate whether there are more than three points with coverage <=6x (False) or not (True). Will only be
                modified if was passed as a parameter (<>None)            
    *******************************************************************************************************************************************"""
        
    # Check whether access is granted
    if(executiongranted<>None):
        executiongranted.acquire()
    
    colours = ['#46a246', '#ff0000', '#00ff00', '#0000ff', '#cc0011', '#007722', '#110066', '#c1c1c1', '#544db1', '#aa5198', '#bbd1e9', '#f1c4ab', '#24687a']
    allx = []
    ally = []
    allzeros = []
    # Calculate graphs for each coveragefile
    for fileidx,coverageFile in enumerate(coverageFiles):
        maxconsecutivelowcovbases = 0
        fdon=file(coverageFile)
        pos=0
        #Calculate size of intervals
        for line in fdon:
            pos=pos+1
        #fix size of intervals
        size=round(pos/npoints*1.0)
        fdon=file(coverageFile)
        
        # initialize variables
        i=0
        x=[]
        y=[]
        xaux=[]
        yaux=[]
        chroms=[]
        chrom=''
        sum=0
        a_point=0
        zeroflag=False
        zeroini=0
        zeroend=0
        chrzero=''
        zeros=[]
        real_pos=0
        ini_inter=0
        warningcounter = 0
        warning = False
        #######################
        
        # Reads and processes each line in the coverage file
        for line in fdon:
            row=line.split('\t')
            
            if int(row[1])==ini_inter:
                real_pos=real_pos+1
            else:
                real_pos=int(row[1])
                ini_inter=int(row[1])
                if zeroflag==True:
                    aux=(chrzero,zeroini,zeroend)
                    zeros.append(aux)
                    zeroflag=False 
                    
            
            
            #split if is a new chromosome
            if chrom!=row[0] and chrom!='':
                chroms.append(chrom)
                chrom=row[0]
                if i!=0:
                    yaux.append(sum/(i*1.0))
                    xaux.append(a_point)
                    x.append(xaux)
                    y.append(yaux)
                a_point=0
                xaux=[]
                yaux=[]
                sum=0
                i=0
                warningcounter = 0
                if zeroflag==True:
                    aux=(chrzero,zeroini,zeroend)
                    zeros.append(aux)
                    zeroflag=False   
            ########################################
            chrom=row[0]    
            #Calculate coverage per interval
            if i<size:
                sum=sum+int(row[3])
                i=i+1
            else:
                yaux.append(sum/(i*1.0))
                xaux.append(a_point)
                i=0
                sum=0
                a_point=a_point+1
            ################################################
            
            #save zeros
            if int(row[3])==0 and zeroflag==False:
                chrzero=row[0]
                zeroini=real_pos
                zeroend=real_pos
                zeroflag=True
            elif int(row[3])==0 and zeroflag==True:
                zeroend=real_pos
            elif int(row[3])!=0 and zeroflag==True:
                aux=(chrzero,zeroini,zeroend)
                zeros.append(aux)
                zeroflag=False   
            ###########
            
#            if (not warning):
            if(int(row[3])<warncoveragethreshold):
                warningcounter += 1
            else:
                if(warningcounter>0): 
                    maxconsecutivelowcovbases = max(maxconsecutivelowcovbases,warningcounter)                    
                warningcounter = 0
            
            if(warningcounter>warnregionsize):
                warning = True
                
            
            
        #Save the last chromosome
        if i!=0:
            yaux.append(sum/(i*1.0))
            xaux.append(a_point)
        chroms.append(chrom)     
        x.append(xaux)
        y.append(yaux)
        ###########################
        
        allx.append(x)
        ally.append(y)
        allzeros.append(zeros)   
            
        if(lowcovbases<>None):
            lowcovbases[fileidx] = maxconsecutivelowcovbases
            
        
    
    pylab.rc('axes', linewidth=3.0)
    #Drawing points
    i=0
    for chr in chroms:
        fig = matplotlib.pyplot.figure(figsize=(25,6))
        ax = fig.add_subplot(111)
        rects = []
        for j in range(len(coverageFiles)):
            rects.append(ax.plot(allx[j][i],ally[j][i],color=colours[j%12],linewidth=3.0))
        ax.axhline(linewidth=3.0)
        ax.axvline(linewidth=3.0)
    
    
        #Axe style###############################
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)   
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.set_xticklabels('')
        ax.set_xlabel(chr)
        ax.set_ylim(bottom=0)    
        ##########################################
        
        for item in ([ax.xaxis.label] + ax.get_yticklabels()):
            item.set_fontsize(30)
            item.set_weight('bold')
        
        if(legend<>None and len(legend)>1):
            
            # Shink current axis by 20%
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width, box.height*0.9])
        
#            ax.legend( tuple([rect[0] for rect in rects]), tuple(legend), loc="upper left", bbox_to_anchor=(1,1) )
            ax.legend( tuple([rect[0] for rect in rects]), tuple(legend), loc="lower left", bbox_to_anchor=(0,1.03) )
#            ax.legend( tuple([rect[0] for rect in rects]), tuple(legend), loc="lower left" )
        
        fig.savefig(outdir+chr+'_Ontarget_Coverage.png')
        matplotlib.pyplot.close(fig)
        i=i+1

    fdw=file(outdir+'NoCoverage.txt', 'w')        
    for i in range(len(coverageFiles)):
        if(legend<>None):
            fdw.write('############################\n') 
            fdw.write('#'+legend[0])
            fdw.write('\n############################\n')
        for region in allzeros[i]:
            fdw.write(region[0]+'\t'+str(region[1])+'\t'+str(region[2])+'\n')
    fdw.close()
            
    if(executiongranted<>None):
        executiongranted.release()
        
    if(status<>None):
        status.value = (not warning)


#testeo    
#ontarget=sys.argv[1]
#out=sys.argv[2]
#print_coverage(ontarget,1000, out)
    
