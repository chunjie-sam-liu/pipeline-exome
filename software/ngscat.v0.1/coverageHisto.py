'''
Created on 06/02/2013

@author: antonior
'''
import numpy
import matplotlib
matplotlib.use('Agg')
import pylab
import xlwt

def histo (alldata,xlab,ylab,outdir,legend=None):

    colours = ['#46a246', '#ff0000', '#00ff00', '#0000ff', '#cc0011', '#007722', '#110066', '#c1c1c1', '#544db1', '#aa5198', '#bbd1e9', '#f1c4ab', '#24687a']
        
    fig=matplotlib.pyplot.figure(figsize=(13,6))
    ax = fig.add_subplot(111)
             
    fig2=matplotlib.pyplot.figure()
    ax2 = fig2.add_subplot(111)

    maximum=0
    for data in alldata:
        maximum = max(maximum, max(data))
        
#    print 'MAXIMUM = '+str(maximum)
    bin = range(1,maximum+1,max(1,maximum/50))  
    histograms = []
    for colouridx,data in enumerate(alldata):   
        histograms.append(ax.hist(data, bins=bin, facecolor=colours[colouridx%12], alpha=0.5)[2])

    #Quitamos los ejes de la derecha y arriba
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    #########################
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)        
    ax.set_xlim(left=1)
    ax.set_xbound(lower=1)
    ax.set_xticks([1]+range(0,maximum+1,maximum/10)[1:])
           
    bp = ax2.boxplot(alldata, sym='+', vert=1, whis=1.5)

    for i in range(len(alldata)):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = zip(boxX,boxY)
        boxPolygon = matplotlib.patches.Polygon(boxCoords, facecolor=colours[i%12])
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []    
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            matplotlib.pyplot.plot(medianX, medianY, 'black')

        # Finally, overplot the sample averages, with horizontal alignment
        # in the center of each box
        ax2.plot([numpy.average(med.get_xdata())], [numpy.average(alldata[i])], color='w', marker='*', markeredgecolor='k')
        ax2.add_patch(boxPolygon)
    
    matplotlib.pyplot.setp(bp['boxes'], color='black')
    matplotlib.pyplot.setp(bp['whiskers'], color='black')
    matplotlib.pyplot.setp(bp['fliers'], color='red', marker='+')
    
    ax2.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',alpha=0.3)
    ax2.set_axisbelow(True)
    #Quitamos los ejes de la derecha y arriba
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.get_xaxis().tick_bottom()
    ax2.get_yaxis().tick_left()
    #########################
    ax2.set_xticklabels('')    
    
    
    if(legend<>None and len(legend)>1):

        if(len(legend[0])>25):
            ax2.set_xticklabels([tag[:25]+'...' for tag in legend])
        else:
            ax2.set_xticklabels(legend)
                
        # Shink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width, box.height*0.9])
    
#            ax.legend( tuple([rect[0] for rect in rects]), tuple(legend), loc="upper left", bbox_to_anchor=(1,1) )
        ax.legend( tuple([histogram[0] for histogram in histograms]), tuple(legend), loc="lower left", bbox_to_anchor=(0,1.03) )
    
    
    fig.savefig(outdir+'Coverage_histo')
    matplotlib.pyplot.close(fig)
    fig2.savefig(outdir+'Coverage_boxp')
    matplotlib.pyplot.close(fig2)





def histo_CV(infiles, outdir, legend=None, executiongranted=None, status=None, meancoverage=None, warnthreshold=40):
    
    if(executiongranted<>None):
        executiongranted.acquire()

    alldata=[]
    for infile in infiles:
        alldata.append([])
        fd=file(infile)
        # Recorremos el archivo y guardamos la ultima columna en un vector
        for line in fd:
            alldata[-1].append(int(line.split("\t")[-1][:-1]))
        fd.close()

#    print 'ALLDATA = '+str(alldata)
    
    #Llamos a la funcion histo para plotear la distribucion (histograma y boxplot)
#    histo(dist,"#46a246",50,"Coverage","Count",outdir)
    histo(alldata,"Coverage","Count",outdir,legend)

    #Imprimimos las estadisticas
    # Initialize the workbook and sheet
    wb = xlwt.Workbook()
    ws = wb.add_sheet('Coverage summary')
    
    # Create header font
    header_style = xlwt.easyxf('font: bold on')
    
    ws.write(0,0,'File',header_style);ws.write(0,1,'# bases coverage 0',header_style);ws.write(0,2,'Q1',header_style);ws.write(0,3,'Q2',header_style);ws.write(0,4,'Q3',header_style);
    ws.write(0,5,'Maximum',header_style);
    ws.write(0,6,'Minimum',header_style);ws.write(0,7,'Mean',header_style)
    
    meanok = True    
    for i,dist in enumerate(alldata):       
        #Sacamos estadisticas
        ndist=numpy.array(dist)
        p25=numpy.percentile(ndist, 25)
        p50=numpy.percentile(ndist, 50)
        p75=numpy.percentile(ndist, 75)
    #    p95=numpy.percentile(ndist, 95)
        maximum=numpy.max(ndist)
        minimum=numpy.min(ndist)
        mean=numpy.average(ndist)                        
        
        if(legend<>None):
            ws.write(i+1,0,legend[i]);
        
        ws.write(i+1,1,len(ndist)-len(ndist.nonzero()));ws.write(i+1,2,p25);ws.write(i+1,3,p50);ws.write(i+1,4,p75);ws.write(i+1,5,maximum);ws.write(i+1,6,minimum);
        ws.write(i+1,7,mean);
        
        meanok = (meanok and mean>=warnthreshold)
        
        if(meancoverage<>None):
            meancoverage[i] = mean

    wb.save(outdir+'/percentile.xls')
    
    if(executiongranted<>None):
        executiongranted.release()
        
    if(status<>None):
        status.value = meanok
        
#histo_CV(['/tmp/test2.coverage', '/tmp/test1.coverage'], '/tmp/')
