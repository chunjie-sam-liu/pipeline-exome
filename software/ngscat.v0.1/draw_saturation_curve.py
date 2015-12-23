#!/usr/bin/python

import optparse
import os
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import string
import xlwt


def draw_saturation_curve(filelist, ylabel, fileout, legend=None, warnthreshold=1e-5):
	"""************************************************************************************************************************************************************
	Task: draws a coverage saturation plot with the data returned by simulated_depth.py
	Inputs:
		filelist: String containing comma separated names of the files that contain the percentaje of covered positions at each depth.
		fileout: String containing the name of the file where the plot will be saved.			  
	Output: a .png image (fileout) containing the saturation plot.					 
	************************************************************************************************************************************************************"""
	
	# samplelist is used to keep the order in which samples appear at filelist, so that we can relate samples with their corresponding entry at the legend list
	samplelist = []
	
	y = {}	
	x = {}
	# Load data stored in each file
	for afile in filelist:
		fd = file(afile)
#		sample = afile.split('_')[-1]
		sample = fd.readline()[:-1]
		parts = fd.readline().split('\t')
		# Check whether there is already an entry for current sample. Load depth and % of covered positions.
		if(sample not in y):
			y[sample] = [string.atof(parts[-1])]
			x[sample] = [string.atof(parts[0])]
			
			# A list of samples is maintained in the order they appear in filelist, which coincides with the order provided at the legend list 
			samplelist.append(sample)
		else:
			y[sample].append(string.atof(parts[-1]))
			x[sample].append(string.atof(parts[0]))
		fd.close()
		   
	fig = pyplot.figure(figsize=(13,6))
	ax = fig.add_subplot(111)
	
	if(legend==None):
		legend = []
		
	rects = []
	# Generate one curve for each sample. Sample ids appear at samplelist in the same order that the legend strings, in other words, sample i in sample list is associated
	# with string i at legend list
	for sample in samplelist:
#		y[sample].sort()
#		y[sample].sort(key=lambda k:x[y[sample].index(k)])		
#		x[sample].sort()
		rects.append(ax.plot(x[sample], y[sample]))
		if(legend==None):
			legend.append(sample)
		
	ax.set_ylabel(ylabel)
	ax.set_xlabel('Depth')
	ax.set_ylim(top=100)
			
	if(len(legend)>1):
		# Shink current axis by 20%
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width, box.height*0.9])
		
		# Add graphic legend
		ax.legend( tuple([rect[0] for rect in rects]), tuple(legend), loc="lower left", bbox_to_anchor=(0,1.03) )
	
	
	fig.savefig(fileout)
	matplotlib.pyplot.close(fig)
	
	# Initialize the workbook and sheet
	wb = xlwt.Workbook()
	ws = wb.add_sheet('Saturation')

	# Create header font
	header_style = xlwt.easyxf('font: bold on')

	status = []
	slopes = []
	for j,sample in enumerate(x):
		ws.write(0, j*2, legend[j], header_style);
		ws.write(1, j*2, 'Depth', header_style); ws.write(1, (j*2)+1, ylabel, header_style);
		for i,value in enumerate(x[sample]): ws.write(i+2,j*2,value); ws.write(i+2,(j*2)+1,y[sample][i])
		
		slopes.append(((y[sample][-1]-y[sample][-2])*1.0/(x[sample][-1]-x[sample][-2])))
		status.append(slopes[-1]<=warnthreshold)
	
	wb.save(os.path.dirname(fileout)+'/values.xls')
	
	return status,slopes
	
	
def main():   
	################################################

	#### Options and arguments #####################

	################################################
	usage="""
	************************************************************************************************************************************************************
	Task: draws a coverage saturation plot with the data returned by simulated_depth.py  
	Output: a .png image containing the saturation plot.					 
	************************************************************************************************************************************************************
	
	
	Usage: %prog --filelist <filelist> --legend <legend> --ylabel <ylabel> --fileout <fileout>
	"""		 

	parser = optparse.OptionParser(usage)
	parser.add_option("--filelist", dest="filelist", help="""String containing comma separated names of the files that contain the percentaje of covered positions at each depth.""")
	parser.add_option("--ylabel", dest="ylabel", help="""String containing the label for the Y axis.""")	
	parser.add_option("--fileout", dest="fileout", help="""String containing the name of the file where the result """)
	
	(options, args) = parser.parse_args()


#	draw_saturation_curve(['/tmp/coverage_3442_1_agilib-agicap', '/tmp/coverage_3442_5_agilib-agicap','/tmp/coverage_3442_10_agilib-agicap','/tmp/coverage_3442_20_agilib-agicap','/tmp/coverage_3442_30_agilib-agicap'], ['Agilib-agicap.', '% covered positions', '/tmp/test.png')
	
	# Check number of arguments	
	if len(sys.argv) < 7:
		parser.print_help()
		sys.exit(1)

	# call core function
	draw_saturation_curve(options.filelist.split(','), options.ylabel, options.fileout)





if __name__=='__main__':
	main()

