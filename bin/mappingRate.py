#!/usr/bin/python
#-*- coding:utf-8 -*-
################################################
#File Name: mappingAndDedup.py
#Author: C.J. Liu
#Mail: samliu@hust.edu.cn
#Created Time: Wed 11 Nov 2015 07:59:51 PM CST
################################################


import os,sys
import argparse
import os.path
import re

root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
software=root + "/software"
data = root + "/data"
index = data + "/hg19/genomeBuild/hg19.fasta"
soft = software + "/ngscat.v0.1/ngscat.py"
cap = data + "/hg19/SureSelect_All_Exon_50mb_with_annotation_hg19_bed_50M.bed"

def usage():
	################################################
	#### Options and arguments #####################
	################################################
	
	description="""	
	Task: Map paired reads to reference genome with bwa-0.7.12, default genome build is hg19.Sort sam files, mark and remove duplicates with picard-1.100
	Output: You will get deduped bam files with its bai index.
	"""
	
	usage = """ %(prog)s -pe1 <fq1> -pe2 <fq2> -i <pwd> -o <pwd>"""
	
	parser = argparse.ArgumentParser(description = description,usage = usage)
	parser.add_argument("-bam", dest="bam", type=str, help="""Required. Input sorted bam file""",required=True)
	parser.add_argument("-i", dest="indir", type=str, help="""Specify input directory. Default is current directory""",default=os.getcwd())
	parser.add_argument("-o", dest="out", type=str, help="""Specify output directory. Default is current directory""",default=os.getcwd())
	parser.add_argument("-idx", dest="index", type=str, help="""GenomeBuild index, default hg19""",default=index)
	parser.add_argument("-cap", dest="cap", type=str, help="""Exome capture file, default hg19 agilent""",default=cap)
	parser.add_argument("-t", dest="nthreads", type=int, help="""Optional. Integer indicating the number of concurrent threads to launch. Default=10.""", default=10)
	parser.add_argument("-s",dest="soft", help="""Specify the location of mapping software""", default=soft)
	parser.add_argument('-v','--version',action='version', version='%(prog)s 1.0')
	args = parser.parse_args()
	
	return args

def mappingRate(bam,indir,out,index=index,cap=cap,soft=soft,t=10):
	## IlluQC_PRLL get quality control
	out = out + os.sep +"MappingRate"
	result = out + os.sep + "mappingrate"
	tmpdir = out + os.sep + "tmp"
	try:
		os.mkdir(out)
		os.mkdir(result)
		print "Create output directory..."
	except:
		print "Directory %s already exists" % out
	try:
		os.mkdir(tmpdir)
		print "Create tmp directory..."
	except:
		print "Directory %s already exists" % tmpdir
	
	inputbam = indir + os.sep + bam
	#ngscat.py --bams bam --bed seqcap.bed --out outdir --reference index --saturation y --depthlist 0.01,0.02,0.025,0.05,0.075,0.1,0.2,0.3,0.4,0.5 --threads 10 --tmp tmpdir
	#Mapping with mem algorithm
	cmd = "python %s --bams %s --bed %s --out %s --reference %s --saturation y --depthlist 0.01,0.02,0.025,0.05,0.075,0.1,0.2,0.3,0.4,0.5 --threads %s --tmp %s" % (soft, inputbam, cap, result, index, t, tmpdir)
	os.system(cmd)
	print "***************mapping Rate done!!!***************"
	
	
def main():
	
	args = usage()	
	mappingRate(args.bam,args.indir,args.out,index = args.index,cap = args.cap,soft = args.soft,t=args.nthreads)

if __name__ == "__main__":
	main()





