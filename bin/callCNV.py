#!/usr/bin/python
#-*- coding:utf-8 -*-
################################################
#File Name: callCNV.py
#Author: C.J. Liu
#Mail: samliu@hust.edu.cn
#Created Time: Wed 11 Nov 2015 08:06:42 PM CST
################################################

import os,sys
import argparse
import os.path
import re

index="/project/liucj/REFDATA/genome_build/hg19-bwa-index/hg19.fasta"
soft="/home/liucj/piplines/resequencing/genome_resequencing/software/CNVnator_v0.3.2"
def usage():
	################################################
	#### Options and arguments #####################
	################################################
	
	description="""	
	Task: Assesses capture performance in terms of sensibility, specificity and uniformity of the coverage.
	Output: An html report will be created at the path indicated with the --out option.
	Before using Calling CNV with cnvnator, be sure that "cern root"-data analysis framework file in the CNVnator path,
	"""
	
	usage = """ %(prog)s -pe1 <fq1> -pe2 <fq2> -i <pwd> -o <pwd>"""
	
	parser = argparse.ArgumentParser(description = description,usage = usage)
	parser.add_argument("-b", dest="bam", type=str, help="""Required. Input read_pe1""",required=True)
	parser.add_argument("-i", dest="indir", type=str, help="""Specify input directory""",default=os.getcwd())
	parser.add_argument("-o", dest="out", type=str, help="""Specify output directory""",default=os.getcwd())
	parser.add_argument("-idx", dest="index", type=str, help="""GenomeBuild index default hg19""",default=index)
	parser.add_argument("-s",dest="softdir", help="""Specify the location of mapping software""", default=soft)
	parser.add_argument('-v','--version',action='version', version='%(prog)s 1.0')
	args = parser.parse_args()
	
	return args

	
def callcnv(bam,indir, out,index=index,soft=soft):
	#Set cernroot environment for cnv.
	#The  cnvnator root out is a cern root data format
	cernRoot = soft + os.sep + "root"
	thisroot = cernRoot + "/bin/thisroot.sh"
	cmd = "source %s" %thisroot
	# print thisroot
	os.system(cmd)
	cnvnator = soft + "/src/cnvnator"
	out = out + os.sep + "CallCNV"
	try:
		os.mkdir(out)
	except:
		print "Direcotry %s already exists" %out
	
	print "cj"
	
def main():
	args = usage()
	callcnv(args.bam, args.indir,args.out,index=args.index,soft=args.softdir)
	


if __name__ == "__main__":
	main()




