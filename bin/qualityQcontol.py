#!/usr/bin/python
#-*- coding:utf-8 -*-
################################################
#File Name: qualityQcontol.py
#Author: C.J. Liu
#Mail: samliu@hust.edu.cn
#Created Time: Wed 11 Nov 2015 07:17:36 PM CST
################################################

import os,sys
import argparse
import os.path

root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
software=root + "/software"
data = root + "/data"

soft= software + "/NGSQCToolkit_v2.3.3/QC/IlluQC_PRLL.pl"

def usage():
	################################################
	#### Options and arguments #####################
	################################################
	
	description="""	
	Task: Trim and evaluate quality of sequencing reads.
	Output: It will create a QualityControl directory containing all result, the html file is the final result.
	"""
	
	usage = """ %(prog)s --pe1 <fq1> --pe2 <fq2>"""
	
	parser = argparse.ArgumentParser(description = description,usage = usage)
	parser.add_argument("-pe1", dest="pe1", type=str, help="""Required. Input read_pe1""",required=True)
	parser.add_argument("-pe2", dest="pe2", type=str, help="""Required. Input reead_pe2""",required=True)
	parser.add_argument("-o", dest="out", type=str, help="""Required. Specify output directory.Default is current directory""",default=os.getcwd())
	parser.add_argument("-i", dest="indir", type=str, help="""Required. Specify input directory.Default is current directory""",default=os.getcwd())
	parser.add_argument("-t", dest="nthreads", type=int, help="""Optional. Integer indicating the number of concurrent threads to launch. Default=10.""", default=10)
	parser.add_argument("-s",dest="soft", help="""Specify the location of NGSQC software""", default=soft)
	parser.add_argument('-v','--version',action='version', version='%(prog)s 1.0')
	args = parser.parse_args()
	
	return args

def run(fq1,fq2,indir,out,soft=soft,t=10):
	## IlluQC_PRLL get quality control
	out = out + os.sep +"QualityControl"
	
	try:
		os.mkdir('%s' % out)
	except:
		print "Directory %s already exists" % out
	# print "perl %s -pe %s %s 2 A -c %s -l 70 -t 2 -o %s -z g" %(soft, fq1, fq2, t,out)
	fq1 = indir + os.sep + fq1
	fq2 = indir + os.sep + fq2
	os.system("perl %s -pe %s %s 2 A -c %s -l 70 -t 2 -o %s -z g" %(soft, fq1, fq2, t,out))
	print "******************Quality control done!!!*******************"

def main():
	
	args = usage()	
	run(args.pe1,args.pe2,args.indir,args.out,args.soft,args.nthreads)
	

if __name__ == "__main__":
	main()




