#!/usr/bin/python
#-*- coding:utf-8 -*-
################################################
#File Name: callSV.py
#Author: C.J. Liu
#Mail: samliu@hust.edu.cn
#Created Time: Wed 11 Nov 2015 08:06:06 PM CST
################################################

#use crest
#Fucking CREST

import os,sys
import argparse
import os.path
import re



def usage():
	################################################
	#### Options and arguments #####################
	################################################
	
	description="""	
	Task: Assesses capture performance in terms of sensibility, specificity and uniformity of the coverage.
	Output: An html report will be created at the path indicated with the --out option.
	"""
	
	usage = """ %(prog)s -pe1 <fq1> -pe2 <fq2> -i <pwd> -o <pwd>"""
	
	parser = argparse.ArgumentParser(description = description,usage = usage)
	parser.add_argument("-b", dest="bam", type=str, help="""Required. Input read_pe1""",required=True)
	parser.add_argument("-i", dest="indir", type=str, help="""Specify input directory""",default=os.getcwd())
	parser.add_argument("-o", dest="out", type=str, help="""Specify output directory""",default=os.getcwd())
	parser.add_argument("-idx", dest="index", type=str, help="""GenomeBuild index default hg19""",default="/project/liucj/REFDATA/genome_build/hg19-bwa-index/hg19.fasta")
	parser.add_argument("-s",dest="softdir", help="""Specify the location of mapping software""", default="/home/liucj/tools/crest-20120209")
	parser.add_argument('-v','--version',action='version', version='%(prog)s 1.0')
	args = parser.parse_args()
	
	return args

def crest(bam,indir,out,index,softdir):
	cmd = "export PATH=%s/:$PATH" % softdir
	os.system(cmd)
	extract = "extractSClip.pl"
	crest = "CREST.pl"
	bam2html = "bam2html.pl"
	cap3="cap3"
	blat="blat"
	# gfserver = softdir + "/gfServer"
	# cmd = "nohup %s start localhost 6666 %s > server.log &" %(gfserver,fa2bit)
	# os.system(cmd)
	#it was running on the server with root
	#gfServer start localhost 6666 hg19.2bit
	inbam = out + os.sep + bam
	# perl extract -i bam --ref_genome index
	cmd = "%s -i %s --ref_genome %s" %(extract, inbam,index)
	os.system(cmd)
	# perl crest -f bam.cover -d bam --ref_genome index --2bitdir os.path.dir(fa2bit) -t fa2bit  --blatserver localhost --blatport 6666 --cap3 cap3 --blat blat
	# perl bam2html -d bam -i predsv --ref_genome index -o html
	
	print 'cj'
	
def main():
	args = usage()
	crest(args.bam,args.indir,args.out,args.index,args.softdir)
	
if __name__ =="__main__":
	main()
	


