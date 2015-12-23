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
soft = software + "/bwa-0.7.12/bwa"
sortsam = software + "/picard-tools-1.100/SortSam.jar"
mark = software + "/picard-tools-1.100/MarkDuplicates.jar"
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
	parser.add_argument("-pe1", dest="pe1", type=str, help="""Required. Input read_pe1""",required=True)
	parser.add_argument("-pe2", dest="pe2", type=str, help="""Required. Input reead_pe2""",required=True)
	parser.add_argument("-i", dest="indir", type=str, help="""Specify input directory. Default is current directory""",default=os.getcwd())
	parser.add_argument("-o", dest="out", type=str, help="""Specify output directory. Default is current directory""",default=os.getcwd())
	parser.add_argument("-idx", dest="index", type=str, help="""GenomeBuild index, default hg19""",default=index)
	parser.add_argument("-t", dest="nthreads", type=int, help="""Optional. Integer indicating the number of concurrent threads to launch. Default=10.""", default=10)
	parser.add_argument("-s",dest="soft", help="""Specify the location of mapping software""", default=soft)
	parser.add_argument("-so",dest="sortsam", help="""Specify the location of mapping software""", default=sortsam)
	parser.add_argument("-m",dest="mark", help="""Specify the location of mapping software""", default=mark)
	parser.add_argument('-v','--version',action='version', version='%(prog)s 1.0')
	args = parser.parse_args()
	
	return args

def mapping(fq1,fq2,indir,out,index,soft=soft,t=10):
	## IlluQC_PRLL get quality control
	out = out + os.sep +"Mapping"
	
	try:
		os.mkdir(out)
	except:
		print "Directory %s already exists" % out
	re1 = indir + os.sep + fq1
	re2 = indir + os.sep + fq2
	sam = re.split(r'\_|\.',fq1)[0] + '.sam'
	
	readgroup = r'"@RG\tID:%s\tLB:%s\tSM:%s\tPL:ILLUMINA"' %(sam, sam, sam)
	#Mapping with mem algorithm
	cmd = soft + " mem -t " + str(t) + " -M -a " + index + " " + re1 + " " + re2 + ' -R ' + readgroup + ">" + out + os.sep + sam
        print cmd
	os.system(cmd)
	print "***************mapping done!!!***************"
	return sam, out

def dedup(sam, out, sortsam=sortsam, mark=mark):
	#sort
	outsam = out + os.sep + sam
	bam = sam + '.bam'
	outbam = out + os.sep + bam
	cmd = "java -Xmx50g -Djava.io.tmpdir=/tmp -jar " + sortsam + " SO=coordinate INPUT=" +  outsam + " OUTPUT=" + outbam + " VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true"
	os.system(cmd)
	#dedup
	dedup = sam + '.dedup.bam'
	outdedup = out + os.sep + dedup
	metrics = out + os.sep + "mark.metrics"
	cmd = "java -Xmx50g -Djava.io.tmpdir=/tmp -jar " + mark + " INPUT=" + outbam + " OUTPUT=" + outdedup + " METRICS_FILE=" + metrics + "  MAX_RECORDS_IN_RAM=5000000 ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true CREATE_INDEX=true"
	os.system(cmd)
	print "***************Dedup done!!!*********************"
	return dedup
	
def main():
	
	args = usage()	
	sam,out = mapping(args.pe1,args.pe2,args.indir,args.out,args.index,soft = args.soft,t=args.nthreads)
	dedup(sam, out, sortsam = args.sortsam, mark = args.mark)

if __name__ == "__main__":
	main()





