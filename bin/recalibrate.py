#!/usr/bin/python
#-*- coding:utf-8 -*-
################################################
#File Name: recalibrate.py
#Author: C.J. Liu
#Mail: samliu@hust.edu.cn
#Created Time: Thu 12 Nov 2015 12:38:16 AM CST
################################################


import os,sys
import argparse

root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
software=root + "/software"
data = root + "/data"

soft=software + "/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar"
mill=data + "/hg19/Mills_and_1000G_gold_standard.indels.hg19.sorted.vcf"
kg=data + "/hg19/1000G_phase1.indels.hg19.sorted.vcf"
dbsnp=data + "/hg19/dbsnp_137.hg19.sorted.vcf"
index = data + "/hg19/genomeBuild/hg19.fasta"
def usage():
	################################################
	#### Options and arguments #####################
	################################################
	
	description="""	
	Task: BQSQ with GATK v3.4-46, some walkers deprecated are not used anymore.
	Output: From deduped file to generate recalibrated file with bai index.
	"""
	
	usage = """ %(prog)s -pe1 <fq1> -pe2 <fq2> -i <pwd> -o <pwd>"""
	
	parser = argparse.ArgumentParser(description = description,usage = usage)
	parser.add_argument("-b", dest="bam", type=str, help="""Required. Input deduped bam file with bai index""",required=True)
	parser.add_argument("-i", dest="indir", type=str, help="""Specify input directory, default direcotry is current directory""",default=os.getcwd())
	parser.add_argument("-o", dest="out", type=str, help="""Specify output directory, default is current directory""",default=os.getcwd())
	parser.add_argument("-idx", dest="index", type=str, help="""GenomeBuild index default hg19""",default=index)
	parser.add_argument("-t", dest="nthreads", type=int, help="""Optional. Integer indicating the number of concurrent threads to launch. Default=10.""", default=10)
	parser.add_argument("-s",dest="soft", help="""Specify the location of mapping software""", default=soft)
	
	parser.add_argument("-m",dest="mill", help="""Specify the location of mapping software""", default=mill)
	parser.add_argument("-k",dest="kg", help="""Specify the location of mapping software""", default=kg)
	parser.add_argument("-d",dest="dbsnp", help="""Specify the location of mapping software""", default=dbsnp)
	parser.add_argument('-v','--version',action='version', version='%(prog)s 1.0')
	args = parser.parse_args()
	
	return args

def realignAndrecal(bam,indir,out,index=index,soft=soft,t=10,mill=mill,kg=kg,dbsnp=dbsnp):
	
	intervals = out + os.sep + bam + '.intervals'
	inbam = out + os.sep + bam
	#CREATE INTERVAL
	cmd = "java -Xmx50g -jar %s -T RealignerTargetCreator -R %s -o %s -I %s -known %s -known %s" %(soft,index,intervals,inbam,mill,kg)
	os.system(cmd)

	#REALIGNMENT
	realign = bam.rstrip('bam') + 'realign.bam'
	outrealign = out + os.sep + realign
	cmd = "java -Xmx50g -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s -known %s -known %s --filter_bases_not_stored" %(soft,index,inbam,intervals,outrealign,mill,kg)
	# print cmd
	os.system(cmd)
	
	#RECALIBRATION
	recal = realign.rstrip("bam") + 'recal.bam'
	grp = out+os.sep + recal + '.grp'
	
	cmd="java -Xmx50g -Djava.io.tmpdir=/tmp -jar %s -T BaseRecalibrator -R %s -I %s -knownSites %s -knownSites %s -knownSites %s -rf BadCigar -o %s" %(soft,index,outrealign,dbsnp,mill,kg,grp)
	# print cmd
	os.system(cmd)
	
	outrecal = out + os.sep + recal
	cmd = "java -Xmx50g -Djava.io.tmpdir=/tmp -jar %s -T PrintReads -R %s -I %s -BQSR %s -o %s" %(soft,index,outrealign,grp,outrecal)
	# print cmd
	os.system(cmd)
	
	#ReduceReads is deprecated in 3.0+, ReduceReads is used for UG calling for low quality variants
	# reduce = recal.rstrip('bam') + 'reduce.bam'
	# outreduce = out + os.sep + reduce
	# cmd = "java -Xmx50g -Djava.io.tmpdir=/tmp/ -jar %s -T ReduceReads -R %s -I %s -o %s" %(soft,index,outrecal,outreduce)
	# print cmd
	os.system(cmd)
	print "***************Mapping recalibration is done!!!************"
	
	return recal
	
def main():
	args = usage()
	realignAndrecal(args.bam,args.indir,args.out,index=args.index,soft = args.soft,t = args.nthreads,mill = args.mill, kg = args.kg,dbsnp=args.dbsnp)

if __name__ == '__main__':
	main()
