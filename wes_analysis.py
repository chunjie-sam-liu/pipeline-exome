#!/usr/bin/python
#-*- coding:utf-8 -*-
################################################
#File Name: genomeResequencing.py
#Author: C.J. Liu
#Mail: samliu@hust.edu.cn
#Created Time: Wed 11 Nov 2015 07:56:27 PM CST
################################################

import __future__
import os,sys
import argparse
import json
import os.path
import pprint
import multiprocessing
root = os.path.dirname(os.path.realpath(__file__))

sys.path.insert(0,root + os.sep + 'bin')
sys.path.insert(0,root + os.sep + 'script')

import qualityQcontol
import mappingAndDedup
import recalibrate
import callVariants
import annotateWithAnnovar
import mappingRate
def usage():
	################################################
	#### Options and arguments #####################
	################################################
	
	description="""	
	Task: Whole genome resequencing data analysis.
	Specification:
			Detailed software version and reference data are list in the config file
	Output:
		The software will create server files cooresponding task and results.
	"""
	
	usage = """ %(prog)s -pe1 <fq1> -pe2 <fq2> -i <pwd> -o <pwd>"""
	
	parser = argparse.ArgumentParser(description = description,usage = usage)
	parser.add_argument("-pe1", dest="pe1", type=str, help="""Required. Input read_pe1""",required=True)
	parser.add_argument("-pe2", dest="pe2", type=str, help="""Required. Input reead_pe2""",required=True)
	parser.add_argument("-i", dest="indir", type=str, help="""Specify input directory""",default=os.getcwd())
	parser.add_argument("-o", dest="out", type=str, help="""Specify output directory""",default=os.getcwd())
	parser.add_argument("-r", dest="reference", type=str, help="""Specify genome assembly hg19|hg38,default=hg19""",default="hg19")
	parser.add_argument("-c", dest="capture", type=str, help="""Specify capture chip sequencing used,default hg19 agilent 50m;optional agilent|illumina|roche""",default="agilent")
	parser.add_argument("-t", dest="nthreads", type=int, help="""Optional. Integer indicating the number of concurrent threads to launch. Default=10.""", default=10)
	parser.add_argument('-v','--version',action='version', version='%(prog)s 1.0')
	args = parser.parse_args()

	return args

def config(build):
	try:
		with open(root + os.sep + 'config','r') as foo:
			print root + os.sep + 'config'
			config = json.load(foo)
			print "loading config file........"
	except :
		print "Can't load config file, config file is required"
		sys.exit(1)
	
	#add root
	#For software
	for key,val in config['software'].items():	
		config['software'][key] = root + os.sep + val
	#for reference
	for key,val in config['reference'][build].items():
		config['reference'][build][key] = root + os.sep + val
	
	#check all file
	
		
	return config

def main():
	args = usage()
	
	conf = config(args.reference)
	soft = conf['software']
	ref = conf['reference'][args.reference]
	script = conf['script']
	# pprint.pprint(conf)
	
	#Test for qualityQcontol
	multiprocessing.Process(target=qualityQcontol.run, args=(args.pe1,args.pe2,args.indir,args.out,soft['IllQC'],args.nthreads)).start()
	
	#Test mapping
	sam,mappingOut=mappingAndDedup.mapping(args.pe1,args.pe2,args.indir,args.out,ref['genomeBuild'],soft=soft['bwa'],t=args.nthreads)
	dedup = mappingAndDedup.dedup(sam, mappingOut, sortsam = soft['sortsam'], mark = soft['mark'])
	
	#Test mappingRate
	bam = sam + '.bam'
	#mappingRate.mappingRate(bam,mappingOut,args.out,index = ref['genomeBuild'],cap = ref[args.capture],soft = soft['ngscat'],t=args.nthreads)
	#In multiprocessing method
	multiprocessing.Process(target=mappingRate.mappingRate,args=(bam,mappingOut,args.out,ref['genomeBuild'],ref[args.capture],soft['ngscat'],args.nthreads)).start()
	
	#Test recalibrate
	recal = recalibrate.realignAndrecal(dedup,mappingOut,mappingOut,index=ref['genomeBuild'],soft = soft['gatk'],t = args.nthreads,mill = ref['mill'], kg = ref['kg'],dbsnp=ref['dbsnp'])
		
	#Test variant calling
	vcf,variantout = callVariants.gatkCall(recal,mappingOut,args.out,index=ref['genomeBuild'],soft = soft['gatk'],t = args.nthreads,mill = ref['mill'], kg = ref['kg'],dbsnp=ref['dbsnp'],omni=ref['omni'],hapmap=ref['hapmap'])
	
	#Test annotation
	annotateWithAnnovar.annotate(vcf,variantout,args.out,hd=ref['humandb'],con=soft['convert2annovar'],tb=soft['table_annovar'])
	
	

if __name__ == '__main__':
	main()
	
	
	
	
	
