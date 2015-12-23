#!/usr/bin/python
#-*- coding:utf-8 -*-
################################################
#File Name: callVariants.py
#Author: C.J. Liu
#Mail: samliu@hust.edu.cn
#Created Time: Wed 11 Nov 2015 08:01:40 PM CST
################################################



#using GATK and lumpy method
import os,sys
import argparse
import re
import multiprocessing

root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
software=root + "/software"
data = root + "/data"

soft=software + "/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar"
mill=data + "/hg19/Mills_and_1000G_gold_standard.indels.hg19.sorted.vcf"
kg=data + "/hg19/1000G_phase1.indels.hg19.sorted.vcf"
dbsnp=data + "/hg19/dbsnp_137.hg19.sorted.vcf"
index = data + "/hg19/genomeBuild/hg19.fasta"

omni=data + "/hg19/1000G_omni2.5.hg19.vcf"
hapmap=data + "/hg19/hapmap_3.3.hg19.vcf"
def usage():
	################################################
	#### Options and arguments #####################
	################################################
	
	description="""	
	Task: Call variants with GATK.
	Output: The script was writted to deal with fewer samples without using complete GATK Best Practice,like VQSR can not run in terms of less than 30 samples.
	"""
	
	usage = """ %(prog)s -pe1 <fq1> -pe2 <fq2> -i <pwd> -o <pwd>"""
	
	parser = argparse.ArgumentParser(description = description,usage = usage)
	parser.add_argument("-b", dest="bam", type=str, help="""Required. Input recalibrated bam file with bai index""",required=True)
	parser.add_argument("-i", dest="indir", type=str, help="""Specify input directory, default is current directory""",default=os.getcwd())
	parser.add_argument("-o", dest="out", type=str, help="""Specify output directory, default is current directory""",default=os.getcwd())
	parser.add_argument("-idx", dest="index", type=str, help="""GenomeBuild index default hg19""",default=index)
	parser.add_argument("-t", dest="nthreads", type=int, help="""Optional. Integer indicating the number of concurrent threads to launch. Default=10.""", default=20)
	parser.add_argument("-s",dest="soft", help="""Specify the path and version of GATK""", default=soft)
	parser.add_argument("-m",dest="mill", help="""Specify mills path for VQSR""", default=mill)
	parser.add_argument("-k",dest="kg", help="""Specify 1000 indels for VQSR""", default=kg)
	parser.add_argument("-d",dest="dbsnp", help="""Specify dbsnp path for VQSR, default version is 137""", default=dbsnp)
	parser.add_argument("-om",dest="omni", help="""Specify omni path""", default=omni)
	parser.add_argument("-hap",dest="hapmap", help="""Specify hapmap path""", default=hapmap)
	parser.add_argument('-v','--version',action='version', version='%(prog)s 1.0')
	args = parser.parse_args()
	
	return args

def gatkCall(bam,indir,out,index=index,soft=soft,t=20,mill=mill,kg=kg,dbsnp=dbsnp,omni=omni,hapmap=hapmap):
	out = out + os.sep +"CallVariants"
	inbam=indir + os.sep + bam
	try:
		os.mkdir(out)
	except:
		print "Directory %s already exists" % out
	
	#Call variants with HaplotypeCaller
	rawvcf = re.split(r'\_|\.', bam)[0] + '.raw_variants.vcf'
	outrawvcf = out + os.sep + rawvcf	
	# is limited to contig chr20, it should be removed when calling genome variants
	cmd = "java -Xmx50g -Djava.io.tmpdir=/tmp  -jar %s -T HaplotypeCaller -R %s -I %s  --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o %s -nct %s" %(soft,index,inbam, outrawvcf, str(t))
	# print cmd
	os.system(cmd)
	###############################################
	#Recalibrate variant quality scores = run VQSR
	###############################################
	####not suitable for small samples
	###############################################
	#snp
	# snprecal = out + os.sep + 'recalibrate_SNP.recal'
	# snptranches = out + os.sep + "recalibrate_SNP.tranches"
	# snprscript = out + os.sep + "recalibrate_SNP_plots.R"
	# cmd = "java -Xmx50g -Djava.io.tmpdir=/tmp  -jar %s -T VariantRecalibrator -R %s -input %s -resource:hapmap,known=false,training=true,truth=true,prior=15.0 %s  -resource:omni,known=false,training=true,truth=true,prior=12.0 %s -resource:1000G,known=false,training=true,truth=false,prior=10.0 %s -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 %s -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum --maxGaussians 4 -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile %s -tranchesFile %s -rscriptFile %s -nt %s" %(soft, index, outrawvcf, hapmap, omni, kg, dbsnp, snprecal, snptranches,snprscript,t)
	# print cmd
	# os.system(cmd)
	
	# recalsnprawindel = out + os.sep + rawvcf.rstrip('vcf') + "recal.snps.raw.indels.vcf"
	# cmd = "java -Xmx50g -Djava.io.tmpdir=/tmp  -jar %s -T ApplyRecalibration -R %s -input %s -mode SNP --ts_filter_level 99.5 -recalFile %s -tranchesFile %s -o %s" %(soft,index,outrawvcf,snprecal,snptranches,recalsnprawindel)
	# os.system(cmd)
	
	#indel
	# indelrecal = out + os.sep + 'recalibrate_INDEL.recal'
	# indeltranches = out + os.sep + "recalibrate_INDEL.tranches"
	# indelrscript = out + os.sep + "recalibrate_INDEL_plots.R"
	# cmd = "java -Xmx50g -Djava.io.tmpdir=/tmp  -jar %s -T VariantRecalibrator -R %s -input %s -resource:mills,known=false,training=true,truth=true,prior=12.0 %s -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 %s -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum  -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile %s -tranchesFile %s -rscriptFile %s -nt %s" %(soft, index, recalsnprawindel, mill, dbsnp,indelrecal, indeltranches, indelrscript, t)
	# os.system(cmd)
	
	# recalvariants = out + os.sep + rawvcf.rstrip('vcf') + "recal.variants.vcf"
	# cmd = "java -Xmx50g -Djava.io.tmpdir=/tmp  -jar %s -T ApplyRecalibration -R %s -input %s -mode INDEL --ts_filter_level 99.0 -recalFile %s -tranchesFile %s -o %s" %(soft,index,recalsnprawindel,indelrecal,indeltranches,recalvariants)
	# os.system(cmd)
	##########
	#Apply hard filters to call set
	##########
	def snprun():
		#For select SNP
		snp = rawvcf.rstrip('vcf')+'SNP.vcf'
		outsnp = out + os.sep + snp
		cmd = "java -Xmx50g -Djava.io.tmpdir=/tmp  -jar %s -T SelectVariants -R %s -V %s  -selectType SNP -o %s" %(soft,index,outrawvcf,outsnp)
		os.system(cmd)
		#For filter SNP
		snpfilter = snp.rstrip('vcf') + 'filter.vcf'
		outsnpfilter = out + os.sep + snpfilter
		cmd = "java -Xmx50g -Djava.io.tmpdir=/tmp  -jar %s -T VariantFiltration -R %s -V %s --filterExpression \"DP < 8 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0\" --filterName GATK_snp_filter -o %s " %(soft,index,outsnp,outsnpfilter)
		os.system(cmd)
		
	def indelrun():
		#For select INDEL
		indel = rawvcf.rstrip('vcf')+'INDEL.vcf'
		outindel = out + os.sep + indel
		cmd = "java -Xmx50g -Djava.io.tmpdir=/tmp  -jar %s -T SelectVariants -R %s -V %s  -selectType INDEL -o %s" %(soft,index,outrawvcf,outindel)
		os.system(cmd)	
		#For filter INDEL
		indelfilter = indel.rstrip('vcf') + 'filter.vcf'
		outindelfilter = out + os.sep + indelfilter
		cmd = "java -Xmx50g -Djava.io.tmpdir=/tmp  -jar %s -T VariantFiltration -R %s -V %s  --filterExpression \"DP < 8 || QD < 3.0 || FS > 200.0 || ReadPosRankSum < -20.0\" --filterName GATK_indel_filter -o %s " %(soft,index,outindel,outindelfilter)
		os.system(cmd)
	
	jobs = []
	
	p = multiprocessing.Process(target=snprun)
	p.start()
	jobs.append(p)
	p = multiprocessing.Process(target=indelrun)
	p.start()
	jobs.append(p)
	for it in jobs:
		it.join()
	
	return rawvcf, out
	print "********************variation calling done!!!*************"
	
def main():
	args = usage()
	print args
	gatkCall(args.bam,args.indir,args.out,index=args.index,soft = args.soft,t = args.nthreads,mill = args.mill, kg = args.kg,dbsnp=args.dbsnp,omni=args.omni,hapmap=args.hapmap)
	
	
if __name__ == "__main__":
	main()










