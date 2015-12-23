#!/usr/bin/python
#-*- coding:utf-8 -*-
################################################
#File Name: annotateWithAnnovar.py
#Author: C.J. Liu
#Mail: samliu@hust.edu.cn
#Created Time: Wed 11 Nov 2015 08:00:59 PM CST
################################################

import os,sys
import os.path
import re
import argparse
import multiprocessing

root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
software=root + "/software"
data = root + "/data"

humandb=data + "/hg19/humandb19"
convert = software + "/annovar/convert2annovar.pl"
tableannovar = software + "/annovar/table_annovar.pl"

def usage():
	################################################
	#### Options and arguments #####################
	################################################
	
	description="""	
	Task: Annotate variants with annovar.
	Output: We have our eye on the multianno.txt file, this is the collected result.
	"""
	
	usage = """ %(prog)s -pe1 <fq1> -pe2 <fq2> -i <pwd> -o <pwd>"""
	
	parser = argparse.ArgumentParser(description = description,usage = usage)
	parser.add_argument("-vcf", dest="vcf", type=str, help="""Required. Input vcf file.""",required=True)
	parser.add_argument("-i", dest="indir", type=str, help="""Specify input directory""",default=os.getcwd())
	parser.add_argument("-o", dest="out", type=str, help="""Specify output directory""",default=os.getcwd())
	parser.add_argument("-hd", dest="humandb", type=str, help="""Annovar annotation files""",default=humandb)
	parser.add_argument("-con", dest="convert", type=str, help="""Annovar convert vcf file to avinput""",default=convert)
	parser.add_argument("-tb", dest="tableannovar", type=str, help="""Annovar table annovar to annotate""",default=tableannovar)
	parser.add_argument('-v','--version',action='version', version='%(prog)s 1.0')
	args = parser.parse_args()
	
	return args

def run(vcf,out,hd,con,tb):
	avinput = os.path.basename(vcf) +  ".avinput"
	outavinput = out + os.sep +  avinput
	cmd = "perl %s -format vcf4 --includeinfo %s > %s" %(con,vcf,outavinput)
	os.system(cmd)
	buildver = 'hg' + str(re.search(r'(\d+)',hd).group(1))
	cmd = "perl %s %s %s -buildver %s -protocol refGene,ensGene,knownGene,ccdsGene,sift,cosmic65,esp6500si_all,snp137 -operation g,g,g,g,f,f,f,f -nastring . -outfile %s" %(tb,outavinput,hd,buildver,outavinput)
	os.system(cmd)
	
	# perl $SCRIPT/maf_stat.pl $ANO_SNP_OUTFILE/$FINAL_NAME.snp.hg19_multianno.txt 

	# perl $SCRIPT/exome_demand.pl -maf $ANO_SNP_OUTFILE/$FINAL_NAME.snp.hg19_multianno.txt.maf -avi $ANO_SNP_OUTFILE/$FINAL_NAME.snp.final.vcf.avinput
	
def annotate(vcf,indir,out,hd=humandb,con=convert,tb=tableannovar):
	out = out + os.sep +"Annotation"
	snp = indir + os.sep + vcf.rstrip('vcf') + "SNP.filter.vcf"
	indel = indir + os.sep + vcf.rstrip('vcf') + "INDEL.filter.vcf"
	snpdir = out + os.sep + "SNP"
	indeldir = out + os.sep + "INDEL"
	try:
		os.mkdir(out)
	except:
		print "Directory %s already exists" % out
	try:
		os.mkdir(snpdir)
	except:
		print "Directory %s already exists" % snpdir
	try:
		os.mkdir(indeldir)
	except:
		print "Directory %s already exists" % indeldir
	# run(snp,snpdir,hd,con,tb)
	multiprocessing.Process(target=run,args=(snp,snpdir,hd,con,tb,)).start()
	multiprocessing.Process(target=run,args=(indel,indeldir,hd,con,tb,)).start()
	print "*********************Anotation done!!!********************"


def main():
	args = usage()
	annotate(vcf,indir,out,hd,con,tb)
	


if __name__ == "__main__":
	main()

