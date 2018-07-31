#########################################################################
# File Name: test.sh
# Author: C.J. Liu
# Mail: samliu@hust.edu.cn
# Created Time: Thu 19 Nov 2015 04:50:58 PM CST
#########################################################################
#!/bin/bash
mkdir process;mkdir result
nohup python ../wes_analysis.pyc -pe1 NA12878.hiseq.wgs_chr20_2mb.30xPE_1.fastq -pe2 NA12878.hiseq.wgs_chr20_2mb.30xPE_2.fastq -i $PWD/sample -o $PWD/process -r hg19 -c agilent -t 20 &
