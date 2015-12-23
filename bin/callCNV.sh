#########################################################################
# File Name: callCNV.sh
# Author: C.J. Liu
# Mail: samliu@hust.edu.cn
# Created Time: Thu 19 Nov 2015 03:41:42 PM CST
#########################################################################
#!/bin/bash
function usage(){
	echo "bash $0 -b <recal_bam> -i <indir> -o <output dir> -r <ref>"
	echo "Description:"
	echo "	Call Copy number variation with CNVnator"
	echo "	CNVnator depends on the cern root"
	echo "	Before using Calling CNV with cnvnator, be sure that \"cern root\"-data analysis framework file in the CNVnator path,"
	echo "	"
	return 1
}

function main(){
	recal=$1
	indir=$2
	out=$3
	ref=$4
	test -d $out/CallCNV || mkdir $out/CallCNV
	out=$out/CallCNV
	
	root=$(cd "$(dirname "$0")";cd ../;pwd)
	data=${root}/data
	software=${root}/software
	
	soft=${software}/CNVnator_v0.3.2
	indexdir=${data}/$ref/genomeBuild
	index=${indexdir}/${ref}.fasta
	test -e "$index" || { echo "File $index doesn't exist"; exit 1; }
	
	#Set cern root environment before running CNVnator
	cernroot=${soft}/root
	test -d $cernroot || { echo "CERN ROOT must be in the $soft"; exit 1; }
	thisroot=${cernroot}/bin/thisroot.sh
	source $thisroot
	
	#Set CNVnator path
	export PATH=${soft}/:${soft}/src/:$PATH
	#cnvnator
	
	test -L ${out}/${recal} || ln -s ${indir}/${recal} ${out}/${recal} 
	#predict with CNV region
	#Extract read mapping from bam/sam files
	cnvnator -genome $index -root ${out}/${recal}.root -tree ${out}/${recal}

	#Generate a histogram
	cnvnator -genome $index -root ${out}/${recal}.root -his 1000
	#Calculate statistics
	cnvnator -genome $index -root ${out}/${recal}.root -stat 1000
	#RD signal partition
	cnvnator -genome $index -root ${out}/${recal}.root -partition 1000 -ngc
	#CNV calling
	cnvnator -genome $index -root ${out}/${recal}.root -call 1000	-ngc > ${out}/${recal}.CNV.txt
	#Convert txt file to vcf file
	cnvnator2VCF.pl ${out}/${recal}.CNV.txt > ${out}/${recal}.CNV.vcf
	
	echo "Copy number variation calling done!!!!"
	
}
i=`pwd`
o=`pwd`
r=hg19

while getopts :i:o:b:r: option
do
	case $option in
		i) i=$OPTARG;;
		o) o=$OPTARG;;
		b) b=$OPTARG;;
		r) r=$OPTARG;;
		*);;
	esac
done

if [ "$i" != "" -a "$o" != "" -a "$r" != "" -a "$b" != "" ]
then
	main $b $i $o $r
else
	 usage
fi


