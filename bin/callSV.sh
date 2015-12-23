#########################################################################
# File Name: callSV.sh
# Author: C.J. Liu
# Mail: samliu@hust.edu.cn
# Created Time: Mon 16 Nov 2015 10:52:07 PM CST
#########################################################################
#!/bin/bash






function usage(){
	echo "bash $0 -b <recal_bam> -i <indir> -o <output dir> -r <ref>"
	echo "Description:"
	echo "	Call strutural variation from recalibration bam file"
	echo "	Use Crest-2012"
	echo "	b for bam argument is required, r is optional default value is hg19"
	echo "	Before running script, you have to be sure gfServer of genome build is on"
	echo "	if not run 'cd $root/data/hg19/genomeBuild/ && $root/software/crest-20120209/gfServer start localhost 6666 hg19.fasta.2bit &' until get Ready"
	return 1
}

function main(){
	recal=$1
	indir=$2
	out=$3
	ref=$4
	test -d $out/CallSV || mkdir $out/CallSV
	out=$out/CallSV
	
	root=$(cd "$(dirname "$0")";cd ../;pwd)
	data=${root}/data
	software=${root}/software
	
	
	soft=${software}/crest-20120209/
	indexdir=${data}/$ref/genomeBuild
	index=${indexdir}/${ref}.fasta
	test -e "$index" || { echo "File $index doesn't exist"; exit 1; }
	#faToTwoBit to end with fa.2bit file. it should be load in the gfServer
	#Go to the path, run "gfServer start localhost 6666 hg19.fasta.2bit & "until you got Query Ready
	index2bit=${index}.2bit
	
	test -e "$index2bit" || { echo "File $index2bit doesn't exist"; exit 1; }

	#Crest.pl use perl module in its own path, be sure the module can be load in the crest
	export PERLLIB=${soft}:$PERLLIB
	#somo of software such as cap3 and blat will be call be crest.pl,but they are not installed in my $HOME path.
	#I put them together in crest path, and export the path so that running crest can find them from export path
	export PATH=${soft}:$PATH
	
	#You will get error "Sorry, the BLAT/iPCR server seems to be down..." if gfServer is off
	#BLAT runs in a stand alone or a Client/Server mode, CREST use Client/Server mode
	#So gfServer must be on, and when using gfClient, you need spycify the host and port the gfServer used
	
	test -L "$out/$recal" || { ln -s $indir/$recal $out/$recal; }
	test -L "$out/${recal}.bai" || ln -s  $indir/${recal%.bam}.bai $out/${recal}.bai
	extractSClip.pl -i $out/$recal --ref_genome $index -o $out
	test -f "${out}/${recal}.cover" || { echo "File ${out}/${recal}.cover doesn't exist" ; exit 1; }
	test -f "${out}/$recal" || { echo "${out}/$recal doesn't exist"; exit 1; }
	
	CREST.pl -f "${out}/${recal}.cover" -d "${out}/$recal" --ref_genome "$index" --2bitdir "$indexdir" -t "$index2bit" --blatserver localhost --blatport 6666 -o $out 
	test -f "${out}/${recal}.predSV.txt" || { echo "${out}/${recal}.predSV.txt doesn't exist"; exit 1; }
	bam2html.pl -d "${out}/${recal}" -i "${out}/${recal}.predSV.txt" --ref_genome "$index" -o "${out}/${recal}.predSV.html" 
	
	echo "************SV calling done***********"
	

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

