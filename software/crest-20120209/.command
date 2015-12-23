ls
l
export PATH=/home/liucj/tools/crest-20120209/CAP3/:/home/liucj/tools/crest-20120209/blatSuite-34/:$PATH
l
echo $HOStNAME
ps xf
kill 78250
l
cd blatSuite-34/
l
chmod -R 777 blatSuite-34/ CAP3/
l
l
cd
ln -s /home/liucj/liucj/REFDATA/genome_build/hg19-bwa-index/hg19.fasta .
l
gfServer start localhost 6666 hg19.fasta.2bit
l
cd blatSuite-34/
l
faToTwoBit hg19.fasta hg19.2bit
gfServer start localhost 6666 hg19.fasta.2bit
gfServer start localhost 6666 hg19.2bit
ps xf
gfServer start localhost 6666 hg19.2bit &
ps xf
free
top 
l
extractSClip.pl -i tumor.bam --ref_genome hg19.fasta 
l
CREST.pl -f tumor.bam.cover -d tumor.bam --ref_genome hg19.fasta --2bitdir $PWD -t $PWD/hg19.2bit --blatserver localhost --blatport 6666
l
cd example/
l
ps xf
CREST.pl -f tumor.bam.cover -d tumor.bam --ref_genome hg19.fasta --2bitdir $PWD -t $PWD/hg19.2bit --blatserver localhost --blatport 6666
l
samtool view tumor.bam
samtools view tumor.bam
l
CREST.pl -f tumor.bam.cover -d tumor.bam --ref_genome hg19.fasta --2bitdir $PWD -t $PWD/hg19.2bit --blatserver localhost --blatport 6666
l
ll
les tumor.bam.predSV.txt
k
CREST.pl -f tumor.bam.cover -d tumor.bam --ref_genome hg19.fasta --2bitdir $PWD -t $PWD/hg19.2bit --blatserver localhost --blatport 6666
export PATH=/home/liucj/tools/samtools-0.1.19/:$PATH
CREST.pl -f tumor.bam.cover -d tumor.bam --ref_genome hg19.fasta --2bitdir $PWD -t $PWD/hg19.2bit --blatserver localhost --blatport 6666
l
satmools view tumor.bam
samtools view tumor.bam
samtools view tumor.bam|les
samtools view tumor.bam|awk '{$2=‘chr' + $2;print $0}'|les
samtools view tumor.bam|awk '{$2="chr" + $2;print $0}'|les
samtools view tumor.bam|awk 'Begin{FS="\t";OFS="\t";}{$2="chr" + $2;print $0}'|les
samtools view tumor.bam|awk 'BEGIN{FS="\t";OFS="\t";}{$2="chr" + $2;print $0}'|les
samtools view tumor.bam|awk 'BEGIN{FS="\t";OFS="\t";}{$3="chr"+$3;print $0}'|les
samtools view tumor.bam|awk 'BEGIN{FS="\t";OFS="\t";}{$3="chr"+$3;print $3}'|les
l
samtools 
l
cd tool-data/
l
man sort
l
cd /home/liucj/piplines/resequencing/genome_resequencing/example/Mapping
ln -s /home/liucj/piplines/resequencing/genome_resequencing/example/Mapping/NA12878.sam.dedup.realign.recal.bam .
ln -s /home/liucj/piplines/resequencing/genome_resequencing/example/Mapping/NA12878.sam.dedup.realign.recal.bai .
l
l
extractSClip.pl -i NA12878.sam.dedup.realign.recal.bam --ref_genome hg19.fasta
l
cd /home/liucj/piplines/resequencing/genome_resequencing/example/Mapping
extractSClip.pl -i NA12878.sam.dedup.realign.recal.bam --ref_genome hg19.fasta
rm NA12878.sam.dedup.realign.recal.bai NA12878.sam.dedup.realign.recal.bam
l
cp /home/liucj/piplines/resequencing/genome_resequencing/example/Mapping/NA12878.sam.dedup.realign.recal.b* .
l
rm NA12878.sam.dedup.realign.recal.bam.grp
extractSClip.pl -i NA12878.sam.dedup.realign.recal.bam --ref_genome hg19.fasta
l
mv NA12878.sam.dedup.realign.recal.bai NA12878.sam.dedup.realign.recal.bam.bai
extractSClip.pl -i NA12878.sam.dedup.realign.recal.bam --ref_genome hg19.fasta
ls
les NA12878.sam.dedup.realign.recal.bam.cover
l
CREST.pl -f NA12878.sam.dedup.realign.recal.bam.cover -d NA12878.sam.dedup.realign.recal.bam --ref_genome hg19.fasta --2bitdir $PWD -t $PWD/hg19.2bit  --blatserver localhost --blatport 666
ps xf
CREST.pl -f NA12878.sam.dedup.realign.recal.bam.cover -d NA12878.sam.dedup.realign.recal.bam --ref_genome hg19.fasta --2bitdir $PWD -t $PWD/hg19.2bit  --blatserver localhost --blatport 666
CREST.pl -f NA12878.sam.dedup.realign.recal.bam.cover -d NA12878.sam.dedup.realign.recal.bam --ref_genome hg19.fasta --2bitdir $PWD -t $PWD/hg19.2bit  --blatserver localhost --blatport 6666
l
les NA12878.sam.dedup.realign.recal.bam.predSV.txt
l
ll
l
bam2html.pl -d NA12878.sam.dedup.realign.recal.bam -i NA12878.sam.dedup.realign.recal.bam.predSV.txt --ref_genome hg19.fasta -o NA12878.sam.dedup.realign.recal.bam.predSV.html
l
rm 4lvkMi7Abr aguMOmuca1 bi972PmDYC H4ts6BnHSq w7PeYAKGxf Y5oU7F6bL0
ls
l
gfServer stop localhost 6666
rm NA12878.sam.dedup.realign.recal.bam*
l
ll
l
rm hg19.*
l
cd ./
cd ../
l
chmod -R CAP3
chmod -R 777 CAP3
chmod -R 777 blatSuite-34/
l
cp CAP3/cap3 .
les README 
l
cp CAP3/* .
cp blatSuite-34/* .
l
ps xf
netstat -a
netstat -n 6666
l
ps xf
kill 27247 27826 28689 28690
ps xf
kill 26289
ll /home/liucj/piplines/resequencing/genome_resequencing/example/Mapping
rm te
l
ll SCValidator.pm
l
les README 
l
les README 
l
rm README 
l
ls
top -u liucj
l
man test
ls
l
l
cd ../
l
cd ../../
