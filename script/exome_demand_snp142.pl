#!/usr/bin/perl -w
################################################
#File Name: exome_demand.pl
#Author: C.J. Liu
#Mail: samliu@hust.edu.cn
#Created Time: Fri 12 Sep 2014 03:18:20 PM CST
################################################

use strict;
#use warnings;
use Getopt::Long;
use File::Basename;

my %opts;
GetOptions(\%opts, 'maf=s', 'avi=s', 'h');
my $usage =<<USAGE;
Description:
	This script is to produce result file as liujy's request;
	I tackle the maf addition file txt.maf to split information into different result file
	
	The input file is txt.maf, the following five result files are:
	txt.maf.exonic	
	txt.maf.updownstream
	txt.maf.utr
	txt.maf.splicing
	txt.maf.ns 

	the last suffix name is the meaning of content of file.
Usage:
	perl exome.demand.pl [option] [file]
	Options:
	maf	input file of maf - txt.maf
	avi	input file of avinput - vcf.avinput
	h	help
USAGE
if(defined $opts{'h'}){print $usage; exit 1;}
if(!defined $opts{'maf'}){print "maf option must be required"; exit 1;}
if(!defined $opts{'avi'}){print "avi option must be required"; exit 1;}




###
#store quality infromation in hash
###
my %hash;
open(AVI, "<$opts{'avi'}") or die "error: cant open file: $opts{'avi'}\t$!\n";
while(<AVI>){
	chomp;
	my @keyArray = (split/\t/)[0,1,2,3,4,14];
	my $key = join("#", @keyArray[0,1,2,3,4]);
	$hash{$key} = [split(/:/, $keyArray[5])];	
}
close AVI;

###
#split file maf
###
my $exonic = $opts{'maf'}.'.exonic';
open(EXONIC, ">$exonic") or die("ERROR: cant open file:$exonic\t$!\n");
my $stream = $opts{'maf'}.'.updownstream';
open(STREAM, ">$stream") or die("ERROR: cant open file: $stream\t$!\n");
my $utr = $opts{'maf'}.'.utr';
open(UTR, ">$utr") or die("ERROR: cant open file: $utr\t$!\n");
my $splicing = $opts{'maf'}.'.splicing';
open(SPLICING, ">$splicing") or die("ERROR: cant open file: $splicing\t$!\n");
my $ns = $opts{'maf'}.'.ns';
open(NS, ">$ns") or die("ERROR: cant open file: $ns\t$!\n");

my @title = qw/chr start end ref alt id dbsnp.ref.count dbsnp.ref.freq dbsnp.alt.count dbsnp.alt.freq genotype alleledepth depth ens.func ens.gene ens.gene.symbol ens.exonic ens.AAchange sift cosmic65 esp6500si_all/;
my $first = join("\t", @title);

print EXONIC "$first\n";
print STREAM "$first\n";
print UTR "$first\n";
print SPLICING "$first\n";
print NS "$first\n";

###
#print output files 
###
open(IN, "<$opts{'maf'}") or die("ERROR: cant open file: $opts{'maf'}\t$!\n");
while(<IN>){
	chomp;
	my @array = (split/\t/)[0,1,2,3,4,24,25,26,27,28,9,10,6,11,12,21,22,23];
	my $func = $array[10];
	my $exon = $array[13];
	
	my $key = join "#", @array[0,1,2,3,4];
	my @qual;
	if(defined $hash{$key}){
		@qual = @{$hash{$key}};
	}
	splice @array, 10, 0, @qual[0, 1, 2];
	my $record = join("\t", @array);
	
	if($func eq "exonic" or $func eq "exonic,splicing"){
		print EXONIC "$record\n";
	}
	elsif($func eq "downstream" or $func eq "upstream" or $func eq "upstream;downstream"){
		print STREAM "$record\n";
	}
	elsif($func eq "UTR3" or $func eq "UTR5" or $func eq "UTR5;UTR3"){
		print UTR "$record\n";
	}
	elsif($func eq "splicing"){
		print SPLICING "$record\n";
	}

	if($exon eq "nonsynonymous SNV" or $exon eq "stopgain SNV" or $exon eq "stoploss SNV" or $exon eq "unknown" or $exon eq "frameshift deletion" or $exon eq "frameshift insertion"){
		print NS "$record\n";
	}

}
close IN;
close EXONIC;
close STREAM;
close UTR;
close SPLICING;
close NS;




