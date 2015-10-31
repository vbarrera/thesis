#! /usr/bin/perl

# @author: Victor Barrera
# Description: This script obtains and prints the nucleotide sequence
# for a list of elements in the form:
# id chr startPosition endPosition


use strict;
use warnings;
use Bio::DB::Fasta;
use Getopt::Long;
use Pod::Usage;


## Declaring variables

my$positionsFile=$ARGV[0];
my$indexpath=$ARGV[1];

my $help=0;
my $man=0;
 
GetOptions ('help|?'=>\$help, 'man'=> \$man) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

## Opening file

unless (open(POSFILE,$positionsFile)){
	print "Error opening file \n";
	pod2usage(1);
	die;
}

## Connection to the index file
my $db=Bio::DB::Fasta->new($indexpath);

## For each element, obtain and print the sequence
while (<POSFILE>){
	chomp $_;
	my@elementsArray=split("\t",$_);
	my$id=$elementsArray[0];
	my$chr=$elementsArray[1];
	if($chr=~/_/){
	next;
	}
	my$startcoord=$elementsArray[2];
	my$endcoord=$elementsArray[3];

	print $id."\t".$chr."\t".$startcoord."\t".$endcoord."\t";
	print $db->seq($chr,$startcoord+1,$endcoord)."\n";
}

exit;

__END__

=head1 NAME

seqgetter

=head1 SYNOPSIS

	seqgetter [Elements FILE] [GENOME INDEX FILE]
	
=head1 DESCRIPTION

This script obtains and prints the nucleotide sequence
for a list of elements in the form:
id chr startPosition endPosition


=cut	


