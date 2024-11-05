#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use Bio::Cigar;  # CPAN this if you want to run this script
use List::Util qw(max);

use vars qw/ %opt /;

# Initialization and usage instructions
sub init() {
	getopts("hv:s:i:g:", \%opt) or usage();
	usage() if $opt{h};
	die "Require path specification for VCF file: -v\n" if not defined $opt{v};
	die "Require path specification for SAM file: -s\n" if not defined $opt{s};
	die "Require individual ID: -i\n" if not defined $opt{i};
	die "Require group ID: -g \n" if not defined $opt{g};
}

sub usage() {
	print STDERR << "EOF";
    This program gathers variant haplotype sites from a SAM alignment file and reports a summary file for those variant sites

    usage: $0 [-h] -v vcf_file -s sam_file -i int -g str
EOF
	exit;
}

init();

#--- Read VCF file and store variants ---
my $vcf;
my $hap;

open VCF, $opt{v};
while (<VCF>) {
	next if /^#/;
	my @line = split "\t";
	next if length($line[3]) > 1;  # Skip indels

	my @snp = split ",", $line[4];
	next if length($snp[0]) > 1;  # Skip complex events

	push @{$vcf->{$line[0]}}, $line[1];
	push @{$vcf->{"ref_".$line[0]}}, $line[3];
	push @{$vcf->{"der_".$line[0]}}, $snp[0];
}
close VCF;

#--- Process SAM file ---
open SAM, $opt{s};
while (<SAM>) {
	next if /^\@/;
	my @lines = split "\t";
	my $id = $lines[2];
	next if $lines[1] >= 256;  # Skip secondary/multi-mapping
	my $st_qpos = $lines[3];
	next if not defined $vcf->{$id};

	next if $lines[5] eq "*";
	my $cigar = Bio::Cigar->new($lines[5]);
	my @qseq = split "", $lines[9];
	my @qseq_qual = split "", $lines[10];
	next if $#qseq < 1;

	my $hapRead = {};
	$hapRead->{"seq"} = "";
	my $has_missing_data = 0;  # Track if there's missing data

	for my $rpos (@{$vcf->{$id}}) {
		my $rpos_adj = $rpos - $st_qpos + 1;

		if ($cigar->reference_length < $rpos_adj || $rpos_adj < 1) {
			$hapRead->{"seq"} .= "N";  # unknown base
			$has_missing_data = 1;
			next;
		}
		my ($qpos, $op) = $cigar->rpos_to_qpos($rpos_adj);

		if (not defined $qpos) {
			$hapRead->{"seq"} .= "X";  # deletion site
			$has_missing_data = 1;
		} else {
			$hapRead->{"seq"} .= $qseq[$qpos - 1];
		}
	}

	# Only proceed if there’s no missing data in hapRead->{"seq"}
	next if $has_missing_data;

	# Populate hap structure with valid haplotypes (no missing data)
	$hap->{$id}->{$hapRead->{"seq"}}->{"ct"}++;
	for my $i (0 .. $#{$vcf->{$id}}) {
		if (!defined ${$hap->{$id}->{$hapRead->{"seq"}}->{"maxC"}}[$i]) {
			${$hap->{$id}->{$hapRead->{"seq"}}->{"maxC"}}[$i] = 0;
			${$hap->{$id}->{$hapRead->{"seq"}}->{"sC"}}[$i] = 0;
		}

		my $q = 10**(-(ord(${$hapRead->{"qual"}}[$i]) - 33) / 10);
		${$hap->{$id}->{$hapRead->{"seq"}}->{"sC"}}[$i] += 1 - $q;
		${$hap->{$id}->{$hapRead->{"seq"}}->{"maxC"}}[$i] = max(1 - $q, ${$hap->{$id}->{$hapRead->{"seq"}}->{"maxC"}}[$i]);
	}
}

#--- Output a haplotype summary file ---
for my $id (keys %{$hap}) {
	for my $h (keys %{$hap->{$id}}) {
		print join "\t", $opt{g}, 
			$opt{i}, 
			$id, 
			$h, 
			$hap->{$id}->{$h}->{"ct"},
			(join ",", @{$hap->{$id}->{$h}->{"sC"}}),
			(join ",", @{$hap->{$id}->{$h}->{"maxC"}})."\n";
	}
}
