open SAM, $opt{s};
while(<SAM>) {
	next if /^\@/; # Skip header lines
	my @lines = split "\t";
	my $id = $lines[2];
	next if $lines[1] >= 256; # skip entries that are secondary alignment or to multiple sites
	my $st_qpos = $lines[3]; # starting query position
	#my $mapq = $lines[4]; # mapping quality score
	# skip if the alignment id is not found in the vcf hash ref
	next if not defined $vcf->{$id};

	next if $lines[5] eq "*";
	my $cigar = Bio::Cigar->new($lines[5]);
	my @qseq = split "", $lines[9];
	my @qseq_qual = split "", $lines[10];

	next if $#qseq < 1;

	my $hapRead = {};
	$hapRead->{"seq"} = "";
	my $ct = 0;
	my $has_missing_data = 0;  # Flag to track if there's any missing data or deletion across positions

	for my $rpos (@{$vcf->{$id}}) {
		my $rpos_adj = $rpos - $st_qpos + 1;
		$ct++;

		# Check for out-of-bound positions or missing reference length
		if ($cigar->reference_length < $rpos_adj || $rpos_adj < 1) {
			$has_missing_data = 1;  # Mark as having missing data
			last;  # Stop processing if any position is missing data
		}
		my ($qpos, $op) = $cigar->rpos_to_qpos($rpos_adj);

		# Check for undefined query position (deletion) and mark missing data
		if (not defined $qpos) {
			$has_missing_data = 1;  # Mark as having missing data
			last;  # Stop processing if any position is missing data
		}

		# Append the base if it's available in the sequence
		$hapRead->{"seq"} .= $qseq[$qpos - 1];
		push @{$hapRead->{"qual"}}, $qseq_qual[$qpos - 1];
	}

	# Skip sequences with any missing or deleted data
	next if $has_missing_data;

	# Populate hap structure with valid haplotypes (fully called SNPs)
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
