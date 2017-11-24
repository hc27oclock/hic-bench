#!/usr/bin/perl
##
# get all possible contact pairs within a certain genomic distance
##

my $pair_file = $ARGV[0];
my $icb_file = $ARGV[1];
my $ld = $ARGV[2]; # lower distance boundry: 1
my $hd = $ARGV[3]; # higher distance boundry: 2000000

my %icb;
my %read;
my @pos;
my %cnt;
my %chrom;

## load stats file
open IN , $icb_file or die $!;
while (<IN>){
	chomp;
	my @lines = split /\t/,$_;
  next if ($lines[1] == 0);
	$icb{$lines[0]}=$lines[1];
	push(@pos, $lines[0]);
}
close IN;

## load pair file
open IN , $pair_file or die $!;
while (<IN>){
	chomp;
	my @lines = split /\t/,$_;
	$cnt{$lines[0]}+=$lines[2];
	$cnt{$lines[1]}+=$lines[2];
	$read{$lines[0]}{$lines[1]}=$lines[2];
}
close IN;


## sep chrom based on pos
for my $i (0..$#pos){
	my @tag = split /\_/, $pos[$i];
	push (@{$chrom{$tag[0]}}, $i);
}

## select position
foreach my $k (keys %chrom){
	for my $ii (0..$#{$chrom{$k}}){
		my $i = $pos[$chrom{$k}[$ii]];
		my @rd1 = split /\_/,$i;
		if (exists $cnt{$i}){ #genomic position in pair file
			for my $jj ($ii..$#{$chrom{$k}}){
				my $j = $pos[$chrom{$k}[$jj]];
				my @rd2 = split /\_/,$j;
				my $dist = abs($rd2[1]-$rd1[1]);
				if ($dist >= $ld and $dist <= $hd){ #dist in range
					if (exists $cnt{$j}){ #genomic position in pair file
						my $icb_2 = $icb{$i}*$icb{$j};
						if (exists $read{$i}{$j}){
							print "$i\t$j\t$read{$i}{$j}\t$dist\t$icb_2\n";
						}elsif (exists $read{$j}{$i}){
							print "$j\t$i\t$read{$j}{$i}\t$dist\t$icb_2\n";
						}else{
							print "$i\t$j\t0\t$dist\t$icb_2\n";
						}
					}
				}
			}
		}
	}
}
