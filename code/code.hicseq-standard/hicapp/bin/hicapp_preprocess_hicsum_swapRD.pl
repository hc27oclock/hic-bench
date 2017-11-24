#!/usr/bin/perl
##
# Swap RD1 and RD2 based on genomic coordinates from hicsum file to make RD1 < RD2
##

## read in hicsum files
open IN  ,  $ARGV[0] or die $!;
while (<IN>){
	chomp;
	my @lines = split /\t/,$_;
	
	my $rd1 = "$lines[0]\t$lines[1]\t$lines[2]\t$lines[3]";
	my $rd2 = "$lines[4]\t$lines[5]\t$lines[6]\t$lines[7]";
	
	my $chr1 = $lines[0];
	$chr1 =~ s/chr//g;
	$chr1 =~ s/X/97/g;
	$chr1 =~ s/Y/98/g;
	$chr1 =~ s/M/99/g;
	my $chr2 = $lines[4];
	$chr2 =~ s/chr//g;
	$chr2 =~ s/X/97/g;
	$chr2 =~ s/Y/98/g;
	$chr2 =~ s/M/99/g;
	
	my $pos1 = $lines[1];
	my $pos2 = $lines[5];
	
	if ($chr1 == $chr2){
		if ($pos1 <= $pos2){
			print "$rd1\t$rd2\n";
		}else{
			print "$rd2\t$rd1\n";
		}
	}elsif ($chr1 < $chr2){
		print "$rd1\t$rd2\n";
	}else{
		print "$rd2\t$rd1\n";
	}
	
}
close IN;
