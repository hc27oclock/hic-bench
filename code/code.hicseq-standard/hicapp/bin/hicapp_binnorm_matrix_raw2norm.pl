#!/usr/bin/perl
####
# get normalized matrix from raw matrix and bias vector
####

use strict;
my %bias;
my @bias;

open IN , $ARGV[1] or die $!;
while (<IN>){
	chomp;
	my @lines = split /\t/,$_;
	$bias{$lines[0]} = $lines[1];
	push(@bias, $lines[1]);
}
close IN;


open IN , $ARGV[0] or die $!;
while (<IN>){
	chomp;
	my @lines = split /\t/,$_;
	if ($. == 1){
		print "$_\n";
	}else{
		my $tag = shift(@lines);
		for my $i (0..$#lines){
      my $b=$bias[$i]*$bias{$tag};
      if ($b==0){
        $lines[$i]=0;
      }else{
			  $lines[$i]=$lines[$i]/$b;
			  $lines[$i]=int($lines[$i]*1000+0.5)/1000;
		  }
    }
		print "$tag\t".join("\t", @lines)."\n";
	}
}

close IN;
