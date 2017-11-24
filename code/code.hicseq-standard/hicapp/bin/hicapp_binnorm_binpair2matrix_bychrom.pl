#!/usr/bin/perl
##
# Generate interaction matrix (by chromosome) from pair file
##

my $pair = $ARGV[0];
my $res = $ARGV[2];
my %itr;
#my @chrom;
my %chrom;
my $matrix = $pair;
my $matrixf = $pair;
$matrix =~ s/pair/matrix/g;
$matrixf =~ s/pair/bychrom.matrix/g;
#print "$matrix\n";
#
if (! -d "$matrixf"){
	mkdir "$matrixf", 0755;
}

## read in 2bed files, and keep intra-chrom pairs
open IN  ,  $pair or die $!;
while (<IN>){
	chomp;
	my @lines = split /\t/,$_;
	my $rd1 = $lines[0];
	my $rd2 = $lines[1];
	my $chr1 = (split /\_/,$rd1)[0];
	my $chr2 = (split /\_/,$rd2)[0];
	if ($chr1 eq $chr2){ #from same chromosome
		$itr{$rd1}{$rd2} = $lines[2];
		$itr{$rd2}{$rd1} = $lines[2];
	}
}
close IN;

## make windows
chomp (my @mat = `bedtools makewindows -g $ARGV[1] -w $res |awk '{print \$1"_"\$2}'`);
#print "pos\t".join("\t", @mat)."\n";

## find chrom
for my $i (0..$#mat){
	my @tag = split /\_/, $mat[$i];
	#push (@chrom, $tag[0]);
	push (@{$chrom{$tag[0]}}, $i);
}

## open all files and add header
foreach my $j (keys %chrom){
	my $fh = "$j.fh";
	open $fh , ">$matrixf/$matrix.$j.$j.mat" or die $!;
	print $fh "pos";
	foreach my $k (@{$chrom{$j}}){
		print $fh "\t$mat[$k]";
	}
	print $fh "\n";

	## write matrix
	foreach my $k (@{$chrom{$j}}){
		print $fh "$mat[$k]";
		foreach my $m (@{$chrom{$j}}){
			if (exists $itr{$mat[$k]}{$mat[$m]}){
				print $fh "\t$itr{$mat[$k]}{$mat[$m]}";
			}else{
				print $fh "\t0";
			}
		}
		print $fh "\n";
	}
	## close file
	close($fh)
}
