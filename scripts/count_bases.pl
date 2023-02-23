#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

# script count_bases.pl
# count base frequency at each position of a list of sequences of same length
# Stéphane Plaisance - VIB-Nucleomics Core - 2023-02-10 v1.00
# visit our Git: https://github.com/Nucleomics-VIB

# bioawk -c fastx '{print $comment}' \
#  V350134580_L01_e1_f91_demux/V350134580_L01_435301ABE103GEX_2.fq.gz \
#  | grep -v AAACCATTCTTCCCACCCTC \
#  |  /opt/biotools/SplitBarcode/dev/count_bases.pl

# store in matrix of 5 rows (bases) and 20 columns (the sequence length in the input data)
my @matrix=();
push @matrix, [(0) x 20] for 0 .. 4;

my @bases=split //,"ATGCN";
my %bases = map { $bases[$_] => $_ } sort 0..$#bases;

# debug
#  while (my ($key, $value) = each (%bases))
#  {
#   $value = $bases{$key};
#   print "  $key <-> $value\n";
#  }

# store max sequence length for final reporting
my $seqlen=0;

foreach my $line ( <STDIN> ) {
    chomp( $line );

	my @nuc=split //, $line;
	$seqlen=( $#nuc > $seqlen ? $#nuc : $seqlen );

	for my $i (0 .. $#nuc) {
	  my $nuc=$nuc[$i];
	  # add 1 to counter for this matrix position
	  $matrix[$bases{$nuc}][$i]+=1;
	}
}

print "max sequence size was ".($seqlen+1)."\n";

print "#\t".(join("\t", 1..20))."\n";
print ($bases[$_].": \t".(join("\t", @{ $matrix[$_] })."\n")) for 0..$#matrix;

# also save to file for gnuplot
my $filename="bases.dat";
open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
print $fh "#\t".(join("\t", 1..20))."\n";
print $fh ($bases[$_].": \t".(join("\t", @{ $matrix[$_] })."\n")) for 0..$#matrix;
close $fh;

# transpose data for gnuplot
system("transpose -t bases.dat > tbases.dat");
# system("gnuplot plotbases.gp");

exit 0;

# bioawk -c fastx '{print $comment}' V350134580_L01_e1_f91_demux/V350134580_L01_435301ABE103GEX_2.fq.gz | grep -v AAACCATTCTTCCCACCCTC |  /opt/biotools/SplitBarcode/dev/count_bases.pl
# max size was 19
#         1       2       3       4       5       6       7       8       9       10      11      12      13      14      15      16      17      18      19      20
# A:      128830  129403  129401  5335    0       122476  8532    0       23      1401    6632    962     522     3758    126025  5378    2466    4306    11433   2048
# T:      574     2       2       0       0       1       120871  129387  9816    94713   117545  1902    2912    2805    2712    3528    2593    3377    115927  4012
# G:      1       0       1       1       1       1       2       0       20      21084   1715    2248    923     2522    545     2831    975     2241    325     14150
# C:      0       0       1       124069  129404  6927    0       18      119545  12204   3507    124293  125047  120320  123     117668  123286  119481  1719    109195
# N:      0       0       0       0       0       0       0       0       1       3       6       0       1       0       0       0       85      0       1       0
# 
# # https://www.xmodulo.com/draw-stacked-histogram-gnuplot.html
# # https://stackoverflow.com/questions/59111654/gnuplot-how-to-plot-bash-array-without-dumping-it-to-a-file
# 
# 
# gnuplot -p -e "set terminal jpeg giant font \"Helvetica\" 16
# 
# set output 'base-substitutions.jpg'
# set key left
# 
# set grid y
# set style data histograms
# set style histogram rowstacked
# set boxwidth 0.5
# set style fill solid 1.0 border -1
# set ytics 10 nomirror
# set yrange [:60]
# set ylabel \"Bases at this position\"
# set ytics 10
# 
# plot \"/dev/stdin\"" < <(
# print "#\t".(join("\t", 1..20))."\n";
# print ($bases[$_].": \t".(join("\t", @{ $matrix[$_] })."\n")) for 0..$#matrix;
# )
# 
# 
