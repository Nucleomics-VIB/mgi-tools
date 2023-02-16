#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use Data::Dumper;
use FileHandle;
# add edit distance function
use Text::Levenshtein::Flexible qw( levenshtein levenshtein_l_all );

# script: SplitDualBarcodes_NC2.pl
# rewritten with edit distance for correction
# Stéphane Plaisance - VIB-Nucleomics Core - 2023-02-14 v1.04
# visit our Git: https://github.com/Nucleomics-VIB
#
# this script runs a modified version of the original 'SplitDualBarcodes.pl'
# (https://github.com/gateswell/SplitBarcode)
# this modified version can be found in our repo 
# (https://github.com/Nucleomics-VIB/mgi-tools)
#
# + small edits:
# -+ replace gzip post-compression through shell scripts by gzip in-pipe compression
# -+ changed headers and file names to match current on-device output
# Text::Levenshtein to correct barcodes

my $version="1.04, 2023_02_16";

my $usage=<<USAGE;
	Usage:
		perl $0 [options]
			*-r1 --read1 <string>		read1.fq.gz
			*-r2 --read2 <string>		read2.fq.gz
			 -e  --errNum <int>		mismatch number for each barcode/index [default: 1]
			*-f  --firstCycle <int>		first cycle of barcode
			*-b  --barcodeList <string>	barcodes list (tab or space-delimited)
			 -rc --revcom	<Y|N>		generate reverse complement of barcode.list [default: N]
			 -c  --compress <Y|N>		compress (.gz) output or not [default: Y]
			 -o  --outdir <string>		output directory [default: ./]
			 -h  --help			print help information and exit
	Example:
		perl $0 -r1 read1.fq.gz -r2 read2.fq.gz -e 2 -f 101 -b barcode.list -o /path/outdir
		
	=======example barcode.list================
	#barcodeName	barcodeSeq1	barcodeSeq2
	barcode_1	ATGCATCTAA	TATAGCCTAG
	barcode_2	AGCTCTGGAC	CTCTATCGTC
	===========================================
	
	(adapted from the original gateswell's SplitDualBarcodes.pl)
	VIB-NC version $version
USAGE

#=============global variants=============
my ($read1,$read2,$errNum,$fc,$bl,$compress,$outdir,$rc,$help);
my (%bchash,$prefix,$ambo1,$ambo2);
#=========================================
GetOptions(
	"read1|r1=s"=>\$read1,
	"read2|r2=s"=>\$read2,
	"errNum|e:i"=>\$errNum,
	"firstCycle|f=i"=>\$fc,
	"barcodeList|b=s"=>\$bl,
	"revcom|rc:s"=>\$rc,
	"compress|c:s"=>\$compress,
	"outdir|o:s"=>\$outdir,
	"help|h:s"=>\$help
);

# default errNum if not provided
#$errNum ||= 1;
if(defined $errNum && $errNum == 0){
	$errNum=0;
}
elsif(!defined $errNum ){
	$errNum ||= 1;
}

$outdir ||= `pwd`;
$compress ||= 'Y';
$rc ||= 'N';

# required arguments or die
if(!$read1 || !$read2 || !$fc || !$bl || $help ){
	die "$usage";
}

#========global variables==========
my (%barhash,%oh,%oribar,%oriname,%correctBar,%correctedBar,%unknownBar,$totalReadsNum);
my (%tagNum,$am1,$am2,@fq,$barcode_len,$barcode_len1,$barcode_len2,%ori_tag,%mis_barNum);
#=========================

my $name=basename($read2);
# the following regexp is matching expected MGI folder naming
# eg V300009631_128A_L01_read_2.fq.gz => V300009631_128A_L01
# prefix is all text before _L0X_read_2 in read2 name
$prefix=$1 if $name=~/(.*)\_(\w+)_2\.fq(.gz)?/;
unless(-d $outdir){
	print STDERR "creating $outdir folder\n";
	`mkdir -p $outdir`;
}
$outdir=abs_path($outdir);
chomp($outdir);

###################################
# create files for undecoded reads
###################################

# create undecoded fastq output file handles
if(uc($compress) eq 'Y'){
	open $am1,"|gzip -9 >$outdir/$prefix\_undecoded_1.fq.gz" or die $!;
	open $am2,"|gzip -9 >$outdir/$prefix\_undecoded_2.fq.gz" or die $!;
}else{
    open $am1,">$outdir/$prefix\_undecoded_1.fq" or die $!;
    open $am2,">$outdir/$prefix\_undecoded_2.fq" or die $!;
}
# add both undecoded filenames to fastq-name array
push @fq,"$outdir/$prefix\_undecoded_1.fq";
push @fq,"$outdir/$prefix\_undecoded_2.fq";

########################################
# load barcode file and store in arrays
########################################

# barcode file handle
open my $fh,$bl or die "$bl barcode list not found !\n$!";

# store barcodes and sample names for use by tryCorrect function
my (@name,@correct,@bc1,@bc2);

# loop through all barcode pairs from provided list
while(<$fh>){	#1	ATGCATCTAA	TATAGCCTAG
	next if /^#/;
	chomp;
	my @tmp=split /\s+/,$_;
	my $barcode_seq;
    # if option rc was set to 'Y'
	if(uc($rc) eq 'Y'){
        # reverse the barcodes
		my $t=reverse(uc($tmp[1]));		#AATCTACGTA
		$tmp[1]=reverse(uc($tmp[2]));	#GATCCGATAT
		$tmp[2]=$t;
		$barcode_seq=$tmp[1].$tmp[2];	#GATCCGATATAATCTACGTA
        # complement the merged sequence
		$barcode_seq=~tr/ATCGN/TAGCN/;	#CTAGGCTATATTAGATGCAT
	}else{
	    # barcodes are read as-is
		$tmp[1]=uc($tmp[1]);			#ATGCATCTAA
		$tmp[2]=uc($tmp[2]);			#TATAGCCTAG
		$barcode_seq=$tmp[1].$tmp[2];	#ATGCATCTAATATAGCCTAG
	}
	
	# add $barcode_name to array of original barcodes
	push(@name,$tmp[0]);
	
	# add barcode name to hash with barcode_seq as key
	$oriname{$barcode_seq} = $tmp[0];	
	
	# add both $barcode sequences to array of original barcodes
	push(@bc1,$tmp[1]);
	push(@bc2,$tmp[2]);

	# add merged to correct barcodes
	push(@correct,$barcode_seq);

	# add merged barcode to hash of original barcodes
	$oribar{$barcode_seq} = 1;
	
	# create demultiplexed fastq output file handles for both reads
	if(uc($compress) eq 'Y'){
 	    open $oh{$tmp[0]}[0],"|gzip -9 >$outdir/$prefix\_$tmp[0]\_1.fq.gz" or die $!;
	    open $oh{$tmp[0]}[1],"|gzip -9 >$outdir/$prefix\_$tmp[0]\_2.fq.gz" or die $!;
    }else{
        open $oh{$tmp[0]}[0],">$outdir/$prefix\_$tmp[0]\_1.fq" or die $!;
	    open $oh{$tmp[0]}[1],">$outdir/$prefix\_$tmp[0]\_2.fq" or die $!;
    }
    # add both demultiplexed read filenames to fastq-name array @fq
	push @fq,"$outdir/$prefix\_$tmp[0]\_1.fq";
	push @fq,"$outdir/$prefix\_$tmp[0]\_2.fq";
}
close $fh;

# measure both barcode lengths from first barcode pair
$barcode_len1=length($bc1[0]);
$barcode_len2=length($bc2[0]);
$barcode_len=$barcode_len1+$barcode_len2;

###############################################
# read from paired raw fastq file and classify
###############################################

my($rd1,$rd2);

# check if read2 is fastq or fastq.gz
if($read2=~/fq$/){
	open $rd1,$read1 or die $!;
	open $rd2,$read2 or die $!;
}elsif($read2=~/fq.gz$/){
	open $rd1,"zcat $read1|" or die $!;
    open $rd2,"zcat $read2|" or die $!;
}else {
	die "read2 file does not have the expected extension (.fq or .fq.gz)!\n";
}

while(<$rd1>){

    # read 4 lines from both reads and store them to fastq fields
	my $head1= $_;
	my $seq1 = <$rd1>;
	my $plus1= <$rd1>;
	my $qual1= <$rd1>;
	my $head2= <$rd2>;
	my $seq2 = <$rd2>;
	my $plus2= <$rd2>;
	my $qual2= <$rd2>;
	
	# chomp all strings
	chomp($head1,$seq1,$plus1,$qual1,$head2,$seq2,$plus2,$qual2);

	# count read pair
	$totalReadsNum ++;
	#$totalReadsNum =~ m/00000$/ and print STDERR "$totalReadsNum\n";
    #$totalReadsNum =~ m/00000$/ and last;
	
	# extract double-barcode sequence from read2 at expected starting position $fc (-f)
	my $barseq=substr($seq2,$fc-1,$barcode_len+1);
	#print STDOUT "# $barseq\n";
	
	# add 1 to this barcode counter for SequenceStat.txt
	# note: can present missmatches with the expected barcode
	$tagNum{$barseq} ++;

    # clip barcode sequence from read2 seq and qual lines
	my $clipseq2=substr($seq2,0,$fc-1); #.substr($seq2,$fc+$barcode_len-1,);
	my $clipqual2=substr($qual2,0,$fc-1); #.substr($qual2,$fc+$barcode_len-1,);	
	
	# add merged barcode sequence to read1 name
	$head1 = $head1." ".$barseq;
	$head2 = $head2." ".$barseq;
	
	#######################################
	# barcode is found in the allowed list
	#######################################
	
	if(exists $oribar{$barseq}){
		# add 1 to counter
    	$correctBar{$barseq} += 1;
   	
    	#print STDOUT "# $barseq is found in list with name $oriname{$barseq}\n";
       	
   		# write reads to files
    	$oh{$oriname{$barseq}}[0]->print("$head1\n$seq1\n$plus1\n$qual1\n");
		$oh{$oriname{$barseq}}[1]->print("$head2\n$clipseq2\n$plus2\n$clipqual2\n");
		# process next read pair
		next;
		}

	#########################################################
	# try find a correction with at most $errNum per barcode
	#########################################################
	
	# extract both barcodes for correction testing
	my $test1 = substr($seq2,$fc-1,$barcode_len1+1);
	my $test2 = substr($seq2,$fc+$barcode_len1-1,$barcode_len2+1);
	
   	my $cor = tryCorrect2($test1,$test2);
   	if ($cor ne "na"){
		# add 1 to counter
   		$correctedBar{$cor} += 1;

   		#print STDOUT "# $name[$cor] $bc1[$cor] $bc2[$cor]\n";

		# write reads to files
    	$oh{$oriname{$cor}}[0]->print("$head1\n$seq1\n$plus1\n$qual1\n");
		$oh{$oriname{$cor}}[1]->print("$head2\n$clipseq2\n$plus2\n$clipqual2\n");
   	} else {
		# add 1 to counter
		$unknownBar{$barseq} += 1;

    	#print STDOUT "# barcode is unexpected, store read pair as undecoded\n";
    	
		# write reads to files
		print $am1 "$head1\n$seq1\n$plus1\n$qual1\n";
		print $am2 "$head2\n$clipseq2\n$plus2\n$clipqual2\n";
	   	}
 	}

# close file handles
close $rd1;close $rd2;
close $am1;close $am2;

###############################
# print summary to BarcodeStat
###############################

open my $BS,">$outdir/BarcodeStat.txt" or die $!;
print $BS "#Barcode\tCorrect\tCorrected\tTotal\tPercentage(%)\n";

my($totalcorrect,$totalcorrected,$totalbarreads,$totalpct);
for my $seq(sort {$oriname{$a} cmp $oriname{$b}} keys %oribar){
	# unlikely but when $correctBar{$barhash{$seq}} is unset => set it to 0
    if (!$correctBar{$seq}) {
        $correctBar{$seq}=0;
    }
    # when $correctedBar{$barhash{$seq}} is unset => set it to 0
    if (!$correctedBar{$seq}) {
        $correctedBar{$seq}=0;
    }
    my $BartotalReads = $correctBar{$seq}+$correctedBar{$seq};
	my $pct = ($BartotalReads/$totalReadsNum)*100;
	$totalcorrect += $correctBar{$seq};
	$totalcorrected += $correctedBar{$seq};
	$totalbarreads += $BartotalReads;
	$totalpct += $pct;
	printf $BS "%s\t%d\t%d\t%d\t%.4f%%\n",$oriname{$seq},$correctBar{$seq},$correctedBar{$seq},$BartotalReads,$pct;
}
printf $BS "Total\t%d\t%d\t%d\t%.4f%%\n",$totalcorrect,$totalcorrected,$totalbarreads,$totalpct;
close $BS;

###############################
# print summary to SequenceStat
###############################

open my $SS,">$outdir/SequenceStat.txt" or die $!;
print $SS "#Sequence\tBarcode\tCount\tPercentage(%)\n";

for my $seq(sort {$tagNum{$b}<=>$tagNum{$a}} keys %tagNum){
	my $pct=($tagNum{$seq}/$totalReadsNum)*100;
	if(exists $oriname{$seq}){
		printf $SS "%s\t%s\t%d\t%.2f%%\n",$seq,$oriname{$seq},$tagNum{$seq},$pct;
	}
	else{
		printf $SS "%s\tundecoded\t%d\t%.2f%%\n",$seq,$tagNum{$seq},$pct;
	}
}
close $SS;

#================================
sub search_idx {
my ($q, @arr) = ( @_ );
#my @matches = grep { /^$q$/ } @arr; # returns elements
my @matches = grep { $arr[$_] =~ /^$q$/ } 0..$#arr;  # returns indices
return @matches;
}

#================================
sub tryCorrect2 {
  my ($test1, $test2) = @_;
  my (@res1,@res2,$merge,$match);
  
  @res1=( map { $_->[0] }
  sort { $a->[1] <=> $b->[1] }
  levenshtein_l_all($errNum, $test1, @bc1) );
     
  # barcode#1 matches or could be corrected
  if ( @res1 ){
    # identify the indices of the first barcode match(es)
    my @idx1 = map {search_idx( $_, @bc1 )} @res1;
    
    # only look at the paired barcode#2 in second list
    @res2=( map { $_->[0] }
    sort { $a->[1] <=> $b->[1] }
    levenshtein_l_all($errNum, $test2, @bc2[@idx1]) );

    # barcode#2 matches or could be corrected
    if ( @res2 ){
      # corrected pair exists (not barcode hopping)
      $merge="@res1"."@res2";
      if ( $oribar{$merge} ){
        return $merge;
      }else{
      return "na";
      }
    } else {
      return "na";
    }
  } else {
    return "na";
  }
}
