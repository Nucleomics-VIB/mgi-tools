#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use Data::Dumper;
use FileHandle;

# SP@NC, 2023-02-09 
# + small edits:
# -+ replace gzip post-compression through shell scripts by gzip in-pipe compression
# + handle undef $correctedBar{$barhash{$barseq}} in case of errNum==0

my $usage=<<USAGE;
	Usage:
		perl $0 [options]
			*-r1 --read1 <string>		read1.fq.gz
			*-r2 --read2 <string>		read2.fq.gz
			 -e  --errNum <int>		mismatch number for each barcode/index [default: 1]
			*-f  --firstCycle <int>		first cylce of barcode
			*-b  --barcodeList <string>	barcodes list
			 -rc --revcom	<Y|N>		generate reverse complement of barcode.list or not
			 -c  --compress <Y|N>		compress (.gz) output or not [default: Y]
			 -o  --outdir <string>		output directory [default: ./]
			 -h  --help			print help information and exit
	Example:
		perl $0 -r1 read1.fq.gz -r2 read2.fq.gz -e 2 -f 101 -b barcode.list -o /path/outdir
		
	============barcode.list===========
	#barcodeNum	barcodeSeq1	barcodeSeq2
	1	ATGCATCTAA	TATAGCCTAG
	2	AGCTCTGGAC	CTCTATCGTC
	===================================
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
#$errNum ||= 1;
if(defined $errNum && $errNum == 0){
	$errNum=0;
}
elsif(!defined $errNum ){
	$errNum ||= 1;
}
$outdir ||= `pwd`;
$compress ||= 'Y';
$rc ||= 'Y';

if(!$read1 || !$read2 || !$fc || !$bl || $help ){
	die "$usage";
}

#========global variables==========
my (%barhash,%oh,%oribar,%correctBar,%correctedBar,%unknownBar,$totalReadsNum);
my (%tagNum,$am1,$am2,@fq,$barcode_len,$barcode_len1,$barcode_len2,%ori_tag,%mis_barNum);
#=========================

my $name=basename($read2);
$prefix=$1 if $name=~/(.*)\_(\w+)_2\.fq(.gz)?/;
unless(-d $outdir){
	print STDERR "$outdir: No such directory, but we will create it\n";
	`mkdir -p $outdir`;
}
$outdir=abs_path($outdir);
chomp($outdir);
open my $fh,$bl or die "$bl No such file, check it !\n$!";

# create unbarcoded output files
if(uc($compress) eq 'Y'){
	open $am1,"|gzip -9 >$outdir/$prefix\_ambiguous_1.fq.gz" or die $!;
	open $am2,"|gzip -9 >$outdir/$prefix\_ambiguous_2.fq.gz" or die $!;
}else{
    open $am1,">$outdir/$prefix\_unbarcoded_1.fq" or die $!;
    open $am2,">$outdir/$prefix\_unbarcoded_2.fq" or die $!;
}
push @fq,"$outdir/$prefix\_unbarcoded_1.fq";
push @fq,"$outdir/$prefix\_unbarcoded_2.fq";

# create summary files 
open my $BS,">$outdir/BarcodeStat.txt" or die $!;
open my $SS,">$outdir/TagStat.txt" or die $!;

# add header to summary files
print $BS "#SpeciesNO\tCorrect\tCorrected\tTotal\tPct\n";
print $SS "#Sequence\tSpeciesNO\treadCount\tPct\n";

# loop through all barcode pairs from provided list
while(<$fh>){	#1	ATGCATCTAA	TATAGCCTAG
	next if /^#/;
	chomp;
	my @tmp=split /\s+/,$_;
	my $barcode_seq;
    # if option rc was set to 'Y'
	if(uc($rc) eq 'Y'){
        # reverse the barcodes
		$tmp[1]=reverse(uc($tmp[1]));	#AATCTACGTA
		$tmp[2]=reverse(uc($tmp[2]));	#GATCCGATAT
		$barcode_seq=$tmp[2].$tmp[1];	#GATCCGATATAATCTACGTA
        # complement the merged sequence
		$barcode_seq=~tr/ATCGN/TAGCN/;	#CTAGGCTATATTAGATGCAT
	}else{
		$tmp[1]=uc($tmp[1]);			#ATGCATCTAA
		$tmp[2]=uc($tmp[2]);			#TATAGCCTAG
		$barcode_seq=$tmp[1].$tmp[2];	#ATGCATCTAATATAGCCTAG
	}
	$oribar{$barcode_seq} = 1;
	$barcode_len1=length($tmp[1]);$barcode_len2=length($tmp[2]);
	&bar_hash($barcode_seq,$tmp[0],$errNum,\%barhash);
	
	# create output handles
	if(uc($compress) eq 'Y'){
 	    open $oh{$barhash{$barcode_seq}}[0],"|gzip>$outdir/$prefix\_$tmp[0]\_1.fq.gz" or die $!;
	    open $oh{$barhash{$barcode_seq}}[1],"|gzip>$outdir/$prefix\_$tmp[0]\_2.fq.gz" or die $!;
    }else{
        open $oh{$barhash{$barcode_seq}}[0],">$outdir/$prefix\_$tmp[0]\_1.fq" or die $!;
	    open $oh{$barhash{$barcode_seq}}[1],">$outdir/$prefix\_$tmp[0]\_2.fq" or die $!;
    }	
	push @fq,"$outdir/$prefix\_$tmp[0]\_1.fq";
	push @fq,"$outdir/$prefix\_$tmp[0]\_2.fq";
}
close $fh;

$barcode_len=$barcode_len1+$barcode_len2;
for my $line_seq(keys %barhash){
		$mis_barNum{$line_seq} ++;
		if ($mis_barNum{$line_seq} >=2){
			print STDERR "Set a smaller mismatch value because of same barcodes\($barhash{$line_seq}\t$line_seq\) are same!\n";
			exit;
			}
}

# read from fastq file(s)
my($rd1,$rd2);
if($read2=~/fq$/){
	open $rd1,$read1 or die $!;
	open $rd2,$read2 or die $!;
}elsif($read2=~/fq.gz$/){
	open $rd1,"zcat $read1|" or die $!;
}   open $rd2,"zcat $read2|" or die $!;

while(<$rd1>){
	my $head1= $_;
	my $seq1 = <$rd1>;
	my $plus1= <$rd1>;
	my $qual1= <$rd1>;
	my $head2= <$rd2>;
	my $seq2 = <$rd2>;
	my $plus2= <$rd2>;
	my $qual2= <$rd2>;
	$totalReadsNum ++;
	chomp($head1,$seq1,$plus1,$qual1,$head2,$seq2,$plus2,$qual2);
	my $barseq=substr($seq2,$fc-1,$barcode_len+1);
	$tagNum{$barseq} ++;
	$ori_tag{$barseq} = substr($seq2,$fc-1,$barcode_len1+1)."-".substr($seq2,$fc+$barcode_len1-1,$barcode_len2+1);
	if(exists $barhash{$barseq}){
	    # trim barcode from seq and qual lines
		my $spitseq2=substr($seq2,0,$fc-1).substr($seq2,$fc+$barcode_len-1,);
		my $spitqual2=substr($qual2,0,$fc-1).substr($qual2,$fc+$barcode_len-1,);
		$oh{$barhash{$barseq}}[0]->print("$head1\n$seq1\n$plus1\n$qual1\n");
		$oh{$barhash{$barseq}}[1]->print("$head2\n$spitseq2\n$plus2\n$spitqual2\n");
		if(exists $oribar{$barseq}){
			$correctBar{$barhash{$barseq}} +=1;
		}else{
			$correctedBar{$barhash{$barseq}} +=1;
		}
	}else{
		my $spitseq2=substr($seq2,0,$fc-1).substr($seq2,$fc+$barcode_len-1,);
		my $spitqual2=substr($qual2,0,$fc-1).substr($qual2,$fc+$barcode_len-1,);
		print $am1 "$head1\n$seq1\n$plus1\n$qual1\n";
		print $am2 "$head2\n$spitseq2\n$plus2\n$spitqual2\n";
		$unknownBar{$barseq} +=1;
	}
}
close $rd1;close $rd2;
close $am1;close $am2;

###############################
# print results to BarcodeStat
###############################

my($totalcorrect,$totalcorrected,$totalbarreads,$totalpct);
for my $seq(sort {$barhash{$a} cmp $barhash{$b}} keys %oribar){
    # when $correctedBar{$barhash{$seq}} is unset => set it to 0
    if (!$correctedBar{$barhash{$seq}}) {
        $correctedBar{$barhash{$seq}}=0;
    }
    my $BartotalReads = $correctBar{$barhash{$seq}}+$correctedBar{$barhash{$seq}};
	my $pct = ($BartotalReads/$totalReadsNum)*100;
	$totalcorrect += $correctBar{$barhash{$seq}};
	$totalcorrected+=$correctedBar{$barhash{$seq}};
	$totalbarreads += $BartotalReads;
	$totalpct += $pct;
	printf $BS "%s\t%d\t%d\t%d\t%.4f%%\n",$barhash{$seq},$correctBar{$barhash{$seq}},$correctedBar{$barhash{$seq}},$BartotalReads,$pct;
}
printf $BS "Total\t%d\t%d\t%d\t%.4f%%\n",$totalcorrect,$totalcorrected,$totalbarreads,$totalpct;
close $BS;

###############################
# print results to TagStat
###############################

for my $seq(sort {$tagNum{$b}<=>$tagNum{$a}} keys %tagNum){
	my $pct=($tagNum{$seq}/$totalReadsNum)*100;
	if(exists $barhash{$seq}){
		printf $SS "%s\t%s\t%d\t%.2f%%\n",$ori_tag{$seq},$barhash{$seq},$tagNum{$seq},$pct;
	}
	else{
		printf $SS "%s\tunknown\t%d\t%.2f%%\n",$ori_tag{$seq},$tagNum{$seq},$pct;
	}
}
close $SS;

#=============subroutine==================
sub bar_hash{
	my ($seq,$name,$errnum,$hash)=@_;
	my ($tmp_seq1,$tmp_seq2);
	my @bases=('A','T','C','G','N');
	if($errnum==0){
		$hash->{$seq} =$name;
		return $hash;
	}else{
		my $seq1 = substr($seq,0,length($seq)/2);	#barcode 1
		my $seq2 = substr($seq,length($seq)/2,);	#barcode 2
		for (my $i=0;$i<length($seq1);$i++){
			for (my $j=0;$j<length($seq2);$j++){
				for (my $k=0;$k<@bases;$k++){
					$tmp_seq1 =substr($seq1,0,$i).$bases[$k].substr($seq1,$i+1,);
					$tmp_seq2 =substr($seq2,0,$j).$bases[$k].substr($seq2,$j+1,);
					my $tmp_seqs = $tmp_seq1.$tmp_seq2;
					if($errnum > 1){
						&bar_hash($tmp_seqs,$name,$errnum-1,$hash);
					}
					else{
						$hash->{$tmp_seqs}=$name;
					}
				}
			}
		}
		return $hash;
	}
}
