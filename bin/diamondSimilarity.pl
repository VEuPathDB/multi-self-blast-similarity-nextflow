#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(first_index);

my ($fasta,$result,$output, $minLen, $minPercent, $minPval, $remMaskedRes,$outputType);

&GetOptions("fasta=s"=> \$fasta,
            "result=s"=> \$result,
            "output=s"=> \$output,
            "minLen=i"=> \$minLen,
            "minPercent=i"=> \$minPercent,
            "minPval=f"=> \$minPval,
            "remMaskedRes=s" => \$remMaskedRes,
            "outputType=s" => \$outputType);

$minPval = $minPval ? $minPval : 1e-5;
$minLen = $minLen ? $minLen : 10;
$minPercent = $minPercent ? $minPercent : 20;
$outputType = $outputType ? $outputType : "both";

my $printSum = 0;
my $printSpan = 0;
if($outputType =~ /sum/i){
    $printSum = 1;
}elsif($outputType =~ /span/i){
    $printSpan = 1;
}elsif($outputType =~ /both/i){
    $printSum = 1;
    $printSpan = 1;
}

my %subjectCountHash = &getSubjectCount($fasta,$result,$minLen,$minPercent,$minPval);
my @deflines = &getSubjectIds($fasta);

open(my $data, '<', $result) || die "Could not open file $result: $!";
open(OUT,">$output");
open(LOG,">diamondSimilarity.log");

my $previousQseqId;
my $previousSseqId = "hello";
my $counter = 0;

while (my $line = <$data>) {
    chomp $line;

    my ($qseqid,$qlen,$sseqid,$slen,$qstart,$qend,$sstart,$send,$evalue,$bitscore,$length,$nident,$pident,$positive,$qframe,$qstrand,$gaps,$qseq) = split(/\t/, $line);
    my $subjectCount = $subjectCountHash{">$qseqid"};

    if ($counter == 0) {
	&addSubjectsNoCount($previousQseqId,$qseqid,1,@deflines);
	$subjectCount == 0 ? print LOG "No HSPS passed requirements for $qseqid\n" : print OUT ">$qseqid ($subjectCount subjects)\n";
	$previousQseqId = $qseqid;
    }

    if ($qseqid eq $previousQseqId) {

	if ($previousSseqId eq $sseqid) {
	    print LOG "Multiple HSPS for $qseqid, $sseqid\n";
	    die;
	}
        if ($remMaskedRes eq 'true') {
	        $qseq =~ s/X//g;
	        length($qseq) >= 10 ? $length = length($qseq) : print LOG "RemoveMaskedResiduesFromLength: error...length is less than 10 following removal of X's\n";
	}
	if ($length >= $minLen && $pident >= $minPercent && $evalue <= $minPval) {
	    my $nonOverlappingPercent;
	    my $roundedPercent;
	    if ($slen < $qlen) {
		$nonOverlappingPercent = ($slen - $gaps)/$slen * 100;
		$roundedPercent = sprintf("%.2f", $nonOverlappingPercent);
	    }
	    if ($slen >= $qlen) {
                $nonOverlappingPercent = ($qlen - $gaps)/$qlen * 100;
		$roundedPercent = sprintf("%.2f", $nonOverlappingPercent);
	    }
            if ($printSum) {
		print OUT "  Sum: $sseqid:$bitscore:$evalue:$sstart:$send:$qstart:$qend:1:$length:$nident:$positive:$qstrand:$qframe:$qlen:$slen:$roundedPercent\n";
	    }
            if ($printSpan) {
		print OUT "   HSP1: $sseqid:$nident:$positive:$length:$bitscore:$evalue:$sstart:$send:$qstart:$qend:$qstrand:$qframe\n";
	    }
	}

    } else {
	  &addSubjectsNoCount($previousQseqId,$qseqid,0,@deflines);
	  $subjectCount == 0 ? print LOG "No HSPS passed requirements for $qseqid\n" : print OUT ">$qseqid ($subjectCount subjects)\n";
	  if ($remMaskedRes eq 'true') {
	      $qseq =~ s/X//g;
	      length($qseq) >= 10 ? $length = length($qseq) : print LOG "RemoveMaskedResiduesFromLength: error...length is less than 10 following removal of X's\n";
	  }
	  if ($length >= $minLen && $pident >= $minPercent && $evalue <= $minPval) {
	      my $nonOverlappingPercent;
	      my $roundedPercent;
	      if ($slen < $qlen) {
		  $nonOverlappingPercent = ($slen - $gaps)/$slen * 100;
		  $roundedPercent = sprintf("%.2f", $nonOverlappingPercent);
	      }
	      if ($slen >= $qlen) {
                  $nonOverlappingPercent = ($qlen - $gaps)/$qlen * 100;
		  $roundedPercent = sprintf("%.2f", $nonOverlappingPercent);
	      }
	      if ($printSum) {
		  print OUT "  Sum: $sseqid:$bitscore:$evalue:$sstart:$send:$qstart:$qend:1:$length:$nident:$positive:$qstrand:$qframe:$qlen:$slen:$roundedPercent\n";
	      }
              if ($printSpan) {
		  print OUT "   HSP1: $sseqid:$nident:$positive:$length:$bitscore:$evalue:$sstart:$send:$qstart:$qend:$qstrand:$qframe\n";
	      }
          }
    }
    $previousQseqId = $qseqid;
    $previousSseqId = $sseqid; 
    $counter += 1;
}

close($data);
close OUT;
close LOG;

#============================================= subroutines ==================================================================================

sub getSubjectIds{
    my ($fasta) = @_;
    my @deflines = `grep "^>" $fasta | awk 'sub(/\\slength=[0-9]+/, "")'`;
    chomp(@deflines);
    return @deflines;
}

sub getSubjectCount{
    my ($fasta,$outFile,$minLen,$minPercent,$minPval) = @_;
    my @deflines = &getSubjectIds($fasta);
    my %hash = map { $_ => 0 } @deflines;
    open(my $data, '<', $outFile) || die "Could not open file $outFile: $!";
    while (my $line = <$data>) {
        chomp $line;
        my ($qseqid,$qlen,$sseqid,$slen,$qstart,$qend,$sstart,$send,$evalue,$bitscore,$length,$nident,$pident,$positive,$qframe,$qstrand,$gaps,$qseq) = split(/\t/, $line);
	if ($length >= $minLen && $pident >= $minPercent && $evalue <= $minPval) {
	    $hash{">$qseqid"} += 1;
	}
    }
    close($data);
    return %hash;
}

sub addSubjectsNoCount {
    my ($prevQSeqId,$qSeqId,$isFirstQSeq,@subjectIds) = @_;
    my $counter = 0;
    if ($isFirstQSeq) {
        until ($subjectIds[$counter] eq ">$qSeqId") {
            print OUT "$subjectIds[$counter] (0 subjects)\n";
	    $counter += 1;
	}
    }
    else {
        my $firstIndex = first_index { $_ eq ">$prevQSeqId" } @subjectIds;
	$firstIndex += 1;
	my $currentIndex = first_index { $_ eq ">$qSeqId" } @subjectIds;
        until ($firstIndex == $currentIndex) {
	    print OUT "$subjectIds[$firstIndex] (0 subjects)\n";
	    $firstIndex += 1;
        }
    }
}
