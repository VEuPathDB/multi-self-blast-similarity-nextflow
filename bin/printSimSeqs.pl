#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(first_index);

my ($result,$output, $minLen, $minPercent, $minPval, $remMaskedRes);

&GetOptions("result=s"=> \$result,
            "output=s"=> \$output,
            "minLen=i"=> \$minLen,
            "minPercent=i"=> \$minPercent,
            "minPval=f"=> \$minPval,
            "remMaskedRes=s" => \$remMaskedRes);

$minPval = $minPval ? $minPval : 1e-5;
$minLen = $minLen ? $minLen : 10;
$minPercent = $minPercent ? $minPercent : 20;

open(my $data, '<', $result) || die "Could not open file $result: $!";
open(OUT,">$output");
open(LOG,">diamondSimilarity.log");

while (my $line = <$data>) {
    chomp $line;
    my ($qseqid,$qlen,$sseqid,$slen,$qstart,$qend,$sstart,$send,$evalue,$bitscore,$length,$nident,$pident,$positive,$qframe,$qstrand,$gaps,$qseq) = split(/\t/, $line);
    next if $sseqid eq $qseqid;
    if ($remMaskedRes eq 'true') {
        $qseq =~ s/X//g;
        length($qseq) >= 10 ? $length = length($qseq) : print LOG "RemoveMaskedResiduesFromLength: error...length is less than 10 following removal of X's\n";
    }
    if ($length >= $minLen && $pident >= $minPercent && $evalue <= $minPval) {
	my $nonOverlappingPercent;
	my $roundedPercent;
	my $qTaxAbbrev = $qseqid;
	$qTaxAbbrev =~ s/\|.*//g;
	my $sTaxAbbrev = $sseqid;
	$sTaxAbbrev =~ s/\|.*//g;	
        if ($slen < $qlen) {
            $nonOverlappingPercent = ($slen - $gaps)/$slen * 100;
            $roundedPercent = sprintf("%.2f", $nonOverlappingPercent);
        }
        if ($slen >= $qlen) {
            $nonOverlappingPercent = ($qlen - $gaps)/$qlen * 100;
            $roundedPercent = sprintf("%.2f", $nonOverlappingPercent);
        }
	my ($mant, $exp);
        if ($evalue =~ /e/) {
            my $pValue = $evalue =~ /^e/ ?  '1' . $evalue  : $evalue ;
            ($mant, $exp) = split(/e/, $pValue);
        } else {
	    $mant = int($evalue);
	    $exp = 0;
        }
	print OUT "$qseqid $sseqid $qTaxAbbrev $sTaxAbbrev $mant $exp $pident $roundedPercent\n";
    }
}	
