#!/usr/bin/env bash

set -euo pipefail
 
diamond blastp \
	-d $database \
	-q $fasta \
	-o out.txt \
	-f 6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length nident pident positive qframe qstrand gaps qseq \
	--comp-based-stats 0 \
	$blastArgs

if [ "$printSimSeqs" = true ]; then

    perl /usr/bin/printSimSeqs.pl \
     --result out.txt \
     --output diamondSimilarity.out \
     --minLen $lengthCutoff \
     --minPercent $percentCutoff \
     --minPval $pValCutoff \
     --remMaskedRes $adjustMatchLength
    
else

    perl /usr/bin/diamondSimilarity.pl \
     --fasta $fasta \
     --result out.txt \
     --output diamondSimilarity.out \
     --minLen $lengthCutoff \
     --minPercent $percentCutoff \
     --minPval $pValCutoff \
     --remMaskedRes $adjustMatchLength \
     --outputType $outputType
    
fi




