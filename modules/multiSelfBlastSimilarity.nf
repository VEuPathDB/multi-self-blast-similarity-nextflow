#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process retrieveFastas {
  input:
    path tarfile

  output:
    path '*'

  script:
    """
    if [[ "$tarfile" == *".gz"* ]];then
        tar xzf $tarfile
    else
        tar xf $tarfile
    fi
    mv **/* .
    rm -r */    
    """
}

process createDatabase {
  input:
    path newdbfasta

  output:
    path 'newdb.dmnd', emit: database
    path newdbfasta, emit: fasta

  script:
    template 'createDatabase.bash'
}

process diamondSimilarity {
  input:
    path fasta
    path database
    val pValCutoff 
    val lengthCutoff 
    val percentCutoff 
    val blastProgram  
    val blastArgs 
    val adjustMatchLength
    val outputType
    val printSimSeqs

  output:
    path 'diamondSimilarity.out', emit: output_file
    path 'diamondSimilarity.log', emit: log_file

  script:
    template 'diamondSimilarity.bash'
}

process sortOutput {
  publishDir params.outputDir, saveAs: {filename->params.dataFile}
  
  input:
    path output
        
  output:
    path 'diamondSimilarity.out'

  script:
    """
    sed 's/^>/\\x00&/' $output  | sort -z | tr -d '\\0' > diamondSimilarity.out    
    """
}

process sortSimSeqs {
  publishDir params.outputDir, saveAs: {filename->params.dataFile}
  
  input:
    path output
        
  output:
    path 'diamondSimilarity.out'

  script:
    """
    cat $output | sort -k 1 > diamondSimilarity.out
    """
}

workflow multiBlastSelf {
  take:
    tarFile

  main:
    seqs = retrieveFastas(tarFile) | flatten()
    database = createDatabase(seqs)
    diamondSimilarityResults = diamondSimilarity(database.fasta, database.database, params.pValCutoff, params.lengthCutoff, params.percentCutoff, params.blastProgram, params.blastArgs, params.adjustMatchLength, params.outputType, params.printSimSeqs)
    if (params.printSimSeqs) {
       diamondSimilarityResults.output_file | collectFile(name: 'similarity.out') | sortSimSeqs
    }
    else {
       diamondSimilarityResults.output_file | collectFile(name: 'similarity.out') | sortOutput
    }
    diamondSimilarityResults.log_file | collectFile(storeDir: params.outputDir, name: params.logFile)
}