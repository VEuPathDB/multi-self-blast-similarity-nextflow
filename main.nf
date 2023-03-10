#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//---------------------------------------------------------------
// Param Checking 
//---------------------------------------------------------------

if(params.tarFile) {
  tarFile = Channel.fromPath( params.tarFile )
}
else {
  throw new Exception("Missing params.tarFile")
}

//--------------------------------------------------------------------------
// Includes
//--------------------------------------------------------------------------

include { multiBlastSelf } from './modules/multiSelfBlastSimilarity.nf'

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------

workflow {
  
  multiBlastSelf(tarFile)

}

