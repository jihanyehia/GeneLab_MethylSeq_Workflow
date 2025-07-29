process TO_PRED {
  // Converts reference gtf into pred 
  storeDir "${ params.derivedStorePath }/BismarkIndices_BT2/${ ref_source }_release${ensemblVersion}/${ meta.organism_sci.capitalize() }"

  input:
    path(genome_gtf)
    val(meta)
    tuple val(ensemblVersion), val(ref_source) // Used for defining storage location 

  output:
    path("${ genome_gtf.baseName }.genepred")

  script: // https://github.com/nextflow-io/nextflow/issues/1359
  """
  gtfToGenePred ${ genome_gtf } ${ genome_gtf.baseName }.genepred
  """
}

process TO_BED {
  // Converts reference genePred into Bed format
  storeDir "${ params.derivedStorePath }/BismarkIndices_BT2/${ ref_source }_release${ensemblVersion}/${ meta.organism_sci.capitalize() }"

  input:
    path(genome_pred)
    val(meta)
    tuple val(ensemblVersion), val(ref_source) // Used for defining storage location 

  output:
    path("${ genome_pred.baseName }.bed")

  script:
  """
  genePredToBed ${ genome_pred } ${ genome_pred.baseName }.bed
  """
}
