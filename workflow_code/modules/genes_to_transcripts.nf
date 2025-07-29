process GENES_TO_TRANSCRIPTS {
  storeDir "${ params.derivedStorePath }/BismarkIndices_BT2/${ ref_source }_release${ensemblVersion}/${ meta.organism_sci.capitalize() }"

  input:
    path(genome_gtf)
    val(meta)
    tuple val(ensemblVersion), val(ref_source) // Used for defining storage location 

  output:
    path("${ genome_gtf.baseName }-gene-to-transcript-map.tsv")

  script:
  """
  awk ' \$3 == \"transcript\" ' ${ genome_gtf } | cut -f 9 | tr -s \";\" \"\t\" | \
    cut -f 1,3 | tr -s \" \" \"\t\" | cut -f 2,4 | tr -d '\"' \
    > ${ genome_gtf.baseName }-gene-to-transcript-map.tsv
  """
}
