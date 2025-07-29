// This ensures DSL2 syntax and process imports
nextflow.enable.dsl=2

include { PARSE_ANNOTATIONS_TABLE; DOWNLOAD_GUNZIP_REFERENCES } from '../modules/annotations.nf'


/**************************************************
* ACTUAL WORKFLOW  ********************************
**************************************************/
workflow REFERENCES {
  take:
    organism_sci

  main:
      // Must run in any approach to find the appropriate annotations database table
      PARSE_ANNOTATIONS_TABLE( params.reference_table, organism_sci)

      if (params.ref_fasta && params.ref_gtf) {
        genome_annotations = Channel.fromPath([params.ref_fasta, params.ref_gtf], checkIfExists: true).toList()
        Channel.value( [params.ensemblVersion, params.ref_source] ) | set { ch_ref_source_version }
      } else {
        // use assets table to find current fasta and gtf urls and associated metadata about those reference files
        
        DOWNLOAD_GUNZIP_REFERENCES( 
          PARSE_ANNOTATIONS_TABLE.out.reference_genome_urls,
          PARSE_ANNOTATIONS_TABLE.out.reference_version_and_source,
        )
        DOWNLOAD_GUNZIP_REFERENCES.out | set{ genome_annotations }
        PARSE_ANNOTATIONS_TABLE.out.reference_version_and_source | set { ch_ref_source_version }
      }

      // use assets table to find current annotations file
      PARSE_ANNOTATIONS_TABLE.out.annotations_db_url | set{ ch_gene_annotations_url }
      PARSE_ANNOTATIONS_TABLE.out.simple_organism_name | set{ ch_simple_organism_name }

  emit:
      genome_annotations = genome_annotations
      gene_annotations = ch_gene_annotations_url
      simple_organism_name = ch_simple_organism_name
      reference_version_and_source = ch_ref_source_version
}
