process PARSE_ANNOTATIONS_TABLE {
  // Extracts data from this kind of table: 
  // https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv

  input:
    val(annotations_csv_url_string)
    val(organism_sci)
  
  output:
    tuple val(organism_sci), val(fasta_url), val(gtf_url), emit: reference_genome_urls
    val(annotations_db_url), emit: annotations_db_url
    tuple val(ensemblVersion), val(ensemblSource), emit: reference_version_and_source
    val(simple_organism_name), emit: simple_organism_name
  
  exec:
    def organisms = [:]
    println "Fetching table from ${annotations_csv_url_string}"
    
    // download data to memory
    annotations_csv_url_string.toURL().splitEachLine(",") {fields ->
          organisms[fields[1]] = fields
    }
    // extract required fields
    organism_key = organism_sci.capitalize().replace("_"," ")
    simple_organism_name = organisms[organism_key][0]
    fasta_url = organisms[organism_key][5]
    gtf_url = organisms[organism_key][6]
    annotations_db_url = organisms[organism_key][10]
    ensemblVersion = organisms[organism_key][3]
    ensemblSource = organisms[organism_key][4]

    println "PARSE_ANNOTATIONS_TABLE:"
    println "Values parsed for '${organism_key}' using process:"
    println "--------------------------------------------------"
    println "- fasta_url: ${fasta_url}"
    println "- gtf_url: ${gtf_url}"
    println "- annotations_db_url: ${annotations_db_url}"
    println "- ensemblVersion: ${ensemblVersion}"
    println "- ensemblSource: ${ensemblSource}"
    println "--------------------------------------------------"
}

process DOWNLOAD_GUNZIP_REFERENCES {
  // Download and decompress genome and annotation files
  tag "Organism: ${ organism_sci }, Ensembl Version: ${ensemblVersion}"
  label 'networkBound'
  storeDir "${ params.derivedStorePath }/BismarkIndices_BT/${ ref_source }_release${ensemblVersion}/${ organism_sci.capitalize() }"

  input:
    tuple val(organism_sci), val(fasta_url), val(gtf_url)
    tuple val(ensemblVersion), val(ref_source)
  
  output:
    path("*.fa")

  script:
  """
  wget ${fasta_url}
  gunzip *.gz
  """
}
