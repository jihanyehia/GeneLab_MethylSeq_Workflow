process BUILD_BISMARK {
  tag "Refs: ${ genomeFasta }, Ensembl Version: ${ensemblVersion}"
  storeDir "${ params.derivedStorePath }/BismarkIndices_BT/${ ref_source }_release${ensemblVersion}"

  label 'maxCPU'
  label 'big_mem'

  input:
    path(genomeFasta)
    val(meta)
    tuple val(ensemblVersion), val(ref_source) // Used for defining storage location 

  output:
    path("${ meta.organism_sci.capitalize() }"), emit: build
    path("${ meta.organism_sci.capitalize() }/genomic_nucleotide_frequencies.txt")
    path("${ meta.organism_sci.capitalize() }/Bisulfite_Genome/")

  script:
    def genome_path = "${ params.derivedStorePath }/BismarkIndices_BT/${ ref_source }_release${ensemblVersion}/${ meta.organism_sci.capitalize() }"
    """
    bam2nuc --genome_folder ${ genome_path } --genomic_composition_only
    bismark_genome_preparation --bowtie2 --parallel ${task.cpus} ${ genome_path }
    """
}

process ALIGN_BISMARK {
  publishDir "${params.outdir}/${params.gldsAccession}/02-Bismark_Alignment",
    mode: params.publish_dir_mode

  tag "Sample: ${ meta.id }"

  input:
    tuple val( meta ), path( reads )
    path(bismark_index_dir)

  output:
    tuple val(meta), path("${ meta.id }/*.bam"), emit: bam
    path("${ meta.id }/*nucleotide_stats.txt"), emit: stats
    path("${ meta.id }/*report.txt"), emit: report
    path("versions.txt"), emit: version

  script:
    def input = meta.paired_end ? "-1 ${ reads[0] } -2 ${ reads[1] }" : "${ reads }"
    """
    bismark --bowtie2 \
      --bam \
      --parallel ${task.cpus} \
      --gzip \
      --non_bs_mm \
      ${ params.non_directional ? '--non_directional' : '' } \
      --nucleotide_coverage \
      --output_dir ${ meta.id }/ \
      --temp_dir /tmp \
      --genome_folder ${ bismark_index_dir } \
      ${ input }

    # remove R1/R2 and _trimmed in filename
    ${ meta.paired_end ? \
      "mv ${ meta.id }/${ meta.id }*_PE_report.txt ${ meta.id }/${ meta.id }_bismark_bt2_PE_report.txt; \
      mv ${ meta.id }/${ meta.id }*.bam ${ meta.id }/${ meta.id }_bismark_bt2_pe.bam" : \
       "mv ${ meta.id }/${ meta.id }*_SE_report.txt ${ meta.id }/${ meta.id }_bismark_bt2_SE_report.txt; \
      mv ${ meta.id }/${ meta.id }*.bam ${ meta.id }/${ meta.id }_bismark_bt2.bam" }

    mv ${ meta.id }/${ meta.id }*.nucleotide_stats.txt ${ meta.id }/${ meta.id }_bismark_bt2.nucleotide_stats.txt

    echo ALIGN_BISMARK_version: `bismark --version` > versions.txt
    """
}

process DEDUPE {
  publishDir "${params.outdir}/${params.gldsAccession}/02-Bismark_Alignment",
    mode: params.publish_dir_mode
  tag "${ meta.id }"

  label 'maxCPU'

  input:
    tuple val(meta), path(bam)

  output:
    tuple val(meta), path("${ meta.id }/*deduplicated.bam"), emit: bam
    path("${ meta.id }/*deduplication_report.txt"), emit: report

  script:
    """
    deduplicate_bismark \
      --output_dir ${ meta.id } \
      ${ bam }
    """
}

process EXTRACT_CALLS {
  publishDir "${params.outdir}/${params.gldsAccession}/03-Bismark_Methylation_Calls",
    mode: params.publish_dir_mode
  tag "${ meta.id }"

  label 'medCPU'
  label 'big_mem'

  input:
    tuple val(meta), path(bam), path(bismark_index_dir)

  output:
    path("${ meta.id }/*M-bias.txt"), emit: bias
    path("${ meta.id }/*splitting_report.txt"), emit: report
    path("${ meta.id }/*bismark.cov.gz"), emit: cov
    path("${ meta.id }/*")

  script:
    // https://www.biostars.org/p/9594498/ - --genome_folder needs to be absolute path
    """
    bismark_methylation_extractor --parallel ${task.cpus} \
      --bedGraph \
      --gzip \
      --comprehensive \
      --output_dir ${ meta.id } \
      --cytosine_report \
      --genome_folder \${PWD}/${bismark_index_dir} \
      ${ bam }
    """
}

process BISMARK_REPORT {
  publishDir "${params.outdir}/${params.gldsAccession}/04-Bismark_Reports/Individual_Sample_Reports",
    mode: params.publish_dir_mode
  tag "${ meta.id }"

  label 'maxCPU'

  input:
    val(meta)
    path(align_report)
    path(align_stats)
    path(call_report)
    path(call_bias)
    path(align_dedupe_report)

  output:
    path("*_report.html"), emit: html

  script:
    def dedup_report = align_dedupe_report.name != "NO_FILE" ? "--dedup_report ${align_dedupe_report}" : ""
    """
    bismark2report --dir . \
      --alignment_report ${ align_report } \
      --nucleotide_report ${ align_stats } \
      ${ dedup_report } \
      --splitting_report ${ call_report } \
      --mbias_report ${ call_bias }
    """
}

process BISMARK_SUMMARY {
  publishDir "${params.outdir}/${params.gldsAccession}/04-Bismark_Reports",
    mode: params.publish_dir_mode
  tag "Dataset-wide"

  label 'maxCPU'

  input:
    path(bam)
    path(align_report)
    path(align_stats)
    path(call_report)
    path(call_bias)
    path(align_dedupe_report)

  output:
    path("bismark_summary_report.txt")
    path("bismark_summary_report.html")

  script:
    """
    bismark2summary *bam
    """
}
