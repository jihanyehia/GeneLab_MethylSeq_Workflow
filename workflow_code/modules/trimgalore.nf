process TRIMGALORE {
  // Stages the raw reads into appropriate publish directory
  publishDir "${params.outdir}/${params.gldsAccession}/01-TrimFilter/Fastq",
      pattern: "*trimmed.fastq.gz",
      mode: params.publish_dir_mode
  publishDir "${params.outdir}/${params.gldsAccession}/01-TrimFilter/Trimming_Reports",
      pattern: "*.txt",
      mode: params.publish_dir_mode
  tag "${ meta.id }"

  input:
    tuple val(meta), path(reads)

  output:
    tuple val(meta), path("${ meta.id }*trimmed.fastq.gz"), emit: reads
    path("${ meta.id }*.txt"), emit: reports
    path("versions.yml"), emit: versions

  script:

    """
    trim_galore --gzip \
    --cores $task.cpus \
    --phred33 \
    ${ meta.paired_end ? '--paired' : '' } \
    ${ params.nugen ? '-a AGATCGGAAGAGC' : '' } \
    ${ meta.rrbs && !params.nugen ? '--rrbs' : '' } \
    ${ meta.rrbs && params.non_directional ? '--non_directional' : '' } \
    --clip_R1 ${ params.clip } --clip_R2 ${ params.clip } --three_prime_clip_R1 ${ params.clip } --three_prime_clip_R2 ${ params.clip } \
    $reads \
    --output_dir .

    # rename with _trimmed suffix
    ${ meta.paired_end ? \
      "mv ${ meta.id }_R1_raw_val_1.fq.gz ${ meta.id }_R1_trimmed.fastq.gz; \
      mv ${ meta.id }_R2_raw_val_2.fq.gz ${ meta.id }_R2_trimmed.fastq.gz" : \
      "mv ${ meta.id }_raw_trimmed.fq.gz ${ meta.id }_trimmed.fastq.gz" }

    echo '"${task.process}":' > versions.yml
    echo "    trimgalore: \$(trim_galore -v | sed -n 's/.*version //p' | head -n 1)" >> versions.yml
    echo "    cutadapt: \$(cutadapt --version)" >> versions.yml
    """
}
