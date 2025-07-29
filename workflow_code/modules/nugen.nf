process NUGEN {
  // Stages the raw reads into appropriate publish directory
  publishDir "${params.outdir}/${params.gldsAccession}/01-TrimFilter/Fastq",
      pattern: "*trimmed.fastq.gz",
      mode: params.publish_dir_mode
  tag "${ meta.id }"
  label 'low_cpu_med_memory'

  input:
    tuple val(meta), path("input/*")
    path("trimRRBSdiversityAdaptCustomers.py")

  output:
    tuple val(meta), path("${ meta.id }*trimmed.fastq.gz"), emit: reads

  script:
    if (meta.paired_end) {
    """
    python trimRRBSdiversityAdaptCustomers.py -1 input/${ meta.id }_R1_trimmed.fastq.gz -2 input/${ meta.id }_R2_trimmed.fastq.gz

    mv input/${ meta.id }_R1_trimmed.fastq_trimmed.fq.gz ${ meta.id }_R1_trimmed.fastq.gz
    mv input/${ meta.id }_R2_trimmed.fastq_trimmed.fq.gz ${ meta.id }_R2_trimmed.fastq.gz
    """
    } else {
    """
    python trimRRBSdiversityAdaptCustomers.py -1 input/${ meta.id }_trimmed.fastq.gz

    mv input/${ meta.id }_trimmed.fastq_trimmed.fq.gz ${ meta.id }_trimmed.fastq.gz
    """
    }
    
}
