process QUALIMAP {
  publishDir "${params.outdir}/${params.gldsAccession}/02-Bismark_Alignment",
    mode: params.publish_dir_mode
  tag "${ meta.id }"
  label 'big_mem'

  input:
    path(genomeGtf)
    tuple val(meta), path(sorted_bam)

  output:
    path("${ meta.id }/${ meta.id }_bismark_bt2_qualimap")

  script:
    """
    qualimap bamqc -bam ${ sorted_bam } \
        -gff ${genomeGtf} \
        -outdir ${ meta.id }/${ meta.id }_bismark_bt2_qualimap/ \
        --collect-overlap-pairs \
        --java-mem-size=8G \
        -nt ${task.cpus}
    """
}
