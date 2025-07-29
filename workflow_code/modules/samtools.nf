process SAMTOOLS_SORT {
  publishDir "${params.outdir}/${params.gldsAccession}/02-Bismark_Alignment",
    mode: params.publish_dir_mode
  tag "${ meta.id }"

  input:
    tuple val(meta), path(bam)

  output:
    tuple val(meta), path("${ meta.id }/*_sorted${ params.suffix }"), emit: sorted_bam
    path("${ meta.id }/*_sorted${ params.suffix }.bai")

  script:
    """
    mkdir -p ${ meta.id } &&
    samtools sort -@${task.cpus} \
        --write-index \
        --no-PG \
        -o ${ meta.id }/${ meta.id }_bismark_bt2_sorted${ params.suffix }##idx##${ meta.id }/${ meta.id }_bismark_bt2_sorted${ params.suffix }.bai \
        ${ bam }
    """
}
