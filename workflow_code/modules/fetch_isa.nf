process FETCH_ISA {

    tag "${params.osdAccession}"

    publishDir "${params.outdir}/${params.gldsAccession}/Metadata",
        mode: params.publish_dir_mode

    output:
        path "*.zip", emit: isa_archive

    script:
    """
    dpt-get-isa-archive --accession ${ params.osdAccession }
    """
}
