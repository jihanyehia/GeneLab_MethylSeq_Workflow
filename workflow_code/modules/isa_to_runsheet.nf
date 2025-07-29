process ISA_TO_RUNSHEET {
    tag "${params.osdAccession}_${params.gldsAccession}"

    publishDir "${params.outdir}/${params.gldsAccession}/Metadata",
        mode: params.publish_dir_mode

    input: 
        path isa_archive
        path dp_tools_plugin

    output:
        path "*.csv", emit: runsheet

    script:
    """
    dpt-isa-to-runsheet --accession ${params.osdAccession} --isa-archive ${isa_archive} --plugin-dir ${dp_tools_plugin}
    """
}
