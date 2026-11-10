process ISA_TO_RUNSHEET {
    tag "${osd_accession}_${glds_accession}"

    input: 
    val(osd_accession)
    val(glds_accession)
    path(isa_archive), stageAs: 'input_isa/*'
    path(dp_tools_plugin)

    output:
    path("*.csv"), emit: runsheet
    path("*.zip"), emit: staged_isa_archive
    path("versions.yml"), emit: version

    script:
    """
    dpt-isa-to-runsheet --accession ${osd_accession} --isa-archive ${isa_archive} --plugin-dir ${dp_tools_plugin}

    # Copy the ISA archive to the output directory
    cp ${isa_archive} .

    echo '"${task.process}":' > versions.yml
    echo "    dp_tools: \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.yml
    """
}
