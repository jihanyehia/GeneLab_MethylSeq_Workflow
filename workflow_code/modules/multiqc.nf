process MULTIQC {
    publishDir "${params.outdir}/${params.gldsAccession}/${params.multiqc_publish_dir}",
      pattern: '*.{html,zip}',
      mode: params.publish_dir_mode
    tag("Dataset-wide")
    
    input:
      path(sample_names)
      path("mqc_in/*") // any number of multiqc compatible files
      path(multiqc_config)
    
    output:
      path("${ params.MQCLabel }_multiqc${ params.assay_suffix }_data.zip"), emit: zipped_data
      path("${ params.MQCLabel }_multiqc${ params.assay_suffix }.html"), emit: html
      path("versions.yml"), emit: versions

    script:
    def config_arg = multiqc_config.name != "NO_FILE" ? "--config ${multiqc_config}" : ""
    """
    multiqc \
            --interactive -o . \
            -n ${ params.MQCLabel }_multiqc${ params.assay_suffix } \
            ${ config_arg } \
            mqc_in
    
    # Clean paths and create zip
    clean_multiqc_paths.py ${ params.MQCLabel }_multiqc${ params.assay_suffix }_data .

    echo "${task.process}:" > versions.yml
    echo "    multiqc: \$(multiqc --version | sed -e "s/multiqc, version //g")" >> versions.yml
    """
}
