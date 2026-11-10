#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

include { GET_ACCESSIONS } from './modules/get_accessions.nf'
include { PARSE_RUNSHEET } from './subworkflows/parse_runsheet.nf'
include { CLEAN_MULTIQC_PATHS } from './modules/clean_multiqc_paths.nf'
include { CLEAN_OUTPUT_PATHS } from './modules/clean_output_paths.nf'
include { PACKAGE_PROCESSING_INFO } from './modules/package_processing_info.nf'
include { GENERATE_RAW_PROTOCOL } from './modules/generate_protocols.nf'
include { GENERATE_README } from './modules/generate_readme.nf'
include { GENERATE_MD5SUMS } from './modules/generate_md5sums.nf'
include { VALIDATE_PROCESSING } from './modules/validate_processing.nf'
include { UPDATE_ASSAY_TABLE } from './modules/update_assay_table.nf'


workflow {
    main:

    // making sure specific expected nextflow version is being used
    if( ! nextflow.version.matches( workflow.manifest.nextflowVersion ) ) {
        println "\n    This workflow requires Nextflow version $workflow.manifest.nextflowVersion, but version $nextflow.version is currently active."
        println "\n    You can set the proper version for this terminal session by running: `export NXF_VER=$workflow.manifest.nextflowVersion`"

        println "\n  Exiting for now.\n"
        exit 1
    }

    if ( params.truncate_to ) {

        println( "\n    WARNING: DEBUG OPTIONS ARE ENABLED!\n" )

        params.truncate_to ? println( "        - Truncating reads to first ${ params.truncate_to } records" ) : null
        params.limit_samples_to ? println( "        - Limiting number of samples to process to ${ params.limit_samples_to }" ) : null
        params.force_single_end ? println( "        - Forcing single-end processing (using read 1 only)" ) : null

        println ""

    }

    osd_accession = channel.value('null')
    glds_accession = channel.value('null')
    ch_directory_name = channel.value('results')
    if ( params.accession ) {
        GET_ACCESSIONS( params.accession, params.api_url )
        osd_accession = GET_ACCESSIONS.out.accessions_txt.map { txt -> txt.readLines()[0].trim() }
        glds_accession = GET_ACCESSIONS.out.accessions_txt.map { txt ->
            def accession = txt.readLines()[1].trim()
            if (!(accession ==~ /GLDS-\d+/)) {
                error "Expected exactly one GLDS accession from GET_ACCESSIONS, but got: ${accession}"
            }
            accession
        }
        ch_directory_name = glds_accession
    }

    ch_root_directory = ch_directory_name.map { name ->
        def directory = file("${launchDir}/${name}")
        if (!directory.isDirectory()) {
            error "Main workflow output directory '${directory}' does not exist or is not a directory.\n" +
                  "Make sure you have run the main workflow first before running post-processing."
        }
        println "Resolved output directory: ${directory}"
        directory
    }

    ch_runsheet = ch_root_directory.flatMap { directory ->
        file("${directory}/Metadata/*_runsheet.csv")
    }

    ch_multiqc_data = ch_root_directory.flatMap { directory ->
        file("${directory}/*/MultiQC_Reports/*_data", type: 'dir')
    }

    ch_unpurged_files = ch_root_directory.flatMap { directory ->
        file("${directory}/Aligned_Sequence_Data/Alignment_Reports/*report.txt")
    }

    ch_processing_info = channel.fromPath("$launchDir/processing_scripts", type: 'dir')

    ch_software_versions = ch_root_directory.flatMap { directory ->
        file("${directory}/GeneLab/software_versions_*.md")
    }

    ch_isa_zip = ch_root_directory.flatMap { directory ->
        file("${directory}/Metadata/*-ISA.zip")
    }

    PARSE_RUNSHEET( ch_runsheet )
    samples = PARSE_RUNSHEET.out.samples
    runsheet_path = PARSE_RUNSHEET.out.runsheet

    // Extract metadata from the first sample and set it as a channel
    samples | first
            | map { meta, _reads -> meta }
            | set { ch_meta }

    // Extract assay suffix from metadata to be used in processes where meta is not passed
    ch_meta | map { meta -> meta.assay_suffix }
            | set { assay_suffix }

    ch_multiqc_data_mapped = ch_multiqc_data.map { file -> [file, file.parent] }
    CLEAN_MULTIQC_PATHS( ch_multiqc_data_mapped )

    ch_unpurged_files_mapped = ch_unpurged_files.map { file -> [file, file.parent] }
    CLEAN_OUTPUT_PATHS( ch_unpurged_files_mapped )
    
    PACKAGE_PROCESSING_INFO( ch_processing_info, assay_suffix )

    // Generate protocol for raw data only when params.post_processing.include_raw_data is true, otherwise skip
    ch_raw_protocol = channel.empty()
    if( params.post_processing.include_raw_data ) {
        GENERATE_RAW_PROTOCOL( ch_meta, ch_software_versions )
        ch_raw_protocol = GENERATE_RAW_PROTOCOL.out.raw_protocol
    }
    
    GENERATE_README( ch_meta, osd_accession )

    ch_ready = GENERATE_README.out.readme                    // README must exist
                .combine( CLEAN_MULTIQC_PATHS.out.collect() )   // Md5sum for _data.zip files from cleaned multiqc paths 
                .combine( CLEAN_OUTPUT_PATHS.out.collect() )   // Md5sum for cleaned alignment report txt files
                .map { true }
    
    GENERATE_MD5SUMS( assay_suffix, ch_root_directory, ch_ready )

    VALIDATE_PROCESSING( ch_meta, runsheet_path, ch_root_directory, glds_accession, ch_ready )

    ch_raw_multiqc_zip = CLEAN_MULTIQC_PATHS.out
                            .flatten()
                            .filter { zip -> zip.name.contains('raw_multiqc') }

    UPDATE_ASSAY_TABLE( ch_meta, runsheet_path, ch_isa_zip, ch_raw_multiqc_zip )

    publish:
    processing_info = ch_root_directory.combine(PACKAGE_PROCESSING_INFO.out.zip)
    raw_protocol = ch_root_directory.combine(ch_raw_protocol)
    readme = ch_root_directory.combine(GENERATE_README.out.readme)
    raw_md5sum = ch_root_directory.combine(GENERATE_MD5SUMS.out.raw_md5sum) // emits nothing if raw md5sums file wasn't generated
    processed_md5sum = ch_root_directory.combine(GENERATE_MD5SUMS.out.processed_md5sum)
    validation_log = ch_root_directory.combine(VALIDATE_PROCESSING.out.validation_log)
    assay_table = ch_root_directory.combine(UPDATE_ASSAY_TABLE.out.assay_table)
}

output {
    processing_info { 
        path {outdir, files -> files >> "${outdir}/GeneLab/"} 
    }

    raw_protocol {
        path {outdir, files -> files >> "${outdir}/GeneLab/"} 
    }

    readme {
        path {outdir, files -> files >> "${outdir}/GeneLab/"} 
    }

    raw_md5sum {
        path {outdir, files -> files >> "${outdir}/GeneLab/"} 
    }

    processed_md5sum {
        path {outdir, files -> files >> "${outdir}/GeneLab/"} 
    }

    validation_log {
        path {outdir, files -> files >> "${outdir}/GeneLab/"} 
    }

    assay_table {
        path {outdir, files -> files >> "${outdir}/GeneLab/updated_curation_tables/"}
    }

}