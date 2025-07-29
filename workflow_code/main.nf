#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

def colorCodes = [
    c_line: "┅" * 70,
    c_back_bright_red: "\u001b[41;1m",
    c_bright_green: "\u001b[32;1m",
    c_blue: "\033[0;34m",
    c_yellow: "\u001b[33;1m",
    c_reset: "\033[0m"
]

// making sure specific expected nextflow version is being used
if( ! nextflow.version.matches( workflow.manifest.nextflowVersion ) ) {
    println "\n    This workflow requires Nextflow version $workflow.manifest.nextflowVersion, but version $nextflow.version is currently active."
    println "\n    You can set the proper version for this terminal session by running: `export NXF_VER=$workflow.manifest.nextflowVersion`"

    println "\n  Exiting for now.\n"
    exit 1
}


////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+

if (params.help) {
    log.info "\n┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅"
    log.info "┇          GeneLab MethylSeq Workflow: $workflow.manifest.version            ┇"
    log.info "┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅\n"
    
    log.info "                        HELP MENU\n"

    log.info "General settings can be configured in the 'nextflow.config' file.\n"

    log.info "  Usage example 1: Processing a GLDS dataset\n"
    log.info "    `nextflow run main.nf --gldsAccession GLDS-397`\n"

    log.info "  Usage example 2: Processing a local dataset (requires user-created runsheet)\n"
    log.info "    `nextflow run main.nf --runsheet my-runsheet.csv`\n\n"

    exit 0
}

////////////////////////////////////////////////////
/* --              DEBUG WARNING               -- */
////////////////////////////////////////////////////
if ( params.truncate_to ) {

    println( "\n    WARNING WARNING: DEBUG OPTIONS ARE ENABLED!\n" )

    params.truncate_to ? println( "        - Truncating reads to first ${ params.truncate_to } records" ) : null

    println ""

}


////////////////////////////////////////////////////
/* --                WORKFLOW                  -- */
////////////////////////////////////////////////////

include { FETCH_ISA } from './modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from './modules/isa_to_runsheet.nf'

include { PARSE_RUNSHEET } from './workflows/parse_runsheet.nf'
include { STAGE_RAW_READS } from './workflows/stage_raw_reads.nf'
include { REFERENCES } from './workflows/references.nf'

include { FASTQC as RAW_FASTQC } from './modules/fastqc.nf' addParams( fastqc_publish_dir: "00-RawData" )
include { FASTQC as TRIMMED_FASTQC } from './modules/fastqc.nf' addParams( fastqc_publish_dir: "01-TrimFilter" )
include { TRIMGALORE } from './modules/trimgalore.nf'
include { NUGEN } from './modules/nugen.nf'

include { BUILD_BISMARK; ALIGN_BISMARK; DEDUPE; EXTRACT_CALLS; BISMARK_REPORT; BISMARK_SUMMARY } from './modules/bismark.nf'
include { QUALIMAP } from './modules/qualimap.nf'
include { DIFF_METHYLATION } from './modules/methylation.nf'

include { SAMTOOLS_SORT as SORT_BAM } from './modules/samtools.nf' addParams( suffix: ".bam" )
include { SAMTOOLS_SORT as SORT_DEDUPED_BAM } from './modules/samtools.nf' addParams( suffix: ".deduplicated.bam" )

include { TO_PRED; TO_BED } from './modules/ucsc.nf'
include { GENES_TO_TRANSCRIPTS } from './modules/genes_to_transcripts.nf'

include { MULTIQC as RAW_MULTIQC } from './modules/multiqc.nf' addParams(MQCLabel:"raw", multiqc_publish_dir: "00-RawData/FastQC_Reports/raw_multiqc_report")
include { MULTIQC as TRIMMED_MULTIQC } from './modules/multiqc.nf' addParams(MQCLabel:"trimmed", multiqc_publish_dir: "01-TrimFilter/Trimmed_Multiqc_Report")
include { MULTIQC as ALIGN_MULTIQC } from './modules/multiqc.nf' addParams(MQCLabel:"align", multiqc_publish_dir: "04-Bismark_Reports")

ch_dp_tools_plugin = params.dp_tools_plugin ? Channel.value(file(params.dp_tools_plugin)) : Channel.value(file("$projectDir/bin/dp_tools__methylseq_dna"))
ch_runsheet = params.runsheet_path ? Channel.fromPath(params.runsheet_path) : null
ch_isa_archive = params.isa_archive_path ? Channel.fromPath(params.isa_archive_path) : null
ch_multiqc_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.fromPath("NO_FILE")

workflow {

    // If no runsheet path, fetch ISA archive from OSDR (if needed) and convert to runsheet
    if (ch_runsheet == null) {
        if (ch_isa_archive == null) {
            FETCH_ISA()
            ch_isa_archive = FETCH_ISA.out.isa_archive
        }
        ISA_TO_RUNSHEET( ch_isa_archive, ch_dp_tools_plugin )
        ch_runsheet = ISA_TO_RUNSHEET.out.runsheet
    }

    // Get sample metadata from runsheet
    PARSE_RUNSHEET( ch_runsheet )

    samples = PARSE_RUNSHEET.out.samples

    // Extract metadata from the first sample and set it as a channel
    samples | first
            | map { meta, reads -> meta }
            | set { ch_meta }

    STAGE_RAW_READS( samples )

    STAGE_RAW_READS.out.raw_reads | RAW_FASTQC

    STAGE_RAW_READS.out.raw_reads | TRIMGALORE

    if (params.nugen) {
        NUGEN(TRIMGALORE.out.reads, "${ projectDir }/bin/trimRRBSdiversityAdaptCustomers.py")
        NUGEN.out.reads | set { ch_trimmed_reads }
    } else {
        TRIMGALORE.out.reads | set { ch_trimmed_reads }
    }

    ch_trimmed_reads | TRIMMED_FASTQC

    samples | map { it[0].id }
        | collectFile(name: "samples.txt", sort: true, newLine: true, storeDir: "${params.outdir}/${params.gldsAccession}/processing_scripts")
        | set { ch_samples_txt }

    RAW_FASTQC.out.zip | map { it -> [ it[1] ] } | flatten | collect | set { raw_mqc_ch }
    RAW_MULTIQC( ch_samples_txt, raw_mqc_ch, ch_multiqc_config )

    TRIMMED_FASTQC.out.zip | map { it -> [ it[1] ] } | flatten | concat( TRIMGALORE.out.reports ) | collect | set { trimmed_mqc_ch }
    TRIMMED_MULTIQC( ch_samples_txt, trimmed_mqc_ch, ch_multiqc_config )

    REFERENCES( ch_meta | map { it.organism_sci } )
    REFERENCES.out.genome_annotations | set { genome_annotations }

    BUILD_BISMARK( 
      genome_annotations,
      ch_meta,
      REFERENCES.out.reference_version_and_source
    )

    // TO_PRED( 
    //   genome_annotations | map { it[1] }, 
    //   ch_meta,
    //   REFERENCES.out.reference_version_and_source
    // )

    // TO_BED( 
    //   TO_PRED.out, 
    //   ch_meta,
    //   REFERENCES.out.reference_version_and_source
    // )

    // GENES_TO_TRANSCRIPTS( 
    //   genome_annotations | map { it[1] }, 
    //   ch_meta,
    //   REFERENCES.out.reference_version_and_source
    // )

    ALIGN_BISMARK( ch_trimmed_reads, BUILD_BISMARK.out.build )

    // ALIGN_BISMARK.out.bam | SORT_BAM

    // QUALIMAP(genome_annotations | map { it[1] }, SORT_BAM.out.sorted_bam)

    // // Dedupe only if not RRBS ch_meta | map { it.rrbs }
    // if ( true ) {
    //     ALIGN_BISMARK.out.bam | set { ch_aligned_reads }
    //     ch_dedupe_report = Channel.fromPath("NO_FILE") | collect
    // } else {
    //     ALIGN_BISMARK.out.bam | DEDUPE
    //     DEDUPE.out.bam | set { ch_aligned_reads }
    //     ch_dedupe_report = DEDUPE.out.report

    //     DEDUPE.out.bam | SORT_DEDUPED_BAM
    // }

    // ch_aligned_reads | combine( BUILD_BISMARK.out.build ) | EXTRACT_CALLS

    // BISMARK_REPORT(
    //     ch_meta,
    //     ALIGN_BISMARK.out.report,
    //     ALIGN_BISMARK.out.stats,
    //     EXTRACT_CALLS.out.report,
    //     EXTRACT_CALLS.out.bias,
    //     ch_dedupe_report
    // )

    // BISMARK_SUMMARY(
    //     ALIGN_BISMARK.out.bam | map { it -> [ it[1] ] } | collect,
    //     ALIGN_BISMARK.out.report | collect,
    //     ALIGN_BISMARK.out.stats | collect,
    //     EXTRACT_CALLS.out.report | collect,
    //     EXTRACT_CALLS.out.bias | collect,
    //     ch_dedupe_report | collect
    // )

    // ALIGN_BISMARK.out.report | concat( ALIGN_BISMARK.out.stats ) 
    //     | concat( EXTRACT_CALLS.out.report ) 
    //     | concat( EXTRACT_CALLS.out.bias )
    //     | concat( QUALIMAP.out )
    //     | collect 
    //     | set { align_mqc_ch }

    // ALIGN_MULTIQC( ch_samples_txt, align_mqc_ch, ch_multiqc_config )

    // DIFF_METHYLATION(
    //     "${ projectDir }/bin/differential_methylation.R",
    //     REFERENCES.out.simple_organism_name,
    //     EXTRACT_CALLS.out.cov | collect,
    //     ch_runsheet,
    //     BUILD_BISMARK.out.build,
    //     params.reference_table,
    //     REFERENCES.out.gene_annotations,
    //     REFERENCES.out.reference_version_and_source
    // )
}
