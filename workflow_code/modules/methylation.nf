process DIFF_METHYLATION {
  publishDir "${params.outdir}/${params.gldsAccession}/05-MethylKit_DM",
    mode: params.publish_dir_mode

  label 'med_mem'

  input:
    path("differential_methylation.R")
    val(simple_organism_name)
    path("bismark_in/*")
    path(runsheet)
    path(bismark_index_dir)
    val(annotations_csv_url_string)
    path(annotation_file)
    tuple val(ensemblVersion), val(ref_source)

  output:
    path("*_GLMethylSeq.csv")
    path("*-percent-methylated.tsv")
    path("*_vs_*/*")

  script:
    def primary_keytype = ref_source == 'ensembl_plants' ? 'TAIR' : 'ENSEMBL'
    """
    Rscript --vanilla differential_methylation.R \
        --bismark_methylation_calls_dir bismark_in \
        --path_to_runsheet ${runsheet} \
        --simple_org_name ${simple_organism_name} \
        --ref_dir ${bismark_index_dir} \
        --methylkit_output_dir . \
        --ref_org_table_link ${annotations_csv_url_string} \
        --ref_annotations_tab_link ${annotation_file} \
        --methRead_mincov ${params.methRead_mincov} \
        --getMethylDiff_difference ${params.MethylDiff_difference} \
        --getMethylDiff_qvalue ${params.MethylDiff_qvalue} \
        --primary_keytype ${primary_keytype} \
        --file_suffix '_GLMethylSeq'
    """
}
