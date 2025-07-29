library(optparse)
## other packages loaded below after help menu might be called

############################################
############ handling arguments ############
############################################

parser <- OptionParser(description = "\n  This is really only intented to be called from within the GeneLab MethylSeq workflow.")

parser <- add_option(parser, c("-v", "--verbose"),
                     action = "store_true",
                     default = FALSE, help = "Print extra output")

parser <- add_option(parser, c("--bismark_methylation_calls_dir"),
                     default = "Bismark_Methylation_calls",
                     help = "Directory holding *bismark.cov.gz files")

parser <- add_option(parser, c("--path_to_runsheet"), help = "Path to the runsheet")

parser <- add_option(parser, c("--simple_org_name"), 
                     help = "Simple organism name (must match 1st column of *annotations.csv ref table)")

parser <- add_option(parser, c("--ref_dir"),
                     default = "Reference_Genome_Files",
                     help = "Directory holding reference genome files (e.g. *.bed and *.gtf files)")

parser <- add_option(parser, c("--methylkit_output_dir"),
                     default = "MethylKit_Outputs",
                     help = "Directory for methylkit output files")

parser <- add_option(parser, c("--limit_samples_to"),
                     default = "all",
                     help = "Limits the number of samples being processed (won't do real factor comparisons if set)")

parser <- add_option(parser, c("--ref_org_table_link"),
                     help = "GeneLab reference annotations.csv file")

parser <- add_option(parser, c("--ref_annotations_tab_link"),
                     help = "Link to reference-genome annotations table")

parser <- add_option(parser, c("--methRead_mincov"), default = 10, type = "integer",
                     help = "Passed to mincov argument of methRead() call")

parser <- add_option(parser, c("--mc_cores"), default = 4, type = "integer",
                     help = "Passed to mc.cores argument of calculateDiffMeth() call")

parser <- add_option(parser, c("--getMethylDiff_difference"), default = 25, type = "integer",
                     help = "Passed to difference argument of getMethylDiff() call")

parser <- add_option(parser, c("--getMethylDiff_qvalue"), default = 0.01, type = "double",
                     help = "Passed to qvalue argument of getMethylDiff() call")

parser <- add_option(parser, c("--primary_keytype"),
                     help = "The keytype to use for mapping annotations (usually 'ENSEMBL' for most things; 'TAIR' for plants)")

parser <- add_option(parser, c("--test"),
                     action = "store_true",
                     default = FALSE,
                     help = "Provide solely this flag to run with a small test dataset that will be downloaded")

parser <- add_option(parser, c("--file_suffix"), 
		     default = "_GLMethylSeq",
		     help = "File suffix to add to output files")

myargs <- parse_args(parser)


############################################
########### for testing purposes ###########
############################################


if ( myargs$test ) {
    
    test_meth_data_link <- "https://figshare.com/ndownloader/files/39246758"
    test_meth_data_tarball <- "MethylSeq-test-meth-call-cov-files.tar"
    test_ref_data_link <- "https://figshare.com/ndownloader/files/38616860"
    test_ref_data_tarball <- "MethylSeq-test-ref-files.tar"
    
    cat("\n  Running test data from:\n\n")
    cat("    - test methylation calls:   ", test_meth_data_link, "\n", sep = "")
    cat("    - test ref files:           ", test_ref_data_link, "\n", sep = "")
    
    cat("\n  NOTICE\n  Any other parameters are ignored when '--test' is provided.\n\n")

    suppressWarnings(suppressMessages(library(curl)))
    
    curl_download(url = test_meth_data_link, destfile = test_meth_data_tarball, quiet = TRUE)
    curl_download(url = test_ref_data_link, destfile = test_ref_data_tarball, quiet = TRUE)
    
    untar(test_meth_data_tarball, restore_times = FALSE)
    untar(test_ref_data_tarball, restore_times = FALSE)
    
    file.remove(test_meth_data_tarball, test_ref_data_tarball)

    myargs$v <- TRUE
    myargs$simple_org_name <- "MOUSE"
    myargs$bismark_methylation_calls_dir <- "test-meth-calls"
    myargs$path_to_runsheet <- "test-meth-calls/test-runsheet.csv"
    myargs$ref_dir <- "test-ref-files"
    myargs$methylkit_output_dir <- "test-MethylKit_Outputs"
    myargs$ref_org_table_link <- "https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv"
    myargs$ref_annotations_tab_link <- "https://figshare.com/ndownloader/files/36597114"
    myargs$methRead_mincov <- 2
    myargs$getMethylDiff_difference <- 1
    myargs$getMethylDiff_qvalue <- 0.5
    myargs$primary_keytype <- "ENSEMBL"

}

if ( myargs$v ) { cat("\n  NOTICE\n  Verbose logging has been specified.\n\n") }

############################################


############################################
########### checking on arguments ##########
############################################

# checking required arguments were set
required_args <- c("path_to_runsheet" = "--path_to_runsheet",
                   "ref_org_table_link" = "--ref_org_table_link",
                   "ref_annotations_tab_link" = "--ref_annotations_tab_link",
                   "primary_keytype" = "--primary_keytype")

for ( arg in names(required_args) ) {

    tryCatch( { get(arg, myargs) }, error = function(e) { 
              
        error_message = paste0("\n  The '", as.character(required_args[arg]),
                            "' argument must be provided.\nCannot proceed.\n\n")
        
        stop(error_message, call. = FALSE)

    })
    
}

# checking primary keytype is what's expected
currently_accepted_keytypes <- c("ENSEMBL", "TAIR")
if ( ! myargs$primary_keytype %in% currently_accepted_keytypes ) {
    
    error_message = paste0("\n  The current potential --primary_keytypes are: ", paste(currently_accepted_keytypes, collapse = ", "),
                           "\nCannot proceed with: ", myargs$primary_keytype, "\n\n")
    
    stop(error_message, call. = FALSE)
    
}

# checking runsheet file exists
if ( ! file.exists(myargs$path_to_runsheet) ) {
    
    stop("\n  The specified --path_to_runsheet does not seem to point to an actual file.\nCannot proceed.\n\n", call. = FALSE)

}

############################################
############# loading packages #############
############################################

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(methylKit)))
suppressWarnings(suppressMessages(library(genomation)))

############################################


############################################
############# helper functions #############
############################################
get_single_file_path <- function(target_dir, search_pattern) {
    
    hits <- list.files(target_dir, pattern = search_pattern, full.names = TRUE)
    
    # checking only one was found matching pattern search
    if ( length(hits) != 1 ) { 
        
        error_message <- paste0("\n  A single file was not found in the ", target_dir,
                             " directory based on the search pattern '",
                             search_pattern, "'.\nCannot proceed.\n\n")
        
        stop(error_message, call. = FALSE)
    }
    
    return(hits)
    
}

order_input_files <- function(sample_names, paths) {
    
    # this function takes in a vector of sample_names and
    # a vector of paths to coverage files, and returns a vector
    # of coverage files in the same order as the sample names
    
    ordered_paths <- c()
    
    for ( sample in sample_names ) {
        
        search_pattern <- paste0(sample, "_bismark")
        hits <- paths[grep(search_pattern, paths)]
        
        # making sure there is exactly one match
        if ( length(hits) != 1 ) {
            
            stop("\n  There was a problem matching up sample names with their coverage files.\nCannot proceed.\n\n", call. = FALSE)
            
        }
        
        ordered_paths <- c(ordered_paths, hits)
        
    }
    
    return(ordered_paths)
    
}

############################################


############################################
############# setting things up ############
############################################

### finding reference bed file
ref_bed_path <- get_single_file_path(myargs$ref_dir, ".*.bed")


### finding reference gene-to-transcript mapping file
ref_gene_transcript_map_path <- get_single_file_path(myargs$ref_dir, ".*-gene-to-transcript-map.tsv")

### reading runsheet
runsheet <- read.csv(myargs$path_to_runsheet)

### getting all factors
# mock factors if wanted for testing more than 1
# runsheet$Factor.Value.Other <- c(rep("Mad", 10), rep("Dog", 6))
# runsheet$Factor.Value.Other2 <- c(rep("LilMac", 6), rep("Fighter", 10))
# runsheet$Factor.Value.Other2 <- c(rep("A", 6), rep("B", 8), rep("C", 2))
factors <- runsheet %>% dplyr::select(starts_with("Factor.Value"))
colnames(factors) <- paste("factor", 1:dim(factors)[2], sep = "_")

# checking if there is more than two unique values in a given factor, as methylkit isn't built for that
    # might not need this after doing combined only way as RNAseq does, need to check later
# for ( factor in colnames(factors)) {

#     curr_unique_entries <- factors[factor] %>% unique() %>% pull()

#     if ( length(curr_unique_entries) > 2 ) {

#         error_message <- paste0("\n  One of the factors has more than two entries:\n\n  ",
#                              paste(curr_unique_entries, collapse = ", "),
#                              "\n\n  Only pairwise comparisons can be done currently. Cannot proceed.\n\n")
        
#         stop(error_message, call. = FALSE)

#     }

# }


# factor dataframe
factor_df <- data.frame(sample_id = runsheet %>% pull("Sample.Name"), factors)

# subsetting runsheet if specified
if ( myargs$limit_samples_to != "all" ) {
    
    myargs$limit_samples_to <- as.integer(myargs$limit_samples_to)

    runsheet <- runsheet[ 1:myargs$limit_samples_to, ]
    
    # making up factors for testing if things were subset (so that we have multiple factor types even if not really in the subset samples)
    len_type_A <- ceiling(myargs$limit_samples_to / 2)
    len_type_B <- myargs$limit_samples_to - len_type_A
    
    factor_1 <- c(rep("A", len_type_A), rep("B", len_type_B))
    
    factor_df <- data.frame(sample_id = runsheet %>% pull("Sample.Name"), factor_1)

}

if ( myargs$v ) { cat("\nFactor df:\n") ; print(factor_df) ; cat("\n") }

# storing just sample names in a vector
sample_names <- runsheet$Sample.Name

# getting list of *.bismark.cov.gz files
bismark_cov_paths <- list.files(myargs$bismark_methylation_calls_dir, pattern = ".*.bismark.cov.gz", full.names = TRUE, recursive = TRUE)

# making sure files list matches length of runsheet
if ( dim(runsheet)[1] != length(bismark_cov_paths) ) {

    error_message <- paste0("\n  The number of '*.bismark.cov.gz' files found in the ", 
                         myargs$bismark_methylation_calls_dir, 
                         " directory\n  does not match the number of samples specified in the runsheet.\nCannot proceed.\n\n")
    
    stop(error_message, call. = FALSE)
    
}

if ( myargs$v ) { cat("\nBismark coverage files:\n") ; print(bismark_cov_paths) ; cat("\n") }

# making output directory
dir.create(myargs$methylkit_output_dir, showWarnings = FALSE)

# making sure coverage-file-paths vector is in the same order as the sample names
bismark_cov_paths <- order_input_files(sample_names, bismark_cov_paths)

# making table with filenames and bismark coverage file path
samples_and_covs_paths_df <- data.frame("sample_id" = sample_names, "coverage_file_path" = bismark_cov_paths)

## GeneLab combines all factors and just runs those contrasts, so setting that up here
# some of this initial structuring comes from Jonathan Oribello's work in the GeneLab RNAseq workflow, thanks Jonathan :)
study_df <- factor_df %>% column_to_rownames("sample_id")

write.table(study_df, file.path(myargs$methylkit_output_dir, paste0("SampleTable",myargs$file_suffix,".csv")), sep = ",", quote = FALSE)

# if there are multiple factors, here we are concatenating them to make one combined one
groups <- apply(study_df, 1, paste, collapse = " & ")
group_names <- paste0("(", groups, ")", sep = "") # human readable group names
safe_group_names <- make.names(groups) # group naming compatible with R models
names(safe_group_names) <- group_names

samples_and_combined_factors_df <- data.frame("sample_id" = sample_names, "combined_factor" = safe_group_names)

if ( myargs$v ) { cat("\nSamples and factors df:\n") ; print(samples_and_combined_factors_df) ; cat("\n") }


##### Format contrasts table, defining pairwise comparisons for all groups #####
contrasts <- combn(levels(factor(safe_group_names)), 2) # generate matrix of pairwise group combinations for comparison
contrast.names <- combn(levels(factor(names(safe_group_names))), 2)


## this way is only doing one-way contrasts
contrast.names <- paste(contrast.names[1,], contrast.names[2,], sep = "v")
colnames(contrasts) <- contrast.names

if ( myargs$v ) { cat("\nContrasts table:\n") ; print(contrasts) ; cat("\n") }

write.table(contrasts, file.path(myargs$methylkit_output_dir, paste0("contrasts",myargs$file_suffix,".csv")), sep = "\t", quote = FALSE, row.names = FALSE)

# making a single table with info needed for methylkit
sample_meth_info_df <- full_join(samples_and_covs_paths_df, samples_and_combined_factors_df, "sample_id")

if ( myargs$v ) { cat("\nPrimary methylkit df:\n") ; print(sample_meth_info_df) ; cat("\n") }

# reading in transcript features
gene.obj <- readTranscriptFeatures(ref_bed_path, up.flank = 1000, 
                                   down.flank = 1000, remove.unusual = TRUE, unique.prom = TRUE)

# reading in gene to transcript mapping file
gene_transcript_map <- 
    read.table(ref_gene_transcript_map_path, sep = "\t", col.names = c("gene_ID", "feature.name"))

# reading in functional annotation table
options(timeout = 600)

functional_annots_tab <- 
    read.table(myargs$ref_annotations_tab_link, sep = "\t", quote = "", header = TRUE)

# reading in reference annotations.csv table
ref_table <- read.csv(myargs$ref_org_table_link)

# at this point just to get a string for the 'assembly' argument of methRead() function
# setting variable with organism and ensembl version number
ensembl_version <- ref_table %>% 
    filter(name == myargs$simple_org_name) %>% pull(ensemblVersion)

org_and_ensembl_version <- paste(myargs$simple_org_name, ensembl_version, sep = "_")

############################################


############################################
########## moving onto methylkit ##########
############################################


### looping through contrasts and running methylkit and creating outputs for each ###
for ( i in 1:dim(contrasts)[2]) { 

    # getting current contrast
    curr_contrasts_vec <- contrasts[, i]
    
    # prefix for output files
    curr_output_prefix <- paste(curr_contrasts_vec, collapse = "_vs_")
    
    # making directory for output files of current contrast
    curr_output_dir <- file.path(myargs$methylkit_output_dir, curr_output_prefix)
    dir.create(curr_output_dir, showWarnings = FALSE)
    
    # getting subset sample info table relevant to current contrast
    curr_sample_info_df <- sample_meth_info_df %>% filter(combined_factor %in% curr_contrasts_vec)

    # getting which samples are relevant to current contrast
    curr_samples_vec <- curr_sample_info_df %>% pull(sample_id)
    
    # getting which files are relevant to current contrast
    curr_coverage_file_paths <- curr_sample_info_df %>% pull(coverage_file_path)
    
    # making binary vector for treatment argument to methRead()
    curr_treatment_vec <- c()
    for ( value in curr_sample_info_df$combined_factor ) {
        
        if ( value == curr_sample_info_df$combined_factor[1] ) {
            
            curr_treatment_vec <- c(curr_treatment_vec, 1)
            
        } else {
            
            curr_treatment_vec <- c(curr_treatment_vec, 0)
            
        }

    }
    

    ### setting up methylkit object
    curr_obj <- methRead(location = as.list(curr_coverage_file_paths),
                         sample.id = as.list(curr_samples_vec),
                         treatment = curr_treatment_vec,
                         pipeline = "bismarkCoverage",
                         assembly = org_and_ensembl_version,
                         header = FALSE,
                         mincov = myargs$methRead_mincov)
    
    ### Individual-base analysis
    # merging samples
    curr_meth <- methylKit::unite(curr_obj)
    
    # calculating differential methylation
    curr_myDiff <- calculateDiffMeth(curr_meth, mc.cores = myargs$mc_cores)
    

    # getting all sig differentially methylated bases (if any)
    curr_myDiff.all_sig <- getMethylDiff(curr_myDiff, difference = myargs$getMethylDiff_difference,
                                         qvalue = myargs$getMethylDiff_qvalue)
    
    # getting just sig hyper-methlated bases (if any)
    curr_myDiff.hyper <- getMethylDiff(curr_myDiff, difference = myargs$getMethylDiff_difference, 
                                       qvalue = myargs$getMethylDiff_qvalue, type = "hyper")

    # getting just sig hypo-methylated bases (if any)
    curr_myDiff.hypo <- getMethylDiff(curr_myDiff, difference = myargs$getMethylDiff_difference,
                                      qvalue = myargs$getMethylDiff_qvalue, type = "hypo")
    
                                    
    # checking that any sig were recovered
    if ( dim(curr_myDiff.all_sig)[1] > 0 ) any_sig <- TRUE else any_sig <- FALSE
    if ( dim(curr_myDiff.hyper)[1] > 0 ) any_sig_hyper <- TRUE else any_sig_hyper <- FALSE
    if ( dim(curr_myDiff.hypo)[1] > 0 ) any_sig_hypo <- TRUE else any_sig_hypo <- FALSE
    
    ### Tile analysis ###
    # tiling
    curr_tiles_obj <- tileMethylCounts(curr_obj, win.size = 1000, step.size = 1000, cov.bases = 10)
    
    # merging tiled samples
    curr_tiles_meth <- methylKit::unite(curr_tiles_obj)
    
    # calculating differential methylation on tiles
    curr_tiles_diff <- calculateDiffMeth(curr_tiles_meth, mc.cores = myargs$mc_cores)

    # getting all sig differentially methylated tiles (if any)
    curr_tiles_myDiff.all_sig <- getMethylDiff(curr_tiles_diff, difference = myargs$getMethylDiff_difference, 
                                               qvalue = myargs$getMethylDiff_qvalue)
    
    # getting just sig hyper-methylated tiles (if any)
    curr_tiles_myDiff.hyper <- getMethylDiff(curr_tiles_diff, difference = myargs$getMethylDiff_difference, 
                                             qvalue = myargs$getMethylDiff_qvalue, type = "hyper")

    # getting just sig hypo-methylated tiles (if any)
    curr_tiles_myDiff.hypo <- getMethylDiff(curr_tiles_diff, difference = myargs$getMethylDiff_difference, 
                                            qvalue = myargs$getMethylDiff_qvalue, type = "hypo")

    # checking that any sig were recovered
    if ( dim(curr_tiles_myDiff.all_sig)[1] > 0 ) any_sig_tiles <- TRUE else any_sig_tiles <- FALSE
    if ( dim(curr_tiles_myDiff.hyper)[1] > 0 ) any_sig_tiles_hyper <- TRUE else any_sig_tiles_hyper <- FALSE
    if ( dim(curr_tiles_myDiff.hypo)[1] > 0 ) any_sig_tiles_hypo <- TRUE else any_sig_tiles_hypo <- FALSE
    
    ### Adding feature information ###
    
    ## adding features to individual-base objects and making df with features added to sig tables (if any sig)
    
    if ( any_sig ) { 
        
        curr_diffAnn <- annotateWithGeneParts(as(curr_myDiff.all_sig, "GRanges"), gene.obj)
        df_list <- list(tibble::rowid_to_column(data.frame(curr_myDiff.all_sig, row.names=NULL)), 
                        tibble::rowid_to_column(getAssociationWithTSS(curr_diffAnn)), 
                        tibble::rowid_to_column(as.data.frame(genomation::getMembers(curr_diffAnn))))
        curr_sig_all_bases_tab_with_features <- df_list %>% purrr::reduce(full_join, by="rowid") %>% 
                                                          .[, !names(.) %in% c("target.row", "rowid")]
        rm(df_list)    
    }
    
    if ( any_sig_hyper ) { 
        
        curr_diffAnn.hyper <- annotateWithGeneParts(as(curr_myDiff.hyper, "GRanges"), gene.obj)
        df_list <- list(tibble::rowid_to_column(data.frame(curr_myDiff.hyper, row.names=NULL)), 
                        tibble::rowid_to_column(getAssociationWithTSS(curr_diffAnn.hyper)), 
                        tibble::rowid_to_column(as.data.frame(genomation::getMembers(curr_diffAnn.hyper))))
        curr_sig_hyper_bases_tab_with_features <- df_list %>% purrr::reduce(full_join, by="rowid") %>%
                                                            .[, !names(.) %in% c("target.row", "rowid")]
        rm(df_list)
    }
    
    if ( any_sig_hypo ) {
        
        curr_diffAnn.hypo <- annotateWithGeneParts(as(curr_myDiff.hypo, "GRanges"), gene.obj)
        df_list <- list(tibble::rowid_to_column(data.frame(curr_myDiff.hypo, row.names=NULL)), 
                        tibble::rowid_to_column(getAssociationWithTSS(curr_diffAnn.hypo)), 
                        tibble::rowid_to_column(as.data.frame(genomation::getMembers(curr_diffAnn.hypo))))
        curr_sig_hypo_bases_tab_with_features <- df_list %>% purrr::reduce(full_join, by="rowid") %>%
                                                            .[, !names(.) %in% c("target.row", "rowid")]
        rm(df_list)
    }


    ## adding features to tiles objects and making df with features added to sig tables (if any sig)
    
    if ( any_sig_tiles ) {
        
        curr_tiles_diffAnn <- annotateWithGeneParts(as(curr_tiles_myDiff.all_sig, "GRanges"), gene.obj)
        df_list <- list(tibble::rowid_to_column(data.frame(curr_tiles_myDiff.all_sig, row.names=NULL)), 
                        tibble::rowid_to_column(getAssociationWithTSS(curr_tiles_diffAnn)), 
                        tibble::rowid_to_column(as.data.frame(genomation::getMembers(curr_tiles_diffAnn))))
        curr_tiles_sig_all_tab_with_features <- df_list %>% purrr::reduce(full_join, by="rowid") %>%
                                                            .[, !names(.) %in% c("target.row", "rowid")]
        rm(df_list)
    }
    
    
    if ( any_sig_tiles_hyper ) { 
        
        curr_tiles_diffAnn.hyper <- annotateWithGeneParts(as(curr_tiles_myDiff.hyper, "GRanges"), gene.obj)
        df_list <- list(tibble::rowid_to_column(data.frame(curr_tiles_myDiff.hyper, row.names=NULL)), 
                        tibble::rowid_to_column(getAssociationWithTSS(curr_tiles_diffAnn.hyper)), 
                        tibble::rowid_to_column(as.data.frame(genomation::getMembers(curr_tiles_diffAnn.hyper))))
        curr_tiles_sig_hyper_tab_with_features <- df_list %>% purrr::reduce(full_join, by="rowid") %>%
                                                            .[, !names(.) %in% c("target.row", "rowid")]
        rm(df_list)
    }

    if ( any_sig_tiles_hypo ) { 

        curr_tiles_diffAnn.hypo <- annotateWithGeneParts(as(curr_tiles_myDiff.hypo, "GRanges"), gene.obj)
        df_list <- list(tibble::rowid_to_column(data.frame(curr_tiles_myDiff.hypo, row.names=NULL)), 
                        tibble::rowid_to_column(getAssociationWithTSS(curr_tiles_diffAnn.hypo)), 
                        tibble::rowid_to_column(as.data.frame(genomation::getMembers(curr_tiles_diffAnn.hypo))))
        curr_tiles_sig_hypo_tab_with_features <- df_list %>% purrr::reduce(full_join, by="rowid") %>%
                                                            .[, !names(.) %in% c("target.row", "rowid")]
        rm(df_list)
    }

    
    ### Adding functional annotations ###
    
    ## for individual-base outputs (if any sig)
    # for each transcript ID in the sig_all_bases_tab_with_features table, getting 
    # its corresponding gene ID and adding that to the table, then adding full annotations
    
    if ( any_sig ) { 
            
        curr_sig_all_bases_tab_with_features_and_gene_IDs <- 
            left_join(curr_sig_all_bases_tab_with_features, gene_transcript_map)

        curr_sig_all_bases_tab_with_features_and_annots <- 
            left_join(curr_sig_all_bases_tab_with_features_and_gene_IDs, 
                      functional_annots_tab, by = c("gene_ID" = myargs$primary_keytype)) %>% arrange(qvalue)

        # renaming "gene_ID" column to be primary keytype
        curr_sig_all_bases_tab_with_features_and_annots <- 
            S4Vectors::rename(curr_sig_all_bases_tab_with_features_and_annots, "gene_ID" = myargs$primary_keytype)
        
    }
    
    if ( any_sig_hyper ) {

        curr_sig_hyper_bases_tab_with_features_and_gene_IDs <- 
            left_join(curr_sig_hyper_bases_tab_with_features, gene_transcript_map)

        curr_sig_hyper_bases_tab_with_features_and_annots <- 
            left_join(curr_sig_hyper_bases_tab_with_features_and_gene_IDs, 
                      functional_annots_tab, by = c("gene_ID" = myargs$primary_keytype)) %>% arrange(qvalue)

        # renaming "gene_ID" column to be primary keytype
        curr_sig_hyper_bases_tab_with_features_and_annots <- 
            S4Vectors::rename(curr_sig_hyper_bases_tab_with_features_and_annots, "gene_ID" = myargs$primary_keytype)

    }

    if ( any_sig_hypo ) {

        curr_sig_hypo_bases_tab_with_features_and_gene_IDs <- 
            left_join(curr_sig_hypo_bases_tab_with_features, gene_transcript_map)

        curr_sig_hypo_bases_tab_with_features_and_annots <- 
            left_join(curr_sig_hypo_bases_tab_with_features_and_gene_IDs, 
                      functional_annots_tab, by = c("gene_ID" = myargs$primary_keytype)) %>% arrange(qvalue)

        # renaming "gene_ID" column to be primary keytype
        curr_sig_hypo_bases_tab_with_features_and_annots <- 
            S4Vectors::rename(curr_sig_hypo_bases_tab_with_features_and_annots, "gene_ID" = myargs$primary_keytype)

    }

    # and writing out (if they didn't have any, still producing the file but it will say "None detected")
    none_detected_message <- "None detected."
    
    # individual bases, all sig (if any)
    curr_sig_all_bases_tab_with_features_and_annots_path <- 
        file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-methylated-bases.tsv"))
    
    if ( any_sig ) { 
        
        write.table(curr_sig_all_bases_tab_with_features_and_annots, curr_sig_all_bases_tab_with_features_and_annots_path, 
                    sep = "\t", quote = FALSE, row.names = FALSE)
    
    } else {
        
        curr_file <- file(curr_sig_all_bases_tab_with_features_and_annots_path)
        writeLines(none_detected_message, curr_file)
        close(curr_file)

        std_err_message <- paste0("\n  NOTICE\n  There were no significantly differentially methylated sites identified in the \n  '", 
                                  curr_output_prefix, "' contrast.\n  You will see an 'In max(i)...' warning about this.\n")
        
        write(std_err_message, stderr())
        
        
    }
    
    # individual bases, hyper sig (if any)
    curr_sig_hyper_bases_tab_with_features_and_annots_path <- 
        file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-hypermethylated-bases.tsv"))
    
    if ( any_sig_hyper ) { 

        write.table(curr_sig_hyper_bases_tab_with_features_and_annots, curr_sig_hyper_bases_tab_with_features_and_annots_path, 
                    sep = "\t", quote = FALSE, row.names = FALSE)
        
    } else {
        
        curr_file <- file(curr_sig_hyper_bases_tab_with_features_and_annots_path)
        writeLines(none_detected_message, curr_file)
        close(curr_file)
        
        std_err_message <- paste0("\n  NOTICE\n  There were no significantly differentially hypermethylated sites identified in the \n  '", 
                                  curr_output_prefix, "' contrast.\n  You will see an 'In max(i)...' warning about this.\n")
        
        write(std_err_message, stderr())
        
    }

    # individual bases, hypo sig (if any)
    curr_sig_hypo_bases_tab_with_features_and_annots_path <- 
        file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-hypomethylated-bases.tsv"))
    
    if ( any_sig_hypo ) { 

        write.table(curr_sig_hypo_bases_tab_with_features_and_annots, curr_sig_hypo_bases_tab_with_features_and_annots_path, 
                    sep = "\t", quote = FALSE, row.names = FALSE)

    } else {
        
        curr_file <- file(curr_sig_hypo_bases_tab_with_features_and_annots_path)
        writeLines(none_detected_message, curr_file)
        close(curr_file)
        
        std_err_message <- paste0("\n  NOTICE\n  There were no significantly differentially hypomethylated sites identified in the \n  '", 
                                  curr_output_prefix, "' contrast.\n  You will see an 'In max(i)...' warning about this.\n")
        
        write(std_err_message, stderr())
        
    }
    
    
    ## for tiles output (if any sig)
    # for each transcript ID in the tiles_sig_all_out_tab_with_features table, getting 
    # its corresponding gene ID and adding that to the table, then adding full annotations

    if ( any_sig_tiles ) {

        curr_sig_all_tiles_tab_with_features_and_gene_IDs <- 
            left_join(curr_tiles_sig_all_tab_with_features, gene_transcript_map)
        
        curr_sig_all_tiles_tab_with_features_and_annots <- 
            left_join(curr_sig_all_tiles_tab_with_features_and_gene_IDs, 
                      functional_annots_tab, by = c("gene_ID" = myargs$primary_keytype)) %>% arrange(qvalue)

        # renaming "gene_ID" column to be primary keytype
        curr_sig_all_tiles_tab_with_features_and_annots <- 
            S4Vectors::rename(curr_sig_all_tiles_tab_with_features_and_annots, "gene_ID" = myargs$primary_keytype)

    }

    if ( any_sig_tiles_hyper ) {

        curr_sig_hyper_tiles_tab_with_features_and_gene_IDs <- 
            left_join(curr_tiles_sig_hyper_tab_with_features, gene_transcript_map)

        curr_sig_hyper_tiles_tab_with_features_and_annots <- 
            left_join(curr_sig_hyper_tiles_tab_with_features_and_gene_IDs, 
                      functional_annots_tab, by = c("gene_ID" = myargs$primary_keytype)) %>% arrange(qvalue)

        # renaming "gene_ID" column to be primary keytype
        curr_sig_hyper_tiles_tab_with_features_and_annots <- 
            S4Vectors::rename(curr_sig_hyper_tiles_tab_with_features_and_annots, "gene_ID" = myargs$primary_keytype)
        
    }    
    
    if ( any_sig_tiles_hypo ) {

        curr_sig_hypo_tiles_tab_with_features_and_gene_IDs <- 
            left_join(curr_tiles_sig_hypo_tab_with_features, gene_transcript_map)
        
        curr_sig_hypo_tiles_tab_with_features_and_annots <- 
            left_join(curr_sig_hypo_tiles_tab_with_features_and_gene_IDs, 
                      functional_annots_tab, by = c("gene_ID" = myargs$primary_keytype)) %>% arrange(qvalue)

        # renaming "gene_ID" column to be primary keytype
        curr_sig_hypo_tiles_tab_with_features_and_annots <- 
            S4Vectors::rename(curr_sig_hypo_tiles_tab_with_features_and_annots, "gene_ID" = myargs$primary_keytype)
        
    }

    
    # and writing out (if they didn't have any, still producing the file but it will say "None detected")

    # tiles, all sig (if any)
    curr_sig_all_tiles_tab_with_features_and_annots_path <- 
        file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-methylated-tiles.tsv"))
    
    if ( any_sig_tiles ) {

        write.table(curr_sig_all_tiles_tab_with_features_and_annots, curr_sig_all_tiles_tab_with_features_and_annots_path, 
                    sep = "\t", quote = FALSE, row.names = FALSE)

    } else {
        
        curr_file <- file(curr_sig_all_tiles_tab_with_features_and_annots_path)
        writeLines(none_detected_message, curr_file)
        close(curr_file)
        
        std_err_message <- paste0("\n  NOTICE\n  There were no significantly differentially methylated tiles identified in the \n  '", 
                                  curr_output_prefix, "' contrast.\n  You will see an 'In max(i)...' warning about this.\n")
        
        write(std_err_message, stderr())
        
    }
    
    # tiles, hyper sig (if any)
    curr_sig_hyper_tiles_tab_with_features_and_annots_path <- 
        file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-hypermethylated-tiles.tsv"))

    if ( any_sig_tiles_hyper ) {

        write.table(curr_sig_hyper_tiles_tab_with_features_and_annots, curr_sig_hyper_tiles_tab_with_features_and_annots_path, 
                    sep = "\t", quote = FALSE, row.names = FALSE)

    } else {
        
        curr_file <- file(curr_sig_hyper_tiles_tab_with_features_and_annots_path)
        writeLines(none_detected_message, curr_file)
        close(curr_file)

        std_err_message <- paste0("\n  NOTICE\n  There were no significantly differentially hypermethylated tiles identified in the \n  '", 
                                  curr_output_prefix, "' contrast.\n  You will see an 'In max(i)...' warning about this.\n")
        
        write(std_err_message, stderr())
        
    }
    

    # tiles, hypo sig (if any)    
    curr_sig_hypo_tiles_tab_with_features_and_annots_path <- 
        file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-hypomethylated-tiles.tsv"))
    
    if ( any_sig_tiles_hypo ) {

        write.table(curr_sig_hypo_tiles_tab_with_features_and_annots, curr_sig_hypo_tiles_tab_with_features_and_annots_path, 
                    sep = "\t", quote = FALSE, row.names = FALSE)
        
    } else {
        
        curr_file <- file(curr_sig_hypo_tiles_tab_with_features_and_annots_path)
        writeLines(none_detected_message, curr_file)
        close(curr_file)
        
        std_err_message <- paste0("\n  NOTICE\n  There were no significantly differentially hypomethylated tiles identified in the \n  '", 
                                  curr_output_prefix, "' contrast.\n  You will see an 'In max(i)...' warning about this.\n")
        
        write(std_err_message, stderr())
        
    }
    
    
    ### Making overview figure of percent diff. methylation across features (if any sig) ### (these pdfs won't exist if there aren't any sig different)
    
    # for bases
    curr_sig_diff_bases_across_features_plot_path <- file.path(curr_output_dir, 
                                                               paste0(curr_output_prefix, "-sig-diff-methylated-bases-across-features.pdf"))
    
    if ( any_sig ) { 
        
        pdf(curr_sig_diff_bases_across_features_plot_path)
        plotTargetAnnotation(curr_diffAnn, precedence = TRUE, main = "% of sig. diff. methylated sites across features")
        dev.off()

    } else {
        
        std_err_message <- paste0("\n  NOTICE\n  There were no significantly differentially methylated sites identified in the \n  '", 
                                  curr_output_prefix, "' contrast, so this pdf overview figure is not being produced:\n\n    ", 
                                  curr_sig_diff_bases_across_features_plot_path, "\n")
        
        write(std_err_message, stderr())

    }

    # for tiles    
    curr_sig_diff_tiles_across_features_plot_path <- file.path(curr_output_dir, 
                                                               paste0(curr_output_prefix, "-sig-diff-methylated-tiles-across-features.pdf"))
    
    if ( any_sig_tiles ) {
    
        pdf(curr_sig_diff_tiles_across_features_plot_path)
        plotTargetAnnotation(curr_tiles_diffAnn, precedence = TRUE, main = "% of sig. diff. methylated tiles across features")
        dev.off()
        
    } else {
        
        std_err_message <- paste0("\n  NOTICE\n  There were no significantly differentially methylated tiles identified in the \n  '", 
                                  curr_output_prefix, "' contrast, so this pdf overview figure is not being produced:\n\n    ", 
                                  curr_sig_diff_tiles_across_features_plot_path, "\n")
        
        write(std_err_message, stderr())
        
    }
    
}

### making and writing out a table of base-level percent methylated (treatment vector doesn't matter here, just making a mock one)
len_1s <- ceiling(length(sample_names) / 2)
len_0s <- length(sample_names) - len_1s
mock_treatment_vec <- c(rep(1, len_1s), rep(0, len_0s))

#setup grouping for output
grouped_sample_meth_info_df <- group_split(sample_meth_info_df %>% group_by(combined_factor))

obj <- methRead(location = as.list(sample_meth_info_df %>% pull(coverage_file_path)),
                sample.id = as.list(sample_meth_info_df %>% pull(sample_id)),
                treatment = mock_treatment_vec,
                pipeline = "bismarkCoverage",
                assembly = org_and_ensembl_version,
                header = FALSE,
                mincov = myargs$methRead_mincov)

meth <- methylKit::unite(obj)

perc.meth <- percMethylation(meth, rowids = TRUE)
perc.meth <- perc.meth %>% data.frame(check.names = FALSE) %>% rownames_to_column("location")

# writing out
for (i in 1:length(grouped_sample_meth_info_df)) {
    curr_group <- grouped_sample_meth_info_df[[i]]$combined_factor[1]
    perc.meth_path <- file.path(myargs$methylkit_output_dir, paste0(curr_group,"_base-level-percent-methylated.tsv"))
    write.table(perc.meth %>% dplyr::select(grouped_sample_meth_info_df[[i]]$sample_id), perc.meth_path,
                sep="\t", quote = FALSE, row.names = FALSE)

}

### making and writing out a table of tile-level percent methylated (contrasts don't matter here)
tiles_obj <- tileMethylCounts(obj, win.size = 1000, step.size = 1000, cov.bases = 10)
tiles_meth <- methylKit::unite(tiles_obj)

tiles_perc.meth <- percMethylation(tiles_meth, rowids = TRUE)
tiles_perc.meth <- tiles_perc.meth %>% data.frame(check.names = FALSE) %>% rownames_to_column("location")

for (i in 1:length(grouped_sample_meth_info_df)) {
    curr_group <- grouped_sample_meth_info_df[[i]]$combined_factor[1]
    tiles_perc.meth_path <- file.path(myargs$methylkit_output_dir, paste0(curr_group,"_tile-level-percent-methylated.tsv"))
    write.table(tiles_perc.meth %>% dplyr::select(grouped_sample_meth_info_df[[i]]$sample_id), tiles_perc.meth_path,
                sep="\t", quote = FALSE, row.names = FALSE)
}
