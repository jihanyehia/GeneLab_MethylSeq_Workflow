process PACKAGE_PROCESSING_INFO { 
    tag "Purging paths of log files and creating zipped processing info"
    label "zip"

    input:
        path(processing_scripts)
        val(assay_suffix)

    output:
        path("processing_info${assay_suffix}.zip"), emit: zip

    script:
    """
    mkdir -p processing_scripts

    for f in nextflow*.txt; do
        echo "Purging file paths from \$f"
        clean_paths.sh "\$f"

        # Add assay suffix to the end of nextflow log files before .txt extension, if not already present
        if [[ "\$f" != *"${assay_suffix}.txt" ]]; then
            mv "\$f" "\${f%.txt}${assay_suffix}.txt"
        fi
    done

    mv nextflow*.txt processing_scripts/ 2>/dev/null || true
        
    # Zip 
    zip -r processing_info${assay_suffix}.zip processing_scripts
    """
}