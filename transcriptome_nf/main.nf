
nextflow.enable.dsl=2

// Workflow definition
workflow {
    // Header

    //Check for output directory
    ensure_output_dirs()

    sam_files = Channel.fromPath("${params.sam_dir}/*.sam")
    parquet_files = sam_to_parquet(sam_files)

}



// Process definitions
process ensure_output_dirs {
    """
    mkdir -p "${params.parquet_dir}"
    """
}


process sam_to_parquet {
    label 'local'

    publishDir "${params.parquet_dir}", mode: 'copy'

    input:
    path sam_file

    output:
    path "${sam_file.baseName}.parquet"

    script:
    """
    python ${baseDir}/scripts/process_sam.py \\
           "${sam_file}" \\
           "${params.primer_file}" \\
           "${sam_file.baseName}.parquet"
    """
}




































//process merge_parquets {
//    label 'local'
//    publishDir "${launchDir}", mode: 'copy'
//    
//    input:
//    // Collect all Parquet files into a single list
//    path parquet_files //"${params.parquet_dir}"
//
//   output:
//    path "${params.output_file}"
//    
//    script:
//    """
//
//    # Debug: see what's staged
//    echo "DEBUG: Checking files in this work dir..."
//    ls -lh
//
//    python ${baseDir}/scripts/merge_parquets.py \\
//           ${parquet_files.join(' ')} \\
//           ${params.output_file}
//    
//   """
//}

//process basecall_dorado {
//    tag "${pod5_file}"
//
//    label 'gpu'
//
//    input:
//    path pod5_file
//
//    output:
//    path "${params.output_dir}/processed_${pod5_file.simpleName}.processed"
//
//    script:
//    """
//    # Run Dorado basecaller
//    dorado basecaller --config ${params.dorado_config} --device cuda \\
//                      --output sam ${pod5_file} > ${params.output_dir}/output_${pod5_file.simpleName}.sam
//
//    # Run your Python script to manipulate the SAM file
//    python process_sam.py ${params.output_dir}/output_${pod5_file.simpleName}.sam \\
//                          ${params.output_dir}/processed_${pod5_file.simpleName}.processed
//
//    # Optional: Cleanup intermediate SAM file
//    rm ${params.output_dir}/output_${pod5_file.simpleName}.sam
//    """
//}

//process recombine_files {
//    tag 'recombine'
//
//    input:
//    path processed_files from basecall_dorado.out.collect()  // Collect all processed files
//
//    output:
//    path "${params.output_dir}/recombined_output.txt"
//
//    script:
//    """
//    # Concatenate all processed files into one
//    cat ${processed_files} > ${params.output_dir}/recombined_output.txt
//    """
//}

