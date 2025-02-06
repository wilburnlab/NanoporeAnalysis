
nextflow.enable.dsl=2

// Workflow definition
workflow {
    // Header

    batch_results = generate_batches("${params.input}")
    read_batches = batch_results.batches | flatten
    processing_results = process_reads(read_batches)
    protein_results = extract_proteins(processing_results.collect())

}

process generate_batches {
    //label 'local'
    
    //scratch true // Write batches to $TMPDIR
    //stageInMode 'copy'

    input:
    path sam_file // can be individual sam or directory of sam files

    output:
    path "${params.read_batch_dir}/*.sam", emit: batches

    script:
    """
    mkdir -p "${params.read_batch_dir}"
    python ${baseDir}/scripts/batch_sam.py \\
           "${sam_file}" \\
           "${params.read_batch_size}" \\
           "${params.read_batch_dir}"
    """

}


process process_reads {
    //label 'local'

    publishDir "${params.proc_read_dir}", mode: 'copy'

    input:
    path read_batch

    output:
    path "${read_batch.baseName}.pqt"

    script:
    """
    mkdir -p "${params.proc_read_dir}"
    python ${baseDir}/scripts/process_sam.py \\
           "${read_batch}" \\
           "${params.primer_file}" \\
           "${read_batch.baseName}.pqt"
    """
}



process extract_proteins {

    publishDir "${launchDir}", mode: 'copy'

    input:
    path pqt_paths

    output:
    path "${params.out_label}_proteins.fasta"
    path "${params.out_label}_protein-counts.tab"

    script:
    """
    python ${baseDir}/scripts/extract_proteins.py \\
           "${params.out_label}" \\
           ${pqt_paths}
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

