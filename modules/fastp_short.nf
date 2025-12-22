// Process: Fastp for adapter trimming and quality control
process FASTP_SHORT {
    container 'https://depot.galaxyproject.org/singularity/fastp:0.24.0--heae3180_1'
    publishDir "${params.output_dir}/short_read/02_preprocessing/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.json", emit: json_report
    path "${sample_id}_fastp.html", emit: html_report

    script:
    def single_end = reads instanceof Path
    // checks if the "reads" input is a single Path object
    // (meaning single end) or multiple (meaning paired end)
    // the result (true or false) is assigned to the 
    // "single_end" variable

    // this allows one process to check run fastp on either
    // single end or paired end reads by checking the number
    // of inputs automatically

    if (single_end) {
        """
        fastp \\
            -i ${reads} \\
            -o ${sample_id}_trimmed_R1.fastq.gz \\
            -j ${sample_id}_fastp.json \\
            -h ${sample_id}_fastp.html \\
            -q 20 \\
            --trim_poly_x \\
            --length_required 50
        """
    //  automatic adapter trimming on single end data

    } else {
        """
        fastp \\
            -i ${reads[0]} \\
            -I ${reads[1]} \\
            -o ${sample_id}_trimmed_R1.fastq.gz \\
            -O ${sample_id}_trimmed_R2.fastq.gz \\
            -j ${sample_id}_fastp.json \\
            -h ${sample_id}_fastp.html \\
            -q 20 \\
            --trim_poly_x \\
            --length_required 50 \\
            --detect_adapter_for_pe \\
            --correction
        """
    }
    // For paired-end reads

        // -i: read1 input

        // -I: read2 input

        // -o: read1 output

        // -O: read2 output

        // -j: json report file

        // -h: html report file

        // -q: sets minimum phred score (20) for base calling; bases below this score will be trimmed

        // --trim_poly_x: enables poly base trimming, default value is 10

        // --length_required: minimum length a read must have after trimming to be kept

        // --detect_adapter_for_pe: automatic detection of adapters for paired end data

        // --correction: enable base correction in overlapped regions (only for PE data)
        //              performs overlap analysis for PE data, trying to find an overlap of each pair of reads.
        //              If a proper overlap is found, it can correct mismatched base pairs in overlapped regions
        //              if one base is high quality while the other is with low quality.
        //              If a base is corrected, the quality of its paired base will be assigned to it.
}
