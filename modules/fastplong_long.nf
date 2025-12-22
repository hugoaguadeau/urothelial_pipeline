// Process: Fastplong for adapter trimming with reports (similar to fastp)
process FASTPLONG_LONG {
    container 'https://depot.galaxyproject.org/singularity/fastplong:0.2.2--heae3180_0'
    publishDir "${params.output_dir}/long_read/02_preprocessing/fastplong", mode: 'copy'

    input:
    tuple val(sample_id), path(input_fastq)

    output:
    tuple val(sample_id), path("${sample_id}-fastplong-output_reads.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastplong.json", emit: json_report
    path "${sample_id}_fastplong.html", emit: html_report

    script:
    """
    fastplong \\
        -i ${input_fastq} \\
        -o ${sample_id}-fastplong-output_reads.fastq.gz \\
        -q 10 \\
        -u 40 \\
        --trim_poly_x \\
        -j ${sample_id}_fastplong.json \\
        -h ${sample_id}_fastplong.html \\
        -R "${sample_id} fastplong report"
    """

    // Key parameters:
    // -i: input file name
    // -o: output file name
    // -q: the quality value that a base is qualified (10 means phred quality >=Q10 is qualified)
    // -u: how many percent of bases are allowed to be unqualified (40 means 40%)
    // --trim_poly_x: enable polyX trimming in 3' ends
    // -j: JSON report output file
    // -h: HTML report output file  
    // -R: report title
}
