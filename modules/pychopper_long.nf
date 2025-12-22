// Process: Pychopper for full-length cDNA read identification
process PYCHOPPER_LONG {
    container 'https://depot.galaxyproject.org/singularity/pychopper:2.7.10--pyhdfd78af_0'
    publishDir "${params.output_dir}/long_read/02_preprocessing/pychopper", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    tuple val(sample_id), path("${sample_id}_full_length.fastq"), emit: full_length_reads
    tuple val(sample_id), path("${sample_id}_unclassified.fastq"), emit: unclassified_reads
    tuple val(sample_id), path("${sample_id}_rescued.fastq"), emit: rescued_reads
    path("${sample_id}_report.pdf"), emit: report
    path("${sample_id}_stats.tsv"), emit: stats

    script:
    """
    pychopper \
        -t ${task.cpus} \
        -r ${sample_id}_report.pdf \
        -S ${sample_id}_stats.tsv \
        -u ${sample_id}_unclassified.fastq \
        -w ${sample_id}_rescued.fastq \
        ${fastq_file} \
        ${sample_id}_full_length.fastq
    """
     
    // flair says to use pychopper for cDNA preprocessing

    // t: threads
   
    // r: Report PDF (pychopper_report.pdf)

    // S: Write statistics to this file

    // u: Write unclassified reads to this file

    // w: Write rescued reads to this file

    // penultimate line: input file

    // final line: output file

}
