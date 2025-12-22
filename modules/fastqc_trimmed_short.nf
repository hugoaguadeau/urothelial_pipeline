// Process: FastQC on Trimmed Reads
process FASTQC_TRIMMED_SHORT {
    container 'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.output_dir}/short_read/01_quality_control/trimmed/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.zip", emit: zip
    path "*.html", emit: html

    script:
    """
    export _JAVA_OPTIONS="-Xmx${task.memory.toGiga()}g"
    fastqc -o . $reads
    """
}
