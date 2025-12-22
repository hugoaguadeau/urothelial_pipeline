// Process: FastQC
process FASTQC_SHORT {
    container 'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.output_dir}/short_read/01_quality_control/raw/fastqc", mode: 'copy'

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
    
    // export _JAVA_OPTIONS="-Xmx${task.memory.toGiga()}g" gives
    // the amount of memory to use. This seems to be needed for
    // this process as it fails without it. I'm not sure why.

    // -o: output directory. If this option is not set then the
    // output file for each sequence file is created in the same
    
    // .: directory as the sequence file which was processed.

    // $reads: the input file/s
}
