// Generate STAR index
process STAR_INDEX_SHORT {
    container 'https://depot.galaxyproject.org/singularity/star:2.7.11b--h5ca1c30_4'
    publishDir "${params.output_dir}/short_read/03_alignment/star_index", mode: 'copy'

    input:
    path reference_genome
    path reference_gtf

    output:
    path "star_index", emit: index

    script:
    """
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${reference_genome} \\
        --sjdbGTFfile ${reference_gtf} \\
        --runThreadN ${task.cpus} \\
        --sjdbOverhang 150
    """
    
    // runMode genomeGenerate: tells STAR to make an index

    // genomdir: output directory for the generated genome index files

    // genomeFastaFiles: path to the reference genome 

    // sjdbGTFfile: allows for the construction of a splice site database and allows for annotation

    // sjbdOverhang: specifies the length of the genomic
    // sequence around the annotated junction to be used in
    // constructing the splice junctions database. Ideally,
    // this length should be equal to the ReadLength-1.
    // In case of reads of varying length, the ideal value
    // is max(ReadLength)-1.
}
