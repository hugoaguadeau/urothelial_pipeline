// Process: StringTie Merge - Create unified transcript reference
process STRINGTIE_MERGE_SHORT {
    container 'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2'
    publishDir "${params.output_dir}/short_read/04_quantification/stringtie", mode: 'copy'

    input:
    path gtf_files          // All individual sample GTFs collected
    path reference_gtf

    output:
    path "merged.gtf", emit: merged_gtf
    path "merged_summary.txt", emit: summary

    script:
    """
    # Create file listing all GTF files for StringTie merge
    ls *.gtf > gtf_list.txt

    # Merge all sample assemblies into unified reference
    stringtie --merge \\
        -G ${reference_gtf} \\
        -o merged.gtf \\
        -p ${task.cpus} \\
        gtf_list.txt

    # Create summary of merged transcripts
    echo "==================================================================" > merged_summary.txt
    echo "           STRINGTIE MERGE SUMMARY" >> merged_summary.txt
    echo "==================================================================" >> merged_summary.txt
    echo "" >> merged_summary.txt
    echo "Number of input GTF files: \$(cat gtf_list.txt | wc -l)" >> merged_summary.txt
    echo "" >> merged_summary.txt
    echo "Merged transcript statistics:" >> merged_summary.txt
    echo "  Total transcripts: \$(grep -v '^#' merged.gtf | grep -c 'transcript')" >> merged_summary.txt
    echo "  Total exons: \$(grep -v '^#' merged.gtf | grep -c 'exon')" >> merged_summary.txt
    echo "" >> merged_summary.txt
    echo "Generated: \$(date)" >> merged_summary.txt
    echo "==================================================================" >> merged_summary.txt
    """

    // --merge: merge mode - combines multiple GTF assemblies
    // -G: reference annotation to guide merging
    // -o: output merged GTF file
    // -p: number of threads
    // gtf_list.txt: file containing paths to all sample GTFs
}
