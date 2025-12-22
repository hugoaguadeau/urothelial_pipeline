// Process: Minimap2 Alignment Summary - Extract mapping statistics
process MINIMAP2_SUMMARY_LONG {
    container 'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0'
    publishDir "${params.output_dir}/long_read/03_alignment/summaries", mode: 'copy'

    input:
    path bam_files  // All the BAM files from minimap2

    output:
    path "minimap2_summary.csv", emit: summary
    path "minimap2_summary.txt", emit: summary_txt

    script:
    """
    #!/bin/bash
    
    # Create CSV header
    echo "sample_id,total_reads,mapped_reads,mapped_percent,unmapped_reads,unmapped_percent,primary_mapped,secondary_alignments,supplementary_alignments" > minimap2_summary.csv

    # Process each BAM file
    for bam in *.sorted.bam; do
        # Extract sample ID from filename
        sample_id=\${bam%.sorted.bam}
        
        # Run samtools flagstat
        flagstat_output=\$(samtools flagstat \$bam)
        
        # Extract statistics
        total_reads=\$(echo "\$flagstat_output" | grep "in total" | awk '{print \$1}')
        mapped_reads=\$(echo "\$flagstat_output" | grep "mapped (" | head -1 | awk '{print \$1}')
        primary_mapped=\$(echo "\$flagstat_output" | grep "primary mapped" | awk '{print \$1}')
        secondary=\$(echo "\$flagstat_output" | grep "secondary" | awk '{print \$1}')
        supplementary=\$(echo "\$flagstat_output" | grep "supplementary" | awk '{print \$1}')
        
        # Calculate values using awk
        result=\$(awk -v total="\$total_reads" -v mapped="\$mapped_reads" 'BEGIN {
            unmapped = total - mapped
            if (total > 0) {
                mapped_pct = (mapped / total) * 100
                unmapped_pct = (unmapped / total) * 100
            } else {
                mapped_pct = 0
                unmapped_pct = 0
            }
            printf "%d,%.2f,%d,%.2f", unmapped, mapped_pct, unmapped, unmapped_pct
        }')
        
        # Split result
        unmapped_reads=\$(echo "\$result" | cut -d',' -f1)
        mapped_percent=\$(echo "\$result" | cut -d',' -f2)
        unmapped_percent=\$(echo "\$result" | cut -d',' -f4)
        
        # Write to CSV
        echo "\$sample_id,\$total_reads,\$mapped_reads,\$mapped_percent,\$unmapped_reads,\$unmapped_percent,\$primary_mapped,\$secondary,\$supplementary" >> minimap2_summary.csv
    done

    # Create text summary
    current_date=\$(date '+%Y-%m-%d %H:%M:%S')
    cat > minimap2_summary.txt << EOF
================================================================================
                     MINIMAP2 ALIGNMENT SUMMARY
                            (LONG READS)
================================================================================

Mapping statistics from Minimap2 aligner
Generated: \$current_date

EOF

    # Add per-sample details
    tail -n +2 minimap2_summary.csv | while IFS=, read -r sample total mapped map_pct unmapped unmap_pct primary secondary supp; do
        cat >> minimap2_summary.txt << SAMPLE
--------------------------------------------------------------------------------
SAMPLE: \$sample
--------------------------------------------------------------------------------
Total reads:              \$(printf "%'d" \$total 2>/dev/null || echo \$total)
Mapped reads:             \$(printf "%'d" \$mapped 2>/dev/null || echo \$mapped) (\${map_pct}%)
Unmapped reads:           \$(printf "%'d" \$unmapped 2>/dev/null || echo \$unmapped) (\${unmap_pct}%)
Primary mapped:           \$(printf "%'d" \$primary 2>/dev/null || echo \$primary)
Secondary alignments:     \$(printf "%'d" \$secondary 2>/dev/null || echo \$secondary)
Supplementary alignments: \$(printf "%'d" \$supp 2>/dev/null || echo \$supp)

SAMPLE
    done

    # Add overall summary
    cat >> minimap2_summary.txt << 'FOOTER'
================================================================================
OVERALL SUMMARY:
--------------------------------------------------------------------------------
FOOTER

    # Calculate overall statistics
    awk -F',' 'NR>1 {
        samples++
        sum_total += \$2
        sum_mapped += \$3
        sum_map_pct += \$4
    }
    END {
        if (samples > 0) {
            printf "Total samples:              %d\\n", samples
            printf "Total reads:                %'"'"'d\\n", sum_total
            printf "Total mapped:               %'"'"'d\\n", sum_mapped
            printf "Mean mapping rate:          %.2f%%\\n", sum_map_pct/samples
        }
    }' minimap2_summary.csv >> minimap2_summary.txt

    echo "================================================================================" >> minimap2_summary.txt
    """
}
