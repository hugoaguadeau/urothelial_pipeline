// Process: Enhanced Volcano Plot - Generate differential expression visualisation
process ENHANCED_VOLCANO_LONG {
    container 'https://depot.galaxyproject.org/singularity/bioconductor-enhancedvolcano:1.24.0--r44hdfd78af_0'
    publishDir "${params.output_dir}/long_read/05_differential_expression/plots/volcano", mode: 'copy'

    input:
    path deseq2_results

    output:
    path "volcano_plot.pdf", emit: pdf
    path "volcano_plot.png", emit: png

    script:
    """
    #!/usr/bin/env Rscript
    library(EnhancedVolcano)
    library(ggplot2)
    
    # Load DESeq2 results table
    res_data <- read.csv("${deseq2_results}", row.names = 1)
    
    # Filter out rows containing NA values
    res_data <- res_data[!is.na(res_data\$padj) & !is.na(res_data\$log2FoldChange), ]
    
    # Calculate axis limits with additional headroom for labels
    max_y <- max(-log10(res_data\$padj[res_data\$padj > 0]), na.rm=TRUE)
    my_ylim <- c(0, max_y * 1.2) 
    
    max_x <- max(abs(res_data\$log2FoldChange), na.rm=TRUE)
    my_xlim <- c(-max_x * 1.1, max_x * 1.1) 
    
    # Define colour scheme for differentially expressed genes
    # Red represents elevated expression in Bladder
    # Blue represents elevated expression in Upper tract (reference condition)
    keyvals <- ifelse(
        res_data\$log2FoldChange > 1.0 & res_data\$padj < 0.05, 'firebrick3',
        ifelse(res_data\$log2FoldChange < -1.0 & res_data\$padj < 0.05, 'royalblue3',
        'grey75'))
    keyvals[is.na(keyvals)] <- 'grey75'
    names(keyvals)[keyvals == 'firebrick3'] <- 'High in Bladder'
    names(keyvals)[keyvals == 'royalblue3'] <- 'High in Uppertract'
    names(keyvals)[keyvals == 'grey75']     <- 'NS'
    
    # Format gene names in italics for display
    lab_italics <- paste0("italic('", rownames(res_data), "')")
    
    # Select top 25 genes by adjusted p-value and absolute fold change
    res_sorted <- res_data[order(res_data\$padj, -abs(res_data\$log2FoldChange)), ]
    top_genes <- head(rownames(res_sorted), 25)
    selectLab_italics <- paste0("italic('", top_genes, "')")
    
    # Generate volcano plot with custom styling
    p <- EnhancedVolcano(res_data,
        lab = lab_italics,
        x = 'log2FoldChange',
        y = 'padj',
        selectLab = selectLab_italics,
        xlab = bquote(~Log[2]~ 'Fold Change'),
        ylab = bquote(~-Log[10]~ 'Adjusted' ~italic(P)),
        title = 'Differential Gene Expression',
        subtitle = 'Urothelial Carcinoma: Bladder vs Upper tract',
        
        # Apply calculated axis limits
        ylim = my_ylim,
        xlim = my_xlim,
        
        # Statistical thresholds
        pCutoff = 0.05,
        FCcutoff = 1.0,
        
        # Visual styling parameters
        pointSize = 3.5,
        labSize = 8.0,
        labCol = 'black',
        labFace = 'bold',
        boxedLabels = TRUE,
        parseLabels = TRUE,
        colAlpha = 0.6,
        
        # Label connector configuration to prevent overlap
        drawConnectors = TRUE,
        widthConnectors = 0.6,
        colConnectors = 'grey30',
        arrowheads = FALSE,
        max.overlaps = Inf,
        
        # Apply custom colour scheme
        colCustom = keyvals,
        
        # Text size parameters
        axisLabSize = 30,
        titleLabSize = 34,
        subtitleLabSize = 28,
        captionLabSize = 24,
        legendLabSize = 24,
        legendIconSize = 8.0,
        
        # Legend positioning
        legendPosition = 'top'
    )
    
    # Apply horizontal legend orientation
    p <- p + theme(legend.direction = "horizontal")
    
    # Export plot in PDF and PNG formats
    pdf("volcano_plot.pdf", width = 14, height = 14)
    print(p)
    dev.off()
    
    png("volcano_plot.png", width = 1400, height = 1400, res = 120)
    print(p)
    dev.off()
    """
}