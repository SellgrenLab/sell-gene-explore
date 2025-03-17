#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://cloud.r-project.org"))


# Load necessary libraries
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
})
# Define command-line options
option_list <- list(
    make_option(c("-g", "--gene"), type = "character", default = NULL, help = "Input gene name", metavar = "CHAR"),
    make_option(c("-a", "--output1"), type = "character", default = "Scatter_bulk_brain_dev.pdf", help = "Scatter plot - Bulk RNAseq", metavar = "FILE"),
    make_option(c("-b", "--output2"), type = "character", default = "Exp_sc_early_dev.pdf", help = "Dotplot & Heatmap - scRNAseq- 1st trimester", metavar = "FILE"),
    make_option(c("-c", "--output3"), type = "character", default = "Exp_sc_all_brain_dev.pdf", help = "Dotplot & Heatmap - scRNAseq- 2nd trimester-Adult", metavar = "FILE"),
    make_option(c("-p", "--palette"), type = "character", default = "YlGnBu", help = "Colour Palette", metavar = "CHAR"),
    make_option(c("-l", "--log"), type = "character", default = "output.log", help = "Log file to save output", metavar = "FILE")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))


# Validate input
if (is.null(opt$gene)) {
  stop("Error: Input file is required. Use -g or --gene to specify the Gene name.", call. = FALSE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_packages <- c("Seurat", "ggplot2", "ggpubr", "dplyr", "ggrepel", 
                                     "cowplot", "RColorBrewer", "grid", "openxlsx", "circlize")

for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

bioc_packages <- c("rhdf5", "ComplexHeatmap")

for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("SeuratExtend", quietly = TRUE)) remotes::install_github("huayc09/SeuratExtend")

suppressPackageStartupMessages({
    require(Seurat)
    require(SeuratExtend)
    require(ggplot2)
    require(ggpubr)
    require(dplyr)
    require(ggrepel)
    require(cowplot)
    require(RColorBrewer)
    require(grid)
    require(rhdf5)
    require(openxlsx)
    require(ComplexHeatmap)
    require(circlize)
})


gene <- opt$gene
pal <- opt$palette

main_plots <- list()

# Create results directory if it doesn't exist
results_dir <- "results"
if (!dir.exists(results_dir)) {
    dir.create(results_dir)
}

# Update output file paths to include the results directory
opt$output1 <- file.path(results_dir, opt$output1)
opt$output2 <- file.path(results_dir, opt$output2)
opt$output3 <- file.path(results_dir, opt$output3)
opt$log <- file.path(results_dir, opt$log)


########################## BRAINSPAN DATA ##########################
expDf <- read.csv("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/External datasets/Bulk Studies/BrainSpan/SCZ_PPI/data/BrainSpan/array_expression_matrix.processed.csv", header=T, stringsAsFactor=F)
geneDf <- read.csv("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/External datasets/Bulk Studies/BrainSpan/SCZ_PPI/data/BrainSpan/rows_metadata.filtered.csv", stringsAsFactor=F)
sampDf <- read.csv("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/External datasets/Bulk Studies/BrainSpan/SCZ_PPI/data/BrainSpan/columns_metadata.csv", stringsAsFactors=F)
ageDf <- read.csv("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/External datasets/Bulk Studies/BrainSpan/SCZ_PPI/data/BrainSpan/age.mapping.csv", stringsAsFactors=F)
sampDf$TimePoint <- plyr::mapvalues(sampDf$age, ageDf$Age, ageDf$Developmental.Period)
sampCounts <- as.data.frame(table(sampDf$TimePoint))
colnames(sampCounts) <- c('TimePoint','NumSamples')

print("Loaded BrainSpan data . . . ")
if (!(gene %in% geneDf$gene_symbol)) {
    stop(paste("Error: Gene", gene, "is not present in the BrainSpan dataset. Check other gene symbol aliases before giving up."), call. = FALSE)
}
tempDf <- expDf[geneDf$row_num[geneDf$gene_symbol == gene],]
missing_samples <- setdiff(sampDf$Identifier, colnames(tempDf))
sampDf <- sampDf[!sampDf$Identifier %in% missing_samples, ]
tempDf <- tempDf[, colnames(tempDf) %in% sampDf$Identifier]
plotDf <- cbind(sampDf, t(tempDf))
colnames(plotDf)[11] <- gene

plotDf <- plotDf %>% mutate(Region = ifelse(structure_acronym %in% c("DFC", "MFC", "PC", "VFC", "M1C", "OFC", "FC"), "NeoCortex",
                                                 ifelse(structure_acronym %in% c("A1C", "ITC", "STC", "TC","TCx"), "NeoCortex",
                                                                ifelse(structure_acronym %in% c("IPC", "S1C", "PC","PCx"), "NeoCortex",
                                                                             ifelse(structure_acronym %in% c("OC", "V1C", "Ocx"), "NeoCortex",
                                                                                            ifelse(structure_acronym %in% c("CB", "CBC", "URL"), "Cerebellum",
                                                                                                         ifelse(structure_acronym %in% c("DTH", "MD", "DIE"), "Thalamus",
                                                                                                                        ifelse(structure_acronym %in% c("CGE", "LGE", "MGE", "STR"), "Striatum",
                                                                                                                                     ifelse(structure_acronym == "HIP", "Hippocampus", 
                                                                                                                                                    ifelse(structure_acronym == "AMY","Amygdala", NA))))))))))

plotDf$age1 <- factor(plotDf$age, levels = unique(plotDf$age[order(gsub("[^a-zA-Z]", "", plotDf$age), as.numeric(gsub("[^0-9]", "", plotDf$age)))]))
plotDf$age1 <- factor(plotDf$age1, levels = c(levels(plotDf$age1)[grep("pcw", levels(plotDf$age1))], levels(plotDf$age1)[grep("mos", levels(plotDf$age1))], levels(plotDf$age1)[grep("yrs", levels(plotDf$age1))]))
plotDf$plot <- as.numeric(plotDf$age1)
plotDf$Region <- factor(plotDf$Region, levels = c("NeoCortex", "Hippocampus", "Striatum", "Thalamus", "Amygdala", "Cerebellum"))

print("BrainSpan data is ready to plot . . . ")

p1 <- ggscatter(plotDf, x = "plot", y = gene, 
            add = "loess",
            color = "Region", 
            palette = "npg",
            xlab = "Age (8 pcw - 40 yrs)", 
            ylab = "Expression Level",
            title = paste("Expression of", gene, "across brain development")) +
                geom_vline(xintercept = which(levels(plotDf$age1) == "26 pcw"), linetype = "dashed", color = "black") +
                theme_minimal() +
                theme(axis.text.x = element_blank())

main_plots[[1]] <- p1

# Save the plot with the gene name
output_file_gene <- file.path(results_dir, paste0(gene, "_Scatter_bulk_brain_dev.pdf"))
ggsave(output_file_gene, plot = main_plots[[1]], width = 8, height = 6)
#ggsave(opt$output1, plot = main_plots[[1]], width = 8, height = 6)


rm(list = c("expDf", "geneDf", "sampDf", "ageDf", "tempDf", "missing_samples", "plotDf"))


######################### First Trimester Single cell Data #########################
print(paste("Loading Fetal First trimester data . . . "))
meta <- read.xlsx("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/External datasets/Single Cell Studies/Linnarsson fetal human data/humanlin/table_S2.xlsx")
meta <- meta[-which(meta$PoolCleanOrder=="--"),]
meta$PoolCleanOrder <- as.numeric(meta$PoolCleanOrder)
meta1 <- dplyr::arrange(meta, PoolCleanOrder)
meta <- meta1
rm(meta1)

#h5file <- H5Fopen("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/External datasets/Single Cell Studies/Linnarsson fetal human data/humanlin/HumanFetalBrainPool.h5")
#lin <- h5read(h5file, name = "shoji/MeanExpression")
#genes_lin <- h5read(h5file, name = "shoji/Gene")
#H5Fclose(h5file) #throws an error, possibly segmentation fault

lin<- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/External datasets/Single Cell Studies/Linnarsson fetal human data/humanlin/Human_Fetal_Brain_Linnarsson_MeanExpression.rds")
genes_lin <- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/External datasets/Single Cell Studies/Linnarsson fetal human data/humanlin/Human_Fetal_Brain_Linnarsson_Gene.rds")
lin <- as.data.frame(lin)
rownames(lin) <- make.unique(genes_lin)
colnames(lin) <- meta$CommonestRegion
if (!(gene %in% rownames(lin))) {
    stop(paste("Error: Gene", gene, "is not present in the Linnarsson dataset. Check other gene symbol aliases before giving up."), call. = FALSE)
} 
lin_sub <- lin[gene,]


print("Preparing fetal first trimester data for plotting . . . ")
unique_regions <- unique(colnames(lin_sub))
lin_grouped <- sapply(unique_regions, function(region) rowMeans(lin_sub[, colnames(lin_sub) == region, drop = FALSE], na.rm = TRUE))
lin_grouped <- as.data.frame(t(lin_grouped))
rownames(lin_grouped) <- rownames(lin_sub)
colnames(lin_grouped) <- gsub("\\..*", "", colnames(lin_grouped))

p2 <- Heatmap(as.matrix(lin_grouped), 
    name = "Expression", 
    col = colorRampPalette(brewer.pal(9, pal))(100),
    cluster_rows = TRUE, 
    cluster_columns = TRUE,
    show_row_names = TRUE, 
    show_column_names = TRUE,
    column_title = "Regions",
    row_title = "Genes"
)

p2_final <- cowplot::plot_grid(grid.grabExpr(draw(p2)))

colnames(lin_sub) <- meta$CommonestRegion
region_indices <- lapply(unique(meta$CommonestRegion), function(region) which(colnames(lin_sub) == region))

max_val <- c()
for(i in seq_along(region_indices)){
 region <- unique(meta$CommonestRegion)[i]
    indices <- region_indices[[i]]
    lin_region <- lin_sub[, indices]
    colnames(lin_region) <- meta$Subclass[indices]
    unique_celltype <- unique(colnames(lin_region))
    lin_grouped <- sapply(unique_celltype, function(celltype) rowMeans(lin_region[, colnames(lin_region) == celltype, drop = FALSE], na.rm = TRUE))
    lin_grouped <- as.data.frame(t(lin_grouped))
    rownames(lin_grouped) <- rownames(lin_region)
    max_val[i] <- max(lin_grouped, na.rm = TRUE)
}

plots_lin <- list()
for (i in seq_along(region_indices)) {
    region <- unique(meta$CommonestRegion)[i]
    indices <- region_indices[[i]]
    lin_region <- lin_sub[, indices]
    colnames(lin_region) <- meta$Subclass[indices]
    unique_celltype <- unique(colnames(lin_region))
    lin_grouped <- sapply(unique_celltype, function(celltype) rowMeans(lin_region[, colnames(lin_region) == celltype, drop = FALSE], na.rm = TRUE))
    lin_grouped <- as.data.frame(t(lin_grouped))
    rownames(lin_grouped) <- rownames(lin_region)
    colnames(lin_grouped) <- gsub("\\..*", "", colnames(lin_grouped))
    
    if (length(unique(as.matrix(lin_grouped))) > 1) {
        if (all(lin_grouped == 0)) {
            message(paste("All values are zero for region:", region))
            next
        }
        plots_lin[[i]] <- Heatmap(as.matrix(lin_grouped),  
            name = "Mean",
            col = colorRamp2(c(0, max(max_val)), colorRampPalette(brewer.pal(9, pal))(2)),
            cluster_rows = TRUE, 
            cluster_columns = TRUE,
            show_row_names = TRUE, 
            show_column_names = TRUE,
            column_title = region,
            row_title = "Expression"
        )
    } else {
        message(paste("Not enough distinct values to plot heatmap for region:", region))
    }
}
print("Fetal first trimester data is ready to plot . . . ")

grobs <- lapply(plots_lin, function(p) grid.grabExpr(draw(p)))
p3 <- cowplot::plot_grid(plotlist = grobs, ncol = 2)

title1 <- ggdraw() + 
  draw_label(
    paste("Expression of", gene, "across brain regions in fetal first trimester"),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
main_plots[[2]] <- cowplot::plot_grid(title1, p2_final, p3, nrow = 3, rel_heights = c(0.1,2, 8))

# Save the plot with the gene name
output_file_gene <- file.path(results_dir, paste0(gene, "_Exp_sc_early_dev.pdf"))
ggsave(output_file_gene, main_plots[[2]], width = 9, height=15)
#ggsave(opt$output2, main_plots[[2]], width = 9, height=15)

rm(list = c("lin","genes_lin", "lin_grouped", "lin_sub", "plots_lin", "grobs"))

######################### Pan Development Single cell Data #########################
print("Loading pan development data . . . ")

velm <- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/External datasets/Single Cell Studies/kriegstein human adolescence/VelmeshevLab_human_brain23.rds")
velm$age <- factor(velm$age, levels = c("2nd trimester", "3rd trimester", "0-1 years", "1-2 years", "2-4 years", "4-10 years", "10-20 years", "Adult"))
velm$plot <- paste0(velm$age, "_", velm$celltype)
Idents(velm) <- "plot"
if (!(gene %in% rownames(velm))) {
    stop(paste("Error: Gene", gene, "is not present in the Velmeshev dataset. Check other gene symbol aliases before giving up."), call. = FALSE)
} 
p5 <- DotPlot2(velm, features = gene, group.by = "age", color_scheme = "YlGnBu") 
p6 <- DotPlot2(velm, features = gene, group.by = "celltype", color_scheme = "YlGnBu") 

toplot2 <- AverageExpression(velm, features = gene, group.by = c("age", "celltype") )$RNA
colnames(toplot2) <- sub("^g", "", colnames(toplot2))
names_of_df <- colnames(toplot2)
cell_groups <- sapply(strsplit(colnames(toplot2), "_"), `[`, 1)
cell_groups <- factor(cell_groups, levels = c("2nd trimester", "3rd trimester", "0-1 years", "1-2 years", "2-4 years", "4-10 years", "10-20 years", "Adult"))
toplot2 <- toplot2[order(cell_groups)]
toplot2 <- as.data.frame(toplot2)
colnames(toplot2) <- gene
rownames(toplot2) <- names_of_df #sapply(strsplit(names_of_df, "_"), `[`, 2)

print("Pan development data is ready to plot . . . ")
p7 <- SeuratExtend::Heatmap(toplot2, features=gene, facet_row = cell_groups, color_scheme="YlGnBu", lab_fill="Expression")

#main_plots[[3]] <- cowplot::plot_grid(p5, p6, p7, nrow = 3, rel_heights = c(1, 1,7))
blank_plot <- ggplot() + theme_void()

title2 <- ggdraw() + 
  draw_label(
    paste("Expression of", gene, "across brain regions throughout brain development"),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
left_col <- cowplot::plot_grid(title2, p5, p6, blank_plot, labels = c('', 'A', 'B', ''), label_size = 10, ncol=1, rel_heights = c(1.5, 1, 1, 1.5))
main_plots[[3]] <- cowplot::plot_grid(left_col,p7, labels = c('', 'C'), label_size = 10, ncol = 2, rel_widths = c(2,1.5),rel_heights=c(1,10))

# Save the plot with the gene name
output_file_gene <- file.path(results_dir, paste0(gene, "_Exp_sc_all_brain_dev.pdf"))
ggsave(output_file_gene, main_plots[[3]], height=12, width=10)
#ggsave(opt$output3, main_plots[[3]], height=12, width=10)

cat("Plots saved:\n", opt$output1, "\n", opt$output2, "\n", opt$output3, "\n")

#devtools::session_info()

