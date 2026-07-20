import os
import sys
import pandas as pd

# This script merges two DEG files, applies significance thresholds, and
# generates a new file with color-coded genes, as well as an R script
# to create UpSetR and scatter plots with dynamic naming.

# Command-line arguments:
# sys.argv[1]: path to first DEG file (e.g., factor1.deg.txt)
# sys.argv[2]: name for factor 1 (e.g., "Week6")
# sys.argv[3]: path to second DEG file (e.g., factor2.deg.txt)
# sys.argv[4]: name for factor 2 (e.g., "Week30_33")
# sys.argv[5]: output prefix (e.g., "my_analysis")
# sys.argv[6]: padj threshold (e.g., 0.05)
# sys.argv[7]: log2FoldChange threshold (e.g., 1.0)
# sys.argv[8]: plot title (e.g., "My Plot Title")

# --- Python Logic ---

# Load first DEG file
factor1DEG = pd.read_csv(sys.argv[1], sep="\t", header=0, index_col=None)
# Filter for Gene, GeneSymbol, and relevant values
factor1DEG = factor1DEG[['Gene', 'GeneSymbol', 'log2FoldChange', 'padj']]
# Use dynamic factor name for column headers
factor1DEG.columns = ['Gene', 'GeneSymbol', sys.argv[2] + "_log2FC", sys.argv[2] + "_Padj"]
# Filter out non-informative gene symbols
factor1DEG = factor1DEG[~factor1DEG['GeneSymbol'].str.contains('^LOC|^Gm[0-9]+|Rik$', na=False)]
factor1DEG = factor1DEG.dropna(subset=['GeneSymbol'])
factor1DEG = factor1DEG.drop_duplicates(subset=['GeneSymbol'], keep='first')

# Load second DEG file
factor2DEG = pd.read_csv(sys.argv[3], sep="\t", header=0, index_col=None)
# Filter for Gene, GeneSymbol, and relevant values
factor2DEG = factor2DEG[['Gene', 'GeneSymbol', 'log2FoldChange', 'padj']]
# Use dynamic factor name for column headers
factor2DEG.columns = ['Gene', 'GeneSymbol', sys.argv[4] + "_log2FC", sys.argv[4] + "_Padj"]
# Filter out non-informative gene symbols
factor2DEG = factor2DEG[~factor2DEG['GeneSymbol'].str.contains('^LOC|^Gm[0-9]+|Rik$', na=False)]
factor2DEG = factor2DEG.dropna(subset=['GeneSymbol'])
factor2DEG = factor2DEG.drop_duplicates(subset=['GeneSymbol'], keep='first')

# Merge the data using 'GeneSymbol' and 'Gene' as the common identifiers
big = pd.merge(factor1DEG, factor2DEG, on=['Gene', 'GeneSymbol'], how="outer")

# Set thresholds from command-line arguments
padj_threshold = float(sys.argv[6])
log2FC_threshold = float(sys.argv[7])

# Define significance for each factor based on BOTH padj and log2FoldChange thresholds
is_significant_factor1 = (big[sys.argv[2] + "_Padj"] <= padj_threshold) & \
                         (big[sys.argv[2] + "_log2FC"].abs() >= log2FC_threshold)

is_significant_factor2 = (big[sys.argv[4] + "_Padj"] <= padj_threshold) & \
                         (big[sys.argv[4] + "_log2FC"].abs() >= log2FC_threshold)

# Assign colors based on these robust significance flags
big['Color'] = 'Grey'  # Default color for non-significant genes
big.loc[is_significant_factor1 & is_significant_factor2, 'Color'] = 'Red'
big.loc[is_significant_factor1 & ~is_significant_factor2, 'Color'] = 'Orange'
big.loc[~is_significant_factor1 & is_significant_factor2, 'Color'] = 'Blue'

# Save the output file with the dynamic prefix
output_file = sys.argv[5]
big.to_csv(output_file + "_DEGs.txt", sep="\t", header=True, index=False)

# --- R Script Generation ---

# Generate the R script with placeholders that will be replaced dynamically.
fw = open("PlotCode33." + sys.argv[5] + ".r", "w")
code = '''
############### UpSetR plot ##################
library("UpSetR")
library(data.table)

args = commandArgs(trailingOnly=TRUE)
# R script arguments:
# args[1]: merged DEG file name (e.g., "my_analysis_DEGs.txt")
# args[2]: output file prefix (e.g., "my_analysis")
# args[3]: plot title (e.g., "My Plot Title")
# args[4]: padj_threshold from Python
# args[5]: log2FC_threshold from Python
# args[6]: factor1_name from Python (e.g., "Week6")
# args[7]: factor2_name from Python (e.g., "Week30_33")

frame <- fread(args[1])
plot_title <- args[3]
padj_cutoff_R = as.numeric(args[4])
log2fc_cutoff_R = as.numeric(args[5])
factor1_name_R = args[6]
factor2_name_R = args[7]

# The grouping is now done by 'GeneSymbol'
frame[,nrow:=1:.N,by=list(GeneSymbol)]

# Create dynamic column names for the UpSetR data frame using factor names
frame[, List1_Up := ifelse(is.na(List1_Padj) | is.na(List1_log2FC), 0,
                           ifelse(List1_Padj <= padj_cutoff_R & List1_log2FC >= log2fc_cutoff_R, 1, 0))]
frame[, List1_Dn := ifelse(is.na(List1_Padj) | is.na(List1_log2FC), 0,
                           ifelse(List1_Padj <= padj_cutoff_R & List1_log2FC <= -log2fc_cutoff_R, 1, 0))]
frame[, List2_Up := ifelse(is.na(List2_Padj) | is.na(List2_log2FC), 0,
                           ifelse(List2_Padj <= padj_cutoff_R & List2_log2FC >= log2fc_cutoff_R, 1, 0))]
frame[, List2_Dn := ifelse(is.na(List2_Padj) | is.na(List2_log2FC), 0,
                           ifelse(List2_Padj <= padj_cutoff_R & List2_log2FC <= -log2fc_cutoff_R, 1, 0))]

############### Convert to a data frame for subsetting ###############
# The R script now correctly uses 'GeneSymbol' instead of 'Gene'
# The dframe is specifically for the UpSetR plot, containing 1s and 0s
dframe <- data.frame(frame[, c('GeneSymbol', 'List1_Up', 'List1_Dn', 'List2_Up', 'List2_Dn')])
dframe[is.na(dframe)] <- 0

############### Save each bar (set/intersection) to separate text files ###############
# Define column names to select for each output file, including 'Gene'
cols_factor1 <- c('Gene', 'GeneSymbol', paste0(factor1_name_R, '_log2FC'))
cols_factor2 <- c('Gene', 'GeneSymbol', paste0(factor2_name_R, '_log2FC'))
cols_both <- c('Gene', 'GeneSymbol', paste0(factor1_name_R, '_log2FC'), paste0(factor2_name_R, '_log2FC'))

# The file names are now dynamically generated using the factor names (args[6] and args[7])
# The subsetting now includes the relevant log2FoldChange columns, as requested.
List1_Up_only <- subset(frame, List1_Up == 1 & List1_Dn == 0 & List2_Up == 0 & List2_Dn == 0)
if (nrow(List1_Up_only) > 0) {
    fwrite(List1_Up_only[, ..cols_factor1], paste0(factor1_name_R, "_Up_only.txt"), sep="\t")
} else {
    fwrite(data.table(matrix(nrow=0, ncol=length(cols_factor1), dimnames=list(NULL, cols_factor1))),
           paste0(factor1_name_R, "_Up_only.txt"), sep="\t")
}

List1_Dn_only <- subset(frame, List1_Dn == 1 & List1_Up == 0 & List2_Up == 0 & List2_Dn == 0)
if (nrow(List1_Dn_only) > 0) {
    fwrite(List1_Dn_only[, ..cols_factor1], paste0(factor1_name_R, "_Dn_only.txt"), sep="\t")
} else {
    fwrite(data.table(matrix(nrow=0, ncol=length(cols_factor1), dimnames=list(NULL, cols_factor1))),
           paste0(factor1_name_R, "_Dn_only.txt"), sep="\t")
}

List2_Up_only <- subset(frame, List2_Up == 1 & List2_Dn == 0 & List1_Up == 0 & List1_Dn == 0)
if (nrow(List2_Up_only) > 0) {
    fwrite(List2_Up_only[, ..cols_factor2], paste0(factor2_name_R, "_Up_only.txt"), sep="\t")
} else {
    fwrite(data.table(matrix(nrow=0, ncol=length(cols_factor2), dimnames=list(NULL, cols_factor2))),
           paste0(factor2_name_R, "_Up_only.txt"), sep="\t")
}

List2_Dn_only <- subset(frame, List2_Dn == 1 & List2_Up == 0 & List1_Up == 0 & List1_Dn == 0)
if (nrow(List2_Dn_only) > 0) {
    fwrite(List2_Dn_only[, ..cols_factor2], paste0(factor2_name_R, "_Dn_only.txt"), sep="\t")
} else {
    fwrite(data.table(matrix(nrow=0, ncol=length(cols_factor2), dimnames=list(NULL, cols_factor2))),
           paste0(factor2_name_R, "_Dn_only.txt"), sep="\t")
}

# Add the new files for all Up/Dn genes for each factor
List1_Up <- subset(frame, List1_Up == 1)
if (nrow(List1_Up) > 0) {
    fwrite(List1_Up[, ..cols_factor1], paste0(factor1_name_R, "_Up.txt"), sep="\t")
} else {
    fwrite(data.table(matrix(nrow=0, ncol=length(cols_factor1), dimnames=list(NULL, cols_factor1))),
           paste0(factor1_name_R, "_Up.txt"), sep="\t")
}

List1_Dn <- subset(frame, List1_Dn == 1)
if (nrow(List1_Dn) > 0) {
    fwrite(List1_Dn[, ..cols_factor1], paste0(factor1_name_R, "_Dn.txt"), sep="\t")
} else {
    fwrite(data.table(matrix(nrow=0, ncol=length(cols_factor1), dimnames=list(NULL, cols_factor1))),
           paste0(factor1_name_R, "_Dn.txt"), sep="\t")
}

List2_Up <- subset(frame, List2_Up == 1)
if (nrow(List2_Up) > 0) {
    fwrite(List2_Up[, ..cols_factor2], paste0(factor2_name_R, "_Up.txt"), sep="\t")
} else {
    fwrite(data.table(matrix(nrow=0, ncol=length(cols_factor2), dimnames=list(NULL, cols_factor2))),
           paste0(factor2_name_R, "_Up.txt"), sep="\t")
}

List2_Dn <- subset(frame, List2_Dn == 1)
if (nrow(List2_Dn) > 0) {
    fwrite(List2_Dn[, ..cols_factor2], paste0(factor2_name_R, "_Dn.txt"), sep="\t")
} else {
    fwrite(data.table(matrix(nrow=0, ncol=length(cols_factor2), dimnames=list(NULL, cols_factor2))),
           paste0(factor2_name_R, "_Dn.txt"), sep="\t")
}

# Intersections (i.e., in multiple sets simultaneously)
List1_Dn_and_List2_Dn <- subset(frame, List1_Dn == 1 & List2_Dn == 1)
if (nrow(List1_Dn_and_List2_Dn) > 0) {
    fwrite(List1_Dn_and_List2_Dn[, ..cols_both], paste0(factor1_name_R, "_Dn_and_", factor2_name_R, "_Dn.txt"), sep="\t")
} else {
    fwrite(data.table(matrix(nrow=0, ncol=length(cols_both), dimnames=list(NULL, cols_both))),
           paste0(factor1_name_R, "_Dn_and_", factor2_name_R, "_Dn.txt"), sep="\t")
}

List1_Dn_and_List2_Up <- subset(frame, List1_Dn == 1 & List2_Up == 1)
if (nrow(List1_Dn_and_List2_Up) > 0) {
    fwrite(List1_Dn_and_List2_Up[, ..cols_both], paste0(factor1_name_R, "_Dn_and_", factor2_name_R, "_Up.txt"), sep="\t")
} else {
    fwrite(data.table(matrix(nrow=0, ncol=length(cols_both), dimnames=list(NULL, cols_both))),
           paste0(factor1_name_R, "_Dn_and_", factor2_name_R, "_Up.txt"), sep="\t")
}

List1_Up_and_List2_Dn <- subset(frame, List1_Up == 1 & List2_Dn == 1)
if (nrow(List1_Up_and_List2_Dn) > 0) {
    fwrite(List1_Up_and_List2_Dn[, ..cols_both], paste0(factor1_name_R, "_Up_and_", factor2_name_R, "_Dn.txt"), sep="\t")
} else {
    fwrite(data.table(matrix(nrow=0, ncol=length(cols_both), dimnames=list(NULL, cols_both))),
           paste0(factor1_name_R, "_Up_and_", factor2_name_R, "_Dn.txt"), sep="\t")
}

List1_Up_and_List2_Up <- subset(frame, List1_Up == 1 & List2_Up == 1)
if (nrow(List1_Up_and_List2_Up) > 0) {
    fwrite(List1_Up_and_List2_Up[, ..cols_both], paste0(factor1_name_R, "_Up_and_", factor2_name_R, "_Up.txt"), sep="\t")
} else {
    fwrite(data.table(matrix(nrow=0, ncol=length(cols_both), dimnames=list(NULL, cols_both))),
           paste0(factor1_name_R, "_Up_and_", factor2_name_R, "_Up.txt"), sep="\t")
}

# Convert to a data frame for the upset plot
dframe <- data.frame(frame[, c('GeneSymbol', 'List1_Up', 'List1_Dn', 'List2_Up', 'List2_Dn')])
dframe[is.na(dframe)] <- 0

# Initialize queries dynamically based on intersections
queries <- list()
if (sum(dframe$List1_Up & dframe$List2_Up) > 0) {
    queries <- append(queries, list(list(query = intersects, params = list("List1_Up", "List2_Up"), color = "red", active = TRUE)))
}
if (sum(dframe$List1_Dn & dframe$List2_Dn) > 0) {
    queries <- append(queries, list(list(query = intersects, params = list("List1_Dn", "List2_Dn"), color = "blue", active = TRUE)))
}

# Generate the upset plot dynamically based on queries
if (length(queries) > 0) {
    p1 <- upset(
        dframe,
        sets = c('List1_Up', 'List1_Dn', 'List2_Up', 'List2_Dn'),
        point.size = 4,
        order.by = "freq",
        line.size = 0.0,
        text.scale = c(3, 1.5, 3, 1.3, 1.5, 2),
        mainbar.y.label = "Overlap size",
        sets.x.label = "Set Size",
        main.bar.color = "grey50",
        sets.bar.color = "black",
        queries = queries
    )
} else {
    # Plot without queries if there are no intersections
    p1 <- upset(
        dframe,
        sets = c('List1_Up', 'List1_Dn', 'List2_Up', 'List2_Dn'),
        point.size = 4,
        order.by = "freq",
        line.size = 0.0,
        text.scale = c(3, 2, 3, 1.5, 2, 2),
        mainbar.y.label = "Overlap size",
        sets.x.label = "Set Size",
        main.bar.color = "grey50",
        sets.bar.color = "black"
    )
}

# The UpSetR plot file name is now dynamically generated
tiff(paste0(args[2],"-",args[3],"_UpSetRPlot.tiff"),width=7, height=5, units='in',res=150)
p1
dev.off()

################## Scatter Density Plot ##################
library(ggplot2)
library(data.table)
library(ggpmisc)

# The R script receives the _DEGs.txt file as args[1]
frame <- fread(args[1])
plot_title <- args[3]

# Calculate the 99th percentile absolute value for setting symmetrical axis limits
max_limit <- max(
  quantile(abs(frame$List1_log2FC), 0.99, na.rm = TRUE),
  quantile(abs(frame$List2_log2FC), 0.99, na.rm = TRUE)
)

# Set the colors
frame$Color <- factor(frame$Color, levels=c("Red", "Orange", "Blue", "Grey"))

# Filter out the grey dots
frame <- frame[Color != "Grey", ]

# Count the number of genes for the legend labels
red_count <- sum(frame$Color == "Red")
orange_count <- sum(frame$Color == "Orange")
blue_count <- sum(frame$Color == "Blue")

# Plot the scatter plot
# The plot file name and axis labels are now dynamically generated
tiff(paste0(args[2], "_ScatterPlot.tiff"), width=6, height=6, units='in', res=150)
ggplot(frame, aes(x=List1_log2FC, y=List2_log2FC, color=Color)) +
  geom_point(size=2, alpha=0.6) +
  scale_color_manual(
    values = c("Red" = "red", "Orange" = "orange", "Blue" = "blue"),
    labels = c(
      paste("Significant in both (", red_count, ")", sep = ""),
      paste("Significant in X-axis (", orange_count, ")", sep = ""),
      paste("Significant in Y-axis (", blue_count, ")", sep = "")
    )
  ) +
  geom_vline(xintercept=0, linetype="dotted") +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_abline(intercept=0, slope=1, color="black", linetype="solid") +
  geom_smooth(method='lm', se=T, color="black") +
  scale_x_continuous(limits = c(-max_limit, max_limit), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-max_limit, max_limit), expand = c(0, 0)) +
  stat_quadrant_counts() +
  labs(
    title = plot_title,
    x = paste0("log2FC of ", factor1_name_R),
    y = paste0("log2FC of ", factor2_name_R)
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.spacing.x = unit(0.2, 'cm'),
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.caption = element_text(hjust = 0.5, size = 10, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(color = guide_legend(
    title = NULL,
    override.aes = list(size = 3)
  ))
dev.off()
'''

# The key is to replace "List1" and "List2" in the R template string with the
# dynamic factor names from the command line. This replacement is done here.
code = code.replace("List1", sys.argv[2])
code = code.replace("List2", sys.argv[4])
fw.write(code)
fw.close()

# The R script is executed with all necessary arguments, including the dynamic names.
cmd = f"Rscript PlotCode33.{sys.argv[5]}.r {output_file}_DEGs.txt {sys.argv[5]} \"{sys.argv[8]}\" {sys.argv[6]} {sys.argv[7]} {sys.argv[2]} {sys.argv[4]}"
os.system(cmd)
