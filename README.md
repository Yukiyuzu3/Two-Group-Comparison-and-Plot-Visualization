# Two-Group-Comparison-and-Plot-Visualization
Analyzing and Visualizing Intron Retention Using IRFinder
# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)

# Clear all data from the current workspace
rm(list=ls())

# Set the working directory for analysis
path <- "/Users/yukiyuzu3/lab/IRFinder_1/txt_fle_dir"
setwd(path)

# Load the DESeq2 package
library(DESeq2)

# Load the script containing functions related to IRFinder analysis
source("/Users/yukiyuzu3/lab/IRFinder_1/bin/DESeq2Constructor.R")

# Read the file containing paths to IRFinder result files
results = read.table("filePaths.txt")

# Save file paths as a vector
paths = as.vector(results$V1)

# Read the file containing experimental conditions
experiment = read.table("experiment.txt",header=T)

# Set experimental conditions as factors (with "Control" as the reference level)
experiment$Condition=factor(experiment$Condition,levels=c("Control","FECD"))

# Remove row names from the experiment data
rownames(experiment)=NULL

# Create a DESeq2 dataset
metaList=DESeqDataSetFromIRFinder(filePaths=paths, designMatrix=experiment, designFormula=~1)

# Get the normalized DESeq2 object
dds = metaList$DESeq2Object

# Display condition data for the dataset
colData(dds)

# Set the analysis model
design(dds) = ~Condition + Condition:IRFinder

# Run DESeq2 analysis
dds = DESeq(dds)

# Get the names of the results
resultsNames(dds)

# Calculate the ratio of IR reads vs. splice reads for the Control condition
res.Control = results(dds, name = "ConditionControl.IRFinderIR")
Control.IR_vs_Splice=2^res.Control$log2FoldChange
IRratio.Control = Control.IR_vs_Splice/(1+Control.IR_vs_Splice)

# Calculate the ratio of IR reads vs. splice reads for the FECD condition
res.FECD = results(dds, name = "ConditionFECD.IRFinderIR")
FECD.IR_vs_Splice=2^res.FECD$log2FoldChange
IRratio.FECD = FECD.IR_vs_Splice/(1+FECD.IR_vs_Splice)

# Test for differences in IR ratio between Control and FECD conditions
res.diff = results(dds, contrast=list("ConditionFECD.IRFinderIR","ConditionControl.IRFinderIR"))
diff_matrix <-na.omit(data.frame(res.diff))

# Filter for differentially expressed genes (DEGs) with IR
DEGsIR <- diff_matrix[diff_matrix$padj<0.05,]
DEGsIR <- na.omit(DEGsIR)

# Separate IR-upregulated and IR-downregulated genes
IRup <- DEGsIR[DEGsIR$log2FoldChange > 1,]
IRdown <- DEGsIR[DEGsIR$log2FoldChange < -1,]

# Prepare IR-downregulated data for CSV export
IRdown_csv <- rownames(IRdown)
IRdown_csv <- cbind(IRdown_csv,IRdown)

# Prepare IR-upregulated data for CSV export
IRup_csv <- rownames(IRup)
IRup_csv <- cbind(IRup_csv,IRup)

# Calculate and plot the change in IR ratio
IR.change = IRratio.FECD - IRratio.Control
plot(IR.change,col=ifelse(res.diff$padj < 0.05 & abs(IR.change)>=0.1, "red", "black"))

#############################################################################################################################

# Combine padj and log2FoldChange into a new data frame
data2 <- cbind(padj = diff_matrix$padj, log2FoldChange = diff_matrix$log2FoldChange)
dim(data2)

# Transform p-values to -log10(p-value) and remove NA rows
data3 <- data2
data3[,1] <- -log10(data2[,1])
data3 <- na.omit(data3)

############################################################
#### Apply conditions based on LFC and P-value
nrow(data3)

# Initialize a 'group' column for classification
group <- as.numeric(matrix(1:nrow(data3), nrow = nrow(data3), ncol=1))
data4 <- as.data.frame(cbind(group,data3))

lfc.threshold <- 1.0 # Log2 Fold Change threshold

x <- 1
while (x <= nrow(data3)) {
  if (data4[x, 2] < -log10(0.05)) { # If p-value is greater than 0.05
    data4[x, 1] <- "A"            # Assign "A"
  } else {
    if (data4[x, 3] > lfc.threshold) {    # If Log2 Fold Change is greater than 1.0
      data4[x, 1] <- "B"          # Assign "B"
    } else {
      if (data4[x, 3] < -lfc.threshold) { # If Log2 Fold Change is less than -1.0
        data4[x, 1] <- "C"        # Assign "C"
      } else {
        data4[x, 1] <- "A"        # Assign "A" if Log2 Fold Change is within the threshold
      }
    }
  }
  x <- x + 1
}

# Check the counts of differentially expressed genes
nrow(data4[data4$group == "A", ])
nrow(data4[data4$group == "B", ])
nrow(data4[data4$group == "C", ])

data4$log2FoldChange <- data4$log2FoldChange
#data4$padj <- -log10(data4$padj) # This line was commented out in the original, so keeping it commented.

# Convert 'group' column to factor type
data4$group <- factor(data4$group, levels = c("A", "B", "C"))

# Volcano plot
library(ggplot2) # Load ggplot2 library
library(reshape2) # Load reshape2 library
library(ggrepel) # Load ggrepel library

x.lab = expression(paste("",{Log[2]}," ","Fold Change", sep="")) # X-axis label
y.lab = expression(paste("-",{Log[10]}," ","(adjusted P-value)", sep="")) # Y-axis label

ggplot(data3, aes(x = log2FoldChange, y = padj, color = group, fill = group)) + # Define data frame for plotting; define data for x and y axes; create a scatterplot object.
  #geom_point() + # This line was commented out in the original.
  geom_point(size = 3.5, shape = 21) + # Define data point style.

  ggtitle("Volcano Plot FECD vs Control") + # Define title.

  labs(x = x.lab, y = y.lab) + # Define labels for x and y axes.

  scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, by = 1.0)) + # Define x limits, add ticks.
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 44, by = 2.0)) + # Define y limits, add ticks.

  theme(
    plot.title = element_text(family = "Arial", size = 11, hjust = 0), # Title size and font.
    plot.subtitle = element_text(family = "Arial", size = 11), # Subtitle size and font.
    axis.text = element_text(family = "Arial", size = 10, colour = "black"), # Size and font of x and y values.
    axis.title = element_text(family = "Arial", size = 10), # Size and font of x and y axes.
    panel.border = element_rect(colour = "black", fill = NA, size = 1), # Black border around the plot area.
    axis.ticks = element_line(colour = "black", size = 1), # Style of x and y ticks.
    legend.position = "none" # Remove legend.
  ) +

  geom_hline(yintercept = 1.30103, colour = "black", linetype = "dashed", size = 0.5) + # Horizontal significance cut-off line (P-value = 0.05).
  geom_vline(xintercept = lfc.threshold, colour = "black", linetype = "dashed", size = 0.5) + # Vertical significance cut-off line (+) for Log2 Fold Change.
  geom_vline (xintercept = -lfc.threshold, colour = "black", linetype = "dashed", size = 0.5) + # Vertical significance cut-off line (-) for Log2 Fold Change.

  annotate("rect", xmin = lfc.threshold, xmax = Inf, ymin = 1.30103, ymax = Inf, alpha = .15) + # Transparent rectangle for significant upregulated region.
  annotate("rect", xmin = -lfc.threshold, xmax = -Inf, ymin = 1.30103, ymax = Inf, alpha = .15) + # Transparent rectangle for significant downregulated region.
  #annotate("text",family = "Arial", x = 3.65, y = 1.45, label = "P-value = 0.05", size = 4, fontface = "plain") + # Label for horizontal cut-off line.
  #annotate("text",family = "Arial", x = 0.39, y = 6.4, label = expression(paste("",{Log[2]}," ","Fold Change = 0.5", sep="")),
  #size = 4, fontface = "plain", srt = 90) + # Label for vertical cut-off line.
  #annotate("text",family = "Arial", x = -0.6, y = 6.4, label = expression(paste("",{Log[2]}," ","Fold Change = -0.5", sep="")),
  #size = 4, fontface = "plain", srt = 90) + # Label for vertical cut-off line.

  scale_colour_manual(values = c("A" = "gray29","B" = "darkred","C" = "blue")) + # Manual color scale for point borders.
  scale_fill_manual(values = c("A" = "gray54","B" = "red","C" = "dodgerblue2")) # Manual fill color scale for points.


h <- geom_label_repel(aes(label = ifelse(group == "B", as.character(Entrez_Gene_ID),'')), # Add repulsive labels for group B (upregulated).
                      label.size = 0.05, # Label border size.
                      nudge_x = 1.0, direction = "y", # Nudge labels along x-axis, prefer y-direction.
                      arrow = arrow(length = unit (0.01, "npc"), type = "closed", ends = "last", angle = 15), # Add arrows from labels to points.
                      fill = "lightskyblue") # Fill color for labels.
# g + h # This line was commented out in the original.

# geom_text(size=1, aes(label = ifelse(group == "B", as.character(Entrez_Gene_ID), "")), hjust = 0, vjust = -0.25) # Add identifiers to a subset of data points.

library(ggplot2) # Load ggplot2 library
library(ggrepel) # Load ggrepel library

# Ensure column names are appropriate before creating the plot
# Assuming 'data4' has 'log2FoldChange', 'padj', and 'group' columns
# Create the volcano plot
g <- ggplot(data4, aes(x = log2FoldChange, y = padj, color = group, fill = group)) +
  geom_point(size = 3.5, shape = 21) +
  ggtitle("Volcano Plot FECD vs Control") +
  labs(x = "Log2 Fold Change", y = "-Log10 adjusted P-value") +

  scale_x_continuous(limits = c(-12, 12), breaks = seq(-12, 12, by = 2.0)) + # Define x limits and add ticks.
  scale_y_continuous(limits = c(0,5), breaks = seq(0, 5, by = 1.0)) + # Define y limits and add ticks.

  theme(
    plot.title = element_text(family = "Arial", size = 11, hjust = 0), # Title size and font.
    plot.subtitle = element_text(family = "Arial", size = 11), # Subtitle size and font.
    axis.text = element_text(family = "Arial", size = 10, colour = "black"), # Size and font of x and y values.
    axis.title = element_text(family = "Arial", size = 10), # Size and font of x and y axes.
    panel.border = element_rect(colour = "black", fill = NA, size = 1), # Black border around the plot area.
    axis.ticks = element_line(colour = "black", size = 1), # Style of x and y ticks.
    legend.position = "none" # Remove legend.
  ) +

  geom_hline(yintercept = 1.30103, colour = "black", linetype = "dashed", size = 0.5) + # Horizontal significance cut-off line.
  geom_vline(xintercept = lfc.threshold, colour = "black", linetype = "dashed", size = 0.5) + # Vertical significance cut-off line (+).
  geom_vline (xintercept = -lfc.threshold, colour = "black", linetype = "dashed", size = 0.5) + # Vertical significance cut-off line (-).

  annotate("rect", xmin = lfc.threshold, xmax = Inf, ymin = 1.30103, ymax = Inf, alpha = .15) + # Transparent rectangle for significant upregulated region.
  annotate("rect", xmin = -lfc.threshold, xmax = -Inf, ymin = 1.30103, ymax = Inf, alpha = .15) + # Transparent rectangle for significant downregulated region.
  theme(legend.position = "none") + # Ensure legend is removed.
  scale_color_manual(values = c("A" = "grey", "B" = "red", "C" = "blue")) + # Manual color scale for point outlines.
  scale_fill_manual(values = c("A" = "grey", "B" = "red", "C" = "blue")) # Manual fill color scale for points.

# Display the plot
print(g)
