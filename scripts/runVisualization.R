#------- Description ----------------------------------------------------------
# This script generates boxplots of scaled gene signature scores for selected 
# patients in the MOHCCN study. It highlights individual patient scores overlaid 
# on the population distribution of each immune-related signature.
#
# Key steps:
#   - Load gene signature scores and scale them to a 0–1 range
#   - Load a list of selected patients and their corresponding signatures
#   - Identify overlapping signatures between selected patients and the full dataset
#   - For each signature, create a boxplot showing:
#       • Distribution of signature scores across all patients
#       • Highlighted points for specific patients of interest
#       • Dashed lines indicating fixed cutoffs at 0.25 and 0.75 (Q1 and Q3)
#
# Output:
#   - A PDF file for each signature boxplot, saved to `data/results/Fig/boxplot/`
#
# Assumptions:
#   - Signature score matrix is named 'MOHCCN_sig_score.rda'
#   - Patient metadata is in 'selected.io.patients.csv'
#   - Scores are already scaled or are scaled within the script
#   - Percentile cutoffs are fixed at 0.25 and 0.75 (not dynamic quantiles)
#------------------------------------------------------------------------------
############################################################
## Load libraries
############################################################
library(ggplot2)

############################################################
## Set up directory
############################################################
dir <- 'data/results'

scaleFun <- function(x){
  
  (x - min(x))/(max(x) - min(x))
  
}
##########################################################
## load signature data and selected IO biomarkers
##########################################################
# MOHCCN signature score data
load(file.path(dir, 'MOHCCN_sig_score.rda'))
sig <- geneSig.score
scaled.sig.data <- t(apply(sig, 1, scaleFun))
rownames(scaled.sig.data) <- rownames(sig)

# MOHCCN selected patients 
res <- read.csv(file.path(dir, 'selected.io.patients.csv'))

sig_name <- intersect(unique(res$signature), rownames(scaled.sig.data))

############################################################
## Visualize
###########################################################

for(k in 1:length(sig_name)){
  
  sig.selected <- data.frame(score = scaled.sig.data[sig_name[k],])
  p.high <- 0.75
  p.low <- 0.25
  
  ## specific patient
  df <- res[res$signature == sig_name[k], ]
  specificScore <- data.frame(score = df[, 'sigScore'])  
  
  pdf(file = file.path(dir, 'Fig/boxplot', paste('boxplot', paste(sig_name[k], ".pdf", sep=""), sep="_")), 
      width = 5, height = 3)
  
  p <- ggplot(sig.selected, aes(x = score, y = "")) +
    geom_boxplot(fill = "#E1E1E1FF", color = "#273649FF", width = 0.05) +
    geom_point(data = specificScore, aes(x = score, y = ""), 
               color = "#8E2322FF", size = 2, shape = 18) + 
    geom_vline(xintercept = p.low, linetype = "dashed", color = "#ABB2A5FF", size = 0.3) +
    geom_vline(xintercept = p.high, linetype = "dashed", color = "#ABB2A5FF", size = 0.3) +
    annotate("text", x = p.low, y = 1.2, label = paste0("Q1: ", round(p.low, 2)), 
             color = "#273649FF", angle = 90, vjust = -0.5, size = 3) +
    annotate("text", x = p.high, y = 1.2, label = paste0("Q3: ", round(p.high, 2)), 
             color = "#273649FF", angle = 90, vjust = -0.5, size = 3) +
    labs(
      title = "",
      x = "",
      y = NULL
    ) +
    theme_void() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=12, face="bold"),
      legend.position="none",
      legend.text = element_text(size = 6, face="bold"),
      legend.title = element_blank()
    ) 
  
  print(p)
  
  dev.off()
  
  
}
