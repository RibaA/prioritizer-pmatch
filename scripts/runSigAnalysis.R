#------------------------------------------------------------------------------
# This script performs immune-related gene signature scoring using expression
# data from cancer immunotherapy patients. It supports multiple signature scoring
# methods including GSVA, ssGSEA, Weighted Mean, and specific published algorithms.
# 
# The script does the following:
#   - Loads raw gene expression data and signature definitions
#   - Applies log2 transformation to the expression matrix (if needed)
#   - Loads signature metadata to determine which method to use per signature
#   - Computes signature scores using the appropriate function:
#       • geneSigGSVA(), geneSigssGSEA(), geneSigMean(), or custom methods
#   - Handles missing values and filters out invalid results
#   - Outputs a final signature score matrix for downstream analysis
#
# Input:
#   - Expression matrix (CSV)
#   - Signature list (.rda files)
#   - Signature metadata (CSV with method annotations)
#
# Output:
#   - Matrix of computed signature scores (RDA file)
#
# Requirements:
#   - R packages: PredictioR, GSVA, dplyr, data.table
#   - Signature scoring functions: geneSigGSVA, geneSigssGSEA, etc.
#------------------------------------------------------------------------------
###########################################################
# Load libraries
###########################################################
library(PredictioR)
library(GSVA)
library(dplyr)
library(data.table)

###########################################################
## Set up directory
###########################################################

dir_in <- 'data/rawdata'
dir_out <- 'data/results'

############################################
## Load data 
############################################
study_icb <- "MOHCCN"
expr <- read.csv(file = file.path(dir_in, "expressions.csv")) # 17741 genes and 452 patients
rownames(expr) <- expr$X ## gene names should be row names, not a column of data frame
expr <- expr[, -1]  
colnames(expr) <- gsub("\\.", "-", colnames(expr))
expr <- log2(expr+1) # log-transformed data

############################################
## Load signature data
############################################
# download signature data from https://github.com/bhklab/SignatureSets
files <- list.files(file.path(dir_in, 'signature')) 

signature <- lapply(1:length(files), function(k){
  
  print(k)
  load(file.path(dir_in, 'signature', files[k]))
  sig

})
names(signature) <- substr(files, 1, nchar(files) - 4)

# signature information
signature_info <- read.csv(file.path(dir_in, "signature_information.csv")) 

##################################################################
## Compute signature score
##################################################################

geneSig.score <- lapply(1:length(signature), function(i){ 
  
  print(paste(i , names(signature)[i], sep="/"))
  sig_name <- names(signature)[i]
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "GSVA"){
    
    geneSig <- geneSigGSVA(dat.icb = expr,
                           sig = signature[[i]],
                           sig.name = sig_name,
                           missing.perc = 0.5,
                           const.int = 1,
                           n.cutoff = 15,
                           sig.perc = 0.8,
                           study = study_icb,
                           gene.annot = "entrez_id")
    
    
    if(sum(!is.na(geneSig)) > 0){
      geneSig <- geneSig[1,]
    }     
    
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Weighted Mean"){
    
    geneSig <- geneSigMean(dat.icb = expr,
                           sig = signature[[i]],
                           sig.name = sig_name,
                           missing.perc = 0.5,
                           const.int = 1,
                           n.cutoff = 15,
                           sig.perc = 0.8,
                           study = study_icb,
                           gene.annot = "entrez_id")
    
  }
  
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "ssGSEA"){
    
    geneSig <- geneSigssGSEA(dat.icb = expr,
                             sig = signature[[i]],
                             sig.name = sig_name,
                             missing.perc = 0.5,
                             const.int = 1,
                             n.cutoff = 15,
                             sig.perc = 0.8,
                             study = study_icb,
                             gene.annot = "entrez_id")
    
    if(sum(!is.na(geneSig)) > 0){
      geneSig <- geneSig[1,]
    }     
    
    
  }
  
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "COX-IS_Bonavita"){
    
    geneSig <- geneSigCOX_IS(dat.icb = expr,
                             sig = signature[[i]],
                             sig.name = sig_name,
                             missing.perc = 0.5,
                             const.int = 1,
                             n.cutoff = 15,
                             sig.perc = 0.8,
                             study = study_icb,
                             gene.annot = "entrez_id")
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPS_Charoentong"){
    
    geneSig <- geneSigIPS(dat.icb = expr,
                          sig = signature[[i]],
                          sig.name = sig_name,
                          missing.perc = 0.5,
                          const.int = 1,
                          n.cutoff = 15,
                          study = study_icb,
                          gene.annot = "entrez_id")
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PredictIO_Bareche"){
    
    geneSig <- geneSigPredictIO(dat.icb = expr,
                                sig = signature[[i]],
                                sig.name = sig_name,
                                missing.perc = 0.5,
                                const.int = 1,
                                n.cutoff = 15,
                                sig.perc = 0.8,
                                study = study_icb,
                                gene.annot = "entrez_id")
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPRES_Hugo"){
    
    geneSig <- geneSigIPRES(dat.icb = expr,
                            sig = signature[[i]],
                            sig.name = sig_name,
                            missing.perc = 0.5,
                            const.int = 1,
                            n.cutoff = 15,
                            sig.perc = 0.8,
                            study = study_icb,
                            gene.annot = "entrez_id")
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PassON_Du"){
    
    geneSig <- geneSigPassON(dat.icb = expr,
                             sig = signature[[i]],
                             sig.name = sig_name,
                             missing.perc = 0.5,
                             const.int = 1,
                             n.cutoff = 15,
                             sig.perc = 0.8,
                             study = study_icb,
                             gene.annot = "entrez_id")
    
  }
  
  if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPSOV_Shen"){
    
    geneSig <- geneSigIPSOV(dat.icb = expr,
                            sig = signature[[i]],
                            sig.name = sig_name,
                            missing.perc = 0.5,
                            const.int = 1,
                            n.cutoff = 15,
                            sig.perc = 0.8,
                            study = study_icb,
                            gene.annot = "entrez_id")
    
  }
  
  
  if(sum(!is.na(geneSig)) > 0){
    
    geneSig <- geneSig
    
  }     
  
  if(sum(!is.na(geneSig)) == 0){
    
    geneSig <- rep(NA, ncol(expr))
    
  }
  
  geneSig
  
})

geneSig.score <- do.call(rbind, geneSig.score)
rownames(geneSig.score) <- names(signature)
remove <- which(is.na(rowSums(geneSig.score)))
if(length(remove) > 0){
  
  geneSig.score <- geneSig.score[-remove, ]
  
}

save(geneSig.score, file=file.path(dir_out, paste(study_icb , "sig_score.rda", sep = "_")))

