###########################################################
## Set up directory
###########################################################
dir <- 'data/results'
dir_selected_sig <- 'data/procdata'

# scale function
scaleFun <- function(x){
  
  (x - min(x))/(max(x) - min(x))
  
}
##########################################################
## Load signature data and selected IO biomarkers
##########################################################
## load signature score named geneSig.score
## Note: This is an example data using subset of MOHCCN. Please load computed signature scores from cBioPortal

# MOHCCN data
load(file.path(dir, 'MOHCCN_sig_score.rda'))
sig <- geneSig.score
scaled.sig.data <- t(apply(sig, 1, scaleFun))
rownames(scaled.sig.data) <- rownames(sig)

# load selected signatures from Annals of Oncology paper
selected_signature <- read.csv(file.path(dir_selected_sig, 'meta_signature_logreg.csv')) 

# head(geneSig.score)
sig_name <- intersect(rownames(scaled.sig.data), selected_signature$signature_name)

########################################################################
## define cut-off for reference population 
########################################################################
## consider the 25th and 75th percentile 

sel.patients <- lapply(1:length(sig_name), function(k){
  
  percentiles <- c(0.25, 0.75)
  
  p.high <- 0.75
  p.low <- 0.25
  
  df <- scaled.sig.data[row.names(scaled.sig.data) == sig_name[k],]
  select.high <- names(df[df >= p.high])
  select.low <- names(df[df <= p.low])

  sigScore.high <- round(df[names(df) %in% select.high], 3)
  sigScore.low <- round(df[names(df) %in% select.low], 3)
  
  data.frame(signature = sig_name[k],
             cutoff = c(rep(p.high, length(select.high)),
                        rep(p.low, length(select.low))),
             sigScore = c(sigScore.high, sigScore.low),
             patientid = c(select.high, select.low))
})

sel.patients <- do.call(rbind, sel.patients)
rownames(sel.patients) <- NULL

write.csv(sel.patients, file = file.path(dir, "selected.io.patients.csv"), row.names = FALSE)

