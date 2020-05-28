library(plyr)
library(edgeR)
library(limma)

preprocess <- function(KO_path) {
  files <- list.files(path=KO_path, pattern="*.txt", full.names=TRUE, recursive=FALSE)
  freq_dfs <- lapply(files, function(x) {
    KO_df <- as.data.frame(read.table(x, header=FALSE, col.names="KO", stringsAsFactors=FALSE))
    new_KO <- plyr::count(KO_df, "KO")
    colnames(new_KO)[2] <- basename(x)
    row.names(new_KO) <- new_KO$KO
    return(new_KO)
  })
  combined <- as.data.frame(freq_dfs[[1]]$KO)
  colnames(combined) <- c("KO")
  for (df in freq_dfs) {
    combined <- merge(combined, df, by=c("KO"), all=TRUE)
  }
  row.names(combined) <- combined$KO
  final_KO <- subset(combined, select=-c(KO))
  colnames(final_KO) <- c("A1HF", "A1LF", "A2HF", "A2LF", "B1HF", "B1LF", "B2HF", "B2LF", "C1HF", "C1LF", "C2HF", "C2LF")
  final_KO[is.na(final_KO)] <- 0
  final_KO <- as.matrix(final_KO)
  return(final_KO)
}

analyze <- function(data, design) {
  y <- DGEList(counts=data,)
  keep <- filterByExpr(y, design=design, min.total.count=100)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- edgeR::calcNormFactors(y)
  y <- estimateDisp(y, design=Design)
  fit <- glmFit(y, design=Design)
  results <- glmLRT(fit, coef=4)
  return(topTags(results, n=Inf))
}

#Create design matrix
Litter <- c(rep("A", times=4),
            rep("B", times=4),
            rep("C", times=4))
Diet <- rep(c("High", "Low"), times=6)
coldata <- cbind(Litter, Diet)
Design <- model.matrix(~0 + Litter + Diet)

#Preprocess KO data and count frequencies of each KO identifier
final_KO <- preprocess(KO_path='./Metagenome/KO_counts')

#Test for enrichment of KO identifiers
results <- analyze(final_KO, design=Design)

#Write analysis outputs
write.table(results, 
            file='/Users/kill337/analysis/manuscript/analysis/local_data/metagenome/output/KO_analysis.txt')








