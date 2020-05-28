library(NormalizeMets)
library(dplyr)
library(imputeLCMD)
library(EnhancedVolcano)

set.seed(1)

preprocess <- function(metabolites_path, annotations_path, annotations=FALSE) {
  metabolites <- data.frame(read.csv(metabolites_path, header=TRUE, row.names='Sample', stringsAsFactors=FALSE))
  log_metabolites <- log2(metabolites)
  log_metabolites[log_metabolites=='-Inf'] <- NaN
  imputed_metabolites <- impute.MinProb(log_metabolites)
  if (annotations==TRUE) {
    annotations <- data.frame(read.csv(annotations_path, header=TRUE, row.names="Metabolite", stringsAsFactors=FALSE))
    return(list(t(imputed_metabolites), annotations))
  }
  else {
    return(list(t(imputed_metabolites)))
  }
}

analyze <- function(design_matrix, data, low_select="dietLow", annotations=FALSE) {
  results <- LinearModelFit(data[[1]], factormat = Design, ruv2=FALSE, moderated=TRUE, padjmethod='BH')
  p_values <- as.data.frame(results$adj.p.value)
  select_low <- subset(p_values, select=low_select)
  select_low$coeffs <- subset(results$coefficients, select=low_select)
  colnames(select_low) <- c("p_value", "coefficient")
  if (annotations==TRUE) {
    merged <- merge(select_low, data[[2]], by="row.names")
    row.names(merged) <- merged$Row.names
    merged <- subset(merged, select = -c(Row.names))
    return(merged)
  }
  else {
    return(select_low)
  }
}

#Import, log2 transform, and impute metabolomics data for both polar and volatile data.
polar_metabolites <- preprocess(metabolites_path="./Metabolomics/inputs/filtered_polar_metabolites.csv", 
                                annotations_path="./Metabolomics/inputs/filtered_polar_annotations.csv",
                                annotations=TRUE)
volatile_metabolites <- preprocess(metabolites_path="./Metabolomics/inputs/filtered_volatile_metabolites.csv", 
                                  annotations_path=NULL,
                                  annotations=FALSE)

#Create design matrix
diet <- c(rep("High", times=6),
          rep("Low", times=6))
litter <- rep(c(rep("A", times=2),
                rep("B", times=2),
                rep("C", times=2)), times=2)
Design <- model.matrix(~0 + litter + diet)

#Compare high and low dietary fiber conditions with LinearModelFit
polar_metabolites_analysis <- analyze(Design, polar_metabolites, annotations=TRUE)
volatile_metabolites_analysis <- analyze(Design, volatile_metabolites, annotations=FALSE)



#Write analysis outputs
Metabolite <- row.names(polar_metabolites_analysis)
polar_table <- cbind(Metabolite, polar_metabolites_analysis)

Metabolite <- row.names(volatile_metabolites_analysis)
volatile_table <- cbind(Metabolite, volatile_metabolites_analysis)

write.table(polar_table, 
            file="./Metabolomics/outputs/polar_metabolites_analysis.tsv",
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
write.table(volatile_table, 
            file="./Metabolomics/outputs/volatile_metabolites_analysis.tsv",
            quote=FALSE,
            sep="\t",
            row.names=FALSE)


### PLOTTING OF DATA ###

plot_data <- function(analysis_df) {
  color_labels <- rep('grey', nrow(analysis_df))
  #names(color_labels) <- 'Not significant'
  
  color_labels[which(analysis_df$p_value < .05 & analysis_df$coefficient > 0)] <- 'green'
  color_labels[which(analysis_df$p_value < .05 & analysis_df$coefficient < 0)] <- 'blue'
  #names(color_labels[which(color_labels == 'green')]) <- 'LF'
  #names(color_labels[which(color_labels == 'blue')]) <- 'HF'
  names(color_labels)[which(color_labels == 'green')] <- 'LF'
  names(color_labels)[which(color_labels == 'blue')] <- 'HF'
  names(color_labels)[which(color_labels == 'grey')] <- 'Not significant'
  
  #names(color_labels[which(color_labels == 'grey')]) <- 'Not significant'
  
  min_coef <- min(floor(analysis_df$coefficient))
  max_coef <- max(ceiling(analysis_df$coefficient))
  max_coef_all <- max(c(abs(min_coef), abs(max_coef)))
  
  max_p <- max(ceiling(-log10(analysis_df$p_value)))

  volcano_plot <- EnhancedVolcano(analysis_df,
                                  lab = rownames(analysis_df),
                                  title = NULL,
                                  subtitle = NULL,
                                  x = 'coefficient',
                                  y = 'p_value',
                                  #xlab = NULL,
                                  #ylab = NULL,

                                  legendVisible=FALSE,
                                  xlab = bquote(~Log[2]~ "fold-change"),
                                  ylab = bquote(~-Log[10]~italic(P)[adj]),
                                  
                                  #xlab = ' ',
                                  #ylab = ' ',
                                  
                                  #legend=c('Not significant','Log (base 2) fold-change ≥ |2|', "Adjusted P-value ≤ 0.05",
                                  #         'Adjusted P-value ≤ 0.05 & Log (base 2) fold-change ≥ |2|'),
                                  FCcutoff = 0,
                                  #cutoffLineWidth=0,
                                  #boxedlabels = TRUE,
                                  #drawConnectors = TRUE,
                                  #widthConnectors = 0.25,
                                  transcriptPointSize = 5,
                                  transcriptLabSize = 6.5,
                                  
                                  #legendLabSize = 10,
                                  #legendIconSize = 3,
                                  #selectLab =  carbohydrateMetabolitesNames,
                                  #selectLab = selectMetabolites,
                                  pCutoff = .05,
                                  #shapeCustom = keyvals.shape,
                                  #colConnectors = 'black',
                                  #xlim = c(-max_coef_all, max_coef_all),
                                  xlim = c(-5, 5),
                                  #ylim = c(0, max_p),
                                  ylim = c(0, 5),
                                  boxedlabels = TRUE,
                                  drawConnectors = TRUE,
                                  widthConnectors = 0.35,
                                  colConnectors = 'grey30',
                                  #lengthConnectors=unit(1, 'npc')
                                  endsConnectors="last",
                                  #col = c("grey30", "forestgreen", "royalblue", "red2"),
                                  colCustom = color_labels,
                                  colAlpha = .65,
                                  #shadeAlpha = 1,
                                  caption=NULL,
                                  legend=NULL,
                                  #legendPosition = "right"
  )
  
  return(volcano_plot)
  #return(color_labels)
}

volatile_plot <- plot_data(volatile_metabolite_analysis)

volatile_out_path <- "./volatiles.png"
metab_plot_nol <- volatile_plot + theme(legend.position="none")
metab_plot_nol



ggsave(volatile_out_path, width=5, scale=2, dpi=600)

#Polar metabolites by class
assign_class <- function(data) {
  data$Metabolite <- row.names(data)
  uni_anns <- unique(data$Class)
  shapes <- seq(from = 1, to = length(uni_anns), by = 1) #
  anns_df <- data.frame(shapes)
  anns_df$Class <- uni_anns
  z <- merge(anns_df, data, by="Class")
  
  z_l <- list()
  for (shape in uni_anns) {
    spec_data <- data[data$Class==shape, ]
    z_l[[shape]] <- spec_data
  }
  return(z_l)
}

df <- assign_class(polar_metabolites_analysis)

for (ds in names(df)) {
  print(ds)
  metab_plot <- plot_data(df[[ds]]) 
  metab_plot_nol <- metab_plot + theme(legend.position="none")
  out_path <- paste("./polar_", ds, "_labeled.png", sep="")
  ggsave(out_path, width=5, scale=2, dpi=600)
}



#metab_plot_nol <- metab_plot + theme(legend.position="none")















