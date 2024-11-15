library(signature.tools.lib)
## define cosmic
cosmicv2_h19  <- read_table("data/COSMIC/COSMIC_v2_SBS_GRCh37.txt") |> column_to_rownames("Type")
cosmicv2_h19 <- cosmicv2_h19[match(rownames(COSMIC30_subs_signatures), rownames(cosmicv2_h19)),]
cosmicv3_h19  <- read_table("data/COSMIC/COSMIC_v3.3.1_SBS_GRCh37.txt") |> column_to_rownames("Type")
cosmicv3_h19 <- cosmicv3_h19[match(rownames(COSMIC30_subs_signatures), rownames(cosmicv3_h19)),]

cosmicv2_h38  <- read_table("data/COSMIC/COSMIC_v2_SBS_GRCh38.txt") |> column_to_rownames("Type")
cosmicv2_h38 <- cosmicv2_h38[match(rownames(COSMIC30_subs_signatures), rownames(cosmicv2_h38)),]
cosmicv3_h38  <- read_table("data/COSMIC/COSMIC_v3.3.1_SBS_GRCh38.txt") |> column_to_rownames("Type")
cosmicv3_h38 <- cosmicv3_h38[match(rownames(COSMIC30_subs_signatures), rownames(cosmicv3_h38)),]

sig_levels <-c("Signature_3","Signature_8",paste("Signature",c(1,2,5,12,13,16,17,18,26,28,30),sep = "_"),"unassigned")
sig_colors <-c(palette.colors(palette = "Set1"), palette.colors(palette = "set3"), palette.colors(palette = "Set2")) |> 
  unique() |> 
  head(length(sig_levels)-1)
sig_color_mapping <- setNames(c(sig_colors,"grey"), sig_levels)

sigsToUse <- c(1,2,3,5,6,8,13,17,18,20,26,30)
Zainal_SBS_fit_vcf <- function(sample, snv_file, signature_to_fit, hg="hg38") {
  print(paste0("Converting vcf ", snv_file, " to SNV catalogue for sample ", sample))
  mut_cat <- vcfToSNVcatalogue(snv_file, genome.v = hg)[["catalogue"]]
  colnames(mut_cat) <- sample
  subs_fit_res <- Fit(catalogues = mut_cat,
                      signatures = signature_to_fit,
                      useBootstrap = TRUE,
                      nboot = 1000,
                      nparallel = 4)
  return(subs_fit_res$exposures |> as.data.frame())
}

Zainal_SBS_fit_from_matrix <- function(sample, matrix, signature_to_fit) {
  mut_cat <- matrix[sample]
  subs_fit_res <- Fit(catalogues = mut_cat,
                      signatures = signature_to_fit,
                      useBootstrap = TRUE,
                      nboot = 1000,
                      nparallel = 4)
  return(subs_fit_res$exposures |> as.data.frame())
}

