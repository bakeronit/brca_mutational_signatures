library(tidyverse)
library(signature.tools.lib)
library(BSgenome.Hsapiens.UCSC.hg19)

#sigsToUse_v2 <- c("Signature_1","Signature_2","Signature_3","Signature_5","Signature_6","Signature_8","Signature_13","Signature_17","Signature_18","Signature_20","Signature_26","Signature_30") ## cosmic v2 signature related to breast cancer cohort
sigsToUse_v2 <- c("Signature_1","Signature_2","Signature_3","Signature_8","Signature_13","Signature_16","Signature_18","Signature_30")
sigsToUse_v3 <- c("SBS1","SBS2","SBS3","SBS5","SBS8","SBS13","SBS18","SBS30","SBS39")

cosmicv2_signatures  <- read_table("data/COSMIC/COSMIC_v2_SBS_GRCh37.txt") |> column_to_rownames("Type")
cosmicv2_signatures <- cosmicv2_signatures[match(rownames(COSMIC30_subs_signatures), rownames(cosmicv2_signatures)),]

cosmicv3_signatures  <- read_table("data/COSMIC/COSMIC_v3.3.1_SBS_GRCh37.txt") |> column_to_rownames("Type")
cosmicv3_signatures <- cosmicv3_signatures[match(rownames(COSMIC30_subs_signatures), rownames(cosmicv3_signatures)),]


# function to load a single vcf, calculate SNV catalogue and fit to reference signatures ---------------------
fit_each_snv_category <- function(fn_snv, sample_id, signatures_to_fit, genome) {
  print(paste0("Converting vcf ",fn_snv, " to SNV catalogue for sample ", sample_id))
  res <- vcfToSNVcatalogue(fn_snv, genome.v = genome)
  colnames(res$catalogue) <- sample_id
  
  print(paste0("Fitting snv catalogue for sample ",sample_id))
  subs_fit_res <- Fit(catalogues = res$catalogue,
                      signatures = signatures_to_fit,
                      useBootstrap = TRUE,
                      nboot = 100,
                      nparallel = 4)
  
  subs_fit_res$exposures |> as_tibble()
}

# New tricks I learned:
# using |> is faster than %>%  :)
# \ is a shortcut to use anonymous function then can use |> like %>%  for multiple placement of .
# imap() can make use of element(.x) and index(.y). use map won't be able to access the name of each element in the function.
# use list_rbind() to combine the list of data frames into a single data frame, and names of list become a column
# Load all vcf files, and get SNV signature fitted results for both cosmic v2 and v3.3.1 ------------------------------------------

snv_sigv2_df <- list.files("data/GRCh37/Nones2019", pattern = "*.snv.vcf.gz$",full.names = TRUE,recursive = TRUE) |>
  ( \(x) set_names(x, str_remove(basename(x), ".snv.vcf.gz")) )() |>  
  imap(~fit_each_snv_category(fn_snv = .x, sample_id = .y, signatures_to_fit = cosmicv2_signatures[sigsToUse_v2], genome="hg19")) |> list_rbind(names_to="sample")

snv_sigv3_df <- list.files("data/GRCh37/Nones2019", pattern = "*.snv.vcf.gz$",full.names = TRUE,recursive = TRUE) |>
  ( \(x) set_names(x, str_remove(basename(x), ".snv.vcf.gz")) )() |>  
  imap(~fit_each_snv_category(fn_snv = .x, sample_id = .y, signatures_to_fit = cosmicv3_signatures[sigsToUse_v3], genome="hg19")) |> list_rbind(names_to="sample")

save(snv_sigv2_df, snv_sigv3_df, file="data/GRCh37_sub_sig.rda")


# Initialize an empty matrix for the HRDetect pipeline ----------------------
col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")

input_matrix_v2 <- matrix(data=NA, nrow = nrow(snv_sigv2_df), ncol = length(col_hrdetect), dimnames = list(snv_sigv2_df$sample,col_hrdetect))
input_matrix_v2[,"SNV3"] <- snv_sigv2_df$Signature_3
input_matrix_v2[,"SNV8"] <- snv_sigv2_df$Signature_8

input_matrix_v3 <- matrix(data=NA, nrow = nrow(snv_sigv3_df), ncol = length(col_hrdetect), dimnames = list(snv_sigv3_df$sample,col_hrdetect))
input_matrix_v3[,"SNV3"] <- snv_sigv3_df$SBS3
input_matrix_v3[,"SNV8"] <- snv_sigv3_df$SBS8

# Load indel, cnv, and sv results -----------------------------------------
indel_files <- list.files("data/GRCh37/Nones2019",pattern = "*.indel.vcf.gz$", full.names = TRUE, recursive = TRUE) |> (\(x) set_names(x, str_remove(basename(x), ".indel.vcf.gz") ))()
cnv_files <- list.files("data/GRCh37/Nones2019", pattern = "*.cnv.txt$", full.names = TRUE, recursive = TRUE) |> (\(x) set_names(x, str_remove(basename(x), ".cnv.txt") ))()
sv_files <- list.files("data/GRCh37/Nones2019",pattern = "*.sv.bedpe$", full.names = TRUE, recursive = TRUE) |> (\(x) set_names(x, str_remove(basename(x), ".sv.bedpe") ))()

res_v2 <- HRDetect_pipeline(input_matrix_v2,
                            genome.v = "hg19",
                            SV_bedpe_files = sv_files,
                            Indels_vcf_files = indel_files,
                            CNV_tab_files = cnv_files,
                            nparallel = 4)

res_v3 <- HRDetect_pipeline(input_matrix_v3,
                            genome.v = "hg19",
                            SV_bedpe_files = sv_files,
                            Indels_vcf_files = indel_files,
                            CNV_tab_files = cnv_files,
                            nparallel = 4)

save(res_v2, res_v3, file="data/GRCh37_hrd_results.rda")
