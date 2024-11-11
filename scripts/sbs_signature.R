## a general function to read a single vcf file and output a matrix
## chr pos ref alt
read_vcf <- function(vcf) {
  df <- read_tsv(vcf, comment = "#", col_names = F) |> select(X1, X2, X4, X5)
  colnames(df) <- c("chr","pos","ref","alt")
  if (startsWith(df$chr, "chr") |> all()) {
    df$chr <- str_remove(df$chr, "chr")
  }
  df <- df |> filter(chr %in% c(1:22,"X")) |> mutate(chr=paste0("chr",chr))
  return(df)
}

## from a name vcf files list to a matrix with 
## sample(row) x weight(col) of each signature
run_deconstructSigs <- function(vcf_list, sigs_to_fit, bsgenome) {
  mut_input <- vcf_list |> map(read_vcf) |> list_rbind(names_to = "Sample") |> 
    as.data.frame() |> 
    mut.to.sigs.input(bsg=bsgenome)
  
  results <- data.frame()
  for (sample in rownames(mut_input)){
    result = whichSignatures(tumor.ref = mut_input,
                           signatures.ref = sigs_to_fit,
                           sample.id = sample,
                           contexts.needed = TRUE,
                           tri.counts.method = "default")
    results <- rbind(results, result$weight)
  }
  results
}
