library(tidyverse)
library(ggmosaic)
result_df <- read_tsv("results/brca_germline_clf_tumor_hrd.add_sv.add_somatic.tsv")
result_df |> select(cohort,sample_id, updated_germline_classification, hr_status) |> 
  filter(updated_germline_classification %in% c("positive_brca1","positive_brca2","negative_brca1/2")) |> 
  mutate(hr_status=ifelse(hr_status=="HR_deficient","HRD","non-HRD"), updated_germline_classification=stringr::str_wrap(group_labeller2[updated_germline_classification],10)) |> 
  ggplot() + geom_mosaic(aes(x=product(updated_germline_classification), fill=hr_status)) + 
  scale_fill_manual(values=c("#0073C2","#EFC000")) + theme_Publication() +
  labs(x="", y="", fill="HR status")


result_df |> select(cohort,sample_id, updated_germline_classification, hrd_type) |> 
  filter(updated_germline_classification %in% c("positive_brca1","positive_brca2","negative_brca1/2"), hrd_type!="cannot_be_determined") |> 
  ggplot() + geom_mosaic(aes(x=product(updated_germline_classification), fill=hrd_type)) + scale_fill_manual(values=c("#0073C2","#EFC000","red")) +
  theme_Publication()



group_color_mapping2 <- c("additional"="#0099B4", "exclude_plp" = "#E5735A","exclude" = "#FDAF91", "negative_brca1/2"="#ADB6B6", "positive_brca1" = "#00468B", "positive_brca2"="#ED0000")
group_labeller2 <- c("VUS BRCA1/2","Exclude P/LP","Exclude VUS","Negative BRCA1/2", "Postive BRCA1", "Positive BRCA2")
names(group_labeller2) <- c("additional","exclude_plp","exclude","negative_brca1/2","positive_brca1","positive_brca2")


exclude_plp_df <- result_df |> select(cohort, sample_id, hrd_type, updated_germline_classification, A, class_A_mutations, pathSV) |> filter(A>0) |> 
  mutate(variant=ifelse(is.na(class_A_mutations), 
                        paste0(str_split_i(pathSV,"_",i=5)," (DEL)"),
                        str_replace(class_A_mutations,"^\\s*(c\\.[^\\[]+)\\[.*?\\]\\(([^)]+)\\).*$","\\2 (\\1)"))
         ) |> select(hrd_type, variant) |> filter(hrd_type!="none") |>  group_by(hrd_type) |> mutate(row=row_number()) |> ungroup() |> 
  mutate(row=ifelse(hrd_type=="BRCA1_type", row+1, row), col=case_when(hrd_type=="BRCA1_type"~11, hrd_type=="BRCA2_type"~8, TRUE~2))

result_df |> select(cohort,sample_id, hrd_type, updated_germline_classification, A, B) |> 
  mutate(germline_classification=ifelse(startsWith(updated_germline_classification,"additional"),"additional",updated_germline_classification)) |> 
  mutate(detail_exclude=ifelse(A>0, "exclude_plp", germline_classification)) |> 
  mutate(detail_exclude=factor(detail_exclude, levels=c("positive_brca1","positive_brca2","exclude_plp","exclude","additional","negative_brca1/2"))) |> 
  group_by(hrd_type, detail_exclude) |> summarise(n=n()) |> ungroup() |> 
  ggplot() + 
  geom_waffle(aes(fill=detail_exclude, values=n),color = "white", size = 1.125, n_rows = 6) + 
  geom_text(data=exclude_plp_df, aes(label=variant, x=col,y=row), hjust = 0, fontface="bold") +
  scale_fill_manual(values=group_color_mapping2, labels=group_labeller2) + scale_x_discrete(expand=c(0.01,0.1)) +
  facet_wrap(~hrd_type, ncol = 1,  strip.position = "left", labeller = labeller(hrd_type=hrdtype_labeller)) + theme_Publication() + labs(x="", y="HR status ",fill="Germline classification") +
  theme( legend.position = "bottom", axis.line.y= element_blank(), axis.ticks = element_blank() ,axis.line.x = element_blank(), axis.text = element_blank()) 
  


layout <-"
AABBB
AABBB
CCCCC
CCCCC
CCCCC
"
pa + pb + pc + plot_layout(design = layout)

pab <- cowplot::plot_grid(pa + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),pb + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), rel_widths = c(0.45, 0.55), align = "hv")
pabc <- cowplot::plot_grid(pab, pc , ncol = 1, rel_heights = c(0.45,0.55))
pabc
