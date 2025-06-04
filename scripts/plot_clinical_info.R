library(ggnewscale)
df2 <- df |> mutate(gender=ifelse(toupper(gender) %in% c("MALE","M"), "Male", ifelse(toupper(gender) %in% c("FEMALE","F"), "Female", "n/a"))) |> 
  mutate(tumor_size_category=case_when(startsWith(tumor_size,"T") ~ tumor_size,
                                       tumor_size == "40+55" ~ "T3",
                                       as.numeric(tumor_size)>1 & as.numeric(tumor_size)<=5 ~ "T1a", 
                                       as.numeric(tumor_size)>5 & as.numeric(tumor_size)<=10 ~ "T1b", 
                                       as.numeric(tumor_size)>10 & as.numeric(tumor_size) <=20 ~ "T1c", 
                                       as.numeric(tumor_size)>20 & as.numeric(tumor_size) <=50 ~ "T2",
                                       as.numeric(tumor_size)>50 ~ "T3",
                                       TRUE ~ "n/a"),
         grade_clean=grade_map[grade] |> unlist(), 
         morphology_clean=morphology_map[morph] |> unlist(),
         sample_id=ifelse(startsWith(sample_id,"MAGIC"), magic_id_conversion[sample_id], sample_id),
         morphology_clean=ifelse(morphology_clean=="Ductal_Carcinoma_In_Situ" & grade=="III", "High_Grade_DCIS", morphology_clean))

hrd_df <- read_tsv("results/supplementary-table4.tsv")
df3 <- hrd_df |> select(sample_id, hr_status, hrd_type, Probability, HRDsum) |> left_join(df2) |> mutate(age_group=case_when(age<40 ~ "<40",
                                  age<50 ~"40-50",
                                  age<60~"50-60",
                                  age<70~"60-70",
                                  TRUE~">70"),
              cohort=case_when(startsWith(sample_id,"FBC")~"Familial",
                               startsWith(sample_id,"TCGA")~ "TCGA",
                               startsWith(sample_id,"MAGIC") ~ "MAGIC",
                               TRUE ~ "Qimprove"),
              hrdetect=ifelse(Probability>=0.7, "HRD","HRP"),
              hrdsum=ifelse(HRDsum>=42, "HRD", "HRP"),
              chord=ifelse(hr_status=="HR_deficient", "HRD","HRP"),
              hrd_type=case_when(hrd_type=="none"~"HRP",
                                 hrd_type=="BRCA1_type"~"BRCA1",
                                 hrd_type=="BRCA2_type"~"BRCA2",
                                 TRUE~"Undetermined"),
              ER=ifelse(ER=="Indeterminate","n/a",ER),
              HER2=ifelse(HER2=="Equivocal","n/a",HER2)) |> 
  select(chord, hrd_type, hrdetect, hrdsum,sample_id, gender, age_group, grade_clean, ER, PR, HER2)

df4 <- df3 |> mutate(tnbc=ifelse(ER=="Negative" & PR=="Negative" & HER2=="Negative", "yes", "no"),
                     tnbc=ifelse(ER=="n/a" | PR=="n/a" | HER2=="n/a", "n/a", tnbc))

sample_id_levels <- df3 |> arrange(chord, hrd_type) |> pull(sample_id)
df3 <- df3 |> pivot_longer(-sample_id,names_to = "field",values_to = "level")
df3$field <- factor(df3$field,levels = c("hrdsum","hrdetect","chord","hrd_type","gender","age_group","grade_clean", "ER","PR","HER2") %>% rev)
df3$sample_id <- factor(df3$sample_id, levels=sample_id_levels)


color_fig <- list(
  hrdsum=c("#0073C2","#EFC000"),
  hrdetect=c("#0073C2","#EFC000"),
  chord=c("#0073C2","#EFC000"),
  hrd_type=c( brewer.pal(3,"Blues") |> rev(),"#EFC000"),
  gender=c("#DEEBF7", "#0073C2"),
  age_group=c("#EEF6FFFF","#B6D6FDFF","#2D6FF0FF", "#1D42B8FF", "#172553FF") |> rev(),
  grade_clean = c( "#ffc857","#e9724c","#c5283d","lightgrey"),
  `ER/PR/HER2`=c("#f77f00","#3772ff", "lightgrey")
)

legend_names <- list(
  chord="HR status",
  hrdetect="HR status",
  hrdsum="HR status",
  hrd_type="HRD type (CHORD)",
  age_group="Age group",
  gender="Gender",
  grade_clean="Grade",
  ER="ER",
  PR="PR",
  HER2="HER2",
)

fact_levels <- list(
  chord=c("HRD", "HRP"),
  hrdsum=c("HRD", "HRP"),
  hrdetect=c("HRD", "HRP"),
  hrd_type=c("BRCA1","BRCA2","Undetermined","HRP"),
  gender=c("Female","Male"),
  age_group=c("<40","40-50","50-60","60-70",">70") |> rev(),
  grade_clean=c("Grade I","Grade II","Grade III","n/a"),
  ER=c("Positive","Negative","n/a"),
  PR=c("Positive","Negative","n/a"),
  HER2=c("Positive","Negative","n/a")
)

generate_row <- function(row_name, df) {
  geom_tile(data=df |> filter(field==row_name) %>% droplevels, 
            mapping=aes(x=sample_id, y=field, fill=factor(level, levels=fact_levels[[row_name]])), 
            color="black")
}

plot <- ggplot() + 
  scale_y_discrete(limits=levels(df3$field), labels=c("HER2","PR","ER","Grade","Age group","Gender","HRD type","CHORD","HRDetect","HRDsum")) +
  theme_minimal(base_size = 14)+
  theme(legend.position="bottom",
        legend.direction  = "vertical",
        panel.grid.major = element_blank(),
        axis.text.x = element_blank()) + 
  labs(x="",y="")


for (i in levels(df3$field)){
  if(i %in% c("HER2","PR","ER")){
    colors <- color_fig[["ER/PR/HER2"]]
    color_name <- "ER/PR/HER2"
  }else{
    colors <-color_fig[[i]]
    color_name <- legend_names[[i]]
  }
  plot <- plot + generate_row(i,df3) +
    scale_fill_manual(values = colors, name=color_name)+ 
    new_scale_fill()
}

plot



get_correlation <- function(df, var_hrd, var_pathology) {
  test_tbl <- df |> select(var_hrd, var_pathology) |> filter(get(var_pathology)!="n/a" & get(var_hrd)!="cannot_be_determined") |> table()
  print(var_pathology)
  print(test_tbl)
  chi_test <- test_tbl |> chisq.test()
  p_value <- chi_test$p.value
  cramer_value <- test_tbl |> cramer_v(test_tbl)
  
  data.frame("x" = var_pathology, "y" = var_hrd, "chi" = p_value, "cramer" = cramer_value)
}

p2 <- c("chord","hrdetect", "hrdsum", "hrd_type") |> map_df(\(y)( c("ER", "PR", "HER2", "grade_clean", "tnbc") |> map_df(\(x) get_correlation(df4, y, x)) ) ) |> 
  mutate(p=case_when(chi<0.0001 ~ "****",
                     chi<0.001 ~ "***",
                     chi<0.01 ~ "**",
                     chi<0.05 ~ "*",
                     TRUE ~ ""),
         font_color=ifelse(cramer<0.3, "white", "black")) |> 
  ggplot(aes(x=x, y=y, fill=cramer)) + geom_tile() + geom_text(aes(label=paste0(round(cramer, 2), p),color=font_color)) + 
  scale_fill_viridis_c(option = "magma") + 
  scale_color_manual(values = c("black","white"), guide="none")  + 
  scale_x_discrete(limits=c("ER","PR","HER2","tnbc", "grade_clean"), labels=c("ER","PR","HER2","TNBC","Grade")) +
  scale_y_discrete(limits=c("hrd_type","chord","hrdetect","hrdsum"), labels=c("HRD type", "CHORD", "HRDetect","HRDsum")) +
  labs(x="", y="", fill="Cramér's V") +
  theme_Publication() +
  theme(legend.text = element_text(size=8.8))


p3 <- hrd_df |> left_join(df4, by="sample_id") |> select(p_hrd, tnbc) |> filter(tnbc!='n/a') |> 
  ggplot(aes(x=tnbc, y=p_hrd)) + geom_boxplot() + geom_point(alpha=0.6) + ggpubr::stat_compare_means(method = "wilcox.test") + theme_Publication() + labs(x="",y="HRD probability (CHORD)") + scale_x_discrete(breaks=c("no", "yes"), labels=c("non-TNBC","TNBC"))

p4 <- hrd_df |> left_join(df4, by="sample_id") |> filter(hr_status=="HR_deficient") |> select(p_BRCA1, tnbc) |> filter(tnbc!='n/a') |> 
  ggplot(aes(x=tnbc, y=p_BRCA1)) + geom_boxplot() + geom_point(alpha=0.6) + ggpubr::stat_compare_means(method = "wilcox.test") + theme_Publication() + labs(x="",y="HRD BRCA1 type probability (CHORD)") + scale_x_discrete(breaks=c("no", "yes"), labels=c("non-TNBC","TNBC"))


plot/(p2+p3+p4 + plot_layout(widths = c(0.7,0.15,0.15))) + plot_layout(heights=c(0.65, 0.35))
