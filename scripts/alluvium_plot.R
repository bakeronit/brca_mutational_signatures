allu_df <-  result_df |> select(cohort, A, B, C, D, C, E, F, sample_id, updated_germline_classification) |>  
  mutate(has_A=ifelse(A>0, "yes", "no") |> factor(levels=c("yes","no")),
         has_B=ifelse(B>0, "yes", "no") |> factor(levels=c("yes","no")),
         has_C=ifelse(C>0, "yes", "no") |> factor(levels=c("yes","no")),
         has_D=ifelse(D>0, "yes", "no") |> factor(levels=c("yes","no")),
         has_E=ifelse(E>0, "yes", "no") |> factor(levels=c("yes","no")),
         has_F=ifelse(F>0, "yes", "no") |> factor(levels=c("yes","no"))) |> 
  select(cohort, has_A, has_B, has_C, has_D, has_C, has_E, has_F, sample_id, updated_germline_classification) |> 
  mutate(germline_classification=ifelse(startsWith(updated_germline_classification,"additional"),"additional",updated_germline_classification)) |> 
  select(-updated_germline_classification) |> 
  group_by(cohort, germline_classification) |> mutate(n=n()) 


need_repel <- c("Q-IMPROvE", "positive_brca1","positive_brca2", "additional")
group_labeller2 <- c("BRCA1/2 VUS","Exclude","Negative BRCA1/2", "Postive BRCA1", "Positive BRCA2", stringr::str_wrap("Familial Breast", 5),"MAGIC","TCGA-BRCA","Q-IMPROvE", "Yes","No")
names(group_labeller2) <- c("additional","exclude","negative_brca1/2","positive_brca1","positive_brca2","Familial Breast","MAGIC","TCGA","Q-IMPROvE", "yes","no")



allu_df$cohort <- factor(allu_df$cohort, levels = cohort_levels)
allu_df$germline_classification <- factor(allu_df$germline_classification, levels=c("positive_brca1","positive_brca2","exclude","additional","negative_brca1/2"))
group_color_mapping <- c("additional"="#0099B4", "exclude" = "#FDAF91", "negative_brca1/2"="#ADB6B6", "positive_brca1" = "#00468B", "positive_brca2"="#ED0000")


p_allu <- ggplot(allu_df |> arrange(cohort),aes(axis6 = cohort, axis5 = has_A, axis4 = has_D, axis3= has_E, axis2=has_B, axis1=germline_classification,
                                                y = n)) +
  scale_x_discrete(limits = c("cohort", "has_A", "has_D", "has_E", "has_B","germline_classification") |> rev() , 
                   labels=c(stringr::str_wrap("Germline classification", 5),"Carries B?","Carries E?","Carries D?","Carries A?","Cohort"), expand = c(.15,.05)) +
  scale_y_continuous(expand = c(0.01, 0.0001)) +  
  geom_alluvium(aes(fill=germline_classification),reverse = FALSE) + scale_fill_manual(values = group_color_mapping) +
  geom_stratum(width = 0.4, reverse = FALSE, aes(fill=ifelse(after_stat(stratum) %in% names(group_color_mapping), as.character(after_stat(stratum)), NA))) + 
  scale_fill_manual(values=group_color_mapping, na.value = "white") +
  geom_text(stat = "stratum", aes(label=ifelse(after_stat(stratum) %in%c("yes","no"), group_labeller2[as.character(after_stat(stratum))], "")), size=2.8, reverse = FALSE)+
  geom_text(stat = "stratum", aes(label=ifelse(after_stat(stratum) %notin% c(need_repel, "yes","no"), group_labeller2[as.character(after_stat(stratum))], "")), size=3.6, reverse = FALSE)+
  geom_text_repel(stat = "stratum", aes(label = ifelse(after_stat(stratum) %in% need_repel & after_stat(x) == 6, as.character(after_stat(stratum)), "")), size = 3.6,
                  direction = "y",
                  nudge_x = .35,
                  reverse = FALSE, segment.color="darkgrey") +
  geom_text_repel(stat = "stratum", aes(label = ifelse(after_stat(stratum) %in% need_repel & after_stat(x) == 1, group_labeller2[as.character(after_stat(stratum))], "")), size = 3.6,
                  direction = "y",
                  nudge_x = -.45,
                  reverse = FALSE, segment.color="darkgrey") +
  theme_void() + coord_flip() + 
  theme(legend.position = "none",
        axis.text.y= element_text())


p_allu_bar<- result_df |> 
  mutate(germline_classification=ifelse(startsWith(updated_germline_classification,"additional"),"additional",updated_germline_classification)) |> 
  count(cohort,germline_classification) |> 
  complete(cohort,germline_classification) |> 
  mutate(n = replace_na(n, 0), cohort=factor(cohort, levels=cohort_levels)) |> 
  mutate(germline_classification = factor(germline_classification, levels=c("positive_brca1","positive_brca2","additional","exclude","negative_brca1/2"))) |> 
  ggplot(aes(x = cohort, y = n, fill = germline_classification)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = n),position = position_dodge(width = 0.7),vjust = -0.5, size=2.5) +
  theme_Publication() +
  scale_x_discrete(limits=cohort_levels, labels=c("Familial Breast","TCGA-BRCA","MAGIC","Q-IMPROvE")) +
  scale_fill_manual( labels=group_labeller2,values = group_color_mapping) + 
  labs(x = "", y = "Count", fill = "") +
  scale_y_continuous(expand = expansion(mult = c(.2, 0.16))) +
  theme(axis.title.y=element_text(vjust=-12))


p_allu /p_allu_bar + plot_layout(heights = c(0.88,0.12))

ggsave("alluvium.pdf", height = 9, width = 8.6)
