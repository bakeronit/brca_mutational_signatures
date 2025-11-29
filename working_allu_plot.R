 p_allu <- ggplot(allu_df |> arrange(cohort),aes(axis6 = cohort, axis5 = has_A, axis4 = has_D, axis3= has_E, axis2=has_B, axis1=germline_classification,
                                                y = n)) +
  scale_x_discrete(limits = c("cohort", "has_A", "has_D", "has_E", "has_B","germline_classification") |> rev() , 
                   labels=c(stringr::str_wrap("Individual Classification", 5),"Carries D","Carries C","Carries B","Carries A","Cohort"), expand = c(.15,.05)) +
  scale_y_continuous(expand = c(0.01, 0.0001)) +  
  geom_alluvium(aes(fill=germline_classification),reverse = FALSE) + scale_fill_manual(values = group_color_mapping) +
  geom_stratum(width = 0.4, reverse = FALSE, aes(fill=ifelse(after_stat(stratum) %in% names(group_color_mapping), as.character(after_stat(stratum)), NA))) + 
  scale_fill_manual(values=group_color_mapping, na.value = "white") +
  geom_text(stat = "stratum", aes(label=ifelse(after_stat(stratum) %in%c("yes","no"), alluvial_labeller[as.character(after_stat(stratum))], "")), size=2.8, reverse = FALSE)+
  geom_text(stat = "stratum", aes(label=ifelse(after_stat(stratum) %notin% c(need_repel, "yes","no"), alluvial_labeller[as.character(after_stat(stratum))], "")), size=3.6, reverse = FALSE)+
  geom_text_repel(stat = "stratum", aes(label = ifelse(after_stat(stratum) %in% need_repel & after_stat(x) == 6, as.character(after_stat(stratum)), "")), size = 3.6,
                  direction = "y",
                  nudge_x = .35,
                  reverse = FALSE, segment.color="darkgrey") +
  geom_text_repel(stat = "stratum", aes(label = ifelse(after_stat(stratum) %in% need_repel & after_stat(x) == 1, alluvial_labeller[as.character(after_stat(stratum))], "")), size = 3.6,
                  direction = "y",
                  nudge_x = -.45,
                  reverse = FALSE, segment.color="darkgrey") +
  theme_void() + coord_flip() + labs(title="Classification of 350 individuals based on variants") +
  theme(legend.position = "none",
        axis.text.y= element_text(),
        plot.title = element_text(hjust = 0.5, vjust=-10))
 
 
p_allu <- ggplot(allu_df |> arrange(cohort),aes(axis6 = cohort, axis5 = has_D, axis4 = has_A, axis3= has_E, axis2=has_B, axis1=germline_classification,
                                                 y = n)) +
   scale_x_discrete(limits = c("cohort", "has_A", "has_D", "has_E", "has_B","germline_classification") |> rev() , 
                    labels=c(stringr::str_wrap("Individual Classification", 5),"Carries D","Carries C","Carries B","Carries A","Cohort"), expand = c(.15,.05)) +
   scale_y_continuous(expand = c(0.01, 0.0001)) +  
   geom_alluvium(aes(fill=germline_classification),reverse = FALSE) + scale_fill_manual(values = group_color_mapping) +
   geom_stratum(width = 0.4, reverse = FALSE, aes(fill=ifelse(after_stat(stratum) %in% names(group_color_mapping), as.character(after_stat(stratum)), NA))) + 
   scale_fill_manual(values=group_color_mapping, na.value = "white") +
   geom_text(stat = "stratum", aes(label=ifelse(after_stat(stratum) %in%c("yes","no"), alluvial_labeller[as.character(after_stat(stratum))], "")), size=2.8, reverse = FALSE)+
   geom_text(stat = "stratum", aes(label=ifelse(after_stat(stratum) %notin% c(need_repel, "yes","no"), alluvial_labeller[as.character(after_stat(stratum))], "")), size=3.6, reverse = FALSE)+
   geom_text_repel(stat = "stratum", aes(label = ifelse(after_stat(stratum) %in% need_repel & after_stat(x) == 6, as.character(after_stat(stratum)), "")), size = 3.6,
                   direction = "y",
                   nudge_x = .35,
                   reverse = FALSE, segment.color="darkgrey") +
   geom_text_repel(stat = "stratum", aes(label = ifelse(after_stat(stratum) %in% need_repel & after_stat(x) == 1, alluvial_labeller[as.character(after_stat(stratum))], "")), size = 3.6,
                   direction = "y",
                   nudge_x = -.45,
                   reverse = FALSE, segment.color="darkgrey") +
   theme_void() + coord_flip() + labs(title="Classification of 350 individuals based on variants") +
   theme(legend.position = "none",
         axis.text.y= element_text(),
         plot.title = element_text(hjust = 0.5, vjust=-10))

 p_allu / p_germline_classification_by_cohort + plot_layout(heights = c(0.88,0.12))
 ggsave("figures/Figure1b_rev.pdf", height = 9, width = 8.6, dpi = 320)
 
 result_df$updated_germline_classification |> table()
 