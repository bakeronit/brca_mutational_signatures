pe <- plot_df |> 
  ggplot(aes(x=tool)) +
  geom_bar(position="fill",width=0.55,aes(fill=updated_germline_classification)) + facet_grid(~hr_status) + scale_fill_manual(values=group_color_mapping, labels=group_labeller) +
  geom_label(data=plot_df |> filter(hr_status!="non-HRD" | updated_germline_classification %notin% c("positive_brca1","positive_brca2")),
             aes(label=..count.., color=type, fill=updated_germline_classification),
             stat="count",position="fill",size=4,hjust=1.1,show.legend = FALSE) + scale_color_manual(values=c("white","black")) +
  geom_label(data=plot_df |> filter(hr_status=="non-HRD" & updated_germline_classification %in% c("positive_brca1","positive_brca2")),
             aes(label=..count.., color=type, fill=updated_germline_classification),
             stat="count",position="fill",size=4, hjust=-.1,show.legend = FALSE) + 
  geom_text(data=label_df,aes(x=tool,y=0.5, label=paste0("Total=",n)),
            position = position_dodge(0.9),vjust=-1.5, fontface="bold") +
  facet_grid(~hr_status, scales = "free_x") + coord_flip() +
  labs(x="",y="Proportion", fill="Germline classification") +
  theme_Publication() +
  theme(strip.text = element_text(size=10), legend.position = "bottom")

pabcd <- plot_grid(pabc, p1, rel_widths = c(0.2,0.8))
plot_grid(pabcd, pp, rel_heights = c(0.6,0.4), nrow = 2, ncol = 1)

layout <- "
ADDDDDD
ADDDDDD
BDDDDDD
BDDDDDD
CDDDDDD
CDDDDDD
EEEEEEE
EEEEEEE
EEEEEEE
"
pa + pb + pc + p1 + pp + plot_layout(design = layout) + theme(plot.margin = margin(0,0,0,0)) 




#### back up
**The number of samples in each group based on the combined HRD prediction between HRDetect and CHORD**
  
  ```{r, germline_classification_combined,  fig.width=9, fig.height=4.8}
result_df |> filter(hr_status_combined!="Inconsistent") |> 
  ggplot(aes(x=hr_status_combined, fill=hr_status_combined)) +
  geom_bar() + 
  geom_text(stat='count', aes(x=hr_status_combined,label=..count..), vjust=-.2, size=4) +
  facet_grid(~updated_germline_classification, labeller = labeller(updated_germline_classification=group_labeller)) +
  scale_fill_manual(values=c("#0073C2","#EFC000")) +
  labs(x="",y="Count", fill="HR status") + ylim(0,150) +
  theme_Publication() +
  theme(strip.text = element_text(size=10))
```







**The number of samples in each group based on CHORD, HRDetect or HRDsum results**
  
  ```{r, germline_classification_chord,  fig.width=9, fig.height=3.2}
result_df |>
  mutate(hr_status=ifelse(hr_status=="HR_deficient","HRD","non-HRD")) |> 
  ggplot(aes(x=hr_status, fill=hr_status)) +
  geom_bar() + 
  geom_text(stat='count', aes(x=hr_status,label=..count..), vjust=-.2, size=4) +
  facet_grid(~updated_germline_classification, labeller = labeller(updated_germline_classification=group_labeller)) +
  scale_fill_manual(values=c("#0073C2","#EFC000")) +
  labs(x="",y="Count", fill="HR status") + ylim(0,150) +
  theme_Publication() +
  theme(strip.text = element_text(size=10), legend.position = "none")
```

```{r, germline_classification_hrdetect,  fig.width=9, fig.height=3.2}
result_df |>
  mutate(hr_status=ifelse(Probability>=0.7,"HRD","non-HRD")) |> 
  ggplot(aes(x=hr_status, fill=hr_status)) +
  geom_bar() + 
  geom_text(stat='count', aes(x=hr_status,label=..count..), vjust=-.2, size=4) +
  facet_grid(~updated_germline_classification, labeller = labeller(updated_germline_classification=group_labeller)) +
  scale_fill_manual(values=c("#0073C2","#EFC000")) +
  labs(x="",y="Count", fill="HR status") + ylim(0,150) +
  theme_Publication() +
  theme(strip.text = element_text(size=10), legend.position = "none")
```

```{r, germline_classification_hrdsum,  fig.width=9, fig.height=3.2}
result_df |>
  mutate(hr_status=ifelse(HRDsum>=42,"HRD","non-HRD")) |> 
  ggplot(aes(x=hr_status, fill=hr_status)) +
  geom_bar() + 
  geom_text(stat='count', aes(x=hr_status,label=..count..), vjust=-.2, size=4) +
  facet_grid(~updated_germline_classification, labeller = labeller(updated_germline_classification=group_labeller)) +
  scale_fill_manual(values=c("#0073C2","#EFC000")) +
  labs(x="",y="Count", fill="HR status") + ylim(0,150) +
  theme_Publication() +
  theme(strip.text = element_text(size=10), legend.position = "none")
```

## plot them together

```{r, germline_classification_all,  fig.width=9, fig.height=7.2}
result_df |>
  mutate(HRDetect=ifelse(Probability>=0.7,"HRD","non-HRD"), 
         HRDsum=ifelse(HRDsum>=42,"HRD","non-HRD"), 
         CHORD=ifelse(hr_status=="HR_deficient","HRD","non-HRD") ) |> 
  select(cohort, sample_id, updated_germline_classification, HRDetect, CHORD, HRDsum) |> 
  pivot_longer(c(HRDetect, CHORD, HRDsum)) |> 
  ggplot(aes(x=value, fill=value)) +
  geom_bar() + 
  geom_text(stat='count', aes(x=value,label=..count..), vjust=-.2, size=4) +
  facet_grid(name~updated_germline_classification, labeller = labeller(updated_germline_classification=group_labeller)) +
  scale_fill_manual(values=c("#0073C2","#EFC000")) +
  labs(x="",y="Count", fill="HR status") + ylim(0,150) +
  theme_Publication() +
  theme(strip.text = element_text(size=10), legend.position = "none")
```