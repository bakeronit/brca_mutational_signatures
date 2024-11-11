.is_na <- function(x) {
  purrr::is_scalar_vector(x) && is.na(x)
}

my_plot_contribution <- function (contribution, signatures = NA, index = NA, coord_flip = FALSE, 
          mode = c("relative", "absolute"), palette = NA) 
{
  mode <- match.arg(mode)
  if (!.is_na(index)) {
    contribution <- contribution[, index, drop = FALSE]
  }
  Sample <- Contribution <- Signature <- NULL
  if (mode == "absolute" & !.is_na(signatures)) {
    total_signatures <- colSums(signatures)
    abs_contribution <- contribution * total_signatures
  }
  #all_sig_level <- c( paste( colnames(cosmic2), "like", sep = "-"), rownames(contribution) |>  grep(pattern="-like", x=_, value=TRUE, invert=TRUE))
  #all_sig_level <- c( paste( colnames(cosmic2), "like", sep = "-"), c("SBSA","SBS3-like","SBSB","SBS18-like"))
  tb <- contribution %>% as.data.frame() %>% tibble::rownames_to_column("Signature") %>% 
    tidyr::pivot_longer(-Signature, names_to = "Sample", 
                        values_to = "Contribution") %>% dplyr::mutate(Sample = factor(Sample, 
                                                                                      levels = unique(Sample)), Signature = factor(Signature, 
                                                                                                                                   levels = unique(Signature))) |>  full_join(meta_data |> select(Donor, BRCA_status), by=c("Sample"="Donor"))
  if (mode == "absolute") {
    bar_geom <- geom_bar(stat = "identity")
    y_lab <- "Absolute contribution \n (no. mutations)"
  }
  else if (mode == "relative") {
    bar_geom <- geom_bar(position = "fill", stat = "identity"
                         )
    y_lab <- "Relative contribution"
  }
  present_sigs <- tb %>% dplyr::filter(Contribution != 0) %>% 
    dplyr::pull(Signature) %>% unique()
  plot <- ggplot(tb, aes(x = Sample, y = Contribution, fill = Signature, color= Signature)) + 
    bar_geom + labs(x = "", y = y_lab) + scale_fill_discrete(breaks = present_sigs) + 
    theme_pubclean(base_size = 12) + theme(panel.grid.minor.x = element_blank(), 
                       panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), 
                       panel.grid.major.y = element_blank())
  if (!.is_na(palette)) {
    plot <- plot + scale_fill_manual(name = "Signature", 
                                     values = palette) +
      scale_color_manual(name= "Signature", values=palette)
  }
  if (coord_flip) {
    #plot <- plot + coord_flip() + xlim(rev(levels(factor(tb$Sample)))) 
    plot <- plot + coord_flip() + facet_wrap(~BRCA_status, scales = "free_x")
  }
  else {
    #plot <- plot + xlim(levels(factor(tb$Sample))) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    plot <- plot + facet_wrap(~BRCA_status, scales = "free_x") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  return(plot)
}

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vetical",
            legend.key.size= unit(0.5, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}
