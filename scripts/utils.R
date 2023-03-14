
## setting for figures

color_fig <- list(
  BRCA_status=c("lightgrey","darkgrey","black"),
  Gender=c("grey","blue3"),
  gender=c("grey","blue3"),
  `Age Group`=rev(c("black","#4C4C4C","#787878","#ACACAC","#D8D3D3")),
  age_group=rev(c("black","#4C4C4C","#787878","#ACACAC","#D8D3D3")),
  Grade=rev(c("lightgrey","blue3","darkgreen","purple3")),
  `ER/PR/HER2`=rev(c("lightgrey","grey","blue3")),
  Morphology=c("lightgrey","purple3","blue3","darkgreen","red3","orange3","yellow4"),
  icd_o_3=c("lightgrey","purple3","blue3","darkgreen","orange3")
)


fact_levels <- list(
  BRCA_status=c("BRCA1","BRCA2","non-BRCA"),
  `Age Group`=c("<40","40-50","50-60","60-70",">70"),
  Gender=c("Female","Male"),
  Grade=c("Grade I","Grade II","Grade III","n/a"),
  Morphology=c("IC NST","MDL","DCIS","ILC","Medullary","Metaplastic","Mucinous"),
  ER=c("Positive","Negative","n/a"),
  PR=c("Positive","Negative","n/a"),
  HER2=c("Positive","Negative","n/a"),
  age_group=c("<40","40-50","50-60","60-70",">70"),
  gender=c("FEMALE","MALE"),
  icd_o_3=c("8500/3","8510/3","8520/3","8523/3","8575/3")
)

legend_names <- list(
  BRCA_status="BRCA status",
  `Age Group`="Age group",
  Gender="Gender",
  gender="Gender",
  Grade="Grade",
  Morphology="Morphology",
  age_group="Age group",
  icd_o_3="Morphology",
  ER="ER",
  PR="PR",
  HER2="HER2"
)

## funcs

generate_row <- function(row_name, df) {
  geom_tile(data=df %>% filter(field==row_name) %>% droplevels, 
            mapping=aes(x=Donor, y=field, fill=factor(level, levels=fact_levels[[row_name]])), 
            color="white")
}
