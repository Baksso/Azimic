library(tidyverse)
library(readxl)

setwd("/Users/bakarysanyang/Desktop/R_workspace/MICA1/")

#Baseline characteristics
baseline <- read_excel("/Users/bakarysanyang/Desktop/R_workspace/MICA1/baseline.xlsx", sheet = 1)


#test continuous variables for normality by Shapiro-wilk test

#maternal age

with(baseline, shapiro.test(maternal_age[arm == "azi"]))

with(baseline, shapiro.test(maternal_age[arm == "placebo"]))

#birth weight 

with(baseline, shapiro.test(birthwgh[arm == "azi"]))

with(baseline, shapiro.test(birthwgh[arm == "placebo"]))


#test for homogeneity of variance in maternal age and birth weight between arms

var.test(maternal_age ~ arm, data = baseline)

var.test(birthwgh ~ arm, data = baseline)


#nonparametric test for difference in maternal age and birth weight between treatment groups

ma_age_stat <- wilcox.test(maternal_age ~ arm, data = baseline, alternative = "two.sided")

birthwgt_stat <- wilcox.test(birthwgh ~ arm, data = baseline, alternative = "two.sided")


#distribution of categorical variables between treatment arms

#ethnicity
ethnicity_stat <- with(baseline, chisq.test(Ethnicity, arm))

#sex
sex_stat <- with(baseline, chisq.test(sex_at_delivery, arm))

#birth season
season_stat <- with(baseline, chisq.test(Season, arm))


#initial QC 

#DNA concentrations

biomass_info <- read.csv(file = "/Users/bakarysanyang/Desktop/R_workspace/MICA1/Biomass_info.csv")

biomass_info$time_point <- factor(biomass_info$time_point, 
                                  levels = c("day 0", "day 6","day 28", "12 months", "Blank"))

ggplot(biomass_info, aes(x=time_point, y=Initial_Conc, color=Arm)) +
  geom_point (position = position_jitterdodge
              (jitter.width = 0.2, jitter.height = 0, dodge.width = 0.8), fill='white', alpha=0.6) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha=0.6, outlier.colour = NA) +
  scale_color_manual(values = c("dark orange", "dark blue", "dark red"), name=NULL,
                     breaks = c("Azithromycin", "Placebo", "Control"),
                     labels = c("Azithromycin", "Placebo", "Control")) +
  theme_classic() +
  labs(y = "Initial DNA concentration (ng/ul)",
       x = "Time Point") +
  theme(axis.title = element_text(size = 10)) 


#identification and removal of contaminants

library(decontam)

otu_data <- read.delim("~/Desktop/R_workspace/MICA1/final.opti_mcc.0.03.filter.shared") %>%
  select(-label)

qc_metadata <- read.csv("~/Desktop/R_workspace/MICA1/qc_metadata.csv")

otu_matrix <- data.matrix(otu_data)

is.matrix(otu_matrix)

is.numeric(otu_matrix)


iscontam_freqprev <- isContaminant(otu_matrix, conc=qc_metadata$conc, neg=qc_metadata$samdf,
                                         method = "either", batch= qc_metadata$Extraction_batch, threshold = c(0.1, 0.5)) %>% filter(contaminant=="TRUE")

#write.csv(iscontam_freqprev, file = "~/Desktop/R_workspace/MICA1/iscontam_freqprev.csv")


#final community analysis

#read in data

metadata <- read_excel("metadata_final.xlsx") %>%
  rename(sample=group) %>%
  mutate(sample=as.character(sample))

taxonomy <- read_tsv(file="/Users/bakarysanyang/Desktop/R_workspace/MICA1/final.clean.cons.pick.taxonomy") %>%
  mutate(Taxonomy=str_replace_all(string = Taxonomy, pattern = "\\(\\d*\\)", replacement="")) %>%
  mutate(Taxonomy=str_replace_all(string = Taxonomy, pattern = ";$", replacement = "")) %>%
  separate(Taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";")


otu_data <- read_tsv("final.clean.subsample.shared", col_types = cols(Group=col_character())) %>%
  select(-label, -numOtus) %>%
  rename(sample=Group) %>%
  pivot_longer(cols = -sample, names_to = "otu", values_to = "count") %>%
  group_by(sample)%>%
  mutate(total_count = sum(count), r_abund=count/total_count)

otu_data %>% group_by(sample) %>% summarise(n=sum(r_abund)) %>% summary()

comb_genus_data <- inner_join(otu_data, taxonomy, by=c('otu' = 'OTU')) %>%
  group_by(sample, genus) %>%
  summarise(agg_r_abund=sum(r_abund)) %>%
  inner_join(., metadata)%>%
  ungroup()

#genus profile

top20_genus <- comb_genus_data %>%
  group_by(genus) %>%
  summarise(mean=mean(agg_r_abund)) %>%
  arrange((desc(mean))) %>% 
  top_n(20, mean) %>%
  pull(genus) 


mycolor <- c("#4C5270", "#EBEBE8","#BCECE0", "#3D550C", "#81B622", "#59981A",
             "#778A35", "#DAD870", "#ECF87F","#D1E2C4", "#31352E", "#8F3B21", 
             "#DD9927", "#FF6501", "#B78326", "#86785C",
             "#ADA7A7", "#FFB85D", "#FF5412", "#FFA384")


top20 <- comb_genus_data %>%
  filter(genus %in% top20_genus)

top20 %>% group_by(sample) %>% summarise(n=sum(agg_r_abund)) %>% summary()

top20$Time_Point <- factor(top20$Time_Point, 
                           levels = c("day 0", "day 6", "day 28", "12 months"))


genus_plot <- top20 %>%
  group_by(genus, Time_Point, Arm) %>%
  summarise(sum_agg_r_abund = sum(agg_r_abund)) %>%
  ggplot(aes(x=Arm, y=sum_agg_r_abund, fill=genus)) +
  geom_col(position = 'fill', width = 0.6) +
  scale_x_discrete (limits = c("Azi", "Placebo")) +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, ) +
  labs(x = " ", y= "Relative Abundance", tag = " ") +
  scale_fill_manual(values = mycolor, aesthetics = c( "color", "fill")) +
  facet_grid(facets = 'Time_Point')
genus_plot


#phylum profile

comb_phylum_data <- inner_join(otu_data, taxonomy, by=c('otu' = 'OTU')) %>%
  group_by(sample, phylum) %>%
  summarise(agg_r_abund=sum(r_abund)) %>%
  inner_join(., metadata)%>%
  ungroup()


mycolor2 <- c("#4C5270", "#EBEBE8","#FFB85D", "#FF5412", "#FFA384","#BCECE0", 
              "#3D550C", "#81B622","#778A35", "#DAD870", "#ECF87F","#D1E2C4")


comb_phylum_data$Time_Point <- factor(comb_phylum_data$Time_Point, 
                                      levels = c("day 0", "day 6", "day 28", "12 months"))

phylum_plot <- ggplot(comb_phylum_data, aes(x=Arm, y=agg_r_abund, fill=phylum)) +
  geom_col(position = 'fill', width = 0.6) +
  scale_x_discrete (limits = c("Azi", "Placebo")) +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, ) +
  labs(x = " ", y= "Relative Abundance") +
  scale_fill_manual(values = mycolor2, aesthetics = c( "color", "fill")) +
  facet_grid(facets = 'Time_Point')

phylum_plot


#alpha diversity analysis

library(readr)
library(lmerTest)
library(emmeans)

alpha <- read_tsv(file="/Users/bakarysanyang/Desktop/R_workspace/MICA1/alpha_diversity/final.clean.subsample.groups.summary") 

alpha$group <- as.character(alpha$group)

metadata_alpha <- inner_join(metadata, alpha, by=c('sample' = 'group'))

#compare alpha diversity between trial arms at each time-point by t-test

day0_alpha <- metadata_alpha %>% group_by(sample, Arm, Time_Point)%>% 
  filter(Time_Point=="day 0")

day0_ttest <- t.test(day0_alpha$npshannon ~ day0_alpha$Arm, alternative = "two.sided",
                     var.equal = TRUE, p.adjust.methods="hochberg")
day0_ttest

day6_alpha <- metadata_alpha %>% group_by(sample, Arm, Time_Point)%>% 
  filter(Time_Point=="day 6")

day6_ttest <- t.test(day6_alpha$npshannon ~ day6_alpha$Arm, alternative = "two.sided",
                     var.equal = TRUE, p.adjust.methods="hochberg")
day6_ttest

day28_alpha <- metadata_alpha %>% group_by(sample, Arm, Time_Point)%>% 
  filter(Time_Point=="day 28")

day28_ttest <- t.test(day28_alpha$npshannon ~ day28_alpha$Arm, alternative = "two.sided",
                      var.equal = TRUE, p.adjust.methods="hochberg")
day28_ttest

m12_alpha <- metadata_alpha %>% group_by(sample, Arm, Time_Point)%>% 
  filter(Time_Point=="12 months")

m12_ttest <- t.test(m12_alpha$npshannon ~ m12_alpha$Arm, alternative = "two.sided",
                    var.equal = TRUE, p.adjust.methods="hochberg")
m12_ttest


#alpha diversity box plot by time_point

ggplot(metadata_alpha, aes(x=Time_Point, y=npshannon, color=Arm)) +
  geom_point (position = 
                position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.8), fill='white', alpha=0.6) +
  scale_color_manual(name=NULL, values = c("dark orange", "dark blue"),
                     breaks = c("Azi", "Placebo"),
                     labels = c("Azithromycin", "Placebo")) +
  geom_boxplot(aes(x=Time_Point, y=npshannon, color=Arm),
               alpha=0.6, outlier.colour = NA, position = position_dodge(width = 0.8)) +
  
  scale_x_discrete(limits= c("day 0", "day 6", "day 28", "12 months")) +
  labs(x="Time points",
       y="Shannon index",
  ) +
  theme_classic() +
  theme(axis.title = element_text(size = 10)) 


#variance in alpha diversity by age

library(lmerTest)
library(emmeans)

metadata_alpha2 <- metadata_alpha %>% rename("maternal_age" = `maternal age`)%>%
  mutate(Time_Point = fct_relevel(Time_Point, "day 0", "day 6", "day 28", "12 months")) 


# Running models by Arm, testing changes in alpha diversity over time

## Azi 

azi_alpha <- metadata_alpha2 %>% filter(Arm=="Azi")

lm_alpha_azi <- lmerTest::lmer(npshannon ~ Time_Point + sexdelvi + maternal_age + (1|Rndo), data = azi_alpha)
lm_alpha_azi %>% broom.mixed::tidy()

time_d0vsall_azi <- emmeans(lm_alpha_azi, trt.vs.ctrlk ~ Time_Point, ref = 1, adjst="Dunnett")
time_d0vsall_azi

time_d6vsall_azi <- emmeans(lm_alpha_azi, trt.vs.ctrlk ~ Time_Point, ref = 2, adjst="Dunnett")
time_d6vsall_azi

time_d28vsall_azi <- emmeans(lm_alpha_azi, trt.vs.ctrlk ~ Time_Point, ref = 3, adjst="Dunnett")
time_d28vsall_azi

time_12mvsall_azi <- emmeans(lm_alpha_azi, trt.vs.ctrlk ~ Time_Point, ref = 4, adjst="Dunnett")
time_12mvsall_azi

##placebo

placebo_alpha <- metadata_alpha2 %>% filter(Arm=="Placebo")

lm_alpha_p <- lmerTest::lmer(npshannon ~ Time_Point + sexdelvi + maternal_age + (1|Rndo), data = placebo_alpha)
lm_alpha_p %>% broom.mixed::tidy()

time_d0vsall_p <- emmeans(lm_alpha_p, trt.vs.ctrlk ~ Time_Point, ref = 1, adjst="Dunnett")
time_d0vsall_p 

time_d6vsall_p <- emmeans(lm_alpha_azi, trt.vs.ctrlk ~ Time_Point, ref = 2, adjst="Dunnett")
time_d6vsall_p

time_d28vsall_p <- emmeans(lm_alpha_azi, trt.vs.ctrlk ~ Time_Point, ref = 3, adjst="Dunnett")
time_d28vsall_p

time_12mvsall_p <- emmeans(lm_alpha_azi, trt.vs.ctrlk ~ Time_Point, ref = 4, adjst="Dunnett")
time_12mvsall_p


#beta diversity analysis

NMDS <- read_tsv(file="/Users/bakarysanyang/Desktop/R_workspace/MICA1/final.clean.subsample.braycurtis.0.03.lt.nmds.axes")

NMDS$group <- as.character(NMDS$group)

# metadata <- read_excel(path="/Users/bakarysanyang/Desktop/R_workspace/MICA1/beta_diversity/metadata_final.xlsx")

metadata_NMDS <- inner_join(metadata, NMDS, by=c('sample' = 'group'))

metadata_NMDS_new <- metadata_NMDS

metadata_NMDS_new$Time_Point <-factor(metadata_NMDS_new$Time_Point,
                                      levels = c("day 0", "day 6", "day 28", "12 months"))

#nmds ordination

nmd_plot <- ggplot(metadata_NMDS_new, aes(x=axis1, y=axis3, color=Time_Point)) +
  geom_point(size=1, alpha=0.6) +
  scale_color_brewer(name=NULL, type = "qualitative", palette = 'Set1',
                     breaks = c("day 0", "day 6", "day 28", "12 months"),
                     labels = c("day 0", "day 6", 
                                "day 28", "12 months")) +
  stat_ellipse(aes(x=axis1, y=axis3, linetype=Arm)) +
  coord_fixed() +
  labs(x="NMDS Axis 1",
       y="NMDS Axis 3") +
  theme_classic()

nmd_plot


#beta diversity stats

library(phyloseq)
library(vegan)

bray_dist <- import_mothur_dist("final.clean.subsample.braycurtis.0.03.lt.dist")
df_bray_dist <- bray_dist%>%as.matrix()%>% as.data.frame()

#Azi vs placebo per time-point

#day0
metadata_d0 <- metadata%>% filter(Time_Point=="day 0")

day0_bray <- import_mothur_dist("final.day0.subsample.braycurtis.0.03.lt.dist")

betadisp_d0 <- betadisper(day0_bray, group = metadata_d0$Arm)
betadisp_d0_perm <- permutest(betadisp_d0, pairwise = TRUE, 
                              control = how(blocks = metadata_d0$Arm, within(type="free")))
betadisp_d0_perm


df_d0_bray <- day0_bray%>%as.matrix()%>% as.data.frame()

betadiv_d0 <- adonis2(df_d0_bray ~ Arm + Ethnicity, data = metadata_d0, 
                      permutations = with(metadata_d0, how(nperm=999),
                                          method = "bray"))
betadiv_d0


#day6
metadata_d6 <- metadata%>% filter(Time_Point=="day 6")

day6_bray <- import_mothur_dist("final.day6.subsample.braycurtis.0.03.lt.dist")

betadisp_d6 <- betadisper(day6_bray, group = metadata_d6$Arm)
betadisp_d6_perm <- permutest(betadisp_d6, pairwise = TRUE, 
                              control = how(blocks = metadata_d6$Arm, within(type="free")))
betadisp_d6_perm

df_d6_bray <- day6_bray%>%as.matrix()%>% as.data.frame()

betadiv_d6 <- adonis2(df_d6_bray ~ Arm + Ethnicity, data = metadata_d6, 
                      permutations = with(metadata_d6, how(nperm=999),
                                          method = "bray"))
betadiv_d6


#day28

metadata_d28 <- metadata%>% filter(Time_Point=="day 28")

day28_bray <- import_mothur_dist("final.day28.subsample.braycurtis.0.03.lt.dist")

betadisp_d28 <- betadisper(day28_bray, group = metadata_d28$Arm)
betadisp_d28_perm <- permutest(betadisp_d28, pairwise = TRUE, 
                               control = how(blocks = metadata_d28$Arm, within(type="free")))
betadisp_d28_perm

df_d28_bray <- day28_bray%>%as.matrix()%>% as.data.frame()

betadiv_d28 <- adonis2(df_d28_bray ~ Arm + Ethnicity, data = metadata_d28, 
                       permutations = with(metadata_d28, how(nperm=999, blocks =Arm),
                                           method = "bray", by="margin"))
betadiv_d28


#12months

metadata_m12 <- metadata%>% filter(Time_Point=="12 months")

oneyr_bray <- import_mothur_dist("final.oneyr.subsample.braycurtis.0.03.lt.dist")

betadisp_1yr <- betadisper(oneyr_bray, group = metadata_m12$Arm)
betadisp_1yr_perm <- permutest(betadisp_1yr, pairwise = TRUE, 
                               control = how(blocks = metadata_m12$Arm, within(type="free")))
betadisp_1yr_perm

df_1yr_bray <- oneyr_bray%>%as.matrix()%>% as.data.frame()

betadiv_1yr <- adonis2(df_1yr_bray ~ Arm + Ethnicity, data = metadata_m12, 
                       permutations = with(metadata_m12, how(nperm=999, blocks = Arm),
                                           method = "bray", by="margin"))
betadiv_1yr

#variance in beta diversity within trial arm

#Azi

metadata_azi <- metadata%>%filter(Arm=="Azi")

azi_dist <- import_mothur_dist("final.Azi.subsample.braycurtis.0.03.lt.dist")

betadisp_azi <- betadisper(azi_dist, group = metadata_azi$Time_Point)
betadisp_azi_perm <- permutest(betadisp_azi, pairwise = TRUE,
                                          control = how(blocks = metadata_azi$Rndo, within(type="free")))
betadisp_azi_perm

df_azi_dist <- as.matrix(azi_dist)%>% data.frame() 

azi_betadiv <- adonis2(df_azi_dist ~  Time_Point, data = metadata_azi,
                       permutations = with(metadata_azi, how(nperm = 999, blocks = Rndo) 
                       ),
                       method = "bray", by="margin")
azi_betadiv


#placebo

# metadata_placebo <- read_tsv(file = "metadata_Placebo.txt")

metadata_placebo <- metadata%>% filter(Arm=="Placebo")

placebo_dist <- import_mothur_dist("final.Placebo.subsample.braycurtis.0.03.lt.dist")

betadisp_p <- betadisper(placebo_dist, group = metadata_placebo$Time_Point)
betadisp_p_perm <- permutest(betadisp_p, pairwise = TRUE,
                             control = how(blocks = metadata_placebo$Rndo, within(type="free")))
betadisp_p_perm

df_placebo_dist <- as.matrix(placebo_dist)%>%as.data.frame()

placebo_betadiv <- adonis2(df_placebo_dist ~  Time_Point, data = metadata_placebo,
                           permutations = with(metadata_placebo, how(nperm = 999, blocks = Rndo)), 
                           method = "bray", by="margin")
placebo_betadiv


#community types

#partition distribution - stacked bar chart

library(colorspace)

partitions <- read_excel(path = "/Users/bakarysanyang/Desktop/R_workspace/MICA1/partitions.xlsx")

stacked_bar <- ggplot(partitions, aes(x=Time_Point, fill=Partition)) +
  geom_bar(position = 'fill') +
  scale_x_discrete(limits=c("d0 A", "d0 P", "d6 A", "d6 P",
                            "d28 A", "d28 P", "12m A", "12m P")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x="Time Points",
       y="Percentage of samples", tag = "a") +
  scale_fill_discrete_sequential(name=NULL, palette = "ag_Sunset") +
  theme_classic() 

stacked_bar

#partitions distribution between arms 

day0_partitions <- partitions %>% filter(Time_Point=="d0 A" | Time_Point=="d0 P")

day6_partitions <- partitions %>% filter(Time_Point=="d6 A" | Time_Point=="d6 P")

day28_partitions <- partitions %>% filter(Time_Point=="d28 A" | Time_Point=="d28 P")

m12_partitions <- partitions %>% filter(Time_Point=="12m A" | Time_Point=="12m P")


#create contingency tables and compute fisher's exact test

#day 0
con_tab_day0 <- table(day0_partitions$Time_Point, day0_partitions$Partition)

fisher_day0 <- fisher.test(con_tab_day0, alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)

fisher_day0

#day 6 
con_tab_day6 <- table(day6_partitions$Time_Point, day6_partitions$Partition)

fisher_day6 <- fisher.test(con_tab_day6, alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)

fisher_day6

#day 28
con_tab_day28 <- table(day28_partitions$Time_Point, day28_partitions$Partition)

fisher_day28 <- fisher.test(con_tab_day28, alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)

fisher_day28

#12months
con_tab_m12 <- table(m12_partitions$Time_Point, m12_partitions$Partition)

fisher_m12 <- fisher.test(con_tab_m12, alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)

fisher_m12


#partitions heatmap

library(microbiomer)
library(ggtext)

theme_set(theme_light())


##load files 
ps <- import_biom("/Users/bakarysanyang/Desktop/R_workspace/MICA1/final.clean.subsample.partition.biom")


# load functions

hc_order <- function(ps) {
  ps.m <- ps %>%
    otu_tab_to_df() %>%
    column_to_rownames("sample_id") %>%
    as.matrix() 
  
  bc <- vegan::vegdist(ps.m, method = "bray")
  bc.m <-as.matrix(bc)
  hc <- hclust(bc, method="average") 
  
  ps_reorder_samples <- function(ps, order_by)  { # order_by is numeric
    sample_data(ps) <- sample_data(ps)[order_by, ]
    otu_table(ps) <- otu_table(ps)[, order_by]
    return(ps)
  }
  
  ps_ord <- ps_reorder_samples(ps, hc$order)
  sample_data(ps_ord)$hc_sample_id <- fct_inorder(rownames(sample_data(ps_ord)))
  return(ps_ord)
}

#heatmap

#data(ps_ord)
ps_ord <- hc_order(ps)

ps_ord <- ps_ord %>% to_RA() # convert to relative abundance

ps_ord %>% meta_to_df() %>% head

#create dummy partition variable

set.seed(187)

## hierarchical clustering

bc <- ps_ord %>%
  otu_table() %>% 
  t %>%
  vegdist(., "bray")

hc <- hclust(bc, "average")

# prepare heatmap

otu_hm <- as(otu_table(ps_ord), "matrix")

otu_hm_order <- otu_hm[order(rowMeans(otu_hm), decreasing = T)[1:20], hc$order] # change if you want to show more or less OTUs
meta_hm_order <- ps_ord %>% meta_to_df() %>% .[hc$order,]

all(colnames(otu_hm_order) == hc$labels[hc$order]) # should be TRUE

prep_hm_data <- function(otu_hm_order, meta_hm_order) {
  hm_data <- otu_hm_order %>%
    data.frame(check.names = F) %>%
    rownames_to_column("OTU") %>%
    pivot_longer(-OTU, names_to = "sample_id", values_to = "RA") %>%
    left_join(
      .,
      meta_hm_order, by = "sample_id") 
  
  return(hm_data)
}

hm_data <- prep_hm_data(otu_hm_order, meta_hm_order)

#annotate OTUs with taxonomy

hm_data_new <- inner_join(hm_data, taxonomy, by="OTU")
hm_data_new$ann_OTU <- base::paste(hm_data_new$OTU, hm_data_new$genus, sep = "_")

#convert variable class

hm_data_new$ann_OTU <- as.factor(hm_data_new$ann_OTU)
hm_data_new$sample_id <- as.factor(hm_data_new$sample_id)

#order variables

hm_data_new$ann_OTU <- fct_inorder(hm_data_new$ann_OTU) %>%
  fct_rev() 
hm_data_new$sample_id <- fct_inorder(hm_data_new$sample_id)

#draw heat map
hm2 <- hm_data_new %>%
  ggplot(aes(x = sample_id, y = ann_OTU, fill = RA)) +
  geom_tile() +
  scale_fill_gradientn(name="Relative abundance", colors = RColorBrewer::brewer.pal(9, "OrRd"), na.value = "grey90") +
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.text.y = element_markdown(size = 8), axis.ticks.x = element_blank(), legend.position = "bottom") +
  xlab("Sample") + ylab("OTU") + labs(tag = "b")

#plot heat map

faceted_heatmap <- hm2 + 
  facet_grid(~ Partition, space = "free", scales = "free_x")


#combine partitions stacked bar and heatmap

library(patchwork)
stacked_bar/faceted_heatmap 


#differential abundance analysis

library(metagenomeSeq)

taxa <- taxonomy %>% select(-Size)%>% column_to_rownames(var = "OTU")%>%as.matrix()
is.matrix(taxa) #should be true

metadata_ps <- metadata %>% column_to_rownames(var = "sample")

otu_shared <- import_mothur(mothur_shared_file = "final.clean.subsample.shared")

row.names(otu_shared)
is.matrix(otu_shared) #should be true
is.numeric(otu_shared) #should be true

ps_obj <- phyloseq(otu_table(otu_shared, taxa_are_rows = TRUE), tax_table(taxa), sample_data(metadata_ps))

taxa_label <- taxonomy %>% select(-Size)

#phylum level
ps_phylum <- tax_glom(ps_obj, taxrank = "phylum")

mrexp_phy <- phyloseq_to_metagenomeSeq(ps_phylum)

#day6 

d6_samp_phy <- which(pData(mrexp_phy)$Time_Point == "day 6")

mrexp_phy_day6 <- mrexp_phy[g_features_to_keep, d6_samp_phy]

mrexp_phy_day6 = filterData(mrexp_phy_day6, depth = 2)

mrexp_phy_day6 <- cumNorm(mrexp_phy_day6, p=0.75)

day6_meta_p <- pData(mrexp_phy_day6) 

day6_meta_p$Arm <- factor(day6_meta_p$Arm, levels = c("Placebo", "Azi"))

day6_model_p <- model.matrix(~ 1 + Arm, data = day6_meta_p)

day6_fitft_p = fitFeatureModel(obj = mrexp_phy_day6, mod = day6_model_p)


head(MRcoefs(day6_fitft_p))

logFC_day6_p <- MRcoefs(day6_fitft_p, number = nrow(mrexp_phy_day6)) %>% rownames_to_column()%>%
  plyr::rename(c("rowname"="OTU"))%>% inner_join(.,taxa_label, by ="OTU")

sig_feat_day6_p <- logFC_day6_p %>% filter(adjPvalues < 0.1) %>% arrange(adjPvalues)

sig_feat_day6_p


#day28

d28_samp_phy <- which(pData(mrexp_phy)$Time_Point == "day 28")

mrexp_phy_day28 <- mrexp_phy[g_features_to_keep, d28_samp_phy]

mrexp_phy_day28 = filterData(mrexp_phy_day28, depth = 2)

mrexp_phy_day28 <- cumNorm(mrexp_phy_day28, p=0.75)

day28_meta_p <- pData(mrexp_phy_day28)

day28_meta_p$Arm <- factor(day28_meta_p$Arm, levels = c("Placebo", "Azi"))

day28_model_p <- model.matrix(~ 1 + Arm, data = day28_meta_p)

day28_fitft_p = fitFeatureModel(obj = mrexp_phy_day28, mod = day28_model_p)


head(MRcoefs(day28_fitft_p))

logFC_day28_p <- MRcoefs(day28_fitft_p, number = nrow(mrexp_phy_day28)) %>% rownames_to_column()%>%
  plyr::rename(c("rowname"="OTU"))%>% inner_join(.,taxa_label, by ="OTU")

sig_feat_day28_p <- logFC_day28_p %>% filter(adjPvalues < 0.1) %>% arrange(adjPvalues)

sig_feat_day28_p


#12months

m12_samp_phy <- which(pData(mrexp_phy)$Time_Point == "12 months")

mrexp_phy_m12 <- mrexp_phy[g_features_to_keep, m12_samp_phy]

mrexp_phy_m12 = filterData(mrexp_phy_m12, depth = 2)

mrexp_phy_m12 <- cumNorm(mrexp_phy_m12, p=0.75)

m12_meta_p <- pData(mrexp_phy_m12)

m12_meta_p$Arm <- factor(m12_meta_p$Arm, levels = c("Placebo", "Azi"))

m12_model_p <- model.matrix(~ 1 + Arm, data = m12_meta_p)

m12_fitft_p = fitFeatureModel(obj = mrexp_phy_m12, mod = m12_model_p)


head(MRcoefs(m12_fitft_p))

logFC_m12_p <- MRcoefs(m12_fitft_p, number = nrow(mrexp_phy_m12)) %>% rownames_to_column()%>%
  plyr::rename(c("rowname"="OTU"))%>% inner_join(.,taxa_label, by ="OTU")

sig_feat_m12_p <- logFC_m12_p %>% filter(adjPvalues < 0.1) %>% arrange(adjPvalues)

sig_feat_m12_p


#genus level
ps_genus <- tax_glom(ps_obj, taxrank = "genus")

mrexp_genus <- phyloseq_to_metagenomeSeq(ps_genus)

#day6 

d6_samp_genus <- which(pData(mrexp_genus)$Time_Point == "day 6")

mrexp_genus_day6 <- mrexp_genus[g_features_to_keep, d6_samp_genus]

mrexp_genus_day6 = filterData(mrexp_genus_day6, depth = 2)

mrexp_genus_day6 <- cumNorm(mrexp_genus_day6, p=0.75)

day6_meta_g <- pData(mrexp_genus_day6) 

day6_meta_g$Arm <- factor(day6_meta_g$Arm, levels = c("Placebo", "Azi"))

day6_model_g <- model.matrix(~ 1 + Arm, data = day6_meta_g)

day6_fitft_g = fitFeatureModel(obj = mrexp_genus_day6, mod = day6_model_g)


head(MRcoefs(day6_fitft_g))

logFC_day6_g <- MRcoefs(day6_fitft_g, number = nrow(mrexp_genus_day6)) %>% rownames_to_column()%>%
  plyr::rename(c("rowname"="OTU"))%>% inner_join(.,taxa_label, by ="OTU")

sig_feat_day6_g <- logFC_day6_g %>% filter(adjPvalues < 0.1) %>% arrange(adjPvalues)

sig_feat_day6_g


#day28

d28_samp_genus <- which(pData(mrexp_genus)$Time_Point == "day 28")

mrexp_genus_day28 <- mrexp_genus[g_features_to_keep, d28_samp_genus]

mrexp_genus_day28 = filterData(mrexp_genus_day28, depth = 2)

mrexp_genus_day28 <- cumNorm(mrexp_genus_day28, p=0.75)

day28_meta_g <- pData(mrexp_genus_day28)

day28_meta_g$Arm <- factor(day28_meta_g$Arm, levels = c("Placebo", "Azi"))

day28_model_g <- model.matrix(~ 1 + Arm, data = day28_meta_g)

day28_fitft_g = fitFeatureModel(obj = mrexp_genus_day28, mod = day28_model_g)


head(MRcoefs(day28_fitft_g))

logFC_day28_g <- MRcoefs(day28_fitft_g, number = nrow(mrexp_genus_day28)) %>% rownames_to_column()%>%
  plyr::rename(c("rowname"="OTU"))%>% inner_join(.,taxa_label, by ="OTU")

sig_feat_day28_g <- logFC_day28_g %>% filter(adjPvalues < 0.1) %>% arrange(adjPvalues)

sig_feat_day28_g


#12months

m12_samp_genus <- which(pData(mrexp_genus)$Time_Point == "12 months")

mrexp_genus_m12 <- mrexp_genus[g_features_to_keep, m12_samp_genus]

mrexp_genus_m12 = filterData(mrexp_genus_m12, depth = 2)

mrexp_genus_m12 <- cumNorm(mrexp_genus_m12, p=0.75)

m12_meta_g <- pData(mrexp_genus_m12)

m12_meta_g$Arm <- factor(m12_meta_g$Arm, levels = c("Placebo", "Azi"))

m12_model_g <- model.matrix(~ 1 + Arm, data = m12_meta_g)

m12_fitft_g = fitFeatureModel(obj = mrexp_genus_m12, mod = m12_model_g)


head(MRcoefs(m12_fitft_g))

logFC_m12_g <- MRcoefs(m12_fitft_g, number = nrow(mrexp_genus_m12)) %>% rownames_to_column()%>%
  plyr::rename(c("rowname"="OTU"))%>% inner_join(.,taxa_label, by ="OTU")

sig_feat_m12_g <- logFC_m12_g %>% filter(adjPvalues < 0.1) %>% arrange(adjPvalues)

sig_feat_m12_g


#differential abundance bar plot

diff_abund <- read_excel("/Users/bakarysanyang/Desktop/R_workspace/MICA1/differential_abund.xlsx") %>%
  filter(., Fold_change != 'NA') 

diff_ordered <- diff_abund

diff_ordered$Time_point <- factor(diff_abund$Time_point, levels = c("day 6", "day 28", "12 months"))

color <- ifelse(diff_ordered$Fold_change < 0, "dark blue", "orange")

diff_plot <- ggplot(diff_ordered, aes(x=reorder(Genus, Fold_change), y=Fold_change)) +
  geom_bar(stat = "identity", fill=color, show.legend = FALSE) +
  geom_hline(yintercept = 0, color = 1) +
  coord_flip() +
  scale_y_continuous() +
  labs(x="Genus", y= "Fold Change") +
  theme(axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4)) +
  theme_classic() +
  facet_wrap(facets =  'Time_point',  scales = "free_x" )

diff_plot

