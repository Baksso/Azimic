library(tidyverse)
library(readxl)
library(patchwork)


here <- getwd()
setwd(here)

#distribution of RVS microbiota type

partitions_mother <- read_excel(path = "metadata_mother.xlsx")

RVS_partition <- partitions_mother %>% 
  select(-group, -Sample_type, -Time_point, -Extraction_batch)%>%
  table(.)
  
partition_distr <- chisq.test(RVS_partition)


partition1_d_plot <- as.data.frame(RVS_partition)%>% 
  filter(partition=='Partition_1')%>%
  ggplot(aes(x=partition, y=Freq, fill=Treatment)) +
  geom_col(position = position_dodge2(), width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("dark orange", "dark blue"), 
                     breaks = c("Azi", "Placebo"),
                     labels = c("Azithromycin", "Placebo"))+
          scale_x_discrete (limits = c("Partition_1")) +
             theme_classic() +
             scale_y_continuous(breaks = seq(0, 40, by=5),
                                sec.axis = sec_axis(trans = ~./58, 
                                                    labels = scales::percent,
                                                    name = "Percentage"))+
             labs(x= " ", y= "Frequency")
  
partition1_d_plot         

partition2_d_plot <- as.data.frame(RVS_partition)%>%
  filter(partition=='Partition_2')%>%
    ggplot(aes(x=partition, y=Freq, fill=Treatment)) +
    geom_col(position = position_dodge2(), width = 0.6) +
  scale_fill_manual(values = c("dark orange", "dark blue"), 
                    breaks = c("Azi", "Placebo"),
                    labels = c("Azithromycin", "Placebo"))+
    scale_x_discrete (limits = c("Partition_2")) +
    theme_classic() +
    scale_y_continuous(breaks = seq(0, 40, by=5),
                       sec.axis = sec_axis(trans = ~./38, 
                                           labels = scales::percent,
                                           name = "Percentage"))+
                         labs(x= " ", y= "Frequency")
                         
partition1_d_plot + partition2_d_plot


#genus profile

metadata_m2 <- read_excel("metadata_mother.xlsx") %>%
  rename(sample=group) %>%
  mutate(sample=as.character(sample))


taxonomy_m <- read_tsv(file="mother.cons.subsample.taxonomy") %>%
  mutate(
    Taxonomy=str_replace_all(string = Taxonomy, pattern = "\\(\\d*\\)", replacement="")) %>%
  mutate(
    Taxonomy=str_replace_all(string = Taxonomy, pattern = ";$", replacement = "")) %>%
  separate(Taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";")

rel_abund <- function(x){
  x/sum(x)
}

otu_data_m <- read_tsv("mother.filter.0.03.subsample.shared", col_types = cols(Group=col_character())) %>%
  select(-label, -numOtus) %>%
  rename(sample=Group) %>%
  pivot_longer(cols = -sample, names_to = "otu", values_to = "count") %>%
  group_by(sample)%>%
  mutate(r_abund=rel_abund(count))

otu_data_m %>% group_by(sample) %>% summarise(n=sum(r_abund)) %>% summary()

comb_genus_data_m <- inner_join(otu_data_m, taxonomy_m, by=c('otu' = 'OTU')) %>%
  group_by(sample, genus) %>%
  summarise(agg_r_abund=sum(r_abund)) %>%
  inner_join(., metadata_m2)%>%
  ungroup()


top20_genus_m <- comb_genus_data_m %>%
  group_by(genus) %>%
  summarise(mean=mean(agg_r_abund)) %>%
  arrange((desc(mean))) %>% 
  top_n(20, mean) %>%
  pull(genus)

#make color palette
mycolor <- c("#4C5270", "#EBEBE8","#BCECE0", "#3D550C", "#81B622", "#59981A",
             "#778A35", "#DAD870", "#ECF87F","#D1E2C4", "#31352E", "#8F3B21", 
             "#DD9927", "#FF6501", "#B78326", "#86785C",
             "#ADA7A7", "#FFB85D", "#FF5412", "#FFA384")

#draw plot

top20_m <- comb_genus_data_m %>%
  filter(genus %in% top20_genus_m)

top20_m %>% group_by(sample) %>% summarise(n=sum(agg_r_abund)) %>% summary()

genus_plot_m <- top20_m %>%
  group_by(partition, genus) %>%
  summarise(sum_agg_r_abund = sum(agg_r_abund)) %>%
  ggplot(aes(x=partition, y=sum_agg_r_abund, fill=genus)) +
  geom_col(position = 'fill', width = 0.6, show.legend = FALSE) +
  scale_x_discrete (limits = c("Partition_1", "Partition_2")) +
  theme_classic() +
  scale_y_continuous(labels = scales::percent, ) +
  labs(x= " ", y= "Relative Abundance", tag = "a") +
  scale_fill_manual(values = mycolor, aesthetics = c( "color", "fill")) 


genus_plot_m

genus_plot2_m <- top20_m %>%
  group_by(partition, sample, genus) %>%
  summarise(sum_agg_r_abund = sum(agg_r_abund)) %>%
  ggplot(aes(x=sample, y=sum_agg_r_abund, fill=genus)) +
  geom_col(position = 'fill', width = 0.6) +
  scale_x_discrete(label=NULL) +
  theme_classic() +
  theme(legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, 'cm'))+
  scale_y_continuous(labels = scales::percent, ) +
  labs(x= " ", y= "Relative Abundance", tag = "b") +
  scale_fill_manual(values = mycolor, aesthetics = c( "color", "fill")) +
  facet_wrap(facets = 'partition', scales = "free", nrow = 2, 
             strip.position = "bottom")

genus_plot2_m

genus_plot_m + genus_plot2_m


#alpha-diversity

#shannon plot

alpha_div_m <- read_tsv(file = "mother.filter.0.03.subsample.groups.summary")%>%
  select(-label, -nseqs) %>% inner_join(., partitions_mother, 
                                        by= c('group'='group') )

shan_plot <- ggplot(alpha_div_m, aes(x=partition, y=npshannon, color=partition)) +
  geom_point (position = 
                position_jitterdodge(
                  jitter.width = 0.2, jitter.height = 0, dodge.width = 0.8), 
              fill = 'white', alpha=0.6, show.legend = FALSE) +
  scale_color_manual(values = c("red", "dark blue")) +
  geom_boxplot(aes(x=partition, y=npshannon, color=partition),
               alpha=0.6, outlier.colour = NA, 
               position = position_dodge(width = 0.8), show.legend = FALSE) +
  labs(x="Partitions",
       y="Shannon index", tag = "a") +
  theme_classic() +
  theme(axis.title = element_text(size = 10)) 

shan_plot

#richness plot

rich_plot <- ggplot(alpha_div_m, aes(x=partition, y=sobs, color=partition)) +
  geom_point (position = 
                position_jitterdodge(
                  jitter.width = 0.2, jitter.height = 0, dodge.width = 0.8), 
              fill = 'white', alpha=0.6, show.legend = FALSE) +
  scale_color_manual(values = c("red", "dark blue")) +
  geom_boxplot(aes(x=partition, y=sobs, color=partition),
               alpha=0.6, outlier.colour = NA,
               position = position_dodge(width = 0.8), show.legend = FALSE) +
  labs(x="Partitions",
       y="OTUs observed", tag = "b") +
  theme_classic() +
  theme(axis.title = element_text(size = 10)) 

shan_plot + rich_plot

#calculate difference in alpha diversity

shapiro.test(alpha_div_m$npshannon) #p-value <0.001

wilcox_alpha_m <- wilcox.test(npshannon ~ partition, data=alpha_div_m, alternative = "two.sided")
wilcox_alpha_m



# community composition ---------------------------------------------------

#ordination of community composition

nmds <- read_tsv(file = "mother.filter.0.03.subsample.braycurtis.0.03.lt.nmds.axes")

metadata_nmds <- inner_join(nmds, partitions_mother, by= c('group'='group'))

nmds_plot <- ggplot(metadata_nmds, aes(x=axis1, y=axis2, color=partition)) +
  geom_point(size=1, alpha=0.6) +
  scale_color_manual(values = c("red", "dark blue")) +
  stat_ellipse(aes(x=axis1, y=axis2), show.legend = FALSE) +
  coord_fixed() +
  labs(x="NMDS Axis 1",
       y="NMDS Axis 2") +
  theme_classic()

nmds_plot


library(phyloseq)
library(vegan)

mother_dist <- 
  import_mothur_dist("mother.filter.0.03.subsample.braycurtis.0.03.lt.dist")

#homogeneity of group dispersion

set.seed(1821965)

dispersion_m <- betadisper(mother_dist, group = metadata_m2$partition)
perm_dis_m <- permutest(dispersion_m)
perm_dis_m

disp_plot_m <- plot(dispersion_m, hull = FALSE, main = " ", sub = " ",
                    col = c("red", "dark blue"), 
                    ellipse = TRUE, label = FALSE,) + 
  legend(x=0.4, y=-0.3,legend = c("Partition_1", "Partition_2", "p = 0.062"), 
         text.col = c("red","dark blue", "black"))


#variance in community composition

mother_dist <- as.matrix(mother_dist) %>% as.data.frame()

permanova_partition_m <- adonis2(mother_dist ~ metadata_m2$partition, 
                                 data = mother_dist, permutations = 999, method = "bray")

permanova_partition_m



# prevalence and abundance of pathogenic taxa in microbiota types ---------


key_pathogens <- data.frame(
  num = c(1:5),
  genus = c("Staphylococcus", "Streptococcus", "Enterococcaceae_unclassified", 
            "Enterobacteriaceae_unclassified", "Pseudomonas"))

patho_genus_data_m <- inner_join(otu_data_m, taxonomy_m, by=c('otu' = 'OTU')) %>%
  group_by(sample, genus) %>%
  summarise(agg_r_abund=sum(r_abund)) %>%
  inner_join(., metadata_m2)%>%
  filter(genus %in% key_pathogens$genus) %>%
  ungroup()


#prevalence in partition_1

patho_genus_prev_m <- patho_genus_data_m %>% group_by(sample) %>% 
  mutate(Present = if_else (agg_r_abund >= 0.001, "yes", "no"))

patho_genus_partition1 <- patho_genus_prev_m %>% filter(partition=="Partition_1")

nsamples_p1 <- length(unique(patho_genus_partition1$sample))

partition1_prev <- table(patho_genus_partition1$genus, 
                         patho_genus_partition1$Present) %>% as.data.frame()%>%
  mutate(prevalence = Freq/nsamples_p1) %>% rename(genus=Var1, Present=Var2)


prev_partition1 <- ggplot(partition1_prev, 
                           aes(x=genus, y=prevalence, fill=Present))+
  geom_col(position = "dodge", width = 0.6) +
  scale_x_discrete () +
  theme_classic() +
  theme(text = element_text(size = 10))+
  coord_flip() +
  scale_y_continuous(labels = scales::percent, ) +
  labs(x= "Genus", y= "Prevalence", tag = "Partition_1") +
  scale_fill_manual(values = c("grey", "black"), aesthetics = c( "color", "fill"))
prev_partition1


#prevalence in partition_1

patho_genus_partition2 <- patho_genus_prev_m %>% filter(partition=="Partition_2")

nsamples_p2 <- length(unique(patho_genus_partition2$sample))

partition2_prev <- table(patho_genus_partition2$genus, 
                          patho_genus_partition2$Present) %>% as.data.frame()%>%
  mutate(prevalence = Freq/nsamples_p2) %>% rename(genus=Var1, Present=Var2)

prev_partition2 <- ggplot(partition2_prev, 
                          aes(x=genus, y=prevalence, fill=Present))+
  geom_col(position = "dodge", width = 0.6) +
  scale_x_discrete () +
  theme_classic() +
  theme(text = element_text(size = 10))+
  coord_flip() +
  scale_y_continuous(labels = scales::percent, ) +
  labs(x= "Genus", y= "Prevalence", tag = "Partition_2") +
  scale_fill_manual(values = c("grey", "black"), aesthetics = c( "color", "fill"))
prev_partition2

patho_prev_partition <- prev_partition1 / prev_partition2 
patho_prev_partition


