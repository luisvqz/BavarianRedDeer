##### rarefied


##  
##   888888ba                 dP    888888ba                             
##   88    `8b                88    88    `8b                            
##  a88aaaa8P' .d8888b. .d888b88    88     88 .d8888b. .d8888b. 88d888b. 
##   88   `8b. 88ooood8 88'  `88    88     88 88ooood8 88ooood8 88'  `88 
##   88     88 88.  ... 88.  .88    88    .8P 88.  ... 88.  ... 88       
##   dP     dP `88888P' `88888P8    8888888P  `88888P' `88888P' dP       
##  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
##                                                                       
##  d8888b.  a8888a  d8888b. d88  
##      `88 d8' ..8b     `88  88  
##  .aaadP' 88 .P 88 .aaadP'  88  
##  88'     88 d' 88 88'      88  
##  88.     Y8'' .8P 88.      88  
##  Y88888P  Y8888P  Y88888P d88P 
##  oooooooooooooooooooooooooooooo
##                                


#install.packages("vegan")
#install.packages("tidyverse")
#install.packages("ggplot2")

#R version 3.6.0
#install.packages("multcompView")
#install.packages("stringr")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("phyloseq")
#install.packages("remotes")
#BiocManager::install("DESeq2", force=T)
#BiocManager::install("genefilter", force=T)
#install.packages("remotes")
#remotes::install_github("twbattaglia/btools")library(btools)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(multcompView)
library(stringr)
library(vegan)
library(btools)
library(tidyverse)
library(rcompanion)

setwd("C:/Users/Luis/Google Drive/ULM CLASSES/Red Deer paper Nov 2021")
#setwd("D:/Gdrive/ULM CLASSES/Red_Deer_2021")
getwd()
red_deer_map <-import_qiime_sample_data ("R-essentials/mapfile_red_deer_2021_bt.txt")
red_deer_tree <- read_tree_greengenes ("R-essentials/tree_red_deer_2021.nwk")
red_deer_biom <-import_biom("R-essentials/red_deer_biome_taxonomy_2021.biom")
red_deer_biom
red_deer_map$year<-as.factor(red_deer_map$year)
#merge into new object <,

red_deer_bay<- merge_phyloseq (red_deer_biom, red_deer_map, red_deer_tree)
red_deer_bay
sample_data(red_deer_bay)

# Filter ASVs present in controls from entire dataset
controls <- red_deer_bay %>%
  subset_samples(condition== "Control") 
controls@sam_data
controls_1 <- prune_taxa(taxa_sums(controls) > 2, controls)
badtaxa<-taxa_names(controls_1)
alltaxa<-taxa_names(red_deer_bay)
alltaxa1 <- alltaxa[!(alltaxa %in% badtaxa)]
red_deer_bay = prune_taxa(alltaxa1, red_deer_bay)
red_deer_bay <- red_deer_bay %>%
  subset_samples(condition!="Control")
red_deer_bay <- red_deer_bay %>%
  subset_samples(condition!="rm")
red_deer_bay <- red_deer_bay %>%
  subset_samples(condition!="lab")

#Change label for taxonomy
tax <- data.frame(tax_table(red_deer_bay))
tax
tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "d__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)
tax.clean[is.na(tax.clean)] <- ""
tax.clean[1:20,1:7]
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}

####### Fille holes in the tax table
tax.clean[is.na(tax.clean)] <- ""

#### re import as matrix into the S4 object
tax_table(red_deer_bay) <- as.matrix(tax.clean)

tax_table(red_deer_bay)[1:100,1:7]

#change "NA" for ""
red_deer_sub = subset_taxa(red_deer_bay, Phylum!="")
red_deer_sub

#remove 2017 data
red_deer_sub<-prune_samples(sample_sums(red_deer_sub)>=12000,red_deer_sub)

red_deer_sub<- red_deer_sub %>%
  subset_samples(year!="2017")

#Dataframe for ASVs counts
sample_sum_df<- data.frame(sum = sample_sums(red_deer_sub))
sample_sum_df

#plot for ASV counts
asvcounts<-ggplot(sample_sum_df, aes(x = sum)) + geom_histogram()
asvcounts+
  ggtitle("Distribution of sample sequencing depth")+
  xlab("Read counts")+
  ylab("Frequency")+theme_classic(base_size = 16)
rank_names(red_deer_bay)

##### estimate alpha####
red_deer_subr<-rarefy_even_depth(red_deer_sub, sample.size = min(sample_sums(red_deer_sub)),
                              rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
Richness<-estimate_richness(red_deer_subr,measures=c("Observed","Shannon"))
FaithsPD<-estimate_pd(red_deer_subr)
Alpha<-sample_data(red_deer_subr)
Alpha$Observed<-Richness$Observed
Alpha$Shannon<-Richness$Shannon
Alpha$FPD<-FaithsPD$PD
head(Alpha)
Alpha<- Alpha %>%
  subset_samples(FPD!="NA")
Alpha
Alpha$condition<-as.factor(Alpha$condition)
Alpha$year<-as.factor(Alpha$year)
Alpha$betaplus2<-as.factor(Alpha$betaplus2)
Alpha$cond_year<-as.factor(Alpha$cond_year)
Alpha$condition<-factor(Alpha$condition, levels = c("FL","WG","G"))
Alpha$locality<-factor(Alpha$locality, levels = c("FL","Neuhüttenwiese","Riedelhäng","Buchenau","Ahornschachten","Scheuereck","Altschönau"), ordered = TRUE)
###### plot themes####

compositionbetaplus<-theme_bw()+theme(axis.title.x = element_blank(),axis.title.y = element_text(face="bold", size=14),
                                   axis.text.x=element_text(colour="black", face="bold", size=12, hjust = 0.5,angle = 90), 
                                   axis.text.y=element_text(colour="black", size = 14),legend.title = element_blank(),legend.text = element_text(size = 14),legend.position='right',
                                   axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())
compositiontheme<-theme_bw()+theme(axis.title.x = element_blank(),axis.title.y = element_text(face="bold", size=16),
                                   axis.text.x=element_text(colour="black", face="bold", size=16, hjust = 0.5), 
                                   axis.text.y=element_text(colour="black", size = 16),legend.title = element_blank(),legend.text = element_text(size = 16),legend.position='right',
                                   axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())
themebloc<-theme_minimal()+theme(axis.title.x = element_text(color = "black", face="bold", size = 16),axis.title.y = element_text(color="black",face="bold", size=16),  
                                 axis.text.x=element_text(colour="black", face="bold", size=14,angle = 0, hjust = 0.5),axis.text.y=element_text(colour="black", size = 14),
                                 title = element_text(color = "black", face="bold", size = 16),legend.text = element_text(size=12),legend.position='none',strip.text.y=element_text(colour = "Black", face="bold", size = 14),
                                 axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())
themebetaplots<-theme_bw() +theme(axis.title.x = element_text(face="bold", size=14),axis.title.y = element_text(face="bold", size=14),axis.text.x=element_text(colour="black", face="bold", size=14, hjust = 0.5),
                                  axis.text.y=element_text(colour="black", size = 14),legend.title = element_blank(),legend.text = element_text(size = 12),
                                  legend.position='top',axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())
themevolcano<-theme_minimal()+theme(text=element_text(face="bold", size=12),axis.title.x = element_text(face="bold", size=14),axis.title.y = element_text(face="bold", size=14), axis.text.x=element_text(colour="black", face="bold", size=14, hjust = 0.5),
                                  axis.text.y=element_text(colour="black", size = 14),strip.text.x=element_text(face="bold", size=14),legend.title = element_blank(),legend.text = element_text(size = 14),legend.position='right',
                                  axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.background  = element_rect())
##### Colors for plots######
familycolors <-c("Bacteroidales_RF16_group" = "#46b3d1",
                 "Lactobacillaceae" = "#ce8e60",
                 "Bacteroidaceae" = "#887530",
                 "Ruminococcaceae" = "#9f4a71",
                 "Prevotellaceae" = "#5db788",
                 "[Eubacterium]_coprostanoligenes_group" = "#6680c7",
                 "Micrococcaceae" = "#e5818a",
                 "Christensenellaceae" = "#a54a47",
                 "Bacillaceae" = "#507d39",
                 "Rikenellaceae" = "#8e549f",
                 "UCG-010" = "#d38dce",
                 "Lachnospiraceae" = "#d89133",
                 "Planococcocaceae" = "#b3b522",
                 "Pseudomonadaceae" = "#d84480",
                 "Oscillospiraceae" = "#ca5629",
                 "Others (< 3% abund.)" = "#b9c4c3")

phylumcolors <-c("Bacteroidota" = "#5db64c",
                 "Proteobacteria" = "#7865d7",
                 "Firmicutes" = "#c751ba",
                 "Others (< 3% abund.)" = "#b9c4c3")

betaplus2colors<-c("#d1502c",
                   "#616ddb",
                   "#96b233",
                   "#a257cc",
                   "#56b84e",
                   "#d154b5",
                   "#458542",
                   "#d93f7f",
                   "#5dc69b",
                   "#d34250",
                   "#4ab2d8",
                   "#dd892f",
                   "#647ec2",
                   "#d0b040",
                   "#8c5199",
                   "#a0b46c",
                   "#c890db",
                   "#686e29",
                   "#a04a72",
                   "#308e78",
                   "#e385a3",
                   "#a48836",
                   "#b35558",
                   "#e1976e",
                   "#9b5c2d",
                   "Others (< 3% abund.)" = "#b9c4c3",
                   "Others (< 5% abund.)" = "#b9c4c3")
locality_colors <-c("Ahornschachten" = "#bf823b",
                    "Riedelhäng" = "#9473c6",
                    "Buchenau" = "#6ca75d",
                    "Neuhüttenwiese" = "#cc546d",
                    "Altschönau" = "#5641f2",
                    "Scheuereck" = "#6a61ad",
                    "FL" = "#d0c67a")
condition_colors<-c("WG"="#ab72a0",
                    "FL"="#ada458", 
                    "G"="#074eaf")
year_colors<-c("2018"="#72a659",
               "2019"="#c55895",
               "2020"="#c77542",
               "2021"="#777acd")
betaplus_colors<-c("Altschönau 2018"="#6d97ff",
            "Altschönau 2019"="#074eaf",
            "Altschönau 2021"="#6270fb",
            "Scheuereck 2018"="#832895",
            "Scheuereck 2019"="#fcaafe",
            "Scheuereck 2021"="#d556da")
ANCOMcolors<-c("Bacteroidaceae"="#c8577b",
               "Bacteroidales_RF16_group"="#ce5d40",
                 "Christensenellaceae"="#956e34",
               "Eggerthellaceae"="#d5a345",
               "Lachnospiraceae"="#82953f",
               "Oscillospiraceae"="#68b846",
               "p-251-o5"="#4fa877",
               "Prevotellaceae"="#4dadcf",
               "Rikenellaceae"="#7975c9",
               "UCG-010"="#bc5bbc",
               "below cutoff"="grey")
cond_year_colors<-c("G 2018"="#c989a5",
                    "G 2019"="#d846ad",
                    "G 2021"="#c24770",
                    "WG 2018"="#cbd356",
                    "WG 2019"="#889b60",
                    "WG 2021"="#70cf54",
                    "FL 2019"="#007aef",
                    "FL 2020"="#254955",
                    "FL 2021"="#84afc2")


#######Alpha plots by condition#######

Shannon_condition<-ggplot(data=Alpha,aes(x=condition, y=Shannon, fill=condition))
Shannon_condition+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Shannon Index")+ggtitle(" Shannon Alpha diversity per Condition")+scale_fill_manual(values=condition_colors)+
  themebloc
ggsave("plots20211031rarefied/Alpha Shannon by condition.png", dpi = 600, width=6,height=4,units=c("in"))
FPD_condition<-ggplot(data=Alpha,aes(x=condition, y=FPD, fill=condition))
FPD_condition+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Faith's PD")+ggtitle("FPD Alpha diversity per Condition")+scale_fill_manual(values=condition_colors)+
  themebloc
ggsave("plots20211031rarefied/Alpha FPD by condition.png", dpi = 600, width=6,height=4,units=c("in"))
Observed_condition<-ggplot(data=Alpha,aes(x=condition, y=Observed, fill=condition))
Observed_condition+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Number of ASVs")+ggtitle(" ASVs Alpha diversity per Condition")+scale_fill_manual(values=condition_colors)+
  themebloc
ggsave("plots20211031rarefied/Alpha ASVs by condition.png", dpi = 600, width=6,height=4,units=c("in"))
#+facet_grid(cols=vars(year),switch = "x")
######Alpha plots by locality#######
Alpha_local_p<-Alpha%>%
  subset_samples(locality!="0")
Alpha_local<- as.data.frame(Alpha_local_p)
Shannon_locality<-ggplot(data=Alpha_local,aes(x=locality, y=Shannon, fill=locality))
Shannon_locality+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Shannon Index")+ggtitle("Shannon Alpha diversity per locality")+scale_fill_manual(values=locality_colors)+
  themebloc
ggsave("plots20211031rarefied/Alpha Shannon by locality.png", dpi = 600, width=14,height=6.21,units=c("in"))
FPD_locality<-ggplot(data=Alpha_local,aes(x=locality, y=FPD, fill=locality))
FPD_locality+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Faith's PD")+ggtitle("FPD Alpha diversity per locality")+scale_fill_manual(values=locality_colors)+
  themebloc
ggsave("plots20211031rarefied/Alpha FPD by locality.png", dpi = 600, width=14,height=6.21,units=c("in"))
Observed_locality<-ggplot(data=Alpha_local,aes(x=locality, y=Observed, fill=locality))
Observed_locality+geom_boxplot()+geom_jitter(width = 0.25)+ylab("Number of ASVs")+ggtitle("ASVs Alpha diversity per locality")+scale_fill_manual(values=locality_colors)+
  themebloc
ggsave("plots20211031rarefied/Alpha ASVs by locality.png", dpi = 600, width=14,height=6.21,units=c("in"))




##### Stats Alpha for Locality######
Alpha_local_df<-as.list(sample_data(Alpha_local))
Alpha_local_df$locality<-as.factor(Alpha_local_df$locality)
neuorder<-c("FL","Neuhüttenwiese","Riedelhäng","Buchenau","Ahornschachten","Scheuereck","Altschönau")
levels(Alpha_local_df$locality)<-neuorder
summary(Alpha_local_df$locality)
summary(Alpha_local_df$FPD)
kw_PD<-kruskal.test(FPD~locality, data= Alpha_local_df)
pw_kw_PD<-pairwise.wilcox.test(Alpha_local_df$FPD, Alpha_local_df$locality,
                               p.adjust.method = "fdr")
PT_pw_kw_PD<-pw_kw_PD$p.value
PT_pw_kw_PD
PT_pw_kw_PD1<-fullPTable(PT_pw_kw_PD)

pdletters<-multcompLetters(PT_pw_kw_PD1,
                           compare="<",
                           threshold=0.05,  # p-value to use as significance threshold
                           Letters=letters,
                           reversed = FALSE)

kw_Shannon<-kruskal.test(Shannon~locality, data=Alpha_local_df)
pw_kw_Shannon<-pairwise.wilcox.test(Alpha_local_df$Shannon, Alpha_local_df$locality,
                                    p.adjust.method = "fdr")
PT_pw_kw_Shannon<-pw_kw_Shannon$p.value
PT_pw_kw_Shannon
PT_pw_kw_Shannon1<-fullPTable(PT_pw_kw_Shannon)

shannonletters<-multcompLetters(PT_pw_kw_Shannon1,
                                compare="<",
                                threshold=0.05,  # p-value to use as significance threshold
                                Letters=letters,
                                reversed = FALSE)

kw_Observed<-kruskal.test(Observed~locality, data=Alpha_local_df)
pw_kw_Observed<-pairwise.wilcox.test(Alpha_local_df$Observed, Alpha_local_df$locality,
                                     p.adjust.method = "fdr")
PT_pw_kw_Observed<-pw_kw_Observed$p.value
PT_pw_kw_Observed
PT_pw_kw_Observed1<-fullPTable(PT_pw_kw_Observed)

observedletters<-multcompLetters(PT_pw_kw_Observed1,
                                 compare="<",
                                 threshold=0.05,  # p-value to use as significance threshold
                                 Letters=letters,
                                 reversed = FALSE)
kw_Shannon
shannonletters
kw_Observed
observedletters
kw_PD
pdletters

Alpha_cond_df<-as.list(sample_data(Alpha))
Alpha_cond_df$condition<-as.factor(Alpha_cond_df$condition)
neuorder_cond<-c("FL","WG","G")
levels(Alpha_cond_df$condition)<-neuorder_cond
summary(Alpha_cond_df$condition)
summary(Alpha_cond_df$FPD)
kw_condition_PD<-kruskal.test(FPD~condition, data= Alpha_cond_df)
pw_kw_condition_PD<-pairwise.wilcox.test(Alpha_cond_df$FPD, Alpha_cond_df$condition,
                                         p.adjust.method = "fdr")
PT_pw_kw_condition_PD<-pw_kw_condition_PD$p.value
PT_pw_kw_condition_PD
PT_pw_kw_condition_PD1<-fullPTable(PT_pw_kw_condition_PD)

pdletters_cond<-multcompLetters(PT_pw_kw_condition_PD1,
                                compare="<",
                                threshold=0.05,  # p-value to use as significance threshold
                                Letters=letters,
                                reversed = FALSE)

kw_condition_Shannon<-kruskal.test(Shannon~condition, data=Alpha_cond_df)
pw_kw_condition_Shannon<-pairwise.wilcox.test(Alpha_cond_df$Shannon, Alpha_cond_df$condition,
                                              p.adjust.method = "fdr")
PT_pw_kw_condition_Shannon<-pw_kw_condition_Shannon$p.value
PT_pw_kw_condition_Shannon
PT_pw_kw_condition_Shannon1<-fullPTable(PT_pw_kw_condition_Shannon)

shannonletters_cond<-multcompLetters(PT_pw_kw_condition_Shannon1,
                                     compare="<",
                                     threshold=0.05,  # p-value to use as significance threshold
                                     Letters=letters,
                                     reversed = FALSE)

kw_condition_Observed<-kruskal.test(Observed~condition, data=Alpha_cond_df)
pw_kw_condition_Observed<-pairwise.wilcox.test(Alpha_cond_df$Observed, Alpha_cond_df$condition,
                                               p.adjust.method = "fdr")
PT_pw_kw_condition_Observed<-pw_kw_condition_Observed$p.value
PT_pw_kw_condition_Observed
PT_pw_kw_condition_Observed1<-fullPTable(PT_pw_kw_condition_Observed)

observedletters_cond<-multcompLetters(PT_pw_kw_condition_Observed1,
                                      compare="<",
                                      threshold=0.05,  # p-value to use as significance threshold
                                      Letters=letters,
                                      reversed = FALSE)
kw_condition_Shannon
shannonletters_cond
kw_condition_Observed
observedletters_cond
kw_condition_PD
pdletters_cond

### Composition by Family locality####

melt_Family_ord <- red_deer_sub %>%
  tax_glom(taxrank = "Family") %>%      
  psmelt() %>%                         
  arrange(Family) 
melt_Family<- melt_Family_ord%>%
  arrange(locality) %>%
  mutate(locality=factor(locality, levels=c("FL","Neuhüttenwiese","Riedelhäng","Buchenau","Ahornschachten","Scheuereck","Altschönau")))                    

melt_Family2 <- aggregate(Abundance ~ Family + locality, 
                          data= melt_Family, 
                          sum)

melt_Family2$Family <- as.character(melt_Family2$Family) 

melt_Family3 <- melt_Family2 %>% 
  dplyr::group_by(locality) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),3)) %>%
  dplyr::filter(rel.freq > 0.03)

#Get Remainers to add in graph
remainers_Family <- melt_Family2 %>% 
  dplyr::group_by(locality) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),3)) %>%
  dplyr::filter(rel.freq < 0.03)

remainers_Family2 <- aggregate(rel.freq ~ locality, 
                               data= remainers_Family, 
                               sum)

remainers_Family2$Family <- "Others (< 3% abund.)"

join_Family <- full_join(melt_Family3,remainers_Family2)

join_Family <- join_Family %>%
  dplyr::mutate(percent = paste0(round(100 * rel.freq, 0), "%"))

join_Family$Family <- as.factor(join_Family$Family)
join_Family$locality<- as.factor(join_Family$locality)
join_Family$Family <- reorder(join_Family$Family, join_Family$Abundance)
family_species <- ggplot(join_Family, aes(x = locality, y = rel.freq, fill = Family)) + 
  geom_bar(stat = "identity", position = "fill") +
  ylab("Relative abundance") + scale_fill_manual(values=familycolors)+
  xlab("locality") +ggtitle("composition by locality")+
  geom_text(aes(label=rel.freq), position=position_fill(vjust=0.5), size=4)+
  compositiontheme
family_species

ggsave("plots20211031rarefied/Taxa composition family level all sites.png", dpi = 600, width=14,height=6.21,units=c("in"))

###Composition by Phylum all site####

melt_Phylum_ord <- red_deer_sub %>%
  tax_glom(taxrank = "Phylum") %>%      
  psmelt() %>%                         
  arrange(Phylum) 
melt_Phylum<- melt_Phylum_ord%>%
  arrange(locality) %>%
  mutate(locality=factor(locality, levels=c("FL","Neuhüttenwiese","Riedelhäng","Buchenau","Ahornschachten","Scheuereck","Altschönau")))

melt_Phylum2 <- aggregate(Abundance ~ Phylum + locality, 
                          data= melt_Phylum, 
                          sum)

melt_Phylum2$Phylum <- as.character(melt_Phylum2$Phylum) 

melt_Phylum3 <- melt_Phylum2 %>% 
  dplyr::group_by(locality) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),3)) %>%
  dplyr::filter(rel.freq > 0.03)

#Get Remainers to add in graph
remainers_Phylum <- melt_Phylum2 %>% 
  dplyr::group_by(locality) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),3)) %>%
  dplyr::filter(rel.freq < 0.03)

remainers_Phylum2 <- aggregate(rel.freq ~ locality, 
                               data= remainers_Phylum, 
                               sum)

remainers_Phylum2$Phylum <- "Others (< 3% abund.)"

join_Phylum <- full_join(melt_Phylum3,remainers_Phylum2)

join_Phylum <- join_Phylum %>%
  dplyr::mutate(percent = paste0(round(100 * rel.freq, 0), "%"))

join_Phylum$Phylum <- as.factor(join_Phylum$Phylum)
join_Phylum$locality <- as.factor(join_Phylum$locality)
join_Phylum$Phylum <- reorder(join_Phylum$Phylum, join_Phylum$Abundance)
phylum_species <- ggplot(join_Phylum, aes(x = locality, y = rel.freq, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  ylab("Relative abundance")+ scale_fill_manual(values = phylumcolors) +
  xlab("locality") +ggtitle("composition by locality")+
  geom_text(aes(label=rel.freq), position=position_fill(vjust=0.5), size=6)+
  compositiontheme
phylum_species
ggsave("plots20211031rarefied/Taxa composition phylum level all sites.png", dpi = 600, width=14,height=6.21,units=c("in"))


####composition plots for condition
red_deer_sub_cond<- red_deer_sub

###Composition by Phylum condition####
melt_Phylum_cond_ord <- red_deer_sub_cond %>%
  tax_glom(taxrank = "Phylum") %>%      
  psmelt() %>%                         
  arrange(Phylum) 
melt_Phylum_cond<- melt_Phylum_cond_ord%>%
  arrange(condition) %>%
  mutate(condition=factor(condition, levels=c("FL","WG","G")))                     
melt_Phylum_cond

melt_Phylum_cond2 <- aggregate(Abundance ~ Phylum + condition, 
                               data= melt_Phylum_cond, 
                               sum)

melt_Phylum_cond2$Phylum <- as.character(melt_Phylum_cond2$Phylum) 

melt_Phylum_cond3 <- melt_Phylum_cond2 %>% 
  dplyr::group_by(condition) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),3)) %>%
  dplyr::filter(rel.freq > 0.03)

#Get Remainers to add in graph
remainers_Phylum_cond <- melt_Phylum_cond2 %>% 
  dplyr::group_by(condition) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),3)) %>%
  dplyr::filter(rel.freq < 0.03)

remainers_Phylum_cond2 <- aggregate(rel.freq ~ condition, 
                                    data= remainers_Phylum_cond, 
                                    sum)

remainers_Phylum_cond2$Phylum <- "Others (< 3% abund.)"

join_Phylum_cond <- full_join(melt_Phylum_cond3,remainers_Phylum_cond2)

join_Phylum_cond <- join_Phylum_cond %>%
  dplyr::mutate(percent = paste0(round(100 * rel.freq, 0), "%"))

join_Phylum_cond$Phylum <- as.factor(join_Phylum_cond$Phylum)
join_Phylum_cond$condition<-as.factor(join_Phylum_cond$condition)
join_Phylum_cond$condition
#levels(join_Phylum_cond$condition)<-c("FL","WG","G")
join_Phylum_cond$Phylum <- reorder(join_Phylum_cond$Phylum, join_Phylum_cond$Abundance)

phylum_species_cond_pt <- ggplot(join_Phylum_cond, aes(x = condition, y = rel.freq, fill = Phylum))+
  geom_bar(stat = "identity", position = "fill")+ylab("Relative abundance")+ scale_fill_manual(values = phylumcolors) +
  xlab("Management condition")+ggtitle("Phylum Composition by Management Condition")+
  geom_text(aes(label=rel.freq), position=position_fill(vjust=0.5), size=6)+compositiontheme
phylum_species_cond_pt
ggsave("plots20211031rarefied/Taxa composition phylum by condition.png", dpi = 600, width=14,height=6.21,units=c("in"))

### Composition by Family Condition####
melt_Family_cond_ord <- red_deer_sub_cond %>%
  tax_glom(taxrank = "Family") %>%      
  psmelt() %>%                         
  arrange(Family) 
melt_Family_cond<- melt_Family_cond_ord%>%
  arrange(condition) %>%
  mutate(condition=factor(condition, levels=c("FL","WG","G")))                     

melt_Family_cond2 <- aggregate(Abundance ~ Family + condition, 
                               data= melt_Family_cond, 
                               sum)

melt_Family_cond2$Family <- as.character(melt_Family_cond2$Family) 

melt_Family_cond3 <- melt_Family_cond2 %>% 
  dplyr::group_by(condition) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),3)) %>%
  dplyr::filter(rel.freq > 0.03)

#Get Remainers to add in graph
remainers_Family_cond <- melt_Family_cond2 %>% 
  dplyr::group_by(condition) %>%
  dplyr::mutate(rel.freq = round(Abundance / sum(Abundance),3)) %>%
  dplyr::filter(rel.freq < 0.03)

remainers_Family_cond2 <- aggregate(rel.freq ~ condition, 
                                    data= remainers_Family_cond, 
                                    sum)

remainers_Family_cond2$Family <- "Others (< 3% abund.)"

join_Family_cond <- full_join(melt_Family_cond3,remainers_Family_cond2)
join_Family_cond
join_Family_cond <- join_Family_cond %>%
  dplyr::mutate(percent = paste0(round(100 * rel.freq, 0), "%"))

join_Family_cond$Family <- as.factor(join_Family_cond$Family)
join_Family_cond$condition<-as.factor(join_Family_cond$condition)
#levels(join_Family_cond$condition)
#levels(join_Family_cond$condition)<-c("FL","WG","G")
join_Family_cond$Family <- reorder(join_Family_cond$Family, join_Family_cond$Abundance)
family_speciesGbt <- ggplot(join_Family_cond, aes(x = condition, y = rel.freq, fill = Family)) + 
  geom_bar(stat = "identity", position = "fill") +
  ylab("Relative abundance") + scale_fill_manual(values=familycolors)+
  xlab("Management Condition") +ggtitle("Family Composition by Management Condition")+
  geom_text(aes(label=rel.freq), position=position_fill(vjust=0.5), size=4)+
  compositiontheme
family_speciesGbt

ggsave("plots20211031rarefied/Taxa composition family by condition.png", dpi = 600, width=14,height=6.21,units=c("in"))


#####Beta Diversity#######

#remove ALL singletons
red_deer_sub1.1<- filter_taxa(red_deer_subr, function (x) {sum(x > 0) >1}, prune=TRUE)
red_deer_sub1.1

otu_table(red_deer_sub1.1)


prevalencedf = apply(X = otu_table(red_deer_sub1.1),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})

prevalencedf = data.frame(Prevalence = prevalencedf,
                          TotalAbundance = taxa_sums(red_deer_sub1.1),
                          tax_table(red_deer_sub1.1))
head(prevalencedf)



nsamples(red_deer_sub1.1)
#so we want the ids to be present in at least 30% of samples, roughtly 200
prevalenceThreshold = 200 #changed to 200

prevalenceThreshold

red_deer_sub1.1<- subset_taxa(red_deer_sub1.1, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
red_deer_sub1.1


plyr::ddply(prevalencedf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})

prevalencedf[1:10,]

prevalencedf1 = subset(prevalencedf, Phylum %in% get_taxa_unique(red_deer_sub1.1, taxonomic.rank = "Phylum"))
prevalencedf1[1:10,]

ggplot(prevalencedf1, aes(TotalAbundance, Prevalence / nsamples(red_deer_sub1.1),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

str(prevalencedf)
str(prevalencedf1)

#Now we use this value to filter out any otu that is not present in at least 30% of all samples
(prevalencedf1$Prevalence >= prevalenceThreshold)

keepTaxa = rownames(prevalencedf1)[(prevalencedf1$Prevalence >= prevalenceThreshold)]
length(keepTaxa)

red_deer_sub1 = prune_taxa(keepTaxa, red_deer_sub1.1)
red_deer_sub1

red_deer_sub1<-prune_samples(sample_sums(red_deer_sub1)>=2000,red_deer_sub1)

# weighted unifrac 
DistW = distance(red_deer_sub1,method="wunifrac")
#unweighted unifrac 
DistUW = distance(red_deer_sub1,method="uunifrac")
DistUW

#create ordination

ordW = ordinate(red_deer_sub1, method = "PCoA", distance = DistW)
ordUW = ordinate(red_deer_sub1, method = "PCoA", distance = DistUW)


#lets visualise how informative each ordination is for each distance matrix
plot_scree(ordW)
plot_scree(ordUW)

#plot ordinations for years
plot_ordination(red_deer_sub1, ordW, color = "year",axes=c(2,1))+  stat_ellipse() +  ggtitle("Weighted Unifrac by year")+
  theme_bw()+themebetaplots+scale_color_manual(values=year_colors)
ggsave("plots20211031rarefied/PCoA Weighted UniFrac by year.png", dpi =600, width=14.23,height=6.21,units=c("in"))
plot_ordination(red_deer_sub1, ordUW, color = "year",axes=c(2,1))+  stat_ellipse() +  ggtitle("UnWeighted UniFrac by year")+
  theme_bw()+themebetaplots+scale_color_manual(values=year_colors)
ggsave("plots20211031rarefied/PCoA Un-Weighted UniFrac by year.png", dpi = 600, width=14.23,height=6.21,units=c("in"))

plot_ordination(red_deer_sub1, ordW, color = "condition",axes=c(2,1))+  stat_ellipse() +  ggtitle("Weighted Unifrac by Management Condition")+
  theme_bw()+themebetaplots+scale_color_manual(values=condition_colors)
ggsave("plots20211031rarefied/PCoA Weighted UniFrac by condition.png", dpi =600, width=14.23,height=6.21,units=c("in"))
plot_ordination(red_deer_sub1, ordUW, color = "condition",axes=c(2,1))+  stat_ellipse() +  ggtitle("UnWeighted UniFrac by Management Condition")+
  theme_bw()+themebetaplots+scale_color_manual(values=condition_colors)
ggsave("plots20211031rarefied/PCoA Un-Weighted UniFrac by condition.png", dpi = 600, width=14.23,height=6.21,units=c("in"))

plot_ordination(red_deer_sub1, ordW, color = "locality",axes=c(2,1))+  stat_ellipse() +  ggtitle("Weighted Unifrac by Management locality")+
  theme_bw()+themebetaplots+scale_color_manual(values=locality_colors)
ggsave("plots20211031rarefied/PCoA Weighted UniFrac by locality.png", dpi =600, width=14.23,height=6.21,units=c("in"))
plot_ordination(red_deer_sub1, ordUW, color = "locality",axes=c(2,1))+  stat_ellipse() +  ggtitle("UnWeighted UniFrac by Management locality")+
  theme_bw()+themebetaplots+scale_color_manual(values=locality_colors)
ggsave("plots20211031rarefied/PCoA Un-Weighted UniFrac by locality.png", dpi = 600, width=14.23,height=6.21,units=c("in"))

####Centroid plots#####

#####Beta for Localities########
#UW Unifrac
metadata<-data.frame(sample_data(red_deer_sub1))

UW_ab<-data.frame(ordUW$vectors[,1],
                 ordUW$vectors[,2])
colnames(UW_ab)[1]<-"MDS1"
colnames(UW_ab)[2]<-"MDS2"

UWcentroid<-cbind(UW_ab,metadata)
centroids_UWy <- as.data.frame(UWcentroid %>% 
                                 dplyr::group_by(locality) %>% # calculate functions below for each group
                                 dplyr::summarise(mean_MDS1=mean(MDS1),
                                                  mean_MDS2=mean(MDS2),
                                                  n_MDS1=length(MDS1),
                                                  n_MDS2=length(MDS2),
                                                  stdv_MDS1=sd(MDS1),
                                                  stdv_MDS2=sd(MDS2),
                                                  se_MDS1=stdv_MDS1/sqrt(n_MDS1),
                                                  se_MDS2=stdv_MDS2/sqrt(n_MDS2)))
UW_unifrac_plot_by_locality <- ggplot(data = centroids_UWy, aes(x=mean_MDS1, y=mean_MDS2, color=locality))+
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mean_MDS2-se_MDS2, 
                    ymax=mean_MDS2+se_MDS2), width=0.02,  size=0.3)+
  geom_errorbarh(aes(xmin=mean_MDS1-se_MDS1, 
                     xmax=mean_MDS1+se_MDS1), height=0.02, size=0.3) +
  geom_point(data = UWcentroid, aes(x=MDS1, y=MDS2, color=locality), alpha=0.5, size=3) +
  scale_color_manual(values = locality_colors)+
  labs(x = "MDS1 [24%]", y="MDS2 [7.7%]")+
  ggtitle("Unweighted UniFrac Locality")+themebetaplots
UW_unifrac_plot_by_locality+stat_ellipse(data=UWcentroid, aes(x=MDS1, y=MDS2, color=locality),inherit.aes = FALSE)
ggsave("plots20211031rarefied/UW_unifrac_plot_by_locality_ctr.png", dpi = 600, width=10.23,height=6.21,units=c("in"))

#WU unifrac
metadata<-data.frame(sample_data(red_deer_sub1))

W_ab<-data.frame(ordW$vectors[,1],
                   ordW$vectors[,2])
colnames(W_ab)[1]<-"MDS1"
colnames(W_ab)[2]<-"MDS2"

Wcentroid<-cbind(W_ab,metadata)

centroids_Wy <- as.data.frame(Wcentroid %>% 
                                dplyr::group_by(locality) %>% # calculate functions below for each group
                                dplyr::summarise(mean_MDS1=mean(MDS1),
                                                 mean_MDS2=mean(MDS2),
                                                 n_MDS1=length(MDS1),
                                                 n_MDS2=length(MDS2),
                                                 stdv_MDS1=sd(MDS1),
                                                 stdv_MDS2=sd(MDS2),
                                                 se_MDS1=stdv_MDS1/sqrt(n_MDS1),
                                                 se_MDS2=stdv_MDS2/sqrt(n_MDS2)))
W_unifrac_plot_by_locality <- ggplot(data = centroids_Wy, aes(x=mean_MDS1, y=mean_MDS2, color=locality))+
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mean_MDS2-se_MDS2, 
                    ymax=mean_MDS2+se_MDS2), width=0.03,  size=0.5)+
  geom_errorbarh(aes(xmin=mean_MDS1-se_MDS1, 
                     xmax=mean_MDS1+se_MDS1), height=0.03, size=0.5) +
  geom_point(data = Wcentroid, aes(x=MDS1, y=MDS2, color=locality), alpha=0.5, size=3) +
  scale_color_manual(values = locality_colors)+
  labs(x = "MDS1 [28.6%]", y="MDS2 [14%]")+
  ggtitle("Weighted UniFrac Locality") +themebetaplots
W_unifrac_plot_by_locality+stat_ellipse(data=Wcentroid, aes(x=MDS1, y=MDS2, color=locality),inherit.aes = FALSE)
ggsave("plots20211031rarefied/W_unifrac_plot_by_locality_ctr.png", dpi =600, width=10.23,height=6.21,units=c("in"))


#####Beta for management########

#UW Unifrac for management conditions

metadata<-data.frame(sample_data(red_deer_sub1))

UW_ab<-data.frame(ordUW$vectors[,1],
                 ordUW$vectors[,2])
colnames(UW_ab)[1]<-"MDS1"
colnames(UW_ab)[2]<-"MDS2"

UWcentroid<-cbind(UW_ab,metadata)
centroids_UWy <- as.data.frame(UWcentroid %>% 
                                 dplyr::group_by(condition) %>% # calculate functions below for each group
                                 dplyr::summarise(mean_MDS1=mean(MDS1),
                                                  mean_MDS2=mean(MDS2),
                                                  n_MDS1=length(MDS1),
                                                  n_MDS2=length(MDS2),
                                                  stdv_MDS1=sd(MDS1),
                                                  stdv_MDS2=sd(MDS2),
                                                  se_MDS1=stdv_MDS1/sqrt(n_MDS1),
                                                  se_MDS2=stdv_MDS2/sqrt(n_MDS2)))
UW_unifrac_plot_by_condition <- ggplot(data = centroids_UWy, aes(x=mean_MDS1, y=mean_MDS2, color=condition))+
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mean_MDS2-se_MDS2, 
                    ymax=mean_MDS2+se_MDS2), width=0.02,  size=0.3)+
  geom_errorbarh(aes(xmin=mean_MDS1-se_MDS1, 
                     xmax=mean_MDS1+se_MDS1), height=0.02, size=0.3) +
  geom_point(data = UWcentroid, aes(x=MDS1, y=MDS2, color=condition), alpha=0.5, size=3) +
  scale_color_manual(values = condition_colors)+
  labs(x = "MDS1 [23%]", y="MDS2 [7.6%]")+
  ggtitle("Unweighted UniFrac by Management Condition") +themebetaplots
UW_unifrac_plot_by_condition+stat_ellipse(data=UWcentroid, aes(x=MDS1, y=MDS2, color=condition),inherit.aes = FALSE)
ggsave("plots20211031rarefied/UW_unifrac_plot_by_condition_ctr.png", dpi = 600, width=10.23,height=6.21,units=c("in"))

#WU unifrac for management conditions
metadata<-data.frame(sample_data(red_deer_sub1))

W_ab<-data.frame(ordW$vectors[,1],
                  ordW$vectors[,2])
colnames(W_ab)[1]<-"MDS1"
colnames(W_ab)[2]<-"MDS2"

Wcentroid<-cbind(W_ab,metadata)

centroids_Wy <- as.data.frame(Wcentroid %>% 
                                dplyr::group_by(condition) %>% # calculate functions below for each group
                                dplyr::summarise(mean_MDS1=mean(MDS1),
                                                 mean_MDS2=mean(MDS2),
                                                 n_MDS1=length(MDS1),
                                                 n_MDS2=length(MDS2),
                                                 stdv_MDS1=sd(MDS1),
                                                 stdv_MDS2=sd(MDS2),
                                                 se_MDS1=stdv_MDS1/sqrt(n_MDS1),
                                                 se_MDS2=stdv_MDS2/sqrt(n_MDS2)))
W_unifrac_plot_by_condition <- ggplot(data = centroids_Wy, aes(x=mean_MDS1, y=mean_MDS2, color=condition))+
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mean_MDS2-se_MDS2, 
                    ymax=mean_MDS2+se_MDS2), width=0.03,  size=0.5)+
  geom_errorbarh(aes(xmin=mean_MDS1-se_MDS1, 
                     xmax=mean_MDS1+se_MDS1), height=0.03, size=0.5) +
  geom_point(data = Wcentroid, aes(x=MDS1, y=MDS2, color=condition), alpha=0.5, size=3) +
  scale_color_manual(values = condition_colors)+
  labs(x = "MDS1 [32.6%]", y="MDS2 [12.3%]")+
  ggtitle("Weighted UniFrac by Management Condition") +themebetaplots
W_unifrac_plot_by_condition +stat_ellipse(data=Wcentroid, aes(x=MDS1, y=MDS2, color=condition),inherit.aes = FALSE)
ggsave("plots20211031rarefied/W_unifrac_plot_by_condition_ctr.png", dpi =600, width=10.23,height=6.21,units=c("in"))

########Beta Diversity Cond_year######
metadata<-data.frame(sample_data(red_deer_sub1))
plot_ordination(red_deer_sub1, ordW, color = "cond_year",axes=c(2,1))+  stat_ellipse() +  ggtitle("Weighted Unifrac by year")+
  theme_bw()+themebetaplots+scale_color_manual(values=cond_year_colors)
ggsave("plots20211031rarefied/PCoA Weighted UniFrac by cond_year.png", dpi =600, width=14.23,height=6.21,units=c("in"))

plot_ordination(red_deer_sub1, ordUW, color = "cond_year",axes=c(2,1))+  stat_ellipse() +  ggtitle("UnWeighted UniFrac by year")+
  theme_bw()+themebetaplots+scale_color_manual(values=cond_year_colors)
ggsave("plots20211031rarefied/PCoA Un-Weighted UniFrac by cond_year.png", dpi = 600, width=14.23,height=6.21,units=c("in"))
cond_year_colors<-c("G 2018"="#c989a5",
                    "G 2019"="#d846ad",
                    "G 2021"="#c24770",
                    "WG 2018"="#cbd356",
                    "WG 2019"="#889b60",
                    "WG 2021"="#70cf54",
                    "FL 2019"="#007aef",
                    "FL 2020"="#254955",
                    "FL 2021"="#84afc2")

#UW unifrac for Cond_year
UW_cond_year_ab<-data.frame(ordUW$vectors[,1],
                            ordUW$vectors[,2])
colnames(UW_cond_year_ab)[1]<-"MDS1"
colnames(UW_cond_year_ab)[2]<-"MDS2"

UWcentroid_cond_year<-cbind(UW_cond_year_ab,metadata)
centroids_UW_cond_year <- as.data.frame(UWcentroid_cond_year %>% 
                                          dplyr::group_by(cond_year) %>% # calculate functions below for each group
                                          dplyr::summarise(mean_MDS1=mean(MDS1),
                                                           mean_MDS2=mean(MDS2),
                                                           n_MDS1=length(MDS1),
                                                           n_MDS2=length(MDS2),
                                                           stdv_MDS1=sd(MDS1),
                                                           stdv_MDS2=sd(MDS2),
                                                           se_MDS1=stdv_MDS1/sqrt(n_MDS1),
                                                           se_MDS2=stdv_MDS2/sqrt(n_MDS2)))
UW_cond_year_unifrac_plot <- ggplot(data = centroids_UW_cond_year, aes(x=mean_MDS1, y=mean_MDS2, color=cond_year))+
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mean_MDS2-se_MDS2, 
                    ymax=mean_MDS2+se_MDS2), width=0.02,  size=0.3)+
  geom_errorbarh(aes(xmin=mean_MDS1-se_MDS1, 
                     xmax=mean_MDS1+se_MDS1), height=0.02, size=0.3) +
  geom_point(data = UWcentroid_cond_year, aes(x=MDS1, y=MDS2, color=cond_year), alpha=0.5, size=3) +
  labs(x = "MDS1 [24%]", y="MDS2 [7%]")+ scale_color_manual(values = cond_year_colors)+
  ggtitle("Unweighted UniFrac cond_year")+themebetaplots
UW_cond_year_unifrac_plot+stat_ellipse(data=UWcentroid_cond_year, aes(x=MDS1, y=MDS2, color=cond_year),inherit.aes = FALSE)
ggsave("plots20211031rarefied/UW_cond_year_unifrac_plot.png", dpi = 600, width=10.23,height=6.21,units=c("in"))

#WU unifrac for Cond_year
metadata<-data.frame(sample_data(red_deer_sub1))

W_cond_year_ab<-data.frame(ordW$vectors[,1],
                           ordW$vectors[,2])
colnames(W_cond_year_ab)[1]<-"MDS1"
colnames(W_cond_year_ab)[2]<-"MDS2"

Wcentroid_cond_year<-cbind(W_cond_year_ab,metadata)
centroids_W_cond_year <- as.data.frame(Wcentroid_cond_year %>% 
                                         dplyr::group_by(cond_year) %>% # calculate functions below for each group
                                         dplyr::summarise(mean_MDS1=mean(MDS1),
                                                          mean_MDS2=mean(MDS2),
                                                          n_MDS1=length(MDS1),
                                                          n_MDS2=length(MDS2),
                                                          stdv_MDS1=sd(MDS1),
                                                          stdv_MDS2=sd(MDS2),
                                                          se_MDS1=stdv_MDS1/sqrt(n_MDS1),
                                                          se_MDS2=stdv_MDS2/sqrt(n_MDS2)))
W_cond_year_unifrac_plot<- ggplot(data = centroids_W_cond_year, aes(x=mean_MDS1, y=mean_MDS2, color=cond_year))+
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mean_MDS2-se_MDS2, 
                    ymax=mean_MDS2+se_MDS2), width=0.03,  size=0.5)+
  geom_errorbarh(aes(xmin=mean_MDS1-se_MDS1, 
                     xmax=mean_MDS1+se_MDS1), height=0.03, size=0.5) +
  geom_point(data = Wcentroid_cond_year, aes(x=MDS1, y=MDS2, color=cond_year), alpha=0.5, size=3)+
  labs(x = "MDS1 [28.6%]", y="MDS2 [14%]")+ scale_color_manual(values = cond_year_colors)+
  ggtitle("Weighted UniFrac cond_year") +themebetaplots
W_cond_year_unifrac_plot+stat_ellipse(data=Wcentroid_cond_year, aes(x=MDS1, y=MDS2, color=cond_year),inherit.aes = FALSE)
ggsave("plots20211031rarefied/W_cond_year_unifrac_plot.png", dpi =600, width=10.23,height=6.21,units=c("in"))





######PERMANOVAS AND BETADISPER CONDITION and YEAR######

metadata_condition<-data.frame(sample_data(red_deer_sub1))
metadata_condition$condition<-as.factor(metadata_condition$condition)
levels(metadata_condition$condition)<-c("FL","WG","G")

#permanova UW
UW_per_condition<-adonis2(DistUW~condition*year, data=metadata_condition, permutations = 9999)
UW_per_condition
#permanova W
W_per_condition<-adonis2(DistW~condition*year, data=metadata_condition, permutations = 9999)
W_per_condition

model_condition_UW <- betadisper(DistUW, metadata_condition$condition, type = "centroid")
permu_model_condition_UW<-permutest(model_condition_UW, permutations = 9999)
anova_model_condition_UW<-anova(model_condition_UW)
TU_model_condition_UW<-TukeyHSD(model_condition_UW)
model_condition_UW
permu_model_condition_UW
anova_model_condition_UW
TU_model_condition_UW
boxplot(model_condition_UW)
png(file="plots20211031rarefied/boxplot_betadispUW.png",
    width=300, height=300)
boxplot(model_condition_UW)
dev.off()
plot(TU_model_condition_UW)
png(file="plots20211031rarefied/betadispUW.png",
    width=300, height=300)
plot(TU_model_condition_UW)
dev.off()



#dispersion for UW
model_condition_W <- betadisper(DistW, metadata_condition$condition, type = "centroid")
permu_model_condition_W<-permutest(model_condition_W, permutations = 9999)
anova_model_condition_W<-anova(model_condition_W)
TU_model_condition_W<-TukeyHSD(model_condition_W)
model_condition_W
permu_model_condition_W
anova_model_condition_W
TU_model_condition_W
boxplot(model_condition_W)
png(file="plots20211031rarefied/boxplot_betadispW.png",
    width=300, height=300)
boxplot(model_condition_W)
dev.off()
plot(TU_model_condition_W)
png(file="plots20211031rarefied/betadispW.png",
    width=300, height=300)
plot(TU_model_condition_W)
dev.off()

#####ANCOM for condition######
library(nlme)
library(tidyverse)
library(ggplot2)
library(compositions)
source("R-essentials/ancom_v2.1.R.txt")
#### for more information on ANCOM see https://github.com/FrederickHuangLin/ANCOM

##removing non-id phylum

red_deer_sub1ANCOM<- subset_taxa(red_deer_sub, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))


#prepare objects for data pre-processing for ANCOM


taxanxom_condition<-otu_table(red_deer_sub1ANCOM)
metanxom_condition<-sample_data(red_deer_sub1ANCOM)

metanxom_condition

##Check for the structure
colnames(taxanxom_condition)
colnames(metanxom_condition)

#pre-processing of data for ANCOM

preprocessancom_condition<-feature_table_pre_process(taxanxom_condition, metanxom_condition, "sample.id", group_var = "condition", out_cut = 0.0, zero_cut = 0.90, 10000, neg_lb = TRUE)

#objects created by the preprocessing

preprocessancom_condition$feature_table
preprocessancom_condition$meta_data
preprocessancom_condition$structure_zeros

#Now we run the ANCOM

resultsANCOM_condition<-ANCOM(preprocessancom_condition$feature_table,preprocessancom_condition$meta_data,preprocessancom_condition$structure_zeros,main_var = "condition",alpha=0.05,adj_formula = NULL,p_adj_method = "bonferroni")

resultsANCOM_condition

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(preprocessancom_condition$structure_zeros), nrow(preprocessancom_condition$feature_table), sum(apply(preprocessancom_condition$structure_zeros, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann_condition = data.frame(x = min(resultsANCOM_condition$fig$data$x), y = cut_off["detected_0.9"], label = "W[0.9]")

fig_condition = resultsANCOM_condition$fig +  
  geom_hline(yintercept = cut_off["detected_0.9"], linetype = "dashed") + 
  geom_text(data = dat_ann_condition, aes(x = x, y = y, label = label), 
            size = 6, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig_condition
full_comparison_condition<-as.data.frame(resultsANCOM_condition$out)
write.table(resultsANCOM_condition$out, "plots20211031rarefied/resultsANCOM_condition.csv")

#extract taxonomy table from phyloseq object
taxnames_condition <- tax_table(red_deer_sub1ANCOM)
#coerce table to a dataframe and keep the rownames because they are the OTU IDs
taxnames_condition_df <- as.data.frame(taxnames_condition@.Data, row.names = row.names(taxnames_condition))
str(taxnames_condition_df)

#change rownames into a proper column
taxnames_condition_df_clean <- tibble::rownames_to_column(taxnames_condition_df, "taxa_id")
werty_condition<-as.data.frame(resultsANCOM_condition$fig$data)
head(werty_condition)
#merge the taxa df with the output of the ancom test:
ANCOMdifferential_condition <- merge(full_comparison_condition, taxnames_condition_df_clean, by = "taxa_id")
ANCOMdifferential_condition <- merge(ANCOMdifferential_condition,werty_condition, by="taxa_id")
head(ANCOMdifferential_condition)
write.table(ANCOMdifferential_condition, "plots20211031rarefied/ANCOM_differential_condition.csv", sep = ",")


head(ANCOMdifferential_condition)
library(ggrepel)


ANCOMdiff_condition2<-ANCOMdifferential_condition
head(ANCOMdiff_condition2)

head(dat_ann_condition)
levels(ANCOMdiff_condition2$Family) <- c(levels(ANCOMdiff_condition2$Family), "below cutoff")
ANCOMdiff_condition2$Family[ANCOMdiff_condition2$W<(dat_ann_condition$y)]<-"below cutoff"
undercut_condition<-ANCOMdiff_condition2%>% filter(W<dat_ann_condition$y)
overcut_condition<-ANCOMdiff_condition2%>% filter(W>dat_ann_condition$y)
undercut_condition
length(undercut_condition)
length(overcut_condition)
ANCOMdiff_condition2
write.table(ANCOMdiff_condition2, "plots20211031rarefied/ANCOM_ww_condition.csv", sep = ",")

ANCOM_ww_condition <- read_delim("plots20211031rarefied/ANCOM_ww_condition.csv",
";", escape_double = FALSE, col_types = cols(W = col_number(),
x = col_number()), trim_ws = TRUE)
 
ANCOMdiff_condition <- as.data.frame(ANCOM_ww_condition)
str(ANCOMdiff_condition)
ANCOMdiff_condition<-ANCOMdiff_condition%>% filter(W!="Inf")
#Volcano plot
ancom_condition_plot<-ggplot(data = ANCOMdiff_condition, aes(y=W, x=x, col=Family))+scale_color_manual(values=ANCOMcolors)+facet_grid(cols =vars(group))+
  geom_point(size=rel(3))+geom_hline(yintercept = cut_off["detected_0.9"], linetype = "dashed") + 
  geom_text(data = dat_ann_condition, aes(x = x, y = y, label = label), 
            size = 5, vjust = -0.5, hjust = 0, color = "Black", parse = TRUE)+theme_linedraw()+
  geom_text_repel(data=ANCOMdiff_condition[which(ANCOMdiff_condition$detected_0.9==TRUE),],aes(y=W,x=x, label=Family), size=rel(3),nudge_y = 1,force=2,segment.size=0,parse=TRUE)+
  xlab(substitute(paste("Differential logs(clr)")))+
  ylab(substitute(paste(italic("W "),"statistic from ANCOM")))+themevolcano+theme(legend.position = "none")
ancom_condition_plot
ggsave("plots20211031rarefied/Ancom-Cond.png", dpi =600, width=18,height=8,units=c("in"))

ancom_condition_plot_genus<-ggplot(data = ANCOMdiff_condition, aes(y=W, x=x, col=Family))+scale_color_manual(values=ANCOMcolors)+facet_grid(cols =vars(group))+
  geom_point(size=rel(3))+geom_hline(yintercept = cut_off["detected_0.9"], linetype = "dashed") + 
  geom_text(data = dat_ann_condition, aes(x = x, y = y, label = label), 
            size = 5, vjust = -0.5, hjust = 0, color = "Black", parse = TRUE)+theme_linedraw()+
  geom_text_repel(data=ANCOMdiff_condition[which(ANCOMdiff_condition$detected_0.9==TRUE),],aes(y=W,x=x, label=Genus), size=rel(3),nudge_y = 1,force=2,segment.size=0,parse=TRUE)+
  xlab(substitute(paste("Differential logs(clr)")))+
  ylab(substitute(paste(italic("W "),"statistic from ANCOM")))+themevolcano+theme(legend.position = "none")
ancom_condition_plot_genus
ggsave("plots20211031rarefied/Ancom-Cond-genus.png", dpi =600, width=18,height=8,units=c("in"))






