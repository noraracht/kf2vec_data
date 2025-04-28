

require(ggplot2); require(scales); require(reshape2); 
#install.packages("dplyr")
require(dplyr)
#require(Hmisc)
library("readxl")
library(RColorBrewer)
library("ggsci")
#install.packages("ggrepel")
library("ggrepel")
library(ggpubr)
library(rstatix)



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd('/Users/admin/Documents/support')
getwd()


library(RColorBrewer)
my_palette = c(brewer.pal(9, "RdBu")[c(1,2, 3, 7, 9)])
my_palette_lst =as.list(strsplit(my_palette, " "))


my_set1 = c(brewer.pal(9, "Set1"))
my_set1_lst =as.list(strsplit(my_set1, " "))


#pl_per_quality=read.csv('/Users/nora/Documents/ml_metagenomics/paper_fig1/kmer_combined_with_quality.txt',sep=" ",header=TRUE)
#pl_per_quality=read.csv('/Users/nora/Documents/ml_metagenomics/depp_queries/kmer_combined_with_quality_v2.txt',sep=" ",header=TRUE)
pl_per_quality=read.csv('kmer_combined_with_quality_v2.txt',sep=" ",header=TRUE)
ncol(pl_per_quality)
nrow(pl_per_quality)

colnames(pl_per_quality)


ggplot(aes(x=cut(error,c(0,4,40),include.lowest =T,right = F), y=reference_representation_score),
       data=pl_per_quality)+
  geom_point(alpha = 0.9)+
  geom_boxplot()+
  theme_classic()+
  stat_compare_means(method = "t.test", label.x = 1.3, label.y = 1)
  #stat_compare_means(method = "anova")


ggplot(aes(x=error, y=n_contigs),
       data=pl_per_quality)+
  geom_point(alpha = 0.9)
  
  stat_summary()+
  stat_summary(geom="line",color="blue",aes(group="1"))+
  stat_summary(geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()+
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="Contig length (KB)", label = c("(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-10000]"))  +
  ylab("Placement error")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.7, hjust=0.5))
  
#label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-10000]" ) )  
#bquote('2560-'(10^4))
#bquote( '(2560-'  *10^4* ']' )

  #coord_cartesian(ylim=c(0,8))
  #data=pl_err)+
  #geom_plot(alpha=0.6)+
  #geom_point(alpha = 0.9, size = 0.5)+
  #stat_smooth(se=F)
  #geom_bar(position = "dodge",stat = "summary",fun = "mean")+
  #geom_boxplot(alpha=0.6)+
  stat_summary(aes(fill=
                     ifelse(cond =="true","Hypothetical",ifelse(cond == "combined","Default","Explore"))),geom="bar",color="black")+
  stat_summary()+
  scale_fill_brewer(name="")+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_abline(color="red")+
  #theme_classic()+
  #scale_fill_discrete(name="Condition") +
  theme_bw()+
  xlab("Condition")+
  ylab("Placement error")+
  scale_x_discrete(labels = c("Shorter\nkmer", "Extra\nlayer", "Global",
                              #  "Global", 
                              "Local", "Default\n(combined)", "True\nclade"),
                   name="")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.7, hjust=0.5),
        legend.position = c(.67,.9),legend.direction = "horizontal")
#theme(axis.line = element_line(colour = "black"),
#      panel.grid.major = element_blank(),
#      panel.grid.minor = element_blank(),
#      panel.border = element_blank(),
#      panel.background = element_blank()) 

#ggsave("wol_dev_queries_pl_perquery.pdf",width=5,height = 4)









#pl_per_pl_per_closest_dist=read.csv('/Users/nora/Documents/ml_metagenomics/paper_fig1/kmer_combined_with_closest_dist.txt',sep=" ",header=TRUE)
#pl_per_pl_per_closest_dist=read.csv('/Users/nora/Documents/ml_metagenomics/depp_queries/kmer_combined_with_closest_dist_v2.txt',sep=" ",header=TRUE)
pl_per_pl_per_closest_dist=read.csv('kmer_combined_with_closest_dist_v2.txt',sep=" ",header=TRUE)
head(pl_per_pl_per_closest_dist)

colnames(pl_per_pl_per_closest_dist)
nrow(pl_per_pl_per_closest_dist)

#pl_per_pl_per_closest_dist_depp=read.csv('/Users/nora/Documents/ml_metagenomics/depp_queries/pl_results_wol_extended_queries_depp_16s_v2.txt',sep=" ",header=TRUE)
pl_per_pl_per_closest_dist_depp=read.csv('pl_results_wol_extended_queries_depp_16s_v2.txt',sep=" ",header=TRUE)
head(pl_per_pl_per_closest_dist_depp)
nrow(pl_per_pl_per_closest_dist_depp)
unique(pl_per_pl_per_closest_dist_depp$cond)

subset_depp <-pl_per_pl_per_closest_dist_depp[pl_per_pl_per_closest_dist_depp$cond %in% c("depp"),]
subset_depp16s <-pl_per_pl_per_closest_dist_depp[pl_per_pl_per_closest_dist_depp$cond %in% c("depp_16s"),]
nrow(subset_depp)

part_depp = merge(x =subset_depp, y=pl_per_pl_per_closest_dist[,c("genome", "neighbor","dist")], by="genome")
part_depp16s = merge(x =subset_depp16s, y=pl_per_pl_per_closest_dist[,c("genome", "neighbor","dist")], by="genome")
part_kf2d = pl_per_pl_per_closest_dist[pl_per_pl_per_closest_dist$genome %in% as.list(part_depp$genome), ]
part_kf2d <- as.data.frame(append(part_kf2d, list(cond = 'kf2d'), after = 4))
head(part_kf2d)

nrow(part_kf2d)

pl_vs_dist_ALL_cond = rbind(part_depp, part_depp16s, part_kf2d)
head(pl_vs_dist_ALL_cond)

summary(aov(dist~pl_error,data=pl_per_pl_per_closest_dist))
cor.test(pl_per_pl_per_closest_dist$dist,pl_per_pl_per_closest_dist$pl_error,method="pearson")

mean(as.data.frame(pl_per_pl_per_closest_dist[pl_per_pl_per_closest_dist$pl_error <=1,])[,"dist"])/100

nrow(pl_per_pl_per_closest_dist[pl_per_pl_per_closest_dist$pl_error <=1,])/nrow(pl_per_pl_per_closest_dist)


nrow(pl_per_pl_per_closest_dist[pl_per_pl_per_closest_dist$pl_error >=7 & pl_per_pl_per_closest_dist$dist <=14.86926,])/nrow(pl_per_pl_per_closest_dist)*100



ggplot(aes(x=cut(pl_error,c(0,1,4,6,40),include.lowest =T,right = T), y=dist/100),
       data=pl_per_pl_per_closest_dist)+
  #geom_point(alpha = 0.9)+
  #geom_boxplot()+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  stat_summary(color="blue")+
  theme_classic()+
  #stat_compare_means(method = "t.test", label.x = 1.3, label.y = 160)
  #stat_compare_means(method = "anova")+
  xlab("Placement error")+
  ylab("Distance to the closest species")+
  annotate(x=3.0,y=1.70,label=bquote('p-value < 2.2*'*10^{"-16"}*', Pearson correlation 0.39'), geom="text")

#expression(paste("p-value = 1.213* ", 10^{-14}))
expression(paste("p-value = 2.2* ", 10^{-16}))
getwd()

#ggsave("wol_dev_queries_closest_dist.pdf",width=5,height = 4)
ggsave("wol_dev_queries_closest_dist_corr.pdf",width=4.8,height = 4)




cor.test(part_depp$dist,part_depp$pl_error,method="pearson")
cor.test(part_depp16s$dist,part_depp16s$pl_error,method="pearson")
cor.test(part_kf2d$dist,part_kf2d$pl_error,method="pearson")

pl_vs_dist_ALL_cond$cond <- factor(pl_vs_dist_ALL_cond$cond, levels = c('kf2d', 'depp', 'depp_16s'))
tail(pl_vs_dist_ALL_cond)

ggplot(aes(x=cut(pl_error,c(0,1,4,6,40),include.lowest =T,right = T), y=dist/100, color = cond),
       data=pl_vs_dist_ALL_cond)+
  #geom_point(alpha = 0.9)+
  #geom_boxplot()+
  #geom_pointrange(aes(ymin=dist/100-sd, ymax=dist/100+sd), position = position_dodge(width = 0.3))+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  stat_summary(position = position_dodge(width = 0.9))+
  theme_classic()+
  scale_color_brewer(name = "", palette = "Dark2")+
  #facet_wrap(.~cond, scales = 'free' )
  #stat_compare_means(method = "t.test", label.x = 1.3, label.y = 160)
  #stat_compare_means(method = "anova")+
  xlab("Placement error")+
  ylab("Distance to closest species")+
  annotate(x=3.0,y=1.70,label=bquote('p-value < 2.2*'*10^{"-16"}*', Pearson correlation 0.37'), geom="text")

ggsave("wol_dev_queries_closest_dist_all_cond.pdf",width=6.2,height = 4)
  

nrow(pl_vs_dist_ALL_cond)/128
unique(pl_vs_dist_ALL_cond$cond)

kf_subset = pl_vs_dist_ALL_cond[pl_vs_dist_ALL_cond$cond=='kf2d',]
nrow(kf_subset)
nrow(kf_subset[kf_subset$pl_error <2 & kf_subset$dist <=20,])/nrow(kf_subset)
nrow(kf_subset[kf_subset$pl_error >=7 & kf_subset$dist <=20,])/nrow(kf_subset)
kf_subset[kf_subset$pl_error <2 ,]
nrow(pl_per_pl_per_closest_dist[pl_per_pl_per_closest_dist$pl_error <2 & pl_per_pl_per_closest_dist$dist <=20,])/nrow(pl_per_pl_per_closest_dist)
nrow(pl_per_pl_per_closest_dist[pl_per_pl_per_closest_dist$pl_error >=7 & pl_per_pl_per_closest_dist$dist <=20,])/nrow(pl_per_pl_per_closest_dist)



ggplot(aes(y=pl_error, x=cut(dist/100,c(0,0.05,0.2,0.4,1.6)), color = cond),
       data=pl_vs_dist_ALL_cond)+
  #geom_point(alpha = 0.9)+
  #geom_boxplot()+
  #geom_pointrange(aes(ymin=dist/100-sd, ymax=dist/100+sd), position = position_dodge(width = 0.3))+
  #geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  stat_summary(position = position_dodge(width = 0.1))+
  theme_classic()+
  scale_color_brewer(name = "", palette = "Dark2", labels = c("kf2vec", "DEPP_381", "DEPP_16S"))+
  #facet_wrap(.~cond, scales = 'free' )
  #stat_compare_means(method = "t.test", label.x = 1.3, label.y = 160)
  #stat_compare_means(method = "anova")+
  ylab("Placement error")+
  xlab("Distance to the closest species")+
  theme(legend.position = c(0.15,0.9), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )
ggsave("wol_dev_queries_closest_dist_all_pointrange.pdf",width=4.8,height = 4)

getwd()
  #annotate(x=3.0,y=1.70,label=bquote('p-value < 2.2*'*10^{"-16"}*', Pearson correlation 0.37'), geom="text")

  #unique(pl_vs_dist_ALL_cond$cond)
  
  
  pl_vs_dist_ALL_cond$cond <- factor(pl_vs_dist_ALL_cond$cond, levels = c('kf2d', 'depp', 'depp_16s'))
  
 
 ggplot(aes(x=pl_error, shape=cut(dist/100,c(0,0.05,0.2,0.4,1.6)), color = cond),
         data=pl_vs_dist_ALL_cond)+
    #geom_point(alpha = 0.9)+
    #geom_boxplot()+
    #geom_pointrange(aes(ymin=dist/100-sd, ymax=dist/100+sd), position = position_dodge(width = 0.3))+
    #geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
    stat_ecdf()+
    theme_classic()+
    scale_color_brewer(name = "", palette = "Dark2", labels = c("kf2vec", "DEPP_381", "DEPP_16S"))+
    facet_wrap(.~cut(dist/100,c(0,0.05,0.2,0.4,1.6)) ,nrow=2)+
    #stat_compare_means(method = "t.test", label.x = 1.3, label.y = 160)
    #stat_compare_means(method = "anova")+
    xlab("Placement error")+
    #scale_x_discrete(name="Placement error", label = c("kf2d", "DEPP_381", "DEPP_16S"))  +
    ylab("ECDF")+
    theme(legend.position = c(0.9,0.18), legend.margin=margin(0,0,0,0),)
          #axis.text.x = element_text(size = 8)
    
    #theme(legend.position = "none")

  ggsave("wol_dev_queries_closest_dist_ALL_COND.pdf",width=6.2,height = 4)

  getwd()
#expression(paste("p-value = 1.213* ", 10^{-14}))
expression(paste("p-value = 2.2* ", 10^{-16}))
getwd()

#ggsave("wol_dev_queries_closest_dist.pdf",width=5,height = 4)
#ggsave("wol_dev_queries_closest_dist_ALL_COND.pdf",width=4.8,height = 4)







bx1
g1 = pl_per_pl_per_closest_dist[pl_per_pl_per_closest_dist$pl_error < 6,]
g2 = pl_per_pl_per_closest_dist[pl_per_pl_per_closest_dist$pl_error >= 6,]

g2$dist
my_test <- t.test(g1$dist, g2$dist)
my_test$p.value

stat.test <- t_test(g1$dist ~ g2$dist) %>%
  add_significance()

#stat.test %>% add_xy_position(x = "Ortamlar")
#stat.test <- stat.test %>% mutate(y.position=5.15)

bx1+stat_pvalue_manual(my_test$p.value, label = "p = {p}")



nrow(df)





##########################################################################
library(dplyr)

head(pl_per_quality)
colnames(pl_per_pl_per_closest_dist)

nrow(pl_per_pl_per_closest_dist)

colnames(pl_per_quality)
ncol(pl_per_quality)

pl_per_quality_num <- pl_per_quality %>%
  type.convert(as.is = TRUE) %>%
  select(where(is.numeric))

ncol(pl_per_quality_num)

g1 = pl_per_quality_num[pl_per_quality_num$error < 5,]
g2 = pl_per_quality_num[pl_per_quality_num$error >= 5,]

nrow(g1)
nrow(g2)

colnames(g1)




#provide column names
df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df) <- c('var', 'p_val')
count = 0

for(i in 4:ncol(g1)) {       # for-loop over columns
  
  
  data1 <- g1[ , i]
  data2 <- g2[ , i]
  
  
  
  #if (unique(data1)!= unique(data2)){
  if  (! setequal(data1, data2)){
    my_test <- t.test(data1, data2)}
  
    if (my_test$p.value < 1)
      {
      print(colnames(g1[i]))
      print(my_test$p.value)
      count=count +1
      
      df[nrow(df) + 1,] = c(colnames(g1[i]), my_test$p.value)
      }
}
df$p_val =p.adjust(as.numeric(df$p_val))


sum(p.adjust(df$p_val)<0.05)

df['var'][df['var'] == 'n_genes_called'] <- 'Number of genes called'
df['var'][df['var'] == 'n_genes_mapped'] <- 'Number of genes mapped'
df['var'][df['var'] == 'n_contigs'] <- 'Number of contigs with mapped genes'
df['var'][df['var'] == 'proportion_genes_retained_in_major_clades'] <- 'Proportion of genes retained in major clades'
df['var'][df['var'] == 'genes_retained_index'] <- 'Genes retained index'
df['var'][df['var'] == 'clade_separation_score'] <- 'Clade separation score'
df['var'][df['var'] == 'contamination_portion'] <- 'Contamination portion'
df['var'][df['var'] == 'n_effective_surplus_clades'] <- 'Number of effective surplus clades'
df['var'][df['var'] == 'mean_hit_identity'] <- 'Mean hit identity'
df['var'][df['var'] == 'reference_representation_score'] <- 'Reference representation score'


df['var'][df['var'] == 'X..genomes'] <- 'Count of genomes'
df['var'][df['var'] == 'X..markers'] <- 'Count of markers'
df['var'][df['var'] == 'X..marker.sets'] <- 'Count of marker sets'
df['var'][df['var'] == 'Strain.heterogeneity'] <- 'Strain heterogeneity'
df['var'][df['var'] == 'Genome.size..bp.'] <- 'Genome size'
df['var'][df['var'] == 'X..ambiguous.bases'] <- 'Count of ambiguous bases'
df['var'][df['var'] == 'X..scaffolds'] <- 'Count of scaffolds'
df['var'][df['var'] == 'X..contigs_x'] <- 'Count of contigs'
df['var'][df['var'] == 'N50..scaffolds.'] <- 'N50 (scaffolds)'
df['var'][df['var'] == 'N50..contigs.'] <- 'N50 (contigs)'
df['var'][df['var'] == 'Mean.scaffold.length..bp.'] <- 'Mean scaffold length'
df['var'][df['var'] == 'Mean.contig.length..bp.'] <- 'Mean contig length'
df['var'][df['var'] == 'Longest.scaffold..bp.'] <- 'Longest scaffold'
df['var'][df['var'] == 'Longest.contig..bp.'] <- 'Longest contig'
df['var'][df['var'] == 'GC.std..scaffolds...1kbp.'] <- 'GC for scaffolds above 1kbp'
df['var'][df['var'] == 'Coding.density'] <- 'Coding density'
df['var'][df['var'] == 'Translation.table'] <- 'Translation table'
df['var'][df['var'] == 'X..predicted.genes'] <- 'Count of predicted genes'
df['var'][df['var'] == 'X0'] <- 'qa0'
df['var'][df['var'] == 'X1'] <- 'qa1'
df['var'][df['var'] == 'X2'] <- 'qa2'
df['var'][df['var'] == 'X3'] <- 'qa3'
df['var'][df['var'] == 'X4'] <- 'qa4'
df['var'][df['var'] == 'X5.'] <- 'qa5+'


df['var'][df['var'] == 'X..contigs.....0.bp.'] <- 'Number of contigs of 0 bp or longer'
df['var'][df['var'] == 'X..contigs.....1000.bp.'] <- 'Number of contigs of 1000 bp or longer'
df['var'][df['var'] == 'X..contigs.....5000.bp.'] <- 'Number of contigs of 5000 bp or longer'
df['var'][df['var'] == 'X..contigs.....10000.bp.'] <- 'Number of contigs of 10000 bp or longer'
df['var'][df['var'] == 'X..contigs.....25000.bp.'] <- 'Number of contigs of 25000 bp or longer'
df['var'][df['var'] == 'X..contigs.....50000.bp.'] <- 'Number of contigs of 50000 bp or longer'

df['var'][df['var'] == 'Total.length.....0.bp.'] <- 'Total length at 0 bp or above'
df['var'][df['var'] == 'Total.length.....1000.bp.'] <- 'Total length at 1000 bp or above'
df['var'][df['var'] == 'Total.length.....5000.bp.'] <- 'Total length at 5000 bp or above'
df['var'][df['var'] == 'Total.length.....10000.bp.'] <- 'Total length at 10000 bp or above'
df['var'][df['var'] == 'Total.length.....25000.bp.'] <- 'Total length at 25000 bp or above'
df['var'][df['var'] == 'Total.length.....50000.bp.'] <- 'Total length at 50000 bp or above'

df['var'][df['var'] == 'X..contigs_y'] <- 'Count of contigs in assembly'
df['var'][df['var'] == 'Largest.contig'] <- 'Largest contig in assembly'
df['var'][df['var'] == 'Total.length'] <- 'Total length of assembly'
df['var'][df['var'] == 'GC....'] <- 'G and C in assembly (%)'
df['var'][df['var'] == 'X..N.s.per.100.kbp'] <- 'Count of N per 100kbp'
df['var'][df['var'] == 'dist'] <- 'Distance to closest species'

nrow(df[df$p_val<=0.05,])

ggplot(aes(y=reorder(var, p_val), x=p_val,shape=p_val<0.05),
       data=df[df$p_val<1,])+
  geom_point(alpha = 0.9)+
  scale_shape_manual(guide="none",values=c(1,16))+
  geom_vline(xintercept = 0.05,color="red",linetype=3)+
  #geom_boxplot()
  theme_classic()+
  #theme(axis.text.x = element_text(angle = 0, vjust = 0.8, hjust=0.5))
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8))+
  xlab("p value")
        #legend.position = c(.82,.88),legend.direction = "vertical", legend.margin=margin(t = -0.5, unit='cm'))
ggsave("quality_per_genome_main.pdf",width=5,height = 4)


ggplot(aes(y=reorder(var, p_val), x=p_val,shape=p_val<0.05),
       data=df[df$p_val<2,])+
  geom_point(alpha = 0.9)+
  scale_shape_manual(guide="none",values=c(1,16))+
  geom_vline(xintercept = 0.05,color="red",linetype=3)+
  #geom_boxplot()
  theme_classic()+
  #theme(axis.text.x = element_text(angle = 0, vjust = 0.8, hjust=0.5))
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 7))+
  xlab("p value")
#legend.position = c(.82,.88),legend.direction = "vertical", legend.margin=margin(t = -0.5, unit='cm'))
ggsave("quality_per_genome_main_sup.pdf",width=5,height = 6)

getwd()

count

df$var

##########################################################################
# FROM HERE

#EXTENDED QUERIES MULTIPLE GENES

#pl_err2=read.csv('/Users/nora/Documents/ml_metagenomics/paper_fig1/pl_results_wol_extended_queries_depp_multgen.txt',sep=" ",header=TRUE)
#pl_err2=read.csv('/Users/nora/Documents/ml_metagenomics/depp_queries/pl_results_wol_extended_queries_depp_multgen_v2.txt',sep=" ",header=TRUE)
pl_err2=read.csv('pl_results_wol_extended_queries_depp_multgen_v2.txt',sep=" ",header=TRUE)

head(pl_err2)
unique(pl_err2$cond)

pl_err2$cond = factor(pl_err2$cond, 
                      levels=c("Chunked_Claded", "Unchunked_Claded", "Chunked_Uncladed", "Unchunked_Uncladed", "depp","depp_380","depp_256","depp_128","depp_64","depp_32","depp_16","depp_4", "depp_1"  ))

cond_to_plot = "Unchunked_Uncladed"

#pl_err2$cond <- factor(pl_err2$cond, levels=c("6kC", "moreLayers", "main", "global", "local", "combined", "true" ))

ggplot(aes(x=cond, y=pl_error),
       data=pl_err2)+
  #data=pl_err)+
  #geom_plot(alpha=0.6)+
  #geom_bar(position = "dodge",stat = "summary",fun = "mean")+
  #geom_boxplot(alpha=0.6)+
  stat_summary()+
  stat_summary(geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_abline(color="red")+
  #theme_classic()+
  #scale_fill_discrete(name="Condition") +
  theme_bw()+
  xlab("Condition")+
  ylab("Placement error")+
  #scale_x_discrete(labels = c("Shorter kmer", "Extra layer", "Main", "Global", "Local", "Combined", "True"))+
  theme(axis.text.x = element_text(angle = 30, vjust = 1.1, hjust=1.2))
#theme(axis.line = element_line(colour = "black"),
#      panel.grid.major = element_blank(),
#      panel.grid.minor = element_blank(),
#      panel.border = element_blank(),
#      panel.background = element_blank()) 

ggsave("dev_queries_multiple_conditions_multiple_genes.pdf",width=5,height = 4)


ggplot(aes(x=as.numeric(as.character(sub(".*_","",cond))), y=pl_error),
       data=pl_err2[grepl("depp_",pl_err2$cond), ])+
  #data=pl_err)+
  #geom_plot(alpha=0.6)+
  #geom_bar(position = "dodge",stat = "summary",fun = "mean")+
  #geom_boxplot(alpha=0.6)+
  stat_summary()+
  stat_summary(geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))
#geom_point()+
#geom_errorbar(stat="summary")+
#geom_line() +
#stat_summary(fun.y = "median", geom = "point", size = 3)+
#geom_boxplot(outlier.shape = NA) +
#geom_abline(color="red")+
#theme_classic()+
#scale_fill_discrete(name="Condition") +
theme_bw()+
  xlab("Condition")+
  ylab("Placement error")+
  #scale_x_discrete(labels = c("Shorter kmer", "Extra layer", "Main", "Global", "Local", "Combined", "True"))+
  theme(axis.text.x = element_text(angle = 30, vjust = 1.1, hjust=1.2))


unique(pl_err2$cond)

ggplot(aes(color=cond, x=pl_error),
       data=pl_err2[pl_err2$cond %in% c("depp","depp_1",cond_to_plot), ])+
  stat_ecdf(alpha=0.5,size=1)+
  theme_classic()
#data=pl_err)+
#geom_plot(alpha=0.6)+
#geom_bar(position = "dodge",stat = "summary",fun = "mean")+
#geom_boxplot(alpha=0.6)+


ggplot(aes(fill=cond, x=cut(pl_error,c(0,1,4,6,11,40),include.lowest =TRUE,right = FALSE)),
       data=pl_err2[pl_err2$cond %in% c("depp","depp_1",cond_to_plot), ])+
  geom_bar(position = position_dodge())+
  
  #geom_text(aes(x=2.5+as.numeric(cond)/2,y=100,label=round(`.`,2),color=cond),
  #          data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c("depp","depp_1",cond_to_plot), ],fun.aggregate = mean,value.var="pl_error"))+
  #
  annotate(x=3.5,y=120,label="Mean error:",geom="text")+
  geom_text(aes(x=3.0, y=100,label=round(`.`,2),color=cond),
            data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c(cond_to_plot ), ],fun.aggregate = mean,value.var="pl_error"))+
  
  geom_text(aes(x=3.5, y=100,label=round(`.`,2),color=cond),
            data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c("depp"), ],fun.aggregate = mean,value.var="pl_error"))+
  
  geom_text(aes(x=4.0, y=100,label=round(`.`,2),color=cond),
            data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c("depp_1") ,],fun.aggregate = mean,value.var="pl_error"))+
  
  
  theme_classic()+
  #scale_fill_manual(name="", values = c("#252525", "#ca0020", "#0571b0" ), 
  #                   labels = c("Simulated", "Corrected", "Uncorrected"))
  #scale_fill_manual(name="", values = c("#2c7bb6", "#fdae61", "#d7191c" ))
  scale_fill_brewer(name = "",palette="Dark2", labels = c("kf2vec", "DEPP_381", "DEPP_1"))+
  scale_color_brewer(name = "",palette="Dark2", labels = c("kf2vec", "DEPP_381", "DEPP_1"))+
  ylab("Count")+
  xlab("Placement error")+
  #theme(axis.title.x = element_blank())+
  theme_classic()+
  theme(legend.position = c(.9,.8), legend.margin=margin(t = -0.5, unit='cm') )

ggsave("kf2d_vs_depp.pdf",width=4.8,height = 4)

getwd()

#geom_point()+
#geom_errorbar(stat="summary")+
#geom_line() +
#stat_summary(fun.y = "median", geom = "point", size = 3)+
#geom_boxplot(outlier.shape = NA) +
#geom_abline(color="red")+
#theme_classic()+
#scale_fill_discrete(name="Condition") +
theme_bw()+
  xlab("Condition")+
  ylab("Placement error")+
  #scale_x_discrete(labels = c("Shorter kmer", "Extra layer", "Main", "Global", "Local", "Combined", "True"))+
  theme(axis.text.x = element_text(angle = 30, vjust = 1.1, hjust=1.2))




##########################################################################
# BEST
# Deft vs Depp 16s queries



#pl_err2=read.csv('/Users/nora/Documents/ml_metagenomics/depp_queries/pl_results_wol_extended_queries_depp_16s_v2.txt',sep=" ",header=TRUE)
pl_err2=read.csv('pl_results_wol_extended_queries_depp_16s_v2.txt',sep=" ",header=TRUE)
pl_err2

cond_to_plot = "Unchunked_Uncladed"


#pl_err2$cond <- factor(pl_err2$cond, levels=c("6kC", "moreLayers", "main", "global", "local", "combined", "true" ))

unique(pl_err2$cond)
#pl_err2$cond <- factor(pl_err2$cond, levels=c("Unchunked_Uncladed", "depp", "depp_16s", "Chunked_Uncladed", "Unchunked_Claded", "Chunked_Claded" ,  ))

pl_err2$cond = factor(pl_err2$cond, 
                      levels=c("Chunked_Claded", "Unchunked_Claded", "Chunked_Uncladed", "Unchunked_Uncladed", "depp","depp_16s"))



ggplot(aes(x=cond, y=pl_error),
       data=pl_err2)+
  #data=pl_err)+
  #geom_plot(alpha=0.6)+
  #geom_bar(position = "dodge",stat = "summary",fun = "mean")+
  #geom_boxplot(alpha=0.6)+
  stat_summary()+
  stat_summary(geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_abline(color="red")+
  #theme_classic()+
  #scale_fill_discrete(name="Condition") +
  theme_bw()+
  xlab("Condition")+
  ylab("Placement error")+
  #scale_x_discrete(labels = c("Shorter kmer", "Extra layer", "Main", "Global", "Local", "Combined", "True"))+
  theme(axis.text.x = element_text(angle = 30, vjust = 1.1, hjust=1.2))

ggsave("dev_queries_multiple_conditions_16s_genes.pdf",width=5,height = 4)


unique(pl_err2$cond)

#nrow(pl_err2[pl_err2$cond %in% c("combined"), ])
#nrow(pl_err2[pl_err2$cond %in% c("combined") & pl_err2$pl_error >4, ])
#nrow(pl_err2[pl_err2$cond %in% c("combined") & pl_err2$pl_error ==0, ])
nrow(pl_err2[pl_err2$cond %in% c("depp") & pl_err2$pl_error ==0, ])
nrow(pl_err2[pl_err2$cond %in% c("depp_16s") & pl_err2$pl_error ==27, ])

#head(pl_err2[pl_err2$cond %in% c("depp","depp_16s","Unchunked_Uncladed"), ])


unique(pl_err2$cond)

w = pl_err2[pl_err2$cond %in% c(cond_to_plot, "depp","depp_16s"), ]
tail(w)
w1 = w[w$cond =='Unchunked_Uncladed',]
nrow(w1[w1$pl_error>4,])/nrow(w1)
tail(w)
w1[order(w1$pl_error, decreasing = TRUE),]

ggplot(aes(fill=cond, x=cut(pl_error,c(0,1,4, 6,11,40),include.lowest =TRUE, right = FALSE)),
       data=pl_err2[pl_err2$cond %in% c(cond_to_plot, "depp","depp_16s"), ])+
  geom_bar(position = position_dodge())+
  annotate(x=4.0,y=40,label="Mean error:",geom="text")+
  #geom_text(aes(x=2.5+as.numeric(cond)/4, y=33,label=round(`.`,2),color=cond),
  geom_text(aes(x=3.2, y=33,label=round(`.`,2),color=cond),
            data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c(cond_to_plot ), ],fun.aggregate = mean,value.var="pl_error"))+
  
  geom_text(aes(x=4.0, y=33,label=round(`.`,2),color=cond),
            data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c("depp"), ],fun.aggregate = mean,value.var="pl_error"))+
  
  geom_text(aes(x=4.7, y=33,label=round(`.`,2),color=cond),
            data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c("depp_16s") ,],fun.aggregate = mean,value.var="pl_error"))+
  
  theme_classic()+
  #scale_fill_manual(name="", values = c("#252525", "#ca0020", "#0571b0" ), 
  #                   labels = c("Simulated", "Corrected", "Uncorrected"))
  #scale_fill_manual(name="", values = c("#2c7bb6", "#fdae61", "#d7191c" ), labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  #scale_color_manual(name="", values = c("#2c7bb6", "#fdae61", "#d7191c" ), labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  scale_fill_brewer(name = "",palette="Dark2", labels = c("kf2vec", "DEPP_381", "DEPP_16S"))+
  scale_color_brewer(name = "",palette="Dark2", labels = c("kf2vec", "DEPP_381", "DEPP_16S"))+
  ylab("Count")+
  xlab("Placement error")+
  #theme(axis.title.x = element_blank())+
  theme_classic()+
  #geom_bar(position = "dodge",stat = "summary",fun = "mean", colour="black", alpha=0.9)+
  theme(legend.position = c(.85,.8), legend.margin=margin(t = -0.5, unit='cm') )

ggsave("kf2d_vs_depp_16s.pdf",width=3.8,height = 4)

getwd()

ggplot(aes(color=cond, x=pl_error),
       data=pl_err2[pl_err2$cond %in% c(cond_to_plot, "depp","depp_16s"), ])+
  stat_ecdf()+
  annotate(x=3.5,y=.40,label="Mean error:",geom="text")+
  #geom_text(aes(x=2.5+as.numeric(cond)/4, y=33,label=round(`.`,2),color=cond),
  geom_text(aes(x=3.0, y=.33,label=round(`.`,2),color=cond),
            data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c(cond_to_plot ), ],fun.aggregate = mean,value.var="pl_error"))+
  
  geom_text(aes(x=3.5, y=.33,label=round(`.`,2),color=cond),
            data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c("depp"), ],fun.aggregate = mean,value.var="pl_error"))+
  
  geom_text(aes(x=4.0, y=.33,label=round(`.`,2),color=cond),
            data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c("depp_16s") ,],fun.aggregate = mean,value.var="pl_error"))+
  
  theme_classic()+
  #scale_fill_manual(name="", values = c("#252525", "#ca0020", "#0571b0" ), 
  #                   labels = c("Simulated", "Corrected", "Uncorrected"))
  #scale_fill_manual(name="", values = c("#2c7bb6", "#fdae61", "#d7191c" ), labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  #scale_color_manual(name="", values = c("#2c7bb6", "#fdae61", "#d7191c" ), labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  scale_fill_brewer(name = "",palette="Dark2", labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  scale_color_brewer(name = "",palette="Dark2", labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  ylab("Count")+
  xlab("Placement error")+
  #theme(axis.title.x = element_blank())+
  theme_classic()+
  #geom_bar(position = "dodge",stat = "summary",fun = "mean", colour="black", alpha=0.9)+
  theme(legend.position = c(.9,.8), legend.margin=margin(t = -0.5, unit='cm') )

getwd()

# For ppt  
#& pl_err2$pl_error <=1
ggplot(aes(fill=cond, x=cut(pl_error,c(0, 1, 40),include.lowest =TRUE,right = FALSE)),
       data=pl_err2[pl_err2$cond %in% c("depp","depp_16s",cond_to_plot), ])+
  geom_bar(position = position_dodge())+
  annotate(x=1.5,y=90,label="Mean error:",geom="text")+
  geom_text(aes(x=1.0+as.numeric(cond)/4,y=84,label=round(`.`,2),color=cond),
            data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c("depp","depp_16s","combined"), ],fun.aggregate = mean,value.var="pl_error"))+
  theme_classic()+
  #scale_fill_manual(name="", values = c("#252525", "#ca0020", "#0571b0" ), 
  #                   labels = c("Simulated", "Corrected", "Uncorrected"))
  #scale_fill_manual(name="", values = c("#2c7bb6", "#fdae61", "#d7191c" ), labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  #scale_color_manual(name="", values = c("#2c7bb6", "#fdae61", "#d7191c" ), labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  scale_fill_brewer(name = "",palette="Dark2", labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  scale_color_brewer(name = "",palette="Dark2", labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  ylab("Count")+
  xlab("Placement error")+
  #theme(axis.title.x = element_blank())+
  theme_classic()+
  theme(legend.position = c(.15,.85), legend.margin=margin(t = -0.5, unit='cm') )


ggsave("kf2d_vs_depp_16s_ppt_v1.pdf",width=5,height = 4)

pl_err2[pl_err2$cond %in% c("depp","depp_16s",cond_to_plot) & (pl_err2$pl_error <=1), ]

ggplot(aes(fill=cond, x=cut(pl_error,c(0, 1, 40),include.lowest =TRUE,right = FALSE)),
       data=pl_err2[pl_err2$cond %in% c("depp","depp_16s",cond_to_plot) & (pl_err2$pl_error <1), ])+
  geom_bar(position = position_dodge())+
  annotate(x=1.0,y=90,label="Mean error:",geom="text")+
  geom_text(aes(x=0.3+as.numeric(cond)/3,y=84,label=round(`.`,2),color=cond),
            data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c("depp","depp_16s",cond_to_plot), ],fun.aggregate = mean,value.var="pl_error"))+
  theme_classic()+
  #scale_fill_manual(name="", values = c("#252525", "#ca0020", "#0571b0" ), 
  #                   labels = c("Simulated", "Corrected", "Uncorrected"))
  #scale_fill_manual(name="", values = c("#2c7bb6", "#fdae61", "#d7191c" ), labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  #scale_color_manual(name="", values = c("#2c7bb6", "#fdae61", "#d7191c" ), labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  scale_fill_brewer(name = "",palette="Dark2", labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  scale_color_brewer(name = "",palette="Dark2", labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  ylab("Count")+
  xlab("Placement error")+
  #theme(axis.title.x = element_blank())+
  theme_classic()+
  theme(legend.position = "none", legend.margin=margin(t = -0.5, unit='cm') )+
  annotate(x=0.7,y=30,label=c("kf2d"), geom="text")+
  annotate(x=1.0,y=30,label=c("DEPP_381"), geom="text")+
  annotate(x=1.3,y=30,label=c("DEPP_16s"), geom="text")

ggsave("kf2d_vs_depp_16s_ppt_v2.pdf",width=5,height = 4)

# Get the means from the un-normed data



ggplot(aes(fill=cond, x = cond, y=pl_error ),
       data=pl_err2[pl_err2$cond %in% c("depp","depp_16s",cond_to_plot), ])+
  geom_bar(stat = "summary", fun = "mean", color = "black", position=position_dodge(), alpha = 0.9)+
  #geom_errorbar(stat='summary', width=.2)+
  geom_pointrange(stat='summary', alpha = 0.9)+
  #stat_summary(geom="range",fun.data = function(x) median_hilow_(x,ci=0.8))+
  #annotate(x=1.0,y=90,label="Mean error:",geom="text")+
  #geom_text(aes(x=0.3+as.numeric(cond)/3,y=84,label=round(`.`,2),color=cond),
  #          data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c("depp","depp_16s","combined"), ],fun.aggregate = mean,value.var="pl_error"))+
  theme_classic()+
  scale_x_discrete(labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  #scale_fill_manual(name="", values = c("#252525", "#ca0020", "#0571b0" ), 
  #                   labels = c("Simulated", "Corrected", "Uncorrected"))
  #scale_fill_manual(name="", values = c("#2c7bb6", "#fdae61", "#d7191c" ), labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  #scale_color_manual(name="", values = c("#2c7bb6", "#fdae61", "#d7191c" ), labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  scale_fill_brewer(name = "",palette="Dark2", labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  #scale_color_brewer(name = "",palette="Dark2", labels = c("kf2d", "DEPP_381", "DEPP_16S"))+
  ylab("Mean placement error")+
  #xlab("")
  #theme(axis.title.x = element_blank())+
  theme_classic()+
  theme(legend.position = "none", axis.title.x = element_blank(),
        legend.margin=margin(t = -0.5, unit='cm') )


ggsave("kf2d_vs_depp_16s_ppt_v3.pdf",width=4,height = 4)

########################################################################## 
# DEPP RANDOM

# Build plots for variabel number of genes with random selection
#pl_err2=read.csv('/Users/nora/Documents/ml_metagenomics/paper_fig1/pl_results_wol_extended_queries_depp_multgen_random.txt',sep=" ",header=TRUE)
#pl_err2=read.csv('/Users/nora/Documents/ml_metagenomics/depp_queries/pl_results_wol_extended_queries_depp_multgen_random_v2.txt',sep=" ",header=TRUE)
pl_err2=read.csv('pl_results_wol_extended_queries_depp_multgen_random_v2.txt',sep=" ",header=TRUE)

head(pl_err2)
nrow(pl_err2)
unique(pl_err2$cond)

pl_err2$cond = factor(pl_err2$cond, 
                      levels=c("Chunked_Claded", "Unchunked_Claded", "Chunked_Uncladed", "Unchunked_Uncladed", "depp","depp_380","depp_256","depp_128","depp_64","depp_32","depp_16","depp_4", "depp_1"  ))

cond_to_plot = "Unchunked_Uncladed"


#pl_err2$cond <- factor(pl_err2$cond, levels=c("6kC", "moreLayers", "main", "global", "local", "combined", "true" ))

ggplot(aes(x=cond, y=pl_error),
       data=pl_err2)+
  #data=pl_err)+
  #geom_plot(alpha=0.6)+
  #geom_bar(position = "dodge",stat = "summary",fun = "mean")+
  #geom_boxplot(alpha=0.6)+
  stat_summary()+
  stat_summary(geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_abline(color="red")+
  #theme_classic()+
  #scale_fill_discrete(name="Condition") +
  theme_bw()+
  xlab("Condition")+
  ylab("Placement error")+
  #scale_x_discrete(labels = c("Shorter kmer", "Extra layer", "Main", "Global", "Local", "Combined", "True"))+
  theme(axis.text.x = element_text(angle = 30, vjust = 1.1, hjust=1.2))
#theme(axis.line = element_line(colour = "black"),
#      panel.grid.major = element_blank(),
#      panel.grid.minor = element_blank(),
#      panel.border = element_blank(),
#      panel.background = element_blank()) 

ggsave("dev_queries_multiple_conditions_multiple_genes_random.pdf",width=5,height = 4)


ggplot(aes(x=as.numeric(as.character(sub(".*_","",cond))), y=pl_error),
       data=pl_err2[grepl("depp_",pl_err2$cond), ])+
  #data=pl_err)+
  #geom_plot(alpha=0.6)+
  #geom_bar(position = "dodge",stat = "summary",fun = "mean")+
  #geom_boxplot(alpha=0.6)+
  stat_summary()+
  stat_summary(geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))
#geom_point()+
#geom_errorbar(stat="summary")+
#geom_line() +
#stat_summary(fun.y = "median", geom = "point", size = 3)+
#geom_boxplot(outlier.shape = NA) +
#geom_abline(color="red")+
#theme_classic()+
#scale_fill_discrete(name="Condition") +
theme_bw()+
  xlab("Condition")+
  ylab("Placement error")+
  #scale_x_discrete(labels = c("Shorter kmer", "Extra layer", "Main", "Global", "Local", "Combined", "True"))+
  theme(axis.text.x = element_text(angle = 30, vjust = 1.1, hjust=1.2))


ggplot(aes(color=cond, x=pl_error),
       data=pl_err2[pl_err2$cond %in% c("depp","depp_1",cond_to_plot), ])+
  stat_ecdf(alpha=0.5,size=1)+
  theme_classic()
#data=pl_err)+
#geom_plot(alpha=0.6)+
#geom_bar(position = "dodge",stat = "summary",fun = "mean")+
#geom_boxplot(alpha=0.6)+

unique(pl_err2$cond)
ggplot(aes(fill=cond, x=cut(pl_error,c(0,1,4,6,11,40),include.lowest =TRUE,right = FALSE)),
       data=pl_err2[pl_err2$cond %in% c("depp","depp_1",cond_to_plot), ])+
  geom_bar(position = position_dodge())+
  
  #geom_text(aes(x=2.5+as.numeric(cond)/4,y=100,label=round(`.`,2),color=cond),
  #          data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c("depp","depp_1",cond_to_plot), ],fun.aggregate = mean,value.var="pl_error"))+
  annotate(x=3.5,y=120,label="Mean error:",geom="text")+
  geom_text(aes(x=3.0, y=100,label=round(`.`,2),color=cond),
            data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c(cond_to_plot ), ],fun.aggregate = mean,value.var="pl_error"))+
  
  geom_text(aes(x=3.5, y=100,label=round(`.`,2),color=cond),
            data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c("depp"), ],fun.aggregate = mean,value.var="pl_error"))+
  
  geom_text(aes(x=4.0, y=100,label=round(`.`,2),color=cond),
            data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c("depp_1") ,],fun.aggregate = mean,value.var="pl_error"))+
  
  theme_classic()+
  #scale_fill_manual(name="", values = c("#252525", "#ca0020", "#0571b0" ), 
  #                   labels = c("Simulated", "Corrected", "Uncorrected"))
  #scale_fill_manual(name="", values = c("#2c7bb6", "#fdae61", "#d7191c" ))
  scale_fill_brewer(name = "",palette="Dark2", labels = c("kf2vec", "DEPP_381", "DEPP_1"))+
  scale_color_brewer(name = "",palette="Dark2", labels = c("kf2vec", "DEPP_381", "DEPP_1"))+
  ylab("Count")+
  xlab("Placement error")+
  #theme(axis.title.x = element_blank())+
  theme_classic()+
  theme(legend.position = c(.9,.8), legend.margin=margin(t = -0.5, unit='cm') )


ggsave("kf2d_vs_depp_random.pdf",width=4.8,height = 4)

getwd()

head(pl_err2[pl_err2$cond %in% c("depp","depp_1",cond_to_plot),])
ggplot(aes(fill=cond, x=cut(pl_error,c(0,1,4,6,11,40),include.lowest =T,right = F)),
       data=pl_err2[pl_err2$cond %in% c("depp","depp_1",cond_to_plot), ])+
  geom_bar(position = position_dodge())+
  theme_classic()+
  #scale_fill_manual(name="", values = c("#252525", "#ca0020", "#0571b0" ), 
  #                   labels = c("Simulated", "Corrected", "Uncorrected"))
  #scale_fill_manual(name="", values = c("#2c7bb6", "#fdae61", "#d7191c" ))
  scale_fill_brewer(name = "",palette="Dark2", labels = c("kf2d", "DEPP_381", "DEPP_1"))+
  ylab("Count")+
  xlab("Placement error")+
  #theme(axis.title.x = element_blank())+
  theme_classic()+
  theme(legend.position = c(.9,.8), legend.margin=margin(t = -0.5, unit='cm') )






##########################################################################  
##########################################################################  
##########################################################################  
##########################################################################  
##########################################################################
