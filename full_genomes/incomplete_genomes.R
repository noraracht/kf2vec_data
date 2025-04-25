

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


library(stringr)



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd('/Users/admin/Documents/support')
getwd()


library(RColorBrewer)
my_palette = c(brewer.pal(9, "RdBu")[c(1,2, 3, 7, 9)])
my_palette_lst =as.list(strsplit(my_palette, " "))


dark2_colors <-brewer.pal(n = 8, name = "Dark2")
dark2_colors

set1_colors <-brewer.pal(n = 8, name = "Set1")
set1_colors

accent_colors <-brewer.pal(n = 8, name = "Accent")
accent_colors

#pl_per_contig=read.csv('v11.3_PerClade_v2_wol_EXTENDED_pl_vs_contig_len.txt',sep=" ",header=TRUE)
pl_per_contig=read.csv('all.pl_err',sep=" ",header=FALSE)
head(pl_per_contig)

pl_per_contig$V1 <- gsub('.pl_err', '', pl_per_contig$V1)
pl_per_contig[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig$V1, '_', 2)
pl_per_contig[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig$greedy, '_', 2)
pl_per_contig[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig$chunks, '_', 2)
head(pl_per_contig)


pl_per_contig$chunks <- gsub('n', '', pl_per_contig$chunks)
pl_per_contig$seed <- gsub('s', '', pl_per_contig$seed)
pl_per_contig$chunk = as.numeric(pl_per_contig$chunk)
head(pl_per_contig)

factor(pl_per_contig$chunks)
pl_per_contig$chunks <- factor(pl_per_contig$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                '10', '20', '30', '40', '50', '100', ''))

factor(pl_per_contig$greedy)

ggplot(aes(x=chunks, y=V4, color = greedy),
       data=pl_per_contig[pl_per_contig$greedy == "ng",])+
  #stat_summary(geom="crossbar")+
  #stat_summary(geom="point")+
  #geom_violin()+
 stat_summary(fun.data = "mean_cl_boot",size = 1, alpha = 0.7)+
  coord_cartesian(ylim=c(0,10))
  #theme_classic()
  

# read k5

pl_per_contig_k5=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k5_v20_8k_s28_complete_q_v2/contigs/all.pl_err',sep=" ",header=FALSE)
head(pl_per_contig_k5)

pl_per_contig_k5$V1 <- gsub('.pl_err', '', pl_per_contig_k5$V1)
pl_per_contig_k5[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k5$V1, '_', 2)
pl_per_contig_k5[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k5$greedy, '_', 2)
pl_per_contig_k5[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k5$chunks, '_', 2)
head(pl_per_contig_k5)


pl_per_contig_k5$chunks <- gsub('n', '', pl_per_contig_k5$chunks)
pl_per_contig_k5$seed <- gsub('s', '', pl_per_contig_k5$seed)
pl_per_contig_k5$chunk = as.numeric(pl_per_contig_k5$chunk)
head(pl_per_contig_k5)

factor(pl_per_contig_k5$chunks)
pl_per_contig_k5$chunks <- factor(pl_per_contig_k5$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                '10', '20', '30', '40', '50', '100', ''))

factor(pl_per_contig_k5$greedy)



# read k6

pl_per_contig_k6=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k6_v20_8k_s28_complete_q_v2/contigs/all.pl_err',sep=" ",header=FALSE)
head(pl_per_contig_k6)

pl_per_contig_k6$V1 <- gsub('.pl_err', '', pl_per_contig_k6$V1)
pl_per_contig_k6[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k6$V1, '_', 2)
pl_per_contig_k6[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k6$greedy, '_', 2)
pl_per_contig_k6[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k6$chunks, '_', 2)
head(pl_per_contig_k6)


pl_per_contig_k6$chunks <- gsub('n', '', pl_per_contig_k6$chunks)
pl_per_contig_k6$seed <- gsub('s', '', pl_per_contig_k6$seed)
pl_per_contig_k6$chunk = as.numeric(pl_per_contig_k6$chunk)
head(pl_per_contig_k6)

factor(pl_per_contig_k6$chunks)
pl_per_contig_k6$chunks <- factor(pl_per_contig_k6$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                      '10', '20', '30', '40', '50', '100', ''))

factor(pl_per_contig_k6$greedy)

# read k8
pl_per_contig_k8=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k8_v20_8k_s28_complete_q_v2/contigs/all.pl_err',sep=" ",header=FALSE)
head(pl_per_contig_k8)

pl_per_contig_k8$V1 <- gsub('.pl_err', '', pl_per_contig_k8$V1)
pl_per_contig_k8[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k8$V1, '_', 2)
pl_per_contig_k8[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k8$greedy, '_', 2)
pl_per_contig_k8[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k8$chunks, '_', 2)
head(pl_per_contig_k8)


pl_per_contig_k8$chunks <- gsub('n', '', pl_per_contig_k8$chunks)
pl_per_contig_k8$seed <- gsub('s', '', pl_per_contig_k8$seed)
pl_per_contig_k8$chunk = as.numeric(pl_per_contig_k8$chunk)
head(pl_per_contig_k8)

factor(pl_per_contig_k8$chunks)
pl_per_contig_k8$chunks <- factor(pl_per_contig_k8$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                      '10', '20', '30', '40', '50', '100', ''))

factor(pl_per_contig_k8$greedy)


pl_per_contig['k'] = 7
pl_per_contig_k5['k'] = 5
pl_per_contig_k6['k'] = 6
pl_per_contig_k8['k'] = 8

pl_per_contig_concat = rbind(pl_per_contig_k5,pl_per_contig_k6, pl_per_contig, pl_per_contig_k8)

factor(pl_per_contig_concat$k)

pl_per_contig_concat


ggplot(aes(x=chunks, y=V4, color = factor(k)),
       data=pl_per_contig_concat[pl_per_contig_concat$greedy == "ng",])+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  #stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()+
  stat_summary(fun.data = "mean_cl_boot",size = 1, alpha = 0.7)+
  coord_cartesian(ylim=c(0,10))+
  scale_x_discrete(name="Chunk length (KB)", label = c("10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  guides(col= guide_legend(title= "k"))+
  #scale_fill_discrete(name = "New Legend Title")
  #guides(fill=guide_legend(title="k"))
  theme(legend.position = c(0.9,0.75), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
        )
  #scale_colour_brewer(palette = "Dark2", name="")

ggsave("incomplete_genomes_all_k.pdf",width=5,height = 4)
  

#theme_classic()



##################
# Multiple k placement error

# Claded Unchunked model with computed clade
df_k3=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k3_v20_8k_s28/all.pl_err',sep=" ",header=FALSE)
df_k4=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k4_v20_8k_s28/all.pl_err',sep=" ",header=FALSE)
df_k5=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k5_v20_8k_s28/all.pl_err',sep=" ",header=FALSE)
df_k6=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k6_v20_8k_s28/all.pl_err',sep=" ",header=FALSE)
#df_k7=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28/all.pl_err',sep=" ",header=FALSE)
df_k7=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24/all.pl_err',sep=" ",header=FALSE)
df_k8=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k8_v20_8k_s28/all.pl_err',sep=" ",header=FALSE)

df_k3['k'] = 3
df_k4['k'] = 4
df_k5['k'] = 5
df_k6['k'] = 6
df_k7['k'] = 7
df_k8['k'] = 8

df_k7$V1 <- gsub('.pl_err', '', df_k7$V1)
df_k7$V1 <- gsub('apples_input_di_mtrx_query_', '', df_k7$V1)
df_k7 <- df_k7[,colSums(is.na(df_k7))<nrow(df_k7)]
colnames(df_k7)[colnames(df_k7)=="V3"] <- "V2"
colnames(df_k7)[colnames(df_k7)=="V4"] <- "V3"
colnames(df_k7)[colnames(df_k7)=="V5"] <- "V4"

colnames(df_k7)
df_k7




colnames(df_k3)
colnames(df_k4)
colnames(df_k5)
colnames(df_k6)
colnames(df_k8)
colnames(df_k7)



df_all_concat = rbind(df_k3, df_k4, df_k5, df_k6, df_k7, df_k8)
means <- aggregate(V3 ~  k, df_all_concat, mean)

head(df_all_concat)

ggplot(aes(x=factor(k), y=V3, fill=factor(k)),
       data=df_all_concat)+
  #geom_boxplot(outlier.shape=NA)+
  geom_boxplot(fill="skyblue1", outlier.shape=NA)+
  geom_text(data = means, aes(label = V3, y = V3 - 0.4))+
  #geom_bar( stat="identity", fill="skyblue", alpha=0.7)
  #stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()+
  #stat_summary(fun.data = "mean",size = 1, alpha = 0.7)+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  coord_cartesian(ylim=c(0,8))+
  scale_x_discrete(name="k", label = c("3","4", "5", "6", "7", "8"))  +
  ylab("Placement error")+
  guides(col= guide_legend(title= "k"))+
  #scale_fill_discrete(name = "New Legend Title")
  #guides(fill=guide_legend(title="k"))
  theme(legend.position = "none")
  
#scale_colour_brewer(palette = "Dark2", name="")

ggsave("var_kmer_len.pdf",width=4.8,height = 4)



getwd()
ggplot(aes(x=factor(k), y=V3),
       data=df_all_concat)+
  #geom_boxplot(outlier.shape=NA)+
  #geom_boxplot(fill="skyblue1", outlier.shape=NA)+
  stat_summary(aes(fill=
                     ifelse(k =="7","Default","Tested")),geom="bar",color="black")+
  geom_text(data = means, aes(label = sprintf("%.1f", round(V3, digits = 2)), y = V3 - 0.45))+

  
  stat_summary()+
  #geom_bar( stat="identity", fill="skyblue", alpha=0.7)
  stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()+
  #stat_summary(fun.data = "mean",size = 1, alpha = 0.7)+
  #stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  coord_cartesian(ylim=c(0,6.5))+
  scale_x_discrete(name="k-mer length", label = c("3","4", "5", "6", "7", "8"))  +
  ylab("Placement error")+
  
  #guides(col= guide_legend(title= "k"))+
  #scale_fill_discrete(name = "New Legend Title")
  
  #scale_colour_brewer(palette = "Set1", name="", labels = c("Test", "Train"))+
  #scale_fill_brewer(palette = "Accent", name="")+
  #scale_fill_brewer(palette = "Dark2", name="")+
  scale_fill_manual(values = c( "#fddbc7", "#d1e5f0", "#e6f5d0"), name="")+
  #scale_fill_manual(values = c(  "#d1e5f0", "#f7f7f7", "#fddbc7"), name="")+
  #scale_fill_manual(values = c(  "#E7298A", "#66A61E", "#A6761D"), name="")+
  #scale_fill_manual(values = c(  "#BEAED4", "#7FC97F", "#FDC086"), name="")+
  #scale_fill_manual(values = c(  dark2_colors[1], dark2_colors[2], "#FDC086"), name="")+
  #guides(fill=guide_legend(title="k"))+
  theme(legend.title = element_blank())+
  
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5), 
        legend.position = c(.82,.88),legend.direction = "vertical", legend.margin=margin(t = -0.5, unit='cm'))

#axis.title.x = element_blank(),

#ggsave("var_kmer_len_v2.pdf",width=3.8,height = 4)
ggsave("var_kmer_len_v2.pdf",width=3.8,height = 4)

getwd()

###################################################### ################## ##################
###################################################### ################## ##################
###################################################### ################## ##################


# Comparison of full genome 500 queries on different models

# Uncladed, unchunked
df_UU =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_single_clade/all.pl_err',sep=" ",header=FALSE)
df_UU ["claded"] = "Uncladed"
df_UU ["chunked"] = "Unchunked"
df_UU ["clsf"] = "na"
df_UU
df_UU = df_UU[!(is.na(df_UU$V2)), ]
df_UU

summary(df_UU $V3)

# Uncladed, chunked MODEL
df_BU =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_single_clade_ChunkModel/all.pl_err',sep=" ",header=FALSE)
df_BU ["claded"] = "Uncladed"
df_BU ["chunked"] = "Chunked"
df_BU ["clsf"] = "na"
nrow(df_BU)

# Claded, unchunked - True clade
df_UC_tr =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_TRUE_CLADE_s24/all.pl_err',sep=" ",header=FALSE)
df_UC_tr ["claded"] = "Claded"
df_UC_tr ["chunked"] = "Unchunked"
df_UC_tr["clsf"] = "tr"
nrow(df_UC_tr)
df_UC_tr

# Claded, unchunked - Computed clade
df_UC_clsf =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24/all.pl_err',sep=" ",header=FALSE)
df_UC_clsf ["claded"] = "Claded"
df_UC_clsf ["chunked"] = "Unchunked"
df_UC_clsf["clsf"] = "cmp"
summary(df_UC_clsf$V4)

# Claded, chunked - True clade
df_BC_tr =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_Chunks_FullGenomQueries_TRUE_CLADE/all.pl_err',sep=" ",header=FALSE)
df_BC_tr ["claded"] = "Claded"
df_BC_tr ["chunked"] = "Chunked"
df_BC_tr["clsf"] = "tr"

# Claded, chunked - Computed clade
df_BC_clsf =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_FullGenomQueries_COMPUTED_CLADE/all.pl_err',sep=" ",header=FALSE)
df_BC_clsf ["claded"] = "Claded"
df_BC_clsf ["chunked"] = "Chunked"
df_BC_clsf["clsf"] = "cmp"

head(df_UU)
nrow(df_UU)



df_inter = rbind(df_BU, df_UC_tr, df_UC_clsf, df_BC_tr, df_BC_clsf)

head(df_inter)

df_inter
# Unify format for all dataframes
df_inter$V1 <- gsub('.pl_err', '', df_inter$V1)
df_inter$V1 <- gsub('apples_input_di_mtrx_query_', '', df_inter$V1)
df_inter <- df_inter[,colSums(is.na(df_inter))<nrow(df_inter)]
colnames(df_inter)[colnames(df_inter)=="V3"] <- "V2"
colnames(df_inter)[colnames(df_inter)=="V4"] <- "V3"
colnames(df_inter)[colnames(df_inter)=="V5"] <- "V4"
head(df_inter)

df_inter = df_inter[!(is.na(df_inter$V2)), ]


nrow(df_inter)
head(df_inter)
df_inter

df_fin = rbind(df_UU, df_inter)
df_fin$condition <- paste(df_fin$claded, "-", df_fin$chunked, "-", df_fin$clsf)

head(df_fin)

q = read.csv('/Users/nora/Documents/ml_metagenomics/tol_quality_scores/quality_comparison_hor.csv')
head(q)

#nov = read.csv('/Users/nora/Documents/ml_metagenomics/pendant.txt',sep="\t",h=F)
#nov$V2 = nov$V2*2/100
#head(nov)

nov = read.csv('/Users/nora/Documents/ml_metagenomics/closest_dev_set_2col.txt',sep="\t",h=FALSE)
#head(nov)
nov$V2 = nov$V2*1/100
head(nov)

names(nov)[2] = "nov"
df_fin = merge(df_fin,nov,by="V1")
head(df_fin)

qdf_fin = merge(q,df_fin, by.x="Assembly", by.y = "V1") 

head(qdf_fin)

ggplot(aes(x=condition, y=V3),
       data=qdf_fin)+
  #data=df_fin[df_fin$chunked %in% c("Unchunked"),])+
  #geom_boxplot(outlier.shape=NA)+
  #geom_boxplot(fill="skyblue1", outlier.shape=NA)+
  stat_summary(aes(fill=
                     ifelse(claded =="Claded","Claded","Uncladed")),geom="bar", )+
 #geom_text(data = df_fin$V3, aes(label = sprintf("%.1f", round(V3, digits = 2)), y = V3 - 0.45))+
  geom_text(aes(label = after_stat(sprintf("%.1f", round(y, digits = 2)))), stat = "summary", fun = "mean", vjust = 2.5, colour = "black" ) +
  
  facet_wrap(.~chunked)
  stat_summary()
  #geom_bar( stat="identity", fill="skyblue", alpha=0.7)
  stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()
  
ggplot(aes(x=condition, y=V3),
         data=qdf_fin[qdf_fin$chunked=="Unchunked" ,])+
    #data=df_fin[df_fin$chunked %in% c("Unchunked"),])+
    #geom_boxplot(outlier.shape=NA)+
    #geom_boxplot(fill="skyblue1", outlier.shape=NA)+
    stat_summary(aes(fill=
                       ifelse(claded =="Claded","Claded","Uncladed")),geom="bar", )+
    #geom_text(data = df_fin$V3, aes(label = sprintf("%.1f", round(V3, digits = 2)), y = V3 - 0.45))+
    geom_text(aes(label = after_stat(sprintf("%.1f", round(y, digits = 2)))), stat = "summary", fun = "mean", vjust = 2.5, colour = "black" ) +
    
    facet_wrap(.~chunked)+
  stat_summary()
  #geom_bar( stat="identity", fill="skyblue", alpha=0.7)
  stat_summary(geom="point")+
    #geom_violin()+
    theme_classic()
  
  quantile(nov$nov,(0:10)/10)
           
  ggplot(aes(color=reorder(paste(claded,ifelse(clsf=="tr","(true)",""),ifelse(chunked=="Chunked","(Chunked)","")),V3), 
                        y=V3,x=cut(nov,breaks=c(0,0.01,0.05,0.15,0.2,0.5,2) )),
         data=qdf_fin[qdf_fin$chunked=="Unchunked"| (qdf_fin$chunked=="Chunked" & qdf_fin$claded=="Claded" & qdf_fin$clsf =="scmp") ,])+
    #data=df_fin[df_fin$chunked %in% c("Unchunked"),])+
    #geom_boxplot(outlier.shape=NA)+
    #geom_boxplot(fill="skyblue1", outlier.shape=NA)+
    stat_summary(aes(group=condition),geom="line" )+
    #stat_summary(aes(group=condition),geom="bar" ,position = position_dodge(0.9),color="black")+
    #stat_summary(position = position_dodge(0.9),color="black")+
    stat_summary()+
    scale_fill_brewer(palette = "Paired",name = "")+
    scale_color_manual(labels = c("Cladded (true)", "Cladded","Uncladded"), name = "",values=c("#a6cee3","#1f78b4","#33a02c","red"))+
    theme_classic()+
    theme(legend.position = c(0.18,0.75))+
    scale_y_continuous("Placement error")+scale_x_discrete("Query novelty")
  
  ggsave("clading-novelty-D1-line.pdf",width = 6.2,height = 4)
  getwd()
  
    #geom_text(data = df_fin$V3, aes(label = sprintf("%.1f", round(V3, digits = 2)), y = V3 - 0.45))+
    #geom_text(aes(label = after_stat(sprintf("%.1f", round(y, digits = 2)))), stat = "summary", 
    #          fun = "mean", vjust = 2.5) +
    
    #facet_wrap(.~chunked)+



head(qdf_fin)
summary(qdf_fin[qdf_fin$condition=="Claded - Chunked - cmp",]$V3)
summary(qdf_fin[qdf_fin$condition=="Claded - Unchunked - cmp",]$V3)

 
 ggplot(aes(color=reorder(paste(claded,ifelse(clsf=="tr","(True)",""),ifelse(chunked=="Chunked","(chunked)","")),V3), 
             y=V3,x=cut(nov,breaks=c(0,0.01,0.05,0.15,0.2,0.5,2) )),
         data=qdf_fin[qdf_fin$clsf =="cmp" ,])+
    #data=df_fin[df_fin$chunked %in% c("Unchunked"),])+
    #geom_boxplot(outlier.shape=NA)+
    #geom_boxplot(fill="skyblue1", outlier.shape=NA)+
    stat_summary(aes(group=condition),geom="line" )+
    #stat_summary(aes(group=condition),geom="bar" ,position = position_dodge(0.9),color="black")+
    #stat_summary(position = position_dodge(0.9),color="black")+
    stat_summary()+
    scale_fill_brewer(palette = "Paired",name = "")+
    scale_color_manual(labels = c("Cladded","Cladded (chunked)"), name = "",values=c("#1f78b4","#e31a1c"))+
    theme_classic()+
    theme(legend.position = c(0.18,0.75))+
    scale_y_continuous("Placement error")+scale_x_discrete("Query novelty")
  ggsave("chunking-D1-line.pdf",width = 4.8,height = 4)
    
 
  ggplot(aes(x=cut(genes_retained_index,breaks=quantile(genes_retained_index,(0:5)/5), include.lowest = TRUE), y=V3,color=condition,group=condition),
         data=qdf_fin[qdf_fin$chunked=="Unchunked",])+
    #data=df_fin[df_fin$chunked %in% c("Unchunked"),])+
    #geom_boxplot(outlier.shape=NA)+
    #geom_boxplot(fill="skyblue1", outlier.shape=NA)+
    stat_summary(geom="line")+
    stat_summary()
  
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7, hjust=0.5), 
        legend.position = c(.82,.88),legend.direction = "vertical", legend.margin=margin(t = -0.5, unit='cm'))+
  coord_cartesian(ylim=c(0,6.5))  
stat_summary(aes(fill=
                     ifelse(k =="7","Default","Tested")),geom="bar",color="black")+
  #geom_text(data = means, aes(label = sprintf("%.1f", round(V3, digits = 2)), y = V3 - 0.45))+
  
  
  stat_summary()
  #geom_bar( stat="identity", fill="skyblue", alpha=0.7)
  stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()+
  #stat_summary(fun.data = "mean",size = 1, alpha = 0.7)+
  #stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  
  scale_x_discrete(name="k-mer length", label = c("3","4", "5", "6", "7", "8"))  +
  ylab("Placement error")+
  
  #guides(col= guide_legend(title= "k"))+
  #scale_fill_discrete(name = "New Legend Title")
  
  #scale_colour_brewer(palette = "Set1", name="", labels = c("Test", "Train"))+
  #scale_fill_brewer(palette = "Accent", name="")+
  #scale_fill_brewer(palette = "Dark2", name="")+
  scale_fill_manual(values = c( "#fddbc7", "#d1e5f0", "#e6f5d0"), name="")+
  #scale_fill_manual(values = c(  "#d1e5f0", "#f7f7f7", "#fddbc7"), name="")+
  #scale_fill_manual(values = c(  "#E7298A", "#66A61E", "#A6761D"), name="")+
  #scale_fill_manual(values = c(  "#BEAED4", "#7FC97F", "#FDC086"), name="")+
  
  #guides(fill=guide_legend(title="k"))+
  theme(legend.title = element_blank())+
  
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5), 
        legend.position = c(.82,.88),legend.direction = "vertical", legend.margin=margin(t = -0.5, unit='cm'))



  
  
# Combine with TRUE clade information  for VIOLIN plots 
df_true_clades =read.csv('/Users/nora/Documents/ml_metagenomics/clade_targets.txt',sep=" ",header=TRUE)
names(df_true_clades)[names(df_true_clades) == 'genome'] <- 'V1'

per_clade_error <- merge(df_fin,df_true_clades, by="V1")

head(per_clade_error)
head(df_fin)
  

a <- per_clade_error[per_clade_error$chunked %in% c("Unchunked") & !per_clade_error$clsf %in% c("na"),]

nrow(a)
tail(a)
ggplot(aes(x=clade, y=V3, color = condition, group = interaction(condition, clade)),
       #data=per_clade_error[per_clade_erro])+
       data=per_clade_error[per_clade_error$chunked %in% c("Unchunked") & ! per_clade_error$clsf %in% c("cmp"),])+
  geom_rect(xmin = 0.5,xmax = 1.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 2.5,xmax = 3.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 4.5,xmax = 5.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 6.5,xmax = 7.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 8.5,xmax = 9.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 10.5,xmax = 11.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 12.5,xmax = 13.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  #scale_color_brewer(palette = my_palette,name="")+
  #scale_color_manual(name="", values = c( "#ca0020", "#0571b0" ), 
  #                   labels = c("Global", "Local"))+
  scale_color_manual(name="", values = c( "#0571b0", "#ca0020" ), 
                     labels = c("Local", "Global"))+
  xlab("True clade")+
  ylab("Placement error")+
  #geom_violin()+
  theme_classic()+
  #geom_boxplot(alpha = 0.9, size = 0.5, outlier.shape=NA)+
  #geom_boxplot( alpha = 0.9, size = 0.5)+
  geom_violin(draw_quantiles = c(0.5), fun.data = function(x) median_hilow_(x,ci=0.8))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  #xlim(-0.1, 14.5)+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = c(.87,.9),legend.direction = "vertical")

ggsave("pl_per_clade_true_violin.pdf",width=6.2,height = 4)

q<-per_clade_error[per_clade_error$chunked %in% c("Unchunked") & ! per_clade_error$clsf %in% c("cmp"),]

for (x in 0:14) {

  q1 = q[q$clade==x,]
  q2 = q1[q1$claded =="Claded",]
  #q3 = q1[q1$claded =="Uncladed",]
  #summary(q2$V3)
  mean_value <- mean(q2$V3)
  median_value <- median(q2$V3)
  #sprintf("Mean: %.3g", mean_value); sprintf("Median: %.3g", median_value)
  print (sprintf("%.3g, %.3g", mean_value, median_value))
}
sd(q2$V3, na.rm = FALSE)
summary(q3$V3)
sd(q3$V3, na.rm = FALSE)
q3[q3$V3>5,]
q2[q2$V3>5,]

summary(q[q$claded =="Claded",]$V3)



ggplot(aes(x=reorder(clade,V3), y=V3, color = condition, group = interaction(condition)),
       #data=per_clade_error)+
  data=per_clade_error[per_clade_error$chunked %in% c("Unchunked") & ! per_clade_error$clsf %in% c("cmp"),])+
  geom_rect(xmin = 0.5,xmax = 1.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 2.5,xmax = 3.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 4.5,xmax = 5.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 6.5,xmax = 7.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 8.5,xmax = 9.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 10.5,xmax = 11.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 12.5,xmax = 13.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  #scale_color_brewer(palette = my_palette,name="")+
  scale_color_manual(name="", values = c(  "#0571b0","#ca0020" ), 
                     labels = c( "Cladded (true)", "Uncladded"))+
  xlab("True clade")+
  ylab("Placement error")+
  #geom_violin()+
  theme_classic()+
  stat_summary(geom="line")+
  stat_summary(fun.data=function(x) mean_ci(x,ci=0.99),
               position=position_dodge(width=0.3))+
  stat_summary(fun.y =function(x) median(x),
               position=position_dodge(width=0.3),geom="point",
               shape=5,size=2)+
  coord_cartesian(ylim=c(0,7))+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = c(.2,.9),legend.direction = "vertical", legend.title = element_blank(),
        legend.margin=margin(0,0,0,0))

ggsave("pl_per_clade_true_violin_v3.pdf",width=6.2,height = 4)
getwd()
?mean_ci

ggplot(aes(group=interaction(clade,condition), x=V3, color = condition),
       #data=per_clade_error)+
       data=per_clade_error[per_clade_error$chunked %in% c("Unchunked") & ! per_clade_error$clsf %in% c("cmp"),])+
  #scale_color_brewer(palette = my_palette,name="")+
  scale_color_manual(name="", values = c(  "#0571b0","#ca0020" ), 
                     labels = c( "Claded (true)", "Uncladed"))+
  xlab("True clade")+
  ylab("Placement error")+
  #geom_violin()+
  theme_classic()+
  stat_ecdf()+
  #coord_cartesian(ylim=c(0,7))+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = c(.2,.9),legend.direction = "vertical", legend.title = element_blank(),
        legend.margin=margin(0,0,0,0))

ggplot(aes(x=V3, color = condition),
       #data=per_clade_error)+
  data=per_clade_error[per_clade_error$chunked %in% c("Unchunked") & ! per_clade_error$clsf %in% c("cmp"),])+
  #scale_color_brewer(palette = my_palette,name="")+
  scale_color_manual(name="", values = c(  "#0571b0", "#ca0020" ), 
                     labels = c( "Unchunked_claded (true clade)", "Unchunked_uncladed"))+
  xlab("True clade")+
  ylab("Placement error")+
  #geom_violin()+
  theme_classic()+
  stat_ecdf(aes(group = interaction(condition,clade)),alpha=0.4,size=0.3)+
  stat_ecdf(aes(group = interaction(condition)),size=1.2)+
  coord_cartesian(ylim=c(0,1))+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = c(.7,.2),legend.direction = "vertical")


ggplot(aes(x=clade, y=V3, color = condition, group = clade),
       #data=per_clade_error)+
  data=per_clade_error[per_clade_error$chunked %in% c("Unchunked") & ! per_clade_error$clsf %in% c("cmp"),])+
  facet_wrap(~condition)+
  scale_color_manual(name="", values = c(  "#0571b0", "#ca0020" ), 
                     labels = c("Unchunked_claded (true clade)", "Unchunked_uncladed"),)+
  geom_violin(alpha = 0.9, size = 0.5, draw_quantiles = c(0.25, 0.5, 0.75),  fun.data = function(x) median_hilow_(x,ci=0.8))+
  #scale_color_brewer(palette = my_palette,name="")+
  theme_classic()+
  #geom_boxplot(alpha = 0.9, size = 0.5, outlier.shape=NA)+
  #geom_violin(draw_quantiles = c(0.5), fun.data = function(x) median_hilow_(x,ci=0.8))+
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  #ylim(NA, 10)+
  #theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5),
  #      legend.position = c(.87,.8),legend.direction = "vertical")
  theme(legend.position = "none")+
  xlab("Clade number")+
  ylab("Placement error")
ggsave("pl_per_clade_true_violin_v2.pdf",width=6.2,height = 4)



# Parameter titration on Claded unchunked model based on BL k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24
# Classification was based on results from k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24


# Claded, unchunked - Computed clade
df_bl =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_bl["cond"] = "BL"
df_bl["cond2"] = "BL"

ggplot(data=df_bl, aes(x=pl_err))+
  stat_ecdf()+
  theme_bw()+
  geom_vline(color="red",xintercept = mean(df_bl$pl_err,na.rm=T))+
  scale_y_continuous(labels = percent, "ECDF")+
  scale_x_continuous(name="Placement error",trans = "sqrt")

mean(df_bl$pl_err<=0)
mean(df_bl$pl_err<=1)
mean(df_bl$pl_err<=2)
mean(df_bl$pl_err>=8)

ggsave("ECDF_pl_err_claded_unchunked_D1.pdf",width=6.5,height = 4.5)
getwd()

my_misclassified = df_bl[df_bl$top_class!=df_bl$true_clade,]
nrow(my_misclassified)

my_misclassified_qdf_fin = merge(my_misclassified,qdf_fin, by.x="genome_x", by.y = "Assembly", all.x = FALSE) 
nrow(my_misclassified_qdf_fin)

my_misclassified_qdf_fin2 = my_misclassified_qdf_fin[my_misclassified_qdf_fin$condition=="Claded - Unchunked - cmp",]


my_misclassified_qdf_fin2$nov

df_param1 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_extra_layer/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_param1["cond"] = "extra_layer"
df_param1["cond2"] = "layer"

df_param2 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_extra_2layer/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_param2["cond"] = "extra_2layer"
df_param2["cond2"] = "layer"


df_param3 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_hilr/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_param3["cond"] = "hilr"
df_param3["cond2"] = "lr"

df_param4 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_hilr2/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_param4["cond"] = "hilr2"
df_param4["cond2"] = "lr"

df_param5 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_hilr_hilrmin/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_param5["cond"] = "hilr_hilrmin"
df_param5["cond2"] = "lrmin"

df_param6 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_hilr2_hilrmin/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_param6["cond"] = "hilr2_hilrmin"
df_param6["cond2"] = "lrmin"

df_param7 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_batch32/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_param7["cond"] = "batch32"
df_param7["cond2"] = "other"

df_param8 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_decay4k/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_param8["cond"] = "decay4k"
df_param8["cond2"] = "other"


df_parameters = rbind(df_bl, df_param1,  df_param2, df_param3, df_param4, df_param5, df_param6, df_param7, df_param8)
head(df_parameters)

df_parameters$cond <- factor(df_parameters$cond, levels = c('BL', 'extra_layer', 'extra_2layer', 'hilr', 'hilr_hilrmin', 'hilr2', 'hilr2_hilrmin', 'decay4k', 'batch32'))



ggplot(aes(x=cond, y=pl_err, fill=cond2, group = cond2),
       data=df_parameters)+
  #data=df_fin[df_fin$chunked %in% c("Unchunked"),])+
  #geom_boxplot(outlier.shape=NA)+
  #geom_boxplot(fill="skyblue1", outlier.shape=NA)+
  stat_summary(geom="bar", alpha = 0.7, color="black",)+
  #geom_text(data = df_fin$V3, aes(label = sprintf("%.1f", round(V3, digits = 2)), y = V3 - 0.45))+
  geom_text(aes(label = after_stat(sprintf("%.1f", round(y, digits = 2)))), stat = "summary", fun = "mean", vjust = 2.5, colour = "black" ) +
  stat_summary()+
  #scale_fill_manual(palette="Dark2")+
  #scale_fill_manual(palette = "Dark2", name="")+
  scale_fill_brewer(name = "", palette="Dark2", labels = c("Default", "Model structure", "Learning rate (lr)", "Learning rate min (lr min)", "Other"))+
  #scale_color_brewer(name = "",palette="Dark2", labels = c("kf2d", "DEPP_381", "DEPP_1"))+
#geom_bar( stat="identity", fill="skyblue", alpha=0.7)
stat_summary(geom="point")+
  #scale_x_discrete(labels = c("Default", "Extra\nlayer", "Extra\ntwo\nlayers", "10x\nlearning\nrate", "100x\nlearning\nrate", "10x\nlearning\nrate min", "100x\nlearning\nrate min", "2x\nbatch\nsize", "2x\nlearning\ndecay"), name = "") +
  scale_x_discrete(labels = c("Default", "Extra layer", "Extra two layers", "10x lr", "10x lr min", "100x lr",  "100x lr min", "2x batch size", "2x lr decay"), name = "") +
  #geom_violin()+
  theme_classic()+
  #theme(legend.position = "none")+
  #xlab("Clade number")+
  ylab("Placement error")+
  #theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 30, vjust = 0.7, hjust=0.7),
      theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 30, vjust = 1.0, hjust=1.0),
        #legend.position = c(.2,.9),
        legend.margin=margin(0,0,0,0),
        legend.key.size = unit(0.5, "cm"),
        legend.direction = "horizontal", legend.position = "bottom")

ggsave("parameters_claded_unchunked_D1_bar.pdf",width=6.5,height = 4.5)

getwd()
##############################################################
# Build train test error plot


err_epoch=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v61_8k_train_model_10KFULL_TOL_Claded_Unchunked_tr_test_err/train_test_err_claded_unchunked_clades_all.txt',sep=" ",header=T)


names(err_epoch)
unique(err_epoch$clade)

ggplot(aes(x=epoch, y=loss_val, color = type),
       data=err_epoch[err_epoch$epoch<= 8000,])+
  #geom_point(alpha = 0.9, size = 0.4)+
  stat_summary(geom="point",alpha = 0.4, size = 0.4)+
  #geom_smooth(se=F,  size=0.5, alpha = 0.9, method = "lm",formula = y ~ poly(x, 12), linetype = "dotted")+
  geom_smooth(aes(group=type),se=F,  size=1, alpha = 1, linetype =1,color="black", formula = y~ poly(x, 8))+
  scale_colour_brewer(palette = "Dark2", name="", labels = c("Test", "Train"))+
  #facet_wrap(.~clade, scales = "free", ncol = 3)+
  theme_bw()+
  xlab("Epoch")+
  ylab("Error")+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = c(.90,.90),legend.direction = "vertical",legend.margin=margin(c(0,0,0,0)),
        legend.title=element_blank())+
  scale_y_log10()+
  scale_x_continuous(trans="sqrt",breaks=c(0,100,500,1000,2000,4000,8000))+
  guides(colour = guide_legend(reverse=TRUE))
  #coord_cartesian(ylim=c(NA, 0.05))

#scale_color_manual(values = "Dark2")
ggsave("dev_queries_err_epoch_claded_unchunked.pdf",width=6.5,height = 4.5)

getwd()

ggplot(aes(x=epoch, y=loss_val, color = type),
       data=err_epoch[err_epoch$epoch<= 8000 & err_epoch$clade=="3",])+
  geom_point(alpha = 0.9, size = 0.4)+
  #geom_smooth(se=F,  size=0.5, alpha = 0.9, method = "lm",formula = y ~ poly(x, 12), linetype = "dotted")+
  #geom_smooth(se=F,  size=0.5, alpha = 0.9, method = "lm", linetype = "dotted")+
  scale_colour_brewer(palette = "Dark2", name="", labels = c("Test", "Train"))+
  #facet_wrap(.~clade, scales = "free", ncol = 3)+
  theme_classic()+
  xlab("Epoch")+
  ylab("Error")+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = c(.9,.9),legend.direction = "vertical",legend.margin=margin(c(0,0,0,0)),
        legend.title=element_blank())+
  scale_y_log10()+
  guides(colour = guide_legend(reverse=TRUE))
#coord_cartesian(ylim=c(NA, 0.05))

#scale_color_manual(values = "Dark2")


ggsave("dev_queries_err_epoch_claded_unchunked_clade3.pdf",width=6.5,height = 4.5)
#ggsave("dev_queries_err_epoch_claded_unchunked_clade3.pdf",width=4.8,height = 4)


# Plot placement error change per epoch
plerror_epoch=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v62_8k_train_model_10K_TOL_pl_error_per_epoch/summary_D1_TOL_pl_error_per_epoch.csv',sep=" ",header=T)
unique(plerror_epoch$epoch)

plerror_epoch$epoch <- factor(plerror_epoch$epoch, levels = c('1', '501', '1001', '1501', '2001', '2501', '3001', '3501', '4001', 
                                                  '4501', '5001', '5501', '6001', '6501','7001', '7501', '8000', 'best'))
head(plerror_epoch)

ggplot(aes(x=epoch, y=pl_err, color = claded, group = claded),
       data=plerror_epoch[!plerror_epoch$epoch=="best", ])+
       #data=plerror_epoch)+
  #geom_line(aes(group=claded))+
  #geom_point() + 
  #geom_smooth(aes(group=claded), method = "lm", formula = y ~ poly(x, 12), linetype = "dotted")+
  #geom_point(alpha = 0.9, size = 0.4)+
  #stat_summary(geom="line")+
  stat_summary(fun.data=function(x) mean_ci(x,ci=0.99),
               position=position_dodge(width=0.3))+
  #geom_smooth(se=F,  size=0.5, alpha = 0.9, method = "lm", linetype = "solid")+
  geom_smooth(aes(group=claded), se=F,  size=0.5, alpha = 0.9, method = "lm",formula = y ~ poly(x, 10), linetype = "dotted")+
  #geom_smooth(se=F,  size=0.5, alpha = 0.9, method = "lm", linetype = "dotted")+
  scale_colour_brewer(palette = "Dark2", name="", labels = c("Claded", "Uncladed" ))+
  scale_color_manual(labels = c("Cladded", "Uncladded"), values = c(  dark2_colors[2], dark2_colors[3]))+
  #facet_wrap(.~clade, scales = "free", ncol = 3)+
  theme_classic()+
  scale_y_continuous("Placement error")+scale_x_discrete("Epoch")+
  scale_x_discrete(name="Epoch", label = c("1", "500", "1000", "1500", "2000", "2500", "3000", "3500", "4000", "4500","5000", "5500", "6000", "6500", "7000", "7500", "8000", "best"))  +
  coord_cartesian(ylim=c(1.2,4.4))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7, hjust=0.5),
        legend.position = c(.9,.9),legend.direction = "vertical",legend.margin=margin(c(0,0,0,0)),
        legend.title=element_blank())+
  #scale_y_log10()+
  guides(colour = guide_legend(reverse=FALSE))
#coord_cartesian(ylim=c(NA, 0.05))

#scale_color_manual(values = "Dark2")


ggsave("dev_queries_placementerr_epoch_D1_unchunked.pdf",width=6.5,height = 4.5)




getwd()


###################################################### ################## ##################
###################################################### ################## ##################
###################################################### ################## ##################


# Comparison v2 and v3 (entire range) of chunked genome queries on different models

# Uncladed, unchunked
df_UU =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_Global_qChunks_v2_3_SINGLE_CLADE/all.pl_err',sep=" ",header=FALSE)
df_UU ["claded"] = "Uncladed"
df_UU ["chunked"] = "Unchunked"
df_UU ["clsf"] = "na"
nrow(df_UU)
#df_UU = df_UU[!(is.na(df_UU$V2)), ]

# Uncladed, chunked MODEL
df_BU =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_Global_Chunks_qChunks_v2_3_SINGLE_CLADE/all.pl_err',sep=" ",header=FALSE)
df_BU ["claded"] = "Uncladed"
df_BU ["chunked"] = "Chunked"
df_BU ["clsf"] = "na"
nrow(df_BU)

# Claded, unchunked - True clade
df_UC_tr =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_Clades_s24_qChunks_v2_3_TRUE_CLADE/all.pl_err',sep=" ",header=FALSE)
df_UC_tr ["claded"] = "Claded"
df_UC_tr ["chunked"] = "Unchunked"
df_UC_tr["clsf"] = "tr"
nrow(df_UC_tr)
df_UC_tr

# Claded, unchunked - Computed clade
df_UC_clsf =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_Clades_s24_qChunks_v2_3_COMPUTED_CLADE/all.pl_err',sep=" ",header=FALSE)
df_UC_clsf ["claded"] = "Claded"
df_UC_clsf ["chunked"] = "Unchunked"
df_UC_clsf["clsf"] = "cmp"


# Claded, chunked - True clade
df_BC_tr =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_qChunks_v2_3_TRUE_CLADE/all.pl_err',sep=" ",header=FALSE)
df_BC_tr ["claded"] = "Claded"
df_BC_tr ["chunked"] = "Chunked"
df_BC_tr["clsf"] = "tr"

# Claded, chunked - Computed clade
df_BC_clsf =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_qChunks_v2_3_COMPUTED_CLADE/all.pl_err',sep=" ",header=FALSE)
df_BC_clsf ["claded"] = "Claded"
df_BC_clsf ["chunked"] = "Chunked"
df_BC_clsf["clsf"] = "cmp"

# Claded, chunked - Computed clade - Previous version of classification without EXP
df_BC_clsf_prev =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_Chunks_qChunks_v2_3_COMPUTED_CLADE/all.pl_err',sep=" ",header=FALSE)
df_BC_clsf_prev ["claded"] = "Claded"
df_BC_clsf_prev ["chunked"] = "Chunked"
df_BC_clsf_prev["clsf"] = "cmpNotExp"


df_inter = rbind(df_UU, df_BU, df_UC_tr, df_UC_clsf, df_BC_tr, df_BC_clsf, df_BC_clsf_prev)

head(df_inter)
nrow(df_inter)

# Unify format for all dataframes
df_inter$V1 <- gsub('.pl_err', '', df_inter$V1)
df_inter$V1 <- gsub('apples_input_di_mtrx_query_', '', df_inter$V1)
df_inter <- df_inter[,colSums(is.na(df_inter))<nrow(df_inter)]
colnames(df_inter)[colnames(df_inter)=="V3"] <- "V2"
colnames(df_inter)[colnames(df_inter)=="V4"] <- "V3"
colnames(df_inter)[colnames(df_inter)=="V5"] <- "V4"
head(df_inter)

df_inter = df_inter[!(is.na(df_inter$V2)), ]
tail (df_inter)

nrow(df_inter)
head(df_inter)
df_inter

df_fin = df_inter
df_fin$cond <- paste(df_fin$claded, "-", df_fin$chunked, "-", df_fin$clsf)

tail(df_fin)



df_fin[c('genome', 'greedy')] <- str_split_fixed(df_fin$V1, '_', 2)
tail(df_fin)
df_fin[c('greedy', 'chunks')] <- str_split_fixed(df_fin$greedy, '_', 2)
tail(df_fin)
tail(df_fin)
df_fin[c('chunks', 'seed')] <- str_split_fixed(df_fin$chunks, '_', 2)
head(df_fin)

df_fin$chunks <- gsub('n', '', df_fin$chunks)
head(df_fin)
df_fin$seed <- gsub('s', '', df_fin$seed)
head(df_fin)

unique(df_fin$chunks)

#df_fin$chunks[df_fin$chunks == '01'] <- 0.1
#df_fin$chunks[df_fin$chunks == '02'] <- 0.2
#df_fin$chunks[df_fin$chunks == '03'] <- 0.3
#df_fin$chunks[df_fin$chunks == '04'] <- 0.4
#df_fin$chunks[df_fin$chunks == '05'] <- 0.5
#df_fin$chunks[df_fin$chunks == '06'] <- 0.6
#df_fin$chunks[df_fin$chunks == '07'] <- 0.7
#df_fin$chunks[df_fin$chunks == '08'] <- 0.8
#df_fin$chunks[df_fin$chunks == '09'] <- 0.9






#df_fin$chunks = as.numeric(df_fin$chunks)
tail(df_fin)

factor(df_fin$chunks)
df_fin$chunks <- factor(df_fin$chunks, levels = c('01', '02', '03', '04', '05', '06', '07', '08', '09', 
                                                  '1', '2', '3', '4', '5','10', '20', '30', '40', '50', '100', ''))
df_fin$chunks
factor(df_fin$greedy)

unique(df_fin$chunks)
head(df_fin)



ggplot(aes(x=chunks, y=V3, color = cond),
       data=df_fin[df_fin$greedy == "ng",])+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  #stat_summary(geom="point")+
  #geom_violin()+
  geom_hline(yintercept = 3.0, col = "grey", linetype='dotted')+
  theme_classic()+
  stat_summary(fun.data = "mean_cl_boot",size = 0.8, alpha = 0.5)+
  coord_cartesian(ylim=c(0,24))+
  scale_x_discrete(name="Chunk length (KB)", label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  guides(col= guide_legend(title= "k"))+
  #scale_fill_discrete(name = "New Legend Title")
  #guides(fill=guide_legend(title="k"))
  theme(legend.position = c(0.8,0.77), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )


e1 = df_fin[df_fin$greedy == "ng" & (df_fin$clsf!="tr" |df_fin$chunked=="Chunked") & !grepl("Exp",df_fin$clsf),]

e2 = e1[e1$chunks %in% c("01", "02", "03", "04", "05", "06", "07", "08", "09", "1"),]
unique(e2$cond)
e3= e2[e2$cond=="Claded - Chunked - cmp",]
e4= e2[e2$cond=="Claded - Chunked - tr",]
summary(e3$V3)
summary(e4$V3)

e5 = e1[!e1$chunks %in% c("01", "02", "03", "04", "05", "06", "07", "08", "09", "1", "full"),]
head(e5)
e6= e5[e5$cond=="Claded - Chunked - cmp",]
summary(e6$V3)


ggplot(aes(x=chunks, y=V3, color = chunked,linetype=paste(claded,ifelse(clsf=="tr","(true)","")),group=cond),
       data=df_fin[df_fin$greedy == "ng" & (df_fin$clsf!="tr" |df_fin$chunked=="Chunked") & !grepl("Exp",df_fin$clsf),])+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  stat_summary(geom="line")+
  #geom_violin()+
  theme_classic()+
  stat_summary( alpha = 0.7,geom="point")+
  stat_summary( alpha = 0.7,geom="errorbar",width=0.1,linetype=1,size=0.3)+
  #coord_cartesian(ylim=c(0,30))+
  scale_x_discrete(name="Fragment length (KB)", label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  scale_colour_brewer(palette = "Dark2", name="")+
  scale_linetype(name="", labels = c("Cladded", "Cladded (true)", "Uncladded"))+
  #facet_wrap(~claded)+
  #guides(col= guide_legend(title= "Conditions"))+
  #scale_fill_discrete(name = "New Legend Title")
  #guides(fill=guide_legend(title="k"))
  theme(legend.position = c(0.75,0.7), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )
#scale_colour_brewer(palette = "Dark2", name="")

ggsave("chunks-D3-line.pdf",width=6.5,height = 4.5)


# df_fin[df_fin$greedy == "ng" & (df_fin$clsf!="tr" |df_fin$chunked=="Chunked") & !grepl("Exp",df_fin$clsf),])
df_fin[df_fin$cond %in% c("Claded - Chunked - cmp", "Uncladed - Unchunked - na", "Uncladed - Chunked - na", "Claded - Unchunked - cmp", "Claded - Chunked - tr","Claded - Unchunked - tr") ,] %>%
  filter(greedy == "ng" ) %>%
  mutate(Claded=paste(claded,ifelse(grepl("- tr",cond),"(true)","") ))  %>%
  select(V3,Claded,chunked,chunks) %>%
  dplyr::group_by(Claded,chunked,chunks) %>%
  dplyr::summarise(merror=mean(V3)) %>%
  #pivot_wider(names_from = Chunked,values_from = merror) %>%
  ggplot(aes(color=Claded, 
             x=chunks, xend=chunks,
             group=interaction(Claded,chunks),
             y=merror))+
  geom_line(arrow = arrow(length=unit(0.2,"cm"), 
                          ends="first", type = "closed"), 
            position = position_dodge(width=0.3))+
  scale_colour_brewer(palette = "Dark2", name="", labels = c("Cladded", "Cladded (true)", "Uncladded"))+
  scale_linetype(name="")+
  #scale_color_manual(name = "",values=c("#1f78b4","#e31a1c"))+
  theme_classic()+
  theme(legend.position = c(0.88,0.75))+
  scale_y_continuous("Placement error")+scale_x_discrete("Contig length (KB)")+
  scale_x_discrete(name="Fragment length (KB)", label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  coord_cartesian(ylim=c(0,25))+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust=1), axis.title.x = element_text(vjust=2.5))+
  theme(plot.margin = margin(5.5,5.5,-1.5,5.5, "pt"))
ggsave("Chunks-D3-arrow.pdf",width=4.8,height = 4.0)


unique(df_fin$cond)
summary(df_fin[df_fin$greedy == "ng" & df_fin$cond %in% c("Claded - Chunked - cmp"),]$V3)
summary(df_fin[df_fin$greedy == "ng" & df_fin$cond %in% c("Uncladed - Chunked - na"),]$V3 )
head(df_fin)

###


# ADD CLASSIFICATION ACCURACY TO THE FRAGMENTED DATA
# Claded, chunked - Computed clade
c1 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_qChunks_v2_3_COMPUTED_CLADE/classes.out',header=TRUE, sep = '\t')
#c1 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_qChunks_v2_3_TRUE_CLADE/assigned_classes.out',header=TRUE, sep = '\t')
#c1 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_qChunks_v2_3_TRUE_CLADE/single_classes.out',header=TRUE, sep = '\t')

d1 = df_fin[df_fin$greedy == "ng" & df_fin$cond %in% c("Claded - Chunked - cmp"),]

#d1[d1$V1=="G000017265_ng_n04_s44",]$V1
#c1[c1$genome=="G000017265_ng_n04_s44",]$genome

t = read.csv('/Users/nora/Documents/ml_metagenomics/clade_targets.txt',header=TRUE, sep = ' ')


cd1 = merge(x = c1, y = d1, by.y = "V1", by.x = "genome", by = NULL)

#cd1[cd1$genome=="G000017265_ng_n04_s44",]


cdt1 = merge(x = cd1, y = t, by.x = "genome.y", by.y = "genome", by = NULL)
head(cdt1)
nrow(cdt1[cdt1$top_class==cdt1$clade,])


# Claded, unchunked - Computed clade
c2 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_Clades_s24_qChunks_v2_3_COMPUTED_CLADE/classes.out',header=TRUE, sep = '\t')
d2 = df_fin[df_fin$greedy == "ng" & df_fin$cond %in% c("Claded - Unchunked - cmp"),]

cd2 = merge(x = c2, y = d2, by.y = "V1", by.x = "genome", by = NULL)
cdt2 = merge(x = cd2, y = t, by.x = "genome.y", by.y = "genome", by = NULL)



cd_all = rbind(cdt1, cdt2)
head(cd_all)

nrow(cd_all[cd_all$top_class==cd_all$clade,]) 

s1=cd_all[cd_all$chunks=='07',]
nrow(s1)
s1c = s1[s1$cond=="Claded - Chunked - cmp",]
nrow(s1c[s1c$top_class==s1c$clade,])/nrow(s1c)


ggplot(aes(x=chunks, y=ifelse(top_class==clade,1,0), color = chunked, group=cond),
       data=cd_all)+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  stat_summary(geom="line")+
  stat_summary( alpha = 0.7,geom="point")+
  stat_summary( alpha = 0.7,geom="errorbar",width=0.1,linetype=1,size=0.3)+
  scale_colour_brewer(palette = "Dark2", name="")+
  scale_linetype(name="")+
  #geom_violin()+
  theme_classic()+
  #coord_cartesian(ylim=c(0,30))+
  scale_x_discrete(name="Fragment length (KB)", label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Classification accuracy")+
  scale_colour_brewer(palette = "Dark2", name="")+
  scale_linetype(name="")+
  #facet_wrap(~claded)+
  #guides(col= guide_legend(title= "Conditions"))+
  #scale_fill_discrete(name = "New Legend Title")
  #guides(fill=guide_legend(title="k"))
  theme(legend.position = c(0.75,0.7), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )


ggsave("Chunks-D4-classified-line.pdf",width=6.5,height = 4.5)


getwd()


unique(df_fin$cond)
unique(df_fin$chunks)
summary(df_fin[df_fin$cond=="Claded - Chunked - cmp" & df_fin$chunks=="06",]$V3)
summary(df_fin[df_fin$cond=="Claded - Unchunked - cmp" & df_fin$chunks=="06",]$V3)

ggplot(aes(color=chunks, x=V3, ),
data=df_fin[df_fin$greedy == "ng" & df_fin$cond =="Claded - Chunked - cmp",])+
  stat_ecdf()+
  theme_classic()+
  theme(legend.position = c(0.68,0.25),legend.direction = "vertical",
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_text( hjust = 0, vjust = 1.4),
        
        #legend.margin=margin(c(0.0,0.0, 0.0,0.0)),
        #legend.box.margin=margin(c(0.0,0.0, -30.0,0.0)),
        legend.text = element_text(size = 10),
        #legend.spacing.y = unit(0.1, "cm"),
  )+
  #scale_color_viridis_d(name="Con.\nlen.\n(KB)")+
  scale_color_viridis_d(name="Contig length (KB)")+
  scale_x_continuous("Placement error")+
  ylab("ECDF")+
  guides(color=guide_legend(ncol=2, byrow = TRUE))




head(df_fin)


###################################################### ################## ##################
###################################################### ################## ##################
###################################################### ################## ##################
# OBSOLETE  # OBSOLETE  # OBSOLETE  # OBSOLETE  # OBSOLETE  # OBSOLETE  # OBSOLETE  # OBSOLETE 

# Training using chunked genomes


# read k7 with default training

pl_per_contig_k7=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_complete_q_v2/contigs/all_subset.pl_err',sep=" ",header=FALSE)
head(pl_per_contig_k7)

pl_per_contig_k7$V1 <- gsub('.pl_err', '', pl_per_contig_k7$V1)
pl_per_contig_k7[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7$V1, '_', 2)
pl_per_contig_k7[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7$greedy, '_', 2)
pl_per_contig_k7[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7$chunks, '_', 2)
tail(pl_per_contig_k7)


pl_per_contig_k7$chunks <- gsub('n', '', pl_per_contig_k7$chunks)
pl_per_contig_k7$seed <- gsub('s', '', pl_per_contig_k7$seed)
pl_per_contig_k7$chunk = as.numeric(pl_per_contig_k7$chunk)
head(pl_per_contig_k7)

factor(pl_per_contig_k7$chunks)
pl_per_contig_k7$chunks <- factor(pl_per_contig_k7$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                      '10', '20', '30', '40', '50', '100', ''))

factor(pl_per_contig_k7$greedy)


# Read k7 trained with chunks for clade 1, 5 and 6


pl_per_contig_k7chunks=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_complete_q_v2_CHUNKS/contigs/all_subset.pl_err',sep=" ",header=FALSE)
head(pl_per_contig_k7chunks)

pl_per_contig_k7chunks$V1 <- gsub('.pl_err', '', pl_per_contig_k7chunks$V1)
pl_per_contig_k7chunks[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7chunks$V1, '_', 2)
pl_per_contig_k7chunks[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7chunks$greedy, '_', 2)
pl_per_contig_k7chunks[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7chunks$chunks, '_', 2)
head(pl_per_contig_k7chunks)


pl_per_contig_k7chunks$chunks <- gsub('n', '', pl_per_contig_k7chunks$chunks)
pl_per_contig_k7chunks$seed <- gsub('s', '', pl_per_contig_k7chunks$seed)
pl_per_contig_k7chunks$chunk = as.numeric(pl_per_contig_k7chunks$chunk)
head(pl_per_contig_k7chunks)

factor(pl_per_contig_k7chunks$chunks)
pl_per_contig_k7chunks$chunks <- factor(pl_per_contig_k7chunks$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                      '10', '20', '30', '40', '50', '100', ''))

factor(pl_per_contig_k7chunks$greedy)





pl_per_contig_k7['cond'] = "default"
pl_per_contig_k7chunks['cond'] = "chunks"

pl_per_contig_concat = rbind(pl_per_contig_k7,pl_per_contig_k7chunks)

#factor(pl_per_contig_concat$k)

pl_per_contig_concat
tail(pl_per_contig_concat)

ggplot(aes(x=chunks, y=V4, color = cond),
       data=pl_per_contig_concat[pl_per_contig_concat$greedy == "ng",])+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  #stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()+
  stat_summary(fun.data = "mean_cl_boot",size = 1, alpha = 0.7)+
  coord_cartesian(ylim=c(0,10))+
  scale_x_discrete(name="Chunk length (KB)", label = c("10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  guides(col= guide_legend(title= "k"))+
  #scale_fill_discrete(name = "New Legend Title")
  #guides(fill=guide_legend(title="k"))
  theme(legend.position = c(0.9,0.75), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )
#scale_colour_brewer(palette = "Dark2", name="")

ggsave("incomplete_genomes_all_k.pdf",width=5,height = 4)





######################### Comparing different stopping criteria #########################


# read k7 with default training

pl_per_contig_k7=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_complete_q_v2/contigs/all_subset.pl_err',sep=" ",header=FALSE)
head(pl_per_contig_k7)

pl_per_contig_k7$V1 <- gsub('.pl_err', '', pl_per_contig_k7$V1)
pl_per_contig_k7[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7$V1, '_', 2)
pl_per_contig_k7[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7$greedy, '_', 2)
pl_per_contig_k7[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7$chunks, '_', 2)
head(pl_per_contig_k7)


pl_per_contig_k7$chunks <- gsub('n', '', pl_per_contig_k7$chunks)
pl_per_contig_k7$seed <- gsub('s', '', pl_per_contig_k7$seed)
pl_per_contig_k7$chunk = as.numeric(pl_per_contig_k7$chunk)
head(pl_per_contig_k7)

factor(pl_per_contig_k7$chunks)
pl_per_contig_k7$chunks <- factor(pl_per_contig_k7$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                      '10', '20', '30', '40', '50', '100', ''))

factor(pl_per_contig_k7$greedy)


# Read chunked training

pl_per_contig_k7best=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v25_8k_s28_clade156_chunks_v1/all_best.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7best$V1 <- gsub('.pl_err', '', pl_per_contig_k7best$V1)
pl_per_contig_k7best[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7best$V1, '_', 2)
pl_per_contig_k7best[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7best$greedy, '_', 2)
pl_per_contig_k7best[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7best$chunks, '_', 2)

pl_per_contig_k7best$chunks <- gsub('n', '', pl_per_contig_k7best$chunks)
pl_per_contig_k7best$seed <- gsub('s', '', pl_per_contig_k7best$seed)
pl_per_contig_k7best$chunk = as.numeric(pl_per_contig_k7best$chunk)
pl_per_contig_k7best$chunks <- factor(pl_per_contig_k7best$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                  '10', '20', '30', '40', '50', '100', ''))


pl_per_contig_k7consec=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v25_8k_s28_clade156_chunks_v1/all_consec.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7consec$V1 <- gsub('.pl_err', '', pl_per_contig_k7consec$V1)
pl_per_contig_k7consec[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7consec$V1, '_', 2)
pl_per_contig_k7consec[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7consec$greedy, '_', 2)
pl_per_contig_k7consec[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7consec$chunks, '_', 2)

pl_per_contig_k7consec$chunks <- gsub('n', '', pl_per_contig_k7consec$chunks)
pl_per_contig_k7consec$seed <- gsub('s', '', pl_per_contig_k7consec$seed)
pl_per_contig_k7consec$chunk = as.numeric(pl_per_contig_k7consec$chunk)
pl_per_contig_k7consec$chunks <- factor(pl_per_contig_k7consec$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                              '10', '20', '30', '40', '50', '100', ''))

pl_per_contig_k7last=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v25_8k_s28_clade156_chunks_v1/all_last.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7last$V1 <- gsub('.pl_err', '', pl_per_contig_k7last$V1)
pl_per_contig_k7last[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7last$V1, '_', 2)
pl_per_contig_k7last[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7last$greedy, '_', 2)
pl_per_contig_k7last[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7last$chunks, '_', 2)

pl_per_contig_k7last$chunks <- gsub('n', '', pl_per_contig_k7last$chunks)
pl_per_contig_k7last$seed <- gsub('s', '', pl_per_contig_k7last$seed)
pl_per_contig_k7last$chunk = as.numeric(pl_per_contig_k7last$chunk)
pl_per_contig_k7last$chunks <- factor(pl_per_contig_k7last$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                              '10', '20', '30', '40', '50', '100', ''))


# Reading 16k data

pl_per_contig_k7last_16k=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v25_16k_s28_clade156_chunks_v1/all_last.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7last_16k$V1 <- gsub('.pl_err', '', pl_per_contig_k7last_16k$V1)
pl_per_contig_k7last_16k[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7last_16k$V1, '_', 2)
pl_per_contig_k7last_16k[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7last_16k$greedy, '_', 2)
pl_per_contig_k7last_16k[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7last_16k$chunks, '_', 2)

pl_per_contig_k7last_16k$chunks <- gsub('n', '', pl_per_contig_k7last_16k$chunks)
pl_per_contig_k7last_16k$seed <- gsub('s', '', pl_per_contig_k7last_16k$seed)
pl_per_contig_k7last_16k$chunk = as.numeric(pl_per_contig_k7last_16k$chunk)
pl_per_contig_k7last_16k$chunks <- factor(pl_per_contig_k7last_16k$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                              '10', '20', '30', '40', '50', '100', ''))


pl_per_contig_k7consec_16k=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v25_16k_s28_clade156_chunks_v1/all_consec.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7consec_16k$V1 <- gsub('.pl_err', '', pl_per_contig_k7consec_16k$V1)
pl_per_contig_k7consec_16k[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7consec_16k$V1, '_', 2)
pl_per_contig_k7consec_16k[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7consec_16k$greedy, '_', 2)
pl_per_contig_k7consec_16k[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7consec_16k$chunks, '_', 2)

pl_per_contig_k7consec_16k$chunks <- gsub('n', '', pl_per_contig_k7consec_16k$chunks)
pl_per_contig_k7consec_16k$seed <- gsub('s', '', pl_per_contig_k7consec_16k$seed)
pl_per_contig_k7consec_16k$chunk = as.numeric(pl_per_contig_k7consec_16k$chunk)
pl_per_contig_k7consec_16k$chunks <- factor(pl_per_contig_k7consec_16k$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                  '10', '20', '30', '40', '50', '100', ''))


pl_per_contig_k7best_16k=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v25_16k_s28_clade156_chunks_v1/all_best.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7best_16k$V1 <- gsub('.pl_err', '', pl_per_contig_k7best_16k$V1)
pl_per_contig_k7best_16k[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7best_16k$V1, '_', 2)
pl_per_contig_k7best_16k[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7best_16k$greedy, '_', 2)
pl_per_contig_k7best_16k[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7best_16k$chunks, '_', 2)

pl_per_contig_k7best_16k$chunks <- gsub('n', '', pl_per_contig_k7best_16k$chunks)
pl_per_contig_k7best_16k$seed <- gsub('s', '', pl_per_contig_k7best_16k$seed)
pl_per_contig_k7best_16k$chunk = as.numeric(pl_per_contig_k7best_16k$chunk)
pl_per_contig_k7best_16k$chunks <- factor(pl_per_contig_k7best_16k$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                              '10', '20', '30', '40', '50', '100', ''))




pl_per_contig_k7['cond'] = "default"
pl_per_contig_k7best['cond'] = "best"
pl_per_contig_k7consec['cond'] = "consec"
pl_per_contig_k7last['cond'] = "last"

pl_per_contig_k7best_16k['cond'] = "best"
pl_per_contig_k7consec_16k['cond'] = "consec"
pl_per_contig_k7last_16k['cond'] = "last"


pl_per_contig_k7['epoch'] = "def"
pl_per_contig_k7best['epoch'] = "8k"
pl_per_contig_k7consec['epoch'] = "8k"
pl_per_contig_k7last['epoch'] = "8k"

pl_per_contig_k7best_16k['epoch'] = "16k"
pl_per_contig_k7consec_16k['epoch'] = "16k"
pl_per_contig_k7last_16k['epoch'] = "16k"

pl_per_contig_concat = rbind(pl_per_contig_k7, pl_per_contig_k7best, pl_per_contig_k7consec, pl_per_contig_k7last)


tail(pl_per_contig_concat)


ggplot(aes(x=chunks, y=V4, fill = cond),
       data=pl_per_contig_concat[pl_per_contig_concat$greedy == "ng",])+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  #stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()+
  stat_summary( alpha = 0.5,position = position_dodge(width=0.75),color="black")+
  stat_summary(alpha = 0.5,geom="bar",position = position_dodge(width=0.75),color="black"
               )+
  coord_cartesian(ylim=c(0,10))+
  scale_x_discrete(name="Chunk length (KB)", label = c("10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  guides(col= guide_legend(title= "k"))+
  scale_fill_brewer(name="",palette = "Spectral")
  #guides(fill=guide_legend(title="k"))
  theme(legend.position = c(0.9,0.75), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )
#scale_colour_brewer(palette = "Dark2", name="")

#ggsave("incomplete_genomes_stopping_cond.pdf",width=5,height = 4)



pl_per_contig_concat2 = rbind(pl_per_contig_k7, pl_per_contig_k7best, pl_per_contig_k7consec, pl_per_contig_k7last, pl_per_contig_k7best_16k, pl_per_contig_k7consec_16k, pl_per_contig_k7last_16k)

ggplot(aes(x=chunks, y=V4, color = epoch),
       data=pl_per_contig_concat2[pl_per_contig_concat2$greedy == "ng" & pl_per_contig_concat2$cond %in% c("best", "default"),])+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  #stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()+
  stat_summary(fun.data = "mean_cl_boot",size = 0.8, alpha = 0.5)+
  coord_cartesian(ylim=c(0,10))+
  scale_x_discrete(name="Chunk length (KB)", label = c("10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  guides(col= guide_legend(title= "k"))+
  #scale_fill_discrete(name = "New Legend Title")
  #guides(fill=guide_legend(title="k"))
  theme(legend.position = c(0.9,0.75), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )




##################################################################

##################################################################

##################################################################
######################### Comparing different stopping criteria #########################


# read k7 with default training

pl_per_contig_k7=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_complete_q_v2/contigs/all_subset.pl_err',sep=" ",header=FALSE)
head(pl_per_contig_k7)

pl_per_contig_k7$V1 <- gsub('.pl_err', '', pl_per_contig_k7$V1)
pl_per_contig_k7[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7$V1, '_', 2)
pl_per_contig_k7[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7$greedy, '_', 2)
pl_per_contig_k7[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7$chunks, '_', 2)
head(pl_per_contig_k7)


pl_per_contig_k7$chunks <- gsub('n', '', pl_per_contig_k7$chunks)
pl_per_contig_k7$seed <- gsub('s', '', pl_per_contig_k7$seed)
pl_per_contig_k7$chunk = as.numeric(pl_per_contig_k7$chunk)
head(pl_per_contig_k7)

factor(pl_per_contig_k7$chunks)
pl_per_contig_k7$chunks <- factor(pl_per_contig_k7$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                      '10', '20', '30', '40', '50', '100', ''))

factor(pl_per_contig_k7$greedy)


# Read chunked training

pl_per_contig_k7best=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v25_8k_s28_clade156_chunks_v1/all_best.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7best$V1 <- gsub('.pl_err', '', pl_per_contig_k7best$V1)
pl_per_contig_k7best[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7best$V1, '_', 2)
pl_per_contig_k7best[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7best$greedy, '_', 2)
pl_per_contig_k7best[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7best$chunks, '_', 2)

pl_per_contig_k7best$chunks <- gsub('n', '', pl_per_contig_k7best$chunks)
pl_per_contig_k7best$seed <- gsub('s', '', pl_per_contig_k7best$seed)
pl_per_contig_k7best$chunk = as.numeric(pl_per_contig_k7best$chunk)
pl_per_contig_k7best$chunks <- factor(pl_per_contig_k7best$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                              '10', '20', '30', '40', '50', '100', ''))


pl_per_contig_k7consec=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v25_8k_s28_clade156_chunks_v1/all_consec.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7consec$V1 <- gsub('.pl_err', '', pl_per_contig_k7consec$V1)
pl_per_contig_k7consec[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7consec$V1, '_', 2)
pl_per_contig_k7consec[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7consec$greedy, '_', 2)
pl_per_contig_k7consec[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7consec$chunks, '_', 2)

pl_per_contig_k7consec$chunks <- gsub('n', '', pl_per_contig_k7consec$chunks)
pl_per_contig_k7consec$seed <- gsub('s', '', pl_per_contig_k7consec$seed)
pl_per_contig_k7consec$chunk = as.numeric(pl_per_contig_k7consec$chunk)
pl_per_contig_k7consec$chunks <- factor(pl_per_contig_k7consec$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                  '10', '20', '30', '40', '50', '100', ''))

pl_per_contig_k7last=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v25_8k_s28_clade156_chunks_v1/all_last.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7last$V1 <- gsub('.pl_err', '', pl_per_contig_k7last$V1)
pl_per_contig_k7last[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7last$V1, '_', 2)
pl_per_contig_k7last[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7last$greedy, '_', 2)
pl_per_contig_k7last[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7last$chunks, '_', 2)

pl_per_contig_k7last$chunks <- gsub('n', '', pl_per_contig_k7last$chunks)
pl_per_contig_k7last$seed <- gsub('s', '', pl_per_contig_k7last$seed)
pl_per_contig_k7last$chunk = as.numeric(pl_per_contig_k7last$chunk)
pl_per_contig_k7last$chunks <- factor(pl_per_contig_k7last$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                              '10', '20', '30', '40', '50', '100', ''))


# Reading 2 row chunked data training v30

pl_per_contig_k7last_16k=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v31_8k_s28_clade156_chunks_v1/all_last.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7last_16k$V1 <- gsub('.pl_err', '', pl_per_contig_k7last_16k$V1)
pl_per_contig_k7last_16k[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7last_16k$V1, '_', 2)
pl_per_contig_k7last_16k[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7last_16k$greedy, '_', 2)
pl_per_contig_k7last_16k[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7last_16k$chunks, '_', 2)

pl_per_contig_k7last_16k$chunks <- gsub('n', '', pl_per_contig_k7last_16k$chunks)
pl_per_contig_k7last_16k$seed <- gsub('s', '', pl_per_contig_k7last_16k$seed)
pl_per_contig_k7last_16k$chunk = as.numeric(pl_per_contig_k7last_16k$chunk)
pl_per_contig_k7last_16k$chunks <- factor(pl_per_contig_k7last_16k$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                      '10', '20', '30', '40', '50', '100', ''))


pl_per_contig_k7consec_16k=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v31_8k_s28_clade156_chunks_v1/all_consec.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7consec_16k$V1 <- gsub('.pl_err', '', pl_per_contig_k7consec_16k$V1)
pl_per_contig_k7consec_16k[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7consec_16k$V1, '_', 2)
pl_per_contig_k7consec_16k[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7consec_16k$greedy, '_', 2)
pl_per_contig_k7consec_16k[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7consec_16k$chunks, '_', 2)

pl_per_contig_k7consec_16k$chunks <- gsub('n', '', pl_per_contig_k7consec_16k$chunks)
pl_per_contig_k7consec_16k$seed <- gsub('s', '', pl_per_contig_k7consec_16k$seed)
pl_per_contig_k7consec_16k$chunk = as.numeric(pl_per_contig_k7consec_16k$chunk)
pl_per_contig_k7consec_16k$chunks <- factor(pl_per_contig_k7consec_16k$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                          '10', '20', '30', '40', '50', '100', ''))


pl_per_contig_k7best_16k=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v31_8k_s28_clade156_chunks_v1/all_best.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7best_16k$V1 <- gsub('.pl_err', '', pl_per_contig_k7best_16k$V1)
pl_per_contig_k7best_16k[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7best_16k$V1, '_', 2)
pl_per_contig_k7best_16k[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7best_16k$greedy, '_', 2)
pl_per_contig_k7best_16k[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7best_16k$chunks, '_', 2)

pl_per_contig_k7best_16k$chunks <- gsub('n', '', pl_per_contig_k7best_16k$chunks)
pl_per_contig_k7best_16k$seed <- gsub('s', '', pl_per_contig_k7best_16k$seed)
pl_per_contig_k7best_16k$chunk = as.numeric(pl_per_contig_k7best_16k$chunk)
pl_per_contig_k7best_16k$chunks <- factor(pl_per_contig_k7best_16k$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                      '10', '20', '30', '40', '50', '100', ''))




pl_per_contig_k7['cond'] = "default"
pl_per_contig_k7best['cond'] = "best"
pl_per_contig_k7consec['cond'] = "consec"
pl_per_contig_k7last['cond'] = "last"

pl_per_contig_k7best_16k['cond'] = "best"
pl_per_contig_k7consec_16k['cond'] = "consec"
pl_per_contig_k7last_16k['cond'] = "last"


pl_per_contig_k7['epoch'] = "def"
pl_per_contig_k7best['epoch'] = "single_row"
pl_per_contig_k7consec['epoch'] = "single_row"
pl_per_contig_k7last['epoch'] = "single_row"

pl_per_contig_k7best_16k['epoch'] = "two_rows"
pl_per_contig_k7consec_16k['epoch'] = "two_rows"
pl_per_contig_k7last_16k['epoch'] = "two_rows"

pl_per_contig_concat = rbind(pl_per_contig_k7, pl_per_contig_k7best, pl_per_contig_k7consec, pl_per_contig_k7last)


tail(pl_per_contig_concat)


ggplot(aes(x=chunks, y=V4, fill = cond),
       data=pl_per_contig_concat[pl_per_contig_concat$greedy == "ng",])+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  #stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()+
  stat_summary( alpha = 0.5,position = position_dodge(width=0.75),color="black")+
  stat_summary(alpha = 0.5,geom="bar",position = position_dodge(width=0.75),color="black"
  )+
  coord_cartesian(ylim=c(0,10))+
  scale_x_discrete(name="Chunk length (KB)", label = c("10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  guides(col= guide_legend(title= "k"))+
  scale_fill_brewer(name="",palette = "Spectral")
#guides(fill=guide_legend(title="k"))
theme(legend.position = c(0.9,0.75), legend.margin=margin(0,0,0,0),
      #axis.text.x = element_text(size = 8)
)
#scale_colour_brewer(palette = "Dark2", name="")

#ggsave("incomplete_genomes_stopping_cond.pdf",width=5,height = 4)



pl_per_contig_concat2 = rbind(pl_per_contig_k7, pl_per_contig_k7best, pl_per_contig_k7consec, pl_per_contig_k7last, pl_per_contig_k7best_16k, pl_per_contig_k7consec_16k, pl_per_contig_k7last_16k)

ggplot(aes(x=chunks, y=V4, color = epoch),
       data=pl_per_contig_concat2[pl_per_contig_concat2$greedy == "ng" & pl_per_contig_concat2$cond %in% c("best", "default"),])+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  #stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()+
  stat_summary(fun.data = "mean_cl_boot",size = 0.8, alpha = 0.5)+
  coord_cartesian(ylim=c(0,10))+
  scale_x_discrete(name="Chunk length (KB)", label = c("10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  guides(col= guide_legend(title= "k"))+
  #scale_fill_discrete(name = "New Legend Title")
  #guides(fill=guide_legend(title="k"))
  theme(legend.position = c(0.9,0.75), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )



#######################################################################################
#######################################################################################
#######################################################################################


# Comparinf dif ways to pull rows together


# read k7 with default training

pl_per_contig_k7=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_complete_q_v2/contigs/all_subset.pl_err',sep=" ",header=FALSE)
head(pl_per_contig_k7)
nrow(pl_per_contig_k7)

pl_per_contig_k7$V1 <- gsub('.pl_err', '', pl_per_contig_k7$V1)
pl_per_contig_k7[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7$V1, '_', 2)
pl_per_contig_k7[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7$greedy, '_', 2)
pl_per_contig_k7[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7$chunks, '_', 2)
head(pl_per_contig_k7)


pl_per_contig_k7$chunks <- gsub('n', '', pl_per_contig_k7$chunks)
pl_per_contig_k7$seed <- gsub('s', '', pl_per_contig_k7$seed)
pl_per_contig_k7$chunk = as.numeric(pl_per_contig_k7$chunk)
head(pl_per_contig_k7)

factor(pl_per_contig_k7$chunks)
pl_per_contig_k7$chunks <- factor(pl_per_contig_k7$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                      '10', '20', '30', '40', '50', '100', ''))

factor(pl_per_contig_k7$greedy)


# Read chunked training

pl_per_contig_k7best=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v25_8k_s28_clade156_chunks_v1/all_best.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7best$V1 <- gsub('.pl_err', '', pl_per_contig_k7best$V1)
pl_per_contig_k7best[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7best$V1, '_', 2)
pl_per_contig_k7best[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7best$greedy, '_', 2)
pl_per_contig_k7best[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7best$chunks, '_', 2)

pl_per_contig_k7best$chunks <- gsub('n', '', pl_per_contig_k7best$chunks)
pl_per_contig_k7best$seed <- gsub('s', '', pl_per_contig_k7best$seed)
pl_per_contig_k7best$chunk = as.numeric(pl_per_contig_k7best$chunk)
pl_per_contig_k7best$chunks <- factor(pl_per_contig_k7best$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                              '10', '20', '30', '40', '50', '100', ''))


pl_per_contig_k7consec=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v25_8k_s28_clade156_chunks_v1/all_consec.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7consec$V1 <- gsub('.pl_err', '', pl_per_contig_k7consec$V1)
pl_per_contig_k7consec[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7consec$V1, '_', 2)
pl_per_contig_k7consec[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7consec$greedy, '_', 2)
pl_per_contig_k7consec[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7consec$chunks, '_', 2)

pl_per_contig_k7consec$chunks <- gsub('n', '', pl_per_contig_k7consec$chunks)
pl_per_contig_k7consec$seed <- gsub('s', '', pl_per_contig_k7consec$seed)
pl_per_contig_k7consec$chunk = as.numeric(pl_per_contig_k7consec$chunk)
pl_per_contig_k7consec$chunks <- factor(pl_per_contig_k7consec$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                  '10', '20', '30', '40', '50', '100', ''))

pl_per_contig_k7last=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v25_8k_s28_clade156_chunks_v1/all_last.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7last$V1 <- gsub('.pl_err', '', pl_per_contig_k7last$V1)
pl_per_contig_k7last[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7last$V1, '_', 2)
pl_per_contig_k7last[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7last$greedy, '_', 2)
pl_per_contig_k7last[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7last$chunks, '_', 2)

pl_per_contig_k7last$chunks <- gsub('n', '', pl_per_contig_k7last$chunks)
pl_per_contig_k7last$seed <- gsub('s', '', pl_per_contig_k7last$seed)
pl_per_contig_k7last$chunk = as.numeric(pl_per_contig_k7last$chunk)
pl_per_contig_k7last$chunks <- factor(pl_per_contig_k7last$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                              '10', '20', '30', '40', '50', '100', ''))


# Reading 2 row chunked data training v34

pl_per_contig_k7last_16k=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v34_8k_s28_clade156_chunks_v1/all_last.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7last_16k$V1 <- gsub('.pl_err', '', pl_per_contig_k7last_16k$V1)
pl_per_contig_k7last_16k[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7last_16k$V1, '_', 2)
pl_per_contig_k7last_16k[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7last_16k$greedy, '_', 2)
pl_per_contig_k7last_16k[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7last_16k$chunks, '_', 2)

pl_per_contig_k7last_16k$chunks <- gsub('n', '', pl_per_contig_k7last_16k$chunks)
pl_per_contig_k7last_16k$seed <- gsub('s', '', pl_per_contig_k7last_16k$seed)
pl_per_contig_k7last_16k$chunk = as.numeric(pl_per_contig_k7last_16k$chunk)
pl_per_contig_k7last_16k$chunks <- factor(pl_per_contig_k7last_16k$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                      '10', '20', '30', '40', '50', '100', ''))


pl_per_contig_k7consec_16k=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v34_8k_s28_clade156_chunks_v1/all_consec.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7consec_16k$V1 <- gsub('.pl_err', '', pl_per_contig_k7consec_16k$V1)
pl_per_contig_k7consec_16k[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7consec_16k$V1, '_', 2)
pl_per_contig_k7consec_16k[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7consec_16k$greedy, '_', 2)
pl_per_contig_k7consec_16k[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7consec_16k$chunks, '_', 2)

pl_per_contig_k7consec_16k$chunks <- gsub('n', '', pl_per_contig_k7consec_16k$chunks)
pl_per_contig_k7consec_16k$seed <- gsub('s', '', pl_per_contig_k7consec_16k$seed)
pl_per_contig_k7consec_16k$chunk = as.numeric(pl_per_contig_k7consec_16k$chunk)
pl_per_contig_k7consec_16k$chunks <- factor(pl_per_contig_k7consec_16k$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                          '10', '20', '30', '40', '50', '100', ''))


pl_per_contig_k7best_16k=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v34_8k_s28_clade156_chunks_v1/all_best.pl_err',sep=" ",header=FALSE)

nrow(pl_per_contig_k7best_16k)

pl_per_contig_k7best_16k$V1 <- gsub('.pl_err', '', pl_per_contig_k7best_16k$V1)
pl_per_contig_k7best_16k[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7best_16k$V1, '_', 2)
pl_per_contig_k7best_16k[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7best_16k$greedy, '_', 2)
pl_per_contig_k7best_16k[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7best_16k$chunks, '_', 2)

pl_per_contig_k7best_16k$chunks <- gsub('n', '', pl_per_contig_k7best_16k$chunks)
pl_per_contig_k7best_16k$seed <- gsub('s', '', pl_per_contig_k7best_16k$seed)
pl_per_contig_k7best_16k$chunk = as.numeric(pl_per_contig_k7best_16k$chunk)
pl_per_contig_k7best_16k$chunks <- factor(pl_per_contig_k7best_16k$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                      '10', '20', '30', '40', '50', '100', ''))




pl_per_contig_k7['cond'] = "default"
pl_per_contig_k7best['cond'] = "best"
pl_per_contig_k7consec['cond'] = "consec"
pl_per_contig_k7last['cond'] = "last"

pl_per_contig_k7best_16k['cond'] = "best"
pl_per_contig_k7consec_16k['cond'] = "consec"
pl_per_contig_k7last_16k['cond'] = "last"


pl_per_contig_k7['epoch'] = "def"
pl_per_contig_k7best['epoch'] = "single_row"
pl_per_contig_k7consec['epoch'] = "single_row"
pl_per_contig_k7last['epoch'] = "single_row"

pl_per_contig_k7best_16k['epoch'] = "two_rows"
pl_per_contig_k7consec_16k['epoch'] = "two_rows"
pl_per_contig_k7last_16k['epoch'] = "two_rows"

pl_per_contig_concat = rbind(pl_per_contig_k7, pl_per_contig_k7best, pl_per_contig_k7consec, pl_per_contig_k7last)


tail(pl_per_contig_concat)


ggplot(aes(x=chunks, y=V4, fill = cond),
       data=pl_per_contig_concat[pl_per_contig_concat$greedy == "ng",])+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  #stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()+
  stat_summary( alpha = 0.5,position = position_dodge(width=0.75),color="black")+
  stat_summary(alpha = 0.5,geom="bar",position = position_dodge(width=0.75),color="black"
  )+
  coord_cartesian(ylim=c(0,10))+
  scale_x_discrete(name="Chunk length (KB)", label = c("10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  guides(col= guide_legend(title= "k"))+
  scale_fill_brewer(name="",palette = "Spectral")
#guides(fill=guide_legend(title="k"))
theme(legend.position = c(0.9,0.75), legend.margin=margin(0,0,0,0),
      #axis.text.x = element_text(size = 8)
)
#scale_colour_brewer(palette = "Dark2", name="")

#ggsave("incomplete_genomes_stopping_cond.pdf",width=5,height = 4)



pl_per_contig_concat2 = rbind(pl_per_contig_k7, pl_per_contig_k7best, pl_per_contig_k7consec, pl_per_contig_k7last, pl_per_contig_k7best_16k, pl_per_contig_k7consec_16k, pl_per_contig_k7last_16k)

ggplot(aes(x=chunks, y=V4, color = epoch),
       data=pl_per_contig_concat2[pl_per_contig_concat2$greedy == "ng" & pl_per_contig_concat2$cond %in% c("best", "default"),])+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  #stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()+
  stat_summary(fun.data = "mean_cl_boot",size = 0.8, alpha = 0.5)+
  coord_cartesian(ylim=c(0,10))+
  scale_x_discrete(name="Chunk length (KB)", label = c("10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  guides(col= guide_legend(title= "k"))+
  #scale_fill_discrete(name = "New Legend Title")
  #guides(fill=guide_legend(title="k"))
  theme(legend.position = c(0.9,0.75), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )



#######################################################################################
#######################################################################################
#######################################################################################



#######################################################################################
#######################################################################################
#######################################################################################


# Comparinf dif ways to pull rows together + trained classiier (still clades 1, 5, and 6 only)


# read k7 with default training

pl_per_contig_k7=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_complete_q_v2/contigs/all_subset.pl_err',sep=" ",header=FALSE)
head(pl_per_contig_k7)

nrow(pl_per_contig_k7)

pl_per_contig_k7$V1 <- gsub('.pl_err', '', pl_per_contig_k7$V1)
pl_per_contig_k7[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7$V1, '_', 2)
pl_per_contig_k7[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7$greedy, '_', 2)
pl_per_contig_k7[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7$chunks, '_', 2)
head(pl_per_contig_k7)


pl_per_contig_k7$chunks <- gsub('n', '', pl_per_contig_k7$chunks)
pl_per_contig_k7$seed <- gsub('s', '', pl_per_contig_k7$seed)
pl_per_contig_k7$chunk = as.numeric(pl_per_contig_k7$chunk)
head(pl_per_contig_k7)

factor(pl_per_contig_k7$chunks)
pl_per_contig_k7$chunks <- factor(pl_per_contig_k7$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                      '10', '20', '30', '40', '50', '100', ''))

factor(pl_per_contig_k7$greedy)


# Read chunked training

pl_per_contig_k7best=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v35_8k_s28_clade156_TrainClassf_complete_q_v2/all_best.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7best$V1 <- gsub('.pl_err', '', pl_per_contig_k7best$V1)
pl_per_contig_k7best[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7best$V1, '_', 2)
pl_per_contig_k7best[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7best$greedy, '_', 2)
pl_per_contig_k7best[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7best$chunks, '_', 2)

pl_per_contig_k7best$chunks <- gsub('n', '', pl_per_contig_k7best$chunks)
pl_per_contig_k7best$seed <- gsub('s', '', pl_per_contig_k7best$seed)
pl_per_contig_k7best$chunk = as.numeric(pl_per_contig_k7best$chunk)
pl_per_contig_k7best$chunks <- factor(pl_per_contig_k7best$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                              '10', '20', '30', '40', '50', '100', ''))


pl_per_contig_k7consec=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v35_8k_s28_clade156_TrainClassf_complete_q_v2/all_consec.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7consec$V1 <- gsub('.pl_err', '', pl_per_contig_k7consec$V1)
pl_per_contig_k7consec[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7consec$V1, '_', 2)
pl_per_contig_k7consec[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7consec$greedy, '_', 2)
pl_per_contig_k7consec[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7consec$chunks, '_', 2)

pl_per_contig_k7consec$chunks <- gsub('n', '', pl_per_contig_k7consec$chunks)
pl_per_contig_k7consec$seed <- gsub('s', '', pl_per_contig_k7consec$seed)
pl_per_contig_k7consec$chunk = as.numeric(pl_per_contig_k7consec$chunk)
pl_per_contig_k7consec$chunks <- factor(pl_per_contig_k7consec$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                  '10', '20', '30', '40', '50', '100', ''))

pl_per_contig_k7last=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v35_8k_s28_clade156_TrainClassf_complete_q_v2/all_last.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7last$V1 <- gsub('.pl_err', '', pl_per_contig_k7last$V1)
pl_per_contig_k7last[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7last$V1, '_', 2)
pl_per_contig_k7last[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7last$greedy, '_', 2)
pl_per_contig_k7last[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7last$chunks, '_', 2)

pl_per_contig_k7last$chunks <- gsub('n', '', pl_per_contig_k7last$chunks)
pl_per_contig_k7last$seed <- gsub('s', '', pl_per_contig_k7last$seed)
pl_per_contig_k7last$chunk = as.numeric(pl_per_contig_k7last$chunk)
pl_per_contig_k7last$chunks <- factor(pl_per_contig_k7last$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                              '10', '20', '30', '40', '50', '100', ''))


# Reading 2 row chunked data training v34

pl_per_contig_k7last_16k=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v34_8k_s28_clade156_chunks_v1/all_last.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7last_16k$V1 <- gsub('.pl_err', '', pl_per_contig_k7last_16k$V1)
pl_per_contig_k7last_16k[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7last_16k$V1, '_', 2)
pl_per_contig_k7last_16k[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7last_16k$greedy, '_', 2)
pl_per_contig_k7last_16k[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7last_16k$chunks, '_', 2)

pl_per_contig_k7last_16k$chunks <- gsub('n', '', pl_per_contig_k7last_16k$chunks)
pl_per_contig_k7last_16k$seed <- gsub('s', '', pl_per_contig_k7last_16k$seed)
pl_per_contig_k7last_16k$chunk = as.numeric(pl_per_contig_k7last_16k$chunk)
pl_per_contig_k7last_16k$chunks <- factor(pl_per_contig_k7last_16k$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                      '10', '20', '30', '40', '50', '100', ''))


pl_per_contig_k7consec_16k=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v34_8k_s28_clade156_chunks_v1/all_consec.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7consec_16k$V1 <- gsub('.pl_err', '', pl_per_contig_k7consec_16k$V1)
pl_per_contig_k7consec_16k[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7consec_16k$V1, '_', 2)
pl_per_contig_k7consec_16k[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7consec_16k$greedy, '_', 2)
pl_per_contig_k7consec_16k[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7consec_16k$chunks, '_', 2)

pl_per_contig_k7consec_16k$chunks <- gsub('n', '', pl_per_contig_k7consec_16k$chunks)
pl_per_contig_k7consec_16k$seed <- gsub('s', '', pl_per_contig_k7consec_16k$seed)
pl_per_contig_k7consec_16k$chunk = as.numeric(pl_per_contig_k7consec_16k$chunk)
pl_per_contig_k7consec_16k$chunks <- factor(pl_per_contig_k7consec_16k$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                          '10', '20', '30', '40', '50', '100', ''))


pl_per_contig_k7best_16k=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v34_8k_s28_clade156_chunks_v1/all_best.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7best_16k$V1 <- gsub('.pl_err', '', pl_per_contig_k7best_16k$V1)
pl_per_contig_k7best_16k[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7best_16k$V1, '_', 2)
pl_per_contig_k7best_16k[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7best_16k$greedy, '_', 2)
pl_per_contig_k7best_16k[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7best_16k$chunks, '_', 2)

pl_per_contig_k7best_16k$chunks <- gsub('n', '', pl_per_contig_k7best_16k$chunks)
pl_per_contig_k7best_16k$seed <- gsub('s', '', pl_per_contig_k7best_16k$seed)
pl_per_contig_k7best_16k$chunk = as.numeric(pl_per_contig_k7best_16k$chunk)
pl_per_contig_k7best_16k$chunks <- factor(pl_per_contig_k7best_16k$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                      '10', '20', '30', '40', '50', '100', ''))




pl_per_contig_k7['cond'] = "default"
pl_per_contig_k7best['cond'] = "best"
pl_per_contig_k7consec['cond'] = "consec"
pl_per_contig_k7last['cond'] = "last"

pl_per_contig_k7best_16k['cond'] = "best"
pl_per_contig_k7consec_16k['cond'] = "consec"
pl_per_contig_k7last_16k['cond'] = "last"


pl_per_contig_k7['epoch'] = "def"
pl_per_contig_k7best['epoch'] = "two_rows_w_classif"
pl_per_contig_k7consec['epoch'] = "two_rows_w_classif"
pl_per_contig_k7last['epoch'] = "two_rows_w_classif"

pl_per_contig_k7best_16k['epoch'] = "two_rows"
pl_per_contig_k7consec_16k['epoch'] = "two_rows"
pl_per_contig_k7last_16k['epoch'] = "two_rows"

pl_per_contig_concat = rbind(pl_per_contig_k7, pl_per_contig_k7best, pl_per_contig_k7consec, pl_per_contig_k7last)


tail(pl_per_contig_concat)


ggplot(aes(x=chunks, y=V4, fill = cond),
       data=pl_per_contig_concat[pl_per_contig_concat$greedy == "ng",])+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  #stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()+
  stat_summary( alpha = 0.5,position = position_dodge(width=0.75),color="black")+
  stat_summary(alpha = 0.5,geom="bar",position = position_dodge(width=0.75),color="black"
  )+
  coord_cartesian(ylim=c(0,10))+
  scale_x_discrete(name="Chunk length (KB)", label = c("10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  guides(col= guide_legend(title= "k"))+
  scale_fill_brewer(name="",palette = "Spectral")
#guides(fill=guide_legend(title="k"))
theme(legend.position = c(0.9,0.75), legend.margin=margin(0,0,0,0),
      #axis.text.x = element_text(size = 8)
)
#scale_colour_brewer(palette = "Dark2", name="")

#ggsave("incomplete_genomes_stopping_cond.pdf",width=5,height = 4)



pl_per_contig_concat2 = rbind(pl_per_contig_k7, pl_per_contig_k7best, pl_per_contig_k7consec, pl_per_contig_k7last, pl_per_contig_k7best_16k, pl_per_contig_k7consec_16k, pl_per_contig_k7last_16k)

ggplot(aes(x=chunks, y=V4, color = epoch),
       data=pl_per_contig_concat2[pl_per_contig_concat2$greedy == "ng" & pl_per_contig_concat2$cond %in% c("best", "default"),])+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  #stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()+
  stat_summary(fun.data = "mean_cl_boot",size = 0.8, alpha = 0.5)+
  coord_cartesian(ylim=c(0,10))+
  scale_x_discrete(name="Chunk length (KB)", label = c("10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  guides(col= guide_legend(title= "k"))+
  #scale_fill_discrete(name = "New Legend Title")
  #guides(fill=guide_legend(title="k"))
  theme(legend.position = c(0.8,0.77), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )

ggsave("incomplete_genomes_two_rows_c156.pdf",width=5,height = 4)

getwd()
# Depp paper (dataset 16s and whole genome data) contact Yuyue
# Green genes 2

#######################################################################################
#######################################################################################
#######################################################################################



############# COMPARING FULL DATASET #############



# read k7 with default training

pl_per_contig_k7=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_complete_q_v2/contigs/all.pl_err',sep=" ",header=FALSE)
head(pl_per_contig_k7)
nrow(pl_per_contig_k7)

pl_per_contig_k7$V1 <- gsub('.pl_err', '', pl_per_contig_k7$V1)
pl_per_contig_k7[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7$V1, '_', 2)
pl_per_contig_k7[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7$greedy, '_', 2)
pl_per_contig_k7[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7$chunks, '_', 2)
head(pl_per_contig_k7)


pl_per_contig_k7$chunks <- gsub('n', '', pl_per_contig_k7$chunks)
pl_per_contig_k7$seed <- gsub('s', '', pl_per_contig_k7$seed)
pl_per_contig_k7$chunk = as.numeric(pl_per_contig_k7$chunk)
head(pl_per_contig_k7)

factor(pl_per_contig_k7$chunks)
pl_per_contig_k7$chunks <- factor(pl_per_contig_k7$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                      '10', '20', '30', '40', '50', '100', ''))

factor(pl_per_contig_k7$greedy)


# Read chunked training

pl_per_contig_k7best=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v35_8k_s28_clade_ALL_chunks_v1/all_best.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7best$V1 <- gsub('.pl_err', '', pl_per_contig_k7best$V1)
pl_per_contig_k7best[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7best$V1, '_', 2)
pl_per_contig_k7best[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7best$greedy, '_', 2)
pl_per_contig_k7best[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7best$chunks, '_', 2)

pl_per_contig_k7best$chunks <- gsub('n', '', pl_per_contig_k7best$chunks)
pl_per_contig_k7best$seed <- gsub('s', '', pl_per_contig_k7best$seed)
pl_per_contig_k7best$chunk = as.numeric(pl_per_contig_k7best$chunk)
pl_per_contig_k7best$chunks <- factor(pl_per_contig_k7best$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                              '10', '20', '30', '40', '50', '100', ''))


pl_per_contig_k7consec=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v35_8k_s28_clade_ALL_chunks_v1/all_consec.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7consec$V1 <- gsub('.pl_err', '', pl_per_contig_k7consec$V1)
pl_per_contig_k7consec[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7consec$V1, '_', 2)
pl_per_contig_k7consec[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7consec$greedy, '_', 2)
pl_per_contig_k7consec[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7consec$chunks, '_', 2)

pl_per_contig_k7consec$chunks <- gsub('n', '', pl_per_contig_k7consec$chunks)
pl_per_contig_k7consec$seed <- gsub('s', '', pl_per_contig_k7consec$seed)
pl_per_contig_k7consec$chunk = as.numeric(pl_per_contig_k7consec$chunk)
pl_per_contig_k7consec$chunks <- factor(pl_per_contig_k7consec$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                  '10', '20', '30', '40', '50', '100', ''))

pl_per_contig_k7last=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v35_8k_s28_clade_ALL_chunks_v1/all_last.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7last$V1 <- gsub('.pl_err', '', pl_per_contig_k7last$V1)
pl_per_contig_k7last[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7last$V1, '_', 2)
pl_per_contig_k7last[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7last$greedy, '_', 2)
pl_per_contig_k7last[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7last$chunks, '_', 2)

pl_per_contig_k7last$chunks <- gsub('n', '', pl_per_contig_k7last$chunks)
pl_per_contig_k7last$seed <- gsub('s', '', pl_per_contig_k7last$seed)
pl_per_contig_k7last$chunk = as.numeric(pl_per_contig_k7last$chunk)
pl_per_contig_k7last$chunks <- factor(pl_per_contig_k7last$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                              '10', '20', '30', '40', '50', '100', ''))


# Reading 2 row chunked data training v34

pl_per_contig_k7last_16k=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v34_8k_s28_clade156_chunks_v1/all_last.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7last_16k$V1 <- gsub('.pl_err', '', pl_per_contig_k7last_16k$V1)
pl_per_contig_k7last_16k[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7last_16k$V1, '_', 2)
pl_per_contig_k7last_16k[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7last_16k$greedy, '_', 2)
pl_per_contig_k7last_16k[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7last_16k$chunks, '_', 2)

pl_per_contig_k7last_16k$chunks <- gsub('n', '', pl_per_contig_k7last_16k$chunks)
pl_per_contig_k7last_16k$seed <- gsub('s', '', pl_per_contig_k7last_16k$seed)
pl_per_contig_k7last_16k$chunk = as.numeric(pl_per_contig_k7last_16k$chunk)
pl_per_contig_k7last_16k$chunks <- factor(pl_per_contig_k7last_16k$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                      '10', '20', '30', '40', '50', '100', ''))


pl_per_contig_k7consec_16k=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v34_8k_s28_clade156_chunks_v1/all_consec.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7consec_16k$V1 <- gsub('.pl_err', '', pl_per_contig_k7consec_16k$V1)
pl_per_contig_k7consec_16k[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7consec_16k$V1, '_', 2)
pl_per_contig_k7consec_16k[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7consec_16k$greedy, '_', 2)
pl_per_contig_k7consec_16k[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7consec_16k$chunks, '_', 2)

pl_per_contig_k7consec_16k$chunks <- gsub('n', '', pl_per_contig_k7consec_16k$chunks)
pl_per_contig_k7consec_16k$seed <- gsub('s', '', pl_per_contig_k7consec_16k$seed)
pl_per_contig_k7consec_16k$chunk = as.numeric(pl_per_contig_k7consec_16k$chunk)
pl_per_contig_k7consec_16k$chunks <- factor(pl_per_contig_k7consec_16k$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                          '10', '20', '30', '40', '50', '100', ''))


pl_per_contig_k7best_16k=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v34_8k_s28_clade156_chunks_v1/all_best.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7best_16k$V1 <- gsub('.pl_err', '', pl_per_contig_k7best_16k$V1)
pl_per_contig_k7best_16k[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7best_16k$V1, '_', 2)
pl_per_contig_k7best_16k[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7best_16k$greedy, '_', 2)
pl_per_contig_k7best_16k[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7best_16k$chunks, '_', 2)

pl_per_contig_k7best_16k$chunks <- gsub('n', '', pl_per_contig_k7best_16k$chunks)
pl_per_contig_k7best_16k$seed <- gsub('s', '', pl_per_contig_k7best_16k$seed)
pl_per_contig_k7best_16k$chunk = as.numeric(pl_per_contig_k7best_16k$chunk)
pl_per_contig_k7best_16k$chunks <- factor(pl_per_contig_k7best_16k$chunks, levels = c('1', '2', '3', '4', '5', 
                                                                                      '10', '20', '30', '40', '50', '100', ''))




pl_per_contig_k7['cond'] = "default"
pl_per_contig_k7best['cond'] = "best"
pl_per_contig_k7consec['cond'] = "consec"
pl_per_contig_k7last['cond'] = "last"

pl_per_contig_k7best_16k['cond'] = "best"
pl_per_contig_k7consec_16k['cond'] = "consec"
pl_per_contig_k7last_16k['cond'] = "last"


pl_per_contig_k7['epoch'] = "def"
pl_per_contig_k7best['epoch'] = "two_rows_w_classif"
pl_per_contig_k7consec['epoch'] = "two_rows_w_classif"
pl_per_contig_k7last['epoch'] = "two_rows_w_classif"

pl_per_contig_k7best_16k['epoch'] = "two_rows"

pl_per_contig_k7consec_16k['epoch'] = "two_rows"
pl_per_contig_k7last_16k['epoch'] = "two_rows"

pl_per_contig_concat = rbind(pl_per_contig_k7, pl_per_contig_k7best, pl_per_contig_k7consec, pl_per_contig_k7last)


tail(pl_per_contig_concat)


ggplot(aes(x=chunks, y=V4, fill = cond),
       data=pl_per_contig_concat[pl_per_contig_concat$greedy == "ng",])+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  #stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()+
  stat_summary( alpha = 0.5,position = position_dodge(width=0.75),color="black")+
  stat_summary(alpha = 0.5,geom="bar",position = position_dodge(width=0.75),color="black"
  )+
  coord_cartesian(ylim=c(0,10))+
  scale_x_discrete(name="Chunk length (KB)", label = c("10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  guides(col= guide_legend(title= "k"))+
  scale_fill_brewer(name="",palette = "Spectral")
#guides(fill=guide_legend(title="k"))
theme(legend.position = c(0.9,0.75), legend.margin=margin(0,0,0,0),
      #axis.text.x = element_text(size = 8)
)
#scale_colour_brewer(palette = "Dark2", name="")

#ggsave("incomplete_genomes_stopping_cond.pdf",width=5,height = 4)



pl_per_contig_concat2 = rbind(pl_per_contig_k7, pl_per_contig_k7best, pl_per_contig_k7consec, pl_per_contig_k7last) #, pl_per_contig_k7best_16k, pl_per_contig_k7consec_16k, pl_per_contig_k7last_16k)

ggplot(aes(x=chunks, y=V4, color = epoch),
       data=pl_per_contig_concat2[pl_per_contig_concat2$greedy == "ng" & pl_per_contig_concat2$cond %in% c("best", "default"),])+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  #stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()+
  stat_summary(fun.data = "mean_cl_boot",size = 0.8, alpha = 0.5)+
  coord_cartesian(ylim=c(0,10))+
  scale_x_discrete(name="Chunk length (KB)", label = c("10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  guides(col= guide_legend(title= "k"))+
  #scale_fill_discrete(name = "New Legend Title")
  #guides(fill=guide_legend(title="k"))
  theme(legend.position = c(0.8,0.77), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )



pl_per_contig_concat2



#############################################################################

# read k7 with default training + RETRAINED CHUNKED MODEL


# Had to exclude clade 14 (excluded species G000018945 14 and G000237805 14)

#pl_per_contig_k7a=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_complete_q_v2/contigs/all_noclade14.pl_err',sep=" ",header=FALSE)
pl_per_contig_k7a=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_complete_q_v2/contigs/all.pl_err',sep=" ",header=FALSE)
head(pl_per_contig_k7a)
nrow(pl_per_contig_k7a)


#pl_per_contig_k7b=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_complete_q_v3/all_best_noclade14.pl_err',sep=" ",header=FALSE)
pl_per_contig_k7b=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_complete_q_v3/all_best.pl_err',sep=" ",header=FALSE)
head(pl_per_contig_k7b)
nrow(pl_per_contig_k7b)

pl_per_contig_k7 = rbind(pl_per_contig_k7a,pl_per_contig_k7b)
nrow(pl_per_contig_k7)

pl_per_contig_k7$V1 <- gsub('.pl_err', '', pl_per_contig_k7$V1)
pl_per_contig_k7[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7$V1, '_', 2)
pl_per_contig_k7[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7$greedy, '_', 2)
pl_per_contig_k7[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7$chunks, '_', 2)
head(pl_per_contig_k7)
tail(pl_per_contig_k7)

pl_per_contig_k7$chunks <- gsub('n', '', pl_per_contig_k7$chunks)
pl_per_contig_k7$seed <- gsub('s', '', pl_per_contig_k7$seed)
pl_per_contig_k7$chunk = as.numeric(pl_per_contig_k7$chunk)
head(pl_per_contig_k7)

factor(pl_per_contig_k7$chunks)
pl_per_contig_k7$chunks <- factor(pl_per_contig_k7$chunks, levels = c('0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', 
                                                                      '1', '2', '3', '4', '5', '10', '20', '30', '40', '50', '100', ''))

factor(pl_per_contig_k7$greedy)




#pl_per_contig_k7best_a=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v35_8k_s28_clade_ALL_chunks_v1/all_best_noclade14.pl_err',sep=" ",header=FALSE)
#pl_per_contig_k7best_b=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v35_8k_s28_clade_All_chunks_complete_q_v3/all_best_noclade14.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7best_a=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_clade_All_chunks_complete_q_v2/all_best.pl_err',sep=" ",header=FALSE)
pl_per_contig_k7best_b=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_clade_All_chunks_complete_q_v3/all_best.pl_err',sep=" ",header=FALSE)

pl_per_contig_k7best <- rbind(pl_per_contig_k7best_a, pl_per_contig_k7best_b)

pl_per_contig_k7best$V1 <- gsub('.pl_err', '', pl_per_contig_k7best$V1)
pl_per_contig_k7best[c('genome', 'greedy')] <- str_split_fixed(pl_per_contig_k7best$V1, '_', 2)
pl_per_contig_k7best[c('greedy', 'chunks')] <- str_split_fixed(pl_per_contig_k7best$greedy, '_', 2)
pl_per_contig_k7best[c('chunks', 'seed')] <- str_split_fixed(pl_per_contig_k7best$chunks, '_', 2)

head(pl_per_contig_k7best)


pl_per_contig_k7best$chunks <- gsub('n', '', pl_per_contig_k7best$chunks)
pl_per_contig_k7best$seed <- gsub('s', '', pl_per_contig_k7best$seed)
pl_per_contig_k7best$chunk = as.numeric(pl_per_contig_k7best$chunk)
pl_per_contig_k7best$chunks <- factor(pl_per_contig_k7best$chunks, levels = c('0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', 
                                                                              '1', '2', '3', '4', '5', '10', '20', '30', '40', '50', '100', ''))



pl_per_contig_k7['cond'] = "default"
pl_per_contig_k7best['cond'] = "best"




pl_per_contig_k7['epoch'] = "def"
pl_per_contig_k7best['epoch'] = "two_rows_w_classif"


pl_per_contig_concat2 = rbind(pl_per_contig_k7, pl_per_contig_k7best)
tail(pl_per_contig_concat2)


ggplot(aes(x=chunks, y=V4, color = epoch),
       data=pl_per_contig_concat2[pl_per_contig_concat2$greedy == "ng",])+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  #stat_summary(geom="point")+
  #geom_violin()+
  geom_hline(yintercept = 3.0, col = "grey", linetype='dotted')+
  theme_classic()+
  stat_summary(fun.data = "mean_cl_boot",size = 0.8, alpha = 0.5)+
  coord_cartesian(ylim=c(0,24))+
  scale_x_discrete(name="Chunk length (KB)", label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  guides(col= guide_legend(title= "k"))+
  #scale_fill_discrete(name = "New Legend Title")
  #guides(fill=guide_legend(title="k"))
  theme(legend.position = c(0.8,0.77), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )



head(pl_per_contig_concat2[pl_per_contig_concat2$greedy == "ng",])


#############################################################################
loss = read.csv("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v24_8k_s28_clade156_chunks_v1/epochs.txt",se=" ",h=F)
lr = read.csv("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v24_8k_s28_clade156_chunks_v1/LR.txt",se=" ",h=F)
loss$clade = c(rep("a",8000),rep("b",8000),rep("c",8000))
lr$clade = c(rep("a",8000),rep("b",8000),rep("c",8000))
loss = merge(loss,lr,by=c("V1","clade"))
nrow(loss)
names(loss)[3]="loss"
names(loss)[5]="LR"
head(loss)

ggplot(aes(x=V1,y=loss),data=loss[loss$loss>0.05,])+geom_point(alpha=0.2)+geom_smooth()+
  facet_wrap(.~clade,scales = "free_y")+
  theme_classic()+
  scale_y_log10()+
  geom_line(aes(y=LR*1000),color="red")

summary(loss[loss$clade=="b" & loss$V1>7500,"loss"])

ggplot(aes(x=V1,y=ifelse(loss>.1,1,0)),data=loss)+geom_smooth()+
  facet_wrap(.~clade,scales = "free_y")+
  theme_classic()+
  #scale_y_log10()+
  geom_line(aes(y=LR*1000),color="red")

###################################################### ################## ##################
###################################################### ################## ##################
###################################################### ################## ##################








data.frame(pl_per_contig[pl_per_contig$kmer_sums <=40000 & pl_per_contig$kmer_sums >=20000,])[,"pl_err_delta"]

res<-quantile(data.frame(pl_per_contig[pl_per_contig$kmer_sums <=640000 & pl_per_contig$kmer_sums >320000,])[,"pl_err_delta"],probs=c(0.2,0.8)) 
res

nrow(pl_per_contig[pl_per_contig$kmer_sums >+1280000,])
nrow(pl_per_contig)

ggplot(aes(x=cut(kmer_sums/10^3,c(2, 4, 8,16,32,64,128,256,1000)*10), y=pl_err_delta),
       data=pl_per_contig)+
  #geom_violin()+
  stat_summary()+
  stat_summary(geom="line",color="blue",aes(group="1"))+
  stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()+
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
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
  #stat_summary(aes(fill=
  #                   ifelse(cond =="true","Hypothetical",ifelse(cond == "combined","Default","Explore"))),geom="bar",color="black")+
  #stat_summary()+
  #scale_fill_brewer(name="")+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_abline(color="red")+
  #theme_classic()+
  #scale_fill_discrete(name="Condition") +
  #theme_bw()+
  #xlab("Condition")+
  #ylab("Placement error")+
  #scale_x_discrete(labels = c("Shorter\nkmer", "Extra\nlayer", "Global",
  #                            #  "Global", 
  #                            "Local", "Default\n(combined)", "True\nclade"),
  #                 name="")+
  #theme(axis.text.x = element_text(angle = 10, vjust = 0.7, hjust=0.5),
  #      legend.position = c(.67,.9),legend.direction = "horizontal")
#theme(axis.line = element_line(colour = "black"),
#      panel.grid.major = element_blank(),
#      panel.grid.minor = element_blank(),
#      panel.border = element_blank(),
#      panel.background = element_blank()) 

ggsave("wol_dev_queries_pl_perquery.pdf",width=5,height = 4)


#For ppt
pl_per_contig
pl_per_contig$cond2 <- ""
pl_per_contig$cond2 <- ifelse(pl_per_contig$kmer_sums <=320000, "Short", ifelse((pl_per_contig$kmer_sums >320000 & pl_per_contig$kmer_sums <=1280000 ), "Optimal", "Chimeric"))
pl_per_contig$cond2 <- factor(pl_per_contig$cond2, levels = c("Short", "Optimal", "Chimeric"))
nrow(pl_per_contig)

unique(pl_per_contig$cond2)
ggplot(aes(x=cut(kmer_sums/10^3,c(2, 4, 8,16,32,64,128,256,1000)*10), y=pl_err_delta, fill = cond2),
       data=pl_per_contig)+
  #geom_violin()+
  stat_summary(geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  stat_summary(geom="line",color="blue",aes(group="1"))+
  stat_summary(alpha = 0.9)+
  theme_classic()+
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="Contig length (KB)", label = c("(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-10000]"))  +
  #scale_colour_brewer(palette = "Set1", name="", labels = c("Test", "Train"))+
  #scale_fill_brewer(palette = "Accent", name="")+
  #scale_fill_brewer(palette = "Set2", name="")+
  scale_fill_manual(values = c(  "#d1e5f0", "#e6f5d0", "#fddbc7"), name="")+
  #scale_fill_manual(values = c(  "#d1e5f0", "#f7f7f7", "#fddbc7"), name="")+
  #scale_fill_manual(values = c(  "#E7298A", "#66A61E", "#A6761D"), name="")+
  #scale_fill_manual(values = c(  "#BEAED4", "#7FC97F", "#FDC086"), name="")+
  ylab("Placement error")+
  theme(axis.text.x = element_text(angle = 15, vjust = 0.7, hjust=0.5))+
  theme(legend.position = c(0.9,0.85), legend.margin=margin(0,0,0,0),legend.direction="vertical", legend.title=element_blank(),
        #axis.text.x = element_text(size = 8)
        )+
  guides(fill = guide_legend(reverse = FALSE))



ggsave("wol_dev_queries_pl_percontig_v2.pdf",width=5,height = 4)






ggplot(aes(x=cut(kmer_sums/10^3,c(2, 4, 8,16,32,64,128,256,1000)*10), y=pl_err_delta, fill = cond2),
       data=pl_per_contig)+
  #geom_violin()+
  #geom_jitter(alpha=0.3)+
  stat_summary(geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  stat_summary(geom="line",color="blue",aes(group="1"))+
  stat_summary(alpha = 0.9)+
  theme_classic()+
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="Contig length (KB)", label = c("(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  #scale_colour_brewer(palette = "Set1", name="", labels = c("Test", "Train"))+
  #scale_fill_brewer(palette = "Accent", name="")+
  #scale_fill_brewer(palette = "Set2", name="")+
  scale_fill_manual(values = c(  "#d1e5f0", "#e6f5d0", "#fddbc7"), name="")+
  #scale_fill_manual(values = c(  "#d1e5f0", "#f7f7f7", "#fddbc7"), name="")+
  #scale_fill_manual(values = c(  "#E7298A", "#66A61E", "#A6761D"), name="")+
  #scale_fill_manual(values = c(  "#BEAED4", "#7FC97F", "#FDC086"), name="")+
  ylab("Placement error")+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5))+
  theme(legend.position = c(0.9,0.85), legend.margin=margin(0,0,0,0),legend.direction="vertical", legend.title=element_blank(),
        #axis.text.x = element_text(size = 8)
  )+
  guides(fill = guide_legend(reverse = FALSE))
ggsave("wol_dev_queries_pl_percontig_v2_ppt.pdf",width=6.2,height = 4)



#pl_err=read.csv('pl_results_dev_queries.txt',sep=" ",header=TRUE)
pl_err=read.csv('pl_results_dev_queries_ROUND2.txt',sep=" ",header=TRUE)

head(pl_err)
mean(data.frame(pl_err[pl_err$cond == '6kC',])[,"pl_error"])
mean(data.frame(pl_err[pl_err$cond == 'moreLayers',])[,"pl_error"])
mean(data.frame(pl_err[pl_err$cond == 'main',])[,"pl_error"])
mean(data.frame(pl_err[pl_err$cond == 'global',])[,"pl_error"])
mean(data.frame(pl_err[pl_err$cond == 'local',])[,"pl_error"])
mean(data.frame(pl_err[pl_err$cond == 'combined',])[,"pl_error"])


#fin_df = rbind(pr_A1_grouped, pr_A01_mean)
#fin_df2 = cbind(pr_A1_grouped, pr_A01_mean)
#total <- merge(pr_A1_grouped, pr_A01_median, by = c("V1"))

# Development queries, multiple conditions

pl_err$cond <- factor(pl_err$cond, levels=c("6kC", "moreLayers", "main", "global", "local", "combined", "true" ))

ggplot(aes(x=cond, y=pl_error),
       data=pl_err[!pl_err$cond %in% c("main"),])+
       #data=pl_err)+
       #data=pl_err)+
  #geom_plot(alpha=0.6)+
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
  theme_classic()+
  #xlab("Condition")+
  ylab("Placement error")+
  scale_x_discrete(labels = c("Shorter\nkmer", "Extra\nlayer", "Global",
                            #  "Global", 
                              "Local", "Default\n(combined)", "True\nclade"),
                   name="")+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5), axis.title.x = element_blank(),
        legend.position = c(.82,.88),legend.direction = "vertical", legend.margin=margin(t = -0.5, unit='cm'))
  #theme(axis.line = element_line(colour = "black"),
  #      panel.grid.major = element_blank(),
  #      panel.grid.minor = element_blank(),
  #      panel.border = element_blank(),
  #      panel.background = element_blank()) 

ggsave("dev_queries_multiple_conditions.pdf",width=3.8,height = 4)


# Plot for ppt
unique(pl_err$cond)
ggplot(aes(x=cond, y=pl_error),
       data=pl_err[!pl_err$cond %in% c("global", "moreLayers", "6kC"),])+
  stat_summary(aes(fill=
                     ifelse(cond =="true","Hypothetical",ifelse(cond == "combined","Default","Explore"))),geom="bar",color="black")+
  stat_summary()+
  scale_fill_brewer(name="")+
  #scale_fill_discrete(name="Condition") +
  theme_classic()+
  #xlab("Condition")+
  ylab("Placement error")+
  scale_x_discrete(labels = c("Global",
                             #  "Global", 
                              "Local", "Default\n(combined)", "True\nclade"),
                   name="")+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5), axis.title.x = element_blank(),
        legend.position = c(.88,.90),legend.direction = "vertical", legend.margin=margin(t = -0.5, unit='cm'))


ggsave("dev_queries_multiple_conditions_for_ppt.pdf",width=5,height = 4)



err_epoch=read.csv('error_per_epoch.txt',sep=" ",header=TRUE)

head(err_epoch)

ggplot(aes(x=epoch, y=error, color = cond),
       data=err_epoch[err_epoch$epoch< 8000,])+
  geom_point(alpha = 0.9, size = 0.4)+
  geom_smooth(se=F,  size=0.5, alpha = 0.9, method = "lm",formula = y ~ poly(x, 12))+
  scale_colour_brewer(palette = "Dark2", name="", labels = c("Test", "Train"))+
  theme_classic()+
  xlab("Epoch")+
  ylab("Error")+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = c(.87,.8),legend.direction = "vertical",
        legend.title=element_blank())
  
  #scale_color_manual(values = "Dark2")
  
 
ggsave("dev_queries_err_epoch.pdf",width=5,height = 4)





##########################################################################  
##########################################################################  
##########################################################################  
##########################################################################  
# OBSOLETE

err_vs_prob=read.csv('/Users/nora/Documents/ml_metagenomics/paper_fig1/err_vs_probability.txt',sep=" ",header=TRUE)
nrow(err_vs_prob)

ggplot(aes(x=top_p, color = cond),
       data=err_vs_prob[err_vs_prob$top_p < 1.0,])+
  stat_ecdf()
  geom_step(aes(x=top_p, color = cond), stat="ecdf")
  #geom_plot(alpha=0.6)+
  #geom_bar(position = "dodge",stat = "summary",fun = "mean")+
  #geom_boxplot(alpha=0.6)+
  #stat_summary()+
  #geom_point()+
  scale_x_continuous()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_abline(color="red")+
  #theme_classic()+
  #scale_fill_discrete(name="Condition") +
  theme_bw()
  xlab("Condition")+
  ylab("Placement error")+
  scale_x_discrete(labels = c("Shorter kmer", "Extra layer", "Main", "Global", "Local", "Combined", "True"))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7, hjust=0.5))
#theme(axis.line = element_line(colour = "black"),
#      panel.grid.major = element_blank(),
#      panel.grid.minor = element_blank(),
#      panel.border = element_blank(),
#      panel.background = element_blank()) 

ggsave("my_plot.pdf",width=5,height = 4)


ggplot(aes(y=top_p, x=cut(pl_error,c(-1,0,1,2,4,8,16,40)), color = cond),
       data=err_vs_prob[,])+
  stat_summary()
#geom_plot(alpha=0.6)+
#geom_bar(position = "dodge",stat = "summary",fun = "mean")+
#geom_boxplot(alpha=0.6)+
#stat_summary()+
#geom_point()+
scale_x_continuous()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_abline(color="red")+
  #theme_classic()+
  #scale_fill_discrete(name="Condition") +
  theme_bw()
xlab("Condition")+
  ylab("Placement error")+
  scale_x_discrete(labels = c("Shorter kmer", "Extra layer", "Main", "Global", "Local", "Combined", "True"))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7, hjust=0.5))







ggplot(aes(x=clade, y=coef, color = method),
  data=corr_stats[corr_stats$corr == "pearson",])+
  #geom_bar(position = "dodge",stat = "identity")+
  geom_point(alpha = 0.9)+
  facet_wrap(~ portion)+
  #geom_boxplot(alpha=0.6)+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  #geom_abline(color="red")+
  theme_classic()
#scale_y_continuous(limits = c(0.0000, 0.001))





ggplot(aes(x=coef,  y=pl_error, color = method),
       data=corr_stats[corr_stats$corr == "pearson" & corr_stats$portion == "queries",])+
  #geom_bar(position = "dodge",stat = "identity")+
  geom_point(alpha = 0.9)+
  #facet_wrap(~ method)+
  theme_classic()+
  geom_smooth(method = lm, se = FALSE)+
  #theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  theme(axis.text.x = element_text(angle = 30, vjust = 1.0, hjust=1))
  #geom_boxplot(alpha=0.6)+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  #geom_abline(color="red")+


  
#scale_y_continuous(limits = c(0.0000, 0.001))


ggplot(aes(x=true/100, y=dist/100, color = method),
       data=dist_df)+
  #geom_bar(position = "dodge",stat = "identity")+
  geom_point(alpha = 0.5, size = 0.5)+
  facet_wrap(~ portion)+
  geom_abline(linetype=2)+
  theme_classic()+
  stat_smooth(method="lm",se=FALSE,size=1.4)+
  scale_color_brewer(palette = "Dark2",name="")
  #geom_boxplot(alpha=0.6)+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  #geom_abline(color="red")+
#scale_y_continuous(limits = c(0.0000, 0.001))

######################
prob_plot=read.csv('/Users/nora/Documents/ml_metagenomics/paper_fig1/output_probabilities.txt',sep=" ",header=FALSE)

prob_plot
ggplot(aes(x=V1, y=V2, color = V3),
       data=prob_plot)+
  #geom_bar(position = "dodge",stat = "identity")+
  geom_line(alpha = 0.9)

#er = median(total$V3.x)
#er
#er = median(total$V3.y)
#er



#################################################################
#################################################################

#################################################################
#################################################################

total
myc = c(3, 3, 3, 4, 4, 4, 4)
mean(myc)
sd(myc)

N <- length(myc)
deviations <- myc - mean(myc)
s <- deviations^2 #
m <- sum(s)/N
sd <- sqrt(m)
print(sd)



f2 <-ggplot(aes(x=pr_A1_grouped,y=pr_A01_mean),data=fin_df[fin_df$V3<0.0005,])+
  geom_line(alpha=0.6)+
  theme_classic()
#scale_y_continuous(limits = c(0.0000, 0.001))
f2




#a=read.csv('kmercounts.txt',sep=" ",header=FALSE)

#c=read.csv('kmer_pl_gpu_shared_500epochs_5000cols.csv',sep=",",header=TRUE)

#data=c[c$epoch>500 & c$epoch<1000,]

# Plot error decay

#f2 <-ggplot(aes(x=epoch,y=log10(val), color=cond),data=c)+
#  geom_line(alpha=0.6)+
#  theme_classic()
#f2

#ggsave("kmer_pl_gpu_shared_1500epochs.pdf",width=5,height = 4)




d=read.csv('mer_pl_gpu_shared_model6.csv',sep=",",header=FALSE)



head(d)
df_train <- subset(d, select = c(1:3))
df_test <- subset(d, select = c(5:6))
df_test <- setNames(df_test, names(df_train))
fin_df = rbind(df_train, df_test)
fin_df


f2 <-ggplot(aes(x=V1,y=log10(V3), color=V2),data=fin_df[fin_df$V3<0.0005,])+
  geom_line(alpha=0.6)+
  theme_classic()
  #scale_y_continuous(limits = c(0.0000, 0.001))
f2

ggsave("mer_pl_gpu_shared_model6.pdf",width=5,height = 4)



######### Plot placement_error ######### 

d=read.csv('pl_error_model_v78.txt',sep=" ",header=FALSE)

#f=read.csv('pl_error_model_16.txt',sep=" ",header=FALSE)
#f=read.csv('pl_error_model_20_7kC_8000.txt',sep=" ",header=FALSE)
#f=read.csv('pl_error_model17_8000_V3_shuffled.txt',sep=" ",header=FALSE)



#f1=read.csv('pl_error_model_12_hiqual.txt',sep=" ",header=FALSE)
f1=read.csv('pl_error_model_v78.txt',sep=" ",header=FALSE)
f2=read.csv('pl_error_model_v77.txt',sep=" ",header=FALSE)




#f2=read.csv('pl_error_model_21_7kC_8000_alpha_v1.txt',sep=" ",header=FALSE)
f3=read.csv('pl_error_model_v76.txt',sep=" ",header=FALSE)
f4=read.csv('pl_error_model_v75.txt',sep=" ",header=FALSE)

f5=read.csv('pl_error_model_21_7kC_8000_lambda_v5.txt',sep=" ",header=FALSE)
#f5b=read.csv('pl_error_model_21_7kC_8000_lambda_v5_excl200.txt',sep=" ",header=FALSE)

f6=read.csv('pl_error_model_21_7kC_8000_lambda_v6.txt',sep=" ",header=FALSE)
f7=read.csv('pl_error_model_21_7kC_8000_lambda_v7.txt',sep=" ",header=FALSE)
f8=read.csv('pl_error_model_21_7kC_8000_lambda_v8.txt',sep=" ",header=FALSE)
f9=read.csv('pl_error_model_21_7kC_8000_lambda_v9.txt',sep=" ",header=FALSE)

f1['condition'] = "v78"
f2['condition'] = "v77"
f3['condition'] = "v76"
f4['condition'] = "v75"
f5['condition'] = "withlambda_v5"
#f5b['condition'] = "withlambda_v5_excl200"

f6['condition'] = "withlambda_v6"
f7['condition'] = "withlambda_v7"
f8['condition'] = "withlambda_v8"
f9['condition'] = "withlambda_v9"

total <- rbind(f1, f2, f3, f4)

# Change colors
p<-ggplot(d, aes(x=V3)) + 
  geom_histogram(color="black", fill="white", binwidth = 1)
p



head(total)
p<-ggplot(total,  aes(x=V3, color=condition )) +
  geom_histogram(fill="white", alpha=0.5, position=position_dodge(width = 0.8), binwidth = 1)+
  #geom_bar(width = 0.8, position = position_dodge(width = 0.9))+
  theme_bw()
  #geom_histogram(color="black", fill="white", binwidth = 1)
p


head(total)

#Compute man placement error)
er = mean(f4$V3)
print(er)



p<-ggplot(f, aes(x=V3)) + 
  stat_ecdf(color="black")
p

ggsave("placement_err_model_12.pdf",width=5,height = 4)

##############################################################################

# Placement erroe before and after pruning based on lambdas


l1=read.csv('pl_error_model_v78.txt',sep=" ",header=FALSE)
l2=read.csv('pl_error_model_0.25_v78.txt',sep=" ",header=FALSE)
l3=read.csv('pl_error_model_0.33_v78.txt',sep=" ",header=FALSE)
l4=read.csv('pl_error_model_0.40_v78.txt',sep=" ",header=FALSE)
l5=read.csv('pl_error_model_0.45_v78.txt',sep=" ",header=FALSE)
l6=read.csv('pl_error_model_lambda_apples_v78.txt',sep=" ",header=FALSE)


l7=read.csv('pl_error_model_lambda_apples_v78_f20_myapples.txt',sep=" ",header=FALSE)
l8=read.csv('pl_error_model_lambda_apples_v78_f20_metinapples.txt',sep=" ",header=FALSE)
l9=read.csv('pl_error_model_v78_f20.txt',sep=" ",header=FALSE)
l10=read.csv('pl_error_model_v78_std_f20.txt',sep=" ",header=FALSE)

total_out = cbind(l1, l4)
total_out = cbind(l1, l2, l3, l4, l5, l6)

total_out = cbind(l7, l8, l9, l10)

write.csv(total_out,'lambda_v78_placement_f_flag.csv')


total <- merge(l1, l4, by = c("V1"))
head(total)


ggplot(aes(x=V3.x, y=V3.y),
            data=total)+
  geom_point(alpha=0.6)+
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  geom_abline(color="red")+
  theme_classic()
#scale_y_continuous(limits = c(0.0000, 0.001))

er = mean(l7$V3)
er
##############################################################################
##############################################################################
##############################################################################


# Correlate lambdas with quality metrics

df_lmbd =read.csv('lambdas_with_quality_v78.csv',sep=",", header=1)

head(df_lmbd)
colnames(df_lmbd)
#df_lmbd = df_lmbd[df_lmbd$lambda_i < 1.0,] # exclude all lambdas that are set to 1.0


nrow(df_lmbd)

res3 <- cor.test(df_lmbd$lambda_i, df_lmbd$contamination_portion, method = "spearman")
res3

res3 <- cor.test(df_lmbd$lambda_i, df_lmbd$clade_separation_score, method = "spearman")
res3



ggplot(aes(y=lambda_i, x=cut(contamination_portion, 20)),
       data=df_lmbd)+
  geom_point(alpha=0.1, size = 1)+
  stat_summary(color = "red")+
  facet_wrap(~clade_separation_score < 0.33)+
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  #geom_abline(color="red")+
  theme_classic()+
  stat_smooth(method = 'lm', aes(colour = 'linear'), se = FALSE)+
  stat_regline_equation(label.y=1.0) +
  stat_cor(aes(label=..rr.label..), label.y=0.9)



ggplot(aes(x=lambda_i, color=contamination_portion < 0.08 | clade_separation_score < 0.45 | reference_representation_score < 0.5),
           data=df_lmbd)+
  stat_ecdf(alpha=0.6)+
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  #geom_abline(color="red")+
  theme_classic()+
  scale_color_discrete(name = "")
  #stat_smooth(method = 'lm', aes(colour = 'linear'), se = FALSE)+
  #stat_regline_equation(label.y=1.0) +
  #stat_cor(aes(label=..rr.label..), label.y=0.9)


ggplot(aes(x=lambda_i, y=reference_representation_score),
       data=df_lmbd)+
  geom_point(alpha=0.6)+
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  #geom_abline(color="red")+
  theme_classic()
  stat_smooth(method = 'lm', aes(colour = 'linear'), se = FALSE)+
  stat_regline_equation(label.y=1.0) +
  stat_cor(aes(label=..rr.label..), label.y=0.9)


ggplot(aes(x=lambda_i, y=Contamination),
       data=df_lmbd)+
  geom_point(alpha=0.6)+
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  #geom_abline(color="red")+
  theme_classic()+
  stat_smooth(method = 'lm', aes(colour = 'linear'), se = FALSE)+
  stat_regline_equation(label.y=1.0) +
  stat_cor(aes(label=..rr.label..), label.y=0.9)


ggplot(aes(x=lambda_i, y=Completeness),
       data=df_lmbd)+
  geom_point(alpha=0.6)+
  facet_wrap(~contamination_portion < 0.08 | clade_separation_score < 0.45 | reference_representation_score < 0.5)+
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  #geom_abline(color="red")+
  theme_classic()+
  stat_smooth(method = 'lm', aes(colour = 'linear'), se = FALSE)+
  stat_regline_equation(label.y=1.0) +
  stat_cor(aes(label=..rr.label..), label.y=12.6)

##############################################################################
##############################################################################
##############################################################################


# Lambda changes per epoch



df_lmbd_per_update =read.csv('lambdas_per_update_v92.csv',sep="\t", header=1)

tail(df_lmbd_per_update)

#df_lmbd_tmp = df_lmbd_per_update[df_lmbd_per_update$epoch == 7800,]
df_lmbd_tmp = df_lmbd_per_update[df_lmbd_per_update$epoch == 7991,]
df_lmbd_tmp_sorted = df_lmbd_tmp[order(df_lmbd_tmp$lambda_i),]
subset = df_lmbd_tmp_sorted[10:40,]




#df_lmbd = df_lmbd[df_lmbd$lambda_i < 1.0,] # exclude all lambdas that are set to 1.0


df_lmbd_clean <- df_lmbd_per_update[df_lmbd_per_update$sample %in% subset$sample, ]

nrow(df_lmbd_clean)

ggplot(aes(x=epoch, y=lambda_i, color = sample),
       data=df_lmbd_clean)+
  geom_line(alpha=0.6)+
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  #geom_abline(color="red")+
  theme_classic()+xlim(7000,8000)+
  scale_color_discrete(guide="none")
  #stat_smooth(method = 'lm', aes(colour = 'linear'), se = FALSE)+
  #stat_regline_equation(label.y=1.0) +
  #stat_cor(aes(label=..rr.label..), label.y=12.6)

ggsave("lambda_per_update_bot_120_130.pdf",width=5,height = 4)


##############################################################################
##############################################################################
##############################################################################

# Lambdas final vs average last 4000 epochs

df_last =read.csv('lambdas_v78.csv',sep="\t", header=1)
df_mean =read.csv('lambdas_per_update_MEAN_v78.csv',sep="\t", header=1)

lambdas_compare <- merge(df_last, df_mean, by = c("sample"))
#lambdas_compare = lambdas_compare[~lambdas_compare$lambda_i.x == 1.0 & lambdas_compare$lambda_i.y == 1.0,]

head(lambdas_compare)


ggplot(aes(x=lambda_i.x, y=lambda_i.y, color = lambda_i_std),
       data=lambdas_compare)+
  geom_point(alpha=0.6)+
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  geom_abline(color="red")+
  theme_classic()+
#scale_y_continuous(limits = c(0.0000, 0.001))


ggsave("lambda_per_update_final_vs_mean_v78.pdf",width=5,height = 4)


##############################################################################
##############################################################################
##############################################################################

# Compare lambdas between trainings with different seeds




##############################################################################
##############################################################################
##############################################################################



df_last =read.csv('pl_error_model_apples_MEAN_v78_f0.txt',sep=" ", header=FALSE)
df_mean =read.csv('pl_error_model_lambda_apples_MEAN_v78_f0_metinsapples.txt',sep=" ", header=FALSE)

lambdas_compare <- merge(df_last, df_mean, by = c("V1"))

lambdas_compare



ggplot(aes(x=V3.x, y=V3.y),
       data=lambdas_compare)+
  geom_point(alpha=0.6)+
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  geom_abline(color="red")+
  theme_classic()
#scale_y_continuous(limits = c(0.0000, 0.001))

write.csv(lambdas_compare,'lambdas_compare_v78_with_mean.csv')


head(lambdas_compare)
##############################################################################
##############################################################################
##############################################################################

# Compare lambda consistency between different seed runs

df1 =read.csv('lambdas_per_update_v91.csv',sep="\t", header=1)
df2 =read.csv('lambdas_per_update_v92.csv',sep="\t", header=1)
target_df1 = df1[df1$epoch == 7991,]
target_df2 = df2[df2$epoch == 7991,]


#df1 =read.csv('lambdas_per_update_v80.csv',sep="\t", header=1)
#df2 =read.csv('lambdas_per_update_v83.csv',sep="\t", header=1)

tail(df1)
tail(df2)


#target_df1 = df1[df1$epoch == 7800,]
#target_df2 = df2[df2$epoch == 7800,]

lambdas_compare <- merge(target_df1, target_df2, by = c("sample"))
#lambdas_compare = lambdas_compare[~lambdas_compare$lambda_i.x == 1.0 & lambdas_compare$lambda_i.y == 1.0,]

tail(lambdas_compare)




ggplot(aes(x=lambda_i.x, y=lambda_i.y),
       data=lambdas_compare)+
  geom_point(alpha=0.6)+
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  geom_abline(color="red")+
  theme_classic()
#scale_y_continuous(limits = c(0.0000, 0.001))


##############################################################################
##############################################################################
##############################################################################

# Compare lambdas vs closest distance


df1 =read.csv('lambdas_vs_closestdist_v78.csv',sep=",", header=1)

head(df1)



lambdas_compare = df1
#lambdas_compare <- merge(target_df1, target_df2, by = c("sample"))
#lambdas_compare = lambdas_compare[~lambdas_compare$lambda_i.x == 1.0 & lambdas_compare$lambda_i.y == 1.0,]






ggplot(aes(x=lambda_i, y=closest_dist),
       #data=lambdas_compare[lambdas_compare$lambda_i<1.0 & lambdas_compare$lambda_i>=0.4,])+
       data=lambdas_compare[lambdas_compare$lambda_i<=1.0 & lambdas_compare$lambda_i>=0.0,])+
  geom_point(alpha=0.6)+
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  #geom_abline(color="red")+
  theme_classic()+
  scale_y_sqrt()+
  stat_smooth(method = 'lm', aes(colour = 'linear'), se = FALSE)+
  stat_regline_equation(label.y=2.0) +
  stat_cor(aes(label=..rr.label..), label.y=1.9)
#scale_y_continuous(limits = c(0.0000, 0.001))



##############################################################################
##############################################################################
##############################################################################

goodbad_df=read.csv('Good_bad_contamination.csv',sep=",",header=TRUE)


colnames(goodbad_df)
goodbad_df


# Change colors
p<-ggplot(goodbad_df, aes(x=Assembly, y=Completeness, color=queires)) + 
  geom_point()
p

p<-ggplot(goodbad_df, aes(x=Assembly, y=Contamination, color=queires)) + 
  geom_point()
p


for(i in colnames(goodbad_df)){
  #print(goodbad_df[[i]])
  #print(i)
  
  if(is.numeric(goodbad_df[[i]])){
   
    y1 = goodbad_df[goodbad_df$queires=='bad',]
    y2 = goodbad_df[goodbad_df$queires=='good',]
    
    #print( t.test(y1[[i]], y2[[i]], paired=FALSE, var.equal = FALSE)$p.value)
    cat ( i, t.test(y1[[i]], y2[[i]], paired=FALSE, var.equal = FALSE)$p.value, "\n")
  }
  
 
  
}

y1 = goodbad_df[goodbad_df$queires=='bad',]
y2 = goodbad_df[goodbad_df$queires=='good',]


t.test(y1$Contamination, y2$Contamination, paired=FALSE, var.equal = FALSE)
t.test(y1$Completeness, y2$Completeness, paired=FALSE, var.equal = FALSE)

  

##############################################################################
##############################################################################
##############################################################################

  geom_vline(xintercept = 0.75, linetype=3,color="gray50")+
  theme_classic()+
  
  
  scale_x_continuous(name="Support",labels=percent)+
  scale_y_continuous(name="ECDF")+
  scale_color_manual(name="", values = c(my_colors[2], my_colors [1], my_colors[8] ), 
                     labels=c("Consensus", "Main","Bin median"))+
  scale_linetype_manual(name="", values = c(1, 2), 
                        labels=c("Incorrect", "Correct"))+
  theme(legend.position = c(.20,.25), 
        legend.margin=margin(t = 0.0, unit='cm') )+
  guides(colour = guide_legend(title = NULL, order = 1, reverse=TRUE, ), 
         linetype = guide_legend(title = NULL, order = 2, reverse=FALSE,))

  
  
# sed '$!N;s/\n/,/' kmer_pl_3layer_lowerLrRate.out | awk 'BEGIN { FS = "[, ]" } ; {print $9, $21}' > kmer_pl_3layer_lowerLrRate.csv
