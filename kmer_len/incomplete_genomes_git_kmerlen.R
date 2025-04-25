

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


##################
# Multiple k placement error

# Claded Unchunked model with computed clade
#df_k3=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k3_v20_8k_s28/all.pl_err',sep=" ",header=FALSE)
#df_k4=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k4_v20_8k_s28/all.pl_err',sep=" ",header=FALSE)
#df_k5=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k5_v20_8k_s28/all.pl_err',sep=" ",header=FALSE)
#df_k6=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k6_v20_8k_s28/all.pl_err',sep=" ",header=FALSE)
##df_k7=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28/all.pl_err',sep=" ",header=FALSE)
#df_k7=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24/all.pl_err',sep=" ",header=FALSE)
#df_k8=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k8_v20_8k_s28/all.pl_err',sep=" ",header=FALSE)

df_k3=read.csv('k3_v20_8k_s28_all.pl_err',sep=" ",header=FALSE)
df_k4=read.csv('k4_v20_8k_s28_all.pl_err',sep=" ",header=FALSE)
df_k5=read.csv('k5_v20_8k_s28_all.pl_err',sep=" ",header=FALSE)
df_k6=read.csv('k6_v20_8k_s28_all.pl_err',sep=" ",header=FALSE)
#df_k7=read.csv('k7_v20_8k_s28_all.pl_err',sep=" ",header=FALSE)
df_k7=read.csv('k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_all.pl_err',sep=" ",header=FALSE)
df_k8=read.csv('k8_v20_8k_s28_all.pl_err',sep=" ",header=FALSE)




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

