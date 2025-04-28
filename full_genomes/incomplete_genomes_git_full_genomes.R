

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



###################################################### ################## ##################
###################################################### ################## ##################
###################################################### ################## ##################


# Comparison of full genome 500 queries on different models

# Uncladed, unchunked
#df_UU =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_single_clade/all.pl_err',sep=" ",header=FALSE)
df_UU =read.csv('k7_v20_8k_s28_single_clade_all.pl_err',sep=" ",header=FALSE)
df_UU ["claded"] = "Uncladed"
df_UU ["chunked"] = "Unchunked"
df_UU ["clsf"] = "na"
df_UU
df_UU = df_UU[!(is.na(df_UU$V2)), ]
df_UU

summary(df_UU $V3)

# Uncladed, chunked MODEL
#df_BU =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_single_clade_ChunkModel/all.pl_err',sep=" ",header=FALSE)
df_BU =read.csv('k7_v20_8k_s28_single_clade_ChunkModel_all.pl_err',sep=" ",header=FALSE)
df_BU ["claded"] = "Uncladed"
df_BU ["chunked"] = "Chunked"
df_BU ["clsf"] = "na"
nrow(df_BU)

# Claded, unchunked - True clade
#df_UC_tr =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_TRUE_CLADE_s24/all.pl_err',sep=" ",header=FALSE)
df_UC_tr =read.csv('k7_v20_8k_s28_FullGenomQueries_TRUE_CLADE_s24_all.pl_err',sep=" ",header=FALSE)
df_UC_tr ["claded"] = "Claded"
df_UC_tr ["chunked"] = "Unchunked"
df_UC_tr["clsf"] = "tr"
nrow(df_UC_tr)
df_UC_tr

# Claded, unchunked - Computed clade
#df_UC_clsf =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24/all.pl_err',sep=" ",header=FALSE)
df_UC_clsf =read.csv('k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_all.pl_err',sep=" ",header=FALSE)
df_UC_clsf ["claded"] = "Claded"
df_UC_clsf ["chunked"] = "Unchunked"
df_UC_clsf["clsf"] = "cmp"
summary(df_UC_clsf$V4)

# Claded, chunked - True clade
#df_BC_tr =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_Chunks_FullGenomQueries_TRUE_CLADE/all.pl_err',sep=" ",header=FALSE)
df_BC_tr =read.csv('k7_v37_8k_s28_TrainClassf_10K_TOL_Chunks_FullGenomQueries_TRUE_CLADE_all.pl_err',sep=" ",header=FALSE)
df_BC_tr ["claded"] = "Claded"
df_BC_tr ["chunked"] = "Chunked"
df_BC_tr["clsf"] = "tr"

# Claded, chunked - Computed clade
#df_BC_clsf =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_FullGenomQueries_COMPUTED_CLADE/all.pl_err',sep=" ",header=FALSE)
df_BC_clsf =read.csv('k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_FullGenomQueries_COMPUTED_CLADE_all.pl_err',sep=" ",header=FALSE)
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

#q = read.csv('/Users/nora/Documents/ml_metagenomics/tol_quality_scores/quality_comparison_hor.csv')
q = read.csv('quality_comparison_hor.csv')
head(q)

#nov = read.csv('/Users/nora/Documents/ml_metagenomics/pendant.txt',sep="\t",h=F)
#nov$V2 = nov$V2*2/100
#head(nov)

#nov = read.csv('/Users/nora/Documents/ml_metagenomics/closest_dev_set_2col.txt',sep="\t",h=FALSE)
nov = read.csv('closest_dev_set_2col.txt',sep="\t",h=FALSE)
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
#df_true_clades =read.csv('/Users/nora/Documents/ml_metagenomics/clade_targets.txt',sep=" ",header=TRUE)
df_true_clades =read.csv('clade_targets.txt',sep=" ",header=TRUE)
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
#df_bl =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
df_bl =read.csv('k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
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

#df_param1 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_extra_layer/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
df_param1 =read.csv('k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_extra_layer/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_param1["cond"] = "extra_layer"
df_param1["cond2"] = "layer"

#df_param2 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_extra_2layer/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
df_param2 =read.csv('k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_extra_2layer/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_param2["cond"] = "extra_2layer"
df_param2["cond2"] = "layer"


#df_param3 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_hilr/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
df_param3 =read.csv('k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_hilr/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_param3["cond"] = "hilr"
df_param3["cond2"] = "lr"

#df_param4 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_hilr2/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
df_param4 =read.csv('k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_hilr2/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_param4["cond"] = "hilr2"
df_param4["cond2"] = "lr"

#df_param5 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_hilr_hilrmin/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
df_param5 =read.csv('k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_hilr_hilrmin/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_param5["cond"] = "hilr_hilrmin"
df_param5["cond2"] = "lrmin"

#df_param6 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_hilr2_hilrmin/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
df_param6 =read.csv('k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_hilr2_hilrmin/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_param6["cond"] = "hilr2_hilrmin"
df_param6["cond2"] = "lrmin"

#df_param7 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_batch32/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
df_param7 =read.csv('k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_batch32/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
#df_UC_clsf ["claded"] = "Claded"
#df_UC_clsf ["chunked"] = "Unchunked"
#df_UC_clsf["clsf"] = "cmp"
df_param7["cond"] = "batch32"
df_param7["cond2"] = "other"

#df_param8 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_decay4k/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
df_param8 =read.csv('k7_v20_8k_s28_FullGenomQueries_COMPUTED_CLADE_s24_v60_decay4k/summary_FULL_GENOME_D1_CladedUnchunkedModel_param_sweep.pl_error',sep=" ",header=TRUE)
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

