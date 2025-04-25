

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
require(plyr)
require(tidyr)



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

################################################################
############ UPDATED MODEL #####################################


#pl_per_contig=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v38_8k_s28_clade_All_TOL_CONTIGS_Chunked_Model/summary_FULL_GENOME_CONTIGS.pl_error', sep=" ",header=T)

#pl_per_contig=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v38_8k_s28_clade_All_TOL_CONTIGS_Full_Model/summary_FULL_GENOME_CONTIGS_FullModel.pl_error', sep=" ",header=T)


#pl_per_contig=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v38_8k_s28_clade_All_TOL_CONTIGS_Global_Model/summary_TOL_full_genome_contigs_six_models.pl_error', sep=" ",header=T)
pl_per_contig=read.csv('summary_TOL_full_genome_contigs_six_models.pl_error', sep=" ",header=T)
#pl_per_contig=read.csv('v11.3_PerClade_v2_wol_EXTENDED_pl_vs_contig_len_total.txt',sep=" ",header=T)
nrow(pl_per_contig)

tail(pl_per_contig)
unique(pl_per_contig$model)
colnames(pl_per_contig)
#names(df_true_clades)[names(df_true_clades) == 'genome'] <- 'V1'

# Uncladed, unchunked - current Global
# Uncladed, chunked - current Global_Chunk

# Claded, unchunked - True clade - Don't have
# Claded, unchunked - current Full
# Claded, chunked - True clade - - Don't have
# Claded, chunked - Computed clade - current ChunkExp



#For ppt
pl_per_contig
pl_per_contig$cond2 <- ""
pl_per_contig$cond2 <- ifelse(pl_per_contig$length_adj <=320000, "Short", ifelse((pl_per_contig$length_adj >320000 & pl_per_contig$length_adj <=1280000 ), "Optimal", "Chimeric"))
pl_per_contig$cond2 <- factor(pl_per_contig$cond2, levels = c("Short", "Optimal", "Chimeric"))
nrow(pl_per_contig)
colnames(pl_per_contig)

# Put all the color values (in hex format) from Dark2 into a vector
myPal <- brewer.pal(8,"Dark2")

# Remove whichever color you don't want
myPal <- myPal[-1]

unique(pl_per_contig$cond2)

# Chunk: remove (uniform)
# ChunkExp: default
# ChunkExp_Global_*: chunk+gloal,local
# Full: No chunking, but clading
# Global: No chunk, no clading
# Global_Chunk: No clading, chunking


#pl_per_contig[pl_per_contig$length_adj >10000,] %>% 
#  model %in% c("ChunkExp","ChunkExp_Global_p_0.9","Full", "Global_Chunk", "ChunkExp_v55",))+
unique(pl_per_contig$model)



pl_per_contig$model[pl_per_contig$model == "ChunkExp"] <- "Claded_chunked"
pl_per_contig$model[pl_per_contig$model == "Global"] <- "Uncladed_unchunked"
pl_per_contig$model[pl_per_contig$model == "Global_Chunk"] <- "Uncladed_chunked"
pl_per_contig$model[pl_per_contig$model == "Full"] <- "Claded_unchunked"
pl_per_contig$model[pl_per_contig$model == "ChunkExp_TrueClade"] <- "Claded_chunked_true"
pl_per_contig$model[pl_per_contig$model == "Full_TrueClade"] <- "Claded_unchunked_true"

unique(pl_per_contig$model)



pl_per_contig$model

  ggplot(aes(x=cut(length_adj/10^3,c(0, 0.5, 1, 2, 4, 8,16,32,64,128,256,2000)*10), y=pl_err),
       #data=pl_per_contig[pl_per_contig$model %in% c("ChunkExp", "Full", "Global", "Global_Chunk") & pl_per_contig$length_adj > 5000,])+
    data=pl_per_contig[pl_per_contig$model %in% c("Claded_chunked", "Uncladed_unchunked", "Uncladed_chunked", "Claded_unchunked", "ChunkExp_TrueClade", "Full_TrueClade", "Chunk") & pl_per_contig$length_adj > 5000,])+
    #data=pl_per_contig)+
    stat_summary(aes(group=model, color = model),geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8), position = position_dodge())+
    #stat_summary(geom="line",color="blue",aes(group="1"))+
    stat_summary(geom="line",aes(group=model, color = model))+
    stat_summary(aes(group=model, color = model),alpha = 0.3)+
    theme_classic()+
    geom_hline(yintercept = 3)+
    coord_cartesian(ylim=c(0, 15))+
    scale_color_brewer(palette = "Paired")
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  #scale_x_discrete(name="Contig length (KB)", label = c("(10-20]", "(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-20000]"))  +
  #scale_colour_brewer(palette = "Set1", name="", labels = c("Test", "Train"))+
  #scale_fill_brewer(palette = "Accent", name="")+
  #scale_fill_brewer(palette = "Set2", name="")+
  #scale_fill_manual(values = c(  "#d1e5f0", "#e6f5d0", "#fddbc7"), name="")+
  #scale_fill_manual(values = c(  "#d1e5f0", "#f7f7f7", "#fddbc7"), name="")+
  #scale_fill_manual(values = c(  "#E7298A", "#66A61E", "#A6761D"), name="")+
  #scale_fill_manual(values = c(  "#BEAED4", "#7FC97F", "#FDC086"), name="")+
  ylab("Placement error")+
  theme(axis.text.x = element_text(angle = 15, vjust = 0.7, hjust=0.5))+
  theme(legend.position = c(0.80,0.80), legend.margin=margin(0,0,0,0),legend.direction="vertical", legend.title=element_blank(),
        legend.box = "horizontal"
        #axis.text.x = element_text(size = 8)
  )+
  coord_cartesian(ylim=c(0, 12))+
  annotate(geom="text",x=0.8,y=12.3,label="24")+
  scale_colour_brewer(palette = "Dark2", name="")+
  scale_fill_brewer(palette = "Set2", name="")
  geom_hline(yintercept = 3)



ggsave("wol_dev_queries_pl_percontig_v2_UPDATED_MODEL.pdf",width=5,height = 4)

require(tidyverse)
pl_per_contig[pl_per_contig$model %in% c("Claded_chunked", "Uncladed_unchunked", "Uncladed_chunked", "Claded_unchunked", "Claded_chunked_true","Claded_unchunked_true") & pl_per_contig$length_adj > 5000,] %>%
  mutate(cl=true_clade==top_class) %>%
  dplyr ::group_by(model,) %>% dplyr::summarise(cla = mean(cl))

head(pl_per_contig)
ggplot(aes(color=reorder(ifelse(grepl("_chunked",model),"Chunked","Unchunked"),pl_err), # reorder(paste(claded,ifelse(clsf=="tr","(True)",""),ifelse(chunked=="Chunked","(Chunked)","")),V3), 
          x=cut(length_adj/10^3,c(0, 0.5, 1, 2, 3, 5, 10, 20,50,200,1500)*10), 
          group=model,
          y=pl_err,
          linetype=paste(sub("_.*","",model),ifelse(grepl("true",model),"(true)","") )),
    #data=pl_per_contig[pl_per_contig$model %in% c("ChunkExp", "Full", "Global", "Global_Chunk") & pl_per_contig$length_adj > 5000,])+
    data=pl_per_contig[pl_per_contig$model %in% c("Claded_chunked", "Uncladed_unchunked", "Uncladed_chunked", "Claded_unchunked", "Claded_chunked_true","Claded_unchunked_true") & pl_per_contig$length_adj > 5000,]
    )+
  #data=df_fin[df_fin$chunked %in% c("Unchunked"),])+
  #geom_boxplot(outlier.shape=NA)+
  #geom_boxplot(fill="skyblue1", outlier.shape=NA)+
  stat_summary(geom="line" )+
  #stat_summary(aes(group=condition),geom="bar" ,position = position_dodge(0.9),color="black")+
  #stat_summary(position = position_dodge(0.9),color="black")+
  stat_summary( alpha = 0.7,geom="point")+
  stat_summary( alpha = 0.7,geom="errorbar",width=0.1,linetype=1,size=0.3)+
  scale_colour_brewer(palette = "Dark2", name="")+
  scale_linetype(name="", labels = c("Cladded", "Cladded (true)", "Uncladded"))+
  #scale_color_manual(name = "",values=c("#1f78b4","#e31a1c"))+
  theme_classic()+
  theme(legend.position = c(0.88,0.75))+
  scale_y_continuous("Placement error")+scale_x_discrete("Contig length (KB)")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1.0, hjust=1), axis.title.x = element_text(vjust=5.5))+
  theme(plot.margin = margin(5.5,5.5,-2.5,5.5, "pt"))
ggsave("Contigs-chunks-D2-line.pdf",width=6.5,height = 4.5)

getwd()

ggplot(aes(color=reorder(ifelse(grepl("_chunked",model),"Chunked","Unchunked"),pl_err), # reorder(paste(claded,ifelse(clsf=="tr","(True)",""),ifelse(chunked=="Chunked","(Chunked)","")),V3), 
           x=cut(length_adj/10^3,c(0, 0.5, 1, 2, 3, 5, 10, 20,50,200,1500)*10), 
           y=pl_err,
           fill=paste(sub("_.*","",model),ifelse(grepl("true",model),"(true)","") )),
       #data=pl_per_contig[pl_per_contig$model %in% c("ChunkExp", "Full", "Global", "Global_Chunk") & pl_per_contig$length_adj > 5000,])+
       data=pl_per_contig[pl_per_contig$model %in% c("Claded_chunked", "Uncladed_unchunked", "Uncladed_chunked", "Claded_unchunked", "Claded_chunked_true","Claded_unchunked_true") & pl_per_contig$length_adj > 5000,]
)+
  geom_boxplot(outlier.alpha = 0.1,outlier.size = 0.4)+
  stat_summary(geom="point",position = position_dodge2(width = 0.75))+
  scale_colour_manual( name="",values=c("black","#3050CC"))+
  scale_linetype(name="")+
  #scale_color_manual(name = "",values=c("#1f78b4","#e31a1c"))+
  theme_classic()+
  scale_fill_brewer(palette = "Pastel1",name="", labels = c("Cladded", "Cladded (true)", "Uncladded"))+
  theme(legend.position = c(0.88,0.75))+
  scale_y_continuous("Placement error")+scale_x_discrete("Contig length (KB)")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1.0, hjust=1), axis.title.x = element_text(vjust=5.5))+
  theme(plot.margin = margin(5.5,5.5,-2.5,5.5, "pt"))
ggsave("Contigs-chunks-D2-box.pdf",width=8.5,height = 5.5)

getwd()


ggplot(aes(color=reorder(ifelse(grepl("_chunked",model),"Chunked","Unchunked"),pl_err), # reorder(paste(claded,ifelse(clsf=="tr","(True)",""),ifelse(chunked=="Chunked","(Chunked)","")),V3), 
           x=cut(length_adj/10^3,c(0, 0.5, 1, 2, 3, 5, 10, 20,50,200,1500)*10), 
           group=model,
           y=ifelse(top_class==true_clade,1,0),
           #linetype=paste(sub("_.*","",model),ifelse(grepl("true",model),"(true)","") )
           ),
       #data=pl_per_contig[pl_per_contig$model %in% c("ChunkExp", "Full", "Global", "Global_Chunk") & pl_per_contig$length_adj > 5000,])+
       data=pl_per_contig[pl_per_contig$model %in% c("Claded_chunked", "Claded_unchunked") & pl_per_contig$length_adj > 5000,]
)+
  #data=df_fin[df_fin$chunked %in% c("Unchunked"),])+
  #geom_boxplot(outlier.shape=NA)+
  #geom_boxplot(fill="skyblue1", outlier.shape=NA)+
  stat_summary(geom="line" )+
  #stat_summary(aes(group=condition),geom="bar" ,position = position_dodge(0.9),color="black")+
  #stat_summary(position = position_dodge(0.9),color="black")+
  stat_summary( alpha = 0.7,geom="point")+
  stat_summary( alpha = 0.7,geom="errorbar",width=0.1,linetype=1,size=0.3)+
  scale_colour_brewer(palette = "Dark2", name="")+
  scale_linetype(name="")+
  #scale_color_manual(name = "",values=c("#1f78b4","#e31a1c"))+
  theme_classic()+
  theme(legend.position = c(0.88,0.25))+
  scale_y_continuous("Classification accuracy")+scale_x_discrete("Contig length (KB)")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1.0, hjust=1), axis.title.x = element_text(vjust=5.5))+
  theme(plot.margin = margin(5.5,5.5,-2.5,5.5, "pt"))

ggsave("Contigs-D3-classified-line.pdf",width=6.5,height = 4.5)


# Chunk vs Chunk_EXP for training classifier

pl_per_contig[pl_per_contig$model %in% c("Claded_chunked", "Uncladed_unchunked", "Uncladed_chunked", "Claded_unchunked", "Claded_chunked_true","Claded_unchunked_true") & pl_per_contig$length_adj > 5000,] %>%
  mutate(Chunked=ifelse(grepl("_chunked",model),"Chunked","Unchunked"),
         Claded=paste(sub("_.*","",model),ifelse(grepl("true",model),"(true)","") ),
         lenbin = cut(length_adj/10^3,c(0, 0.5, 1, 2, 3, 5, 10, 20,50,200,1500)*10))  %>%
  select(pl_err,Claded,Chunked,lenbin) %>%
  dplyr::group_by(Claded,Chunked,lenbin) %>%
  dplyr::summarise(merror=mean(pl_err)) %>%
  #pivot_wider(names_from = Chunked,values_from = merror) %>%
ggplot(aes(color=Claded, 
           x=lenbin, xend=lenbin,
           group=interaction(Claded,lenbin),
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
  #coord_cartesian(ylim=c(0,14))+
  theme(axis.text.x = element_text(angle = 30, vjust = 1.0, hjust=1), axis.title.x = element_text(vjust=5.5))+
  theme(plot.margin = margin(5.5,5.5,-2.5,5.5, "pt"))
#ggsave("Contigs-chunks-D2-arrow.pdf",width=7.3,height = 4.5)
ggsave("Contigs-chunks-D2-arrow.pdf",width=4.8,height = 4.0)
getwd()
#legend.key.size = unit(1.0, "cm")


ggplot(aes(color=cut(length_adj/10^3,c(0, 0.5, 1, 2, 3, 5, 10, 20,50,200,1500)*10), 
           x=pl_err,
           ),
       #data=pl_per_contig[pl_per_contig$model %in% c("ChunkExp", "Full", "Global", "Global_Chunk") & pl_per_contig$length_adj > 5000,])+
       data=pl_per_contig[pl_per_contig$model %in% c("Claded_chunked") & pl_per_contig$length_adj > 5000,]
)+
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
#ggsave("Contigs-chunks-D2-ECDF.pdf",width=6.5,height = 4.5)
ggsave("Contigs-chunks-D2-ECDF.pdf",width=4.8,height = 4.0)


unique(pl_per_contig$model)
r= pl_per_contig[pl_per_contig$model=="Claded_chunked",]
r= pl_per_contig[pl_per_contig$model=="Claded_unchunked",]
nrow(r)
summary(r[r$length_adj>1200000,]$pl_err)
summary(r[r$length_adj>=10000 & r$length_adj<=20000,]$pl_err)
summary(r[r$length_adj<=50000 & r$length_adj>20000,]$pl_err)

data = r[r$length_adj>50000,]$pl_err
quantile(data, probs = c(.25, .5, .90))

nrow(r[r$length_adj>1200000,])/nrow(r)*100

t1= r[r$length_adj>1200000,]$pl_err
sum(t1>4)

nrow(r[r$length_adj>50000,])

chunk_clsf<-pl_per_contig[pl_per_contig$model %in% c("Chunk"),]
chunkExp_clsf<- pl_per_contig[pl_per_contig$model %in% c("Claded_chunked"),]



for(x in c(0, 5000, 10000, 20000, 40000, 80000)){
  print(x)
  
  x1 = nrow(chunk_clsf[chunk_clsf$top_class==chunk_clsf$true_clade & chunk_clsf$length_adj > x,]) /nrow(chunk_clsf[chunk_clsf$length_adj > x,])
  x2 = nrow(chunkExp_clsf[chunkExp_clsf$top_class==chunkExp_clsf$true_clade & chunkExp_clsf$length_adj > x,]) /nrow(chunkExp_clsf[chunkExp_clsf$length_adj > x,])
  #combo_chunk_df = 
  print (x1)
  print (x2)
  chunk_clsf$clsf_acc = ifelse (chunk_clsf$length_adj > x, x1, chunk_clsf$clsf_acc )
  chunkExp_clsf$clsf_acc = ifelse (chunkExp_clsf$length_adj > x, x2, chunkExp_clsf$clsf_acc )
}

unique(chunk_clsf$clsf_acc)

combo_classif <- rbind(chunk_clsf, chunkExp_clsf)


ggplot(aes(x=cut(length_adj/10^3,c(0, 0.5, 1, 2, 4, 8, 2000)*10), y=clsf_acc),
       #data=pl_per_contig[pl_per_contig$model %in% c("ChunkExp", "Full", "Global", "Global_Chunk") & pl_per_contig$length_adj > 5000,])+
       data=combo_classif)+
  #data=pl_per_contig)+
  stat_summary(geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8),aes(group=model, color = model))+
  #stat_summary(geom="line",color="blue",aes(group="1"))+
  stat_summary(geom="line",aes(group=model, color = model))+
  stat_summary(alpha = 0.3,aes(group=model, color = model))+
  theme_classic()+
  geom_hline(yintercept = 3)+
  coord_cartesian(ylim=c(0.7, 1))
#scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
#scale_x_discrete(name="Contig length (KB)", label = c("(10-20]", "(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-20000]"))  +
#scale_colour_brewer(palette = "Set1", name="", labels = c("Test", "Train"))+
#scale_fill_brewer(palette = "Accent", name="")+
#scale_fill_brewer(palette = "Set2", name="")+
#scale_fill_manual(values = c(  "#d1e5f0", "#e6f5d0", "#fddbc7"), name="")+
#scale_fill_manual(values = c(  "#d1e5f0", "#f7f7f7", "#fddbc7"), name="")+
#scale_fill_manual(values = c(  "#E7298A", "#66A61E", "#A6761D"), name="")+
#scale_fill_manual(values = c(  "#BEAED4", "#7FC97F", "#FDC086"), name="")+
ylab("Placement error")+
  theme(axis.text.x = element_text(angle = 15, vjust = 0.7, hjust=0.5))+
  theme(legend.position = c(0.80,0.80), legend.margin=margin(0,0,0,0),legend.direction="vertical", legend.title=element_blank(),
        legend.box = "horizontal"
        #axis.text.x = element_text(size = 8)
  )+
  coord_cartesian(ylim=c(0, 12))+
  annotate(geom="text",x=0.8,y=12.3,label="24")+
  scale_colour_brewer(palette = "Dark2", name="")+
  scale_fill_brewer(palette = "Set2", name="")
geom_hline(yintercept = 3)





ggplot(aes(x=cut(length_adj/10^3,c(0, 0.5, 1, 2, 4, 8,16,32,64,128,256,2000)*10), y=pl_err,
           aes(group=model, color = model)),
       #data=pl_per_contig[pl_per_contig$model %in% c("ChunkExp", "Full", "Global", "Global_Chunk") & pl_per_contig$length_adj > 5000,])+
       data=pl_per_contig[pl_per_contig$model %in% c("Claded_chunked", "Chunk") & pl_per_contig$length_adj > 5000,])+
  #data=pl_per_contig)+
  stat_summary(geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8), position = position_dodge(),aes(group=model, color = model))+
  #stat_summary(geom="line",color="blue",aes(group="1"))+
  stat_summary(geom="line",aes(group=model, color = model))+
  stat_summary(alpha = 0.3,aes(group=model, color = model))+
  theme_classic()+
  geom_hline(yintercept = 3)+
  coord_cartesian(ylim=c(0, 15))
#scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
#scale_x_discrete(name="Contig length (KB)", label = c("(10-20]", "(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-20000]"))  +
#scale_colour_brewer(palette = "Set1", name="", labels = c("Test", "Train"))+
#scale_fill_brewer(palette = "Accent", name="")+
#scale_fill_brewer(palette = "Set2", name="")+
#scale_fill_manual(values = c(  "#d1e5f0", "#e6f5d0", "#fddbc7"), name="")+
#scale_fill_manual(values = c(  "#d1e5f0", "#f7f7f7", "#fddbc7"), name="")+
#scale_fill_manual(values = c(  "#E7298A", "#66A61E", "#A6761D"), name="")+
#scale_fill_manual(values = c(  "#BEAED4", "#7FC97F", "#FDC086"), name="")+
ylab("Placement error")+
  theme(axis.text.x = element_text(angle = 15, vjust = 0.7, hjust=0.5))+
  theme(legend.position = c(0.80,0.80), legend.margin=margin(0,0,0,0),legend.direction="vertical", legend.title=element_blank(),
        legend.box = "horizontal"
        #axis.text.x = element_text(size = 8)
  )+
  coord_cartesian(ylim=c(0, 12))+
  annotate(geom="text",x=0.8,y=12.3,label="24")+
  scale_colour_brewer(palette = "Dark2", name="")+
  scale_fill_brewer(palette = "Set2", name="")
geom_hline(yintercept = 3)








colnames(pl_per_contig)
head(pl_per_contig)

ggplot(aes(x=cut(kmer_sums/10^3,c(2, 4, 8,16,32,64,128,256,1000)*10), y=pl_err),
       data=pl_per_contig)+
  #geom_violin()+
  stat_summary(geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  stat_summary(geom="line", color=cond2, aes(group="1"))+
  stat_summary(alpha = 0.9)+
  theme_classic()
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="Contig length (KB)", label = c("(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-10000]"))  +
  #scale_colour_brewer(palette = "Set1", name="", labels = c("Test", "Train"))+
  #scale_fill_brewer(palette = "Accent", name="")+
  #scale_fill_brewer(palette = "Set2", name="")+
  #scale_fill_manual(values = c(  "#d1e5f0", "#e6f5d0", "#fddbc7"), name="")+
  
  #scale_fill_manual(values = c(  "#d1e5f0", "#f7f7f7", "#fddbc7"), name="")+
  #scale_fill_manual(values = c(  "#E7298A", "#66A61E", "#A6761D"), name="")+
  #scale_fill_manual(values = c(  "#BEAED4", "#7FC97F", "#FDC086"), name="")+
  ylab("Placement error")+
  theme(axis.text.x = element_text(angle = 15, vjust = 0.7, hjust=0.5))+
  theme(legend.position = c(0.9,0.85), legend.margin=margin(0,0,0,0),legend.direction="vertical", legend.title=element_blank(),
        #axis.text.x = element_text(size = 8)
  )+
  guides(fill = guide_legend(reverse = FALSE))
geom_hline(yintercept = 3)




################################################################
# CAMI dataset on 65K TOL CHUNK EXP MODEL


pl_per_contig1=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/cami_gsa_pooled/summary_cami_long_reads_gsa_65KTOL_ChunkExp.pl_error',sep=" ",header=T)
pl_per_contig2=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/cami_gsa_pooled/summary_cami_short_reads_gsa_65KTOL_ChunkExp.pl_error',sep=" ",header=T)
pl_per_contig3=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/cami_gsa_pooled/summary_cami_long_reads_65KTOL_ChunkExp_sample_0.pl_error',sep=" ",header=T)
pl_per_contig3

head(pl_per_contig1)


pl_per_contig1["conditions"] = "gsa_long"
pl_per_contig2["conditions"] = "gsa_short"
pl_per_contig3["conditions"] = "sample_0_long"

#pl_per_contig <-rbind(pl_per_contig1, pl_per_contig2, pl_per_contig3, fill = TRUE)
pl_per_contig <-rbind(pl_per_contig1, pl_per_contig2)
#data.frame(pl_per_contig[pl_per_contig$kmer_sums <=40000 & pl_per_contig$kmer_sums >=20000,])[,"pl_err_delta"]

#pl_per_contig_no_NA <- na.omit(pl_per_contig)

colnames(pl_per_contig)

ggplot(aes(x=cut(length_adj/10^3,c(0.5, 1, 2, 4, 8,16,32,64,128,256,1000)*10), y=pl_err, color = conditions),
       data=pl_per_contig[pl_per_contig$length_adj > 5000,])+
  #geom_violin()+
  stat_summary()+
  stat_summary(geom="line", aes(group=conditions, color = conditions))+
  scale_colour_brewer(palette = "Dark2", name="")+
  #stat_summary(geom="line",color="blue",aes(group="1"))+
  stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()+
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="Contig length (KB)", label = c("(5-10]", "(10-20]", "(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-10000]"))  +
  ylab("Placement error")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.7, hjust=0.5))+
  theme(legend.position = c(0.8,0.8), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )


ggplot(pl_per_contig, aes(x=length_adj))+
  geom_histogram(color="darkblue", fill="lightblue")



# Read gsa short with lineage

pl_per_contig_lin=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/cami_gsa_pooled/summary_cami_short_reads_gsa_65KTOL_ChunkExp_wDist_wLineage.pl_error',sep="\t",header=T)

head(pl_per_contig_lin)

colnames(pl_per_contig_lin)

pl_per_contig_lin["conditions"] = "gsa_short"

ggplot(aes(x=cut(length_adj/10^3,c(0.5, 1, 2, 4, 8,16,32,64,128,256,1000)*10), 
           y=pl_err, color = conditions),
       data=pl_per_contig_lin[pl_per_contig_lin$length_adj > 5000,])+
  #geom_violin()+
  stat_summary()+
  stat_summary(geom="line", aes(group=conditions, color = conditions))+
  scale_colour_brewer(palette = "Dark2", name="")+
  #stat_summary(geom="line",color="blue",aes(group="1"))+
  stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()+
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="Contig length (KB)", label = c("(5-10]", "(10-20]", "(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-10000]"))  +
  ylab("Placement error")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.7, hjust=0.5))+
  theme(legend.position = c(0.8,0.8), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )

summary(pl_per_contig_lin$near_dist)
ggplot(aes(x=cut(near_dist/100,breaks=c(0, 0.01, 0.05,0.15,0.2,0.5, 4), include.lowest = TRUE) , y=pl_err, color = conditions),
       data=pl_per_contig_lin[pl_per_contig_lin$length_adj > 5000 & pl_per_contig_lin$self_rank=="species", ])+
  #geom_violin()+
  stat_summary()+
  #stat_summary(geom="line", aes(group=conditions, color = conditions))+
  scale_colour_brewer(palette = "Dark2", name="")+
  #stat_summary(geom="line",color="blue",aes(group="1"))+
  stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()+
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  #scale_x_discrete(name="Contig length (KB)", label = c("(5-10]", "(10-20]", "(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-10000]"))  +
  ylab("Placement error")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.7, hjust=0.5))+
  theme(legend.position = c(0.8,0.8), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )


pl_per_contig_lin$nbr_rank  <- factor(pl_per_contig_lin$nbr_rank, levels = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"))

#& pl_per_contig_lin$self_rank=="species"

ggplot(aes(x=nbr_rank , y=pl_err, color = conditions),
       data=pl_per_contig_lin[pl_per_contig_lin$length_adj > 5000 & pl_per_contig_lin$self_rank=="species",])+
  #geom_violin()+
  #geom_point()+
  stat_summary()+
  #stat_summary(geom="line", aes(group=conditions, color = conditions))+
  scale_colour_brewer(palette = "Dark2", name="")+
  #stat_summary(geom="line",color="blue",aes(group="1"))+
  stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()+
  facet_wrap(~self_rank)
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  #scale_x_discrete(name="Contig length (KB)", label = c("(5-10]", "(10-20]", "(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-10000]"))  +
  ylab("Placement error")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.7, hjust=0.5))
  theme(legend.position = c(0.8,0.8), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )
  
meta = read.csv('/Users/nora/Documents/ml_metagenomics/cami_short_reads/metadata.tsv',sep="\t")
head(meta)
pl_per_contig_lin = merge(pl_per_contig_lin,meta,by.x="otu_x",by.y="genome_ID")

head(pl_per_contig_lin)
  ggplot(aes(x=cut(length_adj/10^3,c(0.5, 4,32,1000)*10) , y=pl_err, color = novelty_category=="known_strain"),
         data=pl_per_contig_lin[pl_per_contig_lin$length_adj > 5000 ,])+
    #geom_violin()+
    #geom_point(alpha = 0.1)+
    #geom_tile()+
    stat_summary()+
    #facet_wrap(~self_rank)+
    stat_summary(geom="line", aes(group=novelty_category=="known_strain"))+
    scale_colour_brewer(palette = "Dark2", name="")+
    #stat_summary(geom="line",color="blue",aes(group="1"))+
    #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
    theme_classic()+
    #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
    #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
    #scale_x_discrete(name="Contig length (KB)", label = c("(5-10]", "(10-20]", "(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-10000]"))  +
    #ylab("Placement error")+
    theme(axis.text.x = element_text(angle = 10, vjust = 0.7, hjust=0.5))
  theme(legend.position = c(0.8,0.8), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )




################################################################
# Plot AMBER results for gsa short
#folder = "amber_closest_subset_5kfilt"

#folder = "amber_closest_subset_hiqual" # also includes >10k cut off

#folder = "amber_closest_subset_5kfilt_above5k"
folder = "amber_closest_subset_5kfilt_above10k"
#folder = "amber_closest_subset_5kfilt_above20k"

#folder = "amber_closest_distance_long_reads_S0_hiqual"
#folder = "amber_closest_distance_long_reads_S0_above10k"

#folder = "amber_closest_subset_combo"
#folder = "amber_closest_subset_hiqual_ALL"

amber_df=read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "results.tsv"),sep="\t",header=T)


nrow(amber_df)
head(amber_df)

colnames(amber_df)
unique(amber_df$Tool)
head(amber_df$accuracy_bp..unfiltered.ool)

amber_df$rank  <- factor(amber_df$rank, levels = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"))

amber_df$Tool

ggplot(aes(x=rank, y=accuracy_seq, fill = Tool,color = Tool, group = Tool),
       data=amber_df[(!grepl("losest",amber_df$Tool) |grepl("v10",amber_df$Tool)) & 
                       !amber_df$Tool %in% c("Closest_tree_dist_phyl","tax2tree", "Closest_tree_dist_cls", "Closest_tree_dist_ordr", "Closest_tree_dist_fam", "Closest_tree_dist_gen", "Gold standard", "closest_subset_v2", "closest_subset_v3"),])+
      #data=amber_df[amber_df$Tool %in% c("closest_subset"),])+
  #geom_violin()+
  stat_summary()+
  stat_summary(geom="line",  alpha=0.7)+
  #stat_summary(geom="bar" , alpha = 0.7,position = position_dodge())+
  #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()+
  coord_cartesian(ylim=c(0.1, 1.0))+
  #scale_color_brewer(name = "",palette="Set1")
  #scale_color_brewer(name = "",palette="Dark2")
  scale_fill_brewer(name = "",palette="Dark2", labels=c("closest_subset_v10" = "kf2d"))+
  scale_color_brewer(name = "",palette="Dark2", labels=c("closest_subset_v10" = "kf2d"))+
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="", label = c("Superkindom", "Phylum", "Class","Order", "Family", "Genus", "Species"))  +
  ylab("Accuracy")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95))+
theme(legend.position = c(0.3,0.4),legend.margin=margin(0,-2,0,0)
      #,
      #axis.text.x = element_text(size = 8)
)

nrow(amber_df)

head(amber_df)

ggplot(aes(x=rank, y=accuracy_seq, fill = reorder(Tool,accuracy_seq)),
       data=amber_df[(!grepl("losest",amber_df$Tool) |grepl("v10",amber_df$Tool)) & 
                       !amber_df$Tool %in% c("Closest_tree_dist_phyl","tax2tree", "Closest_tree_dist_cls", "Closest_tree_dist_ordr", "Closest_tree_dist_fam", "Closest_tree_dist_gen", "Gold standard", "closest_subset_v2", "closest_subset_v3"),])+
  #stat_summary(geom="line",  alpha=0.7)+
  stat_summary(geom="bar",  position = position_dodge(width=0.8),color="black",width=0.75)+
  theme_classic()+
  coord_cartesian(ylim=c(0.1, 1.0))+
  #scale_color_brewer(name = "",palette="Set1")
  #scale_color_brewer(name = "",palette="Dark2")
  
  #scale_fill_brewer(name = "", palette="Dark2", labels=c("closest_subset_v10" = "kf2d"))+
  #scale_color_brewer(name = "", palette="Dark2", labels=c("closest_subset_v10" = "kf2d"))+
  
  #scale_fill_manual(values = c(dark2_colors[2], dark2_colors[4], dark2_colors[5], dark2_colors[3],
  #                             dark2_colors[1], dark2_colors[6],dark2_colors[7], dark2_colors[8]),
  #                  name = "", labels=c("closest_subset_v10" = "kf2d"))+
  #scale_color_manual(values = c(dark2_colors[2], dark2_colors[4], dark2_colors[5], dark2_colors[3],
  #                              dark2_colors[1], dark2_colors[6],dark2_colors[7], dark2_colors[8]),
  #                  name = "", labels=c("closest_subset_v10" = "kf2d"))+
  
 # scale_color_manual(values = c(my_set1_colors[8], my_set1_colors[4], my_set1_colors[2], my_set1_colors[3], my_set1_colors[1],
#                                my_set1_colors[5], my_set1_colors[6],my_set1_colors[7] ),
#                     labels= c("DIAMOND 0.9.28","MEGAN 6.15.2", "LSHVec cami2", "PhyloPythiaS+ 1.4", "closest_subset_v10" = "kf2d", "Kraken 2.0.8-beta" ), name = "")+
#  scale_fill_manual(values = c(my_set1_colors[8], my_set1_colors[4], my_set1_colors[2], my_set1_colors[3], my_set1_colors[1],
#                                my_set1_colors[5], my_set1_colors[6],my_set1_colors[7] ),
 #                    labels= c("DIAMOND 0.9.28","MEGAN 6.15.2", "LSHVec cami2", "PhyloPythiaS+ 1.4", "closest_subset_v10" = "kf2d", "Kraken 2.0.8-beta" ), name = "")+
 
  scale_color_manual(values = c(#9BA8AEFF", dark2_colors[7], dark2_colors[5], dark2_colors[3], dark2_colors[4],
                                  dark2_colors[2], dark2_colors[6],dark2_colors[8], dark2_colors[5]),
                     labels= c("DIAMOND 0.9.28","MEGAN 6.15.2", "LSHVec cami2", "PhyloPythiaS+ 1.4", "closest_subset_v10" = "kf2vec", "Kraken 2.0.8-beta" ), name = NULL)+
  scale_fill_manual(values = c("#9BA8AEFF", dark2_colors[7], dark2_colors[5], dark2_colors[3], dark2_colors[4],
                                 dark2_colors[2], dark2_colors[6],dark2_colors[8], dark2_colors[5] ),
                    labels= c("DIAMOND 0.9.28","MEGAN 6.15.2", "LSHVec cami2", "PhyloPythiaS+ 1.4", "closest_subset_v10" = "kf2vec", "Kraken 2.0.8-beta" ), name = NULL)+
  
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="", label = c("Superkindom", "Phylum", "Class","Order", "Family", "Genus", "Species"))  +
  ylab("Accuracy")+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.1, hjust=1))+
  theme(#legend.position="bottom",
    legend.position=c(.42,-0.19),
        legend.margin=margin(0, 0,0,0), 
        #legend.text  = element_text(size = 8),
        legend.key.size = unit(0.5, "lines")
        #axis.text.x = element_text(size = 4)
        )+
#guides(color=guide_legend(ncol=2))
guides(fill=guide_legend(ncol=3))


ggsave("CAMI-D5-accuracy-line.pdf",width=3.8,height = 4.0)
getwd()



unique(amber_df$Tool)

stats_data = amber_df[amber_df$Tool=="closest_subset_v10",]
stats_data2 = amber_df[amber_df$Tool=="LSHVec cami2",]
stats_data$accuracy_seq
stats_data2$accuracy_seq
summary(c(stats_data$accuracy_seq[-length(stats_data$accuracy_seq)]))
summary(c(stats_data2$accuracy_seq[-length(stats_data2$accuracy_seq)]))


  ggplot(aes(x=rank, y=accuracy_seq, color = Tool, group = Tool),
         data=amber_df[!amber_df$Tool %in% c("Closest_tree_dist_phyl", "Closest_tree_dist_cls", "Closest_tree_dist_ordr", "Closest_tree_dist_fam", "Closest_tree_dist_gen", "Gold standard"),])+
    #data=amber_df[amber_df$Tool %in% c("closest_subset"),])+
    #geom_violin()+
    stat_summary()+
    #stat_summary(geom="line", aes(group=conditions, color = conditions))
    stat_summary(geom="line" , alpha = 0.7)+
    #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
    theme_classic()+
    coord_cartesian(ylim=c(0.55, 1.0))
  #ylim(0.75, 1.05)
  scale_color_brewer(name = "",palette="Set1")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="Contig length (KB)", label = c("(5-10]", "(10-20]", "(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-10000]"))  +
    ylab("Placement error")+
    theme(axis.text.x = element_text(angle = 10, vjust = 0.7, hjust=0.5))
  

  ggplot(aes(x=rank, y=accuracy_seq..unfiltered., color = Tool, group = Tool),
         data=amber_df[!amber_df$Tool %in% c("Closest_tree_dist_phyl", "Closest_tree_dist_cls", "Closest_tree_dist_ordr", "Closest_tree_dist_fam", "Closest_tree_dist_gen", "Gold standard"),])+
    #data=amber_df[amber_df$Tool %in% c("closest_subset"),])+
    #geom_violin()+
    stat_summary()+
    #stat_summary(geom="line", aes(group=conditions, color = conditions))
    stat_summary(geom="line" , alpha = 0.7)+
    #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
    theme_classic()
    #ylim(0.75, 1.05)
    scale_color_brewer(name = "",palette="Dark2")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="Contig length (KB)", label = c("(5-10]", "(10-20]", "(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-10000]"))  +
    ylab("Placement error")+
    theme(axis.text.x = element_text(angle = 10, vjust = 0.7, hjust=0.5))
  
  
  colnames(amber_df)
  
  ggplot(aes(x=rank, y=f1_score_per_seq..unfiltered., color = Tool, group = Tool),
         data=amber_df[!amber_df$Tool %in% c("Closest_tree_dist_phyl", "Closest_tree_dist_cls", "Closest_tree_dist_ordr", "Closest_tree_dist_fam", "Closest_tree_dist_gen", "Gold standard"),])+
    #data=amber_df[amber_df$Tool %in% c("closest_subset"),])+
    #geom_violin()+
    stat_summary()+
    #stat_summary(geom="line", aes(group=conditions, color = conditions))
    stat_summary(geom="line" , alpha = 0.7)+
    #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
    theme_classic()
    scale_color_brewer(name = "",palette="Dark2")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="Contig length (KB)", label = c("(5-10]", "(10-20]", "(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-10000]"))  +
    ylab("Placement error")+
    theme(axis.text.x = element_text(angle = 10, vjust = 0.7, hjust=0.5))
  


  

  
df1 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/Closest_tree_dist_sp/metrics_per_bin.tsv"),sep="\t",header=T)
df2 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/DIAMOND 0.9.28/metrics_per_bin.tsv"),sep="\t",header=T)
df3 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/Kraken 2.0.8-beta/metrics_per_bin.tsv"),sep="\t",header=T)
df4 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/LSHVec cami2/metrics_per_bin.tsv"),sep="\t",header=T)
df5 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/MEGAN 6.15.2/metrics_per_bin.tsv"),sep="\t",header=T)
df6 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/PhyloPythiaS+ 1.4/metrics_per_bin.tsv"),sep="\t",header=T)
df7 = read.csv( file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset/metrics_per_bin.tsv"),sep="\t",header=T)
df8 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v2/metrics_per_bin.tsv"),sep="\t",header=T)
df9 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/tax2tree/metrics_per_bin.tsv"),sep="\t",header=T)
df10 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v3/metrics_per_bin.tsv"),sep="\t",header=T)
df11 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v4/metrics_per_bin.tsv"),sep="\t",header=T)
df12 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v5/metrics_per_bin.tsv"),sep="\t",header=T)
df13 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v6/metrics_per_bin.tsv"),sep="\t",header=T)
df14 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v7/metrics_per_bin.tsv"),sep="\t",header=T)
df15 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v8/metrics_per_bin.tsv"),sep="\t",header=T)
df16 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v9/metrics_per_bin.tsv"),sep="\t",header=T)
df17 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v10/metrics_per_bin.tsv"),sep="\t",header=T)
#df18 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v10_rep/metrics_per_bin.tsv"),sep="\t",header=T)

df1["Tool"] = "Closest_tree_dist_sp"
df2["Tool"] = "DIAMOND 0.9.28"
df3["Tool"] = "Kraken 2.0.8-beta"
df4["Tool"] = "LSHVec cami2"
df5["Tool"] = "MEGAN 6.15.2"
df6["Tool"] = "PhyloPythiaS+ 1.4"
df7["Tool"] = "closest_subset"
df8["Tool"] = "closest_subset_v2"
df9["Tool"] = "tax2tree"
df10["Tool"] = "closest_subset_v3"
df11["Tool"] = "closest_subset_v4"
df12["Tool"] = "closest_subset_v5"
df13["Tool"] = "closest_subset_v6"
df14["Tool"] = "closest_subset_v7"
df15["Tool"] = "closest_subset_v8"
df16["Tool"] = "closest_subset_v9"
df17["Tool"] = "closest_subset_v10"
#df18["Tool"] = "closest_subset_v10_rep"


require(tidyr)

#install.packages("tidyverse")
require(tidyverse)

metrics_per_bin_df = rbind(df1, df2, df3, df4, df5, df6, df7, df9,  df11, df12,  df13, df14, df15, df16, df17)
colnames (metrics_per_bin_df)
nrow(metrics_per_bin_df)
metrics_per_bin_df$Taxonomic.rank  <- factor(metrics_per_bin_df$Taxonomic.rank, levels = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"))

metrics_per_bin_df%>%group_by(Taxonomic.rank, Tool, sample_id, Filtered) %>% mutate(Sample_Purity=mean(Purity..seq., na.rm = TRUE))

nrow(metrics_per_bin_df)

require(tidyverse)

my_set1_colors = RColorBrewer::brewer.pal(8, "Set1")
my_set1_colors[1]

metrics_per_bin_df %>%
  filter((!grepl("losest",Tool) | grepl("v10",Tool)) & Filtered == "False") %>%
  filter(Tool != "tax2tree") %>%
  filter(Tool != "DIAMOND 0.9.28") %>%
  #mutate(bin=cut(True.size..bp.,c(0,10^4,10^7,10^8.3,10^15))) %>%
  group_by(Taxonomic.rank,Tool) %>%
  dplyr::summarise(comp=sum(Completeness..seq.* True.size..bp.,na.rm = T),
                   purity=sum(Purity..seq.* True.size..bp.,na.rm = T),
                   t=sum( True.size..bp.),
                   n=n()) %>%
ggplot(aes(x=purity/t,y=comp/t, color = reorder(Tool,purity), shape=Taxonomic.rank,
           group = Tool))+
  geom_path(size=0.1)+
  geom_point(size=2,fill=NA)+
  facet_wrap(.~ifelse(Taxonomic.rank=="species","species","higher ranks"),scales="free")+
  #scale_shape_manual(values=c(21,25,5,15,17,16,19), name = "", labels=c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))+
  scale_shape_manual(values=c(21,25,5,15,17,16,19), name = "")+
   #scale_color_brewer(palette = "Dark2")+
  scale_color_manual(values = c(dark2_colors[7], dark2_colors[5], dark2_colors[3], dark2_colors[4],
                                dark2_colors[2], dark2_colors[6],"#9BA8AEFF", dark2_colors[5]),
                     labels= c("MEGAN 6.15.2", "LSHVec cami2", "PhyloPythiaS+ 1.4", "kf2vec", "Kraken 2.0.8-beta" ), name = "")+
  #coord_cartesian(ylim=c(0.,1))+
  #geom_smooth(method=lm)+
  xlab ("Purity")+
  ylab ("Completeness")+
  #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()+
  guides(
    color = guide_legend(title = NULL, byrow = T),
    size = guide_legend(title = NULL, byrow = T),
  ) +
  theme(
    legend.spacing.y = unit(1, "lines"),
    legend.margin = margin(0, 0, 0, 0), legend.box = "vertical", legend.title=element_blank()
  )
  theme(legend.position = c(0.85,0.3), legend.margin=margin(0,0,0,0),  
       axis.text.x = element_text(size = 1)
  )

getwd()
ggsave("CAMI-D5-purity_complt-line.pdf",width=6.2,height = 4.0)



m <- metrics_per_bin_df %>%
  filter((!grepl("losest",Tool) | grepl("v10",Tool)) & Filtered == "False") %>%
  filter(Tool != "tax2tree") %>%
  filter(Tool != "DIAMOND 0.9.28") %>%
  #mutate(bin=cut(True.size..bp.,c(0,10^4,10^7,10^8.3,10^15))) %>%
  group_by(Taxonomic.rank,Tool) %>%
  dplyr::summarise(comp=sum(Completeness..seq.* True.size..bp.,na.rm = T),
                   purity=sum(Purity..seq.* True.size..bp.,na.rm = T),
                   t=sum( True.size..bp.),
                   n=n()) 
m$x = m$purity/m$t
m$y = m$comp/m$t

m[m$Taxonomic.rank=="species",]


metrics_per_bin_df %>%
  filter((!grepl("losest",Tool) | grepl("v10",Tool)) & Filtered == "False") %>%
  filter(Tool != "tax2tree") %>%
  filter(Tool != "DIAMOND 0.9.28") %>%
  #mutate(bin=cut(True.size..bp.,c(0,10^4,10^7,10^8.3,10^15))) %>%
  group_by(Taxonomic.rank,Tool) %>%
  dplyr::summarise(comp=sum(Completeness..seq.* True.size..bp.,na.rm = T),
                   purity=sum(Purity..seq.* True.size..bp.,na.rm = T),
                   t=sum( True.size..bp.),
                   n=n()) %>%
  ggplot(aes(x=purity/t,y=comp/t, color = reorder(Tool,purity), shape=Taxonomic.rank,
             group = Tool))+
  geom_path(size=0.1)+
  geom_point(size=2,fill=NA)+
  facet_wrap(.~ifelse(Taxonomic.rank=="species","species","higher ranks"),scales="free")+
  scale_shape_manual(values=c(21,25,5,15,17,16,19))+
  scale_color_brewer(palette = "Dark2")+
  #coord_cartesian(ylim=c(0.,1))+
  #geom_smooth(method=lm)+
  #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()

ggplot(aes(x=Taxonomic.rank, y=Purity..seq., color = Tool,fill=Tool, group = Tool),
       data=metrics_per_bin_df[
         (!grepl("losest",metrics_per_bin_df$Tool)|grepl("v10",metrics_per_bin_df$Tool)) &
           metrics_per_bin_df$Tool != "tax2tree"
         & metrics_per_bin_df$Filtered == "False",])+
  #data=metrics_per_bin_df)+
  #data=amber_df[amber_df$Tool %in% c("closest_subset"),])+
  #geom_violin()+
  stat_summary()+
  stat_summary(geom="line", alpha=0.7)+
  #stat_summary(geom="bar" , alpha = 0.7, position = position_dodge())+
  coord_cartesian(ylim=c(0.,1))+
  #geom_smooth(method=lm)+
  #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()
  #ylim(0.996, 1.0005)
  #scale_color_brewer(name = "",palette="Dark2")
  
  ggplot(aes(x=Taxonomic.rank, y=Completeness..seq., fill=Tool,color = Tool, group = Tool),
         data=metrics_per_bin_df[ (!grepl("losest",metrics_per_bin_df$Tool)|grepl("v10",metrics_per_bin_df$Tool)) &
                                    metrics_per_bin_df$Tool != "tax2tree" &
                                    metrics_per_bin_df$Filtered == "False",])+
    #data=metrics_per_bin_df)+
    #data=amber_df[amber_df$Tool %in% c("closest_subset"),])+
    #geom_violin()+
    stat_summary()+
    #stat_summary(geom="line", aes(group=conditions, color = conditions))
    #stat_summary(geom="line" , alpha = 0.7)+
    stat_summary(geom="line", alpha=0.7)+
    #stat_summary(geom="bar" , alpha = 0.7, position = position_dodge())+
    #geom_smooth(method=lm)+
    #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
    theme_classic()
  scale_color_brewer(name = "",palette="Dark2")

colnames(metrics_per_bin_df)

metrics_per_bin_df$Purity..seq.

metrics_per_bin_df["gp"] = 0

ggplot(aes(x=Purity..seq., y=Completeness..seq., color = Tool, group = gp),
       data=metrics_per_bin_df[metrics_per_bin_df$Filtered == "False" & metrics_per_bin_df$Taxonomic.rank == "species",])+
  #facet_grid(.~Taxonomic.rank)+
  #data=metrics_per_bin_df)+
  #data=amber_df[amber_df$Tool %in% c("closest_subset"),])+
  #geom_violin()+
  #stat_summary(fun.data = "mean_cl_boot", linewidth = 2, size = 0.1)
  stat_summary(fun = mean, geom = "point")

  #stat_summary(geom="line", aes(group=conditions, color = conditions))
  stat_summary(geom="line" , alpha = 0.7)+
  #geom_smooth(method=lm)+
  #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()+
  scale_color_brewer(name = "",palette="Dark2")



################################################################
################################################################
################################################################

  # OBSOLETE (PREVIOUS VERSION OF PAPER)    

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



#pl_err=read.csv('pl_results_dev_queries.txt',sep=" ",header=T)
pl_err=read.csv('pl_results_dev_queries_ROUND2.txt',sep=" ",header=T)

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



err_epoch=read.csv('error_per_epoch.txt',sep=" ",header=T)

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
#EXTENDED QUERIES

pl_err2=read.csv('/Users/nora/Documents/ml_metagenomics/paper_fig1/pl_results_wol_extended_queries.txt',sep=" ",header=T)

head(pl_err2)




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
  theme(axis.text.x = element_text(angle = 30, vjust = 0.7, hjust=0.5))
#theme(axis.line = element_line(colour = "black"),
#      panel.grid.major = element_blank(),
#      panel.grid.minor = element_blank(),
#      panel.border = element_blank(),
#      panel.background = element_blank()) 

ggsave("dev_queries_multiple_conditions.pdf",width=5,height = 4)

##########################################################################


#EXTENDED QUERIES MULTIPLE GENES

pl_err2=read.csv('/Users/nora/Documents/ml_metagenomics/paper_fig1/pl_results_wol_extended_queries_depp_multgen.txt',sep=" ",header=T)

head(pl_err2)




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

  
ggplot(aes(color=cond, x=pl_error),
         data=pl_err2[pl_err2$cond %in% c("depp","depp_1","combined"), ])+
  stat_ecdf(alpha=0.5,size=1)+
  theme_classic()
    #data=pl_err)+
    #geom_plot(alpha=0.6)+
    #geom_bar(position = "dodge",stat = "summary",fun = "mean")+
    #geom_boxplot(alpha=0.6)+
  

ggplot(aes(fill=cond, x=cut(pl_error,c(0,1,4,6,11,40),include.lowest =TRUE,right = FALSE)),
       data=pl_err2[pl_err2$cond %in% c("depp","depp_1","combined"), ])+
  geom_bar(position = position_dodge())+
  annotate(x=3.5,y=120,label="Mean error:",geom="text")+
  geom_text(aes(x=2.5+as.numeric(cond)/2,y=100,label=round(`.`,2),color=cond),
            data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c("depp","depp_1","combined"), ],fun.aggregate = mean,value.var="pl_error"))+
  theme_classic()+
  #scale_fill_manual(name="", values = c("#252525", "#ca0020", "#0571b0" ), 
  #                   labels = c("Simulated", "Corrected", "Uncorrected"))
  #scale_fill_manual(name="", values = c("#2c7bb6", "#fdae61", "#d7191c" ))
  scale_fill_brewer(name = "",palette="Dark2", labels = c("kf2d", "DEPP_381", "DEPP_1"))+
  scale_color_brewer(name = "",palette="Dark2", labels = c("kf2d", "DEPP_381", "DEPP_1"))+
  ylab("Count")+
  xlab("Placement error")+
  #theme(axis.title.x = element_blank())+
  theme_classic()+
  theme(legend.position = c(.9,.8), legend.margin=margin(t = -0.5, unit='cm') )
  
ggsave("kf2d_vs_depp.pdf",width=5,height = 4)



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
##########################################################################
##########################################################################
# OBSOLETE (PREVIOUS VERSION OF PAPER)  
  
# Deft vs Depp 16s queries



pl_err2=read.csv('/Users/nora/Documents/ml_metagenomics/paper_fig1/pl_results_wol_extended_queries_depp_16s.txt',sep=" ",header=T)
  pl_err2
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
 
ggsave("dev_queries_multiple_conditions_16s_genes.pdf",width=5,height = 4)
  
nrow(pl_err2[pl_err2$cond %in% c("combined"), ])
nrow(pl_err2[pl_err2$cond %in% c("combined") & pl_err2$pl_error >4, ])
nrow(pl_err2[pl_err2$cond %in% c("combined") & pl_err2$pl_error ==0, ])
nrow(pl_err2[pl_err2$cond %in% c("depp") & pl_err2$pl_error ==0, ])
nrow(pl_err2[pl_err2$cond %in% c("depp_16s") & pl_err2$pl_error ==27, ])

  ggplot(aes(fill=cond, x=cut(pl_error,c(0,1,4,6,11,40),include.lowest =TRUE, right = FALSE)),
         data=pl_err2[pl_err2$cond %in% c("depp","depp_16s","combined"), ])+
    geom_bar(position = position_dodge())+
    annotate(x=3.5,y=40,label="Mean error:",geom="text")+
    geom_text(aes(x=2.5+as.numeric(cond)/2,y=33,label=round(`.`,2),color=cond),
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
    theme(legend.position = c(.9,.8), legend.margin=margin(t = -0.5, unit='cm') )
  
  ggsave("kf2d_vs_depp_16s.pdf",width=5,height = 4)
  
  
# For ppt  
#& pl_err2$pl_error <=1
ggplot(aes(fill=cond, x=cut(pl_error,c(0, 1, 40),include.lowest =TRUE,right = FALSE)),
         data=pl_err2[pl_err2$cond %in% c("depp","depp_16s","combined"), ])+
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
  
pl_err2[pl_err2$cond %in% c("depp","depp_16s","combined") & (pl_err2$pl_error <=1), ]

ggplot(aes(fill=cond, x=cut(pl_error,c(0, 1, 40),include.lowest =TRUE,right = FALSE)),
       data=pl_err2[pl_err2$cond %in% c("depp","depp_16s","combined") & (pl_err2$pl_error <1), ])+
  geom_bar(position = position_dodge())+
  annotate(x=1.0,y=90,label="Mean error:",geom="text")+
  geom_text(aes(x=0.3+as.numeric(cond)/3,y=84,label=round(`.`,2),color=cond),
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
  theme(legend.position = "none", legend.margin=margin(t = -0.5, unit='cm') )+
  annotate(x=0.7,y=30,label=c("kf2d"), geom="text")+
  annotate(x=1.0,y=30,label=c("DEPP_381"), geom="text")+
  annotate(x=1.3,y=30,label=c("DEPP_16s"), geom="text")

  ggsave("kf2d_vs_depp_16s_ppt_v2.pdf",width=5,height = 4)
  
  # Get the means from the un-normed data


  
ggplot(aes(fill=cond, x = cond, y=pl_error ),
         data=pl_err2[pl_err2$cond %in% c("depp","depp_16s","combined"), ])+
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

# Build plots for variabel number of genes with random selection
pl_err2=read.csv('/Users/nora/Documents/ml_metagenomics/paper_fig1/pl_results_wol_extended_queries_depp_multgen_random.txt',sep=" ",header=T)
  
head(pl_err2)
  
  
  
  
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
         data=pl_err2[pl_err2$cond %in% c("depp","depp_1","combined"), ])+
    stat_ecdf(alpha=0.5,size=1)+
    theme_classic()
  #data=pl_err)+
  #geom_plot(alpha=0.6)+
  #geom_bar(position = "dodge",stat = "summary",fun = "mean")+
  #geom_boxplot(alpha=0.6)+
  
  
  ggplot(aes(fill=cond, x=cut(pl_error,c(0,1,4,6,11,40),include.lowest =TRUE,right = FALSE)),
         data=pl_err2[pl_err2$cond %in% c("depp","depp_1","combined"), ])+
    geom_bar(position = position_dodge())+
    annotate(x=3.5,y=120,label="Mean error:",geom="text")+
    geom_text(aes(x=2.5+as.numeric(cond)/2,y=100,label=round(`.`,2),color=cond),
              data= dcast(cond~.,data=pl_err2[pl_err2$cond %in% c("depp","depp_1","combined"), ],fun.aggregate = mean,value.var="pl_error"))+
    theme_classic()+
    #scale_fill_manual(name="", values = c("#252525", "#ca0020", "#0571b0" ), 
    #                   labels = c("Simulated", "Corrected", "Uncorrected"))
    #scale_fill_manual(name="", values = c("#2c7bb6", "#fdae61", "#d7191c" ))
    scale_fill_brewer(name = "",palette="Dark2", labels = c("kf2d", "DEPP_381", "DEPP_1"))+
    scale_color_brewer(name = "",palette="Dark2", labels = c("kf2d", "DEPP_381", "DEPP_1"))+
        ylab("Count")+
    xlab("Placement error")+
    #theme(axis.title.x = element_blank())+
    theme_classic()+
    theme(legend.position = c(.9,.8), legend.margin=margin(t = -0.5, unit='cm') )
    
  
  ggsave("kf2d_vs_depp_random.pdf",width=5,height = 4)
  


head(pl_err2[pl_err2$cond %in% c("depp","depp_1","combined"),])
ggplot(aes(fill=cond, x=cut(pl_error,c(0,1,4,6,11,40),include.lowest =T,right = F)),
       data=pl_err2[pl_err2$cond %in% c("depp","depp_1","combined"), ])+
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
  


err_vs_prob=read.csv('err_vs_probability.txt',sep=" ",header=T)
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
prob_plot=read.csv('output_probabilities.txt',sep=" ",header=F)

prob_plot
ggplot(aes(x=V1, y=V2, color = V3),
       data=prob_plot)+
  #geom_bar(position = "dodge",stat = "identity")+
  geom_line(alpha = 0.9)

#er = median(total$V3.x)
#er
#er = median(total$V3.y)
#er

################################################################################
########### Updated software results (v39 training) ###########

pl_global=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_single_clade/all.pl_err',sep=" ",header=F)
pl_claded=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v20_8k_s28_FullGenomQueries_TRUE_CLADE/all.pl_err',sep=" ",header=F)
pl_chunked=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_Chunks_FullGenomQueries_TRUE_CLADE/all.pl_err',sep=" ",header=F)

pl_claded$V1 <- gsub('apples_input_di_mtrx_query_', '', pl_claded$V1)
pl_claded$V1 <- gsub('.pl_err', '', pl_claded$V1)
pl_claded <- pl_claded[,colSums(is.na(pl_claded))<nrow(pl_claded)]
pl_claded <- na.omit(pl_claded)


pl_chunked$V1 <- gsub('apples_input_di_mtrx_query_', '', pl_chunked$V1)
pl_chunked$V1 <- gsub('.pl_err', '', pl_chunked$V1)
pl_chunked <- pl_chunked[,colSums(is.na(pl_chunked))<nrow(pl_chunked)]
pl_chunked <- na.omit(pl_chunked)


colnames(pl_global) <- c("genome", "a", "pl_err", "b")
colnames(pl_claded) <- c("genome", "a", "pl_err", "b")
colnames(pl_chunked) <- c("genome", "a", "pl_err", "b")



pl_global['cond'] = "Global"
pl_claded['cond'] = "Clade"
pl_chunked['cond'] = "Chunk"

true_clade = read.csv('/Users/nora/Documents/ml_metagenomics/clade_targets.txt',sep=" ",header=T)


  
pl_global_merge = merge(pl_global, true_clade, by.x=c("genome"), by.y=c("genome"))
pl_claded_merge = merge(pl_claded, true_clade, by.x=c("genome"), by.y=c("genome"))
pl_chunked_merge = merge(pl_chunked, true_clade, by.x=c("genome"), by.y=c("genome"))


three_combo = rbind(pl_global_merge, pl_claded_merge, pl_chunked_merge)

nrow(three_combo)
colnames(three_combo)

head(three_combo)
dcast(data=three_combo,cond~.,value.var = "pl_err",fun.aggregate = mean)

ggplot(aes(x=reorder(clade,pl_err), y=pl_err, color = cond, group = interaction(cond)),
       data=three_combo[three_combo$cond !="Chunk",])+
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
  geom_rect(xmin = 14.5,xmax = 15.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  #scale_color_brewer(palette = my_palette,name="")+
  #scale_color_manual(name="", values = c( "#252525", "#ca0020", "#0571b0" ), 
                   #  labels = c("Global", "Local (true clade)", "Chunk (true clade)"))+
  #geom_hline(yintercept = 2.0, col = "blue", linetype='dotted', linewidth=0.5)+
  xlab("True clade")+
  ylab("Placement error")+
  #geom_violin()+
  theme_classic()+
  stat_summary(geom="line")+
  stat_summary(fun.data=function(x) mean_ci(x,ci=0.99),
               position=position_dodge(width=0.3))+
  coord_cartesian(ylim=c(0,5))+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = c(.2,.9),legend.direction = "vertical")

ggsave("pl_per_clade_true_violin_v3_updated_training.pdf",width=6.2,height = 4)








################################################################################



ggplot(aes(x=pl_error, color = cond),
       data=per_clade_error)+
  #scale_color_brewer(palette = my_palette,name="")+
  scale_color_manual(name="", values = c( "#ca0020", "#0571b0" ), 
                     labels = c("Global", "Local (true clade)"))+
  xlab("True clade")+
  ylab("Placement error")+
  #geom_violin()+
  theme_classic()+
  stat_ecdf(aes(group = interaction(cond,true_clade)),alpha=0.4,size=0.3)+
  stat_ecdf(aes(group = interaction(cond)),size=1.2)+
  coord_cartesian(ylim=c(0,1))+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = c(.7,.2),legend.direction = "vertical")


ggplot(aes(x=true_clade, y=pl_error, color = cond, group = true_clade),
       data=per_clade_error)+
  facet_wrap(~cond)+
  scale_color_manual(name="", values = c( "#ca0020", "#0571b0" ), 
                     labels = c("Global", "Local"),)+
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




#a=read.csv('kmercounts.txt',sep=" ",header=F)

#c=read.csv('kmer_pl_gpu_shared_500epochs_5000cols.csv',sep=",",header=T)

#data=c[c$epoch>500 & c$epoch<1000,]

# Plot error decay

#f2 <-ggplot(aes(x=epoch,y=log10(val), color=cond),data=c)+
#  geom_line(alpha=0.6)+
#  theme_classic()
#f2

#ggsave("kmer_pl_gpu_shared_1500epochs.pdf",width=5,height = 4)




d=read.csv('mer_pl_gpu_shared_model6.csv',sep=",",header=F)



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

d=read.csv('pl_error_model_v78.txt',sep=" ",header=F)

#f=read.csv('pl_error_model_16.txt',sep=" ",header=F)
#f=read.csv('pl_error_model_20_7kC_8000.txt',sep=" ",header=F)
#f=read.csv('pl_error_model17_8000_V3_shuffled.txt',sep=" ",header=F)



#f1=read.csv('pl_error_model_12_hiqual.txt',sep=" ",header=F)
f1=read.csv('pl_error_model_v78.txt',sep=" ",header=F)
f2=read.csv('pl_error_model_v77.txt',sep=" ",header=F)




#f2=read.csv('pl_error_model_21_7kC_8000_alpha_v1.txt',sep=" ",header=F)
f3=read.csv('pl_error_model_v76.txt',sep=" ",header=F)
f4=read.csv('pl_error_model_v75.txt',sep=" ",header=F)

f5=read.csv('pl_error_model_21_7kC_8000_lambda_v5.txt',sep=" ",header=F)
#f5b=read.csv('pl_error_model_21_7kC_8000_lambda_v5_excl200.txt',sep=" ",header=F)

f6=read.csv('pl_error_model_21_7kC_8000_lambda_v6.txt',sep=" ",header=F)
f7=read.csv('pl_error_model_21_7kC_8000_lambda_v7.txt',sep=" ",header=F)
f8=read.csv('pl_error_model_21_7kC_8000_lambda_v8.txt',sep=" ",header=F)
f9=read.csv('pl_error_model_21_7kC_8000_lambda_v9.txt',sep=" ",header=F)

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


l1=read.csv('pl_error_model_v78.txt',sep=" ",header=F)
l2=read.csv('pl_error_model_0.25_v78.txt',sep=" ",header=F)
l3=read.csv('pl_error_model_0.33_v78.txt',sep=" ",header=F)
l4=read.csv('pl_error_model_0.40_v78.txt',sep=" ",header=F)
l5=read.csv('pl_error_model_0.45_v78.txt',sep=" ",header=F)
l6=read.csv('pl_error_model_lambda_apples_v78.txt',sep=" ",header=F)


l7=read.csv('pl_error_model_lambda_apples_v78_f20_myapples.txt',sep=" ",header=F)
l8=read.csv('pl_error_model_lambda_apples_v78_f20_metinapples.txt',sep=" ",header=F)
l9=read.csv('pl_error_model_v78_f20.txt',sep=" ",header=F)
l10=read.csv('pl_error_model_v78_std_f20.txt',sep=" ",header=F)

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



df_last =read.csv('pl_error_model_apples_MEAN_v78_f0.txt',sep=" ", header=F)
df_mean =read.csv('pl_error_model_lambda_apples_MEAN_v78_f0_metinsapples.txt',sep=" ", header=F)

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

goodbad_df=read.csv('Good_bad_contamination.csv',sep=",",header=T)


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
