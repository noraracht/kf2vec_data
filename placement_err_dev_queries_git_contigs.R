

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


