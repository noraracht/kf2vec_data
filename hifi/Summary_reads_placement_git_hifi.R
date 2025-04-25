
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



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd('/Users/admin/Documents/support')
getwd()


# Summary of HiFi reads
#summary_HiFi=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v40_8k_s28All_HiFi_Reads_p_all/summary_hifi_six_models.pl_error', sep=" ",header=T)

# Splitting datarfame into 2 so it can be uploaded to github with the script
#summary_HiFi=read.csv('summary_hifi_six_models.pl_error', sep=" ",header=T)
#colnames (summary_HiFi)
#tail(summary_HiFi)
#unique(summary_HiFi$model)
#nrow(summary_HiFi)

#part1_models = unique(summary_HiFi$model)[1:5]
#part2_models = unique(summary_HiFi$model)[6:10]

#part1 = summary_HiFi[summary_HiFi$model %in%(part1_models),]
#write.csv(part1, "summary_hifi_six_models.pl_error_part1.csv", row.names = FALSE)

#part2 = summary_HiFi[summary_HiFi$model %in%(part2_models),]
#write.csv(part2, "summary_hifi_six_models.pl_error_part2.csv", row.names = FALSE)

part1_df <- read.csv("summary_hifi_six_models.pl_error_part1.csv")
part2_df <- read.csv("summary_hifi_six_models.pl_error_part2.csv")
#new_summary_HiFi = rbind(part1_df, part2_df)
summary_HiFi== rbind(part1_df, part2_df)

#all.equal(summary_HiFi,new_summary_HiFi)

summary_HiFi$model[summary_HiFi$model == "ChunkExp"] <- "Claded_chunked"
summary_HiFi$model[summary_HiFi$model == "Global"] <- "Uncladed_unchunked"
summary_HiFi$model[summary_HiFi$model == "Global_Chunk"] <- "Uncladed_chunked"
summary_HiFi$model[summary_HiFi$model == "Full"] <- "Claded_unchunked"
summary_HiFi$model[summary_HiFi$model == "ChunkExp_TrueClade"] <- "Claded_chunked_tr"
summary_HiFi$model[summary_HiFi$model == "Full_TrueClade"] <- "Claded_unchunked_tr"


unique(summary_HiFi$model)

summary_HiFi$claded[summary_HiFi$model=="Claded_chunked"] <- "Claded"
summary_HiFi$chunked[summary_HiFi$model=="Claded_chunked"] <- "Chunked"
summary_HiFi$clsf[summary_HiFi$model=="Claded_chunked"] <- "cmp"


summary_HiFi$claded[summary_HiFi$model=="Claded_chunked_tr"] <- "Claded"
summary_HiFi$chunked[summary_HiFi$model=="Claded_chunked_tr"] <- "Chunked"
summary_HiFi$clsf[summary_HiFi$model=="Claded_chunked_tr"] <- "tr"

summary_HiFi$claded[summary_HiFi$model=="Claded_unchunked"] <- "Claded"
summary_HiFi$chunked[summary_HiFi$model=="Claded_unchunked"] <- "Unchunked"
summary_HiFi$clsf[summary_HiFi$model=="Claded_unchunked"] <- "cmp"

summary_HiFi$claded[summary_HiFi$model=="Claded_unchunked_tr"] <- "Claded"
summary_HiFi$chunked[summary_HiFi$model=="Claded_unchunked_tr"] <- "Unchunked"
summary_HiFi$clsf[summary_HiFi$model=="Claded_unchunked_tr"] <- "tr"

summary_HiFi$claded[summary_HiFi$model=="Uncladed_unchunked"] <- "Uncladed"
summary_HiFi$chunked[summary_HiFi$model=="Uncladed_unchunked"] <- "Unchunked"
summary_HiFi$clsf[summary_HiFi$model=="Uncladed_unchunked"] <- "na"

summary_HiFi$claded[summary_HiFi$model=="Uncladed_chunked"] <- "Uncladed"
summary_HiFi$chunked[summary_HiFi$model=="Uncladed_chunked"] <- "Chunked"
summary_HiFi$clsf[summary_HiFi$model=="Uncladed_chunked"] <- "na"


summary_HiFi$cond <- paste(summary_HiFi$claded, "-", summary_HiFi$chunked, "-", summary_HiFi$clsf)




ggplot(aes(x=cut(length_adj/10^3,c(0.1, 0.5, 1, 2, 4, 8,16,32,64, 128, 1000)*10), y=pl_err, color = model),
      # data=summary_HiFi)+
  #data=summary_HiFi[summary_HiFi$length_adj>5000 & summary_HiFi$model %in% c("Chunk", "Full", "Global", "ChunkExp", "Global_Chunk"),])+
  #data=summary_HiFi[summary_HiFi$length_adj>5000 & summary_HiFi$model %in% c("Global", "ChunkExp"),])+
  data=summary_HiFi[summary_HiFi$length_adj>5000 & summary_HiFi$model %in% c("Claded_chunked", "Uncladed_unchunked", "Uncladed_chunked", "Claded_unchunked", "Claded_chunked_tr", "Claded_unchunked_tr"),])+
  #data=summary_HiFi[summary_HiFi$parent_genome %in% chunk_queries & summary_HiFi$top_class==summary_HiFi$true_clade,])+
  #data=summary_HiFi[summary_HiFi$top_class==summary_HiFi$true_clade,])+
  #data=summary_ONT[summary_ONT$parent_genome %in% chunk_queries,])+
  #geom_violin()+
  #facet_wrap(~top_class==true_clade)+
  stat_summary(geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8), position = position_dodge())+
  #stat_summary(geom="line",color="blue",aes(group="1"))+
  stat_summary(geom="line",aes(group=model, color = model))+
  stat_summary(alpha = 0.3)+
  theme_classic()+
  geom_hline(yintercept = 3)+
coord_cartesian(ylim=c(0, 15))
summary_HiFi

ggplot(aes(x=cut(length_adj/10^3,c(0.1, 0.5, 1, 2, 4, 8,16,32,64, 128, 1000)*10), 
           fill = top_class==true_clade),
       data=summary_HiFi)+
  facet_wrap(~model)+
  #stat_summary(geom="line",color="blue",aes(group="1"))+
  geom_bar(aes(y=(..count..)/sum(..count..)),stat="count")+
  theme_classic()



# Chunk vs Chunk_EXP for training classifier

#chunk_clsf<-summary_HiFi[summary_HiFi$model %in% c("Chunk"),]
chunk_clsf<-summary_HiFi[summary_HiFi$model %in% c("Claded_unchunked"),]
chunkExp_clsf<- summary_HiFi[summary_HiFi$model %in% c("Claded_chunked"),]


head(summary_HiFi)
for(x in c(0, 10000, 20000, 30000, 40000, 80000)){
  print(x)
  
  x1 = nrow(chunk_clsf[chunk_clsf$top_class==chunk_clsf$true_clade & chunk_clsf$length_adj >= x,]) /nrow(chunk_clsf[chunk_clsf$length_adj >= x,])
  x2 = nrow(chunkExp_clsf[chunkExp_clsf$top_class==chunkExp_clsf$true_clade & chunkExp_clsf$length_adj >= x,]) /nrow(chunkExp_clsf[chunkExp_clsf$length_adj >= x,])
  #combo_chunk_df = 
  print (x1)
  print (x2)
  chunk_clsf$clsf_acc = ifelse (chunk_clsf$length_adj >= x, x1, chunk_clsf$clsf_acc )
  chunkExp_clsf$clsf_acc = ifelse (chunkExp_clsf$length_adj >= x, x2, chunkExp_clsf$clsf_acc )
}

unique(chunk_clsf$clsf_acc)

combo_classif <- rbind(chunk_clsf, chunkExp_clsf)
head(combo_classif)
r1 = combo_classif[combo_classif$cond=="Claded - Unchunked - cmp",]

nrow(r1[r1$top_class==r1$true_clade,])/nrow(r1)

ggplot(aes(x=cut(length_adj/10^3,c(0, 1, 2, 3, 4, 8, 2000)*10), y=clsf_acc),
       #data=pl_per_contig[pl_per_contig$model %in% c("ChunkExp", "Full", "Global", "Global_Chunk") & pl_per_contig$length_adj > 5000,])+
       data=combo_classif)+
  #data=pl_per_contig)+
  #stat_summary(geom="crossbar", fun.data = function(x) median_hilow_(x,ci=0.8), aes(group=model, color = model))+
  #stat_summary(geom="line",color="blue",aes(group="1"))+
  stat_summary(geom="line",aes(group=model, color = model))+
  stat_summary(alpha = 0.3,aes(group=model, color = model))+
  theme_classic()+
  geom_hline(yintercept = 3)+
  coord_cartesian(ylim=c(0.8, 1))+
  ylab("Classification accuracy")+
  xlab("Read length (KB)")+
  theme(legend.position = c(0.85,0.3), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )+
  scale_colour_brewer(palette = "Dark2", name="", labels = c("Chunked", "Unchunked"),)
ggsave("HiFi-D4-classified-line.pdf",width=4.8,height = 4.0)


ggplot(aes(x=cut(length_adj/10^3,c(0.1, 0.5, 1, 2, 4, 8,16,32,64, 128, 1000)*10), y=pl_err, color = model),
       # data=summary_HiFi)+
       #data=summary_HiFi[summary_HiFi$length_adj>5000 & summary_HiFi$model %in% c("Chunk", "Full", "Global", "ChunkExp", "Global_Chunk"),])+
       #data=summary_HiFi[summary_HiFi$length_adj>5000 & summary_HiFi$model %in% c("Global", "ChunkExp"),])+
       data=summary_HiFi[summary_HiFi$length_adj>5000 & summary_HiFi$model %in% c("Claded_chunked", "Chunk"),])+
  #data=summary_HiFi[summary_HiFi$parent_genome %in% chunk_queries & summary_HiFi$top_class==summary_HiFi$true_clade,])+
  #data=summary_HiFi[summary_HiFi$top_class==summary_HiFi$true_clade,])+
  #data=summary_ONT[summary_ONT$parent_genome %in% chunk_queries,])+
  #geom_violin()+
  #facet_wrap(~top_class==true_clade)+
  stat_summary(geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8), position = position_dodge(),aes(group=model, color = model))+
  #stat_summary(geom="line",color="blue",aes(group="1"))+
  stat_summary(geom="line",aes(group=model, color = model))+
  stat_summary(alpha = 0.3,aes(group=model, color = model))+
  theme_classic()+
  geom_hline(yintercept = 3)+
  coord_cartesian(ylim=c(0, 15))

summary_HiFi

ggplot(aes(x=cut(length_adj/10^3,c(0.1, 0.5, 1, 2, 4, 8,16,32,64, 128, 1000)*10), 
           fill = top_class==true_clade),
       data=summary_HiFi)+
  facet_wrap(~model)+
  #stat_summary(geom="line",color="blue",aes(group="1"))+
  geom_bar(aes(y=(..count..)/sum(..count..)),stat="count")+
  theme_classic()


unique(summary_HiFi$chunked)


head(summary_HiFi[(summary_HiFi$clsf!="tr" |summary_HiFi$chunked=="Chunked") | (summary_HiFi$chunked=="Unchunked"),])

head(summary_HiFi)
ggplot(aes(x=cut(length_adj/10^3,c(0.4, 0.6, 0.8, 1, 1.4, 2, 6)*10), y=pl_err, color = chunked,linetype=paste(claded,ifelse(clsf=="tr","(True)","")),group=cond),
       data=summary_HiFi[!is.na(summary_HiFi$clsf) & !is.na(summary_HiFi$chunked) & !is.na(summary_HiFi$claded)  & summary_HiFi$clsf %in% c("cmp" ,"na")
      ,])+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  stat_summary(geom="line")+
  #geom_violin()+
  theme_classic()+
  stat_summary( alpha = 0.7,geom="point")+
  stat_summary( alpha = 0.7,geom="errorbar",width=0.1,linetype=1,size=0.3)+
  scale_linetype(name="", labels = c("Cladded", "Uncladded"))+
  #coord_cartesian(ylim=c(0,30))+
  #scale_x_discrete(name="Chunk length (KB)", label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  scale_colour_brewer(palette = "Dark2", name="")+
  #facet_wrap(~claded)+
  #guides(col= guide_legend(title= "Conditions"))+
  #scale_fill_discrete(name = "New Legend Title")
  #guides(fill=guide_legend(title="k"))
  theme(legend.position = c(0.8,0.8), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  ) +
  geom_bar (aes(x=cut(length_adj/10^3,c(0.4, 0.6, 0.8, 1, 1.4, 2, 6)*10),y=..count../10000/3.5),
            stat="count", size=.1, fill="#69b3a2", color="black", alpha=.4,inherit.aes = FALSE) + 
  scale_y_continuous(
    
    # Features of the first axis
    name = "Placement error",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*10000*3.5, name="Count",breaks = c(0,50000,100000))
  ) + xlab("Read length (KB)")

ggsave("HiFi-D4-line.pdf",width=4.8,height = 4.0)
getwd()

summary_HiFi %>% 
  mutate(bin = cut(length_adj/10^3,c(0.4, 0.6, 0.8, 1, 1.4, 2, 6)*10)) %>% 
  filter(cond %in% c("Claded - Chunked - cmp","Claded - Unchunked - cmp")) %>%
  group_by(bin,cond) %>%
  dplyr::summarise(n=n(),correct=sum(top_class==true_clade)) %>%
ggplot(aes(x=bin, y=correct/n,
           color=cond,group=cond))+
  geom_line( )+ 
  theme_classic()+ 
  geom_point(alpha=0.5)+
  coord_cartesian(ylim=c(0.7, 1))+
  ylab("Classification accuracy")+
  xlab("Read length (KB)")+
  theme(legend.position = c(0.85,0.3), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )+
  scale_colour_brewer(palette = "Dark2", name="", labels = c("Chunked", "Unchunked"),)
ggsave("HiFi-D4-classified-line.pdf",width=4.8,height = 4.0)
  

#scale_colour_brewer(palette = "Dark2", name="")

#ggsave("HiFi-D4-line.pdf",width=6.5,height = 4.5)

