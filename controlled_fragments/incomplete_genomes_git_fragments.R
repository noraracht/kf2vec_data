

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


# Comparison v2 and v3 (entire range) of chunked genome queries on different models

# Uncladed, unchunked
#df_UU =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_Global_qChunks_v2_3_SINGLE_CLADE/all.pl_err',sep=" ",header=FALSE)
df_UU =read.csv('k7_v37_8k_s28_TrainClassf_10K_TOL_Global_qChunks_v2_3_SINGLE_CLADE_all.pl_err',sep=" ",header=FALSE)
df_UU ["claded"] = "Uncladed"
df_UU ["chunked"] = "Unchunked"
df_UU ["clsf"] = "na"
nrow(df_UU)
#df_UU = df_UU[!(is.na(df_UU$V2)), ]

# Uncladed, chunked MODEL
#df_BU =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_Global_Chunks_qChunks_v2_3_SINGLE_CLADE/all.pl_err',sep=" ",header=FALSE)
df_BU =read.csv('k7_v37_8k_s28_TrainClassf_10K_TOL_Global_Chunks_qChunks_v2_3_SINGLE_CLADE_all.pl_err',sep=" ",header=FALSE)
df_BU ["claded"] = "Uncladed"
df_BU ["chunked"] = "Chunked"
df_BU ["clsf"] = "na"
nrow(df_BU)

# Claded, unchunked - True clade
#df_UC_tr =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_Clades_s24_qChunks_v2_3_TRUE_CLADE/all.pl_err',sep=" ",header=FALSE)
df_UC_tr =read.csv('k7_v37_8k_s28_TrainClassf_10K_TOL_Clades_s24_qChunks_v2_3_TRUE_CLADE_all.pl_err',sep=" ",header=FALSE)
df_UC_tr ["claded"] = "Claded"
df_UC_tr ["chunked"] = "Unchunked"
df_UC_tr["clsf"] = "tr"
nrow(df_UC_tr)
df_UC_tr

# Claded, unchunked - Computed clade
#df_UC_clsf =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_Clades_s24_qChunks_v2_3_COMPUTED_CLADE/all.pl_err',sep=" ",header=FALSE)
df_UC_clsf =read.csv('k7_v37_8k_s28_TrainClassf_10K_TOL_Clades_s24_qChunks_v2_3_COMPUTED_CLADE_all.pl_err',sep=" ",header=FALSE)
df_UC_clsf ["claded"] = "Claded"
df_UC_clsf ["chunked"] = "Unchunked"
df_UC_clsf["clsf"] = "cmp"


# Claded, chunked - True clade
#df_BC_tr =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_qChunks_v2_3_TRUE_CLADE/all.pl_err',sep=" ",header=FALSE)
df_BC_tr =read.csv('k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_qChunks_v2_3_TRUE_CLADE_all.pl_err',sep=" ",header=FALSE)
df_BC_tr ["claded"] = "Claded"
df_BC_tr ["chunked"] = "Chunked"
df_BC_tr["clsf"] = "tr"

# Claded, chunked - Computed clade
#df_BC_clsf =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_qChunks_v2_3_COMPUTED_CLADE/all.pl_err',sep=" ",header=FALSE)
df_BC_clsf =read.csv('k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_qChunks_v2_3_COMPUTED_CLADE_all.pl_err',sep=" ",header=FALSE)
df_BC_clsf ["claded"] = "Claded"
df_BC_clsf ["chunked"] = "Chunked"
df_BC_clsf["clsf"] = "cmp"

# Claded, chunked - Computed clade - Previous version of classification without EXP
#df_BC_clsf_prev =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_Chunks_qChunks_v2_3_COMPUTED_CLADE/all.pl_err',sep=" ",header=FALSE)
df_BC_clsf_prev =read.csv('k7_v37_8k_s28_TrainClassf_10K_TOL_Chunks_qChunks_v2_3_COMPUTED_CLADE_all.pl_err',sep=" ",header=FALSE)
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
#c1 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_qChunks_v2_3_COMPUTED_CLADE/classes.out',header=TRUE, sep = '\t')
#c1 =read.csv('k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_qChunks_v2_3_COMPUTED_CLADE_classes.out',header=TRUE, sep = '\t')

######## SPLIT FILE AND WRITE INTO 4 FILES SO IT CAN FIT TO GITHUB
# Assume you have a large dataframe called `c1`
#n <- nrow(c1)
#split_size <- ceiling(n / 4)

# Split into 4 parts
#df_parts <- split(c1, rep(1:4, each = split_size, length.out = n))

# Write each part to a separate CSV file
#for (i in 1:4) {
#  write.csv(df_parts[[i]], paste0("k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_qChunks_v2_3_COMPUTED_CLADE_classes_part_", i, ".csv"), row.names = FALSE)
#}

# READ SPLIT CHUNKS

# Read all 4 parts
df_list <- lapply(1:4, function(i) read.csv(paste0("k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_qChunks_v2_3_COMPUTED_CLADE_classes_part_", i, ".csv")))

# Combine them back into one dataframe
#new_c1 <- do.call(rbind, df_list)
c1 <- do.call(rbind, df_list)

#all.equal(c1,new_c1)
#nrow(c1)
#nrow(new_c1)

###################################################


#c1 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_qChunks_v2_3_TRUE_CLADE/assigned_classes.out',header=TRUE, sep = '\t')
#c1 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_ChunksExp_qChunks_v2_3_TRUE_CLADE/single_classes.out',header=TRUE, sep = '\t')

d1 = df_fin[df_fin$greedy == "ng" & df_fin$cond %in% c("Claded - Chunked - cmp"),]

#d1[d1$V1=="G000017265_ng_n04_s44",]$V1
#c1[c1$genome=="G000017265_ng_n04_s44",]$genome

#t = read.csv('/Users/nora/Documents/ml_metagenomics/clade_targets.txt',header=TRUE, sep = ' ')
t = read.csv('clade_targets.txt',header=TRUE, sep = ' ')

cd1 = merge(x = c1, y = d1, by.y = "V1", by.x = "genome", by = NULL)

#cd1[cd1$genome=="G000017265_ng_n04_s44",]


cdt1 = merge(x = cd1, y = t, by.x = "genome.y", by.y = "genome", by = NULL)
head(cdt1)
nrow(cdt1[cdt1$top_class==cdt1$clade,])


# Claded, unchunked - Computed clade
#c2 =read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v37_8k_s28_TrainClassf_10K_TOL_Clades_s24_qChunks_v2_3_COMPUTED_CLADE/classes.out',header=TRUE, sep = '\t')
#c2 =read.csv('k7_v37_8k_s28_TrainClassf_10K_TOL_Clades_s24_qChunks_v2_3_COMPUTED_CLADE_classes.out',header=TRUE, sep = '\t')

######## SPLIT FILE AND WRITE INTO 4 FILES SO IT CAN FIT TO GITHUB
# Assume you have a large dataframe called `c2`
#n <- nrow(c2)
#split_size <- ceiling(n / 4)

# Split into 4 parts
#df_parts <- split(c2, rep(1:4, each = split_size, length.out = n))

# Write each part to a separate CSV file
#for (i in 1:4) {
#  write.csv(df_parts[[i]], paste0("k7_v37_8k_s28_TrainClassf_10K_TOL_Clades_s24_qChunks_v2_3_COMPUTED_CLADE_classes_part_", i, ".csv"), row.names = FALSE)
#}

# READ SPLIT CHUNKS

# Read all 4 parts
df_list <- lapply(1:4, function(i) read.csv(paste0("k7_v37_8k_s28_TrainClassf_10K_TOL_Clades_s24_qChunks_v2_3_COMPUTED_CLADE_classes_part_", i, ".csv")))

# Combine them back into one dataframe
#new_c2 <- do.call(rbind, df_list)
c2 <- do.call(rbind, df_list)

#all.equal(c2,new_c2)
#nrow(c2)
#nrow(new_c2)

##############################################






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


