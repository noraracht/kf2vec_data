

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
#install.packages("cowplot")
library("cowplot")

#install.packages("patchwork")
library(patchwork)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd('/Users/admin/Documents/support')
getwd()



pl_err=read.csv('cafe_vs_us_pl_error.txt',sep="\t",header=T)

corr_stats =read.csv('corr_stats_cafe_D2star.txt',sep="\t",header=T)
dist_df =read.csv('distances_subset_cafe_D2star.txt',sep="\t",header=T)
head(pl_err)

pl_err

corr_stats

#fin_df = rbind(pr_A1_grouped, pr_A01_mean)
#fin_df2 = cbind(pr_A1_grouped, pr_A01_mean)
#total <- merge(pr_A1_grouped, pr_A01_median, by = c("V1"))



# Compute mean
tmp1 <- pl_err[pl_err$method == "ours",]
mean(as.numeric(unlist(tmp1["true_clade_pl_err"])))

tmp1 <- pl_err[(pl_err$method == "D2star" & pl_err$cond == "regscaled"),]
mean(as.numeric(unlist(tmp1["true_clade_pl_err"])))

tmp1 <- pl_err[(pl_err$method == "D2shepp" & pl_err$cond == "regscaled"),]
mean(as.numeric(unlist(tmp1["true_clade_pl_err"])))

tmp1 <- pl_err[(pl_err$method == "JS" & pl_err$cond == "regscaled"),]
mean(as.numeric(unlist(tmp1["true_clade_pl_err"])))

tmp1 <- pl_err[(pl_err$method == "Co-phylog" & pl_err$cond == "regscaled"),]
mean(as.numeric(unlist(tmp1["true_clade_pl_err"])))


tmp1 <- pl_err[(pl_err$method == "Cosine" & pl_err$cond == "regscaled"),]
mean(as.numeric(unlist(tmp1["true_clade_pl_err"])))

pl_err2 <- pl_err[! (pl_err$cond == "100scaled"),]


dark2_colors <-brewer.pal(n = 8, name = "Dark2")
dark2_colors

set1_colors <-brewer.pal(n = 8, name = "Set1")
set1_colors

ggplot(aes(x=method, y=true_clade_pl_err, color=cond, fill = cond),
       data=pl_err2[!(pl_err2$method =="Co-phylog" | pl_err2$method =="JS"), ])+
       #data=pl_err)+
  #geom_plot(alpha=0.6)+
  geom_bar(position = "dodge",stat = "summary",fun = "mean", colour="black", alpha=0.9)+
  #geom_boxplot(alpha=0.6)+
  scale_x_discrete(labels = c("Cosine", "CVTree", "D2Shepp", "D2Star", "Eu", "MA", "kf2d"), )+
  #scale_color_manual(values = c( "#ca0020", "#0571b0" ))+
  #scale_fill_manual(values = c( "#ca0020", "#0571b0" ))+
  scale_color_manual(values = c( "#66A61E", "#E6AB02"  ))+
  scale_fill_manual(values = c( "#66A61E","#E6AB02"  ))+
  stat_summary(color = "black")+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  #geom_abline(color="red")+
  theme_classic()+
#scale_y_continuous(limits = c(0.0000, 0.001))+
 #xlab("")+
 ylab("Placement error")+
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 8))
  

ggsave("cafe_metrics_pl_error.pdf",width=3.8,height = 4)


mean(pl_err2[pl_err2$method =="Cosine", ]$true_clade_pl_err)
mean(pl_err2[pl_err2$method =="Eu", ]$true_clade_pl_err)
mean(pl_err2[pl_err2$method =="Ma", ]$true_clade_pl_err)
mean(pl_err2[pl_err2$method =="D2shepp", ]$true_clade_pl_err)
mean(pl_err2[pl_err2$method =="D2star", ]$true_clade_pl_err)
mean(pl_err2[pl_err2$method =="CVtree", ]$true_clade_pl_err)
mean(pl_err2[pl_err2$method =="ours", ]$true_clade_pl_err)





# For ppt
unique(pl_err2$method)
pl_err2

mean(pl_err2[pl_err2$method == "ours",]$true_clade_pl_err)

pl_err2$cond2 <- ""
pl_err2$cond2 <- ifelse(pl_err2$method %in% c("Co-phylog", "JS", 'Cosine', 'Eu', 'Ma'), "Raw", ifelse(pl_err2$method %in% c("CVtree","D2star", "D2shepp"), "Adjusted", "Trained"))
pl_err2$cond2 <- factor(pl_err2$cond2, levels = c("Raw", "Adjusted", "Trained"))
pl_err2$method <- factor(pl_err2$method, levels = c("Co-phylog", "JS", "Cosine", "Eu", "Ma","CVtree", "D2star", "D2shepp",  "ours"))
unique(pl_err2$method)


pl_err2[pl_err2$method =="ours", ]

ggplot(aes(x=method, y=true_clade_pl_err, color = cond2, fill= cond2),
       data=pl_err2[!(pl_err2$method =="Co-phylog" | pl_err2$method =="JS"), ])+
  #data=pl_err)+
  #geom_plot(alpha=0.6)+
  geom_bar(position = "dodge",stat = "summary",fun = "mean", colour="black", alpha=0.9)+
  geom_text(aes(label = after_stat(sprintf("%.1f", round(y, digits = 3)))), stat = "summary", fun = "mean", vjust = 2.5, colour = "black" ) +
  #geom_boxplot(alpha=0.6)+
  #stat_summary(aes(method=
  #                   ifelse(method %in% c( "Eu", "Ma", "Cosine"),"Hypothetical", ifelse(method %in% c( "CVtree","D2star", "D2shepp"),"Hypothetical2", "infection"))))+
  
  scale_x_discrete(labels = c('Cosine', 'Eu', 'Ma',"CVtree", expression(italic("D")[italic("2")]^italic("*")), expression(italic("D")[italic("2")]^italic("S")), 'kf2vec' ))+
  #scale_color_manual(values = c( "#ca0020", "#0571b0" ))+
  #scale_fill_manual(values = c( "#ca0020", "#0571b0" ))+
  scale_color_manual(values = c( "#E6AB02", "#66A61E","#A6761D"  ))+
  scale_fill_manual(values = c(  "#A6761D", "#E6AB02", "#66A61E"))+
  stat_summary(color = "black")+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  #geom_abline(color="red")+
  theme_classic()+
  #scale_y_continuous(limits = c(0.0000, 0.001))+
  #xlab("")+
  ylab("Placement error")+
  #theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 10))
  theme(legend.position = c(0.85,0.9), legend.margin=margin(0,0,0,0),legend.direction="vertical", legend.title=element_blank(),
        axis.title.x = element_blank())+
  guides(fill = guide_legend(reverse = FALSE))
#, axis.text.x = element_text(size = 8)


ggsave("cafe_metrics_pl_error_v2.pdf",width=3.8,height = 4)
#ggsave("cafe_metrics_pl_error_v2.pdf",width=4.8,height = 4.0)
getwd()



ggplot(aes(x=method, y=true_clade_pl_err, color = cond2, fill= cond2),
       data=pl_err2[!(pl_err2$method =="Co-phylog" | pl_err2$method =="JS"), ])+
  #data=pl_err)+
  #geom_plot(alpha=0.6)+
  geom_bar(position = "dodge",stat = "summary",fun = "mean", colour="black", alpha=0.9)+
  #geom_boxplot(alpha=0.6)+
  #stat_summary(aes(method=
  #                   ifelse(method %in% c( "Eu", "Ma", "Cosine"),"Hypothetical", ifelse(method %in% c( "CVtree","D2star", "D2shepp"),"Hypothetical2", "infection"))))+
  
  scale_x_discrete(labels = c('Cosine', 'Eu', 'Ma',"CVtree","D2Star", "D2Shepp", 'kf2d' ))+
  #scale_color_manual(values = c( "#ca0020", "#0571b0" ))+
  #scale_fill_manual(values = c( "#ca0020", "#0571b0" ))+
  scale_color_manual(values = c( "#E6AB02", "#66A61E","#A6761D"  ))+
  scale_fill_manual(values = c(  "#A6761D", "#E6AB02", "#66A61E"))+
  stat_summary(color = "black")+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  #geom_abline(color="red")+
  theme_classic()+
  #scale_y_continuous(limits = c(0.0000, 0.001))+
  #xlab("")+
  ylab("Placement error")+
  #theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 10))
  theme(legend.position = c(0.85,0.9), legend.margin=margin(0,0,0,0),legend.direction="vertical", legend.title=element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 10))+
  guides(fill = guide_legend(reverse = FALSE))


ggsave("cafe_metrics_pl_error_ppt.pdf",width=5,height = 4)




ggplot(aes(x=method, y=true_clade_pl_err, color = cond2, fill= cond2),
       data=pl_err2[!(pl_err2$method =="Co-phylog" | pl_err2$method =="JS"), ])+
  #data=pl_err)+
  #geom_plot(alpha=0.6)+
  geom_bar(position = "dodge",stat = "summary",fun = "mean", colour="black", alpha=0.9)+
  #geom_boxplot(alpha=0.6)+
  #stat_summary(aes(method=
  #                   ifelse(method %in% c( "Eu", "Ma", "Cosine"),"Hypothetical", ifelse(method %in% c( "CVtree","D2star", "D2shepp"),"Hypothetical2", "infection"))))+
  
  scale_x_discrete(labels = c('Cosine', 'Eu', 'Ma',"CVtree","D2Star", "D2Shepp", 'kf2d' ))+
  #scale_color_manual(values = c( "#ca0020", "#0571b0" ))+
  #scale_fill_manual(values = c( "#ca0020", "#0571b0" ))+
  scale_color_manual(values = c( "#E6AB02", "#66A61E","#A6761D"  ))+
  scale_fill_manual(values = c(  "#A6761D", "#E6AB02", "#66A61E"))+
  stat_summary(color = "black")+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  #geom_abline(color="red")+
  theme_classic()+
  #scale_y_continuous(limits = c(0.0000, 0.001))+
  #xlab("")+
  ylab("Placement error")+
  #theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 10))
  theme(legend.position = c(0.85,0.9), legend.margin=margin(0,0,0,0),legend.direction="vertical", legend.title=element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 8.5))+
  guides(fill = guide_legend(reverse = FALSE))

ggsave("cafe_metrics_pl_error_ppt_v2.pdf",width=3.8,height = 4)



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

ggsave("cafe_vs_us_correlation.pdf",width=3.8,height = 4)

#axis.title.y = element_text(size = 10), 
  #facet_wrap(~ portion)+
  #geom_boxplot(alpha=0.6)+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_point(aes(y=v56,color="vs 56"))+
  #geom_point(aes(y=v59,color="vs 59"))+
  #geom_abline(color="red")+
  #theme_classic()

pl_err=read.csv('cafe_vs_us_pl_error.txt',sep="\t",header=T)
corr_stats =read.csv('corr_stats_cafe_D2star.txt',sep="\t",header=T)
dist_df =read.csv('distances_subset_cafe_D2star.txt',sep="\t",header=T)



corr_stats
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

dist_df$portion <- factor(dist_df$portion, levels = c("backbone", "queries"),
                  labels = c("Backbone", "Queries")
)


dist_df$method <- factor(dist_df$method, levels = c("D2star", "us"), 
                         labels = c("D2Star", "kf2d"))


#scale_x_discrete(labels = c('Cosine', 'Eu', 'Ma',"CVtree", expression(italic("D")[italic("2")]^italic("*")), expression(italic("D")[italic("2")]^italic("S")), 'kf2d' ))+

ggplot(aes(x=true/100, y=dist/100, color = method),
       data=dist_df[dist_df$method == "kf2d" | dist_df$method == "D2Star",])+
  #geom_bar(position = "dodge",stat = "identity")+
  geom_point(alpha = 0.4, size = 0.5)+
  facet_wrap(~ portion)+
  geom_abline(linetype=2,size=1.4)+
  theme_classic()+
  stat_smooth(method="lm", se=FALSE,size=1.4)+
  scale_color_brewer(palette = "Dark2",name="", labels = c(expression(italic("D")[italic("2")]^italic("*")), 'kf2vec' ))+
  #geom_boxplot(alpha=0.6)+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_abline(color="red")+
#scale_y_continuous(limits = c(0.0000, 0.001))
 xlab("True distance")+
 ylab("Estimated distance")+
  coord_cartesian(xlim=c(0,4.5),ylim=c(0,4.5))+
theme(legend.position = c(.9,.1), legend.margin=margin(t = -0.5, unit='cm'), )
#ggsave("cafe_vs_us_true_est_distance.pdf",width=5,height = 4)
ggsave("cafe_vs_us_true_est_distance.pdf",width=6.2,height = 4)



#new.lab <- as_labeller(c(kf2d="kf2d", D2Star="italic(D)[italic(2)]^italic('*')" ), label_parsed)
#new.lab <- as_labeller(c(kf2d="kf2d", D2Star=label_bquote("Q") ), label_parsed)



big_plot<-ggplot(aes(x=true/100, y=dist/100, color = portion),
       data=dist_df[dist_df$method == "kf2d" | dist_df$method == "D2Star",])+
  #geom_bar(position = "dodge",stat = "identity")+
  geom_abline(linetype=2,size=0.9)+
  stat_smooth(method="lm",se=FALSE,size=1.4, formula=y~x+0)+
  geom_point(alpha = 0.3, size = 0.5)+
  facet_wrap(~ method, labeller = as_labeller(c(kf2d="kf2vec", D2Star="italic(D)[italic(2)]^italic('*')" ), label_parsed))+
  theme_classic()+
  scale_color_brewer(palette = "Dark2",name="")+
  #scale_color_manual(name="", values = c(dark2_colors[3],dark2_colors[4]))+
  #geom_boxplot(alpha=0.6)+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_abline(color="red")+
  #scale_y_continuous(limits = c(0.0000, 0.001))
  xlab("True distance")+
  ylab("Estimated distance")+
  coord_cartesian(xlim=c(0,4.5),ylim=c(0,4.5))+
  theme(legend.position = c(.1,.9), legend.margin=margin(t = -0.5, unit='cm'), )+
  #stat_cor(aes(label = after_stat(r.label)),label.y = 0.7,label.x = 3.6, p.digits = 2,label.x.npc  = "right",size=3)
stat_cor(aes(label = after_stat(r.label)),method = "pearson", label.y = c(0.4, 1.1) ,label.x = 3.4, p.digits = 2,size=3.3, show.legend=FALSE)+
stat_cor(aes(label = after_stat(p.label)),method = "pearson", label.y = c(0.1, 0.8) ,label.x = 3.4, p.digits = 2,size=3.3, show.legend=FALSE)

  #theme(legend.key = element_blank(), legend.text = element_blank())+
    #guides(color = guide_legend(override.aes = list(linetype = 0, shape = "")))
big_plot
ggsave("cafe_vs_us_true_est_distance_v2.pdf",width=6.2,height = 4)
#ggsave("cafe_vs_us_true_est_distance_v2.pdf",width=4.8,height = 4)

getwd()










corr_stats
corr_stats$method <- factor(corr_stats$method, levels = c("cafe", "us"),
                            labels = c("D2Star", "kf2d"))


corr_stats$portion <- factor(corr_stats$portion, levels = c("backbone", "queries"),
                             labels = c("Backbone", "Queries"))

corr_stats %>% filter(corr == "pearson" & portion =="Queries") %>% group_by(method) %>%
  summarise(m = median(coef))
corr_stats




insert <-ggplot(aes(y=coef, x= method ,color= portion),
                data=corr_stats[corr_stats$corr == "pearson",])+
  #geom_bar(position = "dodge",stat = "identity")+
  geom_boxplot()+
  #ylab("")+
  #scale_color_brewer(name="",palette = "Dark2", guide="none")+
  #scale_color_brewer(name="",palette = "Dark2")+
  #scale_x_discrete(labels=c( expression(italic("D")[italic("2")]^italic("*")),"kf2d"))+
  scale_color_manual(name="", values = c(dark2_colors[1],dark2_colors[2]))+
  scale_x_discrete(labels=c(expression(italic("D")[italic("2")]^italic("*")), 'kf2d' )) +
  ylab("Correlation")+
  theme(
    #legend.position = c(0.45,-0.43), legend.margin=margin(0,0,0,0),legend.direction="horizontal",
     #   legend.text=element_text(size=7),
      #  legend.key.size = unit(0.08, 'cm'), 
    legend.position = "none",
        axis.title.x = element_blank(), axis.title.y = element_text(size=6),
        axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 6),
        plot.margin = unit(c(0,0.0,0,0), "cm"))
#guides(color = guide_legend(reverse = TRUE))
insert


summary(corr_stats[corr_stats$corr == "pearson" & corr_stats$method == "kf2d" & corr_stats$portion=="Backbone",])

plot.with.inset <-
  ggdraw() +
  draw_plot(big_plot) +
  draw_plot(insert, x = 0.82, y = .13, width = .16, height = .3)
plot.with.inset
ggsave(filename = "cafe_vs_us_inset.pdf", 
       plot = plot.with.inset,
       width = 6.2, 
       height = 4)

#ggsave(filename = "cafe_vs_us_inset.pdf", 
#       plot = plot.with.inset,
#       width = 4.8, 
#       height = 4)



#big_plot + inset_element(
#  insert, 
#  left = 0.8, 
#  bottom = 0.02, 
#  right = unit(1, 'npc'), 
#  top = unit(0.3, 'npc')
#)


#ggsave("cafe_vs_us_inset.pdf",width=6.2,height = 4)

#big_plot + coord_cartesian() +
#  annotation_custom(ggplotGrob(insert), xmin = 1, xmax = 3, ymin = 1, ymax = 2)

#er = median(total$V3.x)
#er
#er = median(total$V3.y)
#er


unique(dist_df$method)
unique(dist_df$portion)
a = dist_df[dist_df$method=="kf2d" & dist_df$portion=="Backbone",]$dist
b = dist_df[dist_df$method=="kf2d" & dist_df$portion=="Backbone",]$true

a2 = dist_df[dist_df$method=="kf2d" & dist_df$portion=="Queries",]$dist
b2 = dist_df[dist_df$method=="kf2d" & dist_df$portion=="Queries",]$true

a3 = dist_df[dist_df$method=="D2Star" & dist_df$portion=="Backbone",]$dist
b3 = dist_df[dist_df$method=="D2Star" & dist_df$portion=="Backbone",]$true

a4 = dist_df[dist_df$method=="D2Star" & dist_df$portion=="Queries",]$dist
b4 = dist_df[dist_df$method=="D2Star" & dist_df$portion=="Queries",]$true

a
res <- cor.test(a3, b3, method = 'pearson')
res





#For ppt

ggplot(aes(x=true/100, y=dist/100, color = portion),
                 data=dist_df[(dist_df$method == "kf2d" | dist_df$method == "D2Star") & (dist_df$portion == "Backbone"),])+
  #geom_bar(position = "dodge",stat = "identity")+
  geom_abline(linetype=2,size=0.9)+
  stat_smooth(method="lm",se=FALSE,size=1.4)+
  geom_point(alpha = 0.3, size = 0.5)+
  facet_wrap(~ method)+
  theme_classic()+
  scale_color_brewer(palette = "Dark2",name="")+
  #geom_boxplot(alpha=0.6)+
  #geom_point()+
  #geom_errorbar(stat="summary")+
  #geom_line() +
  #stat_summary(fun.y = "median", geom = "point", size = 3)+
  #geom_boxplot(outlier.shape = NA) +
  #geom_abline(color="red")+
  #scale_y_continuous(limits = c(0.0000, 0.001))
  xlab("True distance")+
  ylab("Estimated distance")+
  coord_cartesian(xlim=c(0,4.5),ylim=c(0,4.5))+
  theme(legend.position = c(.1,.9), legend.margin=margin(t = -0.5, unit='cm'), )

ggsave("cafe_vs_us_true_est_distance_ppt_backbone.pdf",width=6.2,height = 4)



############################################################################
############################################################################
############################################################################

