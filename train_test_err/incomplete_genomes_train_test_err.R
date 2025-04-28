


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

##############################################################
# Build train test error plot


#err_epoch=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v61_8k_train_model_10KFULL_TOL_Claded_Unchunked_tr_test_err/train_test_err_claded_unchunked_clades_all.txt',sep=" ",header=T)
err_epoch=read.csv('train_test_err_claded_unchunked_clades_all.txt',sep=" ",header=T)


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
#plerror_epoch=read.csv('/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls/k7_v62_8k_train_model_10K_TOL_pl_error_per_epoch/summary_D1_TOL_pl_error_per_epoch.csv',sep=" ",header=T)
plerror_epoch=read.csv('summary_D1_TOL_pl_error_per_epoch.csv',sep=" ",header=T)
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


