## 0. Loading packages ##
#############################################################################################
if(!require('data.table')){install.packages('data.table')}
if(!require('ggpubr')){install.packages('ggpubr')}
if(!require('ggplot2')){install.packages('ggplot2')}
if(!require('DTK')){install.packages('DTK')}
if(!require('svglite')){install.packages('svglite')}
library(data.table)
library(ggpubr)
library(ggplot2)
library(DTK)
library(svglite)

###########################
## 1. Preparing the data ##
#############################################################################################
directory <- getwd()
df <- fread(paste(directory, 'Supplementary/runEMIP/w%entries_w%diseases(accumulated)_w%uniprot_w%orpha(processed).csv', sep = '/'),
            select = c(1,2,5,6,7,8,10,11,12))
colnames(df) <- c("stat", "ent", "len", "cair", "int", "mip", "occ", "dis", "age")
df <- df[df$stat=='reviewed',]
df$dis[df$dis != ""] <- 'Disease'
df$dis[df$dis == ""] <- 'Non Disease'
df$age <- factor(df$age, labels = c("Mammalia", "Vertebrata", "Eumetazoa", "Opisthokonta", "Eukaryota", "Euk_Archaea", "Euk+Bac", "Cellular_organisms"),
                        levels = c("Mammalia", "Vertebrata", "Eumetazoa", "Opisthokonta", "Eukaryota", "Euk_Archaea", "Euk+Bac", "Cellular_organisms"), ordered = T)
occ_sbt <- function (x) {
  if (is.na(x)) {
    x <- NA
  }
  else if (x == 0) {
    x <- '3 - No Disease'
  }
  else if (x <= 0.000001) {
    x <- '2 - Extremely rare'
  }
  else if (x > 0.000001) {
    x <- '1 - Rare'
  }
}
df$occ <- sapply(df$occ, occ_sbt)
df$stat <- NULL
df$ent <- NULL
df <- df[!is.na(df$occ),]
df <- df[order(occ),]
df <- df[!is.na(df$age)]

#######################################
## 2. Drawing Errorplots w/ Gene Age ##
#############################################################################################
err_size <- 0.35
ball_size <- 0.1
fs <- 6
fsl <- 6
font_family <- "sans"
color_palette <- "npg"
mip <- ggerrorplot(data = df, x = 'occ', y = 'mip', desc_stat = 'mean_ci',color = "age",
                   error.plot = "errorbar", xlab = "Occurence category", size =  err_size,
                   ylab = "EMIP", palette = color_palette, add = "mean",
                   position = position_dodge(0.5), add.params = list(size = ball_size))+
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  font("legend.text", size = fs, family = font_family)
int <- ggerrorplot(data = df, x = 'occ', y = 'int', desc_stat = 'mean_ci', color = "age",
                   error.plot = "errorbar", xlab = "Occurence category", size = err_size,
                   ylab = "Number of interactions", palette = color_palette, add = "mean",
                   position = position_dodge(0.5), add.params = list(size = ball_size))+
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  font("legend.text", size = fs, family = font_family)
len <- ggerrorplot(data = df, x = 'occ', y = 'len', desc_stat = 'mean_ci', color = "age",
                   error.plot = "errorbar", xlab = "Occurence category", size = err_size,
                   ylab = "Protein length", palette = color_palette, add = "mean",
                   position = position_dodge(0.5), add.params = list(size = ball_size))+
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  font("legend.text", size = fs, family = font_family)
cair <- ggerrorplot(data = df, x = 'occ', y = 'cair', desc_stat = 'mean_ci', color = "age",
                    error.plot = "errorbar", xlab = "Occurence category", size = err_size,
                    ylab = "CAIR", palette = color_palette, add = "mean",
                    position = position_dodge(0.5), add.params = list(size = ball_size)) +
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  font("legend.text", size = fs, family = font_family)
fig <- ggarrange(mip, int, len, cair, ncol = 2, nrow = 2, vjust = 0,
                 common.legend = T ,legend = "bottom", align = 'v',
                 font.label = list(size = fsl, face = "bold", family = font_family))

p_mip <- ggplot_build(mip)
p_int <- ggplot_build(int)
p_len <- ggplot_build(len)
p_cair<- ggplot_build(cair)

#############################################
## 3. Performing Dunnett-Tukey-Kramer test ##
#############################################################################################
f <- gl.unequal(n=3, k = c(nrow(df[df$occ == "1 - Rare"]),
                           nrow(df[df$occ == "2 - Extremely rare"]),
                           nrow(df[df$occ == "3 - No Disease"])))
p_iter <- list(0.05, 0.001, 0.0001)
for (p in p_iter){
  mip_dtk <- as.data.frame(DTK.test(x = df$mip, f = f, a = p)[[2]])
  if (mip_dtk[1,2]*mip_dtk[1,3]>0){
    if (p==0.05){
      mip_sig_2_1 <- "*"
    }
    else if (p==0.001){
      mip_sig_2_1 <- "**"
    }
    else if (p==0.0001){
      mip_sig_2_1 <- "***"
    }
  else {
    mip_sig_2_1 <- "ns"
  }
  }
  if (mip_dtk[2,2]*mip_dtk[2,3]>0){
    if (p==0.05){
      mip_sig_3_1 <- "*"
    }
    else if (p==0.001){
      mip_sig_3_1 <- "**"
    }
    else if (p==0.0001){
      mip_sig_3_1 <- "***"
    }
    else {
      mip_sig_3_1 <- "ns"
    }
  }
  if (mip_dtk[3,2]*mip_dtk[3,3]>0){
    if (p==0.05){
      mip_sig_3_2 <- "*"
    }
    else if (p==0.001){
      mip_sig_3_2 <- "**"
    }
    else if (p==0.0001){
      mip_sig_3_2 <- "***"
    }
    else {
      mip_sig_3_2 <- "ns"
    }
  }
}
mip_sig <- data.frame("group1" = c("2 - Extremely rare", "3 - No Disease", "3 - No Disease"),
                      "group2" = c("1 - Rare", "1 - Rare", "2 - Extremely rare"),
                      "p.signif" = c(mip_sig_2_1, mip_sig_3_1, mip_sig_3_2),
                      "y_pos" = c(max(p_mip[["data"]][[2]][["ymax"]])*0.985,
                                  max(p_mip[["data"]][[2]][["ymax"]])*1.02,
                                  max(p_mip[["data"]][[2]][["ymax"]])*0.97))

for (p in p_iter){
  int_dtk <- as.data.frame(DTK.test(x = df$int, f = f, a = p)[[2]])
  if (int_dtk[1,2]*int_dtk[1,3]>0){
    if (p==0.05){
      int_sig_2_1 <- "*"
    }
    else if (p==0.001){
      int_sig_2_1 <- "**"
    }
    else if (p==0.0001){
      int_sig_2_1 <- "***"
    }
    else {
      int_sig_2_1 <- "ns"
    }
  }
  if (int_dtk[2,2]*int_dtk[2,3]>0){
    if (p==0.05){
      int_sig_3_1 <- "*"
    }
    else if (p==0.001){
      int_sig_3_1 <- "**"
    }
    else if (p==0.0001){
      int_sig_3_1 <- "***"
    }
    else {
      int_sig_3_1 <- "ns"
    }
  }
  if (int_dtk[3,2]*int_dtk[3,3]>0){
    if (p==0.05){
      int_sig_3_2 <- "*"
    }
    else if (p==0.001){
      int_sig_3_2 <- "**"
    }
    else if (p==0.0001){
      int_sig_3_2 <- "***"
    }
    else {
      int_sig_3_2 <- "ns"
    }
  }
}
int_sig <- data.frame("group1" = c("2 - Extremely rare", "3 - No Disease", "3 - No Disease"),
                      "group2" = c("1 - Rare", "1 - Rare", "2 - Extremely rare"),
                      "p.signif" = c(int_sig_2_1, int_sig_3_1, int_sig_3_2),
                      "y_pos" = c(max(p_int[["data"]][[2]][["ymax"]])*0.99,
                                  max(p_int[["data"]][[2]][["ymax"]])*1.03,
                                  max(p_int[["data"]][[2]][["ymax"]])*0.97))

for (p in p_iter){
  len_dtk <- as.data.frame(DTK.test(x = df$len, f = f, a = p)[[2]])
  if (len_dtk[1,2]*len_dtk[1,3]>0){
    if (p==0.05){
      len_sig_2_1 <- "*"
    }
    else if (p==0.001){
      len_sig_2_1 <- "**"
    }
    else if (p==0.0001){
      len_sig_2_1 <- "***"
    }
    else {
      len_sig_2_1 <- "ns"
    }
  }
  if (len_dtk[2,2]*len_dtk[2,3]>0){
    if (p==0.05){
      len_sig_3_1 <- "*"
    }
    else if (p==0.001){
      len_sig_3_1 <- "**"
    }
    else if (p==0.0001){
      len_sig_3_1 <- "***"
    }
    else {
      len_sig_3_1 <- "ns"
    }
  }
  if (len_dtk[3,2]*len_dtk[3,3]>0){
    if (p==0.05){
      len_sig_3_2 <- "*"
    }
    else if (p==0.001){
      len_sig_3_2 <- "**"
    }
    else if (p==0.0001){
      len_sig_3_2 <- "***"
    }
    else {
      len_sig_3_2 <- "ns"
    }
  }
}
len_sig <- data.frame("group1" = c("2 - Extremely rare", "3 - No Disease", "3 - No Disease"),
                      "group2" = c("1 - Rare", "1 - Rare", "2 - Extremely rare"),
                      "p.signif" = c(len_sig_2_1, len_sig_3_1, len_sig_3_2),
                      "y_pos" = c(max(p_len[["data"]][[2]][["ymax"]])*0.99,
                                  max(p_len[["data"]][[2]][["ymax"]])*1.02,
                                  max(p_len[["data"]][[2]][["ymax"]])*0.97))

for (p in p_iter){
  cair_dtk <- as.data.frame(DTK.test(x = df$cair, f = f, a = p)[[2]])
  if (cair_dtk[1,2]*cair_dtk[1,3]>0){
    if (p==0.05){
      cair_sig_2_1 <- "*"
    }
    else if (p==0.001){
      cair_sig_2_1 <- "**"
    }
    else if (p==0.0001){
      cair_sig_2_1 <- "***"
    }
    else {
      cair_sig_2_1 <- "ns"
    }
  }
  if (cair_dtk[2,2]*cair_dtk[2,3]>0){
    if (p==0.05){
      cair_sig_3_1 <- "*"
    }
    else if (p==0.001){
      cair_sig_3_1 <- "**"
    }
    else if (p==0.0001){
      cair_sig_3_1 <- "***"
    }
    else {
      cair_sig_3_1 <- "ns"
    }
  }
  if (cair_dtk[3,2]*cair_dtk[3,3]>0){
    if (p==0.05){
      cair_sig_3_2 <- "*"
    }
    else if (p==0.001){
      cair_sig_3_2 <- "**"
    }
    else if (p==0.0001){
      cair_sig_3_2 <- "***"
    }
    else {
      cair_sig_3_2 <- "ns"
    }
  }
}
cair_sig <- data.frame("group1" = c("2 - Extremely rare", "3 - No Disease", "3 - No Disease"),
                      "group2" = c("1 - Rare", "1 - Rare", "2 - Extremely rare"),
                      "p.signif" = c(cair_sig_2_1, cair_sig_3_1, cair_sig_3_2),
                      "y_pos" = c(max(p_cair[["data"]][[2]][["ymax"]])*0.998,
                                  max(p_cair[["data"]][[2]][["ymax"]])*1.002,
                                  max(p_cair[["data"]][[2]][["ymax"]])*0.996))

################################
## 4. Drawing overall gglines ##
#############################################################################################
l_mip <- ggline(data = df, x = 'occ', y = 'mip', add = 'mean_ci', error.plot = "errorbar",
              ylim = p_mip[["layout"]][["panel_scales_y"]][[1]][["range"]][["range"]],
              ylab = "EMIP", xlab = "Occurence category", plot_type = "l", size = 0.3)+
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  stat_pvalue_manual(data = mip_sig, label = "p.signif",
                     y.position = "y_pos", tip.length = 0.001)

l_int <- ggline(data = df, x = 'occ', y = 'int', add = 'mean_ci', error.plot = "errorbar",
              ylim = p_int[["layout"]][["panel_scales_y"]][[1]][["range"]][["range"]],
              ylab = "Number of interactions", xlab = "Occurence category",
              plot_type = "l", size = 0.3)+
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  stat_pvalue_manual(data = int_sig, label = "p.signif",
                     y.position = "y_pos", tip.length = 0.001)

l_len <- ggline(data = df, x = 'occ', y = 'len', add = 'mean_ci', error.plot = "errorbar",
              ylim = p_len[["layout"]][["panel_scales_y"]][[1]][["range"]][["range"]],
              ylab = "Protein length", xlab = "Occurence category",
              plot_type = "l", size = 0.3)+
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  stat_pvalue_manual(data = len_sig, label = "p.signif",
                     y.position = "y_pos", tip.length = 0.001)

l_cair <- ggline(data = df, x = 'occ', y = 'cair', add = 'mean_ci', error.plot = "errorbar",
               ylim = p_cair[["layout"]][["panel_scales_y"]][[1]][["range"]][["range"]],
               ylab = "CAIR", xlab = "Occurence category", plot_type = "l", size = 0.3)+
  font("xy", size = fs, family = font_family, face = 'bold')+
  font("xy.text", size = fs, family = font_family)+
  stat_pvalue_manual(data = cair_sig, label = "p.signif",
                     y.position = "y_pos", tip.length = 0.004)

overlay_fig <- ggarrange(l_mip, l_int, l_len, l_cair, labels = c("a", "b", "c", "d"),
                         font.label = list(size = fsl, face = "bold", family = font_family),
                         ncol = 2, nrow = 2)

svglite(file = "fig4_legend.svg", width = 7.08, height = 4)
fig
dev.off()

svglite(file = "fig4_ov.svg", width = 7.08, height = 4, bg = "transparent")
overlay_fig
dev.off()
