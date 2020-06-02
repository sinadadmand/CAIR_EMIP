## 0. Installing necessary packages ##
#############################################################################################
if(!require('readxl')){install.packages('readxl')}
if(!require('svglite')){install.packages('svglite')}

## 0.1. Loading packages
########################
library(readxl) # for "read_excel"
library(svglite) # for "svglite"

#############################
## 1. Drawing denisty plots##
#############################################################################################
directory <- getwd()
bw <- 0.0003 # bandwidth of density plots
family <- "sans" # font family of the plots
width <- 3.54 # width of the output plots in inches
height <- 1.82 # height of the output plots in inches
pointsize <- 6 # pointsize of the vector image outputs

## 1.1. Wave of life, i.e. all organisms
########################################
all_organisms <- read_excel(paste(directory, 'Supplementary/runcair/Figure_3_remove_dup.xlsx',
                                  sep = '/'), sheet = "All organisms")
from <- 0.875
to <- 0.945
wave_life_dens <- density(all_organisms$CAIR, from=from, to=to, bw=bw)
svglite(file = 'Fig3_a.svg', pointsize = pointsize, width = width*2, height = height)
plot(wave_life_dens, main = "", col="#E64B35",
     xlab='CAIR of all organisms', xaxt='n',
     ylab='', yaxt='n', family=family)
axis(side=1, at=seq(from = from, to = to, by = 0.005))
title(ylab="CAIR Density", line=0.5, cex.lab=1, family=family)
dev.off()

## 1.2. Proteobacteria
######################
proteobacteria <- read_excel(paste(directory, 'Supplementary/runcair/Figure_3_remove_dup.xlsx',
                                   sep = '/'), sheet = "Proteobacteria")
from <- 0.890
to <- 0.936
proteo_dens <- density(proteobacteria$CAIR, from=from, to=to, bw=bw)
svglite(file = 'Fig3_b.svg', pointsize = pointsize, width = width, height = height)
plot(proteo_dens, main = "", col="#3C5488",
     xlab='CAIR of the Proteobacteria phylum',
     ylab='', yaxt ='n', family=family)
title(ylab="CAIR Density", line=0.5, cex.lab=1, family=family)
dev.off()

## 1.3. Gammaproteobacteria
###########################
gammaproteobacteria <- read_excel(paste(directory, 'Supplementary/runcair/Figure_3_remove_dup.xlsx',
                                        sep = '/'), sheet = "Gammaproteobacteria")
from <- 0.900
to <- 0.936
gammaproteo_dens <- density(gammaproteobacteria$CAIR, from=from, to=to, bw=bw)
svglite(file = 'Fig3_c.svg', pointsize = pointsize, width = width, height = height)
plot(gammaproteo_dens, main = "", col="#00A087",
     xlab='CAIR of the Gammaproteobacteria class',
     ylab='', yaxt ='n', family=family)
title(ylab="CAIR Density", line=0.5, cex.lab=1, family=family)
dev.off()

## 1.4. Enterobacterales
########################
enterobacterales <- read_excel(paste(directory, 'Supplementary/runcair/Figure_3_remove_dup.xlsx',
                                     sep = '/'), sheet = "Enterobacterales")
from <- 0.910
to <- 0.936
entero_dens <- density(enterobacterales$CAIR, from=from, to=to, bw=bw)
svglite(file = 'Fig3_d.svg', pointsize = pointsize, width = width, height = height)
plot(entero_dens, main = "", col="#AE1F63",
     xlab='CAIR of the Enterobacterales order',
     ylab='', yaxt ='n', family=family)
title(ylab="CAIR Density", line=0.5, cex.lab=1, family=family)
dev.off()

## 1.5. Enterobacteriaceae
##########################
enterobacteriaceae <- read_excel(paste(directory, 'Supplementary/runcair/Figure_3_remove_dup.xlsx',
                                       sep = '/'), sheet = "Enterobacteriaceae")
from <- 0.926
to <- 0.935
enterobacter_dens <- density(enterobacteriaceae$CAIR, from=from, to=to, bw=bw)
svglite(file = 'Fig3_e.svg', pointsize = pointsize, width = width, height = height)
plot(enterobacter_dens, main = "", col="#CC9900",
     xlab='CAIR of the Enterobacteriaceae family',
     ylab='', yaxt ='n', family=family)
title(ylab="CAIR Density", line=0.5, cex.lab=1, family=family)
dev.off()

## 1.6. Escherichia
###################
escherichia <- read_excel(paste(directory, 'Supplementary/runcair/Figure_3_remove_dup.xlsx',
                                sep = '/'), sheet = "Escherichia")
from <- 0.930
to <- 0.935
escherichia_dens <- density(escherichia$CAIR, from=from, to=to, bw=bw)
svglite(file = 'Fig3_f.svg', pointsize = pointsize, width = width, height = height)
plot(escherichia_dens, main = "", col="#7E6148",
     xlab='CAIR of the Escherichia genus',
     ylab='', yaxt ='n', family=family)
title(ylab="CAIR Density", line=0.5, cex.lab=1, family=family)
dev.off()

## 1.7. E. coli
###############
e_coli <- read_excel(paste(directory, 'Supplementary/runcair/Figure_3_remove_dup.xlsx',
                                sep = '/'), sheet = "E. coli")
from <- 0.931
to <- 0.935
e_coli_dens <- density(e_coli$CAIR, from=from, to=to, bw=bw)
svglite(file = 'Fig3_g.svg', pointsize = pointsize, width = width, height = height)
plot(e_coli_dens, main = "", col="#1A237E",
     xlab='CAIR of the E. coli organism',
     ylab='', yaxt ='n', family=family)
title(ylab="CAIR Density", line=0.5, cex.lab=1, family=family)
dev.off()
