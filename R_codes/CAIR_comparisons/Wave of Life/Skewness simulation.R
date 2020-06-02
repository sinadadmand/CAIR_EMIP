## 0. Installing necessary packages ##
#############################################################################################
if(!require('data.table')){install.packages('data.table')}
if(!require('dplyr')){install.packages('dplyr')}
if(!require('plyr')){install.packages('plyr')}
if(!require('fGarch')){install.packages('fGarch')}
if(!require('moments')){install.packages('moments')}
if(!require('svglite')){install.packages('svglite')}

## 0.1. Loading packages
########################
library(data.table) # for "fread" and "rbind"
library(dplyr) # for "sample_n"
library(plyr) # for "count"
library(fGarch) # for "rsnorm"
library(moments) # for "skewness"
library(svglite) # for "svglite"

## 0.2. Staring the time
########################
start.time <- Sys.time()


###########################
## 1. Preparing the data ##
#############################################################################################
directory <- getwd()
df <- fread(paste(directory, 'Supplementary/runcair/Complete proteome CAIRs.csv', sep = '/'))
colnames(df) <- c("org_id", "cair", "org", "phyl", "code")
num_org <- plyr::count(df, vars = "code")
num_org <- num_org[-c(1), ] # removing the 0 code which are the unclassified excluded phyla
num_org <- num_org[(num_org$freq > 10), ]
max_samp <- quantile(num_org$freq, c(0.75))
min_samp <- quantile(num_org$freq, c(0.25))
num_org <- num_org[(num_org$freq >= min_samp & num_org$freq <= max_samp), ]
row.names(num_org) <- NULL

################################
## 2. Skewness of each phylum ##
#############################################################################################
skw <- num_org
for (i in num_org$code)
{
  phylum <- df$cair[df$code==i]
  skw$skw[which(num_org$code ==i)] <- skewness(phylum)
}
#View(skw)
print(paste("The mean skewness of interquantile range in tree of life is",
            round(mean(skw$skw), digits = 2)), quote=F)

#########################
## 3. Overall skewness ##
#############################################################################################
iterations <- 1000 # This may take upto 6 hours for calculation
result <- data.frame("index" = c(1:iterations))
for (sk in seq(-4, 4, 0.01)){
  ks <- {}
  if (sk == 0) next
  for (s in 1:iterations){
    set.seed(s)
    rand_samp <- {}
    rand_samp_i <- {}
    for (i in num_org$code)
    {
      rand_samp_i <- sample_n(data.frame(df$cair[df$code==i]), min_samp)
      rand_samp <- rbind(rand_samp_i, rand_samp)
      num_org$mean[which(num_org$code ==i)] <- mean(as.matrix(rand_samp_i))
      num_org$median[which(num_org$code ==i)] <- median(as.matrix(rand_samp_i))
      num_org$sd[which(num_org$code ==i)] <- sd(as.matrix(rand_samp_i))
    }
    simulation <- {}
    for (i in 1:nrow(num_org))
    {
      simulation <- c(rsnorm(n=min_samp, mean = num_org$median[i],
                             sd = num_org$sd[i], xi = sk), simulation)
    }
    ks_prior <- ks
    ks <- ks.test(x = simulation, y = as.matrix(rand_samp))$p.value
    ks <- c(ks, ks_prior)
  }
  result[,paste0(sk)] <- ks
}
mean_col <- as.data.frame(sapply(result[, -c(1)], FUN=mean))
sk <- as.numeric(rownames(mean_col)[which(mean_col %% 1==max(mean_col))])[1]
print(paste("The overall skewness of the tree of life according tp the simulation is",
            sk), quote=F)
# sk (the overall skewness) would come up to be -0.90 after 1000 iterations.

##############################
## 4. Drawing density plots ##
#############################################################################################
set.seed(123)
rand_samp <- {}
rand_samp_i <- {}
for (i in num_org$code)
{
  rand_samp_i <- sample_n(data.frame(df$cair[df$code==i]), min_samp)
  rand_samp <- rbind(rand_samp_i, rand_samp)
  num_org$mean[which(num_org$code ==i)] <- mean(as.matrix(rand_samp_i))
  num_org$median[which(num_org$code ==i)] <- median(as.matrix(rand_samp_i))
  num_org$sd[which(num_org$code ==i)] <- sd(as.matrix(rand_samp_i))
}

simulation <- {}
for (i in 1:nrow(num_org))
{
  simulation <- c(rsnorm(n=min_samp, mean = num_org$median[i],
                         sd = num_org$sd[i], xi = sk), simulation)
}

bw <- 0.0005
from <- 0.88
to <- 0.94
wol_real <- density(as.matrix(rand_samp), from=from, to=to, bw=bw)
svglite(file = 'New_Fig2_b.svg', pointsize = 6, width = 3.46, height = 1.73)
plot(wol_real, main = "", col="#E64B35", xlab='', ylab='', yaxt ='n', family="sans")
title(ylab="CAIR Density", line=0.5, cex.lab=1, family="sans") + rug (num_org$median)
dev.off()

wol_simul <- density(simulation, from=from, to=to, bw=bw)
svglite(file = 'New_Fig2_a.svg', pointsize = 6, width = 3.46, height = 1.73)
plot(wol_simul, main = "", col="#4DBBD5", xlab='', ylab='', yaxt ='n', family="sans")
title(ylab="CAIR Density", line=0.5, cex.lab=1, family="sans") + rug (num_org$median)
dev.off()

svglite(file = 'New_Fig2_a&b.svg', pointsize = 6, width = 3.46, height = 1.73)
bw <- 0.0001
wol_real <- density(as.matrix(rand_samp), from=from, to=to, bw=bw)
wol_simul <- density(simulation, from=from, to=to, bw=bw)
plot(wol_real, main = "", col="#E64B35", xlab='', ylab='', yaxt ='n', family="sans")
lines(wol_simul, main = "", col="#4DBBD5", xlab='', ylab='', yaxt ='n', family="sans")
title(ylab="CAIR Density", line=0.5, cex.lab=1, family="sans") + rug (num_org$median)
dev.off()

print(paste("The p-value of the Kolmogorov-Smirnov test comparing the wave of life",
            "with the simulation is", round(max(mean_col), digits=4)), quote=F)
# the mean of p-values for sk being -0.90 would come up to be 0.4026.

## 4.1. Reporting the time elapsed
##################################
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
