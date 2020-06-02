## 0. Loading packages ##
#############################################################################################
if(!require('data.table')){install.packages('data.table')}
if(!require('svglite')){install.packages('svglite')}
if(!require('igraph')){install.packages('igraph')}
if(!require('tidyverse')){install.packages('tidyverse')}
if(!require('Cairo')){install.packages('Cairo')}
library(data.table)
library(svglite)
library(igraph)
library(tidyverse)
library(Cairo)


###########################
## 1. Preparing the data ##
#############################################################################################
directory <- getwd()
df <- fread(paste(directory, "Supplementary/runEMIP/w%Entries_w%diseases_w%uniprot_wo%orpha.csv", sep = "/"),
            select = c(2, 8))
colnames(df) <- c("ent", "mip")
sizes <- df[, c("ent", "mip")]

# loading the uniprot PPI networks file
links <- fread(paste(directory, "Supplementary/runEMIP/uniprotnetwork_wo_animals.csv", sep = "/"))
colnames(links) <- c("intA", "intB")
links <- links[links$intA != links$intB, ]
links <- unique(links)

############################################
## 2. Customizing the igraph vertex shape ##
#############################################################################################
mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.frame.lty  <- params("vertex", "lty")
  if (length(vertex.frame.lty) != 1 && !is.null(v)) {
    vertex.frame.lty <- vertex.frame.lty[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width, vertex.frame.lty,
         FUN=function(x, y, bg, fg, size, lwd, lty) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd, lty=lty,
                   circles=size, add=TRUE, inches=FALSE)
         })
}
add.vertex.shape("fcircle", clip=igraph.shape.noclip,
                 plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.lty=1,
                                                vertex.frame.width=1))

############################################
## 3. Preparing and drwaing the networks  ##
#############################################################################################
protnet <- function(Label, Entry, Protein_name, Output_file, Title_pos) {
prim_edge <- c()
main_nodes <- c(Entry)
for (i in main_nodes){
  prim_edge = rbind(links[links$intA==i | links$intB==i, ], prim_edge)
}
ents <- c(prim_edge$intA, prim_edge$intB)
ents <- unique(ents)
scnd_edge = c()
for (i in ents){
  scnd_edge = rbind(links[links$intA==i | links$intB==i, ], scnd_edge)
}
scnd_edge <- scnd_edge[, c(1,2)]
sources <- scnd_edge %>%
  distinct(intA) %>%
  rename(ent = intA)
destinations <- scnd_edge %>%
  distinct(intB) %>%
  rename(ent = intB)
nodes <- full_join(sources, destinations, by = "ent")
nodes <- left_join(nodes, sizes, by = "ent")
nodes$mip <- replace_na(nodes$mip, 0)                    
net <- graph.data.frame(scnd_edge, nodes, directed=F)
vrtx <- V(net)
svglite(file = Output_file, width = 1.77, height = 1.77)
par(mar=c(1,1,1,1))
plot(x = net,
     layout=layout.kamada.kawai,
     vertex.shape="fcircle",
     vertex.frame.color="#000000",
     vertex.lty=1,
     vertex.size=0.04*(vrtx$mip)**0.6,
     vertex.frame.width=0.25,
     vertex.label=NA,
     edge.curved=0,
     vertex.color=c("#00A087"),
     edge.color="#464646",
     edge.width=0.12)
title(Protein_name, line = -(7+Title_pos/5), font.main = 1, family="sans", cex.main=0.5)
title(Label, adj=0, line = -0.36, font.main = 2, family="sans", cex.main=0.5)
dev.off()
}

#############################
## 4. Getting output files ##
#############################################################################################
protnet('a)', 'Q8WZ42', 'Connectin', 'network_01.svg', 1)
protnet('b)', 'P05067', 'Amyloid-beta\nprecursor protein', 'network_02.svg', 2)
protnet('c)', 'P0CG48', 'Polyubiquitin-C', 'network_03.svg', 1)
#protnet('d)', 'Q8WXI7', 'Ovarian carcinoma\nantigen CA125', 'network_04.svg', 2)
protnet('e)', 'Q9NRI5', 'Disrupted in schizo-\nphrenia 1 protein', 'network_05.svg', 2)
protnet('f)', 'P04637', 'Cellular tumor\nantigen p53', 'network_06.svg', 2)
protnet('g)', 'Q09472', 'Histone acetyl-\ntransferase p300', 'network_07.svg', 2)
protnet('h)', 'P00533', 'Epidermal growth\nfactor receptor', 'network_08.svg', 2)
protnet('i)', 'P62993', 'Growth factor receptor-\nbound protein 2', 'network_09.svg', 2)
protnet('j)', 'Q8NF91', 'Nesprin-1', 'network_10.svg', 1)
protnet('k)', 'P63104', 'Protein kinase C\ninhibitor protein 1', 'network_11.svg', 2)
protnet('l)', 'P78362', 'SRSF protein kinase 2', 'network_12.svg', 2)
protnet('m)', 'Q03001', 'Dystonin', 'network_13.svg', 1)
protnet('n)', 'Q5VST9', 'Obscurin', 'network_14.svg', 1)
protnet('o)', 'P38398', 'Breast cancer type 1\nsusceptibility protein', 'network_15.svg', 2)
protnet('p)', 'Q15149', 'Plectin', 'network_16.svg', 1)
