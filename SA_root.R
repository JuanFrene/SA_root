packages <- c('edgeR', 'ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq",
              "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4",
              "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire",
              'igraph','Hmisc','VennDiagram',"microbiomeMarker", 'ohchibi','microbiomeSeq')
sapply(packages, require, character.only = TRUE)              

setwd("C:/Users/juanp/Google Drive/labs/Nottingham/SA project/16S analysis/")


###ASV table
asv.tableSA <- readRDS('C:/Users/juanp/Google Drive/labs/Nottingham/SA project/16S analysis/seqtab_final.rds')
asv.tableSA2<- otu_table(asv.tableSA, taxa_are_rows=FALSE)

### Taxonomy Table
taxaSA <- readRDS("C:/Users/juanp/Google Drive/labs/Nottingham/SA project/16S analysis/tax_final.rds")

###Metadata
metaSA <- read.table("Metadata.txt", header = TRUE, row.names = 1)

#Now we can make the phyloseq object
ps1 <- phyloseq(asv.tableSA2, tax_table(taxaSA), sample_data(metaSA))

taxa_names(ps1) <- paste0("ASV", seq(ntaxa(ps1)))

###Only from SA project
ps.2 = subset_samples(ps1, ID !=  "NNN")
ps.pruned <- prune_samples(sample_sums(ps.2)>=500, ps.2)
ps.3 = subset_taxa(ps.pruned, Kingdom == "Bacteria" | Kingdom == "Archaea")
ps.3.pruned <- prune_taxa(taxa_sums(ps.3)>=50, ps.3)
ps.SA  = subset_samples(ps.3.pruned, ID !=  "c-")
ps.SA.perc <- transform_sample_counts(ps.SA, function(x) x / sum(x)) 


###Subset by type of saple
ps.Root = subset_samples(ps.SA, Compartment ==  "Root")
ps.Root <- prune_taxa(taxa_sums(ps.Root)>=1, ps.Root)
ps.Root.perc <- transform_sample_counts(ps.Root, function(x) x / sum(x)) 
ps.Root.H <- transform(ps.Root, "hellinger")


ps.Soil = subset_samples(ps.SA, Compartment ==  "Soil")
ps.Soil <- prune_taxa(taxa_sums(ps.Soil)>=1, ps.Soil)
ps.Soil.perc <- transform_sample_counts(ps.Soil, function(x) x / sum(x)) 
ps.Soil.H <- transform(ps.Soil, "hellinger")


####PCoA
PCoA <- ordinate(ps.SA, "PCoA", "bray")
title ="PCoA - PCoA1 vs PCoA2"
ev1 <- (PCoA$values$Relative_eig[1])*100
ev2 <- (PCoA$values$Relative_eig[2])*100
xlab_text <- paste("PCoA1 (", round(ev1,2), "%)")
ylab_text <- paste("PCoA2 (", round(ev2,2), "%)")

x<-data.frame(PCoA$vectors[, 1:2])
x <- merge(x, ps.SA@sam_data, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL

PCoAmean = aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2)~Mutant,x,mean) #
PCoAmeands = aggregate(cbind(ds.x=Axis.1,ds.y=Axis.2)~Mutant,x,se)
PCoAmean2 =merge(PCoAmean,PCoAmeands,by="Mutant")


paleta_alive <- c("#C200FF",'#8B7500','#FF0000','#8B8B7A',"#FFB919",'#FF7F50',"#00CC7A",'#8B1A1A','#104E8B','#8B1C62',
                  '#1f78b4','#33a02c','#e31a1c','#ff7f00','#6a3d9a','#b15928','#a6cee3','#b2df8a','#fb9a99','#fdbf6f',
                  '#cab2d6','#ffff99','#8dd3c7','#bebada','#fb8072','#80b1d3','black','#fdb462','#d9d9d9','#bc80bd',
                  '#ccebc5','#ffed6f')

PCoAmean2
write.table(x, "PCoA Soil vs Root table.txt")

ggplot(data = x, aes(Axis.1, Axis.2)) +
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  #geom_errorbarh(mapping = aes(xmin =mean.x - ds.x, xmax = mean.x + ds.x),size = 0.02,alpha = 0.2) +
  #geom_errorbar(mapping = aes(ymin =mean.y - ds.y, ymax = mean.y + ds.y),size = 0.02,alpha = 0.2)  + 
  geom_point(size = 4, aes(color = Compartment),stroke = 1) + #
  #theme_ohchibi(size_panel_border = 2) + 
  xlab(label = xlab_text) + ylab(label = ylab_text) +
  scale_fill_manual(values = paleta_alive) +
  scale_color_manual(values = paleta_alive)
#scale_colour_gradientn(colours=c("Red")) #"Blue",

ggplot(data = PCoAmean2, aes(mean.x, mean.y)) +
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_point(size = 4, aes(color = Mutant),stroke = 1) + #
  theme_ohchibi(size_panel_border = 2) + 
  xlab(label = xlab_text) + ylab(label = ylab_text) +
  scale_fill_manual(values = paleta_alive) +
  scale_color_manual(values = paleta_alive)

###Cluster
Table.cluster.Root = cbind(ps.Root@sam_data, ps.Root@otu_table)
Table.cluster.Soil = cbind(ps.Soil@sam_data, ps.Soil@otu_table)
diff.mean <- Table.cluster.Root %>% group_by(Mutant)%>%
  summarise_all(mean)
Log = data.frame(diff.mean)
mydata[1:5,1:5]
rownames(Log) = Log[,1] 
mydata <- scale(Log[,-(1:8)])

row.names(mydata)
fviz_nbclust(mydata, kmeans, method = "gap_stat")

set.seed(123) # for reproducibility
km.res <- kmeans(mydata, 3, nstart = 25)
fviz_cluster(km.res, data = mydata, palette = "jco",
             ggtheme = theme_few(), main="Cluster by bact Frond")
km.res$betweenss

res.hc <- hclust(dist(mydata),  method = "ward.D")
fviz_dend(res.hc, cex = 0.7, k = 2, palette = "jco",rect = TRUE, main="Cluster by Mutant",  horiz = FALSE) 

res.hc$merge

###Barplot
#  Plotting Relative Abundance Bar Charts####
# phylum-level
remotes::install_github("umerijaz/microbiomeSeq")

ps.perc.Soil <- transform_sample_counts(ps.Soil, function(x) x / sum(x)) 
ps.perc.Root <- transform_sample_counts(ps.Root, function(x) x / sum(x)) 

ps.phyla.perc.Soil <- taxa_level(ps.perc.Soil, "Phylum")
ps.phyla.perc.Root <- taxa_level(ps.perc.Root, "Phylum")

# identify the 10 most abundant phylum
phylum.10.Soil <- names(sort(taxa_sums(ps.phyla.perc.Soil), TRUE)[1:10])
phylum.10.Root <- names(sort(taxa_sums(ps.phyla.perc.Soil), TRUE)[1:10])

# now we can use the list of the top phylum and subset those from the phyloseq object
ps.phylum.10.Soil <- prune_taxa(phylum.10.Soil, ps.phyla.perc.Soil)
ps.phylum.10.Root <- prune_taxa(phylum.10.Root, ps.phyla.perc.Root)

melt.phylum.Soil <- psmelt(ps.phylum.10.Soil)
melt.phylum.Root <- psmelt(ps.phylum.10.Root)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#CC6677", "#DDCC77")

#cargar dplyr
library(dplyr)

#Calcular medias por grupos de todas las columnas
melt.phylum.mean.Root <- melt.phylum.Root%>%group_by(Mutant, OTU)%>%
  summarise_all(mean)
melt.phylum.mean.Root=data.frame(melt.phylum.mean.Root[,c(1,2,4)])
unique(melt.phylum.mean.Root[,1])
write.table(melt.phylum.mean.Root, "Taxa Root.txt")


melt.phylum.mean.Soil <- melt.phylum.Soil%>%group_by(Mutant, OTU)%>%
  summarise_all(mean)
melt.phylum.mean.Soil=data.frame(melt.phylum.mean.Soil[,c(1,2,4)])
unique(melt.phylum.mean.Soil[,2])
write.table(melt.phylum.mean.Soil, "Taxa Soil.txt")

###Graficos
melt.phylum.mean.Root$Mutant <- factor(melt.phylum.mean.Root$Mutant,c("Col_0","ANAC029_NAP","ANAC032","ANAC046","ANAC053_NTL4","ANAC074","ANAC102","ANAC2_ATAF1","AtbZIP29",          
                                                                      "AtEBP_ERF72_RAP2.3","AtERF_1","AtERF71_HRE2","AtERF73/HRE1","AtERF73_HRE1","AtNPR1","AtNPR2","AtNPR3",            
                                                                      "AtNPR4","AtRHD6","AtTGA1","AtTGA4","AtZFP1","DAYSLEEPER","gammaMYB2","KUA1","RSL1","WRKY70"))

melt.phylum.mean.Soil$Mutant <- factor(melt.phylum.mean.Soil$Mutant,c('Soil',"Col_0","ANAC029_NAP","ANAC032","ANAC046","ANAC053_NTL4","ANAC074","ANAC102","ANAC2_ATAF1","AtbZIP29",          
                                                                      "AtEBP_ERF72_RAP2.3","AtERF_1","AtERF71_HRE2","AtERF73/HRE1","AtERF73_HRE1","AtNPR1","AtNPR2","AtNPR3",            
                                                                      "AtNPR4","AtRHD6","AtTGA1","AtTGA4","AtZFP1","DAYSLEEPER","gammaMYB2","KUA1","RSL1","WRKY70"))

bact_order = rev(c("Proteobacteria","Bacteroidetes","Actinobacteria", "Firmicutes","Chloroflexi","Cyanobacteria","Verrucomicrobia","Acidobacteria","Nitrospirae","Gemmatimonadetes"))
melt.phylum.mean.Root$OTU <- factor(melt.phylum.mean.Root$OTU,bact_order)
melt.phylum.mean.Soil$OTU <- factor(melt.phylum.mean.Soil$OTU,bact_order)


ggplot(melt.phylum.mean.Root, aes(x = Mutant , y = Abundance, fill = OTU)) + theme_few() +
  geom_bar(stat = "identity") + scale_fill_manual(values = safe_colorblind_palette) + theme(axis.text.x = element_text(angle = 45, hjust=1))

ggplot(melt.phylum.mean.Soil, aes(x = Mutant , y = Abundance, fill = OTU)) + theme_few() +
  geom_bar(stat = "identity") + scale_fill_manual(values = safe_colorblind_palette) + theme(axis.text.x = element_text(angle = 45, hjust=1))

###PERMANOVA
#Define Bray Curtis as the dissimilarity
ps.Soil2 = subset_samples(ps.Soil.perc, Mutant !=  "AtbZIP29")
sample_ps.Soil<- data.frame(sample_data(ps.Soil2))

ps.Root2 = subset_samples(ps.Root.perc, Mutant !=  "AtbZIP29")
sample_ps.Root<- data.frame(sample_data(ps.Root2))

cbn <- combn(x=unique(sample_ps.Soil$Mutant), m = 2)
p <- c()
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(ps.Soil2, Mutant %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = "bray") ~ Mutant, 
                               data = metadata_sub)
  c = data.frame(rbind(cbn[,i]))
  d = data.frame(permanova_pairwise$`Pr(>F)`[1])
  f = cbind(c,d)
  p <- rbind(p, f)
}
colnames(p) = c('col1','col2','Perc')
p_Col3 <- p %>% subset(col1  == 'Col_0') %>% droplevels

Root.anosim.bray
Soil.anosim.bray=cbind(p_Col,p_Col3[,3],p_Col2[,3])
colnames(Soil.anosim.bray) = c('col1','col2','Normal','Perc','H.Trans')

write.table(Soil.anosim.bray, "Soil ANOSIM bray.txt")


#Root and Front vs water
ps.SA
ps.SA_core = core(ps.SA , detection = 0.01, prevalence = 40/100)

tableSA = cbind(ps.SA_core@sam_data, ps.SA_core@otu_table)
tableSA2=tableSyncom[,-c(1:4,7)]
tableSA2$Mutant = as.factor(tableSA2$Mutant)
tableSA2$Compartment = as.factor(tableSA2$Compartment)

melted_tableSA <- tableSA2 %>% melt

ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

melted_tableSA

log2foldchange_SA <- c()
log2foldchange_SA_Complete2 =c()

cbn2 <- combn(x=unique(ps.Soil2@sam_data$Mutant), m = 2)
p <- c()

cbn2
for(Mut in melted_tableSA$Mutant  %>% unique){
  melted_sub2 <- melted_tableSyncom %>% subset(Compartment ==  Comp) %>% droplevels
  for(ASV in melted_sub2$variable  %>% unique){
    melted_sub <- melted_sub2 %>% subset(variable ==  ASV) %>% droplevels #'ASV522'
    
    Front = data.frame(t(melted_sub[melted_sub$Compartment=='Front',][,7]))
    Water = data.frame(t(melted_sub[melted_sub$Compartment=='Water',][,7]))
    
    Front_mean = apply(Front, 1, mean) 
    Water_mean = apply(Water, 1, mean) 
    
    Frontmean_Watermean <- log2(Front_mean+1) - log2(Water_mean+1) 
    
    Frontmean_Watermean_statistic = t.test(Front,Water)
    pvalueFrontmean_Watermean = Frontmean_Watermean_statistic$p.value#Water_mean+1
    
    result = cbind(Frontmean_Watermean,pvalueFrontmean_Watermean)
    row.names(result)= ASV
    log2foldchange_Front <- rbind(log2foldchange_Front,result)
  }
  Compartment <- data.frame(rep('Front', 97))
  Species_ID= data.frame(rep(specie, 97))
  log2foldchange_Front_Complete = cbind(Compartment,Species_ID,row.names(log2foldchange_Front),log2foldchange_Front)
  colnames(log2foldchange_Front_Complete)=c('Compartment','Species_ID','ASV','diff','pvalue')
  
  log2foldchange_Front_Complete2 <- rbind(log2foldchange_Front_Complete2,log2foldchange_Front_Complete)
}

log2foldchange_Root <- c()
log2foldchange_Root_Complete2 =c()

for(specie in melted_tableSyncom$Specie.1  %>% unique){
  melted_sub2 <- melted_tableSyncom %>% subset(Specie.1 ==  specie) %>% droplevels
  for(ASV in melted_sub2$variable  %>% unique){
    melted_sub <- melted_sub2 %>% subset(variable ==  ASV) %>% droplevels #'ASV522'
    
    Root = data.frame(t(melted_sub[melted_sub$Compartment=='Root',][,7]))
    Water = data.frame(t(melted_sub[melted_sub$Compartment=='Water',][,7]))
    
    Root_mean = apply(Root, 1, mean) 
    Water_mean = apply(Water, 1, mean) 
    
    Rootmean_Watermean <- log2(Root_mean+1) - log2(Water_mean+1) 
    
    Rootmean_Watermean_statistic = t.test(Root,Water)
    pvalueRootmean_Watermean = Rootmean_Watermean_statistic$p.value#Water_mean+1
    
    result = cbind(Rootmean_Watermean,pvalueRootmean_Watermean)
    row.names(result)= ASV
    log2foldchange_Root <- rbind(log2foldchange_Root,result)
  }
  Compartment <- data.frame(rep('Root', 97))
  Species_ID= data.frame(rep(specie, 97))
  log2foldchange_Root_Complete = cbind(Compartment,Species_ID,row.names(log2foldchange_Root),log2foldchange_Root)
  colnames(log2foldchange_Root_Complete)=c('Compartment','Species_ID','ASV','diff','pvalue')
  
  log2foldchange_Root_Complete2 <- rbind(log2foldchange_Root_Complete2,log2foldchange_Root_Complete)
}

log2foldchange_Complete = rbind(log2foldchange_Root_Complete2,log2foldchange_Front_Complete2)

log2foldchange_Complete$Significance <- "NoSignificant"
pval_thres <- 0.05
log2foldchange_Complete$Significance[which(log2foldchange_Complete$pvalue < pval_thres)] <- "q < 0.05"
log2foldchange_Complete$Significance <- log2foldchange_Complete$Significance %>% factor

log2foldchange_Complete$Species_ID =factor(log2foldchange_Complete$Species_ID, c('PS','LY5280','SI9227','SP7498','SP9509',
                                                                                 'SI7820','SP5543','SP9192','LM7200','LT9243',
                                                                                 'LP7760','LV7650','LJ9250','LP0049','LM8389','LT9109',
                                                                                 'LA7339','LM7123','LP8539'))

log2foldchange_Completefdffd = rbind(log2foldchange_Complete[log2foldchange_Complete$Species_ID=='LJ9250',], log2foldchange_Complete[log2foldchange_Complete$Species_ID=='LM7123',])

ggplot(data = log2foldchange_Completefdffd, aes(Species_ID,ASV)) +
  geom_raster(aes(fill = diff)) +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.2,width = 1,height = 1) + #
  #geom_text(aes(label = Phylum),color = "black",size = 5)+
  facet_grid(~Compartment, space = "free",scales = "free") +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-0.5,0.5),na.value = "#D9D9D9",oob = squish,name = "Fold Change") +
  scale_color_manual(values = c('grey',"black"),na.value =  "transparent",name = "Significance vs Water") + #Significance Genotype vs Col-0
  theme_ohchibi(size_panel_border = 0.2) +
  theme_few()

#Aggregate the dataframe to display as heatmap

PCoAmean = aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2)~Species*Nutrient,tableSyncom2,mean)
Alpha.mean <- Alfadiv%>%group_by(Species, Compartment)%>%
  summarise_all(mean)


display <- tableSyncom2  %>%
  acast(formula = Species_ID*Compartment*~ASV,fill = 0,
        value.var = "diff") %>% scale

head(log2foldchange_Complete)
display2=data.frame(display)

library(factoextra)
res.display2 <- hclust(dist(t(tableSyncom2[,-(1:5)])),  method = "ward.D")#[c(1,7,13,19,25,31,37),]
fviz_dend(res.display2, cex = 0.35, k = 3, palette = "jco",rect = TRUE, main="Cluster by Samples",  horiz = TRUE) 

label_order_vswater = c('ASV4534','ASV2052','ASV2204','ASV1892','ASV181','ASV2080','ASV3818','ASV99','ASV3925',
                        'ASV3048','ASV3476','ASV572','ASV4398','ASV894','ASV1030','ASV2494','ASV2044','ASV4016',
                        'ASV427','ASV4297','ASV987','ASV662','ASV4235','ASV3323','ASV2755','ASV3613','ASV20',
                        'ASV885','ASV67','ASV1045','ASV582','ASV2307','ASV4061','ASV1579','ASV4273','ASV1348',
                        'ASV3152','ASV522','ASV1762','ASV175','ASV2536','ASV2653','ASV1523','ASV1525','ASV3027',
                        'ASV845','ASV79','ASV3719','ASV3274','ASV1744','ASV421','ASV1081','ASV3096','ASV2089',
                        'ASV1041','ASV3656','ASV3106','ASV3410','ASV385','ASV1322','ASV1637','ASV399','ASV3430',
                        'ASV2877','ASV1246','ASV1982','ASV2522','ASV3477','ASV3731','ASV3045','ASV1190','ASV2448',
                        'ASV3539','ASV1565','ASV4229','ASV1103','ASV4562','ASV3682','ASV18','ASV1321','ASV2209',
                        'ASV154','ASV3570','ASV594','ASV4437','ASV815','ASV3560','ASV1798','ASV2593','ASV3757',
                        'ASV217','ASV3587','ASV2845','ASV1909','ASV3315','ASV2118','ASV2826')

log2foldchange_Complete$ASV=factor(log2foldchange_Complete$ASV, label_order_vswater)
log2foldchange_Complete$Compartment=factor(log2foldchange_Complete$Compartment, c('Root','Front'))

ggplot(data = log2foldchange_Complete, aes(Species_ID,ASV)) +
  geom_raster(aes(fill = diff)) +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.2,width = 1,height = 1) + #
  #geom_text(aes(label = Phylum),color = "black",size = 5)+
  facet_grid(~Compartment, space = "free",scales = "free") +
  theme_few()+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-0.4,0.4),na.value = "#D9D9D9",oob = squish,name = "Fold Change") +
  scale_color_manual(values = c('grey',"black"),na.value =  "transparent",name = "Significance vs Water") + #Significance Genotype vs Col-0
  theme_ohchibi(size_panel_border = 0.2)+
  theme(axis.text.x = element_text(angle = -90,vjust=-0.05,hjust=-0.05, size = 7)) #, 


