library(vegan)
library(phyloseq)
library(ggpubr)

setwd("~/Software_dev/Benchmark_Nygaard/")

load("Nygaard_pipeline_analysis/Nygaard.rdata")

load("SituSeq_Nygaard_dataset/SituSeq_Nygaard.rdata")
SituSeq <- SituSeq.Nygaard

load("NanoASV_Nygaard_dataset/NanoASV_output/nanoasv.out/Results/Rdata/NanoASV.rdata")


sample_names(NanoASV) <- paste0(sample_names(NanoASV), "_NanoASV")
sample_names(Nygaard) <- paste0(sample_names(Nygaard), "_Nygaard")
sample_names(SituSeq) <- paste0(sample_names(SituSeq), "_SituSeq")


MOCK.phy <- merge_phyloseq(NanoASV, SituSeq, Nygaard)


MOCK.phy@sam_data$Pipeline <- sub("^barcode\\d+_", "", sample_names(MOCK.phy))



#Extraction de la table d'OTUs
OTU <- data.frame(MOCK.phy@otu_table)
metadata <- data.frame(MOCK.phy@sam_data)
#Construction des tables de Rochesses et d'abondance
rich <- data.frame(apply(OTU>0,2,sum))

Shan <- data.frame(diversity(t(OTU), index="shannon"))


#Synthese Richesse et DIversite
alpha_tab <- data.frame( "richness" = rich,
                         "shannon" = Shan)
colnames(alpha_tab) <- c("richness", "shannon")

alpha_tab <- cbind(metadata, alpha_tab)



rich_plot <- ggplot(data = alpha_tab, aes(x = ena_accession, y = richness)) +
  geom_col(aes(fill = ena_accession)) + theme_bw() +
  facet_wrap(~ Pipeline) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle=45, colour = "black", vjust=1, hjust = 1, size=5),
        legend.position = "none") + 
  ylab("ASV-level Numerical Richness") + xlab("Nygaard samples ID") + 
  labs(title = "Nygaard dataset bacterial populations Numerical Richness ",
       caption = date())

shannon_plot <- ggplot(data = alpha_tab, aes(x = ena_accession , y = shannon)) +
  geom_col(aes(fill = ena_accession)) + theme_bw() +
  facet_wrap(~ Pipeline) +          
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle=45, colour = "black", vjust=1, hjust = 1, size=5),
        legend.position = "right") + 
  ylab("ASV-level Shannon Index") + xlab("Nygaard samples ID") + 
  labs(title = "Nygaard dataset bacterial populations Shannon Index ",
       caption = date())


pdf("Figs/Pipelines_Alpha_div.pdf", he = 8, wi = 12)
ggarrange(rich_plot, shannon_plot)
dev.off()


alpha_tab_NanoASV <- alpha_tab[alpha_tab$Pipeline == "NanoASV",]
alpha_tab_Nygaard <- alpha_tab[alpha_tab$Pipeline == "Nygaard",]
alpha_tab_SituSeq <- alpha_tab[alpha_tab$Pipeline == "SituSeq",]

mod <- lm(alpha_tab_NanoASV$richness ~ alpha_tab_Nygaard$richness)
summary(mod)
plot(alpha_tab_NanoASV$richness ~ alpha_tab_Nygaard$richness)

mod <- lm(alpha_tab_NanoASV$shannon ~ alpha_tab_Nygaard$shannon)
summary(mod)
plot(alpha_tab_NanoASV$shannon ~ alpha_tab_Nygaard$shannon)

mod <- lm(alpha_tab_NanoASV$richness ~ alpha_tab_SituSeq$richness)
summary(mod)
plot(alpha_tab_NanoASV$richness ~ alpha_tab_SituSeq$richness)

mod <- lm(alpha_tab_NanoASV$shannon ~ alpha_tab_SituSeq$shannon)
summary(mod)
plot(alpha_tab_NanoASV$shannon ~ alpha_tab_SituSeq$shannon)

mod <- lm(alpha_tab_Nygaard$richness ~ alpha_tab_SituSeq$richness)
summary(mod)
plot(alpha_tab_Nygaard$richness ~ alpha_tab_SituSeq$richness)

mod <- lm(alpha_tab_Nygaard$shannon ~ alpha_tab_SituSeq$shannon)
summary(mod)
plot(alpha_tab_Nygaard$shannon ~ alpha_tab_SituSeq$shannon)


library(corrplot)
richness_df <- data.frame(NanoASV_richness = alpha_tab_NanoASV$richness,
                          Nygaard_richness = alpha_tab_Nygaard$richness,
                          SituSeq_richness = alpha_tab_SituSeq$richness)

cor(richness_df)

shannon_df <- data.frame(NanoASV_shannon = alpha_tab_NanoASV$shannon,
                         Nygaard_shannon = alpha_tab_Nygaard$shannon,
                         SituSeq_shannon = alpha_tab_SituSeq$shannon)

cor(shannon_df)



#MOCK.genus <- tax_glom(MOCK.phy, taxrank = "Genus")
MOCK.phylum <- tax_glom(MOCK.phy, taxrank = "Phylum")

couleurs_genus <- c("#FF0000FF","#FF9900FF","#FFCC00FF","#00FF00FF","#6699FFFF",
                    "#CC33FFFF","#99991EFF","#FF00CCFF", "#999999FF","#CC0000FF",
                    "#FFCCCCFF","#FFFF00FF","#CCFF00FF","#358000FF","#0000CCFF",
                    "#99CCFFFF","#00FFFFFF","#CCFFFFFF","#9900CCFF","#CC99FFFF",
                    "#996600FF","#79CC3DFF","#CCCCCCFF","#79CC3DFF","#CCCC99FF","black")

#Genus taxonomy ----

Genus.norm <- transform_sample_counts(MOCK.genus, function(x) x/sum(x))
Genus.melted <- psmelt(Genus.norm)
sub_Genus.melted <- Genus.melted
sub_Genus.melted[sub_Genus.melted$Abundance<0.04,11:16] <- "Z_Others"

#pdf("Figs/Pipelines_Genus_composition.pdf", he = 6, wi = 12)
jpeg("Figs/Pipelines_Genus_composition.jpg", he = 6, wi = 12, units = "in", res = 1200)
ggplot(sub_Genus.melted, aes(x = Sample, y = Abundance, fill = Genus), color = "black") +
  theme_bw() +
  
  geom_bar(stat = "identity", color = "black", size = 0.15, width = 0.85) +
  scale_fill_manual(values = couleurs_genus) + 
  ylab("Relative Abundance") +
  scale_y_continuous(expand = c(0,0)) + #
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.margin = margin(t = 10, r = 50, b = 10, l = 50, unit = "pt")) +  
  labs(title = "Nygaard dataset samples Genus-level composition")
#subtitle = "",
#caption = date()) 
dev.off()

#Phylum taxonomy ----

Phylum.norm <- transform_sample_counts(MOCK.phylum, function(x) x/sum(x))
Phylum.melted <- psmelt(Phylum.norm)
sub_Phylum.melted <- Phylum.melted
sub_Phylum.melted[sub_Phylum.melted$Abundance<0.04,7:8] <- "Z_Others"

#pdf("Figs/Pipelines_Genus_composition.pdf", he = 6, wi = 12)
#jpeg("Figs/Pipelines_Phylum_composition.jpg", he = 6, wi = 12, units = "in", res = 1200)
ggplot(sub_Phylum.melted, aes(x = Sample, y = Abundance, fill = Phylum), color = "black") +
  theme_bw() +
  
  geom_bar(stat = "identity", color = "black", size = 0.15, width = 0.85) +
  scale_fill_manual(values = c(couleurs_genus[1:6], "black")) + 
  ylab("Relative Abundance") +
  scale_y_continuous(expand = c(0,0)) + #
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.margin = margin(t = 10, r = 50, b = 10, l = 50, unit = "pt")) +  
  labs(title = "Nygaard dataset samples Genus-level composition")
#subtitle = "",
#caption = date()) 
dev.off()


# Beta div ----
library(ape)


dist.bray <- vegdist(t(MOCK.genus@otu_table), method = "bray")
pcoa.sub <- pcoa(dist.bray)
beta_tab <- cbind(pcoa.sub$vectors[,1:2], MOCK.genus@sam_data)


#Global Beta diversity ----
pdf("Figs/Beta_div.pdf", he = 8, wi =12)
ggplot(beta_tab) + theme_bw() +
  geom_point(aes(x = Axis.1, y = Axis.2, 
                 color = Pipeline), 
             size = 5) +
  geom_text(data = beta_tab, aes(x = Axis.1, y = Axis.2), label = sub("_.*$", "", rownames(beta_tab)), size = 4) +
  # scale_color_manual(breaks = waiver(), 
  #                    values = c("#69E972", "blue", "brown", "cyan")) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 0, lty = 3) +
  xlab(paste("PCo1 (", round(pcoa.sub$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub$values$Relative_eig[2]*100, 1), "%)")) +
  labs(title = "Nygaard dataset Bacterial population PCoA \nBray-Curtis distance matrix",
       #subtitle = "Plot of length by dose",
       caption = date(), 
       color = "Pipeline",
       shape = "Sample")
dev.off()

#No SituSeq BetDiv ----

NoSituSeq <- subset_samples(MOCK.genus, MOCK.genus@sam_data$Pipeline != "SituSeq")

dist.bray <- vegdist(t(NoSituSeq@otu_table), method = "bray")
pcoa.sub <- pcoa(dist.bray)
beta_tab <- cbind(pcoa.sub$vectors[,1:2], NoSituSeq@sam_data)

pdf("Figs/Beta_div_No_SituSeq.pdf", he = 8, wi =12)
ggplot(beta_tab) + theme_bw() +
  geom_point(aes(x = Axis.1, y = Axis.2, 
                 color = Pipeline),
             size = 5) +
  geom_text(data = beta_tab, aes(x = Axis.1, y = Axis.2), label = sub("_.*$", "", rownames(beta_tab)), size = 4) +
  # scale_color_manual(breaks = waiver(), 
  #                    values = c("#69E972", "blue", "brown", "cyan")) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 0, lty = 3) +
  xlab(paste("PCo1 (", round(pcoa.sub$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub$values$Relative_eig[2]*100, 1), "%)")) +
  labs(title = "Nygaard dataset Bacterial population PCoA \nBray-Curtis distance matrix \nNo SituSeq",
       caption = date(), 
       color = "Pipeline",
       shape = "Sample")
dev.off()


#Individual Beta div ----


NanoASV.beta <- subset_samples(MOCK.genus, MOCK.genus@sam_data$Pipeline == "NanoASV")

dist.bray <- vegdist(t(NanoASV.beta@otu_table), method = "bray")
dist.NanoASV <- dist.bray
pcoa.sub <- pcoa(dist.bray)
beta_tab <- cbind(pcoa.sub$vectors[,1:2], NanoASV.beta@sam_data)

pdf("Figs/NanoASV_Beta_div.pdf", he = 8, wi =12)
ggplot(beta_tab) + theme_bw() +
  geom_point(aes(x = Axis.1, y = Axis.2), 
                 color = "red", 
             size = 5) +
  geom_text(data = beta_tab, aes(x = Axis.1, y = Axis.2), label = sub("_.*$", "", rownames(beta_tab)), size = 4) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 0, lty = 3) +
  xlab(paste("PCo1 (", round(pcoa.sub$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub$values$Relative_eig[2]*100, 1), "%)")) +
  labs(title = "Nygaard dataset Bacterial population PCoA \nBray-Curtis distance matrix \nNanoASV workflow",
       caption = date(), 
       color = "Pipeline",
       shape = "Sample")
dev.off()


Nygaard.beta <- subset_samples(MOCK.genus, MOCK.genus@sam_data$Pipeline == "Nygaard")

dist.bray <- vegdist(t(Nygaard.beta@otu_table), method = "bray")
dist.Nygaard <- dist.bray
pcoa.sub <- pcoa(dist.bray)
beta_tab <- cbind(pcoa.sub$vectors[,1:2], Nygaard.beta@sam_data)

pdf("Figs/Nygaard.beta_Beta_div.pdf", he = 8, wi =12)
ggplot(beta_tab) + theme_bw() +
  geom_point(aes(x = Axis.1, y = Axis.2), 
                 color = "lightgreen", 
             size = 5) +
  geom_text(data = beta_tab, aes(x = Axis.1, y = Axis.2), label = sub("_.*$", "", rownames(beta_tab)), size = 4) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 0, lty = 3) +
  xlab(paste("PCo1 (", round(pcoa.sub$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub$values$Relative_eig[2]*100, 1), "%)")) +
  labs(title = "Nygaard dataset Bacterial population PCoA \nBray-Curtis distance matrix \nNygaard pipeline",
       caption = date(), 
       color = "Pipeline",
       shape = "Sample")
dev.off()


SituSeq.beta <- subset_samples(MOCK.genus, MOCK.genus@sam_data$Pipeline == "SituSeq")

dist.bray <- vegdist(t(SituSeq.beta@otu_table), method = "bray")
dist.SituSeq <- dist.bray
pcoa.sub <- pcoa(dist.bray)
beta_tab <- cbind(pcoa.sub$vectors[,1:2], SituSeq.beta@sam_data)

pdf("Figs/SituSeq.beta_Beta_div.pdf", he = 8, wi =12)
ggplot(beta_tab) + theme_bw() +
  geom_point(aes(x = Axis.1, y = Axis.2), 
             color = "lightblue", 
             size = 5) +
  geom_text(data = beta_tab, aes(x = Axis.1, y = Axis.2), label = sub("_.*$", "", rownames(beta_tab)), size = 4) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 0, lty = 3) +
  xlab(paste("PCo1 (", round(pcoa.sub$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub$values$Relative_eig[2]*100, 1), "%)")) +
  labs(title = "Nygaard dataset Bacterial population PCoA \nBray-Curtis distance matrix \nNygaard pipeline",
       caption = date(), 
       color = "Pipeline",
       shape = "Sample")
dev.off()


#Mantel tests to check if dissimilarity matrices are alike 

library(vegan)

mantel(dist.NanoASV, dist.Nygaard)
mantel(dist.NanoASV, dist.SituSeq)
mantel(dist.Nygaard, dist.SituSeq)
