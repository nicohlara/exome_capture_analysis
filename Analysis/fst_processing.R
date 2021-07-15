#process regional exome capture fst data
#merge into one dataset, plots
#created by N. A. H. Lara
#Last edit: 2021-7-8

#load in libraries
library("tidyverse")
library("reshape2")
library("ggVennDiagram")

#set working directory
setwd("/Users/nico/Desktop/exome_project")

#Set filename variables
regions <- c("Central", "East", "PNWWinter", "PNWSpring","North")
fst_files <- list.files(path="output/imp_fst", pattern = "(.*?)weir.fst")

#set output directory
dir <- "output/imp_fst/"

#function to read in data file for manipulation
readin_fst_files <- function(i) {
  x <- paste(dir,i, sep="")
  z <-read.delim(x, sep = "\t")
  return(z)
}

#create initial dataframe of chromosome+position
t_df <- readin_fst_files(fst_files[1])[1:2]
#add data columns
for (i in fst_files) {
	fst_data <- readin_fst_files(i)
	if (all(fst_data[,"POS"]==t_df[,"POS"])) {
		print(paste("Adding", i))
		j <- "WEIR_AND_COCKERHAM_FST"
		k <- gsub("([_])|[[:punct:]]","\\1",i)
		l <- gsub("weirfst","", k)
		t_df[l] <- as.numeric(unlist(fst_data[j]))
	} else {
		print(paste(i, "failed"))
	}
}
#process data into chromosome files
chromosomes <- unique(t_df$CHROM)
dir.create(paste(dir,"chrom_subsets_fst", sep=""))
for (i in 1:7) {
	filter <- grep(i, chromosomes, value=TRUE)
	chr_subset <- subset(t_df, CHROM %in% filter)
	write.table(chr_subset, file=paste(dir,'chrom_subsets_fst/CHROMS_',i,'_concat_data.tsv', sep=""), quote=FALSE, sep='\t', row.names = FALSE)
}
	
	
dir.create(paste(dir,"plots",sep=""))	
for (i in chromosomes) {
	chrom <- subset(t_df, CHROM %in% i)
	flat_chrom <- melt(chrom, id = "POS", c(colnames(chrom[-1:-2])))
	flat_chrom <- cbind(flat_chrom, data.frame(do.call("rbind", strsplit(as.character(flat_chrom$variable), "_", fixed=TRUE))))
	flat_chrom <- flat_chrom %>% mutate(X2 = factor(X2, levels = regions)) %>% arrange(X1)
	flat_chrom <- flat_chrom %>% mutate(X1 = factor(X1, levels = rev(regions))) %>% arrange(X1)
	fst_sig <- sort(flat_chrom$value, decreasing=TRUE)[length(flat_chrom$value)*.05]
	ggplot(data=flat_chrom, aes(POS, value)) + 
		geom_point(size=.01, col = ifelse(flat_chrom$value > fst_sig, "blue", ifelse(flat_chrom$value < 0, "black", "grey"))) +
		facet_grid(flat_chrom$X1~flat_chrom$X2)
	ggsave(paste(dir,'plots/fst',i,'.png', sep=""))
}

#plot specific region/chromosome combos
fst_plot <- function(chromosome, region, file_name) {
	chrom <- subset(t_df, CHROM %in% chromosome)
	flat_chrom <- melt(chrom, id = "POS", c(colnames(chrom[-1:-2])))
	flat_chrom <- subset(flat_chrom, variable %in% region)
	print(unique(flat_chrom$variable))
	flat_chrom <- cbind(flat_chrom, data.frame(do.call("rbind", strsplit(as.character(flat_chrom$variable), "_", fixed=TRUE))))
	flat_chrom <- flat_chrom %>% mutate(X2 = factor(X2, levels = regions)) %>% arrange(X1)
	flat_chrom <- flat_chrom %>% mutate(X1 = factor(X1, levels = rev(regions))) %>% arrange(X1)
	fst_sig <- sort(flat_chrom$value, decreasing=TRUE)[length(flat_chrom$value)*.05]
	ggplot(data=flat_chrom, aes(POS, value)) + 
		geom_point(size=.01, col = ifelse(flat_chrom$value > fst_sig, "blue", ifelse(flat_chrom$value < 0, "black", "grey"))) +
		facet_grid(flat_chrom$X1~flat_chrom$X2)
	ggsave(paste(dir,'plots/',file_name,'.png', sep=""))
}	
fst_plot("6A", "East", "6A_fst")	



#output truth table of top 5% selection of snps for each region/chromosome combination
venn_df <- data.frame(LOC = c(""))
venn_list <- list()
for (j in regions) {
	loc_df <- data.frame(LOC=c(""))
	for (i in chromosomes) {
		chrom <- subset(t_df, CHROM %in% i)
		chrom <- chrom[, c("CHROM", "POS", j)]
		sig_cut <- sort(chrom[[j]], decreasing=TRUE)[length(chrom[[j]])*.05]
		chrom <- subset(chrom, eval(as.name(paste(j))) >= sig_cut)
		chrom <- data.frame(LOC=paste(chrom$CHROM, chrom$POS, sep="_"))
		loc_df <- rbind(loc_df, chrom)
	}
	loc_df[j] <- c("TRUE")
	print(j)
	venn_list <- c(venn_list, j = list(loc_df$LOC))
	venn_df <- merge(venn_df, loc_df, by = "LOC", all=TRUE)
}
names(venn_list) <- regions
ggVennDiagram(venn_list)
ggsave(paste(dir,'/gene_plots/venn_sig_fst', '.png', sep=""))


#KASP gene subset
g1 <- "TraesCS6A01G091100"
g1b <- "TraesCS6A01G091200"
#this is the gene to focus on
g1c <- "TraesCS6A01G091300"
#position: chr6A:60019880..60023759
###High-confidence annotations v1.1:
#TraesCS6A02G091300
#chr6A:60019299..60023891 (+ strand)

g1_chrom <- "6A"
#overarching region range
#g1_range <- c(59710633,60092778)
#gene range
g1_range <- c(60019299, 60023891)
#subset by location
g1_SNPs <- subset(t_df, POS >= g1_range[1] & POS <= g1_range[2] & CHROM == g1_chrom)

dir.create(paste(dir,"gene_plots",sep=""))	
flat_gene <- melt(g1_SNPs, id = "POS", c(colnames(chrom[-1:-2])))
flat_gene <- cbind(flat_gene, data.frame(do.call("rbind", strsplit(as.character(flat_gene$variable), "_", fixed=TRUE))))
flat_gene <- flat_gene %>% mutate(X2 = factor(X2, levels = regions)) %>% arrange(X1)
flat_gene <- flat_gene %>% mutate(X1 = factor(X1, levels = rev(regions))) %>% arrange(X1)
ggplot(data=flat_gene, aes(POS, value)) + 
	geom_point(size=.01) +
	facet_grid(flat_gene$X1~flat_gene$X2)
	ggsave(paste(dir,'/gene_plots/','gene','.png', sep=""))

