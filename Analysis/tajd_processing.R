#process regional exome capture fst data
#merge into one dataset, plots
#created by N. A. H. Lara
#Last edit: 2021-7-1

library("tidyverse")
#set working directory
setwd("/Users/nico/Desktop/exome_project")
#Set filename variables
#regions <- c("Central", "East", "North", "PNWSpring", "PNWWinter")
regions_reor <- c("Central", "East", "PNWWinter", "PNWSpring","North")
#between_files <- list.files(path="fst", pattern = "between(.*?)weir.fst")
input_files <- list.files(path="output/taj_d", pattern = "(.*?)Tajima.D")
dir <- "output/taj_d/"

#function to read in data
#returns file ready for concatenation
readin_files <- function(i) {
	#reads in file, processes
  x <- paste(dir,i, sep="")
  z <-read.delim(x, sep = "\t")
  z <- z[z$TajimaD != "NaN", ]
  z['Region'] <- sub(".T.*", "", i)
  return(z)
}
#create initial dataframe
t_df <- readin_files(input_files[1])
#add data columns
for (i in input_files[-1]) {
	in_data <- readin_files(i)
	print(head(in_data))
	t_df <- rbind(t_df, in_data)
}

#process data into chromosome files
chromosomes <- unique(t_df$CHROM)
###hard-coded, not ideal
dir.create(paste(dir,"chrom_subsets_tajd", sep=""))
chr_num <- c(1,2,3,4,5,6,7)
for (i in chr_num) {
	#setwd("/Users/nico/Desktop/exome_project/chrom_subsets_fst")
	filter <- grep(i, chromosomes, value=TRUE)
	chr_subset <- subset(t_df, CHROM %in% filter)
	write.table(chr_subset, file=paste(dir,'chrom_subsets_tajd/CHROMS_',i,'_concat_data.tsv', sep=""), quote=FALSE, sep='\t', row.names = FALSE)
}

###Output file to .tsv
#write.table(total_dataframe, file='concat_data.tsv', quote=FALSE, sep='\t', row.names = FALSE)
#print("Done!")

ch1A <-  subset(t_df, CHROM %in% "1A")
#flat_ch1A <- melt(ch1A, id = "POS", c(colnames(ch1A[-1:-2])))
#flat_ch1A <- cbind(flat_ch1A, data.frame(do.call("rbind",strsplit(as.character(flat_ch1A$variable), "_", fixed=TRUE))))
ch1A <- ch1A %>% mutate(Region= factor(Region, levels = regions_reor)) %>% arrange(Region)
#flat_ch1A <- flat_ch1A %>% mutate(X1 = factor(X1, levels = rev(regions_reor))) %>% arrange(X1)
#fst_10_val <- sort(flat_ch1A$value, decreasing=TRUE)[length(flat_ch1A$value)*.05]
#summary(flat_ch1A)
ggplot(data=t_df, aes(BIN_START, TajimaD)) +
	geom_point(size=.01) +
	facet_grid(t_df$Region~t_df$CHROM)
	
ymin <- floor(min(t_df$TajimaD))
ymax <- ceiling(max(t_df$TajimaD))
dir.create(paste(dir,"plots",sep=""))	
for (i in chr_num) {
	filter <- grep(i, chromosomes, value=TRUE)
	chrom <- subset(t_df, CHROM %in% filter)
	chrom <- chrom %>% mutate(Region=factor(Region, levels = regions_reor)) %>% arrange(Region)
	sig_pos <- sort(chrom$TajimaD, decreasing=TRUE)[length(chrom$TajimaD)*.05]
	sig_neg <- sort(chrom$TajimaD, decreasing=FALSE)[length(chrom$TajimaD)*.05]
	ggplot(data=chrom, aes(BIN_START, TajimaD)) + 
		geom_point(size=.01, col = ifelse(chrom$TajimaD > sig_pos, "blue", ifelse(chrom$TajimaD < sig_neg, "orange", "grey"))) + 
		coord_cartesian(ylim = c(ymin,ymax)) +
		facet_grid(chrom$Region~chrom$CHROM) 
	ggsave(paste(dir,'plots/tajD_ch',i,'.png', sep=""))
}
