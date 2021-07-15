#process regional exome capture data
#merge into one dataset, calculate summary stats
#created by N. A. H. Lara
#Last edit: 2021-5-10

#set working directory
setwd("/Users/nico/Desktop/exome_project")
#Set variables
regions <- c("ACROSS", "CENTRAL", "EAST", "NORTH", "PNWSpring", "PNWWinter")
dir <- "output/"

#function to read in data, add data columns
#returns finished file, ready for concatenation
readin_files <- function(i) {
	#reads in file, processes
  x <- paste(dir,"Processed_subsets/",i,"_subsets_stats.txt.gz", sep="")
  z <-read.delim(gzfile(x), sep = "\t")
  return(z)
}

#Calculate Fst at each position using formula: Fst = (2pq-PAa)/2pq
#q = MAF, 2pq = 2*(1-MAF)(MAF)
#PAa = 2pq from ACROSS ###CHECK THIS
Fst <- function(q, a) {
	i <- 2*(1-q)*q
	j <- 2*(1-a)*a
	f <- (i-j)/i
	return(f)
}

#create initial dataframe using total data
t_df <- readin_files(regions[1])[1:5]
#add regional data
for (i in regions) {
	region_data <- readin_files(i)
	if (all(region_data[,"POS"]==t_df[,"POS"])) {
		print(paste("Adding", i))
		regional_variables <- c("F_MISSING", "DP_SUM", "MAF")
		for (j in regional_variables) {
			t_df[paste(i, "_", j, sep="")] <- as.numeric(unlist(region_data[j]))
		}
		#add in Fst
		q <- t_df[[paste(i,"_MAF",sep="")]]
		t_df[paste(i,"_Fst", sep="")] <- Fst(q,aq)
	} else {
		print(paste(i, "failed"))
	}
}


#process data into chromosome files

chromosomes <- unique(total_dataframe$CHROM)
###hard-coded, not ideal
dir.create("chrom_subsets")
chr_num <- c(1,2,3,4,5,6,7)
for (i in chr_num) {
	setwd("/Users/nico/Desktop/exome_project/output/chrom_subsets")
	filter <- grep(i, chromosomes, value=TRUE)
	chr_subset <- subset(total_dataframe, CHROM %in% filter)
	write.table(chr_subset, file=paste('CHROMS_',i,'_concat_data.tsv', sep=""), quote=FALSE, sep='\t', row.names = FALSE)
}

###Output file to .tsv
#write.table(total_dataframe, file='concat_data.tsv', quote=FALSE, sep='\t', row.names = FALSE)
#print("Done!")
