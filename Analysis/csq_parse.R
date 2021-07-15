#process .csq predicted consequences
#created by N. A. H. Lara
#Last edit: 2021-7-2

library("tidyverse")
#set working directory
setwd("/Users/nico/Desktop/exome_project")
#Set filename variables
regions <- c("Central", "East", "PNWWinter", "PNWSpring","North")
#between_files <- list.files(path="fst", pattern = "between(.*?)weir.fst")
dir <- "output/"
#function to read in data



t_df <- read.delim("csq.txt", sep = ";")
col <- strsplit(colnames(t_df)[1], split=";")
col
row1 <- strsplit(t_df[1], split=";")
head(t_df[1:6])
head(t_df[7:13])
head(t_df[14:18])
head(t_df[18])
t_df[2,18]

r18_df <- t_df[18]

r18_df[strsplit(colnames(r18_df), "\\.")] <- strsplit(r18_df[], "\\|")
r18_df <- strsplit(r18_df[1], "\\|")
strsplit(colnames(r18_df), "\\.")
strsplit(r18_df[1,1], "\\|")
strsplit(t_df[1,18], "\\|")
strsplit(t_df[2,18], "\\|")

f <- data.frame(a = c(1:3), b = c(4:6), c=c(7:9))
f[2,3]
f