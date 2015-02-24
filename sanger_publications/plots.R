setwd(dir = "/Users/epigenomics_lab/Desktop/Sanger")
file_names <- dir(path = ".", pattern = "big3_")
files <- lapply(X = file_names, FUN = read.table)
file_names <- gsub(pattern = "big3_", replacement = "", x = file_names)
file_names <- gsub(pattern = ".txt", replacement = "", x = file_names)

# create progress bar
pb <- txtProgressBar(min = 0, max = length(files), style = 3)
#merge over loop
merge.master <- files[[1]]
colnames(merge.master) <- c("name",  file_names[1])
for (i in 2:length(files)) {
  tmp <- files[[i]]
  colnames(tmp) <- c("name",  file_names[i])
  merge.master <- merge(x = merge.master, y = tmp, by = c("name"), all = T)
  setTxtProgressBar(pb, i)
}

merge.master[is.na(merge.master)] <- 0

faculty_names <- read.table(file = "faculty.txt")
merge.master_faculty <- merge.master[grep(paste(faculty_names$V1, collapse='|'), merge.master$name, ignore.case=TRUE),]

library(reshape2)
merge.melt <- melt(merge.master_faculty, value.name = "count", variable.name = "year")
merge.melt <- merge.melt[(order(merge.melt$name)),]
merge.melt$year <- as.integer(as.character(merge.melt$year))
merge.melt_subset <- subset(merge.melt, year >= 2005)

merge.cumsum <- within(merge.melt_subset, {
  cumsum <- ave(x = count, name, FUN = cumsum)
  total <- ave(x = count, name, FUN = sum)
})

merge.cumsum_rank <- transform(unique(subset(merge.cumsum, select = c("name",  "total"))), rank =rank(-total, ties.method = "first"))
merge.cumsum_rank <- merge(x = merge.cumsum_rank, y = merge.cumsum, by=c("name",  "total"))

merge_cumsum_subset <- subset(merge.cumsum_rank, rank <= 10)
library(ggplot2)
ggplot(data = merge_cumsum_subset, aes(x=year, y=cumsum, group=name))+
  geom_line(aes(color=name))+
  #geom_bar(aes(fill=name) stat="identity", position="dodge")+
  xlab("Publishing year")+
  ylab("Cumulative number of publications")+
  ggtitle("Authors at the Sanger Institute with the most papers in Nature, Science or Cell journals")+
  theme_bw()

ggplot(data = merge_cumsum_subset, aes(x=year, y=count, group=name))+
  geom_line(aes(color=name))+
  #geom_bar(aes(fill=name) stat="identity", position="dodge")+
  xlab("Publishing year")+
  ylab("Yearly number of publications")+
  ggtitle("Authors at the Sanger Institute with the most papers in Nature, Science or Cell journals")+
  theme_bw()


