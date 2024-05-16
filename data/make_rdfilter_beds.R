# remove all objects from the workspace
rm(list=ls())

library(tidyverse)

folder_name <- 'variant_filtering/depth'

# read files that have '.regions.bed' extension in folder as a list
files <- list.files(path = folder_name, pattern = '.regions.bed$', full.names = TRUE)

# create a data frame where each column is a file
# and each row is a depth value by reading each file's fourth column
df <- data.frame()

for (file in files) {
  sample <- strsplit(basename(file), '\\.')[[1]][1]
  data <- read.table(file, sep = '\t', header = FALSE)[[4]]
  if (ncol(df) == 0) {
    df <- data.frame(sample = data)
    names(df) <- sample
  } else {
    df[[sample]] <- data
  }
}

# add a column that contains the chromosome, one the start position and one the end position
df$chrom <- read.table(files[1], sep = '\t', header = FALSE)[[1]]
df$start <- read.table(files[1], sep = '\t', header = FALSE)[[2]]
df$end <- read.table(files[1], sep = '\t', header = FALSE)[[3]]

# remove rows that don't have "chr" in the chromosome column
# but not ChrX and ChrY
df <- df[grepl('Chr', df$chrom), ]
df <- df[!grepl('ChrX', df$chrom), ]
df <- df[!grepl('ChrY', df$chrom), ]

# get a new data frame with the sum of each row of columns meeting
# a condition in their name to group them by population:

# ll_ba = bal (Balkan Lynx)
# ll_ca = cau (Caucasian Lynx)
# ll_cr = crp (Carpathian Lynx)
# ll_la = lva (Lynx in Latvia)
# ll_no = nor (Lynx in Norway)
# ll_po = pol (Lynx in NE Poland)
# ll_tu = tva (Lynx in Tuva region, Russia)
# ll_ka or ll_og or ll_to = mng (Lynx in Mongolia)
# ll_ki or ll_ur = wru (Lynx in Western Russia: Kirov and Ural regions)
# ll_ya = yak (Lynx in Yakutia, Russia)
# ll_vl = pyk (Lynx in Primorsky Krai, Russia)

populations <- list(
  bal = names(df)[grepl('ll_ba', names(df))],
  cau = names(df)[grepl('ll_ca', names(df))],
  crp = names(df)[grepl('ll_cr', names(df))],
  mng = names(df)[grepl('ll_ka|ll_og|ll_to', names(df))],
  wru = names(df)[grepl('ll_ki|ll_ur', names(df))],
  lva = names(df)[grepl('ll_la', names(df))],
  nor = names(df)[grepl('ll_no', names(df))],
  pol = names(df)[grepl('ll_po', names(df))],
  tva = names(df)[grepl('ll_tu', names(df))],
  pyk = names(df)[grepl('ll_vl', names(df))],
  yak = names(df)[grepl('ll_ya', names(df))]
)

populations_df <- data.frame()
# for (pop in names(populations)) {
#   populations_df[[pop]] <- rowSums(df[, populations[[pop]]])
# }

for (pop in names(populations)) {
  data <- rowSums(df[, populations[[pop]]])
  if (ncol(populations_df) == 0) {
    populations_df <- data.frame(pop = data)
    names(populations_df) <- pop
  } else {
    populations_df[[pop]] <- data
  }
}

# draw a histogram of each population
# showing only values between and 2 times the mode
# add a vertical line at the mode and one at 1.5 times the mode
for (pop in names(populations_df)) {
  hist(populations_df[[pop]], breaks = 100, main = paste('depth distribution of', pop),
       xlim = c(0, 2.5 * mode(populations_df[[pop]])),
       xlab = 'Depth', col = 'lightblue')
  abline(v = mode(populations_df[[pop]]), col = 'black', lty = 'dashed', lwd = 1)
  abline(v = 1.5 * mode(populations_df[[pop]]), col = 'red', lty = 'dashed', lwd = 1)
  png(filename = paste(folder_name, paste(pop, 'depth_distribution.png', sep = '_'), sep = '/'))
  dev.off()
}

####################################################################################################
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#### this one below works

library(ggplot2)


# Create plot
for (pop in names(populations_df)) {
p <- ggplot(data.frame(Depth = populations_df[[pop]]), aes(x = Depth)) +
  geom_histogram(bins = 100, fill = 'lightblue', color = 'grey80', alpha = 0.5) +
  geom_vline(aes(xintercept = get_mode(populations_df[[pop]])), color = 'deepskyblue2', linetype = 'dashed', linewidth  = 0.4) +
  geom_vline(aes(xintercept = 1.5 * get_mode(populations_df[[pop]])), color = 'firebrick1', linetype = 'dashed', linewidth  = 0.4) +
  xlim(1.3 * get_mode(populations_df[[pop]]), 2.5 * get_mode(populations_df[[pop]])) +
  labs(x = 'Depth', title = paste('Depth distribution of', pop)) +
  theme_minimal()

ggsave(filename = paste(folder_name, paste(pop, 'zoomed_1,3_mode_depth_distribution.png', sep = '_'), sep = '/'), bg = "white", plot = p)

}

# create a new data frame with one row for each region and one column for each population
# plus the chromosome, start and end columns
# the value of each cell is either pass or fail depending on the condition
# that the depth is between 0 and 1.5 times the mode
regions_df <- data.frame(chrom = df$chrom, start = df$start, end = df$end)

# modes <- get_mode(populations_df) * 1.5
# for (pop in names(populations_df)) {
#   regions_df[[pop]] <- ifelse(populations_df[[pop]] <= modes[[pop]], 'pass', 'fail')
# }

for (pop in names(populations_df)) {
  mode_val <- get_mode(populations_df[[pop]]) * 1.5
  regions_df[[pop]] <- ifelse(populations_df[[pop]] <= mode_val, 'pass', 'fail')
}

# write the data frame to a file
write.table(regions_df, file = paste(folder_name, 'regions_depth_filtering.tsv', sep = '/'), sep = '\t', row.names = FALSE)

# write a bed file for each population with the regions that did not pass the filter
for (pop in names(populations_df)) {
  write.table(regions_df[regions_df[[pop]] == 'fail', c('chrom', 'start', 'end')],
              file = paste(folder_name, paste(pop, 'rd_filter.bed', sep = '.'), sep = '/'),
              sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE) # in the original code, i did not have the quote = FALSE, so the chromosomes were in "" and the filtering didn't work. 
}
}


# print the amount of regions that pass and fail in each population
for (pop in names(populations_df)) {
  n_fail <- sum(regions_df[[pop]] == 'fail')
  print(paste(pop, 'fail:', n_fail))
}

# print the amount of regions that have a value of fail in all populations
n_fail_all <- sum(apply(regions_df[names(populations_df)], 1, function(x) all(x == 'fail')))
print(paste('all fail:', n_fail_all))

# print the amount of regions that are fail in 'bal' and pass in the other populations
n_Fbal_Pothers <- sum(regions_df$bal == 'fail' & apply(regions_df[names(populations_df) != 'bal'], 1, function(x) all(x == 'pass')))
print(paste('bal fail and others pass:', n_Fbal_Pothers))
