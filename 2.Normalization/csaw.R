# ============================================ csaw normalization in R ============================================ #
library(csaw)
library(edgeR)
library(dplyr)
library(tidyverse)


# prior to loading bams, bam files were aligned, mapq filtered for reads with mapq =/> 40, duplicates removed with picardtools, and blacklisted regions removed as well. 
# load metadata and create bam.files object. All bams and metata files can be located in the same folder. Bam files must also have respective indexed .bam.bai files in same folder. 
# metadata format should be (in csv):
meta_data <- read.csv("/home/kenza/PROJECTS/MetID21/csaw/MetID21_Metadata.csv", header = TRUE)
meta_data <- as.data.frame(meta_data)
bam.files <- meta_data$Path
bam.files

PolIImetadata <- meta_data[meta_data$IP == "PolII",]
vps15metadata <- meta_data[meta_data$IP == "Vps15",]
H3K4metadata <- meta_data[meta_data$IP == "H3K4",]

########################################################################################
######### PoLII ########################################################################
########################################################################################
# Counting reads into windows
frag.len <- 110
win.width <- 10
param <- readParam(minq=20)

# Estimating average fragment length
max.delay <- 500
dedup.on <- reform(param, dedup=TRUE)
x <- correlateReads(PolIImetadata$Path, max.delay, param=dedup.on)
plot(0:max.delay, x, type="l", ylab="CCF", xlab="Delay (bp)")

# per sample : 
x.list <- lapply(PolIImetadata$Path, function(bam) {
  correlateReads(bam, max.delay, param=dedup.on)
})

plot(0:max.delay, x.list[[1]], type="l", ylim=range(unlist(x.list)),
     ylab="CCF", xlab="Delay (bp)")

for(i in 2:length(x.list)) {
  lines(0:max.delay, x.list[[i]], col=i)
}

frag.len <- mean(sapply(x.list, which.max))

# Independant filtering for count data
#abundances <- aveLogCPM(asDGEList(data))
#summary(abundances)
#keep.simple <- abundances > -1
#filtered.data <- data[keep.simple,]
#summary(keep.simple)

# 1. grands bins (background)
binned <- windowCounts(PolIImetadata$Path, bin=TRUE, width=10000, ext=frag.len, param=param)

# 2. normalisation
y.bin <- asDGEList(binned)
y.bin <- calcNormFactors(y.bin)

colData(binned)$norm.factors <- y.bin$norm.factors 

# 3. petites fenêtres (signal)
data <- windowCounts(PolIImetadata$Path, width=150, ext=frag.len, param=param)

# 4. appliquer normalisation
y <- asDGEList(data)
y <- normOffsets(y, se.out=binned)

#saveRDS(y, "/home/kenza/PROJECTS/MetID21/csaw/results.v2/RangedSummarizedExperiment.rds")
y <- normOffsets(y, se.out=binned)


# visualize normalization (pink line should pass directly through blue dot cloud to represent data was normalized well)
par(mfrow=c(1, 1), mar=c(5, 4, 2, 1.5))
adj.counts <- cpm(asDGEList(binned), log=TRUE)
#normfacs <- filtered.data$norm.factors
normfacs <- y.bin$samples$norm.factors
for (i in seq_len(length(PolIImetadata$Path)-1)) {
  cur.x <- adj.counts[,1]
  cur.y <- adj.counts[,1+i]
  smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y,
                xlab="A", ylab="M", main=paste("1 vs", i+1))
  all.dist <- diff(log2(normfacs[c(i+1, 1)]))
  abline(h=all.dist, col="red")
}

PolIImetadata$norm.factors <- normfacs
write.csv(PolIImetadata, "/home/kenza/PROJECTS/MetID21/csaw/results.v2/PolII.normfact.csv" )

########################################################################################
######### Vps15 ########################################################################
########################################################################################
# Counting reads into windows
frag.len <- 110
win.width <- 10
param <- readParam(minq=20)
data <- windowCounts(vps15metadata$Path, ext=frag.len, width=win.width, param=param)
head(assay(data))

# Estimating average fragment length
max.delay <- 500
dedup.on <- reform(param, dedup=TRUE)
x <- correlateReads(vps15metadata$Path, max.delay, param=dedup.on)
plot(0:max.delay, x, type="l", ylab="CCF", xlab="Delay (bp)")

# per sample : 
x.list <- lapply(vps15metadata$Path, function(bam) {
  correlateReads(bam, max.delay, param=dedup.on)
})

plot(0:max.delay, x.list[[1]], type="l", ylim=range(unlist(x.list)),
     ylab="CCF", xlab="Delay (bp)")

for(i in 2:length(x.list)) {
  lines(0:max.delay, x.list[[i]], col=i)
}

frag.len <- mean(sapply(x.list, which.max))

# Independant filtering for count data
#abundances <- aveLogCPM(asDGEList(data))
#summary(abundances)
#keep.simple <- abundances > -1
#filtered.data <- data[keep.simple,]
#summary(keep.simple)

# 1. grands bins (background)
binned <- windowCounts(vps15metadata$Path, bin=TRUE, width=10000, ext=frag.len, param=param)

# 2. normalisation
y.bin <- asDGEList(binned)
y.bin <- calcNormFactors(y.bin)

colData(binned)$norm.factors <- y.bin$norm.factors 

# 3. petites fenêtres (signal)
data <- windowCounts(vps15metadata$Path, width=150, ext=frag.len, param=param)

# 4. appliquer normalisation
y <- asDGEList(data)
y <- normOffsets(y, se.out=binned)

#saveRDS(y, "/home/kenza/PROJECTS/MetID21/csaw/results.v2/RangedSummarizedExperiment.rds")
y <- normOffsets(y, se.out=binned)


# visualize normalization (pink line should pass directly through blue dot cloud to represent data was normalized well)
par(mfrow=c(1, 1), mar=c(5, 4, 2, 1.5))
adj.counts <- cpm(asDGEList(binned), log=TRUE)
#normfacs <- filtered.data$norm.factors
normfacs <- y.bin$samples$norm.factors
for (i in seq_len(length(vps15metadata$Path)-1)) {
  cur.x <- adj.counts[,1]
  cur.y <- adj.counts[,1+i]
  smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y,
                xlab="A", ylab="M", main=paste("1 vs", i+1))
  all.dist <- diff(log2(normfacs[c(i+1, 1)]))
  abline(h=all.dist, col="red")
}

vps15metadata$norm.factors <- normfacs
write.csv(vps15metadata, "/home/kenza/PROJECTS/MetID21/csaw/results.v2/Vps15.normfact.csv" )


########################################################################################
######### H3K4me3 ######################################################################
########################################################################################
# Counting reads into windows
frag.len <- 110
win.width <- 10
param <- readParam(minq=20)
data <- windowCounts(H3K4metadata$Path, ext=frag.len, width=win.width, param=param)
head(assay(data))

# Estimating average fragment length
max.delay <- 500
dedup.on <- reform(param, dedup=TRUE)
x <- correlateReads(H3K4metadata$Path, max.delay, param=dedup.on)
plot(0:max.delay, x, type="l", ylab="CCF", xlab="Delay (bp)")

# per sample : 
x.list <- lapply(H3K4metadata$Path, function(bam) {
  correlateReads(bam, max.delay, param=dedup.on)
})

plot(0:max.delay, x.list[[1]], type="l", ylim=range(unlist(x.list)),
     ylab="CCF", xlab="Delay (bp)")

for(i in 2:length(x.list)) {
  lines(0:max.delay, x.list[[i]], col=i)
}

frag.len <- mean(sapply(x.list, which.max))

# Independant filtering for count data
#abundances <- aveLogCPM(asDGEList(data))
#summary(abundances)
#keep.simple <- abundances > -1
#filtered.data <- data[keep.simple,]
#summary(keep.simple)

# 1. grands bins (background)
binned <- windowCounts(H3K4metadata$Path, bin=TRUE, width=10000, ext=frag.len, param=param)

# 2. normalisation
y.bin <- asDGEList(binned)
y.bin <- calcNormFactors(y.bin)

colData(binned)$norm.factors <- y.bin$norm.factors 

# 3. petites fenêtres (signal)
data <- windowCounts(H3K4metadata$Path, width=150, ext=frag.len, param=param)

# 4. appliquer normalisation
y <- asDGEList(data)
y <- normOffsets(y, se.out=binned)

#saveRDS(y, "/home/kenza/PROJECTS/MetID21/csaw/results.v2/RangedSummarizedExperiment.rds")
y <- normOffsets(y, se.out=binned)


# visualize normalization (pink line should pass directly through blue dot cloud to represent data was normalized well)
par(mfrow=c(1, 1), mar=c(5, 4, 2, 1.5))
adj.counts <- cpm(asDGEList(binned), log=TRUE)
#normfacs <- filtered.data$norm.factors
normfacs <- y.bin$samples$norm.factors
for (i in seq_len(length(H3K4metadata$Path)-1)) {
  cur.x <- adj.counts[,1]
  cur.y <- adj.counts[,1+i]
  smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y,
                xlab="A", ylab="M", main=paste("1 vs", i+1))
  all.dist <- diff(log2(normfacs[c(i+1, 1)]))
  abline(h=all.dist, col="red")
}

H3K4metadata$norm.factors <- normfacs
write.csv(H3K4metadata, "/home/kenza/PROJECTS/MetID21/csaw/results.v2/H3K4me3.normfact.csv" )












