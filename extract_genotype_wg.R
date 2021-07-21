#!/usr/bin/env Rscript

rm(list = ls())
options(stringsAsFactors = FALSE)
library(vcfR)

args = commandArgs(trailingOnly=TRUE)

chr=args[1]

var.type = args[2]


imputeInDir <- paste("/home/buckleyrm/data/projects/gencove/more_dogs/vcf/chr_split/",chr,"/", sep= "")
imputeOutDir <- "/home/buckleyrm/data/projects/gencove/process_vcfs/vcf_objects/impute/"

wgsInDir <- "/home/buckleyrm/data/projects/gencove/ostrander_vcf/impute_dogs/"
wgsOutDir <- "/home/buckleyrm/data/projects/gencove/process_vcfs/vcf_objects/wgs/"

files <- list.files(imputeInDir)
files <- files[grep("vcf.gz$", files)]
samples <- gsub(".chr.+.vcf.gz", "", files)

gt.impute <- rc.impute <- ac.impute <- gp.impute <- ds.impute <- filt.impute <- NULL
for(i in 1:length(files)){
  vcf <- read.vcfR(paste(imputeInDir,files[i],sep = ""))

  if (var.type == "snp"){
    vcf <- vcf[!is.indel(vcf)]
  } else if (var.type == "indel"){
    vcf <- vcf[is.indel(vcf)]
  } else {
    stop("non valid variant type, use snp or indel")
  }

  vcf <- vcf[getALT(vcf) != "*"]

  if(i == 1){
    makeRowname <- paste(getCHROM(vcf),getPOS(vcf),getREF(vcf),getALT(vcf), sep = "_")
    rowNameDup <- duplicated(makeRowname)
  }
  
  vcf <- vcf[!rowNameDup]

  gt <- extract.gt(vcf)
  rc <- extract.gt(vcf, element = "RC", as.numeric = TRUE) 
  ac <- extract.gt(vcf, element = "AC", as.numeric = TRUE)
  gp <- extract.gt(vcf, element = "GP", as.numeric = FALSE)
  gp <- sapply(lapply(strsplit(gp, ","), as.numeric), max)
  ds <- extract.gt(vcf, element = "DS", as.numeric = TRUE)
  filt <- getFILTER(vcf)
  
  gt.impute <- cbind(gt.impute, gt) 
  rc.impute <- cbind(rc.impute, rc)
  ac.impute <- cbind(ac.impute, ac)
  gp.impute <- cbind(gp.impute, gp)
  ds.impute <- cbind(ds.impute, ds)
  filt.impute <- cbind(filt.impute, filt)
}



rownames(gt.impute) <- rownames(rc.impute) <- rownames(ac.impute) <- rownames(gp.impute) <- rownames(ds.impute) <- rownames(filt.impute) <- paste(getCHROM(vcf),getPOS(vcf),getREF(vcf),getALT(vcf), sep = "_")

colname <- samples

colnames(gt.impute) <- colnames(rc.impute) <- colnames(ac.impute) <- colnames(gp.impute) <- colnames(ds.impute) <- colnames(filt.impute) <- colname


# extract out site.info
site.info <- data.frame(CHR = getCHROM(vcf),
                        POS = getPOS(vcf),
                        REF = getREF(vcf),
                        ALT = getALT(vcf),
                        FILTER = NA,
                        AC = NA,
                        AF = NA,
                        AN = NA,
                        row.names = rownames(gt.impute))



refHomAll <- rowSums(gt.impute != "0/0") < 1 

gt.impute <- gt.impute[!refHomAll,]
rc.impute <- rc.impute[!refHomAll,]
ac.impute <- ac.impute[!refHomAll,]
gp.impute <- gp.impute[!refHomAll,]
ds.impute <- ds.impute[!refHomAll,]
filt.impute <- filt.impute[!refHomAll,]
site.info <- site.info[!refHomAll,]


# fix site.info
site.info$AN <- rowSums(!is.na(gt.impute)) * 2
site.info$AC <- rowSums(gt.impute == "0/1", na.rm = TRUE) + (rowSums(gt.impute == "1/1", na.rm = TRUE) * 2)
site.info$AF <- site.info$AC/site.info$AN

site.info.impute <- site.info


save(gt.impute,rc.impute,ac.impute,gp.impute,ds.impute,filt.impute,site.info.impute, 
     file = paste(imputeOutDir, chr, ".impute.genotypes.",var.type,".RData", sep = ""))



vcf <- read.vcfR(paste(wgsInDir, "impute_wgs.",chr,".vcf", sep = ""))

if (var.type == "snp"){
  vcf <- vcf[!is.indel(vcf)]
} else if (var.type == "indel"){
  vcf <- vcf[is.indel(vcf)]
}

vcf <- vcf[getALT(vcf) != "*"]


makeRowname <- paste(getCHROM(vcf),getPOS(vcf),getREF(vcf),getALT(vcf), sep = "_")
rowNameDup <- duplicated(makeRowname)
vcf <- vcf[!rowNameDup]

pos <- getPOS(vcf)
gt <- extract.gt(vcf)

# get rows with info only in selected samples
gt.names <- gsub("_.*","",colnames(gt))
gt <- gt[,gt.names %in% samples]
keep.rows.count <- rowSums((gt == "0/1" | gt == "1/1" | gt == "0|1" | gt == "1|1") & !is.na(gt)) > 0


# keep only sites with non-zero allele count
vcf <- vcf[keep.rows.count]

rowname <- paste(getCHROM(vcf),getPOS(vcf),getREF(vcf),getALT(vcf), sep = "_")


gt.wgs <- extract.gt(vcf, convertNA = FALSE)[,gt.names %in% samples]
ad.wgs <- extract.gt(vcf, convertNA = FALSE, element = "AD")[,gt.names %in% samples]
dp.wgs <- extract.gt(vcf, convertNA = TRUE, element = "DP", as.numeric = TRUE)[,gt.names %in% samples]
gq.wgs <- extract.gt(vcf, convertNA = TRUE, element = "GQ", as.numeric = TRUE)[,gt.names %in% samples]

rownames(gt.wgs) <- rownames(ad.wgs) <- rownames(dp.wgs) <- rownames(gq.wgs) <- rowname
colnames(gt.wgs) <- colnames(ad.wgs) <- colnames(dp.wgs) <- colnames(gq.wgs) <- gt.names[gt.names %in% samples]

gt.wgs <- (gsub("\\|", "/", gt.wgs))
gt.wgs[gt.wgs == "1/0"] <- "0/1"
gt.wgs[grepl("\\.", gt.wgs)] <- NA

site.info <- data.frame(CHR = getCHROM(vcf), 
                        POS = getPOS(vcf), 
                        REF = getREF(vcf), 
                        ALT = getALT(vcf),
                        FILTER = getFILTER(vcf),
                        AC = NA,
                        AF = NA,
                        AN = NA, 
                        row.names = rowname)

site.info$AN <- rowSums(!is.na(gt.wgs)) * 2
site.info$AC <- rowSums(gt.wgs == "0/1", na.rm = TRUE) + (rowSums(gt.wgs == "1/1", na.rm = TRUE) * 2)
site.info$AF <- site.info$AC/site.info$AN

site.info.wgs <- site.info

gt.wgs <- gt.wgs[,samples]
ad.wgs <- ad.wgs[,samples]
dp.wgs <- dp.wgs[,samples]
gq.wgs <- gq.wgs[,samples]

save(gt.wgs, ad.wgs, dp.wgs, gq.wgs, site.info.wgs, 
     file = paste(wgsOutDir, chr, ".wgs.genotypes.",var.type,".RData", sep = ""))


