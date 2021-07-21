rm(list = ls())
options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly=TRUE)

chr=args[1]
var.type = args[2]

outDir = "/home/buckleyrm/data/projects/gencove/process_vcfs/vcf_objects/venn.filter.out"
inDir = "/home/buckleyrm/data/projects/gencove/process_vcfs/vcf_objects"



load(paste(inDir, "/impute/", chr,".impute.genotypes.", var.type, ".RData", sep = ""))
load(paste(inDir, "/wgs/", chr,".wgs.genotypes.", var.type, ".RData", sep = ""))

name722 = paste("/home/buckleyrm/data/projects/gencove/ostrander_vcf/722_AF/722_info.", chr, ".tsv", sep= "")  

site.info.722 <- read.table(name722, 
                            sep = "\t", header=TRUE, comment.char = "")
rownames(site.info.722) <- paste(site.info.722$X.CHROM,
                                 site.info.722$POS,
                                 site.info.722$REF,
                                 site.info.722$ALT, 
                                 sep = "_")


varGroups <- list(impute = rownames(site.info.impute),
                  wgs = rownames(site.info.wgs),
                  refPop = rownames(site.info.722)[rownames(site.info.722) %in% rownames(site.info.wgs) | rownames(site.info.722) %in% rownames(site.info.impute)])


# count intersects

imp <- varGroups$impute[!(varGroups$impute %in% varGroups$wgs | varGroups$impute %in% varGroups$refPop)]
imp.ref <- varGroups$impute[(varGroups$impute %in% varGroups$refPop) & !(varGroups$impute %in% varGroups$wgs)]
imp.ref.wgs <- varGroups$impute[varGroups$impute %in% varGroups$wgs & varGroups$impute %in% varGroups$refPop]

wgs <- varGroups$wgs[!(varGroups$wgs %in% varGroups$impte | varGroups$wgs %in% varGroups$refPop)]
wgs.ref <- varGroups$wgs[(varGroups$wgs %in% varGroups$refPop) & !(varGroups$wgs %in% varGroups$impute)]
wgs.ref.imp <- varGroups$wgs[varGroups$wgs %in% varGroups$impute & varGroups$wgs %in% varGroups$refPop]

imp.wgs <- varGroups$impute[varGroups$impute %in% varGroups$wgs & !(varGroups$impute %in% varGroups$refPop)]

# create merged tale of all sites and intersects

# create intersect table

venn_count <- data.frame(group = c("wgs", "wgs&ref", "imp", 
                                   "imp&ref", "imp&ref&wgs", "imp&wgs"),
                         count = c(length(wgs), length(wgs.ref), length(imp), 
                                   length(imp.ref),length(imp.ref.wgs), length(imp.wgs))
                         )


allSitesRows <- unique(union(union(varGroups$impute, varGroups$wgs), varGroups$refPop))
allSitesMat <- matrix(data = FALSE, nrow = length(allSitesRows), ncol = 6, 
                      dimnames = list(allSitesRows, gsub("\\&", "\\.", venn_count$group)))
for(i in colnames(allSitesMat)){
  allSitesMat[get(i),i] <- TRUE
}


# count per sample gt per intersect

venn_colSums <- data.frame(
  imp_ref = colSums(gt.impute[imp,] == "0/0", na.rm = TRUE),
  imp.ref_ref = colSums(gt.impute[imp.ref,] == "0/0", na.rm = TRUE),
  imp.wgs.ref_ref = colSums(gt.impute[imp.ref.wgs,] == "0/0", na.rm = TRUE),
  imp.wgs_ref = colSums(gt.impute[imp.wgs,] == "0/0", na.rm = TRUE),
  wgs_ref = colSums(gt.wgs[wgs,] == "0/0", na.rm = TRUE),
  wgs.ref_ref = colSums(gt.wgs[wgs.ref,] == "0/0", na.rm = TRUE),
  wgs.imp.ref_ref = colSums(gt.wgs[wgs.ref.imp,] == "0/0", na.rm = TRUE),
  wgs.imp_ref = colSums(gt.wgs[imp.wgs,] == "0/0", na.rm = TRUE),
  
  imp_het = colSums(gt.impute[imp,] == "0/1", na.rm = TRUE),
  imp.ref_het = colSums(gt.impute[imp.ref,] == "0/1", na.rm = TRUE),
  imp.wgs.ref_het = colSums(gt.impute[imp.ref.wgs,] == "0/1", na.rm = TRUE),
  imp.wgs_het = colSums(gt.impute[imp.wgs,] == "0/1", na.rm = TRUE),
  wgs_het = colSums(gt.wgs[wgs,] == "0/1", na.rm = TRUE),
  wgs.ref_het = colSums(gt.wgs[wgs.ref,] == "0/1", na.rm = TRUE),
  wgs.imp.ref_het = colSums(gt.wgs[wgs.ref.imp,] == "0/1", na.rm = TRUE),
  wgs.imp_het = colSums(gt.wgs[imp.wgs,] == "0/1", na.rm = TRUE),
  
  imp_hom = colSums(gt.impute[imp,] == "1/1", na.rm = TRUE),
  imp.ref_hom = colSums(gt.impute[imp.ref,] == "1/1", na.rm = TRUE),
  imp.wgs.ref_hom = colSums(gt.impute[imp.ref.wgs,] == "1/1", na.rm = TRUE),
  imp.wgs_hom = colSums(gt.impute[imp.wgs,] == "1/1", na.rm = TRUE),
  wgs_hom = colSums(gt.wgs[wgs,] == "1/1", na.rm = TRUE),
  wgs.ref_hom = colSums(gt.wgs[wgs.ref,] == "1/1", na.rm = TRUE),
  wgs.imp.ref_hom = colSums(gt.wgs[wgs.ref.imp,] == "1/1", na.rm = TRUE),
  wgs.imp_hom = colSums(gt.wgs[imp.wgs,] == "1/1", na.rm = TRUE)
)

# count AF per intersect

venn_AF <- data.frame(
  imp = hist(site.info.impute[imp,"AF"] * 100, breaks = 0:100, plot = FALSE)$counts,
  imp.ref = hist(site.info.impute[imp.ref,"AF"] * 100, breaks = 0:100, plot = FALSE)$counts,
  imp.ref.wgs = hist(site.info.impute[imp.ref.wgs,"AF"] * 100, breaks = 0:100, plot = FALSE)$counts,
  imp.wgs = hist(site.info.impute[imp.wgs,"AF"] * 100, breaks = 0:100, plot = FALSE)$counts,
  
  wgs = hist(site.info.wgs[wgs,"AF"] * 100, breaks = 0:100, plot = FALSE)$counts,
  wgs.ref = hist(site.info.wgs[wgs.ref,"AF"] * 100, breaks = 0:100, plot = FALSE)$counts,
  wgs.ref.imp = hist(site.info.wgs[wgs.ref.imp,"AF"] * 100, breaks = 0:100, plot = FALSE)$counts,
  wgs.imp = hist(site.info.wgs[imp.wgs,"AF"] * 100, breaks = 0:100, plot = FALSE)$counts,
  
  ref.imp = hist(site.info.722[imp.ref,"AF"] * 100, breaks = 0:100, plot = FALSE)$counts,
  ref.wgs = hist(site.info.722[wgs.ref,"AF"] * 100, breaks = 0:100, plot = FALSE)$counts,
  ref.wgs.imp = hist(site.info.722[wgs.ref.imp,"AF"] * 100, breaks = 0:100, plot = FALSE)$counts
  )



# write venn table

write.table(venn_count,
            file = paste(outDir, "/","venn.count.",chr,".", var.type, ".tsv", sep = ""),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(venn_colSums,
            file = paste(outDir, "/","venn.count.sample.",chr,".", var.type, ".tsv", sep = ""),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(venn_AF,
            file = paste(outDir, "/","venn.count.af.",chr,".", var.type, ".tsv", sep = ""),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(allSitesMat,
            file = paste(outDir, "/","venn.sites.",chr, ".", var.type, ".tsv", sep = ""),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


#### how do we treat filtering?
# do we ask based on acuracy how do we treat filtering. 
# do we look at incorrectly identifying gt and what filtering level to apply


# filtering cutoff
# enrichment of correct genotypes
df_table <- NULL
correct_gt = as.logical(gt.impute[imp.ref.wgs,] == gt.wgs[imp.ref.wgs,]) 
for(i in seq(.7,1,by = 0.01)){
  df <- data.frame(
  correct_gt = correct_gt,
  threshold = as.logical(gp.impute[imp.ref.wgs,] >= i)
  )
  tab <- table(df)[2:1,2:1]
  df_table <- rbind(df_table,c(i,tab[1:4]))
}
colnames(df_table) <- c("GP.PASS","PASS.TRUE","PASS.FALSE","FAIL.TRUE", "FAIL.FALSE")




df_table_rows <- NULL
below_dist <- NULL

correct_gt = gt.impute[imp.ref.wgs,] == gt.wgs[imp.ref.wgs,] 
for(i in seq(.7,1,by = 0.05)){
  
  below_threshold <- rowSums(!(gp.impute[imp.ref.wgs,] >= i))
  
  below_dist <- cbind(below_dist, 
                      hist(below_threshold, breaks = 0:ncol(gt.impute), plot = FALSE)$count)
  
  for(j in c(1,2,3,4,5,10,15,30)){
    
    threshold_mat <- matrix(TRUE,nrow = nrow(correct_gt), ncol = ncol(correct_gt),
                        dimnames = list(rownames(correct_gt), colnames(correct_gt)))
    
    threshold_mat[below_threshold >= j,] <- FALSE
    
    df <- data.frame(
      correct_gt = as.logical(correct_gt),
      threshold = as.logical(threshold_mat)
    )
    tab <- table(df)
    tab.mat <- matrix(NA, nrow = 2, ncol = 2, 
                      dimnames = list(c("TRUE", "FALSE"), 
                                      c("TRUE", "FALSE"))
                      )
    tab.mat[1:nrow(tab),1:ncol(tab)] <- tab[nrow(tab):1,ncol(tab):1]
    df_table_rows <- rbind(df_table_rows,c(i,j,tab.mat[1:4]))
  }
  
}
colnames(df_table_rows) <- c("GP.PASS","MIN.FAIL","PASS.TRUE","PASS.FALSE","FAIL.TRUE", "FAIL.FALSE")
colnames(below_dist) <- paste("PASS_",seq(.7,1,by = 0.05), sep= "")

write.table(df_table,
            file = paste(outDir, "/", "filter.val.site.", chr, ".", var.type, ".tsv", sep = ""),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(df_table_rows,
            file = paste(outDir, "/", "filter.val.row.", chr, ".", var.type, ".tsv", sep = ""),
            sep = "\t", quote = FALSE, row.names = FALSE)



