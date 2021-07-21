rm(list = ls())
options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly=TRUE)

chr=args[1]
var.type = args[2]

#chr = "chr38"
#var.type = "snp"


outDir = "/home/buckleyrm/data/projects/gencove/process_vcfs/vcf_objects/maf.counts"
inDir = "/home/buckleyrm/data/projects/gencove/process_vcfs/vcf_objects"
sitesDir = "/home/buckleyrm/data/projects/gencove/ostrander_vcf/722_AF"
filterDir = "/home/buckleyrm/data/projects/gencove/process_vcfs/vcf_objects/venn.filter.out"

load(paste(inDir, "/impute/", chr,".impute.genotypes.", var.type, ".RData", sep = ""))
load(paste(inDir, "/wgs/", chr,".wgs.genotypes.", var.type, ".RData", sep = ""))

site.info.722 <- read.table(paste(sitesDir,"/722_info.",chr,".tsv", sep = ""), 
                            sep = "\t", header=TRUE, comment.char = "")
rownames(site.info.722) <- paste(site.info.722$X.CHROM,
                                 site.info.722$POS,
                                 site.info.722$REF,
                                 site.info.722$ALT, 
                                 sep = "_")

site.info.722$MAF <- site.info.722$AF
site.info.722$MAF[site.info.722$AF > 0.5] <- 1 - site.info.722$MAF[site.info.722$AF > 0.5]

filteredSites <- read.table(paste(filterDir,"/venn.sites.",chr,".",var.type,".tsv", sep = ""),
                            sep = "\t", header = TRUE)

# perform filtering

sites.imp.ref.wgs.all <- rownames(filteredSites)[filteredSites$imp.ref.wgs]
sites.imp.ref.wgs.filter <- intersect(rownames(gp.impute)[rowSums(gp.impute < 0.9) < 5], sites.imp.ref.wgs.all)

sites.imp.ref.wgs_imp.wgs.all <- rownames(filteredSites)[filteredSites$imp.ref.wgs | filteredSites$imp.wgs]
sites.imp.ref.wgs_imp.wgs.filter <- intersect(rownames(gp.impute)[rowSums(gp.impute < 0.9) < 5], sites.imp.ref.wgs_imp.wgs.all)




# get MAFs
maf.ref <- site.info.722[sites.imp.ref.wgs.all,"MAF"]

maf.ref.cut <- cut(maf.ref,seq(-0.01,.5,0.01))
names(maf.ref.cut) <- rownames(site.info.722[sites.imp.ref.wgs.all,])


af.imp <- (rowSums(gt.impute == "0/1") + (rowSums(gt.impute == "1/1") * 2)) / (ncol(gt.impute) * 2)
maf.imp <- af.imp
maf.imp[af.imp > .5] <- 1 - maf.imp[af.imp > .5]

maf.imp.cut <- cut(maf.imp,seq(-0.01,.5,0.01))
names(maf.imp.cut) <- names(maf.imp)

af.wgs <- (rowSums(gt.wgs == "0/1", na.rm = TRUE) + (rowSums(gt.wgs == "1/1",na.rm = TRUE) * 2)) / (ncol(gt.wgs) * 2)
maf.wgs <- af.wgs
maf.wgs[af.wgs > .5] <- 1 - maf.wgs[af.wgs > .5]

maf.wgs.cut <- cut(maf.wgs,seq(-0.01,.5,0.01))
names(maf.wgs.cut) <- names(maf.wgs)



# need to compare each gt combination
gt.comb <- c(ref = "0/0", het = "0/1", hom = "1/1")

gt.comb.df <- data.frame(gt1 = rep(gt.comb,3),
                         gt2 = rep(gt.comb,each = 3)
)
rownames(gt.comb.df) <- paste(rep(names(gt.comb),3), rep(names(gt.comb), each =3), sep = "_")

denom <- c("wgs", "imp")

maf.type <- c("ref", "wgs", "imp")

filtered <- c("all", "filter")

con.count <- NULL
var.count <- NULL
name.df <- NULL

for(gt.row in 1:nrow(gt.comb.df)){
  for(d.type in denom){
    for(m.type in maf.type){
      for(f.type in filtered){
        
        gt1 <- gt.comb.df[gt.row,"gt1"]
        gt2 <- gt.comb.df[gt.row,"gt2"]
        
        if(m.type != "ref"){
          sites0 <- get(paste("sites.imp.ref.wgs_imp.wgs.", f.type, sep = ""))
        }else{
          sites0 <- get(paste("sites.imp.ref.wgs.", f.type, sep = ""))
        }
        
        maf0 <- get(paste("maf.", m.type, ".cut", sep= ""))
        maf0 <- maf0[sites0]
        
        wgs <- gt.wgs[sites0,] == gt1
        imp <- gt.impute[sites0,] == gt2
        
	imp[is.na(wgs)] <- NA

        df0 <- data.frame(wgs & imp)
        sp.num <- split(df0, f = maf0)
        sp.num <- sapply(sp.num,colSums,na.rm = TRUE)
        
        den0 <- data.frame(get(d.type))
        sp.den <- split(den0, f = maf0)
        sp.den <- sapply(sp.den,colSums,na.rm = TRUE)
        
        con.count <- c(con.count, list(sp.num))
        var.count <- c(var.count, list(sp.den))
        
        name.df <- rbind(name.df,
                         c(
                           gt.wgs = gt1,
                           gt.imp = gt2,
                           context = d.type,
                           maf = m.type,
                           filter = f.type
                         ))
        
      }
    }
  }
}


rownames(name.df) <- paste(name.df[,"gt.wgs"], name.df[,"gt.imp"], name.df[,"context"], name.df[,"maf"], name.df[,"filter"], sep = "_")
names(con.count) <- rownames(name.df)
names(var.count) <- rownames(name.df)


site.gorups <- c("imp", "imp.ref", "imp.ref.wgs", "imp.wgs")
imp.maf.sites.all <- NULL
imp.maf.sites.filt <- NULL
for(i in site.gorups){
  sites.select.all <- rownames(filteredSites)[filteredSites[,i]]
  sites.select.filt <- intersect(rownames(gp.impute)[rowSums(gp.impute < 0.9) < 5], sites.select.all)
  imp.maf.sites.all<- rbind(imp.maf.sites.all, table(maf.imp.cut[sites.select.all]))
  imp.maf.sites.filt <- rbind(imp.maf.sites.filt, table(maf.imp.cut[sites.select.filt]))
}
rownames(imp.maf.sites.all) <- rownames(imp.maf.sites.filt) <- site.gorups
imp.maf.sites <- list(all = imp.maf.sites.all,filt = imp.maf.sites.filt)



IQS <- NULL
name.IQS <- NULL
for(m.type in maf.type){
  for(f.type in filtered){
    
    if(m.type != "ref"){
      sites0 <- get(paste("sites.imp.ref.wgs_imp.wgs.", f.type, sep = ""))
    }else{
      sites0 <- get(paste("sites.imp.ref.wgs.", f.type, sep = ""))
    }
    
    maf0 <- get(paste("maf.", m.type, ".cut", sep= ""))
    maf0 <- maf0[sites0]
    
    wgs <- gt.wgs[sites0,]
    imp <- gt.impute[sites0,]
    
    n11 <- rowSums(wgs == "0/0" & imp == "0/0", na.rm = TRUE)
    n12 <- rowSums(wgs == "0/1" & imp == "0/0", na.rm = TRUE)
    n13 <- rowSums(wgs == "1/1" & imp == "0/0", na.rm = TRUE)
    
    n21 <- rowSums(wgs == "0/0" & imp == "0/1", na.rm = TRUE)
    n22 <- rowSums(wgs == "0/1" & imp == "0/1", na.rm = TRUE)
    n23 <- rowSums(wgs == "1/1" & imp == "0/1", na.rm = TRUE)
    
    n31 <- rowSums(wgs == "0/0" & imp == "1/1", na.rm = TRUE)
    n32 <- rowSums(wgs == "0/1" & imp == "1/1", na.rm = TRUE)
    n33 <- rowSums(wgs == "1/1" & imp == "1/1", na.rm = TRUE)
    
    pO <- (n11 + n22 + n33)/ncol(wgs)
    
    a <- (n11 + n21 + n31) * (n11 + n12 + n13)
    b <- (n12 + n22 + n32) * (n21 + n22 + n23)
    c <- (n13 + n23 + n33) * (n31 + n32 + n33)
    
    pC <- (a + b + c)/(ncol(wgs)^2)
    
    IQS.val <- (pO - pC)/(1-pC)
    IQS.cut <- cut(IQS.val, seq(-1,1,.01))
    
    agg <- aggregate(IQS.val, 
                     by = list(IQS = IQS.cut, maf = maf0), 
                     FUN = length, drop = FALSE)
    agg$x[is.na(agg$x)] <- 0
    IQS <- c(IQS, 
             list(agg)
    )
    
    name.IQS <- rbind(name.IQS, data.frame(maf = m.type, filter = f.type))
  }
}

rownames(name.IQS) <- paste(name.IQS$maf, name.IQS$filter, sep = "_")
names(IQS) <- rownames(name.IQS)






save(name.df, con.count, var.count, imp.maf.sites, IQS, name.IQS, 
     file = paste(outDir,"/maf_concordance_counts.",chr,".",var.type,".RData", sep= ""))




