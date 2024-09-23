library(rMVP)

MVP.Data(fileVCF="SAP.hecaton.PI.vcf",
         filePhe="Clean_SAP_data.tsv",
         fileKin=FALSE,
         filePC=FALSE,
         out="mvp.cnv")

df <- read.table("mvp.cnv.phe", header=TRUE)

phenotypes <- read.table("mvp.cnv.phe", header=T)
genotypes <- attach.big.matrix("mvp.cnv.geno.desc")
map <- read.table("mvp.cnv.geno.map" , head = TRUE)

# MVP.Data.PC(TRUE, mvp_prefix='mvp.vcf', pcs.keep=3)

# Covariates <- model.matrix(~as.factor(SP), data=df)

# 2:236
for (i in 2:8){
    phenotypes_ph <- phenotypes[,c(1, i)]
    imMVP <- MVP(phe=phenotypes_ph, 
                 geno=genotypes, 
                 map=map, 
                 CV.FarmCPU=df$FL,
                 threshold=0.05, 
                 method=c("MLM", "FarmCPU"))
}

