
source("~/bin/R/functions/stringUtils.R")

bb <- read.table("../data/20080606_1348_B_BMM_LPS_0060_P50.bed",sep="\t",as.is=TRUE,header=TRUE)

# should be 5 plus multiples of 5

for ( i in 2:10 ){
  toks <- split1(bb[i],'\t')  
  ll <- length(toks)
  nfmax <- floor(ll/5)-1
  cat("Ntoks:", ll,nfmax,"\n")
  nf <- 0 
  for ( j in 1:nfmax ){
    feature.type <- toks[5*j+1]
    cat(feature.type,",")
    if (feature.type == "Transcript" | feature.type=="miRNA"){
      nf <- nf+1
    }
  }
  cat("\n")  
  cat("no.features:",nf,"\n")
}

