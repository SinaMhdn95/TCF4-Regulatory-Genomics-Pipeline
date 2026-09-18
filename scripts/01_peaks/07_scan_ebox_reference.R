
#source("http://bioconductor.org/biocLite.R")
#biocLite("BSgenome.HSapiens.UCSC.hg19")

library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens

source("scripts/00_setup.R")
workdir <- file.path(interim_dir, "ebox_reference")
dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
setwd(workdir)

#Prep search for sequence string CANNTG
Ebox <- DNAStringSet(c(	"CAAATG",
				"CAACTG",
				"CAAGTG",
				"CAATTG",
				"CACATG",
				"CACCTG",
				"CACGTG",
				"CACTTG",
				"CAGATG",
				"CAGCTG",
				"CAGGTG",
				"CAGTTG",
				"CATATG",
				"CATCTG",
				"CATGTG",
				"CATTTG" ))
pdict0 <- PDict(Ebox)

#Count instances of string
params <- new("BSParams", X = Hsapiens, FUN = countPDict, simplify=T,
              exclude = c("M","random","hap","Un"))
matches = bsapply(params, pdict = pdict0)
matches
sum(matches)

params@FUN <- matchPDict
mapmatches = bsapply(params, pdict = pdict0)

nchrom <- 24
chr_list <- c(seq(1,nchrom-2), "X", "Y")
#i<-1
for(i in 1:length(chr_list))
{
	tmp <- as.data.frame(eval(parse(text=paste("mapmatches$chr",chr_list[i],sep=""))))
	colnames(tmp) <- c("string", "start", "end", "width")
	tmp1 <- tmp[order(tmp$start),]
	tmp2 <- tmp1[(!duplicated(tmp1)),]
	write.table(tmp2$start, file=paste("CANNTG_",chr_list[i],sep=""), row.names=F,col.names=F,quote=F,append=F, sep=",")
}
