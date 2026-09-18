source("scripts/00_setup.R")
ebox_dir <- file.path(interim_dir, "ebox_reference")
dir.create(ebox_dir, recursive = TRUE, showWarnings = FALSE)
setwd(ebox_dir)

chr = 1;

z = scan(paste0('CANNTG_',chr),what = double())

wdth = 100;

lpos = z - wdth;
rpos = z + wdth + 6;

overlap.ind = (lpos[-1] <= rpos[-length(rpos)]+1);

sum(overlap.ind)

lpos2 = lpos[ c(TRUE, !overlap.ind)];
rpos2 = rpos[ c(!overlap.ind, TRUE)];

head(cbind(lpos ,rpos ),10)
head(cbind(lpos2,rpos2),6)

cat('range =',tail(rpos2,1)-lpos2[1]+1,'\n')
cat('cover =',sum(rpos2-lpos2+1),'\n')
