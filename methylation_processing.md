Do this interactively in a screen
```
screen -S charrmethylR
R
```

Set up R environments

```
library(tidyverse)
library(edgeR)
```

Load methylation data

```
MethylAgeMeta <- read.delim("CharrMethyl_MBB22_23_KL_19.tsv", stringsAsFactors=FALSE)
files <- MethylAgeMeta$File
methylall <- readBismark2DGE(files, sample.names=MethylAgeMeta$ID)
table(methylall$genes$Chr)
```

Compare methylation counts and coverage

```
Methylation <- gl(2,1,ncol(methylall), labels=c("Me","Un"))
Me <- methylall$counts[, Methylation=="Me"]
Un <- methylall$counts[, Methylation=="Un"]
Coverage <- Me + Un
```

Filter for loci with mean coverage ≥ 8

```
MeanCoverage <- rowMeans(Coverage)
HasMeanCoverage <- MeanCoverage >= 8
```

Apply the filter to methylation data

```
highcov <- methylall[HasMeanCoverage,, keep.lib.sizes=FALSE]
```

Calculate M values

```
Me <- highcov$counts[, Methylation=="Me"]
Un <- highcov$counts[, Methylation=="Un"]
M <- data.frame(log2(Me + 2) - log2(Un + 2))
```

Save data

```
write.table(M, "CharrMBBKL.Mmatrix.tsv", sep = "\t", quote = F, row.names = T, col.names = T)
```
