## https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
library(fgsea)
library(data.table)
library(ggplot2)

## use rank file and pathways files
#rnk.file <- system.file("extdata", "naive.vs.th1.rnk", package="fgsea")
#gmt.file <- system.file("extdata", "mouse.reactome.gmt", package="fgsea")
#rnk.file
#gmt.file

## read rank file and pathway files
rnk.file <- "./naive.vs.th1.rnk"
gmt.file <- "./mouse.reactome.gmt"

# Loading ranks

ranks <- read.table(rnk.file, header=TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$t, ranks$ID)
#str(ranks)

# Loading pathways

pathways <- gmtPathways(gmt.file)
#str(head(pathways))

# run fgsea
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500)
head(fgseaRes)

### save result
fwrite(fgseaRes, file="fgseaRes_file.tsv", sep="\t", sep2=c("", " ", ""))

pdf(file="plot2.pdf")

## print plot
plotEnrichment(pathways[["5991130_Programmed_Cell_Death"]],
               ranks) + labs(title="Programmed Cell Death")

plotEnrichment(pathways[["5990980_Cell_Cycle"]],
               ranks) + labs(title="Cell Cycle")
dev.off()

### collapsed pathways
pdf(file="topUpDown.collapsed.pathways2.pdf",height=12,width=18)

collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01],
                                      pathways, ranks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
                         order(-NES), pathway]
plotGseaTable(pathways[mainPathways], ranks, fgseaRes,
              gseaParam = 0.5)

dev.off()
