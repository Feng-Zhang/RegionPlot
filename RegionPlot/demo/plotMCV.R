
data(MCV_pval)
png("MCV.png", width=3200, height=2200,res=300)
plotRegion(MCV_pval)
dev.off()
