setwd("/scratch/vverma3/fm_ukb/magma/")
df <- read.table('clean_GWAS_summary.tab', header = T,  sep = '\t',  stringsAsFactors = F)
names(df)
# qqplot:
library("gap")

qqunif(df$pval,type="unif",logscale=TRUE,base=10,
       col=palette()[4],lcol=palette()[2],ci=T,alpha=0.05)

pdf(file = "qq2.pdf", width = 8, height = 8)
qqunif.plot(df$pval)
dev.off()

write.table(df1, file = "GWAS_sum.tab", append = F, quote = F, sep = "\t",
			eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)

