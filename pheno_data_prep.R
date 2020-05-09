cd /scratch/vivek22/FM_UKB/pheno/
module load gcc/7.3.0 r/3.6.1
R --no-save
# load original data:
df <- readRDS("FM_pheno_cov.RDS")
summary(factor(df$FM))
# 1826 cases / 197050 cntls
table(df$MALE,df$FM)
# only 206 cases are males, drop 'em
df <- df[df$MALE == 0,]
# get whit brit data:
eth <- read.table("/scratch/vivek22/UKB/ukb_sqc_v2.txt")
for(i in 1:40){
		  names(eth)[i + 25] <- paste0("PC",i)
		}

		eth <- eth[,-c(2,4:9,12:18,21,22,25,66:68)]
		names(eth)[1:7] <- c("array","sex","inf_sex","het_miss_outlr","sex_aneup","exces_relatvs","white_brit")
eth <- 

