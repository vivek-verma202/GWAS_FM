#!/bin/bash
mkdir /scratch/vivek22/UKB
mkdir /scratch/vivek22/UKB/geno
mkdir /scratch/vivek22/FM_UKB

# Generic UKB data prep
cd /scratch/vivek22/UKB/geno
## get relevant files:
# http://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=976
wget  -nd  biobank.ctsu.ox.ac.uk/crystal/crystal/docs/UKBioBiLallfreqSNPexclude.dat
# http://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=1967
wget  -nd  biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_imp_mfi.tgz
# http://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=1955
wget  -nd  biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_qc.
# HRC imputation file:
wget http://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

# age (DF:21003)
cd ..
cat ~/projects/def-ldiatc/uk_biobank/ukb8045.r.tab | head -1 | tr $'\t' '\n' | awk 'BEGIN {OFS="\t"}{print $1, NR}' | grep f\.21003
# column 1140 has age information
cat ~/projects/def-ldiatc/uk_biobank/ukb8045.r.tab | cut -d $'\t' -f 1,1140 > age.tsv

# get sample QC file from Marc (ukb_sqc_v2.txt), and clean it using R
# information available on: https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=531
module load gcc/7.3.0 r/3.6.1
R --no-save
		df <- read.table("ukb_sqc_v2.txt", quote="\"", comment.char="")
		for(i in 1:40){
		  names(df)[i + 25] <- paste0("PC",i)
		}

		df <- df[,-c(1,2,4:9,12:18,21,22,25,66:68)]
		names(df)[1:7] <- c("array","sex","inf_sex","het_miss_outlr","sex_aneup","exces_relatvs","white_brit")
		length(which(df$sex != df$inf_sex)) # 378 IDs sex mismatch
		df$sex_discord <- as.numeric(df$sex != df$inf_sex)
		df <- df[,c(1:3,48,4:47)]

		# add ID from fam file:
		IID <- read.table("~/projects/def-ldiatc/uk_biobank/1.fam", quote="\"", comment.char="")
		length(which(IID$V1 != IID$V2)) # 0
		IID <- IID[,1]
		df <- cbind(IID,df)
		sapply(df[,c(5:8)], function(y) summary(as.factor(y)))
		#   sex_discord het_miss_outlr sex_aneup exces_relatvs
		# 0      487999         487409    487725        488189
		# 1             378                968            652               188

		df[5:9] <- lapply(df[5:9], factor)
		sapply(df, function(y) sum(is.na(y))) # no missing values
		saveRDS(df,"UKB_cov.RDS")

		# cleaning: lost 1,999 IDs
		df <- subset(df, sex_discord == 0 & het_miss_outlr == 0 & sex_aneup == 0 & exces_relatvs == 0)
		df <- df[,-c(4:9)]
		df$array <- as.numeric(ifelse(df$array == "UKBL", 1, 0)); names(df)[2] <- "UKBL"
		df$sex <- as.numeric(ifelse(df$sex == "M", 1, 0)); names(df)[3] <- "MALE"
		
		# get age data for clean IDs:
		age <- read.table("age.tsv", quote="\"", comment.char="", header = T)
		names(age) <- c("IID","AGE") 
		df <- merge(age, df, by = "IID", all.y = T)
				
		saveRDS(df,"clean_UKB_cov.RDS")
		write.table(df, file = "clean_UKB_cov.tab", append = F, quote = F, sep = "\t",
					eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
		
		# make clean sample file
		sam <- read.table("~/projects/def-ldiatc/uk_biobank/1.fam", quote="\"", comment.char="")
		sam <- sam[sam$V1 %in% df$IID, ]
		write.table(sam, file = "clean.sample", append = F, quote = F, sep = "\t",
					eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)
		q()

# make FM pheno:
cd /scratch/vivek22/FM_UKB/pheno
cat ~/projects/def-ldiatc/uk_biobank/ukb8045.r.tab | cut -d $'\t' -f 1,576-604,1585-1964,1993-2427 > fm.tsv
awk -f ~/scripts/awk_ukb_pheno fm.tsv > fm.melt
awk '{if ($5 == 1542){print $1}}' fm.melt > fm_ids
cat fm.melt | grep M797 | awk '{print $1}' >> fm_ids
uniq fm_ids > fm_id ; rm fm_ids

cat ~/projects/def-ldiatc/uk_biobank/ukb8045.r.tab | cut -d $'\t' -f 1,522-528 > pain.tsv
awk -f ~/scripts/awk_ukb_pheno pain.tsv > pain.melt
awk '{if ($5 == -7){print $1}}' pain.melt > cntrl_id

# compare first column of fm and cntrl IDs
awk 'NR==FNR{A[$1];next}$1 in A' fm_id cntrl_id | wc -l
# looks like 102 of our controls are FM patients, remove them
sort fm_id > fm; sort cntrl_id > cntrl; rm fm_id cntrl_id
wc -l cntrl fm
comm -23 cntrl fm > f_cntrl
wc -l cntrl f_cntrl fm
mv f_cntrl cntrl

# add FM info as a column:
awk '{print $0, "1"}' fm > fm_id
awk '{print $0, "0"}' cntrl > cntrl_id
cat fm_id cntrl_id > pheno_id

# add covariates
# module load gcc/7.3.0 r/3.6.1
R --no-save
		ph <- read.table("pheno_id")
		names(ph) <- c("IID","FM")
		cov <- readRDS("../UKB/clean_UKB_cov.RDS")
		df <- merge(ph, cov, by = "IID", all.x = T)
		summary(as.factor(df$FM))
		#      0      1
		# 197050   1826
		saveRDS(df,"FM_pheno_cov.RDS")
		write.table(df, file = "FM_pheno_cov.tab", append = F, quote = F, sep = "\t",
					eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
		q()

# export ukb covariates and FM data to Projects (backed -up)
mkdir /project/def-ldiatc/vivek22/UKB/
cp FM_pheno_cov.tab /project/def-ldiatc/vivek22/UKB/
cp ../UKB/clean_UKB_cov.tab /project/def-ldiatc/vivek22/UKB/
cp ../UKB/UKB_cov.RDS /project/def-ldiatc/vivek22/UKB/

#### Geno QC ####

# grossly inspired from: https://github.com/Nealelab/UK_Biobank_GWAS#imputed-v3-variant-qc
# Primary genotype QC parameters for inclusion in GWAS
# Autosomal SNPs
#### SNPs present in HRC imputation file
# pHWE > 1e-10 in QC positive samples
# call rate > 0.95 in QC positive samples
# UKBB INFO score > 0.8
# AF > 0.001 and AF < 0.999 (AF = alternate allele frequency in QC positive samples)

cd /scratch/vivek22/UKB/geno
tar zxvf ukb_imp_mfi.tgz
for i in {1..22}; do
	awk -v chr=$i 'BEGIN {FS="\t"; OFS="\t"} {print chr,$0}' "ukb_mfi_chr${i}_v3.txt" >> ukb_mfi_v3.tsv
done

# make MAF >= 0.001, info >= 0.8 SNP file:
cat ukb_mfi_v3.tsv | awk '{if($7 >= 0.001 && $7 < 0.5 && $9 >= 0.8){print}}' > qc1
# awk '{print $1":"$4}' < qc1 > qc2
# zcat HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz | awk '{print $1":"$2}' > hrc
# comm -12 <(sort qc2) <(sort hrc) > snps #13512648 snps

# steps here are to make chromosome specific filters:
# for i in {1..22}; do
# grep "^$i:" snps > $i.snp
# done


mkdir qced_snps
mv *snp* ./qced_snps
mkdir impute_info
mv *mfi* ./impute_info/

### FM GENO prep ###
cd /scratch/vivek22/FM_UKB
mkdir pheno
mv * pheno
mkdir geno
mkdir /scratch/vivek22/FM_UKB/geno/qc_snp
cd /scratch/vivek22/FM_UKB/geno/qc_snp
for i in {1..22}; do
awk -v j="$i" '{if($1==j) print $3}' /scratch/vivek22/UKB/geno/qc1 > $i.snp
done

# make QCed plink files from bgen files:
for i in {1..22}; do
cat -> geno_qc_${i} << EOF
#!/bin/bash
#SBATCH --account=def-ldiatc
#SBATCH --mail-user=vivek.verma@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=4G
#SBATCH --time=60:00:00

/home/vivek22/qctool/build/release/qctool_v2.0.7 \
-g ~/projects/def-ldiatc/uk_biobank/${i}.bgen \
-s ~/projects/def-ldiatc/uk_biobank/${i}.sample \
-incl-samples /scratch/vivek22/FM_UKB/pheno/pheno_id \
-incl-rsids ${i}.snp \
-og chr_${i} \
-ofiletype binary_ped \
-threads 40
EOF
done

for i in {1..22}; do
sbatch geno_qc_${i}
done


# ***************************** SAIGE STEP1 ********************************* #
mkdir /scratch/vivek22/FM_UKB/saige
cd /scratch/vivek22/FM_UKB/saige

for i in {1..22}; do
cat -> s1_${i} << EOF
#!/bin/bash
#SBATCH --account=def-ldiatc
#SBATCH --mail-user=vivek.verma@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=4G
#SBATCH --time=60:00:00

module load gcc/7.3.0 r/3.6.1
/home/vivek22/R/x86_64-pc-linux-gnu-library/3.6/SAIGE/extdata/step1_fitNULLGLMM1.R \
        --plinkFile=../geno/qc_snp/chr_${i} \
        --phenoFile=../pheno/FM_pheno_cov.tab \
        --phenoCol=FM \
        --covarColList=AGE,UKBL,MALE,PC1,PC2,PC3,PC4,PC5,\
		PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,\
		PC16,PC17,PC18,PC19,PC20,PC21,PC22,PC23,PC24,PC25,PC26,PC27,PC28,PC29,\
		PC30,PC31,PC32,PC33,PC34,PC35,PC36,PC37,PC38,PC39,PC40 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=out1_${i} \
        --nThreads=40 \
        --LOCO=FALSE
EOF
done

for i in {1..22}; do
sbatch s1_${i}
done
########### make vcf and vcf index files ##################:
mkdir /scratch/vivek22/FM_UKB/geno/vcf
cd /scratch/vivek22/FM_UKB/geno/vcf

for i in {1..22}; do
cat -> vcf_${i} << EOF
#!/bin/bash
#SBATCH --account=def-ldiatc
#SBATCH --mail-user=vivek.verma@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:60:00
module load plink/2.00-10252019-avx2
plink2 --bfile /scratch/vivek22/FM_UKB/geno/qc_snp/chr_${i} --recode vcf --out ${i}
EOF
done
for i in {1..22}; do
sbatch vcf_${i}
done

for i in {1..22}; do
cat -> vcf_gz_${i} << EOF
#!/bin/bash
#SBATCH --account=def-ldiatc
#SBATCH --mail-user=vivek.verma@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=40
#SBATCH --ntasks=40
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00

module load bcftools/1.10.2
module load mugqic/tabix
bcftools view -Oz --no-version --threads 40 ${i}.vcf > ${i}.vcf.gz
tabix -p vcf ${i}.vcf.gz
EOF
done

for i in {1..12}; do
sbatch vcf_gz_${i}
done

# make ID file from fam
# confirm if all fam has same number of IDs:
wc -l ../qc_snp/chr_*.fam
cat ../qc_snp/chr_1.fam | awk '{print $2}' > id

# ***************************** SAIGE STEP2 ********************************* #

cd /scratch/vivek22/FM_UKB/saige
for i in {1..22}; do
cat -> s2_${i} << EOF
#!/bin/bash
#SBATCH --account=def-ldiatc
#SBATCH --mail-user=vivek.verma@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=40
#SBATCH --ntasks=40
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
module load gcc/7.3.0 r/3.6.1
/home/vivek22/R/x86_64-pc-linux-gnu-library/3.6/SAIGE/extdata/step2_SPAtests.R \
        --vcfFile=/scratch/vivek22/FM_UKB/geno/vcf/${i}.vcf.gz \ 
        --vcfFileIndex=/scratch/vivek22/FM_UKB/geno/vcf/${i}.vcf.gz.tbi \
        --chrom=${i} \
        --vcfField=GT \
        --sampleFile=/scratch/vivek22/FM_UKB/geno/vcf/id \
        --GMMATmodelFile=/scratch/vivek22/FM_UKB/saige/out1_${i}.rda \
        --varianceRatioFile=/scratch/vivek22/FM_UKB/saige/out1_${i}.varianceRatio.txt \
        --SAIGEOutputFile=/scratch/vivek22/FM_UKB/saige/out1_${i}_30markers.SAIGE.results.txt 
EOF
done
for i in {1..22}; do
sbatch s2_${i}
done

# confirm step2 worked:
wc -l out1_*_30markers.SAIGE.results.txt

# concat:
cat out1_*_30markers.SAIGE.results.txt > saige_fm_ukb_results.txt

# column 8 and 15 have only 1s (imputationInfo, Is.SPA.converge), skip them:
awk '{$8=$15=""; print $0}' < saige_fm_ukb_results.txt | awk '$1 == ($1+0)' > fm_gwas_res.txt
sort -n -k 1,1 -k 2,2 < fm_gwas_res.txt > gwas_res.txt

#module load gcc/7.3.0 r/3.6.1
R --no-save
        df <- read.table("gwas_res.txt", header = F, col.names = c("CHR", "POS", "SNPID", "A1", "A2", "AC_A2", 
        "AF_A2", "N", "BETA", "SE", "Tstat", "p", "p_noSPA", "varT", "varTstar"), stringsAsFactors = F)
        str(df)
        plot(df$p,df$p_noSPA)
        # install.packages("qqman")
        library("qqman")
        pdf(file = "manhattan.pdf", width = 8, height = 8)
        qq(df$p_lam_corrected)
        dev.off()
        pdf(file = "manhattan.pdf", width = 10, height = 5)
        manhattan(df, chr = "CHR", bp = "POS", p = "p", snp = "SNPID",
        col = c("blue4", "orange3"), chrlabs = c(1:22),
        suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
        logp = TRUE, annotatePval = -log10(1e-05),annotateTop = TRUE)
        dev.off()

#manhattan(subset(gwasResults, CHR == 16), highlight = snpsOfInterest, xlim = c(200, 
#    500), main = "Chr 3")

# back up original results:
mkdir /scratch/vivek22/FM_UKB/saige/all_res
cp /scratch/vivek22/FM_UKB/saige/out1_*_30markers.SAIGE.results.txt /scratch/vivek22/FM_UKB/saige/all_res


# ***************   GWAS ONLY FEMALES  *****************
mkdir /scratch/vivek22/FM_UKB/saige/f_only
cd /scratch/vivek22/FM_UKB/saige/f_only
# prepare sample file
awk -F "\t" '{ if ($5 == 0) {print $1}}' < /scratch/vivek22/FM_UKB/pheno/FM_pheno_cov.tab > females
# prepare f_only_FM_pheno_cov.tab
cd /scratch/vivek22/FM_UKB/saige
module load gcc/7.3.0 r/3.6.1
R --no-




for i in {1..22}; do
cat -> s1_${i} << EOF
#!/bin/bash
#SBATCH --account=def-ldiatc
#SBATCH --mail-user=vivek.verma@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=4G
#SBATCH --time=60:00:00

module load gcc/7.3.0 r/3.6.1
/home/vivek22/R/x86_64-pc-linux-gnu-library/3.6/SAIGE/extdata/step1_fitNULLGLMM1.R \
        --plinkFile=../geno/qc_snp/chr_${i} \
        --phenoFile=../pheno/FM_pheno_cov.tab \
        --phenoCol=FM \
        --covarColList=AGE,UKBL,MALE,PC1,PC2,PC3,PC4,PC5,\
		PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,\
		PC16,PC17,PC18,PC19,PC20,PC21,PC22,PC23,PC24,PC25,PC26,PC27,PC28,PC29,\
		PC30,PC31,PC32,PC33,PC34,PC35,PC36,PC37,PC38,PC39,PC40 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=out1_${i} \
        --nThreads=40 \
        --LOCO=FALSE
EOF
done

for i in {1..22}; do
sbatch s1_${i}
done
########### make vcf and vcf index files ##################:
mkdir /scratch/vivek22/FM_UKB/geno/vcf
cd /scratch/vivek22/FM_UKB/geno/vcf

for i in {1..22}; do
cat -> vcf_${i} << EOF
#!/bin/bash
#SBATCH --account=def-ldiatc
#SBATCH --mail-user=vivek.verma@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:60:00
module load plink/2.00-10252019-avx2
plink2 --bfile /scratch/vivek22/FM_UKB/geno/qc_snp/chr_${i} --recode vcf --out ${i}
EOF
done
for i in {1..22}; do
sbatch vcf_${i}
done

for i in {1..22}; do
cat -> vcf_gz_${i} << EOF
#!/bin/bash
#SBATCH --account=def-ldiatc
#SBATCH --mail-user=vivek.verma@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=40
#SBATCH --ntasks=40
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00

module load bcftools/1.10.2
module load mugqic/tabix
bcftools view -Oz --no-version --threads 40 ${i}.vcf > ${i}.vcf.gz
tabix -p vcf ${i}.vcf.gz
EOF
done

for i in {1..12}; do
sbatch vcf_gz_${i}
done






# repeat saige step 2:
for i in {1..22}; do
cat -> s2f_${i} << EOF
#!/bin/bash
#SBATCH --account=def-ldiatc
#SBATCH --mail-user=vivek.verma@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=40
#SBATCH --ntasks=40
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
module load gcc/7.3.0 r/3.6.1
/home/vivek22/R/x86_64-pc-linux-gnu-library/3.6/SAIGE/extdata/step2_SPAtests.R \
        --vcfFile=/scratch/vivek22/FM_UKB/geno/vcf/${i}.vcf.gz \ 
        --vcfFileIndex=/scratch/vivek22/FM_UKB/geno/vcf/${i}.vcf.gz.tbi \
        --chrom=${i} \
        --vcfField=GT \
        --sampleFile=/scratch/vivek22/FM_UKB/saige/f_only/females \
        --GMMATmodelFile=/scratch/vivek22/FM_UKB/saige/out1_${i}.rda \
        --varianceRatioFile=/scratch/vivek22/FM_UKB/saige/out1_${i}.varianceRatio.txt \
        --SAIGEOutputFile=/scratch/vivek22/FM_UKB/saige/f_only/out1_${i}_30markers.SAIGE.results.txt 
EOF
done
for i in {1..22}; do
sbatch s2f_${i}
done

# confirm step2 worked:
wc -l out1_*_30markers.SAIGE.results.txt

# concat:
cat out1_*_30markers.SAIGE.results.txt > saige_fm_ukb_results.txt

# column 8 and 15 have only 1s (imputationInfo, Is.SPA.converge), skip them:
awk '{$8=$15=""; print $0}' < saige_fm_ukb_results.txt | awk '$1 == ($1+0)' > fm_gwas_res.txt
sort -n -k 1,1 -k 2,2 < fm_gwas_res.txt > gwas_res.txt

# local lab server:
mkdir /scratch/vverma3/fm_ukb/magma
R --no-save
        df <- read.table("gwas_res.txt", header = F, col.names = c("CHR", "POS", "SNPID", "A1", "A2", "AC_A2", 
        "AF_A2", "N", "BETA", "SE", "Tstat", "p", "p_noSPA", "varT", "varTstar"), stringsAsFactors = F)
        # fix p values: # https://www.biostars.org/p/43328/
        Z <- qnorm( df[,"p"] / 2 )
        lambda <- median(( Z^2 )/0.454)
        lambda
        # https://www.biostars.org/p/92190/
        min( df[,"p"] )
        chi2 <- qchisq( df[,"p"], lower.tail=F, df=1 )
        df$chi2 <- chi2 / lambda
        df$p_lam_corrected <- pchisq( df$chi2, df=1, lower.tail=F )
        min( df[,"p_lam_corrected"] )
        # check if lambda is 1
        # save the results
        write.table(df, file = "/scratch/vverma3/fm_ukb/corrected_fm_ukb_gwas_res.txt", append = F, quote = F, sep = "\t",
			eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
        saveRDS(df,"/scratch/vverma3/fm_ukb/corrected_fm_ukb_gwas_res.RDS")

        # qqplot: https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
        library(lattice)
        qqunif.plot<-function(pvalues, 
            should.thin=T, thin.obs.places=2, thin.exp.places=2, 
            xlab=expression(paste("Expected (",-log[10], " p-value)")),
            ylab=expression(paste("Observed (",-log[10], " p-value)")), 
            draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
            already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
            par.settings=list(superpose.symbol=list(pch=pch)), ...) {
            #error checking
            if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
            if(!(class(pvalues)=="numeric" || 
                (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
                stop("pvalue vector is not numeric, can't draw plot")
            if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
            if (already.transformed==FALSE) {
                if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
            } else {
                if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
            }
            
            
            grp<-NULL
            n<-1
            exp.x<-c()
            if(is.list(pvalues)) {
                nn<-sapply(pvalues, length)
                rs<-cumsum(nn)
                re<-rs-nn+1
                n<-min(nn)
                if (!is.null(names(pvalues))) {
                    grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
                    names(pvalues)<-NULL
                } else {
                    grp=factor(rep(1:length(pvalues), nn))
                }
                pvo<-pvalues
                pvalues<-numeric(sum(nn))
                exp.x<-numeric(sum(nn))
                for(i in 1:length(pvo)) {
                    if (!already.transformed) {
                        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
                        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
                    } else {
                        pvalues[rs[i]:re[i]] <- pvo[[i]]
                        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
                    }
                }
            } else {
                n <- length(pvalues)+1
                if (!already.transformed) {
                    exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
                    pvalues <- -log10(pvalues)
                } else {
                    exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
                }
            }


            #this is a helper function to draw the confidence interval
            panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
                require(grid)
                conf.points = min(conf.points, n-1);
                mpts<-matrix(nrow=conf.points*2, ncol=2)
                    for(i in seq(from=1, to=conf.points)) {
                            mpts[i,1]<- -log10((i-.5)/n)
                            mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
                            mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
                            mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
                    }
                    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
                }

            #reduce number of points to plot
            if (should.thin==T) {
                if (!is.null(grp)) {
                    thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                        exp.x = round(exp.x, thin.exp.places),
                        grp=grp))
                    grp = thin$grp
                } else {
                    thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                        exp.x = round(exp.x, thin.exp.places)))
                }
                pvalues <- thin$pvalues
                exp.x <- thin$exp.x
            }
            gc()
            
            prepanel.qqunif= function(x,y,...) {
                A = list()
                A$xlim = range(x, y)*1.02
                A$xlim[1]=0
                A$ylim = A$xlim
                return(A)
            }

            #draw the plot
            xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
                prepanel=prepanel, scales=list(axs="i"), pch=pch,
                panel = function(x, y, ...) {
                    if (draw.conf) {
                        panel.qqconf(n, conf.points=conf.points, 
                            conf.col=conf.col, conf.alpha=conf.alpha)
                    };
                    panel.xyplot(x,y, ...);
                    panel.abline(0,1);
                }, par.settings=par.settings, ...
            )
        }
        pdf(file = "/scratch/vverma3/fm_ukb/qq1.pdf", width = 8, height = 8)
        qqunif.plot(df$p)
        dev.off()
        # extract SNP and P for MAGMA:
        df1 <- df[,c(3,17)]
        names(df1) <- c("SNP","P")
        write.table(df1, file = "/scratch/vverma3/fm_ukb/magma/SNP_P.txt", append = F, quote = F, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)

        # extract significant SNPs:
        sig_snps <- as.list(df$SNPID[df$p_lam_corrected < 5e-08])

        # Manhattan plot: https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_Manhattan_Plots_in_R
        library(lattice)
        manhattan.plot<-function(chr, pos, pvalue, 
            sig.level=NA, annotate=NULL, ann.default=list(),
            should.thin=T, thin.pos.places=2, thin.logp.places=2, 
            xlab="Chromosome", ylab=expression(-log[10](p-value)),
            col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {

            if (length(chr)==0) stop("chromosome vector is empty")
            if (length(pos)==0) stop("position vector is empty")
            if (length(pvalue)==0) stop("pvalue vector is empty")

            #make sure we have an ordered factor
            if(!is.ordered(chr)) {
                chr <- ordered(chr)
            } else {
                chr <- chr[,drop=T]
            }

            #make sure positions are in kbp
            if (any(pos>1e6)) pos<-pos/1e6;

            #calculate absolute genomic position
            #from relative chromosomal positions
            posmin <- tapply(pos,chr, min);
            posmax <- tapply(pos,chr, max);
            posshift <- head(c(0,cumsum(posmax)),-1);
            names(posshift) <- levels(chr)
            genpos <- pos + posshift[chr];
            getGenPos<-function(cchr, cpos) {
                p<-posshift[as.character(cchr)]+cpos
                return(p)
            }

            #parse annotations
            grp <- NULL
            ann.settings <- list()
            label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
                col=NULL, fontface=NULL, fontsize=NULL, show=F)
            parse.label<-function(rawval, groupname) {
                r<-list(text=groupname)
                if(is.logical(rawval)) {
                    if(!rawval) {r$show <- F}
                } else if (is.character(rawval) || is.expression(rawval)) {
                    if(nchar(rawval)>=1) {
                        r$text <- rawval
                    }
                } else if (is.list(rawval)) {
                    r <- modifyList(r, rawval)
                }
                return(r)
            }

            if(!is.null(annotate)) {
                if (is.list(annotate)) {
                    grp <- annotate[[1]]
                } else {
                    grp <- annotate
                } 
                if (!is.factor(grp)) {
                    grp <- factor(grp)
                }
            } else {
                grp <- factor(rep(1, times=length(pvalue)))
            }
        
            ann.settings<-vector("list", length(levels(grp)))
            ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)

            if (length(ann.settings)>1) { 
                lcols<-trellis.par.get("superpose.symbol")$col 
                lfills<-trellis.par.get("superpose.symbol")$fill
                for(i in 2:length(levels(grp))) {
                    ann.settings[[i]]<-list(pch=pch, 
                        col=lcols[(i-2) %% length(lcols) +1 ], 
                        fill=lfills[(i-2) %% length(lfills) +1 ], 
                        cex=cex, label=label.default);
                    ann.settings[[i]]$label$show <- T
                }
                names(ann.settings)<-levels(grp)
            }
            for(i in 1:length(ann.settings)) {
                if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
                ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                    parse.label(ann.settings[[i]]$label, levels(grp)[i]))
            }
            if(is.list(annotate) && length(annotate)>1) {
                user.cols <- 2:length(annotate)
                ann.cols <- c()
                if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
                    ann.cols<-match(names(annotate)[-1], names(ann.settings))
                } else {
                    ann.cols<-user.cols-1
                }
                for(i in seq_along(user.cols)) {
                    if(!is.null(annotate[[user.cols[i]]]$label)) {
                        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                            levels(grp)[ann.cols[i]])
                    }
                    ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                        annotate[[user.cols[i]]])
                }
            }
            rm(annotate)

            #reduce number of points plotted
            if(should.thin) {
                thinned <- unique(data.frame(
                    logp=round(-log10(pvalue),thin.logp.places), 
                    pos=round(genpos,thin.pos.places), 
                    chr=chr,
                    grp=grp)
                )
                logp <- thinned$logp
                genpos <- thinned$pos
                chr <- thinned$chr
                grp <- thinned$grp
                rm(thinned)
            } else {
                logp <- -log10(pvalue)
            }
            rm(pos, pvalue)
            gc()

            #custom axis to print chromosome names
            axis.chr <- function(side,...) {
                if(side=="bottom") {
                    panel.axis(side=side, outside=T,
                        at=((posmax+posmin)/2+posshift),
                        labels=levels(chr), 
                        ticks=F, rot=0,
                        check.overlap=F
                    )
                } else if (side=="top" || side=="right") {
                    panel.axis(side=side, draw.labels=F, ticks=F);
                }
                else {
                    axis.default(side=side,...);
                }
            }

            #make sure the y-lim covers the range (plus a bit more to look nice)
            prepanel.chr<-function(x,y,...) { 
                A<-list();
                maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
                A$ylim=c(0,maxy);
                A;
            }

            xyplot(logp~genpos, chr=chr, groups=grp,
                axis=axis.chr, ann.settings=ann.settings, 
                prepanel=prepanel.chr, scales=list(axs="i"),
                panel=function(x, y, ..., getgenpos) {
                    if(!is.na(sig.level)) {
                        #add significance line (if requested)
                        panel.abline(h=-log10(sig.level), lty=2);
                    }
                    panel.superpose(x, y, ..., getgenpos=getgenpos);
                    if(!is.null(panel.extra)) {
                        panel.extra(x,y, getgenpos, ...)
                    }
                },
                panel.groups = function(x,y,..., subscripts, group.number) {
                    A<-list(...)
                    #allow for different annotation settings
                    gs <- ann.settings[[group.number]]
                    A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
                    A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
                    A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
                    A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
                    A$x <- x
                    A$y <- y
                    do.call("panel.xyplot", A)
                    #draw labels (if requested)
                    if(gs$label$show) {
                        gt<-gs$label
                        names(gt)[which(names(gt)=="text")]<-"labels"
                        gt$show<-NULL
                        if(is.character(gt$x) | is.character(gt$y)) {
                            peak = which.max(y)
                            center = mean(range(x))
                            if (is.character(gt$x)) {
                                if(gt$x=="peak") {gt$x<-x[peak]}
                                if(gt$x=="center") {gt$x<-center}
                            }
                            if (is.character(gt$y)) {
                                if(gt$y=="peak") {gt$y<-y[peak]}
                            }
                        }
                        if(is.list(gt$x)) {
                            gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
                        }
                        do.call("panel.text", gt)
                    }
                },
                xlab=xlab, ylab=ylab, 
                panel.extra=panel.extra, getgenpos=getGenPos, ...
            );
        }
        annotateSNPRegions<-function(snps, chr, pos, pvalue, snplist,
            kbaway=0, maxpvalue=1, labels=c(), col=c(), pch=c()) {

            stopifnot(all(length(snps)==length(chr), length(chr)==length(pos),
                length(pos)==length(pvalue)))
            if (length(snplist)==0) stop("snplist vector is empty")

            if(any(pos>1e6)) kbaway<-kbaway*1000

            ann<-rep(0, length(snps))
            for(i in seq_along(snplist)) {
                si<-which(snps==snplist[i])
                ci<-chr[si]
                pi<-pos[si]
                ann[chr==ci & pos >= pi-kbaway & pos <= pi+kbaway & pvalue<=maxpvalue]<-i
            }
            ann<-list(factor(ann, levels=0:length(snplist), labels=c("", snplist)))
            if(length(col)>0 || length(pch)>0 || length(labels)>0) {
                for(i in seq_along(snplist)) {
                    ann[[ snplist[i] ]] = list()
                    if(length(col)>0) { 
                        ann[[ snplist[i] ]]$col = col[ (i-1) %% length(col)+1 ]
                    }
                    if(length(pch)>0) {
                        ann[[ snplist[i] ]]$pch = pch[ (i-1) %% length(pch)+1 ]	
                    }
                                if(length(labels)>0) {
                                        ann[[ snplist[i] ]]$label = labels[ (i-1) %% length(labels)+1 ]
                                }
                }
            }
            return(ann)
        }
        # get significant snps:
        sig <- df[df$p_lam_corrected <= 5e-8,]
        sig <- sig[order(sig$p_lam_corrected),]
        ann<-annotateSNPRegions(df$SNPID, df$CHR, df$POS, df$p_lam_corrected,
            snplist=c("rs12447745","rs11132324","rs1616352","rs559527636", "rs181664805"),col="red", kbaway=50)
        #name=c("USB1","PALLD","MYL12A","PDIA5","intergenic")
                png("mh.png", width=950, height=400)
                print(manhattan.plot(df$CHR, df$POS, df$p_lam_corrected,
                    annotate=ann, ann.default=list(label=list(offset=2)),
                    sig.level=5e-8,
                    key=list(background="white", border=T, padding.text=2, 
                    corner=c(.95, .95), text=list(lab=c("Known","Novel")), 
                    points=list(col=c("red","green"), pch=20))
                    )
                )
        dev.off()

# >>>>    MAGMA    >>>>>>>>>>>>>>>>>>>>>>>>>>> 
# local lab server:
cd /scratch/vverma3/fm_ukb/magma
awk 'min=="" || $2 < min {min=$2} END {print min}' SNP_P.txt

# run gene-based analyses
base=/home/mparis15/GENOME/magma_v1.07b
${base}/magma \
  --bfile      ${base}/g1000_eur \
  --gene-annot ${base}/NCBI37.3.genes.annot \
  --pval SNP_P.txt \
  N=200000 \
  --out magma1


#--------------------------------------------------
# run pathway-based analyses
#  needs the results of gene-based
#
for theGO in BP CC MF; do

# run pathways
${base}/magma \
  --gene-results magma1.genes.raw \
  --set-annot ${base}/human_GO${theGO}_dec2019.set \
  --out magma1_GO${theGO}

# add pathway desc
# add FDR
# pathways with at most 100 genes
grep ^"#" magma1_GO${theGO}.gsa.out \
  > magma1_GO${theGO}_wFDR.gsa.out
echo "VARIABLE TYPE NGENES BETA BETA_STD SE P FDR DESC" \
  | tr ' ' $'\t' \
  >> magma1_GO${theGO}_wFDR.gsa.out
LANG=en_EN join -t $'\t' -1 1 -2 1 \
<( \
cat magma1_GO${theGO}.gsa.out \
  | grep -v ^"#" \
  | tail -n +2 \
  | tr -s ' ' $'\t' \
  | awk '( $3 <= 100 )' \
  | /home/mparis15/MISC_UTILS/add_FDR.exe -t $'\t' -k 7 \
  | LANG=en_EN sort -t $'\t' -k 1,1 \
 ) \
/scratch/mparis15/PATHWAY_2019dec/pathway_desc.txt \
  | sort -t $'\t' -g -k 7,7 \
  >> magma1_GO${theGO}_wFDR.gsa.out
done

# show top pathways
head magma1_GO??_wFDR.gsa.out

# refining GWAS results for MAGMA:
cd /scratch/vverma3/fm_ukb/magma
R --no-save
df <- readRDS("../corrected_fm_ukb_gwas_res.RDS")

# identify effect allele based on freq:
df$Eff_A <- ifelse(df$AF_A2 < 0.5, df$A2, df$A1)
df$Alt_A <- ifelse(df$AF_A2 >= 0.5, df$A2, df$A1)
df$MAF <- ifelse(df$AF_A2 < 0.5, df$AF_A2, 1 - df$AF_A2)

df1 <- subset(df, grepl("rs", SNPID))


# drop ambiguous SNPs (n = 559,787)
df <- df[df$MAF < 0.45,]
# drop very rare SNPs: (n = 54,923)
df <- df[df$MAF > 0.001,]
# drop super - insignificant SNPs: (n = 3,177,754)
df <- df[df$p_lam_corrected < 0.8,]
# total participants = 192,958
# drop useless columns:
df <- df[,c(1:3,18:20,9:10,17)] 
names(df)[9] <- "pval"

write.table(df, file = "clean_GWAS_summary.tab", append = F, quote = F, sep = "\t",
			eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)

write.table(df[,-6], file = "GWAS_summary.tab", append = F, quote = F, sep = "\t",
			eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
q()

gzip -9 GWAS_summary.tab > GWAS_summary.gz