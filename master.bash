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
		# 1             378                968           652               188

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
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00

module load bcftools/1.10.2
module load mugqic/tabix
bcftools view -Oz --no-version --threads 40 ${i}.vcf > ${i}.vcf.gz
tabix -p vcf ${i}.vcf.gz
EOF
done

for i in {1..22}; do
sbatch vcf_gz_${i}
done

# make ID file from fam
# confirm if all fam has same number of IDs:
wc -l ../qc_snp/chr_*.fam

cat ../qc_snp/chr_1.fam | awk '{print $2}' > id

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ***************************** SAIGE STEP2 ********************************* #

cd /scratch/vivek22/FM_UKB/saige
module load gcc/7.3.0 r/3.6.1
time \
/home/vivek22/R/x86_64-pc-linux-gnu-library/3.6/SAIGE/extdata/step2_SPAtests.R \
        --vcfFile=/scratch/vivek22/FM_UKB/geno/vcf/22.vcf.gz \
        --vcfFileIndex=/scratch/vivek22/FM_UKB/geno/vcf/22.vcf.gz.tbi \
        --sampleFile=/scratch/vivek22/FM_UKB/geno/vcf/id \
        --GMMATmodelFile=/scratch/vivek22/FM_UKB/saige/out1_22.rda \
        --varianceRatioFile=/scratch/vivek22/FM_UKB/saige/out1_22.varianceRatio.txt \
        --SAIGEOutputFile=/scratch/vivek22/FM_UKB/saige/out1_22_30markers.SAIGE.results.txt \
        --IsOutputAFinCaseCtrl=TRUE \
        --IsOutputNinCaseCtrl=TRUE \
        --IsOutputPvalueNAinGroupTestforBinary=TRUE 



awk '{ if ($14 < 0.1) { print } }' /scratch/vivek22/FM_UKB/saige/out1_22_30markers.SAIGE.results.txt 

    #1
    #2
    #3









