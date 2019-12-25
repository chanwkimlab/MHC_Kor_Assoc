data_path='data/'



data_genotype_path=data_path+'genotype/'
data_phenotype_path=data_path+'phenotype/'
data_codebook_path=data_path+'codebook/'
data_out_pheno_path=data_path+'out_pheno/'
#data_out_pheno_path=data_path+'out_pheno_sumstats/'
data_out_assoc_path=data_path+'out_assoc/'

imputed_bfile_path=data_genotype_path+'HLA_IMPUTED_Result.KCHIP_HLA.MHC'
imputed_bim_path=data_genotype_path+'HLA_IMPUTED_Result.KCHIP_HLA.MHC.bim'
imputed_fam_path=data_genotype_path+'HLA_IMPUTED_Result.KCHIP_HLA.MHC.fam'
imputed_bed_path=data_genotype_path+'HLA_IMPUTED_Result.KCHIP_HLA.MHC.bed'

codebook_AS_path=data_codebook_path+'codebook_AS.csv'
codebook_CT_path=data_codebook_path+'codebook_CT.csv'
codebook_NC_path=data_codebook_path+'codebook_NC.csv'
codebook_custommerge_path=data_codebook_path+'codebook_custommerge.csv'

phenotype_raw_AS_path=data_phenotype_path+'AS1_PHENO_DATA.txt'
phenotype_raw_CT_path=data_phenotype_path+'CT1_PHENO_DATA.txt'
phenotype_raw_NC_path=data_phenotype_path+'NC1_PHENO_DATA.txt'

conversion_manual_path=data_path+'conversion_manual.tsv'

conversion_manual_codebook_path=data_path+'conversion_manual_codebook.csv'

pheno_file_path=data_out_pheno_path+'{}.pheno'
pheno_sumstatsjpg_file_path=data_out_pheno_path+'{}.jpg'

assoc_file_path=data_out_assoc_path+'{}'