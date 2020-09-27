data_path='data/'

#########Directories############

data_genotype_path=data_path+'genotype/'
data_genotype_1000G=data_genotype_path+'0_1000G/'
data_genotype_unimputed_path=data_genotype_path+'1_KCHIP_unimputed/'
data_genotype_hped_path=data_genotype_path+'2_KCHIP_hped/'
data_genotype_bmarkerphased_path=data_genotype_path+'3_KCHIP_bmarkerphased/'
data_genotype_merge_path=data_genotype_path+'4_merge/'

data_phenotype_path=data_path+'phenotype/'
data_codebook_path=data_path+'codebook/'
data_out_pheno_path=data_path+'out_pheno/'
data_out_assoc_path=data_path+'out_assoc/'
data_out_gcta_path=data_path+'out_gcta/'
data_out_pub_path=data_path+'out_pub/'

#########Files############


PC_path=data_genotype_1000G+'KCHIP_130K_PCA2_annotated2_BHAN.txt'

individual_common_path=data_path+'individual.tsv'

#plink_KCHIP_HLA_AA_SNP_1000G_path=data_genotype_merge_path+'KCHIP_HLA_AA_SNP_1000G'
plink_KCHIP_HLA_AA_SNP_1000G_path=data_genotype_merge_path+'KCHIP_HLA_AA_SNP_1000G'
#plink_KCHIP_HLA_SNP_1000G_path=data_genotype_merge_path+'KCHIP_HLA_SNP_1000G' # for association
plink_KCHIP_SNP_1000G_path=data_genotype_merge_path+'KCHIP_SNP_1000G'
plink_1000G_path=data_genotype_merge_path+'1000G'
#plink_KCHIP_HLA_SNP_path=data_genotype_merge_path+'KCHIP_SNP'
#plink_1000G_path=data_genotype_merge_path+'1000G'

#phased_KCHIP_HLA_AA_path=data_genotype_merge_path+'KCHIP_HLA_AA.bgl.phased'
phased_KCHIP_HLA_AA_SNP_path=data_genotype_merge_path+'KCHIP_HLA_AA_SNP.bgl.phased'


#final_aa_path=data_genotype_merge_path+'KCHIP_HLA.hg18.intersection_HAN.LABELED.NoSameAllele.bMarkers.2field.saveRareAllele.beagle5.1.AGM.bgl.phased.QC.aa'
#final_plink_path=data_genotype_merge_path+'KCHIP_HLA_SNP_1000G_merged'
#final_plink_aa_path=data_genotype_merge_path+'KCHIP_HLA_AA_SNP_1000G_merged'
#final_plink_aa_grm_path=data_genotype_merge_path+'KCHIP_HLA_AA_SNP_1000G_merged.grm'

#################


#bmarker_plink_path=data_genotype_bmarkerphased_path+'KCHIP_HLA.hg18.intersection_HAN.LABELED.NoSameAllele.bMarkers.2field.saveRareAllele'
#bmarker_phased_path=data_genotype_bmarkerphased_path+'KCHIP_HLA.hg18.intersection_HAN.LABELED.NoSameAllele.bMarkers.2field.saveRareAllele.beagle5.1.AGM.bgl.phased'
#bmarker_aa_path=data_genotype_bmarkerphased_path+'KCHIP_HLA.hg18.intersection_HAN.LABELED.NoSameAllele.bMarkers.2field.saveRareAllele.beagle5.1.AGM.bgl.phased.aa'

#merge_plink_path=data_genotype_merge_path+'merge'
#imputed_bfile_path=data_genotype_path+'HLA_IMPUTED_Result.KCHIP_HLA.MHC'
#imputed_bim_path=data_genotype_path+'HLA_IMPUTED_Result.KCHIP_HLA.MHC.bim'
#imputed_fam_path=data_genotype_path+'HLA_IMPUTED_Result.KCHIP_HLA.MHC.fam'
#imputed_bed_path=data_genotype_path+'HLA_IMPUTED_Result.KCHIP_HLA.MHC.bed'

codebook_AS_path=data_codebook_path+'codebook_AS.csv'
codebook_CT_path=data_codebook_path+'codebook_CT.csv'
codebook_NC_path=data_codebook_path+'codebook_NC.csv'
codebook_custommerge_path=data_codebook_path+'codebook_custommerge.csv'

phenotype_raw_AS_path=data_phenotype_path+'AS1_PHENO_DATA.txt'
phenotype_raw_CT_path=data_phenotype_path+'CT1_PHENO_DATA.txt'
phenotype_raw_NC_path=data_phenotype_path+'NC1_PHENO_DATA.txt'

conversion_manual_path=data_path+'conversion_manual.tsv'

conversion_manual_codebook_path=data_path+'conversion_manual_codebook.csv'

pheno_file_path=data_out_pheno_path+'{}.phe'
pheno_all_file_path=data_phenotype_path+'phenotypes_all.tsv'
#pheno_sumstatsjpg_file_path=data_out_pheno_path+'{}.jpg'

#assoc_file_path=data_out_assoc_path+'{}'