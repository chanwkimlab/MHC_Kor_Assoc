# Phenome-wide association study of the major histocompatibility complex region in the Korean population 


## Abstract

In this study, we comprehensively analyzed associations between phenotypes and genetic variants in MHC region for Korean population using the Korean Genome and Epidemiology Study (KoGES) data. This data consists of phenotype and genotype information of ~120,000 Korean individuals11. The phenotypes include answers to an extensive set of questionnaires along with data from physical examinations and clinical investigations collected during the baseline and follow-up phases. The individuals were genotyped using the KCHIP (Korea Biobank Array), a microarray that was designed for Korean population12. 

We wanted to use this new single-population dataset to confirm already reported associations and discover novel associations between HLA variants and phenotypes. To this end, we curated 97 phenotypes from questionnaires in KoGES. We imputed genotypes of 8 classical HLA genes, which comprised a total of 119 alleles. We defined markers for the alleles, amino acid residues at 349 amino acid positions in 8 HLA genes, and 67,091 SNVs in the MHC region. We then conducted PheWAS to identify new phenotypes associated with the region. Our analysis confirmed significant associations identified in previous studies, and found 36 novel association signals that were not reported previously. In addition, we analyzed heritability explained by variants in the MHC region and genetic correlations among phenotypes.

## Sample Figures

## PheWAS
![PheWAS](github_images/phewas.png)

## HLA variant Omnibus Test
![Omnibus](github_images/HLA_omnibus.png)

## heritability explained by variants in the MHC region
![h2](github_images/h2_uni_bivar.png)

# GWAS Manhattan Plot
![wbc](github_images/wbc.step_01.merge.manhattan.png)

![height](github_images/height.step_01.merge.manhattan.png)


# Software used

## Genotype munge
- [Vcftools](http://vcftools.sourceforge.net/)
- [Plink](www.cog-genomics.org/plink/2.0/)

## Imputation
- [CookHLA]()

## Association test
- [GAT (Generic Genome-Wide Association Tool)](https://github.com/ch6845/GAT)


## HLA variant
- [HATK(HLA Analysis Toolkit)](https://github.com/WansonChoi/HATK)

# Citation
(Under preparation)