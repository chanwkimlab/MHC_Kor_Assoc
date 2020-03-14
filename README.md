# HLA and Disease Association study in Korean Population



### Genotype data
SNP: xx

## MHC region genotype imputation
### China asia pan-a
markers=
N=
### CookHLA
citation
### SNP2HLA
citation

include cohort in associationt only if case in cohort>0 (questionnaire exists)


other github repo

Rules


## Association test

### logistic regression

### Amino acid polymorphism omnibus p-value
### gene assignment
https://www.ncbi.nlm.nih.gov/snp/?term=6%5BChromosome%5D+AND+1500000%3A3000000%5BBase+Position%5D
ftp://ftp.ncbi.nih.gov/snp/latest_release

## Plotting
code


## GCTA h2

## GCTA correlation



# Abstract - Introduction
Population specific LD structure ?? -> 다르다??
imputation method, large sample size, reference panel
newly identified loci -> ?
* Large scale cohort
imputing genotype for large numer of individuals is very computationally burdensome
Power
* Imputation high accuracy
CookHLA, HATK
QC method
* Population
Highly polymorphic
reference panel optimized for Korean. Asian
LD structure especially in HLA region is highly dependent on population.
* Findings 1
identification of novel loci in HLA associated with disease.
multiple independent association signal
* Finding 2
Heritability analysis
PheWAS and genentic relationship analysis identifies pleiotropic loci

# Data resource-cohort
We used as a data resource a large prospective cohort study, the Korean genome and epidemiology study (KoGES) initiated by the Korean government (National Research Institute of Health (NIH), Centers for Disease Control and Prevention and the Ministry of Health and Welfare, Korea) (CITE https://academic.oup.com/ije/article/46/2/e20/2622834 http://www.nih.go.kr/menu.es?mid=a50303010100 http://www.nih.go.kr/menu.es?mid=a40504010000). Among the six cohorts included in the study, which can be categorized into population-based and gene-environment model studies, in our analysis we used three population-based cohort studies, city (N=99159, mean age=, sex=),  rural (N=1890, mean age=, sex=), and Ansan,Ansung(region-based ) (N=7607, mean age=, sex=).  Ethical consent~~~~

# Genotype
Genotype data was obtained using the Korea Biobank array (KCHIP) designed by the Korean National Institute of Health (KCDC). The array is a microarray that contains tagging SNPs that maximizes genomic coverage based on the Korean genomic structure  Particularly, it contains functional SNPs, such as, nonsynonymous, expression
quantitative trait loci (eQTL), disease associated SNPs and also contains about 2,000 pharmacogenomic variants.
Detailed QC criteria regarding data from the array is described elsewhere (CITE). 
For the core MHC region(29-34Mb) and flanking MHC region(28-29, 34-35Mb), we used separate imputation procedures.
## HLA region
Accurately genotyping the HLA region located on chromosome 6p21 is difficult due to the following issues (Densely compacted genes, Highly polymorphic, Population-specific LD structure, Long-range haplotype, Densely compacted genes) (CITE). We adopted several approaches to resolve the underlying issues.
### Construction of population-specific HLA imputation reference panel 
First, we used population specific reference panel to impute HLA gene alleles, amino acid polymorphisms, and  SNPs. Pan-China~~~ (CITE). LD structure of HLA region is highly dependent on population. Whole genome sequencing all the individuals using NGS is highly expensive.

```
The imputation accuracy of the constructed HLA imputation reference panel was empirically evaluated by a cross-validation approach12. We randomly split the panel into two data sets (n=560 for each data set). HLA alleles from one of the data sets were masked and then imputed by using another data set as an imputation reference. The concordance between imputed and genotyped HLA allele dosages was calculated separately for each HLA gene and each allele digit. To relatively compare the imputation accuracy among the different reference panels, we evaluated accuracy in the previously reported HLA imputation reference panel of independent Japanese individuals in the same way (n=908)4,12.
```
### HLA imputation and data preparation?
Second, we utilized (HATK / CookHLA), a collection of several tools that helps performing HLA fine-mapping analysis in a straightforward and accurate way. (CITE) The toolkit contains CookHLA, a highly-accurate HLA imputation software which is compatible with individuals >100,000. Also, the toolkit automatically builds a marker panel including SNP markers and variants of HLA region using amino acid and DNA sequence information distributed by the IPD-IMGT/HLA database, a database that provides the official and most detailed information of the HLA region managed by 'WHO Nomenclature Committee For Factors of the HLA System’.
Through these approach, we could type variants in HLA region for 125673 individuals with high accuracy and low cost. In our analysis we extracted 5292 SNPS from the whole SNPs typed using the KCHIP. Then we imputed variants in HLA region (HLA gene allele, Amino acid polymorphism, SNPs) using HATK/CookHLA.
### QC
For variants in the core HLA region, we phased haplotypes using Beagle (CITE). We obtained dosage of 280 four-digit alleles and 396 amino acid polymorphism of 8 HLA genes (A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1), as well as 11,144 SNPs. For the PheWAS, we applied stringent postimputation QC filtering(MAF ≥0.5%) for the binary markers After QC filtering, 135 four-digit alleles and 361 amino acid polymorphisms of HLA genes, and 10,817 SNP remained.
```
MAF>0.5% (for plink file, just used --maf 0.5. for amino acid polymorphism residue in *.aa file set the residue as NA so as to be excluded in omnibus test `is.live.aa=!is.na(newaa)`)
convert to binary markers(HLA allele, AA residue, SNP)(*.bed/bim/fam), phase the markers(*.bgl.phased), and AA call (*.aa)
<bmarkergenerator, Beagle, AllCC_Get_Phased_AA_Calls.R>
```
## flanking HLA region
To obtain variants in flanking HLA region(28-35Mb), we imputed SNPs and indels in the region using ??? SNPs and 1000G reference genome as a reference panel.
`vcf file. plink vcf->plink. `
### QC
For the PheWAS, we applied stringent postimputation QC filtering(MAF≥1% and  Imputation quality score>=0.8). We did not phaed variants in 1000G since there are only 300 SNVs overlapping each other.
## Merge
For the 300 variants coexists in 1000G and HLA imputation reference panel, we adopted genotype from HLA panel.
## final result
```
* (*.bed/bim/fam) unphased, dosage of KCHIP HLA allele, KCHIP SNP, 1000G SNV
* (*.aa) phased KCHIP HLA allele, KCHIP AA polymorphism, KCHIP SNP
```
# Phenotype
The three cohort studies shares core questionnaire and examination items. The identical questionnaires, physical examinations and clinical investigations were mostly used during the baseline and follow-up phases. The participants were questioned by trained interviewers regarding their socio-demographic status. (CITE)

## Phenotype Curation
How we defined phenotypes using questionnaire data can be found in Supplementary Table. The conversion from questionnaires to phenotypes was conducted using a sofware included in our GitHub repository (https://github.com/ch6845/MHC_Korean_assocation), which automatically merge questionnaires from multiple cohorts merge. and manual syntax analyzer. 
### Detailed criteria
```
Case filtering base. do not include
if merge pheno only exist in a cohort, include them.
exclude medication
include treatment  (more explicit??)
hierachy qudstion do not test under if upper is +
exclude male in analysis of breast cancer, cervical cancer. exclude female in analysis of prostate cancer

define individuals affected by diabetes, hyperlipidemia, hypertension, allergic_disease, colon polyps, or rheumatoid_arthritis as 'unhealthy individual'
If 'unhealthy individuals'is in case-> set the individuals as missing
(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6708789/ After adjusting for age, gender, race, diabetes, hyperlipidemia, hypertension and all significant alleles shown by univariate analysis, multivariate analysis showed no independent association with NAFLD.)
(http://dx.doi.org/10.1038/s41588-018-0336-0 For the controls in disease association studies, we constructed a shared control group by excluding individuals affected by diseases known to have associations in the MHC region. )
for quantitative traits, individuals >mean+3*std or <mean-3*std were set as missing
```

## association
### <basic information>
in first step(without conditional analysis) test KCHIP AA by Omnibus test R code and test KCHIP HLA, KCHIP SNP, 1000G SNP by Plink
use Top 4 PC obtained from 1000G as covariates
### for quantitative traits
plink --assoc --linear --covar
Omnibus R code 
null-`glm(pheno.f ~ 1", covarexp, condexp, ", family=gaussian(identity))`
alternative-`glm(pheno.f ~ 1", covarexp, condexp, "+ hap.new.f, family=gaussian(identity),maxit=100)`
`log10pvalue <- pchisq(deviancediff, df=dfdiff, lower.tail=FALSE, log.p=TRUE)/log(10)`
### for binary traits
plink --assoc --logistic
Omnibus R code 
null-`glm(pheno.f ~ 1", covarexp, condexp, "+ hap.new.f, family=binomial(logit),maxit=100)`
alternative-`glm(pheno.f ~ 1", covarexp, condexp, ", family=binomial(logit))`
`log10pvalue <- pchisq(deviancediff, df=dfdiff, lower.tail=FALSE, log.p=TRUE)/log(10)`
    
    
    
### Conditional analysis (new)    
for testing binary allele -> plink
for testing multiallele allele -> multiallelic module   
    
prepare two sets
1. plink (KCHIP HLA,SNP / 1000G SNP)
2. bgl.phased (KCHIP HLA,AA,SNP) -> filter bgl.phased (KCHIP HLA,AA)
    
1. plink 
    1.1 covariate
        KCHIP HLA,SNP 
    1.2 cond
2. multiallele
    2.1 covariate
        
    2.2 cond
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
## conditional analysis (from second step)
(http://dx.doi.org/10.1038/s41588-018-0336-0
When the top-associated variant itself was the HLA gene polymorphism or the SNV and indel in strong LD with any of the HLA gene polymorphisms (r2 ≥ 0.7), we additionally included all the two-digit, four-digit and six-digit alleles and the amino acid polymorphisms of the corresponding HLA gene as covariates in the regression to robustly condition the associations attributable to the HLA gene, as previously described4,18. Otherwise, the top-associated SNV and indels were additionally included as the associated variants.)
### expand top signal
#### case 1 top signal from KCHIP HLA, KCHIP AA, KCHIP SNP inside HLA gene
use all KCHIP HLA,AA,SNP from the corresponding gene as conditional variant
#### case 2 top signal from KCHIP SNP outside HLA gene, 1000G SNP
if r2 with HLA gene polymorphisms(HLA allele, AA, SNP) is larger than 0.7->use all KCHIP HLA,AA,SNP from the corresponding gene as conditional variant
else use the SNP only as covariate
### distribute conditional variants as covariates or conditional variants to Plink, Omnibus R
#### case 1 KCHIP AA
Omnibus-conditional(condvar parameter)-the AA polymorphism position
Plink-covariate-dosage vectors of residues excluding the one most frequent residue to avoid colinearity
#### case 2 KCHIP HLA or KCHIP SNP
Omnibus-conditional(condvar parameter)-HLA or SNP name (condition the variant using conditional parameter not covariate parameter since it is phased)
Plink-conditional(--condition-list)-HLA or SNP name
#### case 3 1000G SNP
Omnibus-covariate-dosage vectors of 1000G SNP
Plink-conditional(--condition-list)-1000G SNP

# post-analysis
## assign top association signals to Gene (rsid->gene using dBSNP or UCSC gene database)
## Manhattan plot generation (plotting SNP, AA p-value together <->plot in HATK)
## 2D association result plot(row:gene-col:phenotype value:association p-value)
to see pleiotropic effect
(http://dx.doi.org/10.1038/s41588-018-0336-0 Fig4 2D Matrix association result plot generation)
## Haseman–Elston regression (momentum based method) heritability, genetic correlation estimation
<GCTA>
### convert plink file of (HLA gene, AA residue, SNP) to GRM(genetic relationship matrix)
### Haseman–Elston regression


