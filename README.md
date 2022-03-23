# Phenome-wide association study of the major histocompatibility complex region in the Korean population identifies novel association signals (Kim et al., 2022)


## Abstract

Human leukocyte antigen (HLA) gene variants in the major histocompatibility complex (MHC) region are associated with numerous complex human diseases and quantitative traits. Previous phenome-wide association studies (PheWAS) for this region demonstrated that HLA association patterns to the phenome have both population-specific and population-shared components. We performed MHC PheWAS in the Korean population by analyzing associations between phenotypes and genetic variants in the MHC region using the Korea Biobank Array project data samples from the Korean Genome and Epidemiology Study cohorts. Using this single-population dataset, we curated and analyzed 82 phenotypes for 125 673 Korean individuals after imputing HLA using CookHLA, a recently developed imputation framework. More than one-third of these phenotypes showed significant associations, confirming 56 known associations and discovering 13 novel association signals that were not reported previously. In addition, we analyzed heritability explained by the variants in the MHC region and genetic correlations among phenotypes based on the MHC variants.

## Software

- For genotype munge,
    - [Vcftools](http://vcftools.sourceforge.net/)
    - [Plink](www.cog-genomics.org/plink/2.0/)

- For HLA Imputation,
    - [CookHLA](https://github.com/WansonChoi/CookHLA)

- For HLA variant association test,
    - [GAT (Generic Genome-Wide Association Tool)](https://github.com/ch6845/GAT)

- For HLA variant munge,
    - [HATK(HLA Analysis Toolkit)](https://github.com/WansonChoi/HATK)

## Download summary statistics
- [Korean Biobank Array Project Website](https://www.koreanchip.org/downloads)

## Citation
If you use any part of this code or our data, please cite our
[paper](https://doi.org/10.1093/hmg/ddac016).
```
@article{kim2022phenome,
  title={Phenome-wide association study of the major histocompatibility complex region in the Korean population identifies novel association signals},
  author={Kim, Chanwoo and Kim, Young Jin and Choi, Wanson and Jang, Hye-Mi and Hwang, Mi Yeong and Jung, Sunwoo and Lim, Hyunjoon and Hong, Sang Bin and Yoon, Kyungheon and Kim, Bong-Jo and Park, Hyun-Young and Han, Buhm},
  journal={Human molecular genetics},
  year={2022}
}
```

## Contact
If you have any inquiries, please feel free to contact
- [Chanwoo Kim](https://chanwoo.kim) (Paul G. Allen School of Computer Science & Engineering @ the University of Washington)
