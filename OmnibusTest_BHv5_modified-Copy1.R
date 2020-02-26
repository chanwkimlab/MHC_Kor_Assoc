###### Omnibus test code
###### Buhm Han, Xinli Hu
###### v1.00 : 2016/12/21
###### v2.00 : 2017/05/11 <- I fixed errors in conditioning.
###### v3.00 : 2017/08/01 <- *.covar file, i.e. collection (group) as covariates.
###### v4.00 : 2018/08/07 <- Now argument "data.prefix" is deprecated and divided into each new argument for those values.

# (Before)
# (1) file prefix
# (2) output prefix
# (3) phenotype name
# (4) rare threshold
# (5) conditional variables

# args <- commandArgs(TRUE)
# data.prefix=args[1] ## prefix of .fam, .aa, .phe, and .covar files
# outfile=args[2]
# phen.name=args[3] ## phenotype name defined in .phe file
# rare.thres=as.numeric(args[4]) ## frequency cut-off to ignore alleles in test
# condvar=args[5] ## comma separated, without spaces. 


# (After)
# (1) output prefix
# (2) .fam
# (3) .aa
# (4) .phe
# (5) pheno-name (2018. 11. 2.)
# (6) .covar
# (7) covar-name (2018. 11. 2.)
# (8) rare threshold
# (9) conditional variables.
install.packages("lmtest")
if (!require("data.table")) install.packages("data.table")
if (!require("tools")) install.packages("tools")
if (!require("lmtest")) install.packages("lmtest")



args <- commandArgs(TRUE)
print('-----------------------------------Argument Information--------------------------------------------')
print(paste('arguments','--out',args[2],'--fam',args[3],'--aa',args[4],'--phe',args[5],'--phename',args[6],'--covar',args[7],'--covarname',args[8],'--condvar',args[9]))
print(paste('assoc_mode:',args[1]))

assoc_mode=args[1]
outfile=args[2]
data.fam = args[3]
data.aa = args[4]
data.phe = args[5]
#data.phe.name=args[6] ## phenotype name defined in .phe file (must be single item.)
data.covar = args[6]
data.covar.name = args[7] # 
# rare.thres=as.numeric(args[8]) ## frequency cut-off to ignore alleles in test
condvar=args[8] ## comma separated, without spaces. 
# If no conditioning, provide "NA"
# Put "_All" for all search results: "AA_DRB1_All,AA_DQA1_All,AA_DQB1_All"  




#print(args)
########## <(1) Preprocessing given arguments.> ##########


# condvar trimming
condvars=unlist(strsplit(condvar,split=','))
if(condvars[1]=="NA") { condvars=c() }

# input data-file name trimming
famfile=paste0(data.fam)
phasedfile=paste0(data.aa) ## Amino acid called file
covfile=paste0(data.covar)
phenfile=paste0(data.phe)

# print(famfile)
# print(phasedfile)
# print(covfile)
# print(phenfile)



########## <(2) Loading Data> ##########

## 1. READ FAM, PHASED, COVAR, PHE
# (1) fam
print('----------------------------------------.fam load------------------------------------------')
fam <- as.matrix(read.table(famfile))
n.fam=dim(fam)[1]
print(paste('Loaded',famfile))
print(paste('individuals:',dim(fam)[1]))

fam=rep(fam, each=2) ## make x2 for haploids



# (2) covar
print('---------------------------------Covariate load----------------------------------------------')

if (file.exists(covfile)) {
  
    # covar <- as.matrix(read.table(covfile, header=T))[,-(1:2)]
    covar <- read.table(covfile, header=T)

    
    ### covar_name processing (2018. 11. 2.)
    if (data.covar.name !="NA"){
      
      covar_target = strsplit(data.covar.name, ",(\\s+)?")[[1]]

      # Subsetting columns
      #print(colnames(covar))
      covar = as.data.frame(covar[,covar_target])
    }
    else{
        covar = as.data.frame(covar[,-(1:2)])
    }
    print(paste('Loaded',covfile))
    print(paste('individual ',dim(covar)[1]))
    print(paste('covar count ',dim(covar)[2]))
    ### Duplicating for .bgl.phased format
    covar=apply(covar, 2, rep, each=2) ## make x2 for haploids
    
    ### Covar Expression
    # covarexp=" + covar.f[,1] + covar.f[,2]"
    covarexp = ""
    
    for (i in 1:dim(covar)[2]) {
      covarexp = paste0(covarexp, " + covar.f[,", i, "]")
    }
    
    is.live.covar=(apply(is.na(covar), 1, sum)==0)

} else {
    print('covariate does not exist')
    covarexp=""
    covar=NA
    is.live.covar=rep(TRUE, 2*n.fam)
}

# (3) pheno
print('-----------------------------------phenotype load----------------------------------------')
pheno <- as.numeric(as.matrix(read.table(phenfile))[, 3])
#print(colnames(read.table(phenfile)))
#print(pheno)
#mode(pheno)="numeric"
#pheno=pheno-1 # 2,1,0 to 1,0,-1
pheno=rep(pheno, each=2) ## make x2 for haploids
 ## who is alive for analysis

if (assoc_mode=="logistic"){
    is.live.phen= ((pheno == 2) | (pheno == 1))
    pheno=pheno-1
}else if (assoc_mode=="linear"){
    is.live.phen= (pheno != -9)
}else{
    stop("invalid assoc_mode(use logistic or linear instead). check assoc_mode parameter")
}

print(paste('Loaded',phenfile))
print(paste('individuals:',length(pheno)/2))
print(paste('individuals without missing phenotypes:',sum(is.live.phen)/2))

# (4) phased
print('-------------------------------------AA load---------------------------------------------')
#phased.in <- as.matrix(read.table(phasedfile))

#if(file_ext(phasedfile)=='aa'){
#    phased_read=fread(phasedfile,header=FALSE) # faster version
#    #save(phased_read,file=paste0(phasedfile,'.RData'))
#    print(paste0("Generated", phasedfile,'.RData', 'Use this file as .aa file parameter to shorten loading time.'))
#}else if(file_ext(phasedfile)=='Rdata'){
#    phased_read=load(file=paste0(phasedfile,'.RData'))
#}else{
#    stop("unsupported aa file")
#}
#if(!file.exists(paste0(phasedfile,'.RData'))){
#    phased_read=fread(phasedfile,header=FALSE) # faster version
#    print(paste("Generating", phasedfile,'.RData', 'Use this file as .aa file parameter to shorten loading time.'))
#    save(phased_read,file=paste0(phasedfile,'.RData'))
#    print(paste("Generated"))
#}else{
#    print(paste("Loading", phasedfile,'.RData','instead of',phasedfile,'to reduce loading time.'))
#    phased_read=load(file=paste0(phasedfile,'.RData')) 
#    print('Loaded')
#}
phased_read=fread(phasedfile,header=FALSE) # faster version
#phased_read=read.table(phasedfile)
phased.in <- as.matrix(phased_read)
variants=phased.in[-(1:5),2]
phased=t(phased.in[-(1:5),-(1:2)]) ## transpose: rows are individuals, cols are markers
colnames(phased)=variants
print(paste('Loaded',phasedfile))
print(paste('Individuals:',dim(phased)[1]/2,"Variants:",dim(phased)[2]))


########## <(Optional) Condvar> ##########
print('-------------------------------------Condvar load---------------------------------------------')
print(paste('condvars:',condvars))
print(paste('num of condvards',length(condvars)))
## 2. DEFINE HAPLOTYPES TO BE CONDITIONED ON
is.live.cond=rep(TRUE, length(pheno)) # who is alive for analysis
condexp=""
hap.cond=NULL
if (length(condvars)>0) {
	condvariants=c()
	for (thiscon in condvars) {
		if (grepl("AA",thiscon) || grepl("HLA",thiscon)) {
			if (grepl("_All",thiscon)) {
				thiscon <- strsplit(thiscon,"_All")[[1]][1] ## remove "_All" before search
            }
		    thiscon <- variants[grepl(thiscon,variants)] ## search variables
		}
		condvariants <- c(condvariants, thiscon)
	}
    #print(condvariants)
    condmatrix=phased[ ,condvariants, drop=FALSE]
    condexp=" + hap.cond.f"
    hap.cond=apply(condmatrix, 1, paste0, collapse='')
    is.live.cond=(apply(is.na(condmatrix), 1, sum)==0)
    
    print(paste('Matched condvar:',condvariants))
}


########## <(3) Regression> ##########
print('-------------------------------------Regression Process---------------------------------------------')

## 3. BEGIN REGRESSION
aavariants <- variants[which(grepl("AA_",variants))] ## test only amino acids
testvariants <- aavariants
#if (!is.null(aacondvariants)) {
#	testvariants <- aavariants[which(!aavariants %in% aacondvariants)]
#}
results <- matrix(ncol=6, nrow=length(testvariants))

print(paste('testvariants:'))
print(t(as.matrix(testvariants)))
print(paste('# of testvariants:',length(testvariants)))
print(paste("association mode",assoc_mode))
print('------------------------------Iteration for each variant starts------------------------------')
#
for (i in 1:length(testvariants)) {
    
    if (i!=2){
        next
        
    }    
    
#for (i in 210:210) {
	#  if (i%%10 == 0) {
	#	  cat("progress:",100*i/length(testvariants),'%','\n')
  	#}
	#variant <- as.character(testvariants[i])
    variant=testvariants[i]
  	vcol <- which(variants==variant)
  	newaa <-phased[,vcol]
    
    is.live.aa=!is.na(newaa)
    is.live=(is.live.phen & is.live.cond & is.live.aa & is.live.covar) ## who is finally alive
    n.is.live=sum(is.live)
    
    print(paste0("[",i,'/',length(testvariants),"] ",variant))
    print(paste0(sum(is.live)/2,'(',(length(is.live)/2),')',' remains after filtering ',
                'phen-',sum(is.live.phen)/2,'(',length(is.live.phen)/2,') ',
                'cond-',sum(is.live.cond)/2,'(',length(is.live.cond)/2,') ',
                'aa-',sum(is.live.aa)/2,'(',length(is.live.aa)/2,') ',
                'covar-',sum(is.live.covar)/2,'(',length(is.live.covar)/2,') '
                )
         )
    #print(length(is.live.phen))
    #print(length(is.live.cond))
    #print(length(is.live.aa))

    ## KEEP ALIVE PEOPLE
    pheno.f=pheno[is.live]
    if (!is.null(hap.cond)) {
        hap.cond.f=hap.cond[is.live]
    } else {
        hap.cond.f=NULL
    }
    newaa.f=newaa[is.live]
    # covar.f=as.matrix(covar[is.live,])  # as.matrix() is introduced in case only one covariate name is given(To prevent automatically converting to vector not matrix.)
    
    if(covarexp!=''){
    print("covar.f was generated")
    covar.f=as.matrix(covar[is.live,])
    }

    ## DEFINE NEW HAPLOTYPES
  	residues <- unique(newaa.f)
  	#hap.new.f = apply(cbind(hap.cond.f,newaa.f),1,paste0,collapse='',sep='')
  	hap.new.f = apply(cbind(newaa.f),1,paste0,collapse='',sep='')
    save.image(file='yoursession2.RData')
  	#ALTERNATIVE MODEL
  	if (length(unique(residues)) == 1) {
    		results[i,] <- c(variant,"NaN","NaN","NaN","NaN",paste0(residues,collapse=','))
  	}
  	if (length(unique(residues)) > 1) {
        
        if (assoc_mode=="logistic"){
            nullexp <- paste0("glm(pheno.f ~ 1", covarexp, condexp, ", family=binomial(logit))")
            altexp <- paste0("glm(pheno.f ~ 1", covarexp, condexp, "+ hap.new.f, family=binomial(logit),maxit=100)")
        }else if (assoc_mode=="linear"){
            nullexp <- paste0("lm(pheno.f ~ 1", covarexp, condexp, ")")
            altexp <- paste0("lm(pheno.f ~ 1", covarexp, condexp, "+ hap.new.f)")

        }else{
            stop("invalid assoc_mode(use logistic or linear instead). check assoc_mode parameter")
        }
        
        glm.null <- eval(parse(text=nullexp))
        glm.alt <- eval(parse(text=altexp))
        
        lrtest_result=lrtest(glm.null,glm.alt)

        
        nulldeviance <- summary(glm.null)$deviance
        nulldf <- summary(glm.null)$df[1]         
        
		#STATISTICS
    	deviancediff <- lrtest_result$Df[2]
    	dfdiff <- lrtest_result$Df[2]
    	log10pvalue <- pchisq(deviancediff, df=dfdiff, lower.tail=FALSE, log.p=TRUE)/log(10)

    	results[i,] <- c(variant, deviancediff, dfdiff, n.is.live, log10pvalue, paste0(residues,collapse=','))
  	}

}


