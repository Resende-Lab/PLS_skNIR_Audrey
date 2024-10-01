############################################################################################
###        Prediction of phytoglycogen content in sweet corn using 
###             Single-kernel Near-Infrared spectroscopy 
###
###   Mahon et al. 2024
###   doi: 
### 
############################################################################################
rm(list=ls())


###-----------------------------------
####---- 1. Packages
###-----------------------------------

library(tidyverse)
library(ggplot2)
library(readr)
library(tidyr)
library("pls")
library(Hmisc)
library(squash)
library(chillR)


###-----------------------------------
####---- 2. Functions
####----          
###-----------------------------------

#SNV function
SNV<-function(nir){ 
  nir <- as.matrix(nir)
  nir<-t(nir)
  NIRsnv<-scale(nir, center=TRUE, scale=TRUE)
  NIRsnv<-t(NIRsnv)
  return(NIRsnv)
  }   # return matrix

# RPD function
RPD <- function(measured, RMSEP){
  sd(measured)/RMSEP
}

###-----------------------------------
####---- 3. Traits
###-----------------------------------

#----------- 3.1 Observation
# Obs.: After carefully analyzing the data, we found that every trait has one
# number of components that better fit the data.
# For the list of traits, one component will be selected.
#


#----------- 3.2 default settings 
# number of k-folds
k_fold = 10

# segments
seg = 10


#----------- 3.3 Specify the trait and ncomp below it

# (KW)
trait = "weight_mg"
ncomp = 11

# # (ST - %)
# trait = "starch_PCT"
# ncomp = 9

# # (ST - mg)
# trait = "starch_mg"
# ncomp = 9

# # (SUC)
# trait = "sucrose_mg"
# ncomp = 4

# # (GC)
# trait = "glucose_PCT"
# ncomp = 8

# # (PG)
# trait = "PG_PCT"
# ncomp = 8

# # (TS)
# trait = "total_sugars_mg"
# ncomp = 9

# # (TC)
# trait = "total_carbohydrates_mg"
# ncomp = 10


n_seed=711



#----------- 3.4 Dataset loading

{


# dataset
if(trait=="weight_mg"){
  dataset = "w"
} else if(trait =="PG_PCT"){
  dataset = "pg"
} else if(trait %in% c("starch_PCT" ,"starch_mg", "total_sugars_mg","glucose_PCT","sucrose_mg", "total_carbohydrates_mg")){
  dataset = "s"
}

  
#----------- 3.5 Spectra pre-treatment

if(trait %in% c("starch_PCT","glucose_PCT", "PG_PCT")){
  spec = "snv"
}else{
  spec= "raw"
}

#----------- 3.6 Threshold for Coefficient of variation (CV = sd/mean). applied only for starch and sugars.
cv_threshold = TRUE  # applied only on dataset "S", no need to change
cv <- paste0("cv_", trait)
cvth = 0.11   # default. !note! weight cv threshold is 0.10
# wavelength: 940-1640nm


#----------- 3.7 read calibration dataset 

# sugars, starch
if(dataset == "s"){
  dat <- read.csv("sugars_WL.NIR.csv")
  if(cv_threshold == TRUE & trait != "total_carbohydrates_mg"){
    dat <- dat %>%  filter(dat[,cv] < cvth)
  }
  
  # weight
}else if(dataset == "w"){
  dat <- read.csv("weight_WL.NIR.csv") 
  d = dat[,grep("weight",colnames(dat))]
  d = cbind.data.frame(d, average = apply(d,1,mean),CV = apply(d,1,function(x) sd(x)/mean(x)))
 
  
  # Phytoglycogen
}else if(dataset == "pg"){
  dat <- read.csv("PG_WL.NIR.csv") 
}


#----------- 3.8 data.frame with genotype + measured + NIR
d <- dat[!(is.na(dat[,trait])), ]
nir <- select(d, num_range("X", 940:1640)) 
NIR <- as.matrix(nir)
data <- data.frame(d[,c("genotype",trait)], I(NIR)) 
data.g <- select(data, genotype, trait) %>% group_by(genotype) %>% summarise_all(mean)


#----------- 3.9 sample stats dataframe
sample <- data.frame(type = c("individual", "genotype"),
                     Size = c(nrow(data), nrow(data.g)),
                     Mean = c(round(mean(data[,trait]),2),round(mean(data[,trait]),2)),
                     SD = c(round(sd(data[,trait]),2),round(sd(data[,trait]),2)),
                     Min = c(round(min(data[,trait]),2),round(min(data[,trait]),2)),
                     Max = c(round(max(data[,trait]),2),round(max(data[,trait]),2)))
#print(sample)

###-----------------------------------
####---- 4. Deploy the model
###-----------------------------------


#----------- 4.1 Create a dataframe to record predicted 
# values. cols=(genotype, trait(measured), pred)

Res = c() # empty for results
set.seed(n_seed)

dat <- data %>%  mutate(set = sample(1:k_fold, nrow(data), replace = TRUE))

#----------- 4.2 Model
for(fold in 1:k_fold){
  #split dataset into test and train
  test = dat %>% filter(set == fold)
  train = dat %>% filter(set !=fold)
  
  # build train.pls model for each fold
  if(spec == "raw"){ # raw NIR
    train.pls <- plsr(get(trait) ~ NIR, ncomp = 25, data = train, validation ="CV", segments = seg)
  } else if(spec == "snv"){ # SNV(NIR)
    train.pls <- plsr(get(trait) ~ SNV(NIR), ncomp = 25, data = train, validation ="CV", segments = seg)
  }
  
  # predict test data (outer-validation) with train.pls
  t <- test %>% 
    select(genotype, all_of(trait)) %>% 
    mutate(predicted = predict(train.pls, ncomp = ncomp, newdata =test) )
  # data
  Res <- rbind(Res, t)   #genotype, trait(measured), pred.                                      
}

##----------- 4.3 Results
X = Res[,trait]
Y = Res[,"predicted"]

# Individual correlation: correlation, RMSEP, RPD
r = cor(X, Y) %>% round(2)
rmse = sqrt(mean((X - Y)^2)) %>% round(2)
rpd = RPD(X, rmse) %>% round(2)
individual <- c(r, rmse, rpd)

# grouped by genotype
Res.geno <- Res %>% group_by(genotype) %>% 
            summarise(across(.col = 1:2, mean)) %>%
            as.data.frame()


X = Res.geno[,trait]
Y = Res.geno[,"predicted"]

# Genotype correlation: correlation, RMSEP, RPD
r.g = cor(X, Y) %>% round(2)
rmse.g = sqrt(mean((X - Y)^2)) %>% round(2)
rpd.g  = RPD(X, rmse.g) %>% round(2)
genotype <- c(r.g, rmse.g, rpd.g)

# Grouping
stats <- data.frame(type =c("individual", "genotype"), 
                    r = c(r, r.g), 
                    RMSEP = c(rmse, rmse.g), 
                    RPD = c(rpd, rpd.g))
cat(trait)
print(stats)

}

