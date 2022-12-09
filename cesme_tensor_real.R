
#setwd("c:/Dropbox/Highd-CESME/simulations/code")

#source("cesme_tensor0729/R/cesme.R")
#source("cesme_tensor0729/R/utilities.R")

#library(cesmetensor)

#install.packages("BiocManager")
#BiocManager::install(version = "3.15")
#BiocManager::install("GEOquery")
#BiocManager::install("Biobase")

#setwd("c:/Dropbox/Highd-CESME/tensor/simulation")
#setwd("c:/dpbox/dropbox/Highd-CESME/tensor/simulation")
setwd("D:/Dropbox/Research/Statisitcs/Highd-CESME/tensor/simulation")
setwd("/Users/wenzhou/Dropbox/Research/Statisitcs/Highd-CESME/tensor/simulation")

install.packages("PARSE_0.1.0.tar.gz", repos=NULL, type="source")
#install.packages("rTensor_1.4.8.tar.gz", repos=NULL, type="source")
rm(list=ls())
library(mvtnorm)
library(mclust)
library(mixtools)
library(DiceKriging)
library(gtools)
library(MASS)
library(msda)
library(sparcl)
library(scvxclustr)
library(cesmetensor)
library(kolmim)
library(kernlab)
library(monreg)
library(coneproj)
library(sphunif)
library(ClustOfVar)
library(rTensor)
library(tensorsparse)
library(PARSE)

library(tidyverse)
library(janitor)
library(flextable)
library(hrbrthemes)
library("oro.nifti")
library(TRES)

error_clus=function(clus, true_clus){
  m=max(c(clus,true_clus))
  Ptmp=permutations(m,m)
  etmp=rep(0,factorial(m))
  for(i in 1:factorial(m)){
    etmp[i]=mean((as.numeric(factor(clus,level=Ptmp[i,]))-true_clus!=0))
  }
  return(min(etmp))
}



#######################################
# ADHD

# data loading

setwd("c:/Dropbox/Highd-CESME/tensor/simulation")
setwd("c:/dpbox/dropbox/Highd-CESME/tensor/simulation")

source("utilities.R")
source("cesme.R")
source("DynTensorCluster.R")
source("tcross.R")

KKI <- read_csv(file = "ADHD200_training_pheno/KKI/KKI_phenotypic.csv") %>% 
  mutate_at(vars(`Full2 IQ`), funs(replace_na(NA))) %>% 
  select(-starts_with("QC_"))


NYU <- read_csv(file = "ADHD200_training_pheno/NYU/NYU_phenotypic.csv") %>% 
  mutate_at(vars(Handedness)
            , funs(case_when(Handedness > 0.1 ~ 1,
                             Handedness < 1 ~ 0))) %>% 
  mutate_at(vars(`Full2 IQ`), funs(replace_na(NA))) %>% 
  select(-starts_with("QC_"))
# replacing all positive-valued Edinburgh Handedness scores with a 1 
# (right-handed) and all negative scores with a 0 (left-handed)
# replace N/A with NA values

OHSU <- read_csv(file = "ADHD200_training_pheno/OHSU/OHSU_phenotypic.csv") %>% 
  mutate_at(vars(`ADHD Index`, `Verbal IQ`, `Performance IQ`, `Full2 IQ`), funs(replace_na(NA))) %>% 
  select(-starts_with("QC_"))

Peking_1 <- read_csv(file = "ADHD200_training_pheno/Peking_1/Peking_1_phenotypic.csv") %>% 
  mutate_at(vars(`Full2 IQ`), funs(replace_na(NA))) %>% 
  select(-starts_with("QC_"))

Peking_2 <- read_csv(file = "ADHD200_training_pheno/Peking_2/Peking_2_phenotypic.csv") %>% 
  mutate_at(vars(`Full2 IQ`), funs(replace_na(NA))) %>% 
  select(-starts_with("QC_"))

Peking_3 <- read_csv(file = "ADHD200_training_pheno/Peking_3/Peking_3_phenotypic.csv") %>% 
  mutate_at(vars(`Full2 IQ`), funs(replace_na(NA))) %>% 
  select(-starts_with("QC_"))

Pitts <- read_csv(file = "ADHD200_training_pheno/Pittsburgh/Pittsburgh_phenotypic.csv") %>% 
  mutate_at(vars(`ADHD Measure`, `ADHD Index`, Inattentive, `Hyper/Impulsive`), funs(replace_na(NA))) %>% 
  select(-starts_with("QC_"))

WashU <- read_csv(file = "ADHD200_training_pheno/WashU/WashU_phenotypic.csv") %>% 
  mutate_at(vars(`ADHD Measure`, `ADHD Index`, Inattentive, `Hyper/Impulsive`), funs(replace_na(NA))) %>% 
  select(-starts_with("QC_"))

NeuroIMAGE <- read_csv(file = "ADHD200_training_pheno/NeuroIMAGE/NeuroIMAGE_phenotypic.csv") %>% 
  mutate_at(vars(`ADHD Measure`, `ADHD Index`, Inattentive, `Hyper/Impulsive`), funs(replace_na(NA))) %>% 
  select(-starts_with("QC_"))

adhd_200_full_join <- bind_rows(KKI, NYU, OHSU, Peking_1,Peking_2,Peking_3, Pitts,WashU,NeuroIMAGE)


adhd_200 <- adhd_200_full_join %>% clean_names(case = "snake") %>% 
  mutate_at(vars(gender, handedness, dx, iq_measure, med_status), funs(factor)) %>% 
  mutate_all(funs(ifelse(.== -999, NA, .)))

adhd_200$gender[which(is.na(adhd_200$gender)==1)] = -100
gender_labels <- c("Missing", "Female", "Male")
adhd_200$gender <- factor(adhd_200$gender, labels = gender_labels)

hand_labels <- c("L", "R", "AmbiD")
adhd_200$handedness <- factor(adhd_200$handedness, labels = hand_labels)

dx_labels <- c("Typical", "ADHD_c", "ADHD_H", "ADHD_I")
adhd_200$dx <- factor(adhd_200$dx, labels = dx_labels)

adhd_200$iq_measure[which(is.na(adhd_200$iq_measure)==1)] = 5
adhd_200$iq_measure[which(adhd_200$iq_measure==4)] = 3

#x <- factor(c(5, 2, 3, 4), exclude = NULL) 
#IQ_labels <- c("Missing", "WASI", "WISCCR", "WASI_sub")

x <- factor(c(5, 2, 3), exclude = NULL) 
IQ_labels <- c("Missing", "WASI", "WISCCR")

adhd_200$iq_measure <- factor(adhd_200$iq_measure, labels = IQ_labels, levels = x) # This dataset only has 1 Level

adhd_200$med_status[which(is.na(adhd_200$med_status)==1)] = 1
x1 <- factor(c(1, 2, 3), exclude = NULL)

md_labels <- c("Missing", "Medicated", "Non-Medicated")
adhd_200$med_status <- factor(adhd_200$med_status, labels = md_labels, levels = x1)

site_labels = c("KKI","NYU","OHSU","Peking","Pittsburgh","WashU","NeuroIMAGE")
adhd_200$site = factor(adhd_200_full_join$Site, labels = site_labels, levels = unique(adhd_200_full_join$Site))

adhd_200_full_join = adhd_200
save(adhd_200_full_join,file = "ADHD_phenotypic.Rdata")

# summarized phenotypic

load(file = "ADHD_phenotypic.Rdata")

pdf("real_proportion.pdf", width = 6, height = 4)
inspectdf::inspect_cat(adhd_200) %>% inspectdf::show_plot(text_labels = TRUE, label_size = 5, col_palette = 0)
dev.off()

# regression coefficient

pred_continuous = adhd_200_full_join[,4]
B1 = adhd_200_full_join[,c(2,3,5,12,17)]
B2 = cbind(y = 0, as.data.frame(B1))
pred = cbind(as.matrix(pred_continuous), model.matrix(y~.,data = B2)[,-1])

save(pred,file = "ADHD_predictor.Rdata")

# summary

NYU.motion <- read.csv(file = "ADHD200_training_motion/NYU_motion.csv")

KKI.motion <- read.csv(file = "ADHD200_training_motion/KKI_motion.csv")
KKI.motion[1,]

OHSU.motion <- read.csv(file = "ADHD200_training_motion/OHSU_motion.csv")
OHSU.motion[,1]

# load data

load(file = "ADHD_phenotypic.Rdata")

filepath = paste0(adhd_200_full_join$site,"_",adhd_200_full_join$`scan_dir_id`,"_1")
filename = paste0("C:/tensor_data/ADHD/",filepath,"/anat_1/NIfTI/rest.nii.gz")
savename = paste0("C:/tensor_data/ADHD_anat",filepath,".Rdata")

filename = paste0("C:/tensor_data/ADHD/",filepath,"/anat_1/NIfTI/rest.nii.gz")
savename = paste0("C:/dropbox/dropbox/tensor_data/ADHD_anat/raw/",filepath,".Rdata")

savename0 = paste0("C:/dropbox/dropbox/tensor_data/ADHD_anat/center/",adhd_200_full_join$site,"_",
                  adhd_200_full_join$`scan_dir_id`,".Rdata")
savename0 = paste0("C:/Users/lozhang/Documents/dropbox/tensor_data/ADHD_anat/center/",adhd_200_full_join$site,"_",
                   adhd_200_full_join$`scan_dir_id`,".Rdata")

d_select = c(121,145,121)

i = 1
while(i<=776){
  load(savename[i])
  dtmp = dim(Mtmp)
  if(any(dtmp<d_select)){
    m = which(dtmp<d_select)
    Mtmp0 = array(0,dim=c(100+dtmp[1],dtmp[-1]))
    Mtmp0[-c(1:50,(dim(Mtmp0)[1]-49):dim(Mtmp0)[1]),,] = Mtmp
    Mtmp = Mtmp0
    dtmp = dim(Mtmp)
  }
  dif = (dtmp-d_select)/2
  ltmp = ceiling(dif)
  utmp = floor(dtmp-dif)
  if(sum(utmp-ltmp+1==d_select)!=3){
    break
  }else{
    Ttmp = Mtmp[ltmp[1]:utmp[1],ltmp[2]:utmp[2],ltmp[3]:utmp[3]]
    save(Ttmp,file = savename0[i])
    i = i+1
  }
  if(i%%10==0){
    print(i)
  }
}

dim.store = c()

save(dim.store,file = "ADHD_anat_dim.Rdata")

load(file = "ADHD_anat_dim.Rdata")
i = dim(dim.store)[1]+1

while(i<=776){
  load(file = "ADHD_anat_dim.Rdata")
  d1 = readNIfTI(filename[i])
  dim.store = rbind(dim.store, dim(d1))
  rownames(dim.store)[i] = filepath[i]
  save(dim.store,file = "ADHD_anat_dim.Rdata")
  i = i+1
}

load(file = "ADHD_anat_dim.Rdata")

dim.min = matrix(0,7,3)
rownames(dim.min) = levels(adhd_200_full_join$Site)

for(i in 1:7){
  dim.min[i,] = apply(dim.store[as.numeric(adhd_200_full_join$Site)==i,],2,min)
}

xtable(dim.min,digits = 0)

apply(dim.store,2,min)

# r=40,20,20
m = c(40,20,20)
matrix.store = array(0,dim = c(776,m))
diff.store = rep(0,776)
approx.store = list(matrix.store = matrix.store, diff.store = diff.store)

save(matrix.store,file = "ADHD_anat_matrix_r40_20_10.Rdata")

i = 1

while(i<=776){
  load(savename0[i])
  Omega = list()
  A = as.tensor(Ttmp)
  A.arm = list()
  for(j in 1:3){
    A.tmp = k_unfold(A,m=j)@data
    s.tmp = apply(A.tmp,1,var)
    Omega[[j]] = sort(order(s.tmp,decreasing = TRUE)[1:m[j]])
    A.arm[[j]] = k_unfold(A,m = j)@data[,Omega[[j]]]
  }
  B = A[Omega[[1]],Omega[[2]],Omega[[3]]]
  approx.store$matrix.store[i,,,] = B@data
  result = tcross_cv(B = B,Arm = A.arm)
  A.approx = B
  for(j in 1:3){
    A.approx = ttm(A.approx,result[[3]][[j]],m=j)
  }
  approx.store$diff.store[i] = sqrt(sum((A-A.approx)@data^2))
  if(i%%50==0){
    print(i)
  }
  i = i+1
}

save(approx.store,file = "ADHD_anat_matrix_r40_20_20.Rdata")

# r=20,40,20
m = c(20,40,20)
matrix.store = array(0,dim = c(776,m))
diff.store = rep(0,776)
approx.store = list(matrix.store = matrix.store, diff.store = diff.store)

save(matrix.store,file = "ADHD_anat_matrix_r20_40_20.Rdata")

i = 1

while(i<=776){
  load(savename0[i])
  Omega = list()
  A = as.tensor(Ttmp)
  A.arm = list()
  for(j in 1:3){
    A.tmp = k_unfold(A,m=j)@data
    s.tmp = apply(A.tmp,1,var)
    Omega[[j]] = sort(order(s.tmp,decreasing = TRUE)[1:m[j]])
    A.arm[[j]] = k_unfold(A,m = j)@data[,Omega[[j]]]
  }
  B = A[Omega[[1]],Omega[[2]],Omega[[3]]]
  approx.store$matrix.store[i,,,] = B@data
  result = tcross_cv(B = B,Arm = A.arm)
  A.approx = B
  for(j in 1:3){
    A.approx = ttm(A.approx,result[[3]][[j]],m=j)
  }
  approx.store$diff.store[i] = sqrt(sum((A-A.approx)@data^2))
  if(i%%50==0){
    print(i)
  }
  i = i+1
}

save(approx.store,file = "ADHD_anat_matrix_r20_40_20.Rdata")

# r=20,20,40
m = c(20,20,40)
matrix.store = array(0,dim = c(776,m))
diff.store = rep(0,776)
approx.store = list(matrix.store = matrix.store, diff.store = diff.store)

save(matrix.store,file = "ADHD_anat_matrix_r20_20_40.Rdata")

i = 1

while(i<=776){
  #load(file = "ADHD_anat_matrix_r10_20_40.Rdata")
  load(savename0[i])
  Omega = list()
  A = as.tensor(Ttmp)
  A.arm = list()
  for(j in 1:3){
    A.tmp = k_unfold(A,m=j)@data
    s.tmp = apply(A.tmp,1,var)
    Omega[[j]] = sort(order(s.tmp,decreasing = TRUE)[1:m[j]])
    A.arm[[j]] = k_unfold(A,m = j)@data[,Omega[[j]]]
  }
  B = A[Omega[[1]],Omega[[2]],Omega[[3]]]
  approx.store$matrix.store[i,,,] = B@data
  result = tcross_cv(B = B,Arm = A.arm)
  A.approx = B
  for(j in 1:3){
    A.approx = ttm(A.approx,result[[3]][[j]],m=j)
  }
  approx.store$diff.store[i] = sqrt(sum((A-A.approx)@data^2))
  #save(approx.store,file = "ADHD_anat_matrix_r10_20_40.Rdata")
  if(i%%50==0){
    print(i)
  }
  i = i+1
}

save(approx.store,file = "ADHD_anat_matrix_r20_20_40.Rdata")

# r=40,20,10
m = c(40,20,10)
matrix.store = array(0,dim = c(776,m))
diff.store = rep(0,776)
approx.store = list(matrix.store = matrix.store, diff.store = diff.store)

save(matrix.store,file = "ADHD_anat_matrix_r40_20_10.Rdata")

i = 1

while(i<=776){
  load(savename0[i])
  Omega = list()
  A = as.tensor(Ttmp)
  A.arm = list()
  for(j in 1:3){
    A.tmp = k_unfold(A,m=j)@data
    s.tmp = apply(A.tmp,1,var)
    Omega[[j]] = sort(order(s.tmp,decreasing = TRUE)[1:m[j]])
    A.arm[[j]] = k_unfold(A,m = j)@data[,Omega[[j]]]
  }
  B = A[Omega[[1]],Omega[[2]],Omega[[3]]]
  approx.store$matrix.store[i,,,] = B@data
  result = tcross_cv(B = B,Arm = A.arm)
  A.approx = B
  for(j in 1:3){
    A.approx = ttm(A.approx,result[[3]][[j]],m=j)
  }
  approx.store$diff.store[i] = sqrt(sum((A-A.approx)@data^2))
  if(i%%50==0){
    print(i)
  }
  i = i+1
}

save(approx.store,file = "ADHD_anat_matrix_r10_40_20.Rdata")

# r=10,40,20
m = c(10,40,20)
matrix.store = array(0,dim = c(776,m))
diff.store = rep(0,776)
approx.store = list(matrix.store = matrix.store, diff.store = diff.store)

save(matrix.store,file = "ADHD_anat_matrix_r10_40_20.Rdata")

i = 1

while(i<=776){
  load(savename0[i])
  Omega = list()
  A = as.tensor(Ttmp)
  A.arm = list()
  for(j in 1:3){
    A.tmp = k_unfold(A,m=j)@data
    s.tmp = apply(A.tmp,1,var)
    Omega[[j]] = sort(order(s.tmp,decreasing = TRUE)[1:m[j]])
    A.arm[[j]] = k_unfold(A,m = j)@data[,Omega[[j]]]
  }
  B = A[Omega[[1]],Omega[[2]],Omega[[3]]]
  approx.store$matrix.store[i,,,] = B@data
  result = tcross_cv(B = B,Arm = A.arm)
  A.approx = B
  for(j in 1:3){
    A.approx = ttm(A.approx,result[[3]][[j]],m=j)
  }
  approx.store$diff.store[i] = sqrt(sum((A-A.approx)@data^2))
  if(i%%50==0){
    print(i)
  }
  i = i+1
}

save(approx.store,file = "ADHD_anat_matrix_r10_40_20.Rdata")

# r=10,20,40
m = c(10,20,40)
matrix.store = array(0,dim = c(776,m))
diff.store = rep(0,776)
approx.store = list(matrix.store = matrix.store, diff.store = diff.store)

save(matrix.store,file = "ADHD_anat_matrix_r10_20_40.Rdata")

i = 1

while(i<=776){
  #load(file = "ADHD_anat_matrix_r10_20_40.Rdata")
  load(savename0[i])
  Omega = list()
  A = as.tensor(Ttmp)
  A.arm = list()
  for(j in 1:3){
    A.tmp = k_unfold(A,m=j)@data
    s.tmp = apply(A.tmp,1,var)
    Omega[[j]] = sort(order(s.tmp,decreasing = TRUE)[1:m[j]])
    A.arm[[j]] = k_unfold(A,m = j)@data[,Omega[[j]]]
  }
  B = A[Omega[[1]],Omega[[2]],Omega[[3]]]
  approx.store$matrix.store[i,,,] = B@data
  result = tcross_cv(B = B,Arm = A.arm)
  A.approx = B
  for(j in 1:3){
    A.approx = ttm(A.approx,result[[3]][[j]],m=j)
  }
  approx.store$diff.store[i] = sqrt(sum((A-A.approx)@data^2))
  #save(approx.store,file = "ADHD_anat_matrix_r10_20_40.Rdata")
  if(i%%50==0){
    print(i)
  }
  i = i+1
}

save(approx.store,file = "ADHD_anat_matrix_r10_20_40.Rdata")

# r=30
m = 10
matrix.store = array(0,dim = c(776,rep(m,3)))
diff.store = rep(0,776)
approx.store = list(matrix.store = matrix.store, diff.store = diff.store)

save(matrix.store,file = "ADHD_anat_matrix_r10.Rdata")

i = 1

while(i<=776){
  load(savename0[i])
  Omega = list()
  A = as.tensor(Ttmp)
  A.arm = list()
  for(j in 1:3){
    A.tmp = k_unfold(A,m=j)@data
    s.tmp = apply(A.tmp,1,var)
    Omega[[j]] = sort(order(s.tmp,decreasing = TRUE)[1:m])
    A.arm[[j]] = k_unfold(A,m = j)@data[,Omega[[j]]]
  }
  B = A[Omega[[1]],Omega[[2]],Omega[[3]]]
  approx.store$matrix.store[i,,,] = B@data
  result = tcross_cv(B = B,Arm = A.arm)
  A.approx = B
  for(j in 1:3){
    A.approx = ttm(A.approx,result[[3]][[j]],m=j)
  }
  approx.store$diff.store[i] = sqrt(sum((A-A.approx)@data^2))
  if(i%%50==0){
    print(i)
  }
  i = i+1
}

save(approx.store,file = "ADHD_anat_matrix_r10.Rdata")


# r=10
m = 10
matrix.store = array(0,dim = c(776,rep(m,3)))

save(matrix.store,file = "ADHD_anat_matrix_r10.Rdata")

load(file = "ADHD_anat_matrix_r10.Rdata")
i = 575

while(i<=776){
  load(file = "ADHD_anat_matrix_r10.Rdata")
  d1 = readNIfTI(filename[i])
  Omega = list()
  A = as.tensor(d1@.Data)
  for(j in 1:3){
    A.tmp = k_unfold(A,m=j)@data
    s.tmp = apply(A.tmp,1,var)
    Omega[[j]] = sort(order(s.tmp,decreasing = TRUE)[1:m])
  }
  matrix.store[i,,,] = A[Omega[[1]],Omega[[2]],Omega[[3]]]@data
  save(matrix.store,file = "ADHD_anat_matrix_r10.Rdata")
  print(i)
  i = i+1
}

load(file = "ADHD_anat_matrix_r10.Rdata")

A1 = as.tensor(matrix.store)
A2 = k_unfold(A1,m=1)
A3 = k_fold(A2,m=4,modes = c(rep(m,3),776))
A4 = A3@data

load(file = "ADHD_predictor.Rdata")

fit_trr = TRR.fit(t(pred),A4)

res.store = fit_trr$residuals

save(res.store,file = "ADHD_anat_matrix_r10_res.Rdata")

# r=20
m = 20
matrix.store = array(0,dim = c(776,rep(m,3)))

save(matrix.store,file = "ADHD_anat_matrix_r20.Rdata")

load(file = "ADHD_anat_matrix_r20.Rdata")
i = 575

while(i<=776){
  load(file = "ADHD_anat_matrix_r20.Rdata")
  d1 = readNIfTI(filename[i])
  Omega = list()
  A = as.tensor(d1@.Data)
  for(j in 1:3){
    A.tmp = k_unfold(A,m=j)@data
    s.tmp = apply(A.tmp,1,var)
    Omega[[j]] = sort(order(s.tmp,decreasing = TRUE)[1:m])
  }
  matrix.store[i,,,] = A[Omega[[1]],Omega[[2]],Omega[[3]]]@data
  save(matrix.store,file = "ADHD_anat_matrix_r20.Rdata")
  print(i)
  i = i+1
}

load(file = "ADHD_predictor.Rdata")

###

m = c(20,40,20)
load(file = "ADHD_anat_matrix_r20_40_20.Rdata")

A1 = as.tensor(log(approx.store$matrix.store+1))
A2 = k_unfold(A1,m=1)
A3 = k_fold(A2,m=4,modes = c(m,776))
A4 = A3@data

fit_trr = TRR.fit(t(pred),A4)

approx.store$res.store = fit_trr$residuals

save(approx.store,file = "ADHD_anat_matrix_r20_40_20.Rdata")

###

m = c(40,20,20)
load(file = "ADHD_anat_matrix_r40_20_20.Rdata")

A1 = as.tensor(log(approx.store$matrix.store+1))
A2 = k_unfold(A1,m=1)
A3 = k_fold(A2,m=4,modes = c(m,776))
A4 = A3@data

fit_trr = TRR.fit(t(pred),A4)

approx.store$res.store = fit_trr$residuals

save(approx.store,file = "ADHD_anat_matrix_r40_20_20.Rdata")

###

m = c(20,20,40)
load(file = "ADHD_anat_matrix_r20_20_40.Rdata")

A1 = as.tensor(log(approx.store$matrix.store+1))
A2 = k_unfold(A1,m=1)
A3 = k_fold(A2,m=4,modes = c(m,776))
A4 = A3@data

fit_trr = TRR.fit(t(pred),A4)

approx.store$res.store = fit_trr$residuals

save(approx.store,file = "ADHD_anat_matrix_r20_20_40.Rdata")


# r=10
m = 10
p.all = rep(m,3)
load(file = "ADHD_anat_matrix_r10_res.Rdata")

# r=20
m = 20
p.all = rep(m,3)
load(file = "ADHD_anat_matrix_r20_res.Rdata")

# r=30
m = 30
p.all = rep(m,3)
load(file = "ADHD_anat_matrix_r30.Rdata")

# r=10,40,20
p.all = c(10,40,20)
load(file = "ADHD_anat_matrix_r10_40_20.Rdata")

# r=40,20,10
p.all = c(40,20,10)
load(file = "ADHD_anat_matrix_r40_20_10.Rdata")

# r=10,20,40
p.all = c(10,20,40)
load(file = "ADHD_anat_matrix_r10_20_40.Rdata")

# r=20,40,20
p.all = c(20,40,20)
load(file = "ADHD_anat_matrix_r20_40_20.Rdata")

# r=40,20,20
p.all = c(40,20,20)
load(file = "ADHD_anat_matrix_r40_20_20.Rdata")

# r=20,20,40
p.all = c(20,20,40)
load(file = "ADHD_anat_matrix_r10_20_40.Rdata")


# eliminate adhd_h
adhd_h = which(adhd_200_full_join$dx=="ADHD_H")
nn = 776-length(adhd_h)

# true label
true_id = as.numeric(adhd_200_full_join$dx)[-adhd_h]
true_id[which(true_id==4)] = 3

A1 = approx.store$res.store
A2 = k_unfold(A1,m=4)
A3 = k_fold(A2,m=1,modes = c(776,p.all))
X = A3[-adhd_h,,,]

iter = 2022
nc = 3
# spectral clustering
set.seed(iter)
sc_clus <- as.numeric(specc(k_unfold(X,m = 1)@data, centers=nc))
sc_err <- error_clus(sc_clus, true_id)


# k-means
set.seed(iter)
km_clus <- kmeans(k_unfold(X,m = 1)@data, centers=nc)$cluster
km_err <- error_clus(km_clus, true_id)


# sparse k-means
set.seed(iter)
skm_clus <- KMeansSparseCluster(k_unfold(X,m = 1)@data, K = nc ,wbounds = 1.1)[[1]]$Cs
skm_err <- error_clus(skm_clus, true_id)

# TBM
set.seed(iter)
tbm.fit = tbmClustering(k_fold(k_unfold(X,m = 1),m=1,modes = c(nn,p.all[1],p.all[2]*p.all[3]))@data,nc,1,1)
tbm_clus = tbm.fit$Cs
tbm_err <- error_clus(tbm_clus, true_id)

# DTC
set.seed(iter)
X.DTC = k_fold(k_unfold(X,m = 1),m=4,modes = c(p.all,nn))
out_stdtruncate = stdtruncate(X.DTC@data,c(p.all/2,nn),c(p.all/2,nn),2)
dtc_clus = mydtc(out_stdtruncate,K = nc)$membership
dtc_err <- error_clus(dtc_clus, true_id)


# CESME
set.seed(iter)
npn.def <- npn.clust.tensor(X, K=nc, N=100, 
                            thre=1/(4*(NROW(X))^0.25 *sqrt(pi*log(NROW(X)))), 
                            C0=1, C2=1, kap=0.8, kmnstart=10, z.initial=sc_clus, g.method=1)
cesme_err <- error_clus(npn.def$z[npn.def$iter+1,,1], true_id)


# DEEM
set.seed(iter)
deem.def <- npn.clust.tensor(X, K=nc, N=100, 
                             thre=1/(4*(NROW(X))^0.25 *sqrt(pi*log(NROW(X)))), 
                             C0=1, C2=1, kap=0.8, kmnstart=10, z.initial=sc_clus, g.method=4)
deem_err <- error_clus(deem.def.def$z[deem.def.def$iter+1,,1], true_id)


xtable(t(as.matrix(c(km_err, skm_err, sc_err,  tbm_err, dtc_err)*100)))

# true label
true_id = as.numeric(adhd_200_full_join$dx)
nn = 776

A1 = approx.store$res.store
A2 = k_unfold(A1,m=4)
A3 = k_fold(A2,m=1,modes = c(776,p.all))
X = A3

iter = 2022
nc = 4
# spectral clustering
set.seed(iter)
sc_clus <- as.numeric(specc(k_unfold(X,m = 1)@data, centers=4))
sc_err <- error_clus(sc_clus, true_id)


# k-means
set.seed(iter)
km_clus <- kmeans(k_unfold(X,m = 1)@data, centers=nc)$cluster
km_err <- error_clus(km_clus, true_id)


# sparse k-means
set.seed(iter)
skm_clus <- KMeansSparseCluster(k_unfold(X,m = 1)@data, K = nc ,wbounds = 1.1)[[1]]$Cs
skm_err <- error_clus(skm_clus, true_id)

# TBM
set.seed(iter)
tbm.fit = tbmClustering(k_fold(k_unfold(X,m = 1),m=1,modes = c(nn,p.all[1],p.all[2]*p.all[3]))@data,nc,1,1)
tbm_clus = tbm.fit$Cs
tbm_err <- error_clus(tbm_clus, true_id)

# DTC
set.seed(iter)
X.DTC = k_fold(k_unfold(X,m = 1),m=4,modes = c(p.all,nn))
out_stdtruncate = stdtruncate(X.DTC@data,c(p.all/2,nn),c(p.all/2,nn),2)
dtc_clus = mydtc(out_stdtruncate,K = nc)$membership
dtc_err <- error_clus(dtc_clus, true_id)


# CESME
set.seed(iter)
npn.def <- npn.clust.tensor(X, K=nc, N=100, 
                            thre=1/(4*(NROW(X))^0.25 *sqrt(pi*log(NROW(X)))), 
                            C0=1, C2=1, kap=0.8, kmnstart=10, z.initial=sc_clus, g.method=1)
cesme_err <- error_clus(npn.def$z[npn.def$iter+1,,1], true_id)


# DEEM
set.seed(iter)
deem.def <- npn.clust.tensor(X, K=nc, N=100, 
                             thre=1/(4*(NROW(X))^0.25 *sqrt(pi*log(NROW(X)))), 
                             C0=1, C2=1, kap=0.8, kmnstart=10, z.initial=sc_clus, g.method=4)
deem_err <- error_clus(deem.def.def$z[deem.def.def$iter+1,,1], true_id)


xtable(t(as.matrix(c(km_err, skm_err, sc_err,  tbm_err, dtc_err)*100)))







