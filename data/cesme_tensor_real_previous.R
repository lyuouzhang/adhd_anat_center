
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


source("utilities.R")
source("cesme.R")
source("DynTensorCluster.R")

error_clus=function(clus, true_clus){
  m=max(c(clus,true_clus))
  Ptmp=permutations(m,m)
  etmp=rep(0,factorial(m))
  for(i in 1:factorial(m)){
    etmp[i]=mean((as.numeric(factor(clus,level=Ptmp[i,]))-true_clus!=0))
  }
  return(min(etmp))
}

##########################################################################

#gds1083 

X = read.table("GDS1083.txt",header = TRUE)[,-(1:2)]

X.tensor_1 = k_fold((as.matrix(X)),m = 3,modes = c(4,27,1124))

X.tensor = k_fold(k_unfold(X.tensor_1,m=2),m = 1,modes = c(27,4,1124))

X.matrix = k_unfold(X.tensor,m=1)

true_id = rep(1:3,9)

# spectral

set.seed(2001)

sc.clust = as.numeric(specc(X.matrix@data,centers=3))

sc.error = error_clus(sc.clust,true_id)

# CESME

set.seed(2001)

npn.def <- npn.clust.tensor(X, K=3, N=100, 
                            thre=1/(4*(27)^0.25 *sqrt(pi*log(27))), 
                            C0=1, C2=1, kap=0.8, kmnstart=10, z.initial=sc.clust, g.method=1)
cesme.err <- error_clus(npn.def$z[npn.def$iter+1,,1], true_id)

# DEEM

set.seed(2001)

deem.def <- npn.clust.tensor(X, K=3, N=100, 
                             thre=1/(4*(27)^0.25 *sqrt(pi*log(27))), 
                             C0=1, C2=1, kap=0.8, kmnstart=10, z.initial=sc.clust, g.method=4)
deem.err <- error_clus(deem.def$z[deem.def$iter+1,,1], true_id)


# K-means

set.seed(2001)

km.clust = kmeans(X.matrix@data,centers=3,nstart = 10)$cluster

km.error = error_clus(km.clust,true_id)

# sparse K-means

set.seed(2001)

skm.clust = KMeansSparseCluster(X.matrix@data,K=3,wbounds = 1.1)[[1]]$Cs

skm.error = error_clus(skm.clust,true_id)

# DTC

set.seed(2001)

X.DTC = k_fold(k_unfold(X.tensor,m=2),m = 2,modes = c(1124,4,27))
out_stdtruncate = stdtruncate(X.DTC@data,sparse_para=c(1124/2,4,27),smooth_para=c(1124/2,4,27),Rank=1,
                              is_first2modes_symmetric = FALSE)
dtc.clus = mydtc(out_stdtruncate,K = 3)$membership
dtc.err <- error_clus(dtc.clus, true_id)

# TBM

set.seed(2001)

tbm.fit = tbmClustering(X.tensor@data,3,1,1)
tbm.clus = tbm.fit$Cs
tbm.err <- error_clus(tbm.clus, true_id)

# APFP

set.seed(2001)

apfp.fit = apfp(K = 3, lambda = 1:3, y = X.matrix@data, kms.nstart = 10)
apfp.clus = apfp.fit$s.hat.best
apfp.err <- error_clus(apfp.clus, true_id)

# EM

set.seed(2001)

em.clus = Mclust(X.matrix@data,G=3,modelNames = "EII")$classification
em.err <- error_clus(em.clus, true_id)




#################################################

# GSE69761

X = read.table("GSE69761_36188g-148c-fpkm.txt",header = TRUE)

X.name = X[,1]
X.data = X[,-1]

dim(X.data)
X.name[10000:10010]


##################################################

# spatial transcript
library(dplyr)

setwd("c:/Dropbox/Highd-CESME/tensor/simulation/FIST-master")
setwd("C:/Users/lozhang/Documents/dropbox/tensor_data/FIST-master")

# Package required to run this script
packages <- c("argparser", "Matrix", "R.matlab","data.table")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages],repos = "http://cran.us.r-project.org")
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# arguments parser
p <- arg_parser("Convert spatial transcriptomics matrix data into tensor format")
p <- add_argument(p, "--input", help="Path to input folder", type="character")
p <- add_argument(p, "--output", help="Path to output folder", type="character")
argv <- parse_args(p)

argv$input = "C:/Users/lozhang/Documents/dropbox/tensor_data/FIST-master/FIST_Tutorial"
argv$output = "C:/Users/lozhang/Documents/dropbox/tensor_data/FIST_Tutorial_Output"

# List all tissues for convertion
tissue_names <- list.dirs(paste0(argv$input, "/"), full.names = FALSE, recursive = FALSE)

#tn = tissue_names[4]

for(tn in tissue_names){
  
  
  # Set path to filtered feature-barcode matrix data and spatial coordinates
  matrix_dir <- paste0(argv$input, "/", tn, "/filtered_feature_bc_matrix/")
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  spatial_info_dir <- paste0(argv$input, "/", tn, "/spatial/")
  barcode_position.path <- paste0(spatial_info_dir, "tissue_positions_list.csv")
  
  #cluster.path = paste0(argv$input, "/", tn, "/analysis/clustering/graphclust/clusters.csv")
  
  # Read filtered feature-barcode matrix data
  mat <- readMM(file = matrix.path)
  feature.names <- read.delim(features.path,
                              header = FALSE,
                              stringsAsFactors = FALSE)
  barcode.names <- read.delim(barcode.path,
                              header = FALSE,
                              stringsAsFactors = FALSE)
  colnames(mat) <- barcode.names$V1
  rownames(mat) <- feature.names$V1
  
  #clus = read.csv(cluster.path)
  
  #total_counts = colSums(mat)
  
  # Read spatial coordinates
  barcode_position <- read.delim(barcode_position.path,
                                 header = FALSE,
                                 sep = ",",
                                 stringsAsFactors = FALSE)  
  
  #X = data.frame(x_coords = barcode_position$V3[match(colnames(mat), barcode_position$V1)]+1,
  #               y_coords = barcode_position$V4[match(colnames(mat), barcode_position$V1)]+1,
  #               total_counts = total_counts,
  #               clus = clus$Cluster)
  
  #save(X,file = paste0(tn,".Rdata"))

  # Align spatial coordinates between rows
  t <- data.frame(x_coords = barcode_position$V3[match(colnames(mat)[mat@j+1], barcode_position$V1)]+1,
                  y_coords = barcode_position$V4[match(colnames(mat)[mat@j+1], barcode_position$V1)]+1,
                  variable = mat@i+1,
                  value = mat@x)
  
  #length(unique(t$x_coords))
  #length(unique(t$y_coords))
  
  t$x_aligned_coords <- t$x_coords
  t$y_aligned_coords <- t$y_coords
  
  t$y_aligned_coords[which(t$x_coords%%2==1)] <- t$y_coords[which(t$x_coords%%2==1)]%/%2+1
  t$y_aligned_coords[which(t$x_coords%%2==0)] <- t$y_coords[which(t$x_coords%%2==0)]%/%2
  
  t <- t[order(t$x_aligned_coords, t$y_aligned_coords), ]
  
  t$x_aligned_coords = t$x_aligned_coords - min(t$x_aligned_coords) + 1
  t$y_aligned_coords = t$y_aligned_coords - min(t$y_aligned_coords) + 1
  
  unique.value = unique(t$variable)
  
  X = array(0,dim = c(length(unique.value),max(unique(t$y_aligned_coords)),max(unique(t$x_aligned_coords))))
  

  for(j in 1:max(t$y_aligned_coords)){
    for(k in 1:max(t$x_aligned_coords)){
      A = t[intersect(which(t$y_aligned_coords==j), which(t$x_aligned_coords==k)),c(3,4)]
      if(dim(A)[1]==0){
        next
      }
      B = aggregate(value~variable,data=A,sum)
      X[match(B[,1],unique.value),j,k] = B[,2]
    }
  }
  
  zero = apply(X<3,1,sum)
  length(which(zero==0))
  
  
  #sum(t$value>3)
  length(unique(t$x_aligned_coords))
  table(t$x_aligned_coords)
  length(unique(t$y_aligned_coords))
  table(t$y_aligned_coords)
  length(unique(t$variable))
  t$x_aligned_coords[which(t$variable==1)]
  
  
  

  
  
  # Write gene names of tensor slice for each tissue to output folder
  fwrite(feature.names, paste0(argv$output, "/", tn, "_gene.csv"))
  
  # Write gene names of tensor slice for each tissue to output folder
  fwrite(clus, paste0(argv$output, "/", tn, "_clus.csv"))
  
  
  
  
  
  
  # Write tensor data to output folder
  writeMat(paste0(argv$output, "/", tn, ".mat"), V = t, X=78, Y=64, Z=nrow(mat))
}

###############################################################################


setwd("c:/Dropbox/Highd-CESME/tensor/simulation")

load("V1_Breast_Cancer_Block_A_Section_1.Rdata")

dim(X)[1]
length(unique(X$x_coords))
length(unique(X$y_coords))
length(unique(X$clus))


load("V1_Breast_Cancer_Block_A_Section_2.Rdata")

dim(X)[1]
length(unique(X$x_coords))
length(unique(X$y_coords))
length(unique(X$clus))


load("V1_Human_Heart.Rdata")

dim(X)[1]
length(unique(X$x_coords))
length(unique(X$y_coords))
length(unique(X$clus))


load("V1_Human_Lymph_Node.Rdata")

dim(X)[1]
length(unique(X$x_coords))
length(unique(X$y_coords))
length(unique(X$clus))


load("V1_Mouse_Kidney.Rdata")

dim(X)[1]
length(unique(X$x_coords))
length(unique(X$y_coords))
length(unique(X$clus))


load("V1_Adult_Mouse_Brain.Rdata")

dim(X)[1]
length(unique(X$x_coords))
length(unique(X$y_coords))
length(unique(X$clus))


load("V1_Mouse_Brain_Sagittal_Posterior.Rdata")

dim(X)[1]
length(unique(X$x_coords))
length(unique(X$y_coords))
length(unique(X$clus))


load("V1_Mouse_Brain_Sagittal_Posterior_Section_2.Rdata")

dim(X)[1]
length(unique(X$x_coords))
length(unique(X$y_coords))
length(unique(X$clus))


load("V1_Mouse_Brain_Sagittal_Anterior.Rdata")

dim(X)[1]
length(unique(X$x_coords))
length(unique(X$y_coords))
length(unique(X$clus))


load("V1_Mouse_Brain_Sagittal_Anterior_Section_2.Rdata")

dim(X)[1]
length(unique(X$x_coords))
length(unique(X$y_coords))
length(unique(X$clus))

#########################################################

# SPARK

setwd("c:/Dropbox/Highd-CESME/tensor/simulation/SPARK-Analysis-master")

source("funcs/funcs.R")

# MERFISH
data = read.csv("processed_data/MERFISH_Animal18_Bregma0.11_countdata.csv")
info = read.csv("processed_data/MERFISH_Animal18_Bregma0.11_info.csv")

total_counts = apply(data[,-1],1,sum)
X = cbind.data.frame(total_counts = total_counts,x = info$x,y = info$y,clus = info$Cell_class)
save(X,file="spark_merfish.Rdata")

# seqFISH
data = read.csv("processed_data/seqFISH_field43_countdata.csv")
info = read.csv("processed_data/seqFISH_field43_info.csv")

total_counts = apply(data[,-1],1,sum)
X = cbind.data.frame(total_counts = total_counts,x = info$x,y = info$y,clus = info$z)
save(X,file="spark_seqFISH.Rdata")

# MOB
data = read.csv("processed_data/sim_MOB_pattern2_fc3_tau35_count_power1.csv")
info = read.csv("processed_data/Rep11_MOB_info_spark.csv")

X = info[,c(4,2,3)]
save(X,file="spark_MOB.Rdata")

# BC
setwd("c:/Dropbox/Highd-CESME/tensor/simulation/SPARK-master")
load("data/Layer2_BC_Count.rds")

info <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",1)),
                         y=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",2)),
                         total_counts=apply(rawcount,2,sum))
rownames(info) <- colnames(rawcount)

X = info
save(X,file="spark_BC.Rdata")

################################################
setwd("c:/Dropbox/Highd-CESME/tensor/simulation")

##########################################
# MERFISH
load(file = "spark_merfish.Rdata")

dim(X)[1]
length(unique(X$x))
length(unique(X$y))
length(unique(X$clus))


##########################################
# seqFISH
load(file = "spark_seqFISH.Rdata")

dim(X)[1]
length(unique(X$x))
length(unique(X$y))
length(unique(X$clus))

##########################################
# MOB
load(file = "spark_MOB.Rdata")

dim(X)[1]
length(unique(X$x))
length(unique(X$y))
length(unique(X$clus))


##########################################
# BC
load(file = "spark_BC.Rdata")

dim(X)[1]
length(unique(X$x))
length(unique(X$y))
length(unique(X$clus))



# Alleoscope

Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
BiocManager::install("rtracklayer")
devtools::install_github("seasoncloud/Alleloscope") # install
library(Alleloscope) # load

#memory.limit(4000)

setwd("C:/Users/lozhang/Documents/dropbox/tensor_data/Alleloscope-main") # set path to the github folder
setwd("C:/dropbox/dropbox/tensor_data/Alleloscope-main") # set path to the github folder

dir_path <- "c:/Dropbox/Highd-CESME/tensor/simulation/Alleoscope/output/"
dir_path <- "c:/dpbox/dropbox/Highd-CESME/tensor/simulation/Alleoscope/output/"
dir.create(dir_path) # set up output directory

data(centromere.GRCh38)
data(telomere.GRCh38)
size=read.table("data-raw/sizes.cellranger-GRCh38-1.0.0.txt", stringsAsFactors = F)

# SNP by cell matrices for ref and alt alleles
barcodes=read.table("data-raw/SNU601/scDNA/barcodes_sub.tsv", sep='\t', stringsAsFactors = F, header=F)
alt_all=readMM("data-raw/SNU601/scDNA/alt_all_sub.mtx")
ref_all=readMM("data-raw/SNU601/scDNA/ref_all_sub.mtx")
var_all=read.table("data-raw/SNU601/scDNA/var_all_sub.vcf", header = F, sep='\t', stringsAsFactors = F)

# bin by cell matrices for tumor and normal for segmentation
raw_counts=read.table("data-raw/SNU601/scDNA/tumor_sub.txt", sep='\t', header=T, row.names = 1,stringsAsFactors = F)
ref_counts=read.table("data-raw/SNU601/scDNA/normal_sub.txt", sep='\t', header=T, row.names = 1,stringsAsFactors = F) # Normal sample from patient 6198 was used for the cell line.

Obj=Createobj(alt_all =alt_all, ref_all = ref_all, var_all = var_all ,samplename='Sample', genome_assembly="GRCh38", dir_path=dir_path, barcodes=barcodes, size=size, assay='scDNAseq')

Obj_filtered=Matrix_filter(Obj=Obj, cell_filter=1000, SNP_filter=20, min_vaf = 0.1, max_vaf = 0.9, centro=centromere.GRCh38, telo=telomere.GRCh38) 

Obj_filtered=Segmentation(Obj_filtered=Obj_filtered, 
                          raw_counts=raw_counts, # from matched DNA sequencing (bulk/single)
                          ref_counts=ref_counts, # from matched DNA sequencing (bulk/single)
                          plot_seg = TRUE)

Obj_filtered=Segments_filter(Obj_filtered=Obj_filtered, nSNP=2000)

save(Obj_filtered,file = "tmp_result.rds")

load(file = "tmp_result.rds")

Obj_filtered=Est_regions(Obj_filtered = Obj_filtered, max_nSNP = 3000, plot_stat = T,cont = TRUE)

Obj_filtered=Select_normal(Obj_filtered = Obj_filtered, raw_counts=raw_counts, plot_theta = TRUE)

print(Obj_filtered$select_normal$region_normal)
Obj_filtered$ref=Obj_filtered$select_normal$region_normal[1] # choose one normal region

Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='tumor', raw_counts=raw_counts)  # for tumor
Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='cellline', raw_counts=raw_counts, ref_counts = ref_counts ) # for cell line without normal cells in the tumor sample.

Obj_filtered=Genotype(Obj_filtered = Obj_filtered)

linplot=Lineage_plot(Obj_filtered = Obj_filtered, nSNP = 2000,  nclust = 10 )

tmp = load("c:/Dropbox/Highd-CESME/tensor/simulation/Alleoscope/output/rds/seg_table_all_Sample.rds")
tmp = load("c:/Dropbox/Highd-CESME/tensor/simulation/Alleoscope/output/rds/EMresults/chr3")

# TAP-Seq

install.packages("RCurl",type = "binary")
install.packages("XML",type = "binary")

BiocManager::install("TAPseq")
library(TAPseq)

setwd("c:/Dropbox/Highd-CESME/tensor/simulation")

vignette("tapseq_primer_design", package = "TAPseq")

library(TAPseq)
library(GenomicRanges)
library(BiocParallel)

# gene annotations for chromosome 11 genomic region
data("chr11_genes")

# convert to GRangesList containing annotated exons per gene. for the sake of time we only use
# a small subset of the chr11 genes. use the full set for a more realistic example
target_genes <- split(chr11_genes, f = chr11_genes$gene_name)
target_genes <- target_genes[18:27]

# K562 Drop-seq read data (this is just a small example file within the R package)
dropseq_bam <- system.file("extdata", "chr11_k562_dropseq.bam", package = "TAPseq")

# register backend for parallelization
register(MulticoreParam(workers = 5))

# infer polyA sites from Drop-seq data
polyA_sites <- inferPolyASites(target_genes, bam = dropseq_bam, polyA_downstream = 50,
                               wdsize = 100, min_cvrg = 1, parallel = TRUE)




####################


table(dim.store[,1])
table(dim.store[,2])
table(dim.store[,3])


d1=readNIfTI("rest.nii.gz")
d1=readNIfTI("anat_rest.nii.gz")
dim(d1)
d1[,,,1]


# low dimensional approximation 

source("tcross.R")

p.all = rep(50,3)
A = as.tensor(array(rnorm(prod(p.all)),dim = p.all))

# random
m = 10
Omega = list()
for(i in 1:3){
  Omega[[i]] = sort(sample(p.all[i],m))
}

g = 10
Sigma = list()
for(i in 1:3){
  Sigma[[i]] = sort(sample(m^2,g))
}

# largest marginal variance

m = 10
Omega = list()
for(i in 1:3){
  A.tmp = k_unfold(A,m=i)@data
  s.tmp = apply(A.tmp,1,var)
  Omega[[i]] = sort(order(s.tmp,decreasing = TRUE)[1:m])
}

g = 10
Sigma = list()

A.tmp = k_unfold(A[,Omega[[2]],Omega[[3]]],m=1)@data
s.tmp = apply(A.tmp,2,var)
Sigma[[1]] = sort(order(s.tmp,decreasing = TRUE)[1:g])

A.tmp = k_unfold(A[Omega[[1]],,Omega[[3]]],m=2)@data
s.tmp = apply(A.tmp,2,var)
Sigma[[2]] = sort(order(s.tmp,decreasing = TRUE)[1:g])

A.tmp = k_unfold(A[Omega[[1]],Omega[[2]],],m=3)@data
s.tmp = apply(A.tmp,2,var)
Sigma[[3]] = sort(order(s.tmp,decreasing = TRUE)[1:g])


##############################
A.arm = list()
for(i in 1:3){
  A.arm[[i]] = k_unfold(A,m = i)@data[,Sigma[[i]]]
}

B = A[Omega[[1]],Omega[[2]],Omega[[3]]]

result = tcross_cv(B = B,Arm = A.arm)






