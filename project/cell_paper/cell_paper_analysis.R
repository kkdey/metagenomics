
## Testing the CountClust package on the new metagenomics data

library(fossil) # for Chao1
library(vegan)     # for Simpson
library(metagenomeSeq) # data structured
library(breakaway) # for breakaway
library(devtools)
#install_github('kkdey/CountClust')
library(CountClust)
library(maptpx)

setwd('/Users/kushal/Documents/metagenomics/project/cell_paper/')

load("cellPaperRaw.rdata")
# data is of the MRexperiment format:
# to grab the raw counts it's just MRcounts(cellPaper)
# this is of course the raw OTU data...

pheno_data <- pData(cellPaper);
age = pData(cellPaper)$Age_at_Collection
status = pData(cellPaper)$Case_Control
gender =pData(cellPaper)$gender;
counts = MRcounts(cellPaper)

h = diversity(t(counts))  ## determining the species diversity


## determining the species richness (Chao 1984,1987)

c1= sapply(1:ncol(cellPaper),function(i){
  chao1(counts[,i])
})

names(c1) = names(h) = names(status) = names(age) = colnames(cellPaper)=names(gender)

plot(c1~age,col=status,bg=status,pch=21,bty="l",ylab="Chao1",xlab="Age at collection")

###  applying the CountClust functions 

counts <- t(counts);

## handling NA if any in the data
counts_preprocess <- handleNA(counts)$data;

## Removing the sparse OTUs
counts_filtered <- RemoveSparseFeatures(counts_preprocess,0.99)$data; 

## the OTUs that are very sparse 
otu_removed <- RemoveSparseFeatures(counts_preprocess,0.99)$sparse_features;

## Performing the Structure plots for 3 choices (k=2,3,4) with Age, gender and status as metadata 

nclus_vec <- 5:6;

if(!dir.exists("Structure")) dir.create("Structure")
if(!dir.exists("Structure/batch_uncorrected")) dir.create("Structure/batch_uncorrected")

bayesfac <- array(0,length(nclus_vec));

samp_metadata <- cbind.data.frame(age,status,gender);
colnames(samp_metadata) = c("Age","Status","Gender");

for(num in 1:length(nclus_vec))
{
  if(!dir.exists(paste0("Structure/batch_uncorrected/clus_",nclus_vec[num]))) dir.create(paste0("Structure/batch_uncorrected/clus_",nclus_vec[num]))
  obj <- StructureObj(counts_filtered,nclus_vec[num],samp_metadata = samp_metadata, tol=0.001, batch_lab = NULL, path=paste0("Structure/batch_uncorrected/clus_",nclus_vec[num]),partition=c('FALSE','TRUE','TRUE'));
}




