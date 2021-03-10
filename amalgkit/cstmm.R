library(limma)
library(tidyr)
library(dplyr)

if (length(commandArgs(trailingOnly=TRUE))==0) {
  mode = "debug"
} else {
  mode = "batch"
}

if (mode=="debug") {

  dir_work = '/Users/s229181/MSN/'
  dir_ortho = paste0(dir_work, "OrthoFinder/Results_Feb09_2/Orthogroups")
  dir_count = paste0(dir_work, "counts/")
  
  file_singlecopy = paste0(dir_ortho, '/Orthogroups_SingleCopyOrthologues.txt')
  file_orthogroup = paste0(dir_ortho, '/Orthogroups.tsv')
  single_orthogroups = read.table(file_singlecopy, header=FALSE)$V1
  df_og = read.table(file_orthogroup, header=TRUE, sep='\t', row.names=1)
  df_singleog = df_og[(rownames(df_og) %in% single_orthogroups),]
  
  
} else if (mode=="batch") {
  args = commandArgs(trailingOnly=TRUE)
  dir_count = args[1]
  dir_orthofinder_results = args[2]
  dir_work = args[3]
  
  file_singlecopy = paste0(dir_orthofinder_results, '/Orthogroups_SingleCopyOrthologues.txt')
  file_orthogroup = paste0(dir_orthofinder_results, '/Orthogroups.tsv')
  single_orthogroups = read.table(file_singlecopy, header=FALSE)$V1
  df_og = read.table(file_orthogroup, header=TRUE, sep='\t', row.names=1)
  df_singleog = df_og[(rownames(df_og) %in% single_orthogroups),]
  
}

if (!endsWith(dir_count, '/')) {
  dir_count = paste0(dir_count, '/')
}

if (!endsWith(dir_work, '/')) {
  dir_work = paste0(dir_work, '/')
}

# set directory
if (!file.exists(dir_work)) {
  dir.create(dir_work)
}
setwd(dir_work)


# set output directory and create if not already there
dir_tmm = paste0(dir_work, 'cross_species_tmm_normalized_counts/')
if (!file.exists(dir_tmm)) {

  dir.create(dir_tmm)
}




 # for (col in 1:ncol(df_singleog)) {
 #   df_singleog[,col] = sub('.*_', '', df_singleog[,col])
 # }
df_og = NULL
spp_filled = colnames(df_singleog)
spp = sub('_', ' ', spp_filled)
num_selected_spp = length(spp)
cat('num selected species:', num_selected_spp, '\n')
cat('selected species:', spp, '\n')

uncorrected = list()
df_sog = df_singleog

for (sp in spp_filled) {
  
  infile = list.files(path = dir_count, pattern=paste0(sp,".*count.\\.tsv"))
  if (length(infile)> 1){
    cat("Multiple files found for: ", sp ,"\n")
    cat("Trying: ",infile[1], ".\n")
  }
  infile_path = paste0(dir_count, infile[1])
  if (file.exists(infile_path)) {
    cat('Input file found, reading:', infile[1], '\n')
    dat = read.delim(infile_path,header = T,row.names=1, sep='\t')
    dat = dat[,(colnames(dat)!='length')]
    genus_cap = tolower(substring(sp,1,1))
    spec_name = strsplit(sp,'_')[[1]][2]
    colnames(dat) = paste(sp, colnames(dat), sep='_')
    uncorrected[[sp]] = dat
    df_sog = merge(df_sog, dat, by.x=sp, by.y="row.names", all.x=TRUE, all.y=FALSE, sort=FALSE)
  } else {
    cat('Input file not found:', infile, '\n')
  }
}


df_sog = df_sog[,-(1:num_selected_spp)]
#df_sog = apply(df_sog, 2, as.integer)
rownames(df_sog) = rownames(df_singleog)

df_nonzero = df_sog[,apply(df_sog, 2, sum)!=0]
cnf_in = edgeR::DGEList(counts = df_nonzero)

cnf_out1 = edgeR::calcNormFactors(cnf_in, method='TMM', refColumn=NULL)
x = cnf_out1[[2]][['norm.factors']]
cat('Round 1: Median normalization factor =', median(x), '\n')
median_value = sort(x)[ceiling(length(x)/2)]
median_index = (1:length(x))[x==median_value]

cnf_out2 = edgeR::calcNormFactors(cnf_in, method='TMM', refColumn=median_index)
cat('Round 2: Median normalization factor =', median(cnf_out2[[2]][['norm.factors']]), '\n')

df_nf = cnf_out2[[2]]
df_nf[['sample']] = rownames(df_nf)
df_nf = df_nf[,c('sample','group','lib.size','norm.factors')]
write.table(df_nf, 'normalization_factor.tsv', row.names=FALSE, sep='\t')

xlim = c(-2,2)
bins = seq(-2,2,0.1)
file_name='normalization_factor_histogram.pdf'
pdf(file_name, height=3.3, width=7.2) # full figure size = 9.7 x 7.2
x = log2(cnf_out2[[2]][['norm.factors']])
x[x>xlim[2]] = xlim[2]
x[x<xlim[1]] = xlim[1]
hist(x, xlab='log2(TMM normalization factor)', ylab='Count', main=NULL, col='black', xlim=xlim, breaks=bins, las=1)
graphics.off()

for (sp in names(uncorrected)) {
  infile = list.files(path = dir_count, pattern=paste0(sp,".*count.\\.tsv"))
  if (length(infile) > 1){
    cat("Multiple files found for: ", sp , "\n")
    cat("Trying: ",infile[1], ".\n")
  }
  infile_path = paste0(dir_count, infile[1])
  if (file.exists(infile_path)) {
    cat('Input file found, reading:', infile[1], '\n')
    dat = read.delim(infile_path, header=TRUE, row.names=1, sep='\t')
    dat = dat[,(colnames(dat)!='length')]
    df_nf_sp = cnf_out2[[2]][startsWith(rownames(cnf_out2[[2]]),sp),]
    
    for (i in 1:length(df_nf_sp[,1])){
        # check if dat colnames start with the species name
        # and modify SRR variable to match the format of dat colnames
      
        if(all(startsWith(colnames(dat), sp))){
          SRR = as.character(row.names(df_nf_sp[i]))
        }
        else{
          SRR = as.character(strsplit(row.names(df_nf_sp[i,]),'_')[[1]][3])
        }
        # manually apply normfactor
          dat[,SRR]<-dat[,SRR]/as.double(df_nf_sp[i,"norm.factors"])
          
        }
  
   write.table(dat,paste0(dir_tmm,sp,"_cstmm_counts.tsv") , sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
  }


}
  
  