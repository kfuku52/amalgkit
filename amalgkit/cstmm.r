suppressPackageStartupMessages(library(edgeR, quietly=TRUE))

mode = ifelse(length(commandArgs(trailingOnly=TRUE))==1, 'debug', 'batch')

if (mode=="debug") {
  dir_work = '/Users/s229181/MSN/'
  dir_ortho = paste0(dir_work, "OrthoFinder/Results_Jun22/Orthogroups")
  dir_count = paste0(dir_work, "counts/")
  mode_tmm = 'multi_species'
  #dir_work = '/Users/kef74yk/Dropbox_p/collaborators/Ken Naito/20210509_Vigna/gfe_data'
  #dir_ortho = "/Users/kef74yk/Dropbox_p/collaborators/Ken Naito/20210509_Vigna/gfe_data/Orthogroups"
  #dir_count = "/Users/kef74yk/Dropbox_p/collaborators/Ken Naito/20210509_Vigna/gfe_data/merge"
  setwd(dir_work)
} else if (mode=="batch") {
  args = commandArgs(trailingOnly=TRUE)
  dir_count = args[1]
  dir_ortho = args[2]
  dir_work = args[3]
  mode_tmm = args[4]
}

get_spp_filled = function(dir_count, df_gc=NA) {
  sciname_dirs = list.files(dir_count)
  spp_filled = c()
  for (sciname_dir in sciname_dirs) {
    count_files = list.files(path = file.path(dir_count, sciname_dir), pattern = ".*est_counts\\.tsv")
    if (length(count_files)==1) {
      spp_filled = c(spp_filled, count_files)
    }
  }
  spp_filled = sub('_', '|', spp_filled)
  spp_filled = sub('_.*', '', spp_filled)
  spp_filled = sub('\\|', '_', spp_filled)
  if ('data.frame' %in% class(df_gc)) {
    is_missing_in_genecount = (!spp_filled %in% colnames(df_gc))
    if (sum(is_missing_in_genecount)) {
      for (sp in spp_filled[is_missing_in_genecount]) {
        warning(paste0('Species excluded. Not found in OrthoFinder\'s GeneCount table: ', sp))
      }
    }
    spp_filled = spp_filled[!is_missing_in_genecount]
  }
  cat('Detected species:', spp_filled, '\n')
  return(spp_filled)
}

get_singlecopy_og = function(df_gc, spp_filled) {
  is_singlecopy = TRUE
  for (sp in spp_filled) {
    is_singlecopy = is_singlecopy & (df_gc[,sp]==1)
  }
  sc_og = df_gc[is_singlecopy,'Orthogroup']
  cat(length(sc_og), 'single-copy orthogroups were detected for the', length(spp_filled), 'species.\n')
  return(sc_og)
}

read_est_counts = function(dir_count, sp) {
  sciname_path = file.path(dir_count, sp)
  infile = list.files(path=sciname_path, pattern=".*est_counts\\.tsv")
  if (length(infile)> 1){
    stop(paste0("Multiple *count.tsv files found for: ", sp ,"\n"))
  } else if (length(infile)==0) {
    warning(paste0("Skipping. No *est_counts.tsv files found for: ", sp ,"\n"))
    return(NULL)
  }
  infile_path = file.path(sciname_path, infile[1])
  cat('Input file found, reading:', infile[1], '\n')
  dat = read.delim(infile_path, header = T, row.names=1, sep='\t', check.names=FALSE)
  dat = dat[,(colnames(dat)!='length')]
  colnames(dat) = paste(sp, colnames(dat), sep='_')
  return(dat)
}

get_uncorrected = function(dir_count, file_genecount=NA) {
  if (is.na(file_genecount)) {
    df_gc = NA
  } else {
    df_gc = read.table(file_genecount, header=TRUE, sep='\t', check.names=FALSE)
  }
  spp_filled = get_spp_filled(dir_count, df_gc)
  uncorrected = list()
  for (sp in spp_filled) {
    dat = read_est_counts(dir_count, sp)
    if (is.null(dat)) {
      next
    }
    uncorrected[[sp]] = dat
  }
  return(uncorrected)
}

get_df_sog = function(file_genecount, file_orthogroup, dir_count, uncorrected) {
  df_gc = read.table(file_genecount, header=TRUE, sep='\t', check.names=FALSE)
  spp_filled = get_spp_filled(dir_count, df_gc)
  single_orthogroups = get_singlecopy_og(df_gc, spp_filled)
  df_og = read.table(file_orthogroup, header=TRUE, sep='\t', row.names=1, check.names=FALSE)
  df_singleog = df_og[(rownames(df_og) %in% single_orthogroups), spp_filled, drop=FALSE]
  df_sog = df_singleog
  for (sp in spp_filled) {
    if (!sp %in% names(uncorrected)) {
      next
    }
    df_sog = merge(df_sog, uncorrected[[sp]], by.x=sp, by.y="row.names", all.x=TRUE, all.y=FALSE, sort=FALSE)
  }
  df_sog = df_sog[,-(1:length(spp_filled))]
  rownames(df_sog) = rownames(df_singleog)
  return(df_sog)
}

get_df_nonzero = function(df_sog) {
  is_na_containing_row = apply(df_sog, 1, function(x){any(is.na(x))})
  cat('Removing', sum(is_na_containing_row), 'out of', nrow(df_sog), 'orthogroups because missing values are observed in at least 1 species.\n')
  is_no_count_col = apply(df_sog, 2, function(x){sum(x, na.rm=TRUE)==0})
  cat('Removing', sum(is_no_count_col), 'out of', ncol(df_sog), 'RNA-seq samples because read mapping values are all zero.\n')
  df_nonzero = df_sog[!is_na_containing_row,!is_no_count_col]
  return(df_nonzero)
}

create_eff_length_symlink = function(dir_count, dir_cstmm, sp) {
  path_sp = file.path(dir_count, sp)
  eff_length_files = list.files(path = path_sp, pattern = ".*eff_length\\.tsv")
  if (length(eff_length_files)==1) {
    path_target = file.path(path_sp, eff_length_files[1])
    path_link = file.path(dir_cstmm, sp, eff_length_files[1])
    if (file.exists(path_link)) {
      file.remove(path_link)
    }
    file.symlink(from=path_target, to=path_link)
  } else {
    warning(paste0('No eff_length.tsv file found: ', path_sp))
  }
}

# set directory
if (!file.exists(dir_work)) {
  dir.create(dir_work)
}
setwd(dir_work)

dir_cstmm = file.path(dir_work,'cstmm')
if (!file.exists(dir_cstmm)) {
  dir.create(dir_cstmm)
}

file_genecount = file.path(dir_ortho, 'Orthogroups.GeneCount.tsv')
file_orthogroup = file.path(dir_ortho, 'Orthogroups.tsv')

if (mode_tmm=='single_species') {
  uncorrected = get_uncorrected(dir_count=dir_count, file_genecount=NA)
  stopifnot(length(names(uncorrected))==1)
  sp = names(uncorrected)[[1]]
  df_sog = uncorrected[[sp]]
  df_nonzero = get_df_nonzero(df_sog)
} else if (mode_tmm=='multi_species') {
  uncorrected = get_uncorrected(dir_count=dir_count, file_genecount=file_genecount)
  df_sog = get_df_sog(file_genecount, file_orthogroup, dir_count, uncorrected)
  df_nonzero = get_df_nonzero(df_sog)
}

cnf_in = edgeR::DGEList(counts = df_nonzero)

cat('Round 1: Performing TMM normalization to determine the appropriate baseline.\n')
cnf_out1 = edgeR::calcNormFactors(cnf_in, method='TMM', refColumn=NULL)
x = cnf_out1[[2]][['norm.factors']]
cat('Round 1: Median TMM normalization factor =', median(x), '\n')
median_value = sort(x)[ceiling(length(x)/2)]
median_index = (1:length(x))[x==median_value]

cat('Round 2: Performing TMM normalization for output.\n')
cnf_out2 = edgeR::calcNormFactors(cnf_in, method='TMM', refColumn=median_index)
cat('Round 2: Median TMM normalization factor =', median(cnf_out2[[2]][['norm.factors']]), '\n')

df_nf = cnf_out2[[2]]
df_nf[['sample']] = rownames(df_nf)
df_nf = df_nf[,c('sample','group','lib.size','norm.factors')]
write.table(df_nf, 'cstmm/normalization_factor.tsv', row.names=FALSE, sep='\t')

xlim = c(-2,2)
bins = seq(-2,2,0.1)
file_name='cstmm/normalization_factor_histogram.pdf'
pdf(file_name, height=3.3, width=7.2) # full figure size = 9.7 x 7.2
x = log2(cnf_out2[[2]][['norm.factors']])
x[x>xlim[2]] = xlim[2]
x[x<xlim[1]] = xlim[1]
hist(x, xlab='log2(TMM normalization factor)', ylab='Count', main=NULL, col='black', xlim=xlim, breaks=bins, las=1)
graphics.off()

for (sp in names(uncorrected)) {
  cat('Applying TMM normalization factors:', sp, '\n')
  dat = uncorrected[[sp]]
  df_nf_sp = cnf_out2[[2]][startsWith(rownames(cnf_out2[[2]]),sp),]

  for (i in 1:length(df_nf_sp[,1])){
    # check if dat colnames start with the species name
    # and modify SRR variable to match the format of dat colnames
    # TODO: is this necessary? Is there any situation where the species name prefix is missing?
    if(all(startsWith(colnames(dat), sp))){
      SRR = as.character(row.names(df_nf_sp[i,]))
    }
    else{
      SRR = as.character(strsplit(row.names(df_nf_sp[i,]),'_')[[1]][3])
    }
    # manually apply normfactor
    tmm_normalization_factor = as.double(df_nf_sp[i,"norm.factors"])
    dat[,SRR] = dat[,SRR]/tmm_normalization_factor

  }
  dat_out = cbind(target_id=rownames(dat), dat)
  rownames(dat_out) = NULL
  colnames(dat_out) = sub(paste0(sp, '_'), '', colnames(dat_out))
  dir_cstmm_sp = file.path(dir_cstmm, sp)
  if (!file.exists(dir_cstmm_sp)) {
    dir.create(dir_cstmm_sp)
  }
  create_eff_length_symlink(dir_count, dir_cstmm, sp)
  file_path = file.path(dir_cstmm_sp, paste0(sp, "_cstmm_counts.tsv"))
  write.table(dat_out, file_path, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)

}
cat('cstmm.r completed!\n')