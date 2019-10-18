
if (length(commandArgs(trailingOnly=TRUE))==1) {
    mode = "debug"    
} else {
    mode = "batch"
}

if (mode=="debug") {
    infile = '/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/kallisto_summary/tpm.masked.kallisto.tsv/Canis_lupus.tpm.masked.kallisto.tsv'
    outfile = '/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/kallisto_summary/tpm.masked.kallisto.gene.log.tsv/Canis_lupus.tpm.masked.kallisto.gene.log.tsv'
    idfile = '/Users/kf/Dropbox/kfdata/02_Data/my_db/Ensembl/release-91/id_mapping/clupus_gene_ensembl.tsv'
    transcript2gene = 1
    log_transf = 1
} else if (mode=="batch") {
    args = commandArgs(trailingOnly=TRUE)
    infile = args[1]
    outfile = args[2]
    idfile = args[3]
    transcript2gene = as.integer(args[4])
    log_transf  = as.integer(args[5])
}

read_transcriptome = function(infile) {
    tc = read.table(infile, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    rownames(tc) = tc$target_id
    tc = data.frame(tc[,3:ncol(tc)], stringsAsFactors=FALSE)
    return(tc)
}

transcript2gene_transformation = function(tc, idfile) {
    sp_mapping = read.table(idfile, quote="", sep="\t", header=TRUE, stringsAsFactors=FALSE)
    tc$ensembl_transcript_id = sub("\\..*", "", rownames(tc))
    tc = merge(tc, sp_mapping[,c("Gene.stable.ID","Transcript.stable.ID")], by.x="ensembl_transcript_id", by.y="Transcript.stable.ID", sort=FALSE)
    tc = tc[,-1]
    tc = aggregate(x=tc[,-ncol(tc)], by=list(tc$Gene.stable.ID), FUN=sum)
    rownames(tc) = tc$Group.1
    tc = tc[,-1]
    return(tc)
}

transcriptome_transformation = function(infile, idfile, log_transf, transcript2gene) {
    tc = read_transcriptome(infile)
    if (transcript2gene) {
        cat('Transcript IDs are converted to gene IDs. Total expression value per gene are returned.\n')
        tc = transcript2gene_transformation(tc, idfile)
    }
    if (log_transf) {
        cat('Data are log-transformed.\n')
        tc = log2(tc+1)
    }
    is_allna = apply(is.na(tc), 2, all)
    if (any(is_allna)) {
        cat('Removed due to all NA:', '\n')
        print(colnames(tc)[is_allna])
        tc = tc[, !is_allna]
    }
    if (any(is.na(tc))) {
        cat('NA was detected in the output table:', '\n')
    }
    return(tc)
}


tc = transcriptome_transformation(infile, idfile, log_transf, transcript2gene)
write.table(format(tc, nsmall=5), file=outfile, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
