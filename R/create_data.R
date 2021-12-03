################# 1. GET DATA PER GENE ################# 
gff_to_gene_data <- function(gff, gene_data){
  gene_id = substring(strsplit(gff$ATTRIBUTES, ';')[[1]][1], 4)
  gene_contig = gff$CONTIG
  gene_start = gff$START
  gene_end = gff$END
  gene_strand = gff$STRAND
  
  cols = colnames(gene_data)[!(colnames(gene_data) %in% c('POS','CONTIG','rn'))]
  for (j in cols) set(gene_data, j = j, value = as.vector(lapply(strsplit(gene_data[[j]], ','), as.integer)))
  
  return(list(gene_id = gene_id,
              gene_contig = gene_contig,
              gene_strand = gene_strand,
              gene_start = gene_start,
              gene_end = gene_end,
              gene_length = gene_end - gene_start + 1,
              gene_data = gene_data))}

#' GetGenesData
#'
#' Metapoly's import function, takes as input a gff file open via read.delim (see
#' example on github), and a VCF file loaded with vcfR (read.vcfR). Returns the 
#' polymorphism data in the format used for all MetaPoly's analyses.
#'
#' @export
GetGenesData <- function(gff, vcf){
  t0 = Sys.time()
  # Load gff and vcf, convert to data.table, get pos of the VCF variants
  print('Launching - MetaPoly GetGenesData: loading VCF and GFF data...')

  print('   - Formatting GFF file...')
  colnames(gff) = c('CONTIG','ANNOT','TYPE','START','END','SCORE','STRAND','PHASE','ATTRIBUTES')
  gff$GENE = vapply(gff$ATTRIBUTES, function(x) strsplit(strsplit(x, ';')[[1]][1], 'ID=')[[1]][2], FUN.VALUE = character(1))
  print(Sys.time()-t0)
  
  print('   - Extracting allele depths per site...')
  vcf <- vcfR::extract.indels(vcf, return.indels = FALSE)
  samples = colnames(vcf@gt)[colnames(vcf@gt) != 'FORMAT']
  data <- data.table::as.data.table(vcfR::extract.gt(vcf, element = 'AD'), keep.rownames=TRUE)
  data[,POS := vapply(data[['rn']], function(x) as.integer(strsplit(x,'_')[[1]][3]), FUN.VALUE = integer(1))]
  data[,CONTIG := vapply(data[['rn']], function(x) sub("_[^_]+$", "",x), FUN.VALUE = character(1))]

  print(Sys.time()-t0)
  
  # Create data per gene
  print('   - Gathering data per gene...')
  n_genes <- nrow(gff)
  genes_data = lapply(1:n_genes, function(i) gff_to_gene_data(gff[i,], data[which((CONTIG==gff$CONTIG[i]) & (POS > gff$START[i]) & (POS < gff$END[i]))]))
  print(Sys.time()-t0)
  print(' - Data loaded!')
  return(genes_data)}

############### 2. POLYMORPHISM SUMMARY ################ 
CalcPi <- function(ac){return(1-sum((ac/sum(ac))**2))}
CalcMAJF <- function(ac){return(max(ac / sum(ac)))}
  
#if(sum(ac > 0) > 1){
#     af = ac / sum(ac)
#     
#     #return(mean(sample(1:length(af), 10000, prob = af, replace = T) != sample(1:length(af), 10000, prob = af, replace = T)))}
# else{return(0)}}

GetSnpData <- function(gene_data){
  depth = colMeans(apply(gene_data, c(1,2), function(ac) sum(unlist(ac))))
  snp_n = colSums(apply(gene_data, c(1,2), function(ac) ifelse(length(ac[[1]][ac[[1]] > 0]) > 1, 1, 0)))
  majf = colMeans(apply(gene_data, c(1,2), function(ac) CalcMAJF(as.matrix(ac)[1][[1]])), na.rm = T)
  npi = colMeans(apply(gene_data, c(1,2), function(ac) CalcPi(as.matrix(ac)[1][[1]])), na.rm = T)
  return(list(depth=depth,snp_n=snp_n,majf=majf,npi=npi))}

SumSnpData <- function(data, samp_vec){
  if (nrow(data$gene_data) > 0){snp_data = GetSnpData(data$gene_data[,..cols])
                                       out_df = data.frame(gene_id = rep(data$gene_id,length(samp_vec)),
                                                sample = as.vector(samp_vec),
                                                variable = as.numeric(names(samp_vec)),
                                                SNP_N = snp_data$snp_n,
                                                DEPTH = snp_data$depth,
                                                MAJF = snp_data$majf,
                                                PI = snp_data$npi,
                                                gene_length = rep(data$gene_length,length(samp_vec)))
                                      return(out_df)}
  else{snp_data = GetSnpData(data$gene_data[,..cols])
       out_df = data.frame(gene_id = rep(data$gene_id,length(samp_vec)),
                           sample = as.vector(samp_vec),
                           variable = as.numeric(names(samp_vec)),
                           SNP_N = NA,
                           DEPTH = NA,
                           MAJF = NA,
                           PI = NA,
                           gene_length = rep(data$gene_length,length(samp_vec)))
  return(out_df)}}

#' PolySummary
#'
#' Summarises the polymorphism in the different genes, across samples. The number
#' of SNP (SNP_N), sequencing depth (DEPTH), major allele frequency (MAJF), and 
#' allele frequencies evenness (EVENNESS) are measured and returned in a large df.
#'
#' @export
PolySummary <- function(data, samp_vec, n_cores){
  print('Launching - MetaPoly PolySummary: summarising the polymorphism data...')
  t0 = Sys.time()
  
  print(' - Computing metrics')
  PolyDf = data.frame()
  registerDoMC(n_cores)
  PolyDf = foreach(i=1:length(data), .combine = rbind) %dopar% SumSnpData(data[[i]], samp_vec)
    
  cat(' genes done\n')
  PolyDf$Cons_index = ((PolyDf$gene_length - PolyDf$SNP_N)/PolyDf$gene_length) + ((PolyDf$SNP_N/PolyDf$gene_length)*PolyDf$MAJF)
  print(Sys.time()-t0)
  return(PolyDf)}

#' SummariseSamples
#'
#' Summarises the polymorphism data generated by PolySummary by sample, to get genome-
#' wide data per sample
#'
#' @export
SummariseSamples <- function(poly_summary, val_depth){
  genes_to_keep = c()
  for (gene in unique(poly_summary$gene_id)){
    if (min(poly_summary$DEPTH[poly_summary$gene_id == gene]) > val_depth){
      genes_to_keep = c(genes_to_keep, gene)}}
  poly_summary = poly_summary[poly_summary$gene_id %in% genes_to_keep,]
  sample_df = data.frame()
  for (sample in unique(poly_summary$sample)){
    mean_snp_den =  weighted.mean(poly_summary$SNP_N[poly_summary$sample == sample] / poly_summary$gene_length[poly_summary$sample == sample], poly_summary$gene_length[poly_summary$sample == sample], na.rm=T)
    mean_pi =  weighted.mean(poly_summary$PI[poly_summary$sample == sample], poly_summary$gene_length[poly_summary$sample == sample], na.rm=T)
    mean_depth =  weighted.mean(poly_summary$DEPTH[poly_summary$sample == sample], poly_summary$gene_length[poly_summary$sample == sample], na.rm=T)
    mean_majf =  weighted.mean(poly_summary$MAJF[poly_summary$sample == sample], poly_summary$gene_length[poly_summary$sample == sample], na.rm=T)
    mean_cons =  weighted.mean(poly_summary$Cons_index[poly_summary$sample == sample], poly_summary$gene_length[poly_summary$sample == sample], na.rm=T)
    sample_df = rbind(sample_df, data.frame(sample=sample,MEAN_PI=mean_pi,MEAN_CONS_I=mean_cons,MEAN_SNP_DEN=mean_snp_den,MEAN_DEPTH=mean_depth,MEAN_MAJF=mean_majf))}
  return(list(table=sample_df,n_genes=length(genes_to_keep)))}