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
CalcMAJF <- function(ac){
  if(unique(ac) == 0){return(NA)}
  af = ac / sum(ac)
  majf = max(af)
  if (majf == 1){majf = NA}
  return(majf)}

CalcEven <- function(ac){
  af = ac / sum(ac)
  evenness = sum(-vapply(af, function(x) x * log(x), FUN.VALUE = numeric(1)))/log(length(af))
  return(evenness)}

CalcPi <- function(ac){
  af = ac / sum(ac)
  pi = mean(sample(1:length(af), 1000, prob = af, replace = T) != sample(1:length(af), 1000, prob = af, replace = T))
  return(pi)}

GetSnpData <- function(gene_data){
  depth = colMeans(apply(gene_data, c(1,2), function(ac) sum(unlist(ac))))
  snp_n = colSums(apply(gene_data, c(1,2), function(ac) ifelse(length(ac[[1]][ac[[1]] > 0]) > 1, 1, 0)))
  evenness = colMeans(apply(gene_data, c(1,2), function(ac) CalcEven(as.matrix(ac)[1][[1]])), na.rm = T)
  majf = colMeans(apply(gene_data, c(1,2), function(ac) CalcMAJF(as.matrix(ac)[1][[1]])), na.rm = T)
  pi = colMeans(apply(gene_data, c(1,2), function(ac) CalcPi(as.matrix(ac)[1][[1]])), na.rm = T)
  return(list(depth=depth,snp_n=snp_n,evenness=evenness,majf=majf))}

#' PolySummary
#'
#' Summarises the polymorphism in the different genes, across samples. The number
#' of SNP (SNP_N), sequencing depth (DEPTH), major allele frequency (MAJF), and 
#' allele frequencies evenness (EVENNESS) are measured and returned in a large df.
#'
#' @export
PolySummary <- function(data, samp_vec){
  print('Launching - MetaPoly PolySummary: summarising the polymorphism data...')
  t0 = Sys.time()
  
  print(' - Computing metrics')
  PolyDf = data.frame()
  cols = as.vector(samp_vec)
  count=0
  for (i in 1:length(data)){
    count = count +  1
    cat("\r",count)
    if (nrow(data[[i]]$gene_data) > 0){snp_data = GetSnpData(data[[i]]$gene_data[,..cols])
                                       PolyDf = rbind(PolyDf , data.frame(gene_id = rep(data[[i]]$gene_id,length(samp_vec)),
                                                                          sample = as.vector(samp_vec),
                                                                          variable = as.numeric(names(samp_vec)),
                                                                          SNP_N = snp_data$snp_n,
                                                                          DEPTH = snp_data$depth,
                                                                          MAJF = snp_data$majf,
                                                                          EVENNESS = snp_data$evenness,
                                                                          PI = snp_data$pi,
                                                                          gene_length = rep(data[[i]]$gene_length,length(samp_vec))))}}
  cat(' genes done\n')
  print(Sys.time()-t0)
  return(PolyDf)}

