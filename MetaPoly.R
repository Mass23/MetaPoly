library(vcfR)
library(ape)
library(stringr)
library(data.table)

################# GET DATA PER GENE ################# 
gff_to_gene_data <- function(gff, gene_data){
  gene_id = substring(strsplit(gff$ATTRIBUTES, ';')[[1]][1], 4)
  gene_contig = gff$CONTIG
  gene_start = gff$START
  gene_end = gff$END
  
  cols = colnames(gene_data)[!(colnames(gene_data) %in% c('POS','CONTIG','rn'))]
  for (j in cols) set(gene_data, j = j, value = lapply(strsplit(gene_data[[j]], ','), as.integer))
  
  return(list(gene_id = gene_id,
              gene_contig = gene_contig,
              gene_start = gene_start,
              gene_end = gene_end,
              gene_length = gene_end - gene_start,
              gene_data = gene_data))}

GetGenesData <- function(gff, vcf){
  t0 = Sys.time()
  # Load gff and vcf, convert to data.table, get pos of the VCF variants
  print('Launching - MetaPoly GetGeneData: loading VCF and GFF data...')

  print('   - Formatting GFF file...')
  colnames(gff) = c('CONTIG','ANNOT','TYPE','START','END','SCORE','STRAND','PHASE','ATTRIBUTES')
  gff$GENE = vapply(gff$ATTRIBUTES, function(x) strsplit(strsplit(x, ';')[[1]][1], 'ID=')[[1]][2], FUN.VALUE = character(1))
  print(Sys.time()-t0)
  
  print('   - Extracting allele depths per site...')
  vcf <- extract.indels(vcf, return.indels = FALSE)
  samples = colnames(vcf@gt)[colnames(vcf@gt) != 'FORMAT']
  data <- as.data.table(extract.gt(vcf, element = 'AD'), keep.rownames=TRUE)
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

################# COMPUTE SNP density ################# 
GetDepth <- function(gene_data){
  return(colMeans(apply(gene_data, c(1,2), function(x) sum(unlist(x)))))}

GetSnpN <- function(gene_data){
  return(colSums(apply(gene_data, c(1,2), function(x) ifelse(length(x[[1]][x[[1]] > 0]), 1, 0))))}

fit_cor_gene <- function(gene_data, gene, min_samp_per_group, samp_vec){
  cor_test = cor.test(gene_data$res_m, gene_data$variable)
  p = as.numeric(cor_test$p.value)
  cor = as.numeric(cor_test$estimate)
  mean_coef = mean(gene_data$res_m, na.rm=T)
  low_ci = cor_test$conf.int[1]
  high_ci = cor_test$conf.int[2]
  return(data.frame(gene_id=gene, cor=cor, p=p, mean_coef=mean_coef, low_ci=low_ci, high_ci=high_ci))}


PiCorr <- function(data, min_samp_per_group, samp_vec){
  print('Launching - MetaPoly PiCorr: a polymorphism-variable correlation tool for metagenomic data')
  t0 = Sys.time()
  
  print(' - Computing SNP density and depth...')
  model_df = data.frame()
  cols = as.vector(samp_vec)
  count=0
  for (i in 1:length(data)){
      count = count +  1
      cat("\r",count)
      if (nrow(data[[i]]$gene_data) > 0){model_df = rbind(model_df , data.frame(gene_id = rep(data[[i]]$gene_id,length(samp_vec)),
                                                                                sample = as.vector(samp_vec),
                                                                                variable = as.integer(names(samp_vec)),
                                                                                snp_den = as.vector(GetSnpN(data[[i]]$gene_data[,..cols])),
                                                                                depth = as.vector(GetDepth(data[[i]]$gene_data[,..cols])),
                                                                                gene_length = rep(data[[i]]$gene_length,length(samp_vec))))}}
  cat(' genes done\n')
  print(model_df)
  print(Sys.time() - t0)
  
  print(' - Fitting the poisson model on data...')
  model = glm(data = model_df, family = poisson(), formula = snp_den ~ depth + gene_length + sample, control = list(maxit = 100))
  summary(model)
  model_df$res_m = model$residuals
  print(Sys.time() - t0)
  
  print(' - Computing sample coefficients...')
  coefs = coefficients(model)
  coefs = coefs[startsWith(names(coefs),'sample')]
  coefs_df = as.data.frame(coefs)
  rownames(coefs_df) = vapply(rownames(coefs_df), function(x) strsplit(x, 'sample')[[1]][2], FUN.VALUE = character(1))
  coefs_df$type = vapply(rownames(coefs_df), function(x) names(samp_vec)[samp_vec == x], character(1))
  print(Sys.time() - t0)
  
  print(' - Computing correlations of polymorphism with the variables of interest per gene...')
  corr_df = do.call(rbind, lapply(unique(model_df$gene_id), function(gene) fit_cor_gene(model_df[model_df$gene_id == gene,], gene, min_samp_per_group, samp_vec)))
  corr_df$padj = p.adjust(corr_df$p, method = 'holm')
  sign_genes = corr_df[corr_df$padj < 0.05,]
  pos_genes = as.vector(na.omit(sign_genes$gene_id[sign_genes$cor > 0]))
  neg_genes = as.vector(na.omit(sign_genes$gene_id[sign_genes$cor < 0]))
  print(Sys.time() - t0)
  
  print(' - Analysis done!')
  return(list(pi_corr_res = corr_df, pos_genes = pos_genes, neg_genes = neg_genes, coefs = coefs_df))}



