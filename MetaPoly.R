library(vcfR)
library(ape)
library(stringr)

################# GET DATA PER GENE ################# 
gff_to_gene_data <- function(gff, alleles){
  gene_id = substring(strsplit(gff$V9, ';')[[1]][1], 4)
  gene_contig = gff$V1
  gene_start = gff$V4
  gene_end = gff$V5
  return(list(gene_id = gene_id,
              gene_contig = gene_contig,
              gene_start = gene_start,
              gene_end = gene_end,
              gene_length = gene_end - gene_start,
              gene_data = alleles[(as.character(alleles$contig) == gene_contig) & (as.integer(alleles$pos) > gene_start) & (as.integer(alleles$pos) < gene_end),]))}

get_genes_data <- function(gff, alleles){
  genes = lapply(1:nrow(gff), function(i) gff_to_gene_data(gff[i,], alleles))
  return(genes)}


################# COMPUTE PI ################# 
is_polymorphic <- function(char){
  n_all = as.integer(unique(strsplit(char,',')[[1]]))
  n_all = n_all[n_all > 0]
  if (length(n_all) > 1){return(TRUE)}
  else {return(FALSE)}}

mean_depth <- function(char){
  return(mean(as.integer(strsplit(char,',')[[1]])))}

calc_pi_gene <- function(gene_data, sample, group, gene_id, gene_length){
  pi = sum(vapply(gene_data[,sample], function(x) is_polymorphic(x), logical(1)))
  depth = sum(vapply(gene_data[,sample], function(x) mean_depth(x), double(1)))
  pi_df = data.frame(pi=pi, depth=depth, sample=sample, group=group, gene_id=gene_id, gene_length=gene_length)
  return(pi_df)}

fit_cor_gene <- function(gene_df, gene, min_samp_per_group, group1, group2){
  gene_data = full_pi_df[full_pi_df$gene_id == gene,]
  gene_data$type[gene_data$type == group1] = '1'
  gene_data$type[gene_data$type == group2] = '0'
  gene_data$type = as.integer(as.numeric(gene_data$type))
  if ((sum(gene_data$type == 1) > min_samp_per_group) & (sum(gene_data$type == 0) > min_samp_per_group)){
    cor_test = cor.test(gene_data$res_m, gene_data$type)
    p = as.numeric(cor_test$p.value)
    cor = as.numeric(cor_test$estimate)
    mean_coef = mean(gene_data$res_m, na.rm=T)}
  else{p = NA
  cor = NA
  mean_coef = NA}
  return(data.frame(gene_id=gene, cor=cor, p=p, mean_coef=mean_coef))}

PiCorr <- function(alleles_genes, min_samp_per_group, samp_vec, group1, group2){
  print('Launching - MetaPoly PiCorr: a polymorphism-variable correlation tool for metagenomic data')
  t0 = Sys.time()
  
  print(' - Computing nuc. diversity (pi) per gene, per sample...')
  pi_iter = expand.grid(1:length(alleles_genes), samp_vec)
  pi_df = as.data.frame(do.call(rbind,lapply(1:nrow(pi_iter), function(i) calc_pi_gene(alleles_genes[[pi_iter$Var1[i]]]$gene_data, pi_iter$Var2[i], names(samp_vec)[samp_vec == pi_iter$Var2[i]], alleles_genes[[pi_iter$Var1[i]]]$gene_id, alleles_genes[[pi_iter$Var1[i]]]$gene_length))))
  print(pi_df)
  print(Sys.time() - t0)
  
  print(' - Fitting the poisson model on data...')
  model = glm(data = pi_df, family = poisson(), formula = pi ~ depth + gene_length + sample, control = list(maxit = 10))
  summary(model)
  pi_df$res_m = model$residuals
  print(Sys.time() - t0)
  
  # samples coefficients
  print(' - Gathering coefficients per samples...')
  coefs = coefficients(model)
  coefs = coefs[startsWith(names(coefs),'sample')]
  coefs_df = as.data.frame(coefs)
  rownames(coefs_df) = vapply(rownames(coefs_df), function(x) strsplit(x, 'sample')[[1]][2], FUN.VALUE = character(1))
  coefs_df$type = vapply(rownames(coefs_df), function(x) names(samp_vec)[samp_vec == x], character(1))
  print(Sys.time() - t0)
  
  print(' - Computing correlations of polymorphism with the variables of interest per gene...')
  corr_df = na.omit(do.call(rbind, lapply(unique(pi_df$gene_id), function(gene) fit_cor_gene(pi_df, gene, min_samp_per_group, group1, group2))))
  corr_df$padj = p.adjust(corr_df$p, method = 'holm')
  sign_genes = na.omit(corr_df[corr_df$padj < 0.05,])
  pos_genes = sign_genes$gene_id[sign_genes$cor > 0]
  neg_genes = sign_genes$gene_id[sign_genes$cor < 0]
  print(Sys.time() - t0)
  
  print(' - Analysis done!')
  return(list(pi_corr_res = corr_df, pos_genes = pos_genes, neg_genes = neg_genes, coefs = coefs_df))}