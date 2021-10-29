cog_functions = c('J'='Translation, ribosomal structure and biogenesis',
                  'A'='RNA processing and modification',
                  'K'='Transcription',
                  'L'='Replication, recombination and repair',
                  'B'='Chromatin structure and dynamics',
                  'D'='Cell cycle control, cell division, chromosome partitioning',
                  'Y'='Nuclear structure',
                  'V'='Defense mechanisms',
                  'T'='Signal transduction mechanisms',
                  'M'='Cell wall/membrane/envelope biogenesis',
                  'N'='Cell motility',
                  'Z'='Cytoskeleton',
                  'W'='Extracellular structures',
                  'U'='Intracellular trafficking, secretion, and vesicular transport',
                  'O'='Posttranslational modification, protein turnover, chaperones',
                  'X'='Mobilome: prophages, transposons',
                  'C'='Energy production and conversion',
                  'G'='Carbohydrate transport and metabolism',
                  'E'='Amino acid transport and metabolism',
                  'F'='Nucleotide transport and metabolism',
                  'H'='Coenzyme transport and metabolism',
                  'I'='Lipid transport and metabolism',
                  'P'='Inorganic ion transport and metabolism',                     
                  'Q'='Secondary metabolites biosynthesis, transport and catabolism',
                  'R'='General function prediction only',
                  'S'='Function unknown')


################# DELTALLELE ################# 

#CalcDiffAFS <- function(data, samp1, samp2){
#  return(vapply(1:length(samp1),))}

#GetAFSDiff <- function(gene_data, pop1, pop2, comb_n){
#  gene_data = apply(gene_data, c(1,2), function(x) clr(x[[1]]))
#  combinations = combn(c(pop1,pop2),2, function(x) )
#  
#  return()}

deltAllele <- function(data, min_samp_per_group, samp_vec){
  print('Launching - MetaPoly DeltAllele: an allele frequency difference tool for metagenomic data')
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
                                                                              variable = as.numeric(names(samp_vec)),
                                                                              AFS_diff = as.vector(GetAFSDiff(data[[i]]$gene_data[,..cols])),
                                                                              depth = as.vector(GetDepth(data[[i]]$gene_data[,..cols])),
                                                                              gene_length = rep(data[[i]]$gene_length,length(samp_vec))))}}
  cat(' genes done\n')
  model_df = model_df[model_df$depth>0,]
  print(Sys.time() - t0)
  
  print(' - Fitting the poisson model on data...')
  model = glm(data = model_df, family = poisson(), formula = AFS_diff ~ depth + gene_length + sample, control = list(maxit = 100))
  print(summary(model))
  model_df$res_m = model$residuals
  print(Sys.time() - t0)
  
  print(' - Computing sample coefficients...')
  coefs = coefficients(model)
  coefs = coefs[startsWith(names(coefs),'sample')]
  coefs_df = as.data.frame(coefs)
  rownames(coefs_df) = vapply(rownames(coefs_df), function(x) strsplit(x, 'sample')[[1]][2], FUN.VALUE = character(1))
  coefs_df$type = vapply(rownames(coefs_df), function(x) as.numeric(names(samp_vec)[samp_vec == x]), numeric(1))
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
