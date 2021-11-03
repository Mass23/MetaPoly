################# POLYCORR ################# 
fit_cor_gene <- function(gene_data, gene, min_samp, samp_vec){
  gene_data = na.omit(gene_data[gene_data$DEPTH > 9,])
  if (length(unique(gene_data$sample)) > min_samp){
  cor_test = cor.test(gene_data$res_m, gene_data$variable)
  p = as.numeric(cor_test$p.value)
  cor = as.numeric(cor_test$estimate)
  mean_coef = mean(gene_data$res_m, na.rm=T)
  low_ci = cor_test$conf.int[1]
  high_ci = cor_test$conf.int[2]
  return(data.frame(gene_id=gene, cor=cor, p=p, mean_coef=mean_coef, low_ci=low_ci, high_ci=high_ci))}
  else{return(data.frame(gene_id=gene, cor=NA, p=NA, mean_coef=NA, low_ci=NA, high_ci=NA))}}

#' PolyCorr
#'
#' This functions computes the gene-wise correlation between the number of SNPs, and a target 
#' variable than can be either continuous, either discrete (handles two treatments). It creates
#' a zero-inflated poisson model of the SNP count, taking the sample, log(sequencing depth), 
#' and the gene length as fixed effects. The idea is to identify genes showing signs of selective
#' sweeps associated with the target variable. There are three inputs:
#'   - data: the polymorphism data created by MetaPoly's "GetGenesData" function
#'   - min_samp: the minimal number of sample having the required depth to run the correlation test.
#'   - samp_vec: a vector that assigns each sample (vector values) to each population (vector
#'     names). The method handles two populations encoded as 0 and 1's.
#'
#' @export
PolyCorr <- function(data, min_samp, samp_vec){
  print('Launching - MetaPoly PolyCorr: a polymorphism-variable correlation tool for metagenomic data')
  t0 = Sys.time()
  
  print(' - Fitting a zero-inflated poisson model on the data...')
  data = as.data.table(data[,colnames(data) %in% c('SNP_N', 'DEPTH', 'gene_length', 'sample', 'gene_id')])
  model_df = data[data$DEPTH>0,]
  model = pscl::zeroinfl(SNP_N ~ log(DEPTH) + gene_length + as.factor(sample) | log(DEPTH), data = model_df)
  print(summary(model))
  model_df$res_m = model$residuals
  model_df$variable = vapply(model_df$sample, function(x) as.numeric(names(samp_vec)[samp_vec == x]), numeric(1))
  print(Sys.time() - t0)
  
  print(' - Computing sample coefficients...')
  coefs = coefficients(model)
  coefs = coefs[startsWith(names(coefs),'sample')]
  coefs_df = as.data.frame(coefs)
  rownames(coefs_df) = vapply(rownames(coefs_df), function(x) strsplit(x, '(sample)')[[1]][2], FUN.VALUE = character(1))
  coefs_df$type = vapply(rownames(coefs_df), function(x) as.numeric(names(samp_vec)[samp_vec == x]), numeric(1))
  print(Sys.time() - t0)
  
  print(' - Computing correlations of polymorphism with the variable of interest per gene...')
  corr_df = do.call(rbind, lapply(unique(model_df$gene_id), function(gene) fit_cor_gene(model_df[model_df$gene_id == gene,], gene, min_samp, samp_vec)))
  corr_df$padj = p.adjust(corr_df$p, method = 'holm')
  sign_genes = corr_df[corr_df$padj < 0.05,]
  pos_genes = as.vector(na.omit(sign_genes$gene_id[sign_genes$cor > 0]))
  neg_genes = as.vector(na.omit(sign_genes$gene_id[sign_genes$cor < 0]))
  print(Sys.time() - t0)
  
  print(' - Analysis done!')
  return(list(pi_corr_res = corr_df, pos_genes = pos_genes, neg_genes = neg_genes, coefs = coefs_df))}

#' @export
PlotPolyCorr <- function(res_df, coefs_df, plots_name, boolean_var = TRUE){
  dir.create(paste0(plots_name,'_res'))
  
  res_df$Significant = 'No'
  res_df$Significant[res_df$padj < 0.05] = 'Yes'
  
  # Plot the distribution of correlations
  ggplot2::ggplot(res_df, aggplot2::es(x=cor,fill=Significant)) + ggplot2::geom_histogram() + ggplot2::theme_minimal() +
    ggplot2::geom_vline(ggplot2::aes(xintercept=mean(res_df$cor, na.rm = T)), linetype="dashed") + 
    ggplot2::geom_vline(ggplot2::aes(xintercept=0), color='darkgrey', linetype="dashed") + 
    ggplot2::xlab('Correlation coef.') + ggsci::scale_fill_jco() 
  ggplot2::ggsave(paste0(plots_name,'_res','/',plots_name,'_dist.pdf'), width = 4, height = 4)
  
  # Plot the distribution of p-values
  ggplot2::ggplot(res_df, ggplot2::aes(x=p,fill=Significant)) + ggplot2::geom_histogram() + ggplot2::theme_minimal() + 
    xlab('p') + scale_fill_jco() 
  ggplot2::ggsave(paste0(plots_name,'_res','/',plots_name,'_p.pdf'), width = 4, height = 4)
  
  # Plot the coefficients
  if (boolean_var == TRUE){
    ggplot2::ggplot(coefs_df, ggplot2::aes(x=as.factor(type),y=coefs,color=as.factor(type))) + ggplot2::geom_boxplot() + ggplot2::theme_minimal() + 
      ggpubr::stat_compare_means() + ggsci::scale_color_jco() + ggplot2::xlab('Variable') + ggplot2::ylab('Sample coefs.') + ggplot2::theme(legend.position = 'none')
    ggplot2::ggsave(paste0(plots_name,'_res','/',plots_name,'_coefficients.pdf'), width = 3, height = 4)}
  
  else {ggplot2::ggplot(coefs_df, ggplot2::aes(x=type,y=coefs)) + ggplot2::geom_point() + ggplot2::theme_minimal() + ggplot2::geom_smooth() +
        ggsci::scale_color_jco() + ggplot2::xlab('Variable') + ggplot2::ylab('Sample coefs.') + ggplot2::theme(legend.position = 'none')
    ggplot2::ggsave(paste0(plots_name,'_res','/',plots_name,'_coefficients.pdf'), width = 3, height = 4)}}

################ ENRICH ########
CalcEnrichment <- function(gff, gene_list, cog_table){
  gff$gene = vapply(gff$V9, function(x) strsplit(strsplit(x,';')[[1]][1],'ID=')[[1]][2], FUN.VALUE = character(1))
  gff$cog = vapply(gff$V9, function(x) strsplit(strsplit(x,'COG:')[[1]][2],';')[[1]][1], FUN.VALUE = character(1))
  gff$cog_f = vapply(gff$cog, function(x) ifelse(is.na(x), 'NoCOG', substr(cog_func$V2[cog_func$V1 == x],1,1)), FUN.VALUE = character(1))
  
  gff$sign = 'b_No'
  gff$sign[gff$gene %in% gene_list] = 'a_Yes'

  full_con_tab = table(gff$cog_f, gff$sign)
  enrich_df = data.frame()
  for (cog_func in rownames(full_con_tab)){
    if(cog_func != 'NoCOG'){
      func_vals = full_con_tab[rownames(full_con_tab) == cog_func,]
      other_vals =  colSums(full_con_tab[rownames(full_con_tab) != cog_func,])
      cont_table = as.matrix(rbind(func_vals, other_vals))
      fish_test = fisher.test(cont_table, 'greater')
      p = fish_test$p.value
      or = fish_test$estimate
      l_or = fish_test$conf.int[1]
      h_or = fish_test$conf.int[2]
      enrich_df = rbind(enrich_df, data.frame(Function=cog_func, p=p, OR=or, low_CI=l_or, high_CI=h_or))}}
  enrich_df$padj = p.adjust(enrich_df$p, method = 'holm')
  enrich_df$Function_long = vapply(enrich_df$Function, function(x) as.character(cog_functions[x]), character(1))
  return(list(enrich=enrich_df,cogs=as.vector(na.omit(unique(gff$cog[gff$gene %in% gene_list])))))}

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
