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

################# GET DATA PER GENE ################# 
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

GetGenesData <- function(gff, vcf){
  t0 = Sys.time()
  # Load gff and vcf, convert to data.table, get pos of the VCF variants
  print('Launching - MetaPoly GetGenesData: loading VCF and GFF data...')

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

############### POLYMORPHISM SUMMARY ################ 
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

GetSnpData <- function(gene_data){
  depth = colMeans(apply(gene_data, c(1,2), function(ac) sum(unlist(ac))))
  snp_n = colSums(apply(gene_data, c(1,2), function(ac) ifelse(length(ac[[1]][ac[[1]] > 0]) > 1, 1, 0)))
  evenness = colMeans(apply(gene_data, c(1,2), function(ac) CalcEven(as.matrix(ac)[1][[1]])), na.rm = T)
  majf = colMeans(apply(gene_data, c(1,2), function(ac) CalcMAJF(as.matrix(ac)[1][[1]])), na.rm = T)
  return(list(depth=depth,snp_n=snp_n,evenness=evenness,majf=majf))}

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
                                                                          gene_length = rep(data[[i]]$gene_length,length(samp_vec))))}}
  cat(' genes done\n')
  return(PolyDf)}







################# POLYCORR ################# 


fit_cor_gene <- function(gene_data, gene, min_samp, samp_vec){
  gene_data = na.omit(gene_data[gene_data$depth > 9,])
  if (length(unique(gene_data$sample)) > min_samp){
  cor_test = cor.test(gene_data$res_m, gene_data$variable)
  p = as.numeric(cor_test$p.value)
  cor = as.numeric(cor_test$estimate)
  mean_coef = mean(gene_data$res_m, na.rm=T)
  low_ci = cor_test$conf.int[1]
  high_ci = cor_test$conf.int[2]
  return(data.frame(gene_id=gene, cor=cor, p=p, mean_coef=mean_coef, low_ci=low_ci, high_ci=high_ci))}
  else{return(data.frame(gene_id=gene, cor=NA, p=NA, mean_coef=NA, low_ci=NA, high_ci=NA))}}

PolyCorr <- function(data, min_samp, samp_vec){
  print('Launching - MetaPoly PolyCorr: a polymorphism-variable correlation tool for metagenomic data')
  t0 = Sys.time()
  
  print(' - Computing SNP density and depth...')
  model_df = data.frame()
  cols = as.vector(samp_vec)
  count=0
  for (i in 1:length(data)){
      count = count +  1
      cat("\r",count)
      if (nrow(data[[i]]$gene_data) > 0){snp_res = GetSnpN(data[[i]]$gene_data[,..cols])
                                         model_df = rbind(model_df , data.frame(gene_id = rep(data[[i]]$gene_id,length(samp_vec)),
                                                                                sample = as.vector(samp_vec),
                                                                                variable = as.numeric(names(samp_vec)),
                                                                                snp_n = snp_res$snp_n,
                                                                                depth = snp_res$depth,
                                                                                gene_length = rep(data[[i]]$gene_length,length(samp_vec))))}}
  cat(' genes done\n')
  model_df = model_df[model_df$depth>0,]
  print(Sys.time() - t0)
  
  print(' - Fitting the poisson model on data...')
  model = glm(data = model_df, family = poisson(), formula = snp_n ~ log(depth) + gene_length + sample, control = list(maxit = 100))
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
  
  print(' - Computing correlations of polymorphism with the variable of interest per gene...')
  corr_df = do.call(rbind, lapply(unique(model_df$gene_id), function(gene) fit_cor_gene(model_df[model_df$gene_id == gene,], gene, min_samp, samp_vec)))
  corr_df$padj = p.adjust(corr_df$p, method = 'holm')
  sign_genes = corr_df[corr_df$padj < 0.05,]
  pos_genes = as.vector(na.omit(sign_genes$gene_id[sign_genes$cor > 0]))
  neg_genes = as.vector(na.omit(sign_genes$gene_id[sign_genes$cor < 0]))
  print(Sys.time() - t0)
  
  print(' - Analysis done!')
  return(list(pi_corr_res = corr_df, pos_genes = pos_genes, neg_genes = neg_genes, coefs = coefs_df))}

PlotPolyCorr <- function(res_df, coefs_df, plots_name, boolean_var = TRUE){
  dir.create(paste0(plots_name,'_res'))
  
  res_df$Significant = 'No'
  res_df$Significant[res_df$padj < 0.05] = 'Yes'
  
  # Plot the distribution of correlations
  ggplot(res_df, aes(x=cor,fill=Significant)) + geom_histogram() + theme_minimal() +
    geom_vline(aes(xintercept=mean(res_df$cor, na.rm = T)), linetype="dashed") + 
    geom_vline(aes(xintercept=0), color='darkgrey', linetype="dashed") + 
    xlab('Correlation coef.') + scale_fill_jco() 
  ggsave(paste0(plots_name,'_res','/',plots_name,'_dist.pdf'), width = 4, height = 4)
  
  # Plot the distribution of p-values
  ggplot(res_df, aes(x=p,fill=Significant)) + geom_histogram() + theme_minimal() + 
    xlab('p') + scale_fill_jco() 
  ggsave(paste0(plots_name,'_res','/',plots_name,'_p.pdf'), width = 4, height = 4)
  
  # Plot the coefficients
  if (boolean_var == TRUE){
    ggplot(coefs_df, aes(x=as.factor(type),y=coefs,color=as.factor(type))) + geom_boxplot() + theme_minimal() + 
      stat_compare_means() + scale_color_jco() + xlab('Variable') + ylab('Sample coefs.') + theme(legend.position = 'none')
    ggsave(paste0(plots_name,'_res','/',plots_name,'_coefficients.pdf'), width = 3, height = 4)}
  
  else {ggplot(coefs_df, aes(x=type,y=coefs)) + geom_point() + theme_minimal() + geom_smooth() +
        scale_color_jco() + xlab('Variable') + ylab('Sample coefs.') + theme(legend.position = 'none')
    ggsave(paste0(plots_name,'_res','/',plots_name,'_coefficients.pdf'), width = 3, height = 4)}}

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

