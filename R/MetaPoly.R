library(vcfR)
library(ape)
library(stringr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(compositions)
library(seqinr)

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
  for (j in cols) set(gene_data, j = j, value = lapply(strsplit(gene_data[[j]], ','), as.integer))
  
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

################# POLYCORR ################# 
GetDepth <- function(gene_data){
  return(colMeans(apply(gene_data, c(1,2), function(x) sum(unlist(x)))))}

GetSnpN <- function(gene_data){
  return(colSums(apply(gene_data, c(1,2), function(x) ifelse(length(x[[1]][x[[1]] > 0]), 1, 0))))}

fit_cor_gene <- function(gene_data, gene, min_samp, samp_vec){
  gene_data = gene_data[gene_data$depth > 9,]
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
      if (nrow(data[[i]]$gene_data) > 0){model_df = rbind(model_df , data.frame(gene_id = rep(data[[i]]$gene_id,length(samp_vec)),
                                                                                sample = as.vector(samp_vec),
                                                                                variable = as.numeric(names(samp_vec)),
                                                                                snp_den = as.vector(GetSnpN(data[[i]]$gene_data[,..cols])),
                                                                                depth = as.vector(GetDepth(data[[i]]$gene_data[,..cols])),
                                                                                gene_length = rep(data[[i]]$gene_length,length(samp_vec))))}}
  cat(' genes done\n')
  model_df = model_df[model_df$depth>0,]
  print(Sys.time() - t0)
  
  print(' - Fitting the poisson model on data...')
  model = glm(data = model_df, family = poisson(), formula = snp_den ~ log(depth) + gene_length + sample, control = list(maxit = 100))
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

################# MKTEST ################# 
calc_dn_ds_pn_ps <- function(gene_data, gene_seq, gene_strand, gene_start, gene_end, min_depth, samp_vec, alleles_data){
  for (i in 1:nrow(gene_data)){
    snp_data = gene_data[i,]
    
    # get the codon position and snp pos in codon
    codon_n = floor((snp_data[['POS']] - gene_start) / 3)
    codon_start = codon_n*3
    codon_end = (codon_n*3)+2
    pos_snp = snp_data[['POS']] - gene_start
    pos_in_codon = (3*((pos_snp/3) - floor(pos_snp/3)) ) + 1
    ref_codon = gene_seq[codon_start:codon_end]
    if (gene_strand == '-'){ref_codon = rev(ref_codon)}
    
    # get the amino acids for each allele
    alleles = lapply(alleles_data$alleles[(alleles_data$CHROM == snp_data$CONTIG) & (alleles_data$POS == snp_data$POS)], tolower)
    amino_acids = unlist(lapply(alleles[[1]], function(x)translate(replace(ref_codon,pos_in_codon,x))))
    print(amino_acids)
    
  }
}

MKTest <- function(data, vcf, min_depth = 10, samp_vec, fasta){
  print('Launching - MetaPoly MKTest: a McDonald-Kreitman test tool for metagenomic data')
  t0 = Sys.time()
  
  print(' - Loading fasta and GFF data')
  fasta_loaded = read.fasta(fasta)
  alleles_data = as.data.frame(vcf@fix)[,c('CHROM','POS','REF','ALT')]
  alleles_data$POS = as.integer(alleles_data$POS)
  
  print(' - Computing dn, ds, pn ,ps...')
  mk_df = data.frame()
  cols = as.vector(samp_vec)
  count=0
  for (i in 1:10){
    count = count +  1
    cat("\r",count)
    if (nrow(data[[i]]$gene_data) > 0){
      # Subset snps in the gene
      alleles_data_sub = alleles_data[(alleles_data$CHROM == data[[i]]$gene_contig) & (alleles_data$POS > data[[i]]$gene_start) & (alleles_data$POS < data[[i]]$gene_end),]
      alleles_data_sub$snp = vapply(1:nrow(alleles_data_sub), function(i) paste(alleles_data_sub$CHROM[i], alleles_data_sub$POS[i], sep = '_'), FUN.VALUE = character(1))
      alleles_data_sub$alleles = strsplit(paste0(alleles_data_sub$REF, alleles_data_sub$ALT),'')
      alleles_data_sub = alleles_data_sub[alleles_data_sub$snp %in% data[[i]]$gene_data[['rn']],]

      gene_seq = getSequence(fasta_loaded[data[[i]]$gene_contig])[[1]][data[[i]]$gene_start:data[[i]]$gene_end]
      gene_res = calc_dn_ds_pn_ps(data[[i]]$gene_data, gene_seq, data[[i]]$gene_strand, data[[i]]$gene_start, data[[i]]$gene_end, min_depth, samp_vec, alleles_data_sub)
      mk_df = rbind(mk_df , data.frame(gene_id = data[[i]]$gene_id,
                                       dn = gene_res$dn,
                                       ds = gene_res$ds,
                                       pn = gene_res$pn,
                                       ps = gene_res$ps,
                                       mk = (gene_res$dn / gene_res$ds) / (gene_res$pn / gene_res$ps),
                                       gene_length = data[[i]]$gene_length))}}
  cat(' genes done\n')
  print(' - Analysis done!')
  return(mk_df)}


#data_mt = GetGenesData(gff, vcf)
#MKTest(data_mt, vcf, min_depth = 10, samples_vec, 'Bio17-1_NCBI.fa')






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

################# EVENNESS ################# 

CalcEvenness <- function(AC){
  AF = AC / sum(AC)
  return(sum(-vapply(AF, function(x) x * log(x), FUN.VALUE = numeric(1)))/log(length(AF)))}

GetEvenness <- function(gene_data){
  even_data = apply(gene_data, c(1,2), function(x) CalcEvenness(as.matrix(x)[1][[1]]))
  even_data[even_data == 0] = NA
  return(colMeans(even_data, na.rm = T))}

EvenCalc <- function(data, samp_vec){
  print('Launching - MetaPoly EVENNESS: a polymorphism quantification tool for metagenomic data')
  t0 = Sys.time()
  
  print(' - Computing Evenness and SNP density...')
  evenness_df = data.frame()
  cols = as.vector(samp_vec)
  count=0
  for (i in 1:length(data)){
    count = count +  1
    cat("\r",count)
    if (nrow(data[[i]]$gene_data) > 0){
      evenness_df = rbind(evenness_df, data.frame(gene_id = rep(data[[i]]$gene_id,length(samp_vec)),
                                                  sample = as.vector(samp_vec),
                                                  variable = as.numeric(names(samp_vec)),
                                                  evenness_poly = as.vector(GetEvenness(data[[i]]$gene_data[,..cols])),
                                                  snp_n = as.vector(GetSnpN(data[[i]]$gene_data[,..cols])),
                                                  depth = as.vector(GetDepth(data[[i]]$gene_data[,..cols])),
                                                  gene_length = rep(data[[i]]$gene_length,length(samp_vec))))}}
  cat(' genes done\n')
  evenness_df$prop_poly = evenness_df$snp_n / evenness_df$gene_length
  evenness_df$evenness_all = vapply(1:nrow(evenness_df), function(i) weighted.mean(c(evenness_df$evenness_poly[i], 0), c(evenness_df$prop_poly[i], 1-evenness_df$prop_poly[i])), FUN.VALUE=numeric(1))
  print(' - Analysis done!')
  return(evenness_df)}

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

