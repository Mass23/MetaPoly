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


################# PolyDiff ################# 

PairFst <- function(data, samp1, samp2, pos){
  ac1 = data[[samp1]][[pos]]
  ac2 = data[[samp2]][[pos]]
  depth1 = sum(ac1)
  depth2 = sum(ac2)
  fst = mean(sample(1:length(ac1), 100, prob = ac1, replace = T) != sample(1:length(ac2), 100, prob = ac2, replace = T))
  return(data.frame(fst = fst, depth_mean = exp(mean(log(c(depth1,depth2))))))}

CalcFst <- function(data, n_snp, samp_vec){
  df_within = rbind(expand.grid(samples_vec[names(samples_vec) == 0], samples_vec[names(samples_vec) == 0],1:n_snp),
                    expand.grid(samples_vec[names(samples_vec) == 1], samples_vec[names(samples_vec) == 1],1:n_snp))
  df_between = expand.grid(names(samp_vec)[samp_vec == 0], names(samp_vec)[samp_vec == 1], 1:n_snp)
  # remove comps of the same sample
  df_within = df_within[df_within$V1 != df_within$V2]
  df_between = df_between[df_between$V1 != df_between$V2]
  # if more than 100 comparisons, sub sample 
  if(nrow(df_within) > 100){
    df_within = df_within[sample(1:norw(df_within),100),]}
  if(nrow(df_between) > 100){
    df_between = df_between[sample(1:norw(df_between),100),]}
  # Calc Fst
  within_res = do.call(rbind, 1:nrow(df_within), function(i) PairFst(data, df_within$V1[i], df_within$V2[i], df_within$V3[i]))
  between_res = do.call(rbind, 1:nrow(df_between), function(i) PairFst(data, df_between$V1[i], df_between$V2[i], df_between$V3[i]))
  return(list(mean_fst = mean(between_res$fst / within_res$fst), mean_depth = mean(c(within_res$depth_mean, between_res$depth_mean))))
}

#' @export
PolyDiff <- function(data, min_samp_per_group, samp_vec){
  print('Launching - MetaPoly PolyDiff: an Fst calculatuon tool for metagenomic data')
  t0 = Sys.time()
  
  print(' - Computing Fst across genes...')
  fst_df = data.frame()
  cols = as.vector(samp_vec)
  count=0
  for (i in 1:length(data)){
    count = count +  1
    cat("\r",count)
    if (nrow(data[[i]]$gene_data) > 0){gene_res = CalcFst(data[[i]]$gene_data, samp_vec)
                                       fst_df = rbind(fst_df , data.frame(gene_id = rep(data[[i]]$gene_id,length(samp_vec)),
                                                                              sample = as.vector(samp_vec),
                                                                              variable = as.numeric(names(samp_vec)),
                                                                              fst = gene_res$mean_fst,
                                                                              depth = gene_res$mean_depth,
                                                                              gene_length = rep(data[[i]]$gene_length,length(samp_vec)))}}
  cat(' genes done\n')
  fst_df = fst_df[fst_df$depth>0,]
  print(Sys.time() - t0)
  
  print(' - Analysis done!')
  return(fst_df}
