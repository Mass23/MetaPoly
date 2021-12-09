################# PolyDiv ################# 

CalcFst <- function(data, samp_vec){
  fst_df = data.frame()
  g1 = as.vector(samp_vec[names(samp_vec) == 0])
  g2 = as.vector(samp_vec[names(samp_vec) == 1])
  for (i in 1:nrow(data)){
    snp_data = data[i,]
    depth1 = unlist(lapply(snp_data[,..g1], function(x) sum(x[[1]])))
    depth2 = unlist(lapply(snp_data[,..g2], function(x) sum(x[[1]])))
    af1 = as.data.frame(lapply(snp_data[,..g1], function(x) x[[1]]/sum(x[[1]])))
    af2 = as.data.frame(lapply(snp_data[,..g2], function(x) x[[1]]/sum(x[[1]])))
    af1_filtered = af1[colnames(af1) %in% names(depth1)[depth1 > 9]]
    af2_filtered = af2[colnames(af2) %in% names(depth2)[depth2 > 9]]
    n_af1 = ncol(af1_filtered)
    n_af2 = ncol(af2_filtered)
    af1 = rowMeans(af1_filtered)
    af2 = rowMeans(af2_filtered)
    
    pi_between = 1-sum(af1 * af2)
    pi1 = 1-sum(af1**2)
    pi2 = 1-sum(af2**2)
    pi_within = mean(c(pi1,pi2))
    fst = (pi_between - pi_within)/pi_between
    fst_df = rbind(fst_df, data.frame(SNP=as.character(i),
                                      MIN_DEPTH=min(c(depth1,depth2)),
                                      MEAN_DEPTH=exp(mean(log(c(depth1,depth2)))),
                                      PI_B=pi_between,
                                      PI_W=pi_within,
                                      FST=fst))}
  return(list(DEPTH=mean(fst_df$MEAN_DEPTH),FST=median(fst_df$FST),SNP_N=nrow(fst_df[is.finite(fst_df$FST),])))}

#' PolyDiv
#'
#' This functions computes Fst at the gene level, between two sets of metagenomic samples. It 
#' compares the average allele frequencies across samples for each polymorphic site, and then 
#' outputs the median Fst value. The inputs are:
#'   - data: the polymorphism data created by MetaPoly's "GetGenesData" function
#'   - samp_vec: a vector that assigns each sample (vector values) to each population (vector
#'     names). The method handles two populations encoded as 0 and 1's
#'
#' @export
PolyDiv <- function(data, samp_vec){
  print('Launching - MetaPoly PolyDiv: an Fst calculation tool for metagenomic data')
  t0 = Sys.time()
  
  print(' - Computing Fst across genes...')
  fst_df = data.frame()
  cols = as.vector(samp_vec)
  count=0
  for (i in 1:length(data)){
    count = count +  1
    cat("\r",count)
    if (nrow(data[[i]]$gene_data) > 0){gene_res = CalcFst(na.omit(data[[i]]$gene_data), samp_vec)
                                       fst_df = rbind(fst_df , data.frame(gene_id = data[[i]]$gene_id,
                                                                          FST = gene_res$FST,
                                                                          DEPTH = gene_res$DEPTH,
                                                                          SNP_N = gene_res$SNP_N))}}
  cat(' genes done\n')
  print(Sys.time() - t0)
  print(' - Analysis done!')
  return(fst_df)}



CalcFstNet <- function(data, samp_vec){
  fst_df = data.frame()
  g1 = as.vector(samp_vec[names(samp_vec) == 0])
  g2 = as.vector(samp_vec[names(samp_vec) == 1])
  rows_sel = 1:nrow(data)
  if (nrow(data) > 10){
    rows_sel = sample(1:nrow(data), 10)}
  for (i in rows_sel){
    snp_data = data[i,..samp_vec]
    depth1 = unlist(lapply(snp_data[,..g1], function(x) sum(x[[1]])))
    depth2 = unlist(lapply(snp_data[,..g2], function(x) sum(x[[1]])))
    af1 = as.data.frame(lapply(snp_data[,..g1], function(x) x[[1]]/sum(x[[1]])))
    af2 = as.data.frame(lapply(snp_data[,..g2], function(x) x[[1]]/sum(x[[1]])))
    n_af1 = ncol(af1)
    n_af2 = ncol(af2)
    af1 = rowMeans(af1)
    af2 = rowMeans(af2)
    
    pi_between = 1-sum(af1 * af2)
    pi1 = 1-sum(af1**2)
    pi2 = 1-sum(af2**2)
    pi_within = mean(c(pi1,pi2))
    fst = (pi_between - pi_within)/pi_between
    fst_df = rbind(fst_df, data.frame(SNP=as.character(i),
                                      MIN_DEPTH=min(c(depth1,depth2)),
                                      MEAN_DEPTH=exp(mean(log(c(depth1,depth2)))),
                                      PI_B=pi_between,
                                      PI_W=pi_within,
                                      FST=fst))}
  return(list(DEPTH=mean(fst_df$MIN_DEPTH),FST=median(fst_df$FST),SNP_N=nrow(fst_df[is.finite(fst_df$FST),])))}


#' PolyNet
#'
#' to write
#'
#' @export
PolyNet <- function(data, samp_vec, poly_summary, depth_thr, snp_n_thr){
  print('Launching - MetaPoly PolyNet: an Fst network creation tool for metagenomic data')
  t0 = Sys.time()
  
  print(' - Listing genes with higher min. depth than threshold...')
  genes_to_keep = c()
  for (gene in unique(poly_summary$gene_id)){
    if (sum(is.na(poly_summary$DEPTH[poly_summary$gene_id == gene])) == 0){
    if (min(poly_summary$DEPTH[poly_summary$gene_id == gene]) >= depth_thr){
    if (min(poly_summary$SNP_N[poly_summary$gene_id == gene]) >= snp_n_thr){
      genes_to_keep = c(genes_to_keep, gene)}}}}
  print(paste0('   - Genes kept N:', length(genes_to_keep)))
  
  print(' - Computing Fst across genes...')
  fst_df = data.frame()
  cols = as.vector(samp_vec)
  count=0
  for (i in 1:length(data)){
    if (nrow(data[[i]]$gene_data[,..cols]) > 0){
      if (data[[i]]$gene_id %in% genes_to_keep){
        count = count +  1
        cat("\r",count)
        pair_comb = as.data.frame(combinations(n = length(samp_vec), r = 2, v = samp_vec, repeats.allowed = FALSE))
        for (j in 1:nrow(pair_comb)){
          comb_vec = c(as.character(pair_comb$V1[j]),as.character(pair_comb$V2[j]))
          names(comb_vec) = c(0,1)
          comb_res = CalcFstNet(data[[i]]$gene_data[,..comb_vec], comb_vec)
          fst_df = rbind(fst_df , data.frame(gene_id = data[[i]]$gene_id,
                                             s1 = pair_comb$V1[j],
                                             s2 = pair_comb$V2[j],
                                             FST = comb_res$FST,
                                             SNP_N = comb_res$SNP_N))}}}}
  cat(' genes done\n')
  print(Sys.time() - t0)
  print(' - Analysis done!')
  return(fst_df)}

