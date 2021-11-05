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
    af2_filtered = af1[colnames(af1) %in% names(depth1)[depth1 > 9]]
    n_af1 = ncol(af1_filtered)
    n_af2 = ncol(af2_filtered)
    af1 = rowMeans(af1_filtered)
    af2 = rowMeans(af2_filtered)
    
    if ((n_af1 > 4) & (n_af2 > 4)){
         pi_between = mean(sample(1:length(af1), 1000, prob = af1, replace = T) != sample(1:length(af2), 1000, prob = af2, replace = T))
         pi1 = mean(sample(1:length(af1), 1000, prob = af1, replace = T) != sample(1:length(af1), 1000, prob = af1, replace = T))
         pi2 = mean(sample(1:length(af2), 1000, prob = af2, replace = T) != sample(1:length(af2), 1000, prob = af2, replace = T))
         pi_within = mean(c(pi1,pi2))
         fst = (pi_between - pi_within)/pi_between}
    else{pi_between = NA
         pi_within = NA
         fst=NA}
    fst_df = rbind(fst_df, data.frame(SNP=as.character(i),
                                      DEPTH1=mean(depth1, na.rm=T),
                                      DEPTH2=mean(depth2, na.rm=T),
                                      MEAN_DEPTH=mean(c(depth1,depth2)),
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
    if (nrow(data[[i]]$gene_data) > 0){gene_res = CalcFst(data[[i]]$gene_data, samp_vec)
                                       fst_df = rbind(fst_df , data.frame(gene_id = data[[i]]$gene_id,
                                                                          FST = gene_res$FST,
                                                                          DEPTH = gene_res$DEPTH,
                                                                          SNP_N = gene_res$SNP_N))}}
  cat(' genes done\n')
  print(Sys.time() - t0)
  print(' - Analysis done!')
  return(fst_df)}



