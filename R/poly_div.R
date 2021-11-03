################# PolyDiff ################# 
geom_mean <- function(x,na.rm=TRUE){return(exp(mean(log(x),na.rm=na.rm)))}


PairPI <- function(data, s1, s2){
  ac1 = as.integer(data[[as.character(s1)]][[1]])
  depth1 = sum(ac1)
  af1 = ac1/depth1
  ac2 = as.integer(data[[as.character(s2)]][[1]])
  depth2 = sum(ac2)
  af2 = ac2/depth2
  pi_between = mean(sample(1:length(af1), 1000, prob = af1, replace = T) != sample(1:length(af2), 1000, prob = af2, replace = T))
  pi_within = mean(sample(1:length(af1), 1000, prob = colMeans(rbind(af1,af2)), replace = T) != sample(1:length(af1), 1000, prob = colMeans(rbind(af1,af2)), replace = T))
  return(data.frame(DEPTH1=depth1, DEPTH2=depth2, DEPTH_GM=geom_mean(depth1, depth2), PI_b=pi_between, PI_w=pi_within))}

CombType <- function(samp_vec, samp1, samp2){
  if (names(samp_vec)[samp_vec == samp1] == names(samp_vec)[samp_vec == samp2]){return('within')}
  else{return('between')}}

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
         pi_within = mean(sample(1:length(af1), 1000, prob = colMeans(rbind(af1,af2)), replace = T) != sample(1:length(af1), 1000, prob = colMeans(rbind(af1,af2)), replace = T))
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
  return(list(DEPTH=mean(fst_df$MEAN_DEPTH),FST=median(fst_df$FST)))}

#' @export
PolyDiv <- function(data, samp_vec){
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
                                       fst_df = rbind(fst_df , data.frame(gene_id = data[[i]]$gene_id,
                                                                          FST = gene_res$FST,
                                                                          DEPTH = gene_res$DEPTH))}}
  cat(' genes done\n')
  print(Sys.time() - t0)
  print(' - Analysis done!')
  return(fst_df)}

##### Other funcs
# computing fst for every pair...
CalcFstPerSample <- function(data, samp_vec){
  combn_df = as.data.frame(expand.grid(samp_vec, samp_vec))
  # remove comps of the same sample
  combn_df = combn_df[as.character(combn_df$Var1) != as.character(combn_df$Var2),]
  
  full_df = data.frame()
  for (i in 1:nrow(data)){
    snp_data = data[i,]
    fst_combn = expand.grid(samples_vec, samples_vec)
    fst_combn$SNP = as.character(i)
    fst_combn = fst_combn[fst_combn$Var1 != fst_combn$Var2,]
    res_df = lapply(as.list(1:nrow(fst_combn)), function(i) PairPI(snp_data, fst_combn$Var1[[i]][1], fst_combn$Var2[[i]][1]))
    res_df = do.call(rbind, res_df)
    fst_combn = cbind(fst_combn, res_df)
    full_df = rbind(full_df, fst_combn)}
  return(full_df)}
