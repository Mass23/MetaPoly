setwd('~/Documents/PhD/Others/MetaPoly/MetaPoly')
library(devtools)
install_github('https://github.com/Mass23/MetaPoly', force = T, upgrade = 'never')
library(MetaPoly)

library(ggplot2)
library(ggpubr)
library(ggsci)

vcf = vcfR::read.vcfR('data/WWTP/concatenated_nonMP_Bio17-1_filtered.bcf.gz')
genome = ape::read.dna('data/WWTP/Bio17-1_NCBI.fa', format = "fasta")
gff <- read.delim("data/WWTP/Bio17-1.gff", header=F, comment.char="#", sep='\t', quote = '')
gff = gff[gff$V3 == 'CDS',]

metadata = read.csv('data/WWTP/WWTP_samples.csv')
metadata = metadata[grepl('D',metadata$Sample),]
samples_vec = metadata$Sample
metadata$Date = as.Date(metadata$Date, format = "%Y-%m-%d")
metadata$time_diff = abs(metadata$Date - metadata$Date[metadata$Sample=='D15'])
names(samples_vec) = as.numeric(vapply(samples_vec, function(x) metadata$Date[metadata$Sample == x], FUN.VALUE = numeric(1)))

# Load data, get polymorphism summary
data_mt = GetGenesData(gff, vcf)

######################################################################################################
# POLYSUMMARY
# 1. data generation
mt_poly = PolySummary(data_mt, samples_vec[grepl('D',samples_vec)], 6)
write.csv(mt_poly, file = 'data/WWTP/WWTP_PolySummary.csv', row.names = F, quote = F)
mt_poly = read.csv('data/WWTP/WWTP_PolySummary.csv')
# generate sample summaries
mt_samples = SummariseSamples(mt_poly, 5)

# 2. data visualisation
mt_samples$table$season = vapply(mt_samples$table$sample, function(x) metadata$Season[metadata$Sample == x], FUN.VALUE = character(1))
mt_samples$table$date = as.Date(vapply(mt_samples$table$sample, function(x) as.character(metadata$Date[metadata$Sample == x]), FUN.VALUE = character(1)))
mt_samples$table$time_diff = vapply(mt_samples$table$sample, function(x) metadata$time_diff[metadata$Sample == x], FUN.VALUE = numeric(1))
mt_samples$table$group = 'Baseline'
mt_samples$table$group[mt_samples$table$sample %in% c('D05','D15')] = 'Shift'

p1 = ggplot() + geom_point(mt_samples$table, mapping = aes(x=date,y=MODEL_INT,color=season,shape=group), size=5, alpha=0.7) + 
  xlab('') + ylab('Model Intercept') + scale_color_jco() + theme_minimal() + theme(axis.text.x = element_blank()) 
p2 = ggplot() + geom_point(mt_samples$table, mapping = aes(x=date,y=MODEL_SLOPE,color=season,shape=group), size=5, alpha=0.7) + 
  xlab('Date') + ylab('Model Slope') + scale_color_jco() + theme_minimal() 
p3 = ggplot(mt_samples$table, mapping = aes(x=log(MEAN_DEPTH),y=MEAN_SNP_DEN,color=season,shape=group)) + geom_point(size=5, alpha=0.7) + 
  xlab('log(Seq. depth)') + ylab('SNP density') + geom_smooth(method='lm',se=F) + scale_color_jco() + theme_minimal() 

ggarrange(p1,p2,p3, nrow = 3,ncol=1,align='h', common.legend = TRUE, legend = 'right', heights = c(2, 2, 4))
ggsave('figures/fig1_WWTP_poly_summary.jpg', width=7,height = 8)

# 3. data analysis
mod_snp_n_all = lm(data=mt_samples$table, MEAN_SNP_DEN ~ log(MEAN_DEPTH):season -1)
summary(mod_snp_n_all) # all seasons:log(DEPTH) interactions p < 0.0001, r2 = 0.7588
mod_snp_n_shift_autumn = lm(data=mt_samples$table[mt_samples$table$season == 'Autumn',], MEAN_SNP_DEN ~ log(MEAN_DEPTH) + group)
summary(mod_snp_n_shift_autumn) # shift p = 4.06e-08   , adj. r2 = 0.9995

mod_ndiv_all = lm(data=mt_samples$table, MEAN_PIP ~ log(MEAN_DEPTH):season)
summary(mod_ndiv_all) # all seasons:log(DEPTH) interactions p < 0.0001, r2 = 0.9261
mod_ndiv_shift = lm(data=mt_samples$table[mt_samples$table$season %in% c('Autumn'),], MEAN_PIP ~ log(MEAN_DEPTH) + group)
summary(mod_ndiv_shift) # shift p = 0.09277 , adj. r2 = 0.09277 

######################################################################################################
# POLYCORR
# 1. data generation
metadata$time_diff = log(abs(as.integer(metadata$Date - as.Date(metadata$Date[metadata$Sample == 'D15']))+1))
samples_vec = metadata$Sample
names(samples_vec) = metadata$time_diff

mt_polycorr = PolyCorr(mt_poly, 0, samples_vec, 'holm')
PlotPolyCorr(mt_polycorr$pi_corr_res, mt_polycorr$coefs, 'fig_polycorr', boolean_var = F)

cog_func = read.csv('data/cog-20.def.tab',sep = '\t',header = F)
gff$gene = vapply(gff$V9, function(x) strsplit(strsplit(x,';')[[1]][1],'ID=')[[1]][2], FUN.VALUE = character(1))
gff_to_test = gff[gff$gene %in% mt_polycorr$pi_corr_res$gene_id,]
pos_enrich = CalcEnrichment(gff_to_test, mt_polycorr$pos_genes)
neg_enrich = CalcEnrichment(gff_to_test, mt_polycorr$neg_genes)

# 2. data analysis
pos_enrich$enrich[(pos_enrich$enrich$padj < 0.05) & (pos_enrich$enrich$OR > 1),]
neg_enrich$enrich[(neg_enrich$enrich$padj < 0.05) & (neg_enrich$enrich$OR > 1),]
#             Function           p       OR   low_CI  high_CI                                                            cogs        padj                  Function_long
# odds ratio7        I 0.000462648 3.780909 1.771067 7.357895 COG2230,COG1960,COG1012,COG1804,COG1022,COG0743,COG1028,COG0761  0.01017826 Lipid transport and metabolism
# COG2230 - Cyclopropane fatty-acyl-phospholipid synthase and related methyltransferases
# COG1960 - Acyl-CoA dehydrogenase related to the alkylation response protein (Fatty acid biosynthesis)
# COG1012 - Acyl-CoA reductase or other NAD-dependent aldehyde dehydrogenase (Proline degradation)
# COG1804 - Crotonobetainyl-CoA:carnitine CoA-transferase CaiB and related acyl-CoA transferases
# COG1022 - Long-chain acyl-CoA synthetase (AMP-forming)
# COG0743 - 1-deoxy-D-xylulose 5-phosphate reductoisomerase (Isoprenoid biosynthesis)
# COG1028 - NAD(P)-dependent dehydrogenase, short-chain alcohol dehydrogenase family (Fatty acid biosynthesis)
# COG0761 - 4-Hydroxy-3-methylbut-2-enyl diphosphate reductase (Isoprenoid biosynthesis)

# 3. data visualisation
ggplot(neg_enrich$enrich[(neg_enrich$enrich$OR > 1),]) + geom_point(aes(x=OR,y=Function_long,color=Function_long),size=10) + geom_errorbarh(aes(xmin=low_CI,xmax=high_CI,y=Function_long,color=Function_long),size=2) + 
  xlab('Odds ratio') + ylab('') + theme_minimal() + scale_color_jco() + theme(legend.position = 'none') + geom_vline(xintercept = 1, linetype='dashed')
ggsave('figures/Neg_enrich.pdf', width=8, height = 5)


######################################################################################################
# POLYCORR
# 1. data generation
samples_vec = metadata$Sample
names(samples_vec) = as.integer(metadata$Date > as.Date(metadata$Date[metadata$Sample == 'D15']))

mt_div = PolyDiv(data_mt, samples_vec)
mt_div_norm = na.omit(mt_div[mt_div$SNP_N >= 5,])
fst_outliers = mt_div_norm$gene_id[mt_div_norm$FST > quantile(mt_div_norm$FST, probs = 0.95)]
mt_div_norm$start = vapply(mt_div_norm$gene_id, function(x) gff$V4[gff$gene == x], FUN.VALUE = numeric(1))
mt_div_norm$contig = vapply(mt_div_norm$gene_id, function(x) gff$V1[gff$gene == x], FUN.VALUE = character(1))

ggplot(mt_div_norm, aes(x=start,y=FST)) + facet_grid(~contig, space = 'free_x', scales = 'free_x') + theme_minimal() + geom_point()

mt_div_norm[mt_div_norm$FST > 0.05,]

gff[gff$gene == 'HAEMFFGE_00989',]
gff[gff$gene == 'HAEMFFGE_00990',]
gff[gff$gene == 'HAEMFFGE_00991',]
gff[gff$gene == 'HAEMFFGE_00993',]
gff[gff$gene == 'HAEMFFGE_00996',]

gff[gff$gene == 'HAEMFFGE_01277',]

gff[gff$gene == 'HAEMFFGE_02900',]

gff[gff$gene == 'HAEMFFGE_03314',]

# 2. data analysis
gff$gene = vapply(gff$V9, function(x) strsplit(strsplit(x,';')[[1]][1],'ID=')[[1]][2], FUN.VALUE = character(1))
gff_to_test = gff[gff$gene %in% mt_div_norm$gene_id,]
fst_enrichment = CalcEnrichment(gff_to_test, fst_outliers)
fst_enrichment$enrich[fst_enrichment$enrich$padj < 0.05,]
# odds ratio8        J 0.00161255 2.355162 1.373683 3.89203 0.0354761 Translation, ribosomal structure and biogenesis







################# META TRANSCRIPTOMICS ################# 
library(DESeq2)
library(IHW)
TPM <- function(counts, lengths){
  rpk = counts / (lengths / 1000)
  scaling_factor = sum(rpk) / 1000000
  return(rpk / scaling_factor)}

# 1. data generation
files = list.files(path = 'data/WWTP/metaT/', pattern = '.tsv$')[3:53]
first_file = read.csv(paste0('data/WWTP/metaT/', files[1]), sep='\t', comment.char = '#')
trans_tab = data.frame(GeneID = first_file$Geneid,Length = first_file$Length)
for (sample in metadata$Sample){
  file = paste(c('data/WWTP/metaT/Bio17-1_NCBI_',sample,'.counts.tsv'), collapse = '')
  count_file = read.csv(file, sep='\t', comment.char = '#')
  counts = as.numeric(count_file[,endsWith(colnames(count_file), paste(c(sample,'.sorted.bam'), collapse = ''))])
  #trans_tab[sample] = TPM(counts, trans_tab$Length)}
  trans_tab[sample] = counts}
write.csv(trans_tab, file = 'data/WWTP/WWTP_metaT.csv', row.names = F, quote = F)
trans_tab = read.csv('data/WWTP/WWTP_metaT.csv')

counts = trans_tab[,colnames(trans_tab) %in% metadata$Sample]
rownames(counts) = trans_tab$GeneID
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ time_diff_norm)
dds <- DESeq(dds)
res <- as.data.frame(results(dds, filterFun=ihw))
res$Group = 'NotDEG'
res$Group[(res$log2FoldChange > 0) & (res$padj < 0.05)] = 'PosDEG'
res$Group[(res$log2FoldChange < 0) & (res$padj < 0.05)] = 'NegDEG'

ggplot(res, aes(x=log2FoldChange,y=-log(pvalue),color=Group,size=baseMean)) + geom_point() + scale_color_jco() + theme_minimal()
ggsave('figures/deseq2_DEG.pdf', width=7,height = 5)

pos_deg_enrich = CalcEnrichment(gff, rownames(res)[res$Group=='PosDEG'])
neg_deg_enrich = CalcEnrichment(gff, rownames(res)[res$Group=='NegDEG'])

res$PolyCorr_cor = NA
res$PolyCorr_padj = NA
res$PolyCorr_group = NA
for (gene in rownames(res)){
  if(gene %in% mt_polycorr$pi_corr_res$gene_id){
  pi_corr_coef = mt_polycorr$pi_corr_res$cor[mt_polycorr$pi_corr_res$gene_id == gene]
  pi_corr_padj = mt_polycorr$pi_corr_res$padj[mt_polycorr$pi_corr_res$gene_id == gene]
  if (!(is.na(pi_corr_padj))){
  res$PolyCorr_cor[rownames(res) == gene] = pi_corr_coef
  res$PolyCorr_padj[rownames(res) == gene] = pi_corr_padj
  if (pi_corr_padj < 0.05){
    if (pi_corr_coef > 0){res$PolyCorr_group[rownames(res) == gene] = 'PosAssoc'}
    else {res$PolyCorr_group[rownames(res) == gene] = 'NegAssoc'}}
  else{res$PolyCorr_group[rownames(res) == gene] = 'NotAssoc'}}}}

table(res$Group,res$PolyCorr_group)
#         NegAssoc  NotAssoc  PosAssoc
# NegDEG         6        8       13
# NotDEG       920      686     1888
# PosDEG        29       20       93


ggplot(res, aes(x=Group,y=PolyCorr_cor, color=Group)) + geom_boxplot() + stat_compare_means(label.y = 2.3) + 
  stat_compare_means(comparisons = list(c('PosDEG','NegDEG'), c('PosDEG', 'NotDEG'), c('NegDEG', 'NotDEG')) ,label.y = c(1.6,1.8,2))  + 
  scale_color_jco() + theme_minimal() + xlab('') + ylab('PolyCorr') + theme(legend.position = 'none')
ggsave('figures/metaT_polycorr_comp.pdf',width = 3, height = 4)


################# POPULATION DIVERGENCE ################# 
# time-series analysis
library(gtools)
library(igraph)
library(lubridate)
samples_vec = samples_vec[grepl('D',samples_vec)]
poly_net = PolyNet(data_mt, samples_vec, mt_poly, 6, 3)

library(dplyr)
net_data = as.data.frame(poly_net %>% group_by(s1,s2) %>% summarise(fst=mean(FST, na.rm=T)))
metadata$yday = yday(metadata$Date)

#net_data$s1_yday = vapply(net_data$s1, function(x) metadata$yday[metadata$Sample==x], FUN.VALUE = double(1))
#net_data$s2_yday = vapply(net_data$s2, function(x) metadata$yday[metadata$Sample==x], FUN.VALUE = double(1))
#net_data$diff_yday = abs(net_data$s1_yday - net_data$s2_yday)
#net_data$shift = 'No'
#net_data$shift[net_data$s1 %in% c('D20','D05','D15')] = 'yes'
#net_data$shift[net_data$s2 %in% c('D20','D05','D15')] = 'yes'
#ggplot(net_data, aes(x=diff_yday, y=fst, color=shift)) + geom_point() + geom_smooth()

graph_mat <- graph.data.frame(net_data, directed=FALSE)
graph_mat = as_adjacency_matrix(graph_mat, attr="fst", sparse=FALSE)

# graph creation
graph = graph_from_adjacency_matrix(graph_mat, mode="undirected", weighted=TRUE, diag=F)
hist(E(graph)$weight)
E(graph)$weight = 1-E(graph)$weight
hist(E(graph)$weight)

# filtering
#hist(E(graph)$weight)
graph=delete.edges(graph, which(E(graph)$weight < 0.99))

# colours
V(graph)$period <- vapply(V(graph)$name, function(x) metadata$Variable[metadata$Sample == x], FUN.VALUE = character(1))
V(graph)$color = 'darkgrey'
V(graph)$color[V(graph)$period == 'during'] = 'steelblue'
V(graph)$color[V(graph)$period == 'after'] = 'turquoise'

# shapes
V(graph)$season <- vapply(V(graph)$name, function(x) metadata$Season[metadata$Sample == x], FUN.VALUE = character(1))
V(graph)$shape <- "crectangle"
V(graph)$shape[V(graph)$season == 'Autumn'] <- "circle"
V(graph)$shape[V(graph)$season == 'Winter'] <- "csquare"
V(graph)$shape[V(graph)$season == 'Spring'] <- "square"

E(graph)$weight = (E(graph)$weight-mean(E(graph)$weight))/sd(E(graph)$weight)
E(graph)$weight = log(E(graph)$weight + abs(min(E(graph)$weight)) + 1)
# plot with season as shape and colour as period
plot(graph, edge.width = E(graph)$weight, vertex.shape = V(graph)$shape, layout = layout_nicely, vertex.size = 10)


####### CLUSTERING
c1 = cluster_leading_eigen(graph)
modularity(c1)
V(graph)$cluster = membership(c1)


V(graph)$season <- vapply(V(graph)$name, function(x) metadata$Season[metadata$Sample == x], FUN.VALUE = character(1))
V(graph)$color <- "tomato"
V(graph)$color[V(graph)$cluster == 1] <- "steelblue"
V(graph)$color[V(graph)$cluster == 2] <- "lightgreen"
V(graph)$color[V(graph)$cluster == 3] <- "forestgreen"

pdf('figures/fst_net.pdf')
plot(graph, edge.width = E(graph)$weight, vertex.color=V(graph)$color, layout = layout_nicely, vertex.size = 10)
dev.off()

metadata$cluster_graph = vapply(metadata$Sample, function(x) ifelse(x %in% c('D05','D15'), '', as.character(membership(c1)[names(membership(c1)) == x])), FUN.VALUE = character(1))
metadata$cluster_graph[metadata$cluster_graph == ''] = NA
ggplot(na.omit(metadata), aes(x=Date,y=cluster_graph,color=cluster_graph)) + geom_point(size=5) + xlab('') + ylab('Cluster') +
  scale_color_manual(values=c('steelblue','lightgreen','forestgreen','tomato')) + theme_minimal() + theme(legend.position = 'none')
ggsave('figures/clusters_seasons.pdf', width = 4,height = 2)




















#####################
merged_fst = data.frame(gene_id=intersect(mt_fst_ba_norm$gene_id,mt_fst_bd_norm$gene_id))
merged_fst$FST_ba = vapply(merged_fst$gene_id, function(x) mt_fst_ba_norm$FST[mt_fst_ba_norm$gene_id == x], FUN.VALUE = numeric(1))
merged_fst$FST_bd = vapply(merged_fst$gene_id, function(x) mt_fst_bd_norm$FST[mt_fst_bd_norm$gene_id == x], FUN.VALUE = numeric(1))
merged_fst$FST_da = vapply(merged_fst$gene_id, function(x) mt_fst_da_norm$FST[mt_fst_da_norm$gene_id == x], FUN.VALUE = numeric(1))

merged_fst$Log2FC = vapply(merged_fst$gene_id, function(x) res$log2FoldChange[rownames(res) == x], FUN.VALUE = numeric(1))
merged_fst$BaseMean = vapply(merged_fst$gene_id, function(x) res$baseMean[rownames(res) == x], FUN.VALUE = numeric(1))
merged_fst$padj_DE = vapply(merged_fst$gene_id, function(x) res$padj[rownames(res) == x], FUN.VALUE = numeric(1))
merged_fst$PolyCorr = vapply(merged_fst$gene_id, function(x) res$PolyCorr_cor[rownames(res) == x], FUN.VALUE = numeric(1))
merged_fst$PolyP = vapply(merged_fst$gene_id, function(x) res$PolyCorr_padj[rownames(res) == x], FUN.VALUE = numeric(1))
merged_fst$Description = vapply(merged_fst$gene_id, function(x) gff$V9[gff$gene == x], FUN.VALUE = character(1))

merged_fst$col_bd = 'no'
merged_fst$col_bd[merged_fst$FST_bd > quantile(merged_fst$FST_bd, probs=0.9)] = 'yes'
merged_fst$col_ba = 'no'
merged_fst$col_ba[merged_fst$FST_ba > quantile(merged_fst$FST_ba, probs=0.9)] = 'yes'
merged_fst$col_da = 'no'
merged_fst$col_da[merged_fst$FST_da > quantile(merged_fst$FST_da, probs=0.9)] = 'yes'

merged_fst$POS = vapply(merged_fst$gene_id, function(x) gff$V4[gff$gene == x], FUN.VALUE = numeric(1))
merged_fst$CONTIG = vapply(merged_fst$gene_id, function(x) gff$V1[gff$gene == x], FUN.VALUE = character(1))



p1 = ggplot(merged_fst, aes(x=FST_bd,fill=col_bd)) + geom_histogram() + theme_minimal() + 
  scale_fill_jco() +  theme(legend.position = 'none', strip.text.x = element_blank()) + xlab('') + ylab('Divergence (Fst)')
p3 = ggplot(merged_fst, aes(x=POS,y=PolyCorr)) + geom_point(aes(color = col_bd),alpha=0.8) + geom_hline(yintercept = quantile(merged_fst$FST_bd, probs=0.9), color='darkgrey', linetype='dashed') + 
  facet_grid(~CONTIG,scales = 'free_x', space = 'free_x') + scale_color_jco() + 
  theme_minimal() + theme(legend.position = 'none', strip.text.x = element_blank()) + xlab('') + ylab('PolyCorr')
p2 = ggplot(merged_fst, aes(x=POS,y=FST_bd)) + geom_point(aes(color = col_bd),alpha=0.8) + geom_hline(yintercept = quantile(merged_fst$FST_bd, probs=0.9), color='darkgrey', linetype='dashed') + 
  facet_grid(~CONTIG,scales = 'free_x', space = 'free_x') + scale_color_jco() + 
  theme_minimal() + theme(legend.position = 'none', strip.text.x = element_blank()) + xlab('') + ylab('Divergence (Fst)')
ggarrange(p1,p2,p3,nrow = 3)
ggsave('figures/poly_corr_fst_manhattan.pdf', width = 10, height = 10)

merged_fst$col_type = 'no'
merged_fst$col_type[merged_fst$gene_id %in% merged_fst$gene_id[merged_fst$FST_ba > 0.1]] = 'yes'

p1 = ggplot(merged_fst, aes(x=PolyCorr,y=FST_bd)) + geom_point(aes(color=col_type)) + geom_smooth(color='black') + scale_color_jco() + theme_minimal() + 
  xlab('Polycorr') + ylab('Divergence (Fst)') + theme(legend.position = 'none')
p2 = ggplot(merged_fst, aes(x=PolyCorr,y=FST_ba)) + geom_point(aes(color=col_type)) + geom_smooth(color='black') + scale_color_jco() + theme_minimal() + 
  xlab('Polycorr') + ylab('Divergence (Fst)') + theme(legend.position = 'none')
p3 = ggplot(merged_fst, aes(x=FST_bd,y=FST_ba)) + geom_point(aes(color=col_type)) + geom_smooth(color='black') + scale_color_jco() + theme_minimal() + 
  xlab('Divergence during (Fst)') + ylab('Divergence after (Fst)') + theme(legend.position = 'none') + geom_abline(intercept = 0, slope = 1, color = 'darkgrey', linetype = 'dashed')
ggarrange(p1,p2,p3,labels = c('Before/During', 'Before/After', 'During/After comparison'))
ggsave('figures/poly_corr_divergence.pdf', width = 5, height = 4)

# link with DE analysis
merged_fst$DEG = 'NotDEG'
merged_fst$DEG[(merged_fst$Log2FC > 0) & (merged_fst$padj_DE < 0.05)] = 'PosDEG'
merged_fst$DEG[(merged_fst$Log2FC < 0) & (merged_fst$padj_DE < 0.05)] = 'NegDEG'
ggplot(merged_fst, aes(x=DEG,y=FST_bd, color=DEG)) + geom_boxplot() + stat_compare_means() + 
  stat_compare_means(comparisons = list(c('PosDEG','NegDEG'), c('PosDEG', 'NotDEG'), c('NegDEG', 'NotDEG')))  + 
  scale_color_jco() + theme_minimal() + xlab('') + ylab('Divergence(Fst)') + theme(legend.position = 'none')
ggsave('figures/metaT_fst_comp.pdf',width = 3, height = 4)

merged_fst[merged_fst$FST_bd > 0.14,]


plot_dendrogram(c1)
S = cor(graph_mat, method="pearson")
D = 1-S
d = as.dist(D)
cc = hclust(d, method = "average")
plot(cc)
clusters.list = rect.hclust(cc, k = 4, border="blue")
