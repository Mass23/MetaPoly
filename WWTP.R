setwd('~/Documents/PhD/Others/MetaPoly/MetaPoly')

library(devtools)
install_github('https://github.com/Mass23/MetaPoly')
library(MetaPoly)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(IHW)

vcf = vcfR::read.vcfR('data/WWTP/Bio17-1_NCBI_filtered.bcf.gz')
genome = ape::read.dna('data/WWTP/Bio17-1_NCBI.fa', format = "fasta")
gff <- read.delim("data/WWTP/Bio17-1.gff", header=F, comment.char="#", sep='\t', quote = '')
gff = gff[gff$V3 == 'CDS',]

metadata = read.csv('data/WWTP/WWTP_samples.csv')
metadata = metadata[grepl('D',metadata$Sample),]
samples_vec = metadata$Sample
metadata$Date = as.Date(metadata$Date, format = "%d-%m-%Y")
metadata$time_diff = abs(difftime(metadata$Date, as.Date('2011-11-23', tryFormats = "%Y-%m-%d"), units='days'))
names(samples_vec) = log(as.numeric(metadata$time_diff)+1)

# Load data, get polymorphism summary
data_mt = GetGenesData(gff, vcf)

################# POLYSUMMARY ################# 
# 1. data generation
mt_poly = PolySummary(data_mt, samples_vec[grepl('D',samples_vec)])
write.csv(mt_poly, file = 'data/WWTP/WWTP_PolySummary.csv', row.names = F, quote = F)
mt_poly = read.csv('data/WWTP/WWTP_PolySummary.csv')
# generate sample summaries
mt_samples = SummariseSamples(mt_poly, 5)

# 2. data visualisation
mt_samples$table$season = vapply(mt_samples$table$sample, function(x) metadata$Season[metadata$Sample == x], FUN.VALUE = character(1))
mt_samples$table$date = as.Date(vapply(mt_samples$table$sample, function(x) as.character(metadata$Date[metadata$Sample == x]), FUN.VALUE = character(1)))

mt_samples$table$group = 'Baseline'
mt_samples$table$group[mt_samples$table$sample %in% c('D05','D15')] = 'Shift'

p1 = ggplot() + geom_point(mt_samples$table, mapping = aes(x=log(MEAN_DEPTH),y=MEAN_SNP_DEN,color=season,shape=group), size=5, alpha=0.7) + 
  geom_smooth(mt_samples$table[mt_samples$table$group == 'Baseline',], mapping = aes(x=log(MEAN_DEPTH),y=MEAN_SNP_DEN,color=season), method='lm',se=F,fullrange = T) + 
  xlab('') + ylab('SNP Den.') + scale_color_jco() + theme_minimal() + theme(axis.text.x = element_blank())

p2 = ggplot() + geom_point(mt_samples$table, mapping = aes(x=log(MEAN_DEPTH),y=MEAN_PI,color=season,shape=group), size=5, alpha=0.7) + 
  geom_smooth(mt_samples$table[mt_samples$table$group == 'Baseline',], mapping = aes(x=log(MEAN_DEPTH),y=MEAN_PI,color=season), method='lm',se=F,fullrange = T) + 
  xlab('Log Depth') + ylab('Nuc. Diversity') + scale_color_jco() + theme_minimal()
ggarrange(p1,p2,nrow = 2,align='v',  common.legend = TRUE, legend = 'right')
ggsave('WWTP_poly_summary.pdf', width=4,height = 6)

# 3. data analysis
mod_snp_n_all = lm(data=mt_samples$table, MEAN_SNP_DEN ~ log(MEAN_DEPTH):season)
summary(mod_snp_n_all) # all seasons:log(DEPTH) interactions p < 0.0001, r2 = 0.6721
mod_snp_n_shift = lm(data=mt_samples$table[mt_samples$table$season == 'Autumn',], MEAN_SNP_DEN ~ log(MEAN_DEPTH) + group)
summary(mod_snp_n_shift) # shift p = 3.31e-06, r2 = 0.9983
mod_ndiv_all = lm(data=mt_samples$table, MEAN_PI ~ log(MEAN_DEPTH):season)
summary(mod_ndiv_all) # all seasons:log(DEPTH) interactions p < 0.0001, r2 = 0.6721
mod_ndiv_shift = lm(data=mt_samples$table[mt_samples$table$season == 'Autumn',], MEAN_PI ~ log(MEAN_DEPTH) + group)
summary(mod_ndiv_shift) # shift p = 0.6, r2 = 0.9631

################# POLYCORR ################# 
# 1. data generation
mt_polycorr = PolyCorr(mt_poly, 9, samples_vec, 'fdr')
PlotPolyCorr(mt_polycorr$pi_corr_res, mt_polycorr$coefs, 'Microthrix_WWTP', boolean_var = F)

cog_func = read.csv('data/cog-20.def.tab',sep = '\t',header = F)
gff$gene = vapply(gff$V9, function(x) strsplit(strsplit(x,';')[[1]][1],'ID=')[[1]][2], FUN.VALUE = character(1))
gff_to_test = gff[gff$gene %in% mt_polycorr$pi_corr_res$gene_id,]
pos_enrich = CalcEnrichment(gff_to_test, mt_polycorr$pos_genes)
neg_enrich = CalcEnrichment(gff_to_test, mt_polycorr$neg_genes)

mt_polycorr$pi_corr_res$sign = 'No'
mt_polycorr$pi_corr_res$sign[mt_polycorr$pi_corr_res$padj < 0.05] = 'Yes'
ggplot(mt_polycorr$pi_corr_res, aes(x=cor,y=-log(p),color=sign)) + geom_point() + theme_minimal() + scale_color_jco()

# 2. data analysis
pos_enrich$enrich[(pos_enrich$enrich$padj < 0.05),]

neg_enrich$enrich[(neg_enrich$enrich$padj < 0.05) & (neg_enrich$enrich$OR > 1),]
#              Function            p       OR   low_CI   high_CI         padj                                                Function_long
# odds ratio1         C 6.697720e-04 1.967523 1.321758 2.898294 1.138612e-02                             Energy production and conversion
# odds ratio3         E 1.515654e-04 2.004957 1.384084 2.879182 2.728178e-03                          Amino acid transport and metabolism
# odds ratio4         F 7.402115e-07 4.004115 2.250095 7.170449 1.554444e-05                          Nucleotide transport and metabolism
# odds ratio7         I 6.161104e-09 2.970190 2.045444 4.301015 1.355443e-07                               Lipid transport and metabolism
# odds ratio8         J 4.360373e-05 2.091103 1.464532 2.963444 8.720747e-04              Translation, ribosomal structure and biogenesis
# odds ratio10        L 7.988219e-04 2.067540 1.342677 3.147488 1.278115e-02                        Replication, recombination and repair
# odds ratio15        Q 5.924271e-05 3.803353 1.900858 7.665006 1.125611e-03 Secondary metabolites biosynthesis, transport and catabolism

# 3. data visualisation
ggplot(neg_enrich$enrich[(neg_enrich$enrich$padj < 0.05) & (neg_enrich$enrich$OR > 1),]) + geom_point(aes(x=OR,y=Function_long,color=Function_long),size=10) + geom_errorbarh(aes(xmin=low_CI,xmax=high_CI,y=Function_long,color=Function_long),size=2) + 
  xlab('Odds ratio') + ylab('') + theme_minimal() + scale_color_jco() + theme(legend.position = 'none') + geom_vline(xintercept = 1, linetype='dashed')
ggsave('Microthrix_WWTP_res/Neg_enrich.pdf', width=6, height = 5)

out_sign_neg = cog_func[(cog_func$V1 %in% neg_enrich$cogs) & (cog_func$V2 %in% neg_enrich$enrich$Function[(neg_enrich$enrich$padj < 0.05) & (neg_enrich$enrich$OR > 1)]),]
write.csv(out_sign_neg, file='WWTP_mt_neg_sign_funcs_cogs.csv')



################# META TRANSCRIPTOMICS ################# 
library(DESeq2)
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

counts = trans_tab[,3:53]
rownames(counts) = trans_tab$GeneID
metadata$log_time_diff = log(as.integer(metadata$time_diff) + 1)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ log_time_diff)
dds <- DESeq(dds)
res <- as.data.frame(results(dds, filterFun=ihw))
res$Group = 'NotDEG'
res$Group[(res$log2FoldChange > 0) & (res$padj < 0.05)] = 'PosDEG'
res$Group[(res$log2FoldChange < 0) & (res$padj < 0.05)] = 'NegDEG'

ggplot(res, aes(x=log2FoldChange,y=-log(pvalue),color=Group,size=baseMean)) + geom_point() + scale_color_jco() + theme_minimal()

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

ggplot(res, aes(x=PolyCorr_cor,y=log2FoldChange)) + geom_point() + geom_smooth()
table(res$Group,res$PolyCorr_group)

ggplot(res, aes(x=Group,y=PolyCorr_cor, color=Group)) + geom_boxplot() + stat_compare_means(label.y = 2.2) + 
  stat_compare_means(comparisons = list(c('PosDEG','NegDEG'), c('PosDEG', 'NotDEG'), c('NegDEG', 'NotDEG')) ,label.y = c(1.6,1.8,2))  + 
  scale_color_jco() + theme_minimal()

# 2. data analysis

genes_trans$Group = 'Not assoc.'
genes_trans$Group[genes_trans$gene_id %in% mt_polycorr$pos_genes] = 'Pos. assoc.'
genes_trans$Group[genes_trans$gene_id %in% mt_polycorr$neg_genes] = 'Neg. assoc.'

genes_trans$Group_trans = 'Not DEG'
genes_trans$Group_trans[(genes_trans$TPM_padj < 0.05) & (genes_trans$TPM_cor > 0)] = 'Pos. DEG'
genes_trans$Group_trans[(genes_trans$TPM_padj < 0.05) & (genes_trans$TPM_cor < 0)] = 'Neg. DEG'

ggplot(genes_trans, aes(x=Group,y=TPM_mean, color=Group)) + geom_boxplot() + stat_compare_means(label.y = 14) + 
  stat_compare_means(comparisons = list(c('Pos. assoc.','Neg. assoc.'), c('Pos. assoc.', 'Not assoc.'), c('Neg. assoc.', 'Not assoc.')) ,label.y = c(10,11,12))  + 
  scale_color_jco() + theme_minimal()

ggplot(genes_trans, aes(x=Group,y=TPM_cor, color=Group)) + geom_boxplot() + stat_compare_means(label.y = 1.1) + 
  stat_compare_means(comparisons = list(c('Pos. assoc.','Neg. assoc.'), c('Pos. assoc.', 'Not assoc.'), c('Neg. assoc.', 'Not assoc.')) ,label.y = c(0.65,0.8,0.9))  + 
  scale_color_jco() + theme_minimal()

na.omit(res[(res$PolyCorr_padj < 0.0001) & (res$PolyCorr_cor < 0) & (res$log2FoldChange < 0),])

neg_deg = CalcEnrichment(gff, rownames(res)[res$Group=='NegDEG'])

pos_deg = CalcEnrichment(gff, rownames(res)[res$Group=='PosDEG'])











###################################################################################################
# Compare before - during
samples_vec = metadata$Sample[metadata$Variable %in% c('before','during')]
names(samples_vec) = rep(0,length(samples_vec))
names(samples_vec)[samples_vec %in% metadata$Sample[metadata$Variable == 'during']] = 1

mt_fst_bd = PolyDiv(data_mt, samples_vec)
mt_fst_bd_norm = na.omit(mt_fst_bd[mt_fst_bd$SNP_N >= 10,])
wilcox.test(na.omit(mt_fst_bd_norm$FST), rep(0,nrow(na.omit(mt_fst_bd_norm))), alternative = 'greater')
# no evidence for higher PI between shift and before populations, than within before and shift

###################################################################################################
# Compare before - after
samples_vec = metadata$Sample[metadata$Season %in% c('Spring','Winter')]
names(samples_vec) = rep(0,length(samples_vec))
names(samples_vec)[samples_vec %in% metadata$Sample[metadata$Season == 'Winter']] = 1

mt_fst_ba = PolyDiv(data_mt, samples_vec)
mt_fst_ba_norm = na.omit(mt_fst_ba[mt_fst_ba$SNP_N >= 10,])
wilcox.test(na.omit(mt_fst_ba_norm$FST), rep(0,nrow(na.omit(mt_fst_ba_norm))), alternative = 'greater')
