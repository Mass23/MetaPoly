setwd('~/Documents/PhD/Others/MetaPoly/MetaPoly')

library(devtools)
install_github('https://github.com/Mass23/MetaPoly')
library(MetaPoly)
library(ggplot2)
library(ggsci)
library(ggpubr)

vcf = vcfR::read.vcfR('data/GFS/metabat_res.836_filtered.bcf.gz')
genome = ape::read.dna('data/GFS/metabat_res.836.fa', format = "fasta")
gff <- read.delim("data/GFS/PROKKA_09242021.gff", header=F, comment.char="#", sep='\t', quote = '')
gff = gff[gff$V3 == 'CDS',]

metadata = read.csv('data/GFS/samples.csv')
samples_vec_n = colnames(vcf@gt)
samples_vec_n = samples_vec_n[samples_vec_n != 'FORMAT']
format_name <- function(name){
  name = gsub('GL_R.._','',name)
  name = gsub('GL_R._','',name)
  name = gsub('_.$','',name)
  name = gsub('B$','',name)
  name = gsub('UP','Up',name)
  name = gsub('DN','Down',name)
  name[!(grepl('GL_',name))] = gsub('GL','GL_',name[!(grepl('GL_',name))])
  return(name)}
names_samp_vec = vapply(samples_vec_n, function(x) log(metadata$snout_dist[metadata$location == x]), FUN.VALUE = numeric(1))

samples_vec = colnames(vcf@gt)[colnames(vcf@gt) != 'FORMAT']
names(samples_vec) = names_samp_vec

# Load data, get polymorphism summary
data_mt = GetGenesData(gff, vcf)

################# POLYSUMMARY ################# 
mt_poly = PolySummary(data_mt, samples_vec)
mt_poly$habitat = 'Sediment'
mt_poly$habitat[startsWith(mt_poly$sample, 'GL_R')] = 'Rock'

min_depth_comp = data.frame()
for (i in 1:10){
  mt_samples = SummariseSamples(mt_poly, i)
  if(nrow(mt_samples$table) > 0){
  mt_samples$table$MIN_DEPTH = i
  mt_samples$table$N_genes = mt_samples$n_genes
  min_depth_comp = rbind(min_depth_comp, mt_samples$table)}}
min_depth_comp$habitat = 'Sediment'
min_depth_comp$habitat[startsWith(min_depth_comp$sample, 'GL_R')] = 'Rock'
min_depth_comp$sample_clean = vapply(min_depth_comp$sample, format_name, FUN.VALUE = character(1))
min_depth_comp$expedition = vapply(min_depth_comp$sample_clean, function(x) metadata$expedition[metadata$location == x], FUN.VALUE = character(1))

ggplot(min_depth_comp,aes(x=MIN_DEPTH,y=log(MEAN_PI),group=sample,color=expedition,shape=habitat)) + geom_line() + geom_point(size=3) + geom_text(aes(label=sample_clean),nudge_y = 0.01)


GI_samples = SummariseSamples(mt_poly, 10)
GI_samples$table$habitat = 'Sediment'
GI_samples$table$habitat[startsWith(GI_samples$table$sample, 'GL_R')] = 'Rock'
GI_samples$table$sample_clean = vapply(GI_samples$table$sample, format_name, FUN.VALUE = character(1))
GI_samples$table$snout_dist = as.numeric(vapply(GI_samples$table$sample_clean, function(x) metadata$snout_dist[metadata$location == x], FUN.VALUE = numeric(1)))
GI_samples$table$gl_a = as.numeric(vapply(GI_samples$table$sample_clean, function(x) metadata$gl_a[metadata$location == x], FUN.VALUE = numeric(1)))
GI_samples$table$log_GI = log(sqrt(GI_samples$table$gl_a) / (sqrt(GI_samples$table$gl_a) + GI_samples$table$snout_dist))
GI_samples$table$expedition = vapply(GI_samples$table$sample_clean, function(x) metadata$expedition[metadata$location == x], FUN.VALUE = character(1))

mt_samples$turb = vapply(mt_samples$sample_clean, function(x) metadata$turb[metadata$location == x], FUN.VALUE = numeric(1))
mt_samples$water_temp = vapply(mt_samples$sample_clean, function(x) metadata$water_temp[metadata$location == x], FUN.VALUE = numeric(1))
mt_samples$glacier = vapply(mt_samples$sample_clean, function(x) as.character(strsplit(x,split='_')[[1]][2]), FUN.VALUE = character(1))



ggplot(GI_samples$table, aes(y=MEAN_SNP_DEN,x=log_GI,color=expedition,size=log(MEAN_DEPTH))) + geom_point() + geom_smooth(method='lm') + facet_grid(~habitat)
summary(lm(MEAN_PI ~ log(MEAN_DEPTH) + habitat, data=mt_samples))

ggplot(GI_samples$table, aes(y=MEAN_PI,x=log_GI,color=expedition,size=log(MEAN_DEPTH))) + geom_point() + geom_smooth(method='lm') + facet_grid(~habitat)
summary(lm(MEAN_PI ~ log(MEAN_DEPTH) + habitat, data=mt_samples))

ggplot(GI_samples$table, aes(y=MEAN_PI,x=MEAN_SNP_DEN,color=habitat,size=log(MEAN_DEPTH))) + geom_point() + geom_smooth(method='lm')



# models
mod_snp_n_all = lm(data=mt_samples, snp_den ~ log(depth):season)
summary(mod_snp_n_all)
mod_snp_n_shift = lm(data=mt_samples[mt_samples$season == 'Autumn',], snp_den ~ log(depth) + group)
summary(mod_snp_n_shift)

mod_ndiv_all = lm(data=mt_samples, pi ~ log(depth):season)
summary(mod_ndiv_all)
mod_ndiv_shift = lm(data=mt_samples[mt_samples$season == 'Autumn',], pi ~ log(depth) + group)
summary(mod_ndiv_shift)

################# POLYCORR ################# 
mt_polycorr = PolyCorr(mt_poly, 9, samples_vec, 'fdr')
PlotPolyCorr(mt_polycorr$pi_corr_res, mt_polycorr$coefs, 'Microthrix_WWTP', boolean_var = F)



gff$gene = vapply(gff$V9, function(x) strsplit(strsplit(x,';')[[1]][1],'ID=')[[1]][2], FUN.VALUE = character(1))
gff_to_test = gff[gff$gene %in% mt_polycorr$pi_corr_res$gene_id,]
pos_enrich = CalcEnrichment(gff_to_test, mt_polycorr$pos_genes, cog_func)
neg_enrich = CalcEnrichment(gff_to_test, mt_polycorr$neg_genes, cog_func)

mt_polycorr$pi_corr_res$sign = 'No'
mt_polycorr$pi_corr_res$sign[mt_polycorr$pi_corr_res$padj < 0.05] = 'Yes'
ggplot(mt_polycorr$pi_corr_res, aes(x=cor,y=-log(p),color=sign)) + geom_point() + theme_minimal() + scale_color_jco()

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
ggplot(neg_enrich$enrich[(neg_enrich$enrich$padj < 0.05) & (neg_enrich$enrich$OR > 1),]) + geom_point(aes(x=OR,y=Function_long,color=Function_long),size=10) + geom_errorbarh(aes(xmin=low_CI,xmax=high_CI,y=Function_long,color=Function_long),size=2) + 
  xlab('Odds ratio') + ylab('') + theme_minimal() + scale_color_jco() + theme(legend.position = 'none') + geom_vline(xintercept = 1, linetype='dashed')
ggsave('Microthrix_WWTP_res/Neg_enrich.pdf', width=6, height = 5)

out_sign_neg = cog_func[(cog_func$V1 %in% neg_enrich$cogs) & (cog_func$V2 %in% neg_enrich$enrich$Function[(neg_enrich$enrich$padj < 0.05) & (neg_enrich$enrich$OR > 1)]),]
write.csv(out_sign_neg, file='WWTP_mt_neg_sign_funcs_cogs.csv')

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











