library(devtools)
install_github('https://github.com/Mass23/MetaPoly')
library(MetaPoly)
library(ggplot2)
library(ggsci)
library(ggpubr)

setwd('~/Documents/PhD/Others/MetaPoly/MetaPoly')
vcf = vcfR::read.vcfR('data/WWTP/Bio17-1_NCBI_filtered.bcf.gz')
genome = ape::read.dna('data/WWTP/Bio17-1_NCBI.fa', format = "fasta")
gff <- read.delim("data/WWTP/Bio17-1_NCBI_prokka.gff", header=F, comment.char="#", sep='\t', quote = '')
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
mt_poly = PolySummary(data_mt, samples_vec[grepl('D',samples_vec)])
mt_poly$CONS_INDEX = ((mt_poly$gene_length - mt_poly$SNP_N)/mt_poly$gene_length) + ((mt_poly$SNP_N/mt_poly$gene_length)*mt_poly$MAJF)
sample_df = data.frame()
for (sample in unique(mt_poly$sample)){
  mean_cons = weighted.mean(mt_poly$CONS_INDEX[mt_poly$sample == sample], mt_poly$gene_length[mt_poly$sample == sample], na.rm=T)
  mean_snp_den =  weighted.mean(mt_poly$SNP_N[mt_poly$sample == sample] / mt_poly$gene_length[mt_poly$sample == sample], mt_poly$gene_length[mt_poly$sample == sample], na.rm=T)
  mean_even =  weighted.mean(mt_poly$EVENNESS[mt_poly$sample == sample], mt_poly$gene_length[mt_poly$sample == sample], na.rm=T)
  mean_depth=  weighted.mean(mt_poly$DEPTH[mt_poly$sample == sample], mt_poly$gene_length[mt_poly$sample == sample], na.rm=T)
  sample_df = rbind(sample_df, data.frame(sample=sample,cons=mean_cons,snp_den=mean_snp_den,even=mean_even,depth=mean_depth))}

sample_df$time = as.Date(vapply(sample_df$sample, function(x) as.character(metadata$Date[metadata$Sample == x]), FUN.VALUE = character(1)))
sample_df$season = vapply(sample_df$sample, function(x) metadata$Season[metadata$Sample == x], FUN.VALUE = character(1))
sample_df$group = 'Baseline'
sample_df$group[sample_df$sample %in% c('D05','D15')] = 'Shift'

p1 = ggplot() + geom_point(sample_df, mapping = aes(x=log(depth),y=cons,color=season,shape=group), size=5) + 
  geom_smooth(sample_df[sample_df$group == 'Baseline',], mapping = aes(x=log(depth),y=cons,color=season), method='lm',se=F,fullrange = T) + 
  xlab('') + ylab('Cons. Index') + scale_color_jco() + theme_minimal() + theme(axis.text.x = element_blank())

p2 = ggplot() + geom_point(sample_df, mapping = aes(x=log(depth),y=snp_den,color=season,shape=group), size=5) + 
  geom_smooth(sample_df[sample_df$group == 'Baseline',], mapping = aes(x=log(depth),y=snp_den,color=season), method='lm',se=F,fullrange = T) + 
  xlab('') + ylab('SNP Den.') + scale_color_jco() + theme_minimal() + theme(axis.text.x = element_blank())

p3 = ggplot() + geom_point(sample_df, mapping = aes(x=log(depth),y=even,color=season,shape=group), size=5) + 
  geom_smooth(sample_df[sample_df$group == 'Baseline',], mapping = aes(x=log(depth),y=even,color=season), method='lm',se=F,fullrange = T) + 
  xlab('Log Depth') + ylab('AF Evenness') + scale_color_jco() + theme_minimal()
ggarrange(p1,p2,p3,nrow = 4,align='v',  common.legend = TRUE, legend = 'right')
ggsave('WWTP_poly_summary.pdf', width=5,height = 10)

################# POLYCORR ################# 
mt_polycorr = PolyCorr(mt_poly, 9, samples_vec, 'fdr')
PlotPolyCorr(mt_polycorr$pi_corr_res, mt_polycorr$coefs, 'Microthrix_WWTP', boolean_var = F)
cog_func = read.table('data/cog-20.def.tab', sep='\t')
pos_enrich = CalcEnrichment(gff, mt_polycorr$pos_genes, cog_func)
neg_enrich = CalcEnrichment(gff, mt_polycorr$neg_genes, cog_func)

pos_enrich$enrich[pos_enrich$enrich$padj < 0.05,]
#             Function            p        OR    low_CI   high_CI         padj                       Function_long
# odds ratio1        C 1.250414e-05 0.3466012 0.1976789 0.5834878 0.0002625870    Energy production and conversion
# odds ratio3        E 2.364640e-04 0.4127144 0.2427603 0.6795412 0.0044928168 Amino acid transport and metabolism
# odds ratio4        F 2.736090e-03 0.2913363 0.1068114 0.6881908 0.0492496135 Nucleotide transport and metabolism
# odds ratio7        I 4.105441e-05 0.4041271 0.2485946 0.6395266 0.0008210881      Lipid transport and metabolism
ggplot(pos_enrich$enrich[pos_enrich$enrich$padj < 0.05,]) + geom_point(aes(x=OR,y=Function_long,color=Function_long),size=10) + geom_errorbarh(aes(xmin=low_CI,xmax=high_CI,y=Function_long,color=Function_long),size=2) + 
  xlab('Odds ratio') + ylab('') + theme_minimal() + scale_color_jco() + theme(legend.position = 'none') + geom_vline(xintercept = 1, linetype='dashed')
ggsave('Microthrix_WWTP_res/Pos_enrich.pdf', width=5, height = 4)

neg_enrich$enrich[neg_enrich$enrich$padj < 0.05,]
#              Function            p       OR   low_CI   high_CI         padj                                                Function_long
# odds ratio1         C 7.416405e-04 2.281996 1.397382  3.664675 1.409117e-02                             Energy production and conversion
# odds ratio4         F 1.801926e-05 4.599847 2.200722  9.699622 3.603852e-04                          Nucleotide transport and metabolism
# odds ratio7         I 9.842083e-10 3.672821 2.405202  5.591781 2.066837e-08                               Lipid transport and metabolism
# odds ratio8         J 2.905802e-03 2.123257 1.283893  3.442975 4.939864e-02              Translation, ribosomal structure and biogenesis
# odds ratio15        Q 2.497434e-03 4.550399 1.551631 13.612134 4.495382e-02 Secondary metabolites biosynthesis, transport and catabolism
ggplot(neg_enrich$enrich[neg_enrich$enrich$padj < 0.05,]) + geom_point(aes(x=OR,y=Function_long,color=Function_long),size=10) + geom_errorbarh(aes(xmin=low_CI,xmax=high_CI,y=Function_long,color=Function_long),size=2) + 
  xlab('Odds ratio') + ylab('') + theme_minimal() + scale_color_jco() + theme(legend.position = 'none') + geom_vline(xintercept = 1, linetype='dashed')
ggsave('Microthrix_WWTP_res/Neg_enrich.pdf', width=6, height = 5)

samples_vec = metadata$Sample
names(samples_vec) = rep(0,length(samples_vec))
names(samples_vec)[samples_vec %in% c('D05','D15')] = 1

mt_fst = PolyDiv(data_mt, samples_vec)
mt_fst_norm = na.omit(mt_fst[mt_fst$SNP_N >= 10,])
wilcox.test(na.omit(mt_fst_norm$FST), rep(0,nrow(na.omit(mt_fst_norm))), alternative = 'greater')
# no evidence for higher PI between shift and baseline populations, than within baselines and shift














