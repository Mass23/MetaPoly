library(devtools)
install_github('https://github.com/Mass23/MetaPoly')
library(MetaPoly)

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
mt_poly = PolySummary(data_mt, samples_vec[grepl('D',samples_vec)])
mt_polycorr = PolyCorr(mt_poly, 4, samples_vec)




mt_poly_f$CONS_INDEX = ((mt_poly_f$gene_length - mt_poly_f$SNP_N)/mt_poly_f$gene_length) + (mt_poly_f$SNP_N/mt_poly_f$gene_length*mt_poly_f$MAJF)


library(bayestestR)
library(rstanarm)
library(ggplot2)
mt_poly_f$EVENNESS = mt_poly_f$EVENNESS * 0.99999

df_mod_lm = data.frame()
for (gene in unique(mt_poly_f$gene_id)){
gene_data = na.omit(mt_poly_f[mt_poly_f$gene_id==gene,])
if(nrow(gene_data) > 9){
mod = lm(EVENNESS ~ log(DEPTH) + variable, weights = sqrt(SNP_N), data=gene_data)
mod_sum = summary(mod)
coef = mod_sum$coefficients[3,1]
p_val = mod_sum$coefficients[3,4]
df_mod_lm = rbind(df_mod_lm, data.frame(gene=gene, coef=coef, p=p_val))}}

hist(df_mod_lm$coef)
hist(df_mod_lm$p)


coefs = mod$coefficients
coefs = coefs[startsWith(names(coefs),'sample')]
names(coefs) = vapply(names(coefs), function(x) strsplit(as.character(x), split='mple')[[1]][2], FUN.VALUE = character(1))
metadata = metadata[metadata$Sample %in% names(coefs),]
metadata$Even_mod = vapply(metadata$Sample, function(x) coefs[names(coefs) == x], FUN.VALUE = numeric(1))

ggplot(metadata, aes(x=Date,y=Even_mod,color=Season)) + geom_point()


library(pscl)
library(ggplot2)
snp_n_mod <- 
summary(snp_n_mod)

coefs = snp_n_mod$coefficients$count
coefs = coefs[startsWith(names(coefs),'as.factor(sample)')]
names(coefs) = vapply(names(coefs), function(x) strsplit(as.character(x), split=')')[[1]][2], FUN.VALUE = character(1))
metadata = metadata[metadata$Sample %in% names(coefs),]
metadata$SNP_N_COEFS = vapply(metadata$Sample, function(x) coefs[names(coefs) == x], FUN.VALUE = numeric(1))

ggplot(metadata, aes(x=Date,y=SNP_N_COEFS,color=Season)) + geom_point()


ggplot(mt_poly_f, aes(x=DEPTH,y=EVENNESS)) + geom_point() + geom_smooth()
