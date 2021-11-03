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

samples_vec = metadata$Sample
names(samples_vec) = metadata$Test

mt_poly = PolySummary(data_mt, samples_vec[grepl('D',samples_vec)])
mt_polycorr = PolyCorr(mt_poly, 4, samples_vec)
mt_fst = PolyDiff(data_mt, samples_vec)

gff$gene = vapply(gff$V9, function(x) strsplit(strsplit(x,';')[[1]][1],'ID=')[[1]][2], FUN.VALUE = character(1))
gff$cog = vapply(gff$V9, function(x) strsplit(strsplit(x,'COG:')[[1]][2],';')[[1]][1], FUN.VALUE = character(1))
cog_func = read.table('data/cog-20.def.tab', sep='\t')
gff$cog_f = vapply(gff$cog, function(x) ifelse(is.na(x), 'NoCOG', substr(cog_func$V2[cog_func$V1 == x],1,1)), FUN.VALUE = character(1))

mt_fst$COG = vapply(mt_fst$gene_id, function(x) gff$cog_f[gff$gene == x], FUN.VALUE = character(1))
mt_fst$COG[is.na(mt_fst$COG)] = 'NoCOG'
mt_fst$COG_long = vapply(mt_fst$COG, function(x) ifelse(x == 'NoCOG', 'NoCOG', cog_functions[names(cog_functions) == x]), FUN.VALUE =  character(1))
mt_fst$fst_norm = mt_fst$FST
mt_fst$fst_norm[mt_fst$fst_norm < 0] = 0
ggplot(mt_fst, aes(x=fst_norm,y=COG_long)) + geom_boxplot()

mod = zoib::zoib(model=fst_norm ~ COG_long|1|COG_long|1|1, data=mt_fst)

pvals = c()
for (sample in unique(test$Var1)){
  for (snp in unique(test$SNP)){
    pvals = c(pvals, cor.test(test$PI[(test$Var2 == sample) & (test$SNP == snp)], test$DEPTH[(test$Var2 == sample)  & (test$SNP == snp)])$p.value)
  }
}

x = glm(formula = PI ~ as.factor(COMP) + as.factor(SNP) + as.factor(Var1) + as.factor(Var2) + DEPTH, data=test, family = quasibinomial(), na.action = na.omit)



hist(mt_fst$fst)

mt_poly$CONS_INDEX = ((mt_poly$gene_length - mt_poly$SNP_N)/mt_poly$gene_length) + (mt_poly$SNP_N/mt_poly$gene_length*mt_poly$MAJF)


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
